"""modelo_pinn.py (v2) -- Marco 4: Modelo Substituto Preditivo Baseado em Fisica (PINN) da planta termica do Marco 3.
Verdade-base: modelo_termico.py (mesma pasta). Janela de predicao W = 300 s. NAO e Gemeo Digital (sem enlace com hardware/telemetria).
Mudancas em relacao a v1 (todas com flag; checkpoints v1 continuam carregaveis e avaliaveis com as metricas novas):
  1. Perdas em K-EQUIVALENTE (cada residuo dividido pela sensibilidade fisica a 1 K), todas com peso 1 por padrao [HIPOTESE DECLARADA].
     Valores de L NAO sao comparaveis com v1.
  2. hard_ic = True por padrao (condicao inicial exata por construcao).
  3. Fase L-BFGS (lote fixo) apos o Adam; --init continua de um checkpoint.
  4. Amostragem correlacionada (T_j0 - T_raiz0 em [5, 40] K) para o TREINO. Os casos de avaliacao R00..R19 sao os MESMOS de v1 (semente 1234);
     cada caso informa se esta dentro do dominio de treino.
  5. Metricas em K: rms_K, max_abs_C, erro no pico de T_j (negativo = subestima = inseguro), vies, resumo por dominio. rel_L2_var foi removida.
  6. Estimador de Delta_ponta (ajuste por minimos quadrados ponderados sobre uma grade, usando a propria PINN como modelo direto) + avaliacao.
Etiquetas: [HIPOTESE DECLARADA] = escolha minha sem fonte. Nao validado experimentalmente. Obs. (auditoria 22/11/2026): os 14 testes pytest foram executados pelo autor; este docstring estava desatualizado.
"""
import argparse
import json
import math
import time

import numpy as np
import torch
from torch import nn

import modelo_termico as mt

torch.set_default_dtype(torch.float32)

T_REF, DTK, W_S = 320.0, 50.0, 300.0     # escalas de normalizacao [HIPOTESE DECLARADA]
A_PL = 49.89754                          # m2, area de projeto do Marco 3 [USUARIO]
PNAMES = ("f", "q_abs", "q_rest", "tj0", "tr0", "dtip")
LO = torch.tensor([300.0, 90.0, 800.0, 280.0, 265.0, 0.0])      # escala de NORMALIZACAO da entrada (identica a v1)
HI = torch.tensor([1410.0, 145.0, 7500.0, 360.0, 330.0, 15.0])
RNG_ = {"f": (300.0, 1410.0), "q_abs": (90.0, 145.0), "q_rest": (800.0, 7500.0), "tr0": (265.0, 330.0), "dtip": (0.0, 15.0)}
DJR = (5.0, 40.0)                        # T_j0 - T_raiz0 [K] no treino [HIPOTESE DECLARADA; S1 medido: ~12.6 a ~32 K]
W_DEF = (1.0, 1.0, 1.0, 1.0)             # pesos (pde, bc, mod, ic) em K^2


class Phys:
    """Constantes fisicas da planta, lidas de modelo_termico (uma unica fonte)."""

    def __init__(self, a_pl=A_PL, par=None):
        p = mt._merge(par)
        self.p, self.a_pl, self.L = p, a_pl, p["L"]
        self.kd, self.ca = p["k"] * p["delta"], p["c_area"]
        self.es, self.ts = 2.0 * mt.optics(p)[1] * mt.SIGMA, p["t_space"]
        self.R, self.cmod, self.cfl = p["r_int"], p["c_mod"], p["c_fluid"]
        self.ell = a_pl / self.L
        self.s, self.kap = p["n_mod"] / self.ell, p["kappa_dab"]
        self.h = 4.0 * self.es * T_REF ** 3          # W/(m2 K): sensibilidade do PDE a 1 K (radiacao linearizada)


def p_mod_t(f, tj_c, util=1.0):
    """P(f, Tj) por modulo [W]; copia em torch de modelo_termico.ia_p (testada contra ela)."""
    ia = mt.IA
    fv = ia["fv"] * ia["f_nom"]
    v = ia["v_min"] + (ia["v_nom"] - ia["v_min"]) * torch.clamp((f - fv) / (ia["f_nom"] - fv), 0.0, 1.0)
    p_dyn = util * (1.0 - ia["fix"] - ia["leak"]) * ia["tdp"] * (v / ia["v_nom"]) ** 2 * (f / ia["f_nom"])
    i_leak = ia["leak"] * ia["tdp"] / ia["v_nom"] * (v / ia["v_nom"]) ** (ia["kv"] - 1.0) * torch.exp(ia["beta"] * (tj_c - ia["t_ref"]))
    return ia["fix"] * ia["tdp"] + p_dyn + v * i_leak


def norm_p(P):
    return 2.0 * (P - LO) / (HI - LO) - 1.0


def ic_T(xi, P):
    """Perfil inicial do campo: quadratico, inclinacao nula na ponta (medido no solver: desvio <= 0.03 K em regime periodico)."""
    tr = P[:, 4:5]
    tt = tr - P[:, 5:6]
    return tt + (tr - tt) * (1.0 - xi) ** 2


class PINN(nn.Module):
    def __init__(self, width=64, depth=4, hard_ic=True):
        super().__init__()

        def mlp(nin):
            layers, d = [], nin
            for _ in range(depth):
                layers += [nn.Linear(d, width), nn.Tanh()]
                d = width
            layers.append(nn.Linear(d, 1))
            return nn.Sequential(*layers)
        self.fin, self.mod, self.hard_ic = mlp(2 + len(PNAMES)), mlp(1 + len(PNAMES)), hard_ic

    def field(self, xi, tau, P):
        """T(xi, tau) em K. xi, tau: [n,1]; P: [n,6] fisico, na ordem PNAMES."""
        raw = self.fin(torch.cat([2.0 * xi - 1.0, 2.0 * tau - 1.0, norm_p(P)], 1))
        return (ic_T(xi, P) + tau * DTK * raw) if self.hard_ic else (T_REF + DTK * raw)

    def tj(self, tau, P):
        raw = self.mod(torch.cat([2.0 * tau - 1.0, norm_p(P)], 1))
        return (P[:, 3:4] + tau * DTK * raw) if self.hard_ic else (T_REF + DTK * raw)


def grad(y, x):
    return torch.autograd.grad(y, x, grad_outputs=torch.ones_like(y), create_graph=True)[0]


# ---------------- residuos em K-EQUIVALENTE
def res_pde(field, ph, xi, tau, P):
    """(c_a T_t - k delta T_xx + 2 eps sigma (T^4 - Ts^4) - q_abs) / h_rad  [K]."""
    xi = xi.clone().requires_grad_(True)
    tau = tau.clone().requires_grad_(True)
    T = field(xi, tau, P)
    Tt = grad(T, tau) / W_S
    Tx = grad(T, xi) / ph.L
    Txx = grad(Tx, xi) / ph.L
    return (ph.ca * Tt - ph.kd * Txx + ph.es * (T ** 4 - ph.ts ** 4) - P[:, 1:2]) / ph.h


def res_tip(field, ph, tau, P):
    """dT/dxi em xi = 1 [K por unidade de xi]."""
    xi = torch.ones_like(tau).requires_grad_(True)
    return grad(field(xi, tau, P), xi)


def res_root_mod(field, tjfun, ph, tau, P):
    """Raiz: balanco de fluxo dividido por (s/R) [K]. Modulo: c_mod dTj/dt - (P - (Tj - T0)/R), multiplicado por R [K]."""
    tau = tau.clone().requires_grad_(True)
    xi = torch.zeros_like(tau).requires_grad_(True)
    T0 = field(xi, tau, P)
    Tx0, Tt0 = grad(T0, xi) / ph.L, grad(T0, tau) / W_S
    Tj = tjfun(tau, P)
    Tjt = grad(Tj, tau) / W_S
    Pm = p_mod_t(P[:, 0:1], Tj - 273.15)
    r_root = (ph.cfl / ph.ell * Tt0 - ph.kd * Tx0 - ph.s * (Tj - T0) / ph.R - ph.s * ph.kap * Pm - P[:, 2:3] / ph.ell) / (ph.s / ph.R)
    r_mod = (ph.cmod * Tjt - (Pm - (Tj - T0) / ph.R)) * ph.R
    return r_root, r_mod


def res_ic(field, tjfun, xi, P):
    tau0 = torch.zeros_like(xi)
    return field(xi, tau0, P) - ic_T(xi, P), tjfun(tau0, P) - P[:, 3:4]


# ---------------- amostragem e perdas
def rep(P, k):
    return P.repeat_interleave(k, 0)


def sample_P(n, gen):
    """Amostragem de TREINO: T_j0 = T_raiz0 + U(DJR) [correlacionada]."""
    u = torch.rand(n, len(PNAMES), generator=gen)

    def lerp(k, c):
        return RNG_[k][0] + u[:, c:c + 1] * (RNG_[k][1] - RNG_[k][0])
    tr0 = lerp("tr0", 4)
    tj0 = tr0 + DJR[0] + u[:, 3:4] * (DJR[1] - DJR[0])
    return torch.cat([lerp("f", 0), lerp("q_abs", 1), lerp("q_rest", 2), tj0, tr0, lerp("dtip", 5)], 1)


def sample_P_legacy(n, gen):
    """Amostragem de v1 (independente). Usada SO nos casos R00..R19 de avaliacao, para comparabilidade com E1-E3."""
    return LO + torch.rand(n, len(PNAMES), generator=gen) * (HI - LO)


def make_batch(ns, gen, nint=96, nbc=32, nic=32, sampling="corr"):
    """sampling='corr' (padrao v2): T_j0 correlacionada com T_raiz0. sampling='indep': amostragem independente de v1 (comparavel com E1-E3)."""
    if sampling not in ("corr", "indep"):
        raise ValueError("sampling deve ser 'corr' ou 'indep'")
    P = sample_P(ns, gen) if sampling == "corr" else sample_P_legacy(ns, gen)
    return {"P": P, "nint": nint, "nbc": nbc, "nic": nic,
            "xi_i": torch.rand(ns * nint, 1, generator=gen), "ta_i": torch.rand(ns * nint, 1, generator=gen),
            "ta_b": torch.rand(ns * nbc, 1, generator=gen), "xi_c": torch.rand(ns * nic, 1, generator=gen)}


def losses_from_batch(model, ph, b):
    """Retorna (L_pde, L_bc, L_mod, L_ic) em K^2: media de quadrados de residuos em K-equivalente."""
    P = b["P"]
    L_pde = res_pde(model.field, ph, b["xi_i"], b["ta_i"], rep(P, b["nint"])).pow(2).mean()
    Pb = rep(P, b["nbc"])
    r_root, r_mod = res_root_mod(model.field, model.tj, ph, b["ta_b"], Pb)
    L_bc = res_tip(model.field, ph, b["ta_b"], Pb).pow(2).mean() + r_root.pow(2).mean()
    r_f, r_j = res_ic(model.field, model.tj, b["xi_c"], rep(P, b["nic"]))
    return L_pde, L_bc, r_mod.pow(2).mean(), r_f.pow(2).mean() + r_j.pow(2).mean()


def losses(model, ph, ns, gen, nint=96, nbc=32, nic=32):
    return losses_from_batch(model, ph, make_batch(ns, gen, nint, nbc, nic))


def total(parts, w=W_DEF):
    return w[0] * parts[0] + w[1] * parts[1] + w[2] * parts[2] + w[3] * parts[3]


def _rec(fase, it, parts, loss):
    r = {"fase": fase, "it": it, "loss": loss.item()}
    for n, p in zip(("L_pde", "L_bc", "L_mod", "L_ic"), parts):
        r[n] = p.item()
        r["rms_K_" + n[2:]] = math.sqrt(p.item())
    return r


def lbfgs_phase(model, ph, fixed, steps, w=W_DEF, log_every=10, hist=None):
    """L-BFGS (Wolfe forte) em lote FIXO de colocacao. Cada .step() faz ate 20 iteracoes."""
    opt = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=20, history_size=50, line_search_fn="strong_wolfe",
                            tolerance_grad=1e-12, tolerance_change=1e-15)

    def closure():
        opt.zero_grad()
        loss = total(losses_from_batch(model, ph, fixed), w)
        loss.backward()
        return loss
    hist = [] if hist is None else hist
    for k in range(1, steps + 1):
        opt.step(closure)
        if k % log_every == 0 or k == 1 or k == steps:
            parts = losses_from_batch(model, ph, fixed)
            hist.append(_rec("lbfgs", k, parts, total(parts, w)))
            print(json.dumps(hist[-1]), flush=True)
    return hist


def train(epochs=20000, ns=48, lr=1e-3, seed=0, w=W_DEF, width=64, depth=4, hard_ic=True, out="pinn.pt", log_every=250,
          lbfgs_steps=0, ns_lbfgs=96, init=None, sampling="corr"):
    """Adam + cosseno em lotes novos a cada passo, depois L-BFGS opcional em lote fixo. Sementes fixas."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    gen = torch.Generator().manual_seed(seed)
    ph, model = Phys(), PINN(width, depth, hard_ic)
    if init:
        ck = torch.load(init, map_location="cpu")
        width, depth, hard_ic = ck["width"], ck["depth"], ck["hard_ic"]
        model = PINN(width, depth, hard_ic)
        model.load_state_dict(ck["state"])
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1), eta_min=lr * 0.01)
    hist = []
    for it in range(1, epochs + 1):
        opt.zero_grad()
        parts = losses_from_batch(model, ph, make_batch(ns, gen, sampling=sampling))
        loss = total(parts, w)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        opt.step()
        sched.step()
        if it % log_every == 0 or it == 1:
            hist.append(_rec("adam", it, parts, loss))
            print(json.dumps(hist[-1]), flush=True)
    if lbfgs_steps > 0:
        fixed = make_batch(ns_lbfgs, torch.Generator().manual_seed(seed + 1000))
        lbfgs_phase(model, ph, fixed, lbfgs_steps, w, hist=hist)
    if out:
        torch.save({"state": model.state_dict(), "width": width, "depth": depth, "hard_ic": hard_ic, "seed": seed, "version": 2, "hist": hist,
                    "planta_sha256": _planta_sha()}, out)      # impressao digital da planta usada no treino (registrar em contrato_parametros.CKPT_PLANTA_SHA256)
    return model, hist


def _planta_sha():
    from contrato_parametros import plant_fingerprint
    return plant_fingerprint()


def load(path):
    ck = torch.load(path, map_location="cpu")
    m = PINN(ck["width"], ck["depth"], ck["hard_ic"])
    m.load_state_dict(ck["state"])
    m.eval()
    return m


# ---------------- referencia (solver do Marco 3) e avaliacao
def ic_np(xi, tt, tr):
    return tt + (tr - tt) * (1.0 - np.asarray(xi)) ** 2


def solver_window(case, a_pl=A_PL, N=160, dt=1.0, W=W_S):
    """Janela de W s com forcante estatica (f, q_abs, q_rest); 2 meios passos de Euler implicito no inicio (Rannacher), depois CN."""
    f, qa, qr, tj0, tr0, dtip = case
    pl = mt.Plant(a_pl, {"n_nodes": N}, mt.static_drive(f, qa, qr))
    y = np.concatenate(([tj0], ic_np(np.linspace(0.0, 1.0, N + 1), tr0 - dtip, tr0)))
    Y = [y.copy()]
    for k in range(int(round(W / dt))):
        ta, tb = k * dt, (k + 1) * dt
        if k == 0:
            y, _ = pl.step(y, ta, 0.5 * (ta + tb), 1.0)
            y, _ = pl.step(y, 0.5 * (ta + tb), tb, 1.0)
        else:
            y, _ = pl.step(y, ta, tb, 0.5)
        Y.append(y.copy())
    return np.arange(len(Y)) * dt, np.array(Y)


@torch.no_grad()
def eval_grid(model, case, N=160, nt=301):
    xi, tau = torch.linspace(0.0, 1.0, N + 1), torch.linspace(0.0, 1.0, nt)
    XI, TA = torch.meshgrid(xi, tau, indexing="ij")
    P = torch.tensor(case, dtype=torch.float32).repeat(XI.numel(), 1)
    T = model.field(XI.reshape(-1, 1), TA.reshape(-1, 1), P).reshape(N + 1, nt).T
    tj = model.tj(tau.reshape(-1, 1), P[:nt]).reshape(-1)
    return tau.numpy().astype(float) * W_S, tj.numpy().astype(float), T.numpy().astype(float)


def _trap(y, x, axis=-1):
    y = np.moveaxis(np.asarray(y, float), axis, -1)
    return (0.5 * (y[..., 1:] + y[..., :-1]) * np.diff(x)).sum(-1)


def energy_resid(ph, case, t, Tj, T):
    """Residuo relativo do balanco de energia integrado na janela: (dE_armazenada - int(q_in - q_out) dt) / int(q_in + q_out) dt."""
    f, qa, qr = case[0], case[1], case[2]
    xs = np.linspace(0.0, ph.L, T.shape[1])
    P = np.array([mt.ia_p(float(f), float(x) - 273.15)[0] for x in Tj])
    stored = ph.p["n_mod"] * ph.cmod * Tj + ph.ell * (ph.ca * _trap(T, xs, 1) + (ph.cfl / ph.ell) * T[:, 0])
    q_in = ph.p["n_mod"] * P * (1.0 + ph.kap) + qr + ph.ell * qa * ph.L
    q_out = ph.ell * ph.es * _trap(T ** 4 - ph.ts ** 4, xs, 1)
    return float((stored[-1] - stored[0] - _trap(q_in - q_out, t)) / _trap(q_in + q_out, t))


def in_domain(case):
    """True se o caso esta dentro do dominio de TREINO da v2 (faixas de RNG_ e T_j0 - T_raiz0 em DJR)."""
    f, qa, qr, tj0, tr0, dtip = case
    ok = all(RNG_[k][0] <= v <= RNG_[k][1] for k, v in (("f", f), ("q_abs", qa), ("q_rest", qr), ("tr0", tr0), ("dtip", dtip)))
    return bool(ok and DJR[0] <= tj0 - tr0 <= DJR[1])


def compare(model, case, ph, N=160):
    t, Y = solver_window(case, ph.a_pl, N)
    _, tjp, Tp = eval_grid(model, case, N, len(t))
    Tr, Tjr = Y[:, 1:], Y[:, 0]
    dT, dJ = Tp - Tr, tjp - Tjr
    diff = np.concatenate([dT.ravel(), dJ])
    span = float(max(Tr.max(), Tjr.max()) - min(Tr.min(), Tjr.min()))
    return {"in_domain": in_domain(case), "rms_K": float(np.sqrt(np.mean(diff ** 2))), "max_abs_C": float(np.abs(diff).max()),
            "bias_K": float(diff.mean()), "tj_max_abs_C": float(np.abs(dJ).max()), "tj_final_err_C": float(dJ[-1]),
            "tj_peak_err_C": float(tjp.max() - Tjr.max()), "tj_final_ref_C": float(Tjr[-1] - 273.15),
            "max_abs_rel_faixa": float(np.abs(diff).max() / max(span, 1e-9)),
            "energy_resid_pinn": energy_resid(ph, case, t, tjp, Tp), "energy_resid_solver": energy_resid(ph, case, t, Tjr, Tr)}


@torch.no_grad()
def predict_state(model, case, taus):
    """Interface para o Marco 5: (Tj, T_raiz, T_ponta)(tau) para tau em [0,1]; case: tensor [6] fisico (PNAMES)."""
    tau = taus.reshape(-1, 1)
    n = tau.shape[0]
    P = case.reshape(1, -1).repeat(n, 1)
    T = model.field(torch.cat([torch.zeros(n, 1), torch.ones(n, 1)], 0), torch.cat([tau, tau], 0), torch.cat([P, P], 0))
    return model.tj(tau, P).reshape(-1), T[:n].reshape(-1), T[n:].reshape(-1)


def measure_latency(model, n=2000, nt=16, threads=1):
    """Latencia de predict_state (float32, CPU, 1 thread) com pesos quaisquer: nao depende do treino."""
    torch.set_num_threads(threads)
    case = torch.tensor([1000.0, 120.0, 4000.0, 330.0, 310.0, 8.0])
    taus = torch.linspace(1.0 / nt, 1.0, nt)
    for _ in range(100):
        predict_state(model, case, taus)
    ts = []
    for _ in range(n):
        t0 = time.perf_counter()
        predict_state(model, case, taus)
        ts.append((time.perf_counter() - t0) * 1e3)
    ts = np.array(ts)
    return {"mean_ms": float(ts.mean()), "p50_ms": float(np.percentile(ts, 50)), "p99_ms": float(np.percentile(ts, 99)),
            "max_ms": float(ts.max()), "n": n, "nt": nt, "threads": threads}


def cost(model, nt=16):
    def c(net):
        lins = [m for m in net if isinstance(m, nn.Linear)]
        return sum(m.in_features * m.out_features for m in lins), sum(m.out_features for m in lins[:-1])
    mf, tf = c(model.fin)
    mm, tm = c(model.mod)
    return {"macs_predict_state": nt * (mm + 2 * mf), "tanh_predict_state": nt * (tm + 2 * tf), "parametros": sum(p.numel() for p in model.parameters())}


def random_cases(n, seed):
    """MESMOS casos R00.. de v1 (amostragem independente, semente 1234): comparabilidade com E1-E3."""
    return [tuple(float(x) for x in row) for row in sample_P_legacy(n, torch.Generator().manual_seed(seed))]


def _summ(rows):
    v = [r["max_abs_C"] for r in rows]
    if not v:
        return {"n": 0}
    return {"n": len(v), "pior_max_abs_C": float(max(v)), "mediana_max_abs_C": float(np.median(v)), "n_le_0p5C": int(sum(x <= 0.5 for x in v)),
            "rms_K_mediana": float(np.median([r["rms_K"] for r in rows])), "pior_subestimacao_pico_Tj_C": float(min(r["tj_peak_err_C"] for r in rows)),
            "pior_energia_pinn": float(max(abs(r["energy_resid_pinn"]) for r in rows))}


def evaluate(model, ncases=20, seed=1234, N=160):
    ph = Phys()
    cases = {"A_sol_quente": (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0), "B_degrau_frio": (1353.752, 135.2, 5551.99, 285.0, 273.0, 5.0),
             "C_eclipse": (508.07, 99.3, 1510.72, 351.0, 321.0, 6.0)}
    for i, c in enumerate(random_cases(ncases, seed)):
        cases["R%02d" % i] = c
    rows = {k: compare(model, c, ph, N) for k, c in cases.items()}
    parts = losses(model, ph, 64, torch.Generator().manual_seed(999))
    vals = list(rows.values())
    return {"losses_validacao_K2": {n: p.item() for n, p in zip(("L_pde", "L_bc", "L_mod", "L_ic"), parts)},
            "losses_validacao_rms_K": {n: math.sqrt(p.item()) for n, p in zip(("L_pde", "L_bc", "L_mod", "L_ic"), parts)},
            "resumo": {"todos": _summ(vals), "no_dominio_treino_v2": _summ([r for r in vals if r["in_domain"]]),
                       "fora_dominio_treino_v2": _summ([r for r in vals if not r["in_domain"]])},
            "casos": rows, "latencia_cpu": measure_latency(model), "custo": cost(model)}


# ---------------- estimador de Delta_ponta (modelo direto = a propria PINN)
@torch.no_grad()
def estimate_dtip(model, case5, taus, tj_obs, tr_obs, sj=0.5, sr=0.05, grid=None):
    """Minimos quadrados ponderados sobre uma grade de Delta_ponta. case5 = (f, q_abs, q_rest, tj0, tr0); taus em [0,1];
    tj_obs, tr_obs [K]. Retorna {'dtip', 'tj_end_K' (T_j em tau = 1 predita com o Delta estimado), 'cost'}.
    LIMITE: o modelo direto e a PINN, entao o estimador herda o erro dela."""
    grid = np.linspace(0.0, 15.0, 31) if grid is None else np.asarray(grid, float)
    taus = np.asarray(taus, float)
    G, K = len(grid), len(taus)
    tau = torch.tensor(np.tile(taus, G), dtype=torch.float32).reshape(-1, 1)
    P = torch.tensor([[*case5, float(d)] for d in grid for _ in range(K)], dtype=torch.float32)
    tjp = model.tj(tau, P).reshape(G, K).numpy().astype(float)
    trp = model.field(torch.zeros_like(tau), tau, P).reshape(G, K).numpy().astype(float)
    c = ((np.asarray(tj_obs)[None, :] - tjp) ** 2).sum(1) / sj ** 2 + ((np.asarray(tr_obs)[None, :] - trp) ** 2).sum(1) / sr ** 2
    i = int(np.argmin(c))
    P1 = torch.tensor([[*case5, float(grid[i])]], dtype=torch.float32)
    return {"dtip": float(grid[i]), "tj_end_K": float(model.tj(torch.ones(1, 1), P1)), "cost": float(c[i])}


def _tj_end(model, case5, dtip):
    with torch.no_grad():
        return float(model.tj(torch.ones(1, 1), torch.tensor([[*case5, float(dtip)]], dtype=torch.float32)))


def evaluate_estimator(model, ncases=20, seed=77, wobs=60.0, sj=0.5, sr=0.05, bj=1.5, br=0.25, dt_obs=2.0, grid=None):
    """Observa T_j e T_raiz por wobs s (ruido sj, sr; vies uniforme +-bj em T_j ou +-br no RTD), estima Delta_ponta e compara T_j(300 s) com o solver.
    Casos sorteados DENTRO do dominio de treino (amostragem correlacionada). Sensores: valores [HIPOTESE DECLARADA] do usuario."""
    grid = np.linspace(0.0, 15.0, 61) if grid is None else np.asarray(grid, float)
    ph = Phys()
    rng = np.random.default_rng(seed)
    P = sample_P(ncases, torch.Generator().manual_seed(seed)).numpy().astype(float)
    k = int(wobs / dt_obs) + 1
    idx = (np.arange(k) * dt_obs).astype(int)
    taus = idx / W_S
    err = {m: [] for m in ("prior", "pinn_dtip_verdadeiro", "ruido", "vies_Tj", "vies_RTD")}
    derr = {m: [] for m in ("ruido", "vies_Tj", "vies_RTD")}
    for row in P:
        case = tuple(float(x) for x in row)
        case5, dtrue = case[:5], case[5]
        _, Y = solver_window(case, ph.a_pl, 160)
        tj_ref = float(Y[-1, 0])
        z = rng.standard_normal((2, k))
        yj, yr = Y[idx, 0] + sj * z[0], Y[idx, 1] + sr * z[1]
        variants = {"ruido": (yj, yr), "vies_Tj": (yj + rng.uniform(-bj, bj), yr), "vies_RTD": (yj, yr + rng.uniform(-br, br))}
        err["prior"].append(_tj_end(model, case5, 7.5) - tj_ref)
        err["pinn_dtip_verdadeiro"].append(_tj_end(model, case5, dtrue) - tj_ref)
        for nm, (a, b) in variants.items():
            est = estimate_dtip(model, case5, taus, a, b, sj, sr, grid)
            err[nm].append(est["tj_end_K"] - tj_ref)
            derr[nm].append(est["dtip"] - dtrue)
    rms = lambda v: float(np.sqrt(np.mean(np.square(v))))
    case5 = tuple(float(x) for x in P[0][:5])
    taus16 = np.linspace(0.0, wobs / W_S, 16)
    obs = np.full(16, 330.0)
    g31 = np.linspace(0.0, 15.0, 31)
    for _ in range(20):
        estimate_dtip(model, case5, taus16, obs, obs - 20.0, sj, sr, g31)
    ts = []
    for _ in range(200):
        t0 = time.perf_counter()
        estimate_dtip(model, case5, taus16, obs, obs - 20.0, sj, sr, g31)
        ts.append((time.perf_counter() - t0) * 1e3)
    return {"n": int(ncases), "wobs_s": wobs, "rms_erro_Tj300_K": {m: rms(v) for m, v in err.items()},
            "max_abs_erro_Tj300_K": {m: float(np.max(np.abs(v))) for m, v in err.items()},
            "rms_erro_dtip_K": {m: rms(v) for m, v in derr.items()},
            "latencia_estimador_cpu_ms_(grade31_k16)": {"mean": float(np.mean(ts)), "p99": float(np.percentile(ts, 99))}}


# ---------------- testes (pytest)
def test_p_torch_igual_numpy():
    f = torch.linspace(300.0, 1410.0, 25, dtype=torch.float64)
    for tj in (40.0, 80.0, 100.0):
        a = p_mod_t(f, torch.full_like(f, tj)).numpy()
        b = np.array([mt.ia_p(float(x), tj)[0] for x in f])
        assert np.max(np.abs(a - b) / b) < 1e-9


def test_ic_inclinacao_nula_na_ponta():
    xi = torch.ones(5, 1, dtype=torch.float64, requires_grad=True)
    P = torch.tensor([[1000.0, 120.0, 4000.0, 330.0, 310.0, 8.0]], dtype=torch.float64).repeat(5, 1)
    assert float(grad(ic_T(xi, P), xi).abs().max()) < 1e-12


def test_amostragem_de_treino_respeita_faixas_e_correlacao():
    P = sample_P(2000, torch.Generator().manual_seed(0))
    d = P[:, 3] - P[:, 4]
    assert float(d.min()) >= DJR[0] - 1e-4 and float(d.max()) <= DJR[1] + 1e-4
    for c, k in ((0, "f"), (1, "q_abs"), (2, "q_rest"), (4, "tr0"), (5, "dtip")):
        assert float(P[:, c].min()) >= RNG_[k][0] - 1e-4 and float(P[:, c].max()) <= RNG_[k][1] + 1e-4
    assert make_batch(4, torch.Generator().manual_seed(1), sampling="indep")["P"].shape == (4, 6)


def test_residuo_pde_esta_em_kelvin_equivalente():
    ph = Phys()
    q = ph.es * (320.0 ** 4 - ph.ts ** 4) + ph.h * 2.0
    P = torch.tensor([[1000.0, q, 4000.0, 330.0, 310.0, 8.0]], dtype=torch.float64).repeat(7, 1)
    xi, tau = torch.rand(7, 1, dtype=torch.float64), torch.rand(7, 1, dtype=torch.float64)

    def field(x, t, p):
        return torch.zeros_like(x) + 320.0 + 0.0 * t + 0.0 * x ** 2
    r = res_pde(field, ph, xi, tau, P)
    assert float((r + 2.0).abs().max()) < 1e-6          # radiacao em excesso de 2 K-equivalente -> residuo = -2


def test_residuo_continuo_em_perfil_estacionario_do_solver():
    # tolerancias escolhidas por mim (nao calibradas por execucao): diagnostico, nao prova.
    ph, case = Phys(), (1353.752, 135.2, 5551.99)
    y = mt.Plant(A_PL, {"n_nodes": 160}, mt.static_drive(*case)).steady()
    assert y is not None
    c = np.polynomial.polynomial.polyfit(np.linspace(0.0, 1.0, 161), y[1:] - T_REF, 12)
    tj_s = float(y[0])

    def field(xi, tau, P):
        out = torch.zeros_like(xi) + T_REF
        for k, ck in enumerate(c):
            out = out + float(ck) * xi ** k
        return out + 0.0 * tau

    def tjf(tau, P):
        return torch.zeros_like(tau) + tj_s + 0.0 * tau
    g, n = torch.Generator().manual_seed(5), 200
    P = torch.tensor([case[0], case[1], case[2], 0.0, 0.0, 0.0], dtype=torch.float64).repeat(n, 1)
    xi = (0.02 + 0.96 * torch.rand(n, 1, generator=g)).double()
    tau = torch.rand(n, 1, generator=g).double()
    r_root, r_mod = res_root_mod(field, tjf, ph, tau, P)
    # tolerancias = as de v1 (que passaram localmente) convertidas para K-equivalente: PDE 1e-2*1000/h = 0.75 K; ponta 1e-2*DTK = 0.5 K;
    # raiz 5e-2*(k delta DTK/L)/(s/R) = 5.9 K; modulo 1e-2*400 W*R = 0.32 K. Diagnostico, nao prova.
    assert float(res_pde(field, ph, xi, tau, P).abs().max()) < 0.75
    assert float(res_tip(field, ph, tau, P).abs().max()) < 0.5
    assert float(r_root.abs().max()) < 5.9 and float(r_mod.abs().max()) < 0.32


def test_solver_window_conserva_energia():
    ph, case = Phys(), (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0)
    t, Y = solver_window(case, N=80)
    assert abs(energy_resid(ph, case, t, Y[:, 0], Y[:, 1:])) < 1e-5


def test_forward_formas_e_finitude():
    m, P = PINN(), sample_P(7, torch.Generator().manual_seed(1))
    xi, ta = torch.rand(7, 1), torch.rand(7, 1)
    assert m.field(xi, ta, P).shape == (7, 1) and m.tj(ta, P).shape == (7, 1)
    assert bool(torch.isfinite(m.field(xi, ta, P)).all()) and bool(torch.isfinite(m.tj(ta, P)).all())
    z = torch.zeros(7, 1)
    assert float((m.field(xi, z, P) - ic_T(xi, P)).abs().max()) < 1e-5 and float((m.tj(z, P) - P[:, 3:4]).abs().max()) < 1e-5   # hard_ic


def test_perdas_finitas_e_gradientes():
    m, ph = PINN(width=16, depth=2), Phys()
    parts = losses(m, ph, 8, torch.Generator().manual_seed(2), nint=8, nbc=4, nic=4)
    assert all(bool(torch.isfinite(x)) for x in parts)
    total(parts).backward()
    assert all(bool(torch.isfinite(p.grad).all()) for p in m.parameters() if p.grad is not None)


def test_treino_curto_reduz_perda():
    _, hist = train(epochs=150, ns=16, seed=3, out=None, log_every=75)
    assert hist[-1]["loss"] < hist[0]["loss"]


def test_reprodutibilidade_semente():
    a, _ = train(epochs=3, ns=8, seed=7, out=None, log_every=100)
    b, _ = train(epochs=3, ns=8, seed=7, out=None, log_every=100)
    assert all(torch.allclose(x, y, atol=1e-7) for x, y in zip(a.parameters(), b.parameters()))


def test_lbfgs_nao_aumenta_a_perda():
    torch.manual_seed(0)
    m, ph = PINN(width=16, depth=2), Phys()
    fixed = make_batch(8, torch.Generator().manual_seed(4), nint=8, nbc=4, nic=4)
    before = total(losses_from_batch(m, ph, fixed)).item()
    h = lbfgs_phase(m, ph, fixed, 2, log_every=1)
    assert h[-1]["loss"] <= before + 1e-9, (before, h[-1]["loss"])


def test_pipeline_de_avaliacao_roda():
    r = compare(PINN(width=16, depth=2), (1000.0, 120.0, 4000.0, 330.0, 310.0, 8.0), Phys(), N=40)
    assert all(math.isfinite(float(v)) for v in r.values())


def test_estimador_recupera_dtip_com_modelo_direto_conhecido():
    class Fake:
        def tj(self, tau, P):
            return P[:, 3:4] + tau * (2.0 - 0.1 * P[:, 5:6])

        def field(self, xi, tau, P):
            return P[:, 4:5] - tau * 0.4 * P[:, 5:6]
    taus = np.linspace(0.0, 0.2, 16)
    tj_obs = 330.0 + taus * (2.0 - 0.1 * 6.0)
    tr_obs = 310.0 - taus * 0.4 * 6.0
    est = estimate_dtip(Fake(), (1000.0, 120.0, 4000.0, 330.0, 310.0), taus, tj_obs, tr_obs, grid=np.linspace(0.0, 15.0, 31))
    assert abs(est["dtip"] - 6.0) < 1e-9 and abs(est["tj_end_K"] - (330.0 + (2.0 - 0.6))) < 1e-4, est


def test_latencia_cpu_menor_que_5ms():
    lat = measure_latency(PINN(), n=300)
    assert lat["mean_ms"] < 5.0 and lat["p99_ms"] < 5.0, lat


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train")
    t.add_argument("--epochs", type=int, default=20000)
    t.add_argument("--ns", type=int, default=48)
    t.add_argument("--lr", type=float, default=1e-3)
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--width", type=int, default=64)
    t.add_argument("--depth", type=int, default=4)
    t.add_argument("--no_hard_ic", action="store_true")
    t.add_argument("--sampling", choices=("corr", "indep"), default="corr")
    t.add_argument("--lbfgs_steps", type=int, default=0)
    t.add_argument("--ns_lbfgs", type=int, default=96)
    t.add_argument("--init", default=None)
    for nm in ("w_pde", "w_bc", "w_mod", "w_ic"):
        t.add_argument("--" + nm, type=float, default=1.0)
    t.add_argument("--out", default="pinn.pt")
    e = sub.add_parser("eval")
    e.add_argument("--ckpt", default="pinn.pt")
    e.add_argument("--ncases", type=int, default=20)
    e.add_argument("--reduced", action="store_true")
    s = sub.add_parser("est")
    s.add_argument("--ckpt", default="pinn.pt")
    s.add_argument("--ncases", type=int, default=20)
    s.add_argument("--wobs", type=float, default=60.0)
    a = ap.parse_args()
    if a.cmd == "train":
        train(a.epochs, a.ns, a.lr, a.seed, (a.w_pde, a.w_bc, a.w_mod, a.w_ic), a.width, a.depth, not a.no_hard_ic, a.out,
              lbfgs_steps=a.lbfgs_steps, ns_lbfgs=a.ns_lbfgs, init=a.init, sampling=a.sampling)
    elif a.cmd == "eval":
        r = evaluate(load(a.ckpt), a.ncases)
        if a.reduced:
            r = {k: r[k] for k in ("losses_validacao_K2", "losses_validacao_rms_K", "resumo", "latencia_cpu", "custo")} | \
                {"casos": {k: v for k, v in r["casos"].items() if k in ("A_sol_quente", "B_degrau_frio", "R00", "R04")}}
        print(json.dumps(r, indent=1, ensure_ascii=False))
    else:
        print(json.dumps(evaluate_estimator(load(a.ckpt), a.ncases, wobs=a.wobs), indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()

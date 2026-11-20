"""modelo_pinn.py -- Marco 4: Modelo Substituto Preditivo Baseado em Fisica (PINN) da planta termica do Marco 3.
Verdade-base: modelo_termico.py (mesma pasta). Janela de predicao W = 300 s. NAO e Gemeo Digital (sem enlace com hardware/telemetria).
Etiquetas: [HIPOTESE DECLARADA] = escolha minha sem fonte; [CALCULADO]; [USUARIO] = confirmado localmente.
Nao validado experimentalmente. Codigo NAO executado pelo autor (sem torch no sandbox).
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

# ---- escalas e dominio de treino [HIPOTESE DECLARADA]
T_REF, DTK, W_S = 320.0, 50.0, 300.0
A_PL = 49.89754                      # m2, area de projeto do Marco 3 [USUARIO]
PNAMES = ("f", "q_abs", "q_rest", "tj0", "tr0", "dtip")
_RANGE = {"f": (300.0, 1410.0), "q_abs": (90.0, 145.0), "q_rest": (800.0, 7500.0),
          "tj0": (280.0, 360.0), "tr0": (265.0, 330.0), "dtip": (0.0, 15.0)}
LO = torch.tensor([_RANGE[k][0] for k in PNAMES])
HI = torch.tensor([_RANGE[k][1] for k in PNAMES])


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
    """Perfil inicial do campo: quadratico, inclinacao nula na ponta."""
    tr = P[:, 4:5]
    tt = tr - P[:, 5:6]
    return tt + (tr - tt) * (1.0 - xi) ** 2


class PINN(nn.Module):
    def __init__(self, width=64, depth=4, hard_ic=False):
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
        return ic_T(xi, P) + tau * DTK * raw if self.hard_ic else T_REF + DTK * raw

    def tj(self, tau, P):
        raw = self.mod(torch.cat([2.0 * tau - 1.0, norm_p(P)], 1))
        return P[:, 3:4] + tau * DTK * raw if self.hard_ic else T_REF + DTK * raw


def grad(y, x):
    return torch.autograd.grad(y, x, grad_outputs=torch.ones_like(y), create_graph=True)[0]


def res_pde(field, ph, xi, tau, P):
    xi, tau = xi.clone().requires_grad_(True), tau.clone().requires_grad_(True)
    T = field(xi, tau, P)
    Tt = grad(T, tau) / W_S
    Tx = grad(T, xi) / ph.L
    Txx = grad(Tx, xi) / ph.L
    return (ph.ca * Tt - ph.kd * Txx + ph.es * (T ** 4 - ph.ts ** 4) - P[:, 1:2]) / 1000.0


def res_tip(field, ph, tau, P):
    xi = torch.ones_like(tau).requires_grad_(True)
    return grad(field(xi, tau, P), xi) / DTK


def res_root_mod(field, tjfun, ph, tau, P):
    tau = tau.clone().requires_grad_(True)
    xi = torch.zeros_like(tau).requires_grad_(True)
    T0 = field(xi, tau, P)
    Tx0, Tt0 = grad(T0, xi) / ph.L, grad(T0, tau) / W_S
    Tj = tjfun(tau, P)
    Tjt = grad(Tj, tau) / W_S
    Pm = p_mod_t(P[:, 0:1], Tj - 273.15)
    r_root = (ph.cfl / ph.ell * Tt0 - ph.kd * Tx0 - ph.s * (Tj - T0) / ph.R - ph.s * ph.kap * Pm - P[:, 2:3] / ph.ell) / (ph.kd * DTK / ph.L)
    r_mod = (ph.cmod * Tjt - (Pm - (Tj - T0) / ph.R)) / 400.0
    return r_root, r_mod


def res_ic(field, tjfun, xi, P):
    tau0 = torch.zeros_like(xi)
    return (field(xi, tau0, P) - ic_T(xi, P)) / DTK, (tjfun(tau0, P) - P[:, 3:4]) / DTK


def sample_P(n, gen):
    return LO + torch.rand(n, len(PNAMES), generator=gen) * (HI - LO)


def rep(P, k):
    return P.repeat_interleave(k, 0)


def losses(model, ph, ns, gen, nint=96, nbc=32, nic=32):
    """Retorna (L_pde, L_bc, L_mod, L_ic), cada uma media de quadrados de residuos normalizados."""
    P = sample_P(ns, gen)
    xi, ta = torch.rand(ns * nint, 1, generator=gen), torch.rand(ns * nint, 1, generator=gen)
    L_pde = res_pde(model.field, ph, xi, ta, rep(P, nint)).pow(2).mean()
    Pb, tb = rep(P, nbc), torch.rand(ns * nbc, 1, generator=gen)
    r_root, r_mod = res_root_mod(model.field, model.tj, ph, tb, Pb)
    L_bc = res_tip(model.field, ph, tb, Pb).pow(2).mean() + r_root.pow(2).mean()
    xc = torch.rand(ns * nic, 1, generator=gen)
    r_f, r_j = res_ic(model.field, model.tj, xc, rep(P, nic))
    return L_pde, L_bc, r_mod.pow(2).mean(), r_f.pow(2).mean() + r_j.pow(2).mean()


def train(epochs=6000, ns=48, lr=1e-3, seed=0, w=(1.0, 1.0, 1.0, 10.0), width=64, depth=4, hard_ic=False,
          out="pinn.pt", log_every=250):
    """w = (w_pde, w_bc, w_mod, w_ic) [HIPOTESE DECLARADA]. Adam + cosseno; sementes fixas."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    gen = torch.Generator().manual_seed(seed)
    ph, model = Phys(), PINN(width, depth, hard_ic)
    opt = torch.optim.Adam(model.parameters(), lr=lr)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=max(epochs, 1), eta_min=lr * 0.01)
    hist = []
    for it in range(1, epochs + 1):
        opt.zero_grad()
        Lp, Lb, Lm, Li = losses(model, ph, ns, gen)
        loss = w[0] * Lp + w[1] * Lb + w[2] * Lm + w[3] * Li
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        opt.step()
        sched.step()
        if it % log_every == 0 or it == 1:
            hist.append({"it": it, "L_pde": float(Lp), "L_bc": float(Lb), "L_mod": float(Lm), "L_ic": float(Li), "loss": float(loss)})
            print(json.dumps(hist[-1]), flush=True)
    if out:
        torch.save({"state": model.state_dict(), "width": width, "depth": depth, "hard_ic": hard_ic, "seed": seed, "hist": hist}, out)
    return model, hist


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


def compare(model, case, ph, N=160):
    t, Y = solver_window(case, ph.a_pl, N)
    _, tjp, Tp = eval_grid(model, case, N, len(t))
    Tr, Tjr = Y[:, 1:], Y[:, 0]
    dT, dJ = Tp - Tr, tjp - Tjr
    diff = np.concatenate([dT.ravel(), dJ])
    ref = np.concatenate([Tr.ravel(), Tjr])
    ref0 = np.concatenate([np.broadcast_to(Tr[0], Tr.shape).ravel(), np.full(Tjr.shape, Tjr[0])])
    return {"rel_L2_C": float(np.linalg.norm(diff) / np.linalg.norm(ref - 273.15)),
            "rel_L2_var": float(np.linalg.norm(diff) / max(np.linalg.norm(ref - ref0), 1e-12)),
            "max_abs_C": float(np.abs(diff).max()), "tj_max_abs_C": float(np.abs(dJ).max()),
            "tj_final_err_C": float(dJ[-1]), "tj_final_ref_C": float(Tjr[-1] - 273.15),
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
    return [tuple(float(x) for x in row) for row in sample_P(n, torch.Generator().manual_seed(seed))]


def evaluate(model, ncases=20, seed=1234, N=160):
    ph = Phys()
    cases = {"A_sol_quente": (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0), "B_degrau_frio": (1353.752, 135.2, 5551.99, 285.0, 273.0, 5.0),
             "C_eclipse": (508.07, 99.3, 1510.72, 351.0, 321.0, 6.0)}
    for i, c in enumerate(random_cases(ncases, seed)):
        cases["R%02d" % i] = c
    rows = {k: compare(model, c, ph, N) for k, c in cases.items()}
    Lp, Lb, Lm, Li = losses(model, ph, 64, torch.Generator().manual_seed(999))
    worst = {m: max(abs(r[m]) for r in rows.values()) for m in ("rel_L2_C", "rel_L2_var", "max_abs_C", "tj_max_abs_C", "energy_resid_pinn", "energy_resid_solver")}
    return {"losses_validacao": {"L_pde": float(Lp), "L_bc": float(Lb), "L_mod": float(Lm), "L_ic": float(Li)},
            "pior_caso": worst, "casos": rows, "latencia_cpu": measure_latency(model), "custo": cost(model)}


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


def test_residuo_continuo_em_perfil_estacionario_do_solver():
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
    assert float(res_pde(field, ph, xi, tau, P).abs().max()) < 1e-2
    assert float(res_tip(field, ph, tau, P).abs().max()) < 1e-2
    assert float(r_root.abs().max()) < 5e-2 and float(r_mod.abs().max()) < 1e-2


def test_solver_window_conserva_energia():
    ph, case = Phys(), (1353.752, 135.2, 5551.99, 351.0, 321.0, 6.0)
    t, Y = solver_window(case, N=80)
    assert abs(energy_resid(ph, case, t, Y[:, 0], Y[:, 1:])) < 1e-5


def test_forward_formas_e_finitude():
    m, P = PINN(), sample_P(7, torch.Generator().manual_seed(1))
    xi, ta = torch.rand(7, 1), torch.rand(7, 1)
    assert m.field(xi, ta, P).shape == (7, 1) and m.tj(ta, P).shape == (7, 1)
    assert bool(torch.isfinite(m.field(xi, ta, P)).all()) and bool(torch.isfinite(m.tj(ta, P)).all())


def test_perdas_finitas_e_gradientes():
    m, ph = PINN(width=16, depth=2), Phys()
    parts = losses(m, ph, 8, torch.Generator().manual_seed(2), nint=8, nbc=4, nic=4)
    assert all(bool(torch.isfinite(x)) for x in parts)
    sum(parts).backward()
    assert all(bool(torch.isfinite(p.grad).all()) for p in m.parameters() if p.grad is not None)


def test_treino_curto_reduz_perda():
    _, hist = train(epochs=150, ns=16, seed=3, out=None, log_every=75)
    assert hist[-1]["loss"] < hist[0]["loss"]


def test_reprodutibilidade_semente():
    a, _ = train(epochs=3, ns=8, seed=7, out=None, log_every=100)
    b, _ = train(epochs=3, ns=8, seed=7, out=None, log_every=100)
    assert all(torch.allclose(x, y, atol=1e-7) for x, y in zip(a.parameters(), b.parameters()))


def test_pipeline_de_avaliacao_roda():
    r = compare(PINN(width=16, depth=2), (1000.0, 120.0, 4000.0, 330.0, 310.0, 8.0), Phys(), N=40)
    assert all(math.isfinite(v) for v in r.values())


def test_latencia_cpu_menor_que_5ms():
    lat = measure_latency(PINN(), n=300)
    assert lat["mean_ms"] < 5.0 and lat["p99_ms"] < 5.0, lat


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("train")
    t.add_argument("--epochs", type=int, default=6000)
    t.add_argument("--ns", type=int, default=48)
    t.add_argument("--lr", type=float, default=1e-3)
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--width", type=int, default=64)
    t.add_argument("--depth", type=int, default=4)
    t.add_argument("--hard_ic", action="store_true")
    for nm, df in (("w_pde", 1.0), ("w_bc", 1.0), ("w_mod", 1.0), ("w_ic", 10.0)):
        t.add_argument("--" + nm, type=float, default=df)
    t.add_argument("--out", default="pinn.pt")
    e = sub.add_parser("eval")
    e.add_argument("--ckpt", default="pinn.pt")
    e.add_argument("--ncases", type=int, default=20)
    a = ap.parse_args()
    if a.cmd == "train":
        train(a.epochs, a.ns, a.lr, a.seed, (a.w_pde, a.w_bc, a.w_mod, a.w_ic), a.width, a.depth, a.hard_ic, a.out)
    else:
        print(json.dumps(evaluate(load(a.ckpt), a.ncases), indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()

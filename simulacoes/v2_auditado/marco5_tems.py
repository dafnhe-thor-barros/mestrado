"""marco5_tems.py -- Marco 5: malha fechada TEMS (DVFS-only) sobre a planta do Marco 3.
Controladores: S0 (40 kW continuo), S1 (estatico energeticamente viavel), R (PI reativo), TEMS com preditor RC3 (Euler / Radau), SOLVER (oraculo) ou PINN.
Planta "verdadeira" = Plant (modelo_termico.py, N = 40). O controlador so ve sensores (T_j, T_raiz) com ruido/vies do contrato.
Requer: numpy, scipy (so RC3 Radau), modelo_termico.py e contrato_parametros.py na mesma pasta. torch + modelo_pinn.py SO para o preditor PINN.
Uso:  python marco5_tems.py --scenarios nominal burst --ctrls S0 S1 R TEMS_RC3 --seeds 0 1 --sensor sem_calibracao --norb 3 --out saida.json
"""
import argparse
import json
import math
import time

import numpy as np

from modelo_termico import (T0K, SIGMA, T_ORB, T_ECL, IA, Plant, static_drive, _merge, optics,
                            ia_p, ia_f_for_power, size_bess, analytic_orbit, orbit_breaks)
from contrato_parametros import (A_PL_PROJETO, PROJETO_1P0, LEGADO_2B_0P90, M_PRED_K, M_OOD_MAX_K, OOD_CORTES, DTIP_PRIOR, DTIP_RANGE,
                                 SENSOR_MARGEM_K, sensor_cfg, contract_hash, sha256_file, verify_checkpoint)


def _now():
    return time.perf_counter()


def _r(o, nd=7):
    if isinstance(o, (float, np.floating)):
        o = float(o)
        return float("%.*g" % (nd, o)) if math.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, dict):
        return {k: _r(v, nd) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_r(v, nd) for v in o]
    if isinstance(o, np.ndarray):
        return _r(o.tolist(), nd)
    return o


def ceil_div(a, b):
    return max(1, int(math.ceil(a / b - 1e-9)))


def bisect_max(ok, lo, hi, n=7):
    """Maior f em [lo, hi] com ok(f) verdadeiro, assumindo ok monotono decrescente em f e ok(lo) verdadeiro."""
    if ok(hi):
        return hi
    for _ in range(n):
        mid = 0.5 * (lo + hi)
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return lo


def plant_par(c, extra=None):
    d = {"n_mod": c["n_mod"], "kappa_dab": 1.0 / (c["eta_dab"] * c["eta_pol"]) - 1.0}
    if extra:
        d.update(extra)
    return d


def _bess_kw(c):
    kw = {k: c[k] for k in ("eta_buck", "eta_dab", "eta_pol", "eta_bd", "p_aux", "eta_ch", "eta_dis", "fade", "n_ser", "e_cell", "soc_min", "soc_max")}
    kw["p_pv"] = c["p_pv"]
    return kw


def battery_for(c):
    s = size_bess(T_ECL, T_ORB, 18000.0, **_bess_kw(c))
    return s, s["n_cells"] * c["e_cell"] * (1.0 - c["fade"])


def sun_cap(c):
    return analytic_orbit(T_ECL, T_ORB, 18000.0, **_bess_kw(c))["P_TI_sun_sust_W"]


def env_eclipse(t, p):
    T, te = p["t_orb"], p["t_ecl"]
    return abs((t % T) - 0.5 * T) < 0.5 * te


def env_q_abs(t, p):
    """Absorvido por area de planta (2 faces): albedo (alpha G a F mu) + IR terrestre; opcao 1 (phi = 90 graus, sem Sol direto)."""
    alpha, eps = optics(p)
    T = p["t_orb"]
    mu = max(0.0, math.cos(math.radians(p["beta_deg"])) * math.cos(2.0 * math.pi * (t % T) / T))
    return 2.0 * (alpha * p["g_sun"] * p["albedo"] * p["f_view"] * mu + eps * p["q_ir"] * p["f_view"])


def orbit_break_times(p, t0, t1):
    nb = int(t1 / p["t_orb"]) + 2
    return [b for b in orbit_breaks(p, nb) if t0 + 1e-6 < b < t1 - 1e-6]


def bus_step(c, pv_out, p_ti, dt, E, ecap, floor):
    """Balanco do barramento num passo dt [s]. E em Wh quimicos. floor = fracao minima de SoC que a protecao respeita.
    Retorna a TI efetivamente servida (reduzida se o BESS nao cobre), o calor externo ao no (q_rest) e o residuo de balanco."""
    chain = c["eta_dab"] * c["eta_pol"]
    l_req = p_ti / chain + c["p_aux"]
    pv_bus = c["eta_buck"] * pv_out
    net = pv_bus - l_req
    e_max, e_min = c["soc_max"] * ecap, floor * ecap
    ti, chem_in, chem_out, q_bd, q_bat = p_ti, 0.0, 0.0, 0.0, 0.0
    if net >= 0.0:
        head = max(0.0, e_max - E) * 3600.0 / dt
        to_bat = min(net, head / (c["eta_ch"] * c["eta_bd"]))
        s = c["eta_bd"] * to_bat
        chem_in = c["eta_ch"] * s
        q_bd, q_bat = to_bat - s, s - chem_in
        used_bus = l_req + to_bat
    else:
        avail = max(0.0, E - e_min) * 3600.0 / dt
        deficit = -net
        dmax = avail * c["eta_dis"] * c["eta_bd"]
        if deficit > dmax:
            deficit = dmax
            ti = max(0.0, (pv_bus + deficit - c["p_aux"]) * chain)
        s = deficit / c["eta_bd"]
        chem_out = s / c["eta_dis"]
        q_bd, q_bat = s - deficit, chem_out - s
        used_bus = pv_bus
    used_out = used_bus / c["eta_buck"]
    q_buck = used_out - used_bus
    q_rest = c["q_aux"] + q_buck + q_bd + (q_bat if c["bat_losses_in_qrest"] else 0.0)
    resid = used_out - (q_buck + ti / chain + c["p_aux"] + q_bd + q_bat + chem_in - chem_out)
    return {"ti": ti, "E": E + (chem_in - chem_out) * dt / 3600.0, "q_rest": q_rest, "used_out": used_out,
            "curtail": pv_out - used_out, "resid": resid}


def q_rest_model(c, ecap, pv_on, p_ti):
    """Modelo do controlador para o calor externo ao no: BESS suposto capaz de cobrir o balanco."""
    return bus_step(c, c["p_pv"] if pv_on else 0.0, p_ti, 60.0, 0.825 * ecap, ecap, 0.0)["q_rest"]


class Sensors:
    def __init__(self, cfg, seed):
        self.c = dict(cfg)
        self.rng = np.random.default_rng(10007 + int(seed))

    def read(self, tj_k, tr_k):
        c = self.c
        z = self.rng.standard_normal(2)
        tj = tj_k + c["bias_sign_tj"] * c["bias_tj"] + c["sigma_tj"] * z[0]
        tr = tr_k + c["bias_sign_rtd"] * c["bias_rtd"] + c["sigma_rtd"] * z[1]
        qm = c.get("quant_mode")
        if qm:
            x = (tj - T0K) / c["quant_step"]
            tj = (np.round(x) if qm == "round" else np.floor(x)) * c["quant_step"] + T0K
        return float(tj), float(tr)


def quad_profile(N, tr, dtip):
    xi = np.linspace(0.0, 1.0, N + 1)
    tt = tr - dtip
    return tt + (tr - tt) * (1.0 - xi) ** 2


# ---------------- preditores
class PredictorBase:
    uses_truth = False
    name = "BASE"

    def domain_distance(self, case):
        return 0.0

    def estimate(self, obs, f, qa, qr, tj0, tr0, sj, sr, grid=None):
        """Minimos quadrados ponderados de Delta_ponta no inicio da janela; o modelo direto e o proprio preditor."""
        grid = np.linspace(0.0, 15.0, 31) if grid is None else np.asarray(grid, float)
        dur = float(obs["t"][-1])
        n = len(obs["t"]) - 1
        best = None
        for d in grid:
            w = self.window(self.init_state(tj0, tr0, float(d), None), f, qa, qr, dur, dur / n)
            m = min(len(w["tj"]), n + 1)
            cost = np.sum((obs["tj"][:m] - w["tj"][:m]) ** 2) / sj ** 2 + np.sum((obs["tr"][:m] - w["tr"][:m]) ** 2) / sr ** 2
            if best is None or cost < best[0]:
                best = (cost, float(d))
        return best[1]


class RC3(PredictorBase):
    """RC de 3 nos (juncao, coletor/raiz, aleta com temperatura MEDIA de perfil parabolico). Parametros = constantes diretas do Marco 3, sem ajuste.
    G_fin = ell k delta / (L/3). Tm = T_raiz - (2/3) Delta_ponta."""

    def __init__(self, par, a_pl, method="euler"):
        p = _merge(par)
        alpha, eps = optics(p)
        ell = a_pl / p["L"]
        self.method, self.name = method, "RC3_" + method
        self.n, self.kap, self.ts, self.a_pl = p["n_mod"], p["kappa_dab"], p["t_space"], a_pl
        self.cj, self.gj = p["n_mod"] * p["c_mod"], p["n_mod"] / p["r_int"]
        self.c0, self.cm = p["c_fluid"], p["c_area"] * a_pl
        self.gf, self.esa = 3.0 * ell * p["k"] * p["delta"] / p["L"], 2.0 * eps * SIGMA * a_pl

    def init_state(self, tj, tr, dtip, y):
        return np.array([tj, tr, tr - (2.0 / 3.0) * dtip], float)

    def dtip_of(self, s):
        return float(1.5 * (s[1] - s[2]))

    def _rhs(self, x, f, qa, qr):
        tj, t0, tm = x
        pm, pl = ia_p(f, tj - T0K)
        P = self.n * pm
        g = self.gj * (tj - t0)
        q = self.gf * (t0 - tm)
        v = np.array([(P - g) / self.cj, (g + self.kap * P + qr - q) / self.c0,
                      (q - self.esa * (tm ** 4 - self.ts ** 4) + qa * self.a_pl) / self.cm])
        return v, self.n * IA["beta"] * pl

    def _step(self, x0, h, f, qa, qr):
        x = x0.copy()
        for _ in range(30):
            fv, dP = self._rhs(x, f, qa, qr)
            J = np.eye(3) - h * np.array([
                [(dP - self.gj) / self.cj, self.gj / self.cj, 0.0],
                [(self.gj + self.kap * dP) / self.c0, (-self.gj - self.gf) / self.c0, self.gf / self.c0],
                [0.0, self.gf / self.cm, (-self.gf - 4.0 * self.esa * x[2] ** 3) / self.cm]])
            dx = np.linalg.solve(J, -(x - x0 - h * fv))
            x = x + dx
            if float(np.max(np.abs(dx))) < 1e-9:
                return x
        raise RuntimeError("RC3: Newton nao convergiu")

    def window(self, state, f, qa, qr, dur, dt_max):
        x0 = np.array(state, float)
        n = ceil_div(dur, dt_max)
        if self.method == "radau":
            from scipy.integrate import solve_ivp
            sol = solve_ivp(lambda t, x: self._rhs(x, f, qa, qr)[0], (0.0, dur), x0, method="Radau",
                            t_eval=np.linspace(0.0, dur, n + 1), rtol=1e-8, atol=1e-6)
            Y = sol.y.T
        else:
            h, x, xs = dur / n, x0, [x0]
            for _ in range(n):
                x = self._step(x, h, f, qa, qr)
                xs.append(x)
            Y = np.array(xs)
        return {"tj": Y[:, 0], "tr": Y[:, 1], "end": Y[-1].copy()}


class SolverPred(PredictorBase):
    """Oraculo: o proprio solver do Marco 3 (N = 40) com passo dt_max; recebe o estado COMPLETO da planta."""
    uses_truth = True
    name = "SOLVER"

    def __init__(self, par, a_pl, n_nodes=40):
        p = dict(par)
        p["n_nodes"] = n_nodes
        self.pl = Plant(a_pl, p, static_drive(1000.0, 0.0, 0.0))
        self.N = self.pl.N

    def init_state(self, tj, tr, dtip, y):
        if y is not None:
            return np.array(y, float)
        return np.concatenate(([tj], quad_profile(self.N, tr, dtip)))

    def dtip_of(self, y):
        return float(y[1] - y[-1])

    def window(self, state, f, qa, qr, dur, dt_max):
        self.pl.drive = static_drive(f, qa, qr)
        n = ceil_div(dur, dt_max)
        h = dur / n
        y = np.array(state, float)
        tj, tr = [y[0]], [y[1]]
        for k in range(n):
            ta, tb = k * h, (k + 1) * h
            if k == 0:
                y, _ = self.pl.step(y, ta, 0.5 * (ta + tb), 1.0)
                y, _ = self.pl.step(y, 0.5 * (ta + tb), tb, 1.0)
            else:
                y, _ = self.pl.step(y, ta, tb, 0.5)
            tj.append(y[0])
            tr.append(y[1])
        return {"tj": np.array(tj), "tr": np.array(tr), "end": y}


class PINNPred(PredictorBase):
    """PINN v2 (modelo_pinn.py) como preditor.
    ATENCAO [NAO VERIFICADO]: este codigo NUNCA foi executado (o sandbox do assistente nao tem torch) e os NOMES da API de modelo_pinn usados abaixo
    (load, W_S, predict_state, estimate_dtip, RNG_, DJR) foram assumidos, nao lidos. test_marco5.py::test_pinn_api_presente aponta qualquer divergencia.
    Estado = (T_j, T_raiz, Delta_ponta)."""
    name = "PINN"

    def __init__(self, ckpt, verify=True):
        import torch
        import modelo_pinn
        self.torch, self.mp = torch, modelo_pinn
        if verify:
            verify_checkpoint(ckpt)
        self.model = self.mp.load(ckpt)
        self.W = self.mp.W_S

    def init_state(self, tj, tr, dtip, y):
        return np.array([tj, tr, min(max(dtip, 0.0), 15.0)], float)

    def dtip_of(self, s):
        return float(s[2])

    def window(self, state, f, qa, qr, dur, dt_max):
        if dur > self.W + 1e-9:
            raise ValueError("janela maior que W")
        n = min(max(ceil_div(dur, dt_max), 4), 32)
        taus = self.torch.linspace(0.0, dur / self.W, n + 1)
        case = self.torch.tensor([f, qa, qr, state[0], state[1], state[2]], dtype=self.torch.float32)
        tj, tr, tt = self.mp.predict_state(self.model, case, taus)
        tj, tr, tt = tj.numpy().astype(float), tr.numpy().astype(float), tt.numpy().astype(float)
        end = np.array([tj[-1], tr[-1], min(max(tr[-1] - tt[-1], 0.0), 15.0)])
        return {"tj": tj, "tr": tr, "end": end}

    def estimate(self, obs, f, qa, qr, tj0, tr0, sj, sr, grid=None):
        taus = np.asarray(obs["t"], float) / self.W
        return float(self.mp.estimate_dtip(self.model, (f, qa, qr, tj0, tr0), taus, obs["tj"], obs["tr"], sj, sr,
                                           np.linspace(0.0, 15.0, 61) if grid is None else grid)["dtip"])

    def domain_distance(self, case):
        """Distancia L-infinito normalizada a caixa de treino da v2. [HIPOTESE DECLARADA] reconstrucao minha: nao li a definicao do arquivo de auditoria."""
        f, qa, qr, tj0, tr0, dtip = case
        mp = self.mp
        ex = [max(0.0, lo - v, v - hi) / (hi - lo) for v, (lo, hi) in
              ((f, mp.RNG_["f"]), (qa, mp.RNG_["q_abs"]), (qr, mp.RNG_["q_rest"]), (tr0, mp.RNG_["tr0"]), (dtip, mp.RNG_["dtip"]), (tj0 - tr0, mp.DJR))]
        return float(max(ex))


# ---------------- controlador
class Controller:
    def __init__(self, cfg, c, p_nom, a_pl, ecap):
        self.cfg, self.c, self.p, self.a_pl, self.ecap = cfg, c, p_nom, a_pl, ecap
        self.kind = cfg.get("ctrl", "TEMS")
        self.sup = cfg.get("supervisor", self.kind in ("R", "TEMS"))
        self.floor = c["soc_min"] if self.sup else 0.0
        self.dtc, self.W = c["dt_ctrl"], c["w_s"]
        self.pred, self.fb = cfg.get("pred"), cfg.get("fallback")
        self.m_pred, self.m_sensor, self.m_total = cfg.get("m_pred", M_PRED_K), cfg.get("m_sensor", 0.0), cfg.get("m_total")
        self.boost, self.dt_pred = cfg.get("boost", False), cfg.get("dt_pred", 15.0)
        scfg = cfg.get("sensors") or sensor_cfg("ideal")
        self.sj, self.sr = max(scfg["sigma_tj"], 0.05), max(scfg["sigma_rtd"], 0.01)
        self.sun = sun_cap(c)
        n = c["n_mod"]
        self.f_s1 = (ia_f_for_power(self.sun / n, c["tj_lim_c"])[0], ia_f_for_power(18000.0 / n, c["tj_lim_c"])[0])
        self.kp, self.ki, self.r_margin = cfg.get("r_kp", 20.0), cfg.get("r_ki", 0.5), cfg.get("r_margin", 1.0)
        self.f_r, self.e_prev = c["f_max"], None
        self.win, self.f_prev_applied, self.dtip = [], None, DTIP_PRIOR
        self.pending, self.pred_err = None, []
        self.pred_err_t, self.infeasible_t = [], []
        self.n_dec = self.n_infeasible = self.n_fallback = 0
        self.decision_ms = []

    def push(self, t, tj, tr):
        self.win.append((t, tj, tr))

    def power_cap(self, t, E):
        if not self.sup:
            return 1e12
        c, p = self.c, self.p
        if env_eclipse(t + 1e-6, p):
            tau = t % p["t_orb"]
            t_rem = max(0.5 * p["t_orb"] + 0.5 * p["t_ecl"] - tau, 1.0)
            avail = max(0.0, E - c["soc_min"] * self.ecap) * 3600.0 / t_rem
            return max(0.0, (avail * c["eta_dis"] * c["eta_bd"] - c["p_aux"]) * c["eta_dab"] * c["eta_pol"])
        return self.sun

    def decide(self, t, y_true, E, f_prev, f_dem, is_ctrl):
        c, n = self.c, self.c["n_mod"]
        if self.kind == "S0":
            return c["f_max"]
        if self.kind == "S1":
            return self.f_s1[1] if env_eclipse(t + 1e-6, self.p) else self.f_s1[0]
        if not is_ctrl:
            return f_prev
        t0 = _now()
        tj_k, tr_k = self.win[-1][1], self.win[-1][2]
        tj_c = tj_k - T0K
        if self.pending is not None:
            self.pred_err.append(self.pending - float(y_true[0]))   # > 0: superestimou T_j (lado seguro); < 0: subestimou
            self.pred_err_t.append(float(t))
            self.pending = None
        cap = self.power_cap(t, E)
        f_cap = ia_f_for_power(cap / n, tj_c)[0] if cap < 1e11 else c["f_max"]
        f_hi = float(max(c["f_min"], min(c["f_max"], f_dem, f_cap)))
        if self.kind == "R":
            e = (c["tj_lim_c"] - self.r_margin) - tj_c
            df = self.ki * self.dtc * e + (self.kp * (e - self.e_prev) if self.e_prev is not None else 0.0)
            f = float(min(max(self.f_r + df, c["f_min"]), f_hi))
            self.f_r, self.e_prev = f, e
        else:
            f = self.tems(t, y_true, E, f_hi, tj_k, tr_k)
        self.n_dec += 1
        self.f_prev_applied = min(f, f_dem)
        self.win = [self.win[-1]]
        self.decision_ms.append((_now() - t0) * 1e3)
        return f

    def estimate_now(self, pr, t):
        c = self.c
        if pr.uses_truth or self.f_prev_applied is None or len(self.win) < 3:
            return self.dtip
        t0w, tj0, tr0 = self.win[0]
        if orbit_break_times(self.p, t0w, t):
            return self.dtip
        obs = {"t": np.array([s[0] - t0w for s in self.win]), "tj": np.array([s[1] for s in self.win]), "tr": np.array([s[2] for s in self.win])}
        tm = 0.5 * (t0w + t)
        qa = env_q_abs(tm, self.p)
        qr = q_rest_model(c, self.ecap, not env_eclipse(tm, self.p), c["n_mod"] * ia_p(self.f_prev_applied, tj0 - T0K)[0])
        d0 = pr.estimate(obs, self.f_prev_applied, qa, qr, tj0, tr0, self.sj, self.sr)
        w = pr.window(pr.init_state(tj0, tr0, d0, None), self.f_prev_applied, qa, qr, t - t0w, (t - t0w) / (len(obs["t"]) - 1))
        self.dtip = float(min(max(pr.dtip_of(w["end"]), DTIP_RANGE[0]), DTIP_RANGE[1]))
        return self.dtip

    def horizon(self, pr, state, t0, f1, f2, tj_c):
        c, p, W = self.c, self.p, self.W
        d1 = min(self.dtc, W)
        cuts = sorted(set([t0, t0 + W, t0 + d1] + orbit_break_times(p, t0, t0 + W)))
        peak, tjn = float(state[0]), None
        for a, b in zip(cuts[:-1], cuts[1:]):
            if b - a < 1e-6:
                continue
            f = f1 if a < t0 + d1 - 1e-6 else f2
            tm = 0.5 * (a + b)
            qr = q_rest_model(c, self.ecap, not env_eclipse(tm, p), c["n_mod"] * ia_p(f, tj_c)[0])
            w = pr.window(state, f, env_q_abs(tm, p), qr, b - a, self.dt_pred)
            peak = max(peak, float(np.max(w["tj"])))
            state = w["end"]
            if abs(b - (t0 + d1)) < 1e-6:
                tjn = float(w["tj"][-1])
        return peak, tjn

    def tems(self, t, y_true, E, f_hi, tj_k, tr_k):
        c, pr = self.c, self.pred
        f_lo, tj_c = c["f_min"], tj_k - T0K
        if f_hi <= f_lo + 1e-9:
            return f_lo
        dtip = self.estimate_now(pr, t)
        tm = t + 0.5 * self.dtc
        qa = env_q_abs(tm, self.p)
        qr = q_rest_model(c, self.ecap, not env_eclipse(tm, self.p), c["n_mod"] * ia_p(f_hi, tj_c)[0])
        d = pr.domain_distance((f_hi, qa, qr, tj_k, tr_k, dtip))
        m_pred = self.m_pred
        if d > OOD_CORTES[1] and self.fb is not None:
            pr = self.fb
            self.n_fallback += 1
            dtip = self.estimate_now(pr, t)
        elif d > OOD_CORTES[0]:
            m_pred += M_OOD_MAX_K * (d - OOD_CORTES[0]) / (OOD_CORTES[1] - OOD_CORTES[0])
        state = pr.init_state(tj_k, tr_k, dtip, y_true)
        m_tot = self.m_total if self.m_total is not None else m_pred + self.m_sensor
        limit = T0K + c["tj_lim_c"] - m_tot

        def hp(f1, f2):
            return self.horizon(pr, state, t, f1, f2, tj_c)
        if hp(f_lo, f_lo)[0] > limit:
            self.n_infeasible += 1
            self.infeasible_t.append(float(t))
            return f_lo
        f1 = bisect_max(lambda f: hp(f, f)[0] <= limit, f_lo, f_hi)
        if self.boost and f1 < f_hi:
            best = (f1 * self.dtc + f1 * (self.W - self.dtc), f1)
            for fc in np.linspace(f1, f_hi, 4)[1:]:
                if hp(fc, f_lo)[0] > limit:
                    continue
                f2c = bisect_max(lambda f: hp(fc, f)[0] <= limit, f_lo, f_hi, 5)
                if fc * self.dtc + f2c * (self.W - self.dtc) > best[0]:
                    best = (fc * self.dtc + f2c * (self.W - self.dtc), fc)
            f1 = best[1]
        self.pending = hp(f1, f1)[1]
        return float(f1)


# ---------------- malha fechada
def advance(pl, y, ta, tb, d, theta, first):
    pl.drive = static_drive(d["f"], d["q_abs"], d["q_rest"], 1.0)
    subs = [(ta, 0.5 * (ta + tb), 1.0), (0.5 * (ta + tb), tb, 1.0)] if first else [(ta, tb, theta)]
    res = fabs = 0.0
    its = 0
    for sa, sb, th in subs:
        y1, it = pl.step(y, sa, sb, th)
        qi0, qo0, _ = pl._flows(y, sa + 1e-7)
        qi1, qo1, _ = pl._flows(y1, sb - 1e-7)
        flow = (sb - sa) * (th * (qi1 - qo1) + (1.0 - th) * (qi0 - qo0))
        res += (pl._stored(y1) - pl._stored(y)) - flow
        fabs += abs(flow)
        y, its = y1, its + it
    return y, res, fabs, its


def run_closed_loop(cfg):
    c = cfg.get("contrato", PROJETO_1P0)
    par_true, par_nom = plant_par(c, cfg.get("plant_par")), plant_par(c)
    p_true, p_nom = _merge(par_true), _merge(par_nom)
    a_pl = cfg.get("a_pl", A_PL_PROJETO)
    n, dt, dtc, T = c["n_mod"], c["dt_sim"], c["dt_ctrl"], p_true["t_orb"]
    t_end = float(cfg.get("t_end", cfg.get("n_orb", 3) * T))
    s_bess, ecap = battery_for(c)
    pl = Plant(a_pl, par_true, static_drive(c["f_max"], 0.0, 0.0))
    y = np.array(cfg["y0"], float) if cfg.get("y0") is not None else np.concatenate(([320.0], np.full(pl.N + 1, 300.0)))
    soc0 = cfg.get("soc0", c["soc_max"])
    E = soc0 * ecap
    sens = Sensors(cfg.get("sensors") or sensor_cfg("ideal"), cfg.get("seed", 0))
    ctrl = Controller(cfg, c, p_nom, a_pl, ecap)
    dem = cfg.get("demand")
    ev = [0.0, t_end] + [x for x in np.arange(0.0, t_end, dtc)] + [b for b in orbit_breaks(p_true, int(t_end / T) + 2) if b < t_end] + [x for x in cfg.get("event_times", ()) if 0.0 < x < t_end]
    ev = sorted(set(round(float(x), 6) for x in ev))
    ctrl.push(0.0, *sens.read(y[0], y[1]))
    rec = {k: [] for k in ("t", "h", "tj", "tr", "tt", "soc", "f", "p_ti", "pv_used", "shed")}
    f_cmd, max_res, sum_abs, steps = c["f_max"], 0.0, 0.0, 0
    bres = 0.0
    for a, b in zip(ev[:-1], ev[1:]):
        is_ctrl = abs(a / dtc - round(a / dtc)) < 1e-9
        f_dem = float(dem(a + 1e-6)) if dem else c["f_max"]
        f_cmd = ctrl.decide(a, y, E, f_cmd, f_dem, is_ctrl)
        ns = ceil_div(b - a, dt)
        h = (b - a) / ns
        for i in range(ns):
            ta, tb = a + i * h, a + (i + 1) * h
            tm = 0.5 * (ta + tb)
            f_dem = float(dem(tm)) if dem else c["f_max"]
            f_set = float(min(max(min(f_cmd, f_dem), c["f_min"]), c["f_max"]))
            ecl = env_eclipse(tm, p_true)
            tjc = y[0] - T0K
            p_req = n * ia_p(f_set, tjc)[0]
            bs = bus_step(c, 0.0 if ecl else c["p_pv"], p_req, h, E, ecap, ctrl.floor)
            bres = max(bres, abs(bs["resid"]))
            f_eff = f_set if bs["ti"] >= p_req - 1e-6 else float(ia_f_for_power(bs["ti"] / n, tjc)[0])
            d = {"f": f_eff, "q_abs": env_q_abs(tm, p_true), "q_rest": bs["q_rest"]}
            y, res, fabs, its = advance(pl, y, ta, tb, d, 0.5, i == 0)
            E = bs["E"]
            max_res, sum_abs, steps = max(max_res, abs(res)), sum_abs + fabs, steps + 1
            ctrl.push(tb, *sens.read(y[0], y[1]))
            for k, v in (("t", tb), ("h", h), ("tj", y[0]), ("tr", y[1]), ("tt", y[-1]), ("soc", E / ecap), ("f", f_eff), ("p_ti", bs["ti"]),
                         ("pv_used", bs["used_out"]), ("shed", max(0.0, p_req - bs["ti"]))):
                rec[k].append(float(v))
    out = summarize(c, rec, ecap, T, t_end, ctrl, soc0, max_res, sum_abs / max(steps, 1))
    out["n_cells"], out["ecap_Wh"], out["bus_resid_W"] = s_bess["n_cells"], ecap, bres
    if cfg.get("series"):
        out["series"] = rec
    return out


def summarize(c, rec, ecap, T, t_end, ctrl, soc0, max_res, mean_flow):
    t, h = np.array(rec["t"]), np.array(rec["h"])
    tj, soc, f = np.array(rec["tj"]) - T0K, np.array(rec["soc"]), np.array(rec["f"])
    last = t > t_end - T + 1e-9 if t_end >= T else np.ones_like(t, bool)
    i0 = int(np.argmax(last))
    s0 = soc[i0 - 1] if i0 > 0 else soc0
    lim = c["tj_lim_c"]
    ti_wh = float(np.sum(np.array(rec["p_ti"])[last] * h[last]) / 3600.0)
    pv_wh = float(np.sum(np.array(rec["pv_used"])[last] * h[last]) / 3600.0)
    de = float((soc[-1] - s0) * ecap)
    pe = np.array(ctrl.pred_err) if ctrl.pred_err else np.zeros(1)
    pet = np.array(ctrl.pred_err_t) if ctrl.pred_err_t else np.zeros(0)
    pe_ult = pe[pet > t_end - T + 1e-9] if (len(pet) and len(pet) == len(ctrl.pred_err)) else np.zeros(0)   # erro de predicao so na ULTIMA orbita (fora do transitorio de partida)
    def janela(k):
        """PUE e variacao de energia do BESS nas k ultimas orbitas. Exige t_end >= (k+1) T: a janela precisa deixar >= 1 orbita de aquecimento ANTES dela,
        senao o transitorio de partida (SoC 0,95 -> ~0,84) entra na conta e o criterio nao mede regime."""
        if t_end < (k + 1) * T - 1e-6:
            return None, None, None
        m = t > t_end - k * T + 1e-9
        j0 = int(np.argmax(m))
        sa = soc[j0 - 1] if j0 > 0 else soc0
        ti = float(np.sum(np.array(rec["p_ti"])[m] * h[m]) / 3600.0)
        pv = float(np.sum(np.array(rec["pv_used"])[m] * h[m]) / 3600.0)
        d = float((soc[-1] - sa) * ecap)
        return (pv / ti if ti > 0 else None), d, (abs(d) / ti if ti > 0 else None)
    pue3, de3, rel3 = janela(3)
    rel1 = abs(de) / ti_wh if ti_wh > 0 else None
    return {"tj_max_C": float(tj.max()), "tj_max_last_C": float(tj[last].max()), "t_above_80_s": float(np.sum(h[tj > lim])), "t_above_85_s": float(np.sum(h[tj > 85.0])),
            "viol_int_Ks": float(np.sum(np.maximum(tj - lim, 0.0) * h)),
            "soc_min": float(soc.min()), "soc_max": float(soc.max()), "soc_min_last": float(soc[last].min()),
            "t_soc_below_min_s": float(np.sum(h[soc < c["soc_min"] - 1e-6])), "shed_Wh": float(np.sum(np.array(rec["shed"]) * h) / 3600.0),
            "tops_frac_mean": float(np.sum(f * h) / np.sum(h) / c["f_max"]), "tops_frac_last": float(np.sum(f[last] * h[last]) / np.sum(h[last]) / c["f_max"]),
            "pue_last": (pv_wh / ti_wh) if ti_wh > 0 else None, "dE_bess_last_Wh": de, "pue_periodico": bool(abs(de) / ecap < 1e-3),
            "energy_resid_rel": float(max_res / max(mean_flow, 1e-300)), "pred_err_rms_K": float(np.sqrt(np.mean(pe ** 2))), "pred_err_max_K": float(np.max(np.abs(pe))),
            "pred_err_min_signed_K": float(np.min(pe)), "pred_err_max_signed_K": float(np.max(pe)),
            "pred_err_rms_last_K": (float(np.sqrt(np.mean(pe_ult ** 2))) if len(pe_ult) else None), "pred_err_min_last_K": (float(pe_ult.min()) if len(pe_ult) else None),
            "pred_err_t_worst_s": (float(ctrl.pred_err_t[int(np.argmax(np.abs(pe)))]) if ctrl.pred_err_t else None), "infeasible_t_s": list(ctrl.infeasible_t),
            "dE_over_ti_last": rel1, "pue_last3": pue3, "dE_bess_last3_Wh": de3, "dE_over_ti_last3": rel3,
            "n_dec": ctrl.n_dec, "n_infeasible": ctrl.n_infeasible, "n_fallback": ctrl.n_fallback,
            "decision_ms_mean": float(np.mean(ctrl.decision_ms)) if ctrl.decision_ms else 0.0,
            "decision_ms_p99": float(np.percentile(ctrl.decision_ms, 99)) if ctrl.decision_ms else 0.0,
            "f_min_MHz": float(f.min()), "f_max_MHz": float(f.max())}


def make_pred(kind, c, extra_true, a_pl, ckpt=None):
    par_nom, par_true = plant_par(c), plant_par(c, extra_true)
    if kind == "RC3":
        return RC3(par_nom, a_pl, "euler")
    if kind == "RC3R":
        return RC3(par_nom, a_pl, "radau")
    if kind == "SOLVER":
        return SolverPred(par_true, a_pl)
    if kind == "PINN":
        if ckpt is None:
            raise ValueError("preditor PINN exige --ckpt (checkpoint com hash travado no contrato)")
        return PINNPred(ckpt)
    raise ValueError("preditor desconhecido: %s" % kind)


CONTROLADORES = {"S0": ("S0", None), "S1": ("S1", None), "R": ("R", None), "R_SINT": ("R", None), "TEMS_PINN": ("TEMS", "PINN"), "TEMS_RC3": ("TEMS", "RC3"),
                 "TEMS_RC3R": ("TEMS", "RC3R"), "TEMS_SOLVER": ("TEMS", "SOLVER")}


def f_for_ti(c, p_ti, tj_c=80.0):
    return float(ia_f_for_power(p_ti / c["n_mod"], tj_c)[0])


def make_demand(c, seed, t_end, per_orbit=4):
    rng = np.random.default_rng(1000 + int(seed))
    starts = np.sort(rng.uniform(0.0, max(t_end - 400.0, 1.0), size=max(1, int(per_orbit * t_end / T_ORB))))
    f_hi, f_lo = c["f_max"], f_for_ti(c, 20000.0)
    evs = [float(s) for s in starts] + [float(s) + 300.0 for s in starts]
    return (lambda t: f_hi if any(s <= t < s + 300.0 for s in starts) else f_lo), evs


PERTURBACOES = {"r_int_alto": {"r_int": 0.12}, "r_int_muito_alto": {"r_int": 0.16}, "c_mod_menos30": {"c_mod": 350.0}, "c_mod_mais30": {"c_mod": 650.0},
                "eps_baixo": {"alpha": 0.36, "eps": 0.80}, "combinada": {"r_int": 0.12, "c_mod": 650.0, "alpha": 0.36, "eps": 0.80}}   # niveis [HIPOTESE DECLARADA]


def scenario_cfg(nome, seed, c, n_orb, sensor="ideal", burst_per_orbit=4, **sens_over):
    t_end = n_orb * T_ORB
    cfg = {"contrato": c, "n_orb": n_orb, "seed": seed, "sensors": sensor_cfg(sensor, **sens_over), "m_sensor": SENSOR_MARGEM_K[sensor]}
    if nome == "burst":
        cfg["demand"], cfg["event_times"] = make_demand(c, seed, t_end, burst_per_orbit)
    elif nome in PERTURBACOES:
        cfg["plant_par"] = PERTURBACOES[nome]
    elif nome == "ood":
        cfg["y0"] = np.concatenate(([345.0], np.full(41, 290.0)))
    elif nome != "nominal":
        raise ValueError("cenario desconhecido: %s" % nome)
    return cfg


def run_scenario(nome, ctrl_label, seed, c=PROJETO_1P0, n_orb=3, sensor="ideal", ckpt=None, m_total=None, blocks_boost=False, cfg_over=None, burst_per_orbit=4, series=False, **sens_over):
    """cfg_over: sobrescreve chaves do cfg (ex.: {"r_kp":..., "r_ki":..., "r_margin":...} para o reativo)."""
    cfg = scenario_cfg(nome, seed, c, n_orb, sensor, burst_per_orbit=burst_per_orbit, **sens_over)
    if series:
        cfg["series"] = True
    kind, pk = CONTROLADORES[ctrl_label]
    cfg["ctrl"], cfg["boost"], cfg["m_total"] = kind, blocks_boost, m_total
    if cfg_over:
        cfg.update(cfg_over)
    if pk:
        cfg["pred"] = make_pred(pk, c, cfg.get("plant_par") if pk == "SOLVER" else None, cfg.get("a_pl", A_PL_PROJETO), ckpt)
        if pk == "PINN":
            cfg["fallback"] = make_pred("RC3", c, None, A_PL_PROJETO)
    return run_closed_loop(cfg)


def ultima_orbita(s, t_end):
    """Series da ULTIMA orbita (para as figuras do Capitulo 4): t_s, tj_C, tr_C, f_MHz, soc."""
    t = np.array(s["t"])
    m = t > t_end - T_ORB + 1e-9
    return {"t_s": t[m].tolist(), "tj_C": (np.array(s["tj"])[m] - T0K).tolist(), "tr_C": (np.array(s["tr"])[m] - T0K).tolist(),
            "f_MHz": np.array(s["f"])[m].tolist(), "soc": np.array(s["soc"])[m].tolist()}


def chave_sintonia(sensor, quant="none"):
    """Chave da sintonia do reativo: o R e sintonizado POR (sensor, quantizacao), porque a margem otima depende de ambos."""
    return "%s|%s" % (sensor, quant)


def carrega_sintonia(path, seeds_teste, sensor, quant="none"):
    """Le o JSON da sintonia do reativo e valida (ver valida_sintonia)."""
    with open(path, encoding="utf-8") as fh:
        t = json.load(fh)
    return valida_sintonia(t, path, seeds_teste, sensor, quant)


def valida_sintonia(t, path, seeds_teste, sensor, quant="none"):
    """RECUSA: contrato diferente, (sensor, quant) ausente/inviavel, ou sementes de teste que coincidam com as de calibracao."""
    if t.get("contrato_sha256") != contract_hash(PROJETO_1P0):
        raise RuntimeError("sintonia do R feita com outro contrato (hash diferente): refazer sintonia_reativo.py")
    comum = set(int(s) for s in seeds_teste) & set(int(s) for s in t["seeds_cal"])
    if comum:
        raise RuntimeError("sementes de teste coincidem com as de calibracao do R: %s" % sorted(comum))
    k = chave_sintonia(sensor, quant)
    s = t.get("sensores", {}).get(k)
    if s is None or not s.get("viavel"):
        raise RuntimeError("sem sintonia viavel do R para '%s' em %s" % (k, path))
    return {"r_kp": s["kp"], "r_ki": s["ki"], "r_margin": s["margem"]}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Marco 5: malha fechada TEMS (DVFS-only)")
    ap.add_argument("--scenarios", nargs="+", default=["nominal"])
    ap.add_argument("--ctrls", nargs="+", default=["S0", "S1", "R", "TEMS_RC3"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[0])
    ap.add_argument("--sensor", default="sem_calibracao", choices=sorted(SENSOR_MARGEM_K))
    ap.add_argument("--norb", type=int, default=3)
    ap.add_argument("--ckpt", default=None)
    ap.add_argument("--quant", default="none", choices=["none", "floor", "trunc", "round"], help="quantizacao inteira (1 C) da leitura de T_j (NVML); trunc = floor")
    ap.add_argument("--rtune", default=None, help="JSON de sintonia_reativo.py (obrigatorio para o controlador R_SINT)")
    ap.add_argument("--mtotal", "--m_t", dest="mtotal", type=float, default=None, help="margem TOTAL de T_j [K] do TEMS (substitui m_pred + m_sensor); usado na varredura Pareto do cenario 6")
    ap.add_argument("--dt_ctrl", type=float, default=60.0, choices=[60.0, 300.0], help="periodo de decisao [s]; 300 e uma VARIANTE do contrato (hash diferente)")
    ap.add_argument("--burst_per_orbit", type=int, default=4, help="rajadas de 300 s por orbita no cenario burst")
    ap.add_argument("--save_traces", default=None, help="JSON com as series da ULTIMA orbita (t_s, tj_C, tr_C, f_MHz, soc)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    c_run = PROJETO_1P0 if a.dt_ctrl == 60.0 else dict(PROJETO_1P0, nome="PROJETO_1P0_dt%d" % int(a.dt_ctrl), dt_ctrl=a.dt_ctrl)
    if "R_SINT" in a.ctrls and a.dt_ctrl != 60.0:
        raise SystemExit("R_SINT so tem sintonia para dt_ctrl = 60 s")
    over = None
    if "R_SINT" in a.ctrls:
        if not a.rtune:
            raise SystemExit("R_SINT exige --rtune sintonia_R.json (gerado por sintonia_reativo.py)")
        over = carrega_sintonia(a.rtune, a.seeds, a.sensor, a.quant)
    qkw = {} if a.quant == "none" else {"quant_mode": a.quant}
    res, traces = [], []
    for sc in a.scenarios:
        for ct in a.ctrls:
            for sd in a.seeds:
                r = run_scenario(sc, ct, sd, c=c_run, n_orb=a.norb, sensor=a.sensor, ckpt=a.ckpt, m_total=a.mtotal, cfg_over=over if ct == "R_SINT" else None,
                                 burst_per_orbit=a.burst_per_orbit, series=bool(a.save_traces), **qkw)
                if a.save_traces:
                    traces.append(dict(cenario=sc, ctrl=ct, seed=sd, sensor=a.sensor, quant=a.quant, **ultima_orbita(r.pop("series"), a.norb * T_ORB)))
                res.append(dict(cenario=sc, ctrl=ct, seed=sd, sensor=a.sensor, quant=a.quant, mtotal=a.mtotal, dt_ctrl=a.dt_ctrl, burst_per_orbit=a.burst_per_orbit, **r))
    payload = {"contrato": c_run["nome"], "contrato_sha256": contract_hash(c_run),
               "ckpt_sha256": sha256_file(a.ckpt) if a.ckpt else None,
               "rtune_sha256": sha256_file(a.rtune) if a.rtune else None, "rtune_usado": over, "resultados": res}
    txt = json.dumps(_r(payload), ensure_ascii=False)
    if a.out:
        with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(txt)
    if a.save_traces:
        with open(a.save_traces, "w", encoding="utf-8", newline="\n") as fh:
            json.dump(_r({"contrato": c_run["nome"], "contrato_sha256": contract_hash(c_run), "traces": traces}), fh, ensure_ascii=False)
    print(txt)
    return payload


if __name__ == "__main__":
    main()

"""modelo_termico.py -- Marco 3: planta termica de referencia no vacuo.
Planta = (i) no concentrado de MODULO (juncao, C_mod, R_int ate o coletor) acoplado a (ii) coletor/fluido (C_fluid) na raiz de
(iii) ALETA DE RADIADOR 1D (EDP do calor), meia-celula entre tubo (x=0) e plano de simetria adiabatico (x=L), duas faces radiando:
    c_a dT/dt = d/dx( k delta dT/dx ) - 2 eps sigma (T^4 - Ts^4) + q_abs(t)         0 < x < L     [W/m2 de planta]
    x=0: no do coletor (fluido + modulos);   x=L: dT/dx = 0
O contorno 3D  -k grad T . n = eps sigma (T^4 - Ts^4) - q_amb  nas faces e integrado na espessura (aleta fina, Bi << 1) e vira o sorvedouro/fonte acima.
Numerica: volumes finitos centrados em vertices (meias celulas nas pontas) + metodo theta (theta=0.5 Crank-Nicolson) + Newton, solver tridiagonal proprio.
Etiquetas: [FONTE-RESUMO] visto so em resumo de busca (PDF NAO aberto) | [CALCULADO] | [HIPOTESE DECLARADA] (com sensibilidade) | [PENDENTE] sem fonte
| [NAO VERIFICADO] memoria | [LEGADO] valor do legado, sem fonte | [USUARIO] confirmado localmente pelo usuario.
Autocontido (numpy). Se modelo_ia.py / modelo_barramento.py estiverem na pasta, run_selftests() faz conferencia cruzada. Nao validado em TVAC/orbita.
"""
import math
import json
import numpy as np

SIGMA = 5.670374419e-8      # W/(m2 K4)  [NAO VERIFICADO] memoria (CODATA 2018)
T0K = 273.15
_EPS = 1.0e-7               # s: avaliacao das forcantes DENTRO do intervalo (limites laterais nas descontinuidades)

# ---- Marco 1 [USUARIO] (LTAN 18:00, 21/06/2026)
T_ORB, T_ECL, BETA_DEG = 5926.379071134441, 1105.3908, 58.37416
F_VIEW_90, G_SUN, ALBEDO, Q_IR = 0.232871441130663, 1361.0, 0.30, 237.0   # F(90 deg) [USUARIO]; G_SUN IAU; albedo e Q_IR NASA SP-8067 [FONTE-RESUMO]

# ---- revestimentos do radiador (alpha_s, eps_IR)
COATINGS = {
    "Z93_EOL30": (0.36, 0.90),   # [FONTE-RESUMO] NTRS AIAA 92-2168 (reg. 19920048673): Z-93, estimativa EOL de 30 anos em LEO. USADO como limite conservador p/ 10 anos [HIPOTESE DECLARADA]
    "Z93_BOL": (0.15, 0.92),     # [FONTE-RESUMO] idem, BOL
    "AgFEP_BOL": (0.085, 0.76),  # [FONTE-RESUMO] NTRS 19740009173 (1974) Teflon FEP prateado comercial; EOL em LEO NAO tem numero lido (LDEF: aumento local de alpha ate ~0.25)
    "LEGADO": (0.14, 0.88),      # [LEGADO] 4_1, sem fonte
}

PAR = {
    "coating": "Z93_EOL30", "alpha": None, "eps": None,
    "k": 167.0,           # W/(m K) [NAO VERIFICADO] aluminio, memoria
    "delta": 2.0e-3,      # m [LEGADO 4_1 fin_efficiency_diagnostic]
    "L": 0.075,           # m [LEGADO 4_1] meia distancia entre tubos
    "c_area": 4050.0,     # J/(m2 K) [LEGADO 4_1: 4.5 kg/m2 x 900 J/kgK]
    "n_nodes": 40,        # [HIPOTESE DECLARADA] malha; convergencia medida em report_conv
    "r_int": 0.08,        # K/W por modulo, juncao -> coletor [HIPOTESE DECLARADA] sem fonte; sensibilidade 0.04..0.16
    "c_mod": 500.0,       # J/K por modulo (juncao+placa fria) [HIPOTESE DECLARADA]
    "n_mod": 100,         # 100 x 400 W = 40 kW [USUARIO]
    "c_fluid": 5.0e4,     # J/K, fluido + coletor [HIPOTESE DECLARADA]
    "kappa_dab": 1.0 / 0.9655 - 1.0,   # perda do DAB / P_TI (eta_pol,ext = 1) [HIPOTESE DECLARADA, 2b]
    "t_space": 3.0,       # K [enunciado do projeto]
    "t_orb": T_ORB, "t_ecl": T_ECL, "beta_deg": BETA_DEG,
    "q_ir": Q_IR, "albedo": ALBEDO, "g_sun": G_SUN, "f_view": F_VIEW_90,
    "linear_h": 0.0, "t_inf": 0.0,    # modo linear SO para testes analiticos
}

# ---- modulo de IA (copia de modelo_ia.py DEFAULTS; conferida por teste: P(1410 MHz, 80 C) = 400 W)
IA = {"tdp": 400.0, "f_nom": 1410.0, "f_min": 300.0, "v_nom": 0.80, "v_min": 0.65, "fv": 0.40, "fix": 0.20, "leak": 0.20, "t_ref": 80.0, "beta": 0.02, "kv": 3.0}


def ia_v(f):
    fv = IA["fv"] * IA["f_nom"]
    return IA["v_min"] if f <= fv else IA["v_min"] + (IA["v_nom"] - IA["v_min"]) * (f - fv) / (IA["f_nom"] - fv)


def ia_p(f, tj_c, util=1.0):
    """P = P_fix + P_dyn + V I_leak(Tj). Retorna (P, P_leak) em W por modulo."""
    if not (IA["f_min"] - 1e-9 <= f <= IA["f_nom"] + 1e-9) or not (0.0 <= util <= 1.0):
        raise ValueError("f fora de [f_min, f_nom] ou util fora de [0, 1]")
    v = ia_v(f)
    p_dyn = util * (1.0 - IA["fix"] - IA["leak"]) * IA["tdp"] * (v / IA["v_nom"]) ** 2 * (f / IA["f_nom"])
    i_leak = IA["leak"] * IA["tdp"] / IA["v_nom"] * (v / IA["v_nom"]) ** (IA["kv"] - 1.0) * math.exp(IA["beta"] * (tj_c - IA["t_ref"]))
    return IA["fix"] * IA["tdp"] + p_dyn + v * i_leak, v * i_leak


def ia_f_for_power(p_t, tj_c, util=1.0):
    lo, hi = IA["f_min"], IA["f_nom"]
    if p_t <= ia_p(lo, tj_c, util)[0]:
        return lo, "piso"
    if p_t >= ia_p(hi, tj_c, util)[0]:
        return hi, "teto"
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if ia_p(mid, tj_c, util)[0] < p_t:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi), "ok"


# ---- barramento e BESS (copia fiel das equacoes de modelo_barramento.py e modelo_bess.py; eta_pol,ext = 1.0 por decisao do usuario)
BUS = {"eta_buck": 0.91, "eta_dab": 0.9655, "eta_pol": 1.0, "eta_bd": 0.97, "p_aux": 1000.0, "p_pv": 50000.0}
BESS = {"eta_ch": 0.995, "eta_dis": 0.97, "fade": 0.20, "n_ser": 108, "e_cell": 16.2, "m_cell": 0.115, "soc_min": 0.70, "soc_max": 0.95}
AUX_HEAT_FRAC = 903.2 / 1000.0   # [LEGADO 4_1] fracao de P_aux que vira calor
Q_OISL = 343.2                   # W [LEGADO 4_1]


def bus_heat(p_pv_out, p_ti, **kw):
    """Calor no no (pior caso: arranjo processa toda a potencia; BESS aceita o excedente). Mesma estrutura de bus_flow() do 2b."""
    b = dict(BUS)
    b.update(kw)
    for k in ("eta_buck", "eta_dab", "eta_pol", "eta_bd"):
        if not (0.0 < b[k] <= 1.0):
            raise ValueError("%s em (0, 1]" % k)
    if not (math.isfinite(p_pv_out) and p_pv_out >= 0.0 and math.isfinite(p_ti) and p_ti >= 0.0):
        raise ValueError("potencias finitas e >= 0")
    p48 = p_ti / b["eta_pol"]
    p400 = p48 / b["eta_dab"]
    l_bus = p400 + b["p_aux"]
    pv_bus = b["eta_buck"] * p_pv_out
    net = pv_bus - l_bus
    q_bd = net * (1.0 - b["eta_bd"]) if net >= 0.0 else (-net) / b["eta_bd"] - (-net)
    q_buck = p_pv_out - pv_bus
    heat_aux = b["p_aux"] * AUX_HEAT_FRAC
    return {"p_bus": p400, "q_pol": p48 - p_ti, "q_dab": p400 - p48, "heat_aux": heat_aux, "q_buck": q_buck, "q_bd": q_bd, "l_bus": l_bus, "net_bus": net,
            "q_node": p_ti + (p48 - p_ti) + (p400 - p48) + heat_aux + q_buck + q_bd}


def analytic_orbit(t_ecl, T, p_ti_ecl, **kw):
    q = dict(BUS, **BESS)
    bad = set(kw) - set(q)
    if bad:
        raise ValueError("parametros desconhecidos: %s" % sorted(bad))
    q.update(kw)
    if not (0.0 <= t_ecl < T) or not (math.isfinite(p_ti_ecl) and p_ti_ecl >= 0.0):
        raise ValueError("t_ecl em [0, T) e p_ti_ecl >= 0")
    chain = q["eta_dab"] * q["eta_pol"]
    l_ecl = p_ti_ecl / chain + q["p_aux"]
    s_out = l_ecl / q["eta_bd"]
    e_chem = s_out * t_ecl / q["eta_dis"] / 3600.0
    surplus = l_ecl * t_ecl / (q["eta_ch"] * q["eta_dis"] * q["eta_bd"] ** 2 * (T - t_ecl))
    l_sun = q["eta_buck"] * q["p_pv"] - surplus
    return {"L_ecl_bus_W": l_ecl, "E_chem_out_Wh": e_chem, "P_TI_sun_sust_W": max(0.0, (l_sun - q["p_aux"]) * chain), "surplus_bus_W": surplus}


def size_bess(t_ecl, T, p_ti_ecl, **kw):
    q = dict(BUS, **BESS)
    q.update({k: v for k, v in kw.items() if k in q})
    if not (0.0 <= q["fade"] < 1.0):
        raise ValueError("fade em [0, 1)")
    a = analytic_orbit(t_ecl, T, p_ti_ecl, **kw)
    e_eol = a["E_chem_out_Wh"] / (q["soc_max"] - q["soc_min"])
    e_bol = e_eol / (1.0 - q["fade"])
    n_par_real = e_bol / (q["n_ser"] * q["e_cell"])
    n_par = max(1, math.ceil(n_par_real - 1e-12))
    n_cells = q["n_ser"] * n_par
    return {"n_par": n_par, "n_cells": n_cells, "mass_cells_kg": n_cells * q["m_cell"], "E_inst_BOL_Wh": n_cells * q["e_cell"], "E_req_EOL_Wh": e_eol,
            "dod_eol": a["E_chem_out_Wh"] / (n_cells * q["e_cell"] * (1.0 - q["fade"])), "P_TI_sun_sust_W": a["P_TI_sun_sust_W"], "L_ecl_bus_W": a["L_ecl_bus_W"]}


# ---- parametros
def _merge(par=None):
    p = dict(PAR)
    if par:
        bad = set(par) - set(PAR)
        if bad:
            raise ValueError("parametros desconhecidos: %s" % sorted(bad))
        p.update(par)
    for k in ("k", "delta", "L", "c_area", "r_int", "c_mod", "g_sun"):
        if not (math.isfinite(p[k]) and p[k] > 0.0):
            raise ValueError("%s deve ser finito e > 0" % k)
    if not (isinstance(p["n_nodes"], int) and p["n_nodes"] >= 2 and isinstance(p["n_mod"], int) and p["n_mod"] >= 0):
        raise ValueError("n_nodes >= 2 (int) e n_mod >= 0 (int)")
    if not (p["c_fluid"] >= 0.0 and p["t_space"] >= 0.0 and p["kappa_dab"] >= 0.0 and 0.0 < p["t_ecl"] < 0.5 * p["t_orb"] and 0.0 <= p["beta_deg"] < 90.0):
        raise ValueError("c_fluid, t_space, kappa_dab >= 0; 0 < t_ecl < t_orb/2; 0 <= beta < 90")
    if p["alpha"] is None or p["eps"] is None:
        if p["coating"] not in COATINGS:
            raise ValueError("coating desconhecido")
    return p


def optics(p):
    return (p["alpha"], p["eps"]) if (p["alpha"] is not None and p["eps"] is not None) else COATINGS[p["coating"]]


def static_drive(f, q_abs, q_rest, util=1.0):
    d = {"f": f, "util": util, "q_abs": q_abs, "q_rest": q_rest, "ecl": False}
    return lambda t: d


def orbit_drive(p, f_sun, f_ecl, q_rest_sun, q_rest_ecl, util=1.0):
    """Eclipse centrada em theta = pi (duracao t_ecl do Marco 1). q_abs = 2 faces x [alpha G a F max(0, cos(beta) cos(theta)) + eps Q_IR F]  (opcao 1: phi = 90 deg, sem Sol direto)."""
    alpha, eps = optics(p)
    T, te, cb = p["t_orb"], p["t_ecl"], math.cos(math.radians(p["beta_deg"]))
    q_ir = eps * p["q_ir"] * p["f_view"]

    def drive(t):
        tau = t % T
        ecl = abs(tau - 0.5 * T) < 0.5 * te
        mu = max(0.0, cb * math.cos(2.0 * math.pi * tau / T))
        return {"f": f_ecl if ecl else f_sun, "util": util, "q_abs": 2.0 * (alpha * p["g_sun"] * p["albedo"] * p["f_view"] * mu + q_ir),
                "q_rest": q_rest_ecl if ecl else q_rest_sun, "ecl": ecl}
    return drive


def orbit_breaks(p, n_orb):
    T, te = p["t_orb"], p["t_ecl"]
    return [k * T + fr for k in range(n_orb) for fr in (0.25 * T, 0.5 * T - 0.5 * te, 0.5 * T + 0.5 * te, 0.75 * T, T)]


def _thomas(sub, diag, sup, rhs):
    n = len(diag)
    sub, diag, sup, rhs = sub.tolist(), diag.tolist(), sup.tolist(), rhs.tolist()
    cp, dp = [0.0] * n, [0.0] * n
    cp[0], dp[0] = sup[0] / diag[0], rhs[0] / diag[0]
    for k in range(1, n):
        m = diag[k] - sub[k - 1] * cp[k - 1]
        cp[k] = sup[k] / m if k < n - 1 else 0.0
        dp[k] = (rhs[k] - sub[k - 1] * dp[k - 1]) / m
    x = [0.0] * n
    x[-1] = dp[-1]
    for k in range(n - 2, -1, -1):
        x[k] = dp[k] - cp[k] * x[k + 1]
    return np.array(x)


class Plant:
    """Estado y = [Tj, T_0, ..., T_N] em K (Tj = juncao de UM modulo representativo; T_0 = coletor/raiz; T_N = ponta da aleta)."""

    def __init__(self, a_pl, par=None, drive=None):
        p = _merge(par)
        if not (isinstance(a_pl, (int, float)) and math.isfinite(a_pl) and a_pl > 0.0) or drive is None:
            raise ValueError("a_pl finito > 0 e drive obrigatorio")
        self.p, self.drive, self.a_pl = p, drive, float(a_pl)
        N = self.N = p["n_nodes"]
        self.L, self.dx = p["L"], p["L"] / N
        self.G = p["k"] * p["delta"] / self.dx
        w = np.full(N + 1, self.dx)
        w[0] = w[-1] = 0.5 * self.dx
        self.w = w
        self.ell = self.a_pl / self.L                 # comprimento total de meias-celulas [m]: planta = ell x L
        self.alpha, self.eps = optics(p)
        self.eps_t = 2.0 * self.eps
        cap = np.empty(N + 2)
        cap[0] = p["c_mod"]
        cap[1:] = p["c_area"] * w
        cap[1] += p["c_fluid"] / self.ell
        self.cap = cap

    def _F(self, y, t):
        p, N, G, w = self.p, self.N, self.G, self.w
        d = self.drive(t)
        tj, T = y[0], y[1:]
        P, pleak = ia_p(d["f"], tj - T0K, d["util"])
        dP = IA["beta"] * pleak
        R, s, kap = p["r_int"], p["n_mod"] / self.ell, p["kappa_dab"]
        if p["linear_h"] > 0.0:
            rad, drad, ab = w * p["linear_h"] * (T - p["t_inf"]), w * p["linear_h"], np.zeros(N + 1)
        else:
            rad = w * self.eps_t * SIGMA * (T ** 4 - p["t_space"] ** 4)
            drad = w * self.eps_t * SIGMA * 4.0 * T ** 3
            ab = w * d["q_abs"]
        q = G * (T[:-1] - T[1:])
        net = np.zeros(N + 1)
        net[:-1] -= q
        net[1:] += q
        F = np.empty(N + 2)
        F[0] = P - (tj - T[0]) / R
        F[1:] = net - rad + ab
        F[1] += s * (tj - T[0]) / R + s * kap * P + d["q_rest"] / self.ell
        diag = np.empty(N + 2)
        sub, sup = np.full(N + 1, G), np.full(N + 1, G)
        diag[0], sup[0], sub[0] = dP - 1.0 / R, 1.0 / R, s / R + s * kap * dP
        diag[1] = -G - drad[0] - s / R
        diag[2:N + 1] = -2.0 * G - drad[1:N]
        diag[N + 1] = -G - drad[N]
        return F, diag, sub, sup

    def _flows(self, y, t):
        p, d, T, w = self.p, self.drive(t), y[1:], self.w
        P, _ = ia_p(d["f"], y[0] - T0K, d["util"])
        if p["linear_h"] > 0.0:
            rad, ab = w * p["linear_h"] * (T - p["t_inf"]), np.zeros_like(T)
        else:
            rad, ab = w * self.eps_t * SIGMA * (T ** 4 - p["t_space"] ** 4), w * d["q_abs"]
        return p["n_mod"] * P * (1.0 + p["kappa_dab"]) + d["q_rest"] + self.ell * float(ab.sum()), self.ell * float(rad.sum()), P

    def _stored(self, y):
        return self.p["n_mod"] * self.cap[0] * y[0] + self.ell * float((self.cap[1:] * y[1:]).sum())

    def step(self, y_old, t_old, t_new, theta=0.5, steady=False, tol=1e-9, itmax=60, maxstep=30.0):
        dt = t_new - t_old
        to, tn = (t_new, t_new) if steady else (t_old + _EPS, t_new - _EPS)
        fo = None if (steady or theta >= 1.0) else self._F(y_old, to)[0]
        y = np.array(y_old, float)
        for it in range(itmax):
            F, dg, sb, sp = self._F(y, tn)
            if steady:
                res, jd, js, ju = F, dg, sb, sp
            else:
                res = self.cap * (y - y_old) / dt - theta * F - ((1.0 - theta) * fo if fo is not None else 0.0)
                jd, js, ju = self.cap / dt - theta * dg, -theta * sb, -theta * sp
            dlt = _thomas(js, jd, ju, -res)
            mx = float(np.max(np.abs(dlt)))
            if not math.isfinite(mx):
                raise RuntimeError("Newton divergiu (nao finito)")
            y = y + (dlt * (maxstep / mx) if mx > maxstep else dlt)
            if mx < tol:
                return y, it + 1
        raise RuntimeError("Newton nao convergiu em %d iteracoes" % itmax)

    def steady(self, y0=None, t=0.0):
        """Estado estacionario da forcante estatica em t. Retorna None se nao houver solucao estavel (fuga termica / divergencia)."""
        y = np.array(y0, float) if y0 is not None else np.concatenate(([340.0], np.full(self.N + 1, 320.0)))

        def ok(v):
            return v is not None and bool(np.all(np.isfinite(v))) and v.max() < 800.0 and v.min() >= self.p["t_space"] - 1e-9
        try:
            y1, _ = self.step(y, t, t, steady=True, itmax=80)
        except (RuntimeError, ValueError, OverflowError):
            y1 = None
        if ok(y1):
            return y1
        yy, dtp = y.copy(), 5.0
        try:
            for _ in range(150):
                yn, _ = self.step(yy, 0.0, dtp, theta=1.0)
                done = float(np.max(np.abs(yn - yy))) < 1e-8
                yy, dtp = yn, dtp * 1.6
                if not ok(yy):
                    return None
                if done:
                    break
            y2, _ = self.step(yy, t, t, steady=True, itmax=80)
        except (RuntimeError, ValueError, OverflowError):
            return None
        return y2 if ok(y2) else None

    def simulate(self, y0, t_end, dt, theta=0.5, breaks=(), rannacher=2):
        """Integra em [0, t_end]. Crank-Nicolson NAO amortece modos rigidos (fator -> -1 quando dt >> tau): os nos da aleta tem tau ~ 0.01 s, dt ~ 25-100 s.
        Por isso os 'rannacher' primeiros passos de CADA segmento (inicio e apos cada descontinuidade da forcante) sao feitos com Euler implicito
        em dois meios passos (amortecimento forte), depois theta. rannacher=0 reproduz o CN puro (usado no teste analitico concentrado)."""
        if not (t_end > 0.0 and dt > 0.0 and 0.0 <= theta <= 1.0 and isinstance(rannacher, int) and rannacher >= 0):
            raise ValueError("t_end, dt > 0; theta em [0, 1]; rannacher inteiro >= 0")
        edges = [0.0]
        for b in sorted(breaks):
            if 1e-6 < b < t_end - 1e-6 and b - edges[-1] > 1e-6:
                edges.append(float(b))
        edges.append(float(t_end))
        y = np.array(y0, float)
        ts, tj, t0, tn = [0.0], [y[0]], [y[1]], [y[-1]]
        edge_y, max_res, sum_res, sum_flow, steps, nit = [], 0.0, 0.0, 0.0, 0, 0
        for a, b in zip(edges[:-1], edges[1:]):
            n = max(1, int(round((b - a) / dt)))
            h = (b - a) / n
            for i in range(n):
                ta, tb = a + i * h, a + (i + 1) * h
                subs = [(ta, 0.5 * (ta + tb), 1.0), (0.5 * (ta + tb), tb, 1.0)] if i < rannacher else [(ta, tb, theta)]
                for sa, sb, th in subs:
                    y1, it = self.step(y, sa, sb, th)
                    qi0, qo0, _ = self._flows(y, sa + _EPS)
                    qi1, qo1, _ = self._flows(y1, sb - _EPS)
                    flow = (sb - sa) * (th * (qi1 - qo1) + (1.0 - th) * (qi0 - qo0))
                    res = (self._stored(y1) - self._stored(y)) - flow
                    max_res, sum_res, sum_flow = max(max_res, abs(res)), sum_res + res, sum_flow + abs(flow)
                    y, steps, nit = y1, steps + 1, nit + it
                    ts.append(sb); tj.append(y[0]); t0.append(y[1]); tn.append(y[-1])
            edge_y.append(y.copy())
        return {"t": np.array(ts), "tj": np.array(tj), "t0": np.array(t0), "tN": np.array(tn), "y": y, "edge_y": edge_y, "max_step_resid_J": max_res,
                "sum_resid_J": sum_res, "rel_resid": max_res / max(sum_flow / max(steps, 1), 1e-300), "steps": steps, "newton_its": nit}


def fin_first_integral(t_root, t_tip, k_delta, eps_t):
    """Aleta radiante estacionaria, Ts = 0, ponta adiabatica: (T')^2 = c (T^5 - T_tip^5), c = 2 eps_t sigma/(5 k delta).
    L = int_{T_tip}^{T_root} dT / sqrt(c (T^5 - T_tip^5)) (substituicao T = T_tip + u^2, integrando suave);  q_raiz = k delta sqrt(c (T_root^5 - T_tip^5))."""
    c = 2.0 * eps_t * SIGMA / (5.0 * k_delta)
    xg, wg = np.polynomial.legendre.leggauss(400)
    umax = math.sqrt(t_root - t_tip)
    u, wu = 0.5 * umax * (xg + 1.0), 0.5 * umax * wg
    a, T = t_tip, t_tip + u * u
    return float((wu * 2.0 / np.sqrt(c * (T ** 4 + T ** 3 * a + T ** 2 * a * a + T * a ** 3 + a ** 4))).sum()), k_delta * math.sqrt(c * (t_root ** 5 - t_tip ** 5))


def describe(pl, y):
    d, p, T, w = pl.drive(0.0), pl.p, y[1:], pl.w
    ts = p["t_space"]
    net = float((w * (pl.eps_t * SIGMA * (T ** 4 - ts ** 4) - d["q_abs"])).sum())
    iso = pl.L * (pl.eps_t * SIGMA * (T[0] ** 4 - ts ** 4) - d["q_abs"])
    P, _ = ia_p(d["f"], y[0] - T0K, d["util"])
    q_in = p["n_mod"] * P * (1.0 + p["kappa_dab"]) + d["q_rest"]
    t_eq = (d["q_abs"] / (pl.eps_t * SIGMA) + ts ** 4) ** 0.25
    return {"a_pl_m2": pl.a_pl, "a_emit_2faces_m2": 2.0 * pl.a_pl, "tj_C": y[0] - T0K, "t_root_C": T[0] - T0K, "t_tip_C": T[-1] - T0K, "eta_fin": net / iso,
            "q_rej_W": q_in, "q_rej_Wm2_planta": q_in / pl.a_pl, "balanco_resid_W": pl.ell * net - q_in, "p_mod_W": P, "t_eq_C": t_eq - T0K,
            "R_mod_eq_K_W": (y[0] - t_eq) / P, "R_rad_eq_K_W_por_modulo": (T[0] - t_eq) / (q_in / max(p["n_mod"], 1)),
            "bi_face": 4.0 * pl.eps * SIGMA * T[0] ** 3 * 0.5 * p["delta"] / p["k"]}


def design_area(p_ti, par=None, tj_lim_c=80.0, bounds=(5.0, 600.0), p_pv=None):
    """Menor area de planta A_pl tal que Tj <= tj_lim em regime estacionario no caso quente (PV ligado, albedo de pico, f ajustada a Tj = limite)."""
    p = _merge(par)
    n = p["n_mod"]
    if n < 1 or not (math.isfinite(p_ti) and p_ti > 0.0):
        raise ValueError("n_mod >= 1 e p_ti > 0")
    f, flag = ia_f_for_power(p_ti / n, tj_lim_c)
    bh = bus_heat(BUS["p_pv"] if p_pv is None else p_pv, p_ti)
    alpha, eps = optics(p)
    q_abs = 2.0 * (alpha * p["g_sun"] * p["albedo"] * p["f_view"] * math.cos(math.radians(p["beta_deg"])) + eps * p["q_ir"] * p["f_view"])
    drive = static_drive(f, q_abs, bh["heat_aux"] + bh["q_buck"] + bh["q_bd"])
    last = {}

    def tj_at(a):
        pl = Plant(a, p, drive)
        y = pl.steady(last.get("y"))
        if y is not None:
            last["y"] = y
        return None if y is None else y[0] - T0K, pl, y
    lo, hi = bounds
    t_hi, pl_hi, y_hi = tj_at(hi)
    if t_hi is None or t_hi > tj_lim_c:
        raise RuntimeError("limite superior de area insuficiente")
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        t_mid, pl_mid, y_mid = tj_at(mid)
        if t_mid is None or t_mid > tj_lim_c:
            lo = mid
        else:
            hi, pl_hi, y_hi = mid, pl_mid, y_mid
        if hi - lo < 1e-5 * hi:
            break
    out = describe(pl_hi, y_hi)
    out.update({"p_ti_W": p_ti, "f_MHz": f, "flag_f": flag, "a_pl": hi})
    return out


# ================= testes =================
def _orders(errs):
    return [math.log(max(errs[i], 1e-300) / max(errs[i + 1], 1e-300), 2.0) for i in range(len(errs) - 1)]


def _lumped_errors(dts, theta, t0=400.0, t_end=2000.0):
    par = {"n_nodes": 2, "k": 1e4, "n_mod": 0, "c_fluid": 0.0, "t_space": 0.0}
    pm = _merge(par)
    exact = (t0 ** -3 + 3.0 * (2.0 * optics(pm)[1]) * SIGMA / pm["c_area"] * t_end) ** (-1.0 / 3.0)
    errs = []
    for dt in dts:
        pl = Plant(10.0, par, static_drive(300.0, 0.0, 0.0, 0.0))
        errs.append(abs(pl.simulate(np.array([300.0] + [t0] * 3), t_end, dt, theta, rannacher=0)["y"][1] - exact))
    return errs


def _linear_fin_errors(ns, qp=10.0, h=11.0, tinf=280.0):
    kd, L = PAR["k"] * PAR["delta"], PAR["L"]
    m = math.sqrt(h / kd)
    th0, thn = qp / (kd * m) / math.tanh(m * L), qp / (kd * m) / math.sinh(m * L)
    out = []
    for N in ns:
        pl = Plant(10.0, {"n_nodes": N, "n_mod": 0, "c_fluid": 0.0, "linear_h": h, "t_inf": tinf}, static_drive(300.0, 0.0, qp * (10.0 / L), 0.0))
        y = pl.steady()
        out.append((abs(y[1] - (tinf + th0)), abs(y[-1] - (tinf + thn))))
    return out


def _nonlinear_fin_errors(ns, t_root=330.0, t_tip=300.0):
    L, q = fin_first_integral(t_root, t_tip, PAR["k"] * PAR["delta"], 2.0 * optics(_merge())[1])
    out = []
    for N in ns:
        pl = Plant(10.0, {"L": L, "n_nodes": N, "n_mod": 0, "c_fluid": 0.0, "t_space": 0.0}, static_drive(300.0, 0.0, q * (10.0 / L), 0.0))
        y = pl.steady()
        out.append((abs(y[1] - t_root), abs(y[-1] - t_tip)))
    return out


def test_regressao_2b_2c_contra_saidas_locais_do_usuario():
    h = bus_heat(50000.0, 40000.0, eta_pol=0.90)
    assert abs(h["q_node"] - 51483.17) < 0.02, ("q_node", h["q_node"])
    assert abs((0.91 * 50000.0 - 1000.0) * 0.9655 * 0.90 - 38668.275) < 1e-3
    a = analytic_orbit(1105.39, T_ORB, 18000.0, eta_pol=0.90, eta_ch=1.0, eta_dis=1.0)
    assert abs(a["P_TI_sun_sust_W"] - 34070.12) < 0.05, ("2b P_TI_sol", a["P_TI_sun_sust_W"])
    s0 = size_bess(1105.39, T_ORB, 18000.0, eta_pol=0.90, eta_ch=1.0, eta_dis=1.0, fade=0.0)
    assert abs(s0["E_req_EOL_Wh"] - 27495.03) < 0.05, ("2b E_nom", s0["E_req_EOL_Wh"])
    s = size_bess(T_ECL, T_ORB, 18000.0, eta_pol=0.90)
    assert s["n_par"] == 21 and s["n_cells"] == 2268, ("2c n", s["n_par"], s["n_cells"])
    assert abs(s["mass_cells_kg"] - 260.82) < 0.01 and abs(s["E_inst_BOL_Wh"] - 36741.6) < 0.05 and abs(s["E_req_EOL_Wh"] - 28345.41) < 0.05, s
    assert abs(s["dod_eol"] - 0.2411) < 5e-5 and abs(s["P_TI_sun_sust_W"] - 33904.08) < 0.05, s


def test_modulo_calibrado_em_400_W():
    assert abs(ia_p(1410.0, 80.0)[0] - 400.0) < 1e-9
    assert ia_p(900.0, 90.0)[0] > ia_p(900.0, 70.0)[0] > 0.0 and abs(ia_f_for_power(ia_p(800.0, 70.0)[0], 70.0)[0] - 800.0) < 1e-6


def test_lumped_resfriamento_radiativo_analitico_ordem_temporal():
    dts = (100.0, 50.0, 25.0, 12.5)
    cn, be = _lumped_errors(dts, 0.5), _lumped_errors(dts, 1.0)
    oc, ob = _orders(cn), _orders(be)
    assert cn[-1] < 5e-3 and 1.8 <= oc[-1] <= 2.3, ("CN", cn, oc)
    assert 0.8 <= ob[-1] <= 1.3, ("BE", be, ob)


def test_aleta_linear_analitica_ordem_espacial():
    e = _linear_fin_errors((10, 20, 40, 80))
    o0, on = _orders([x[0] for x in e]), _orders([x[1] for x in e])
    assert e[-1][0] < 5e-4 and e[-1][1] < 5e-4 and o0[-1] >= 1.8 and on[-1] >= 1.8, (e, o0, on)


def test_aleta_radiante_nao_linear_vs_integral_primeira():
    e = _nonlinear_fin_errors((20, 40, 80, 160))
    o0, on = _orders([x[0] for x in e]), _orders([x[1] for x in e])
    assert e[-1][0] < 2e-3 and e[-1][1] < 2e-3 and o0[-1] >= 1.8 and on[-1] >= 1.8, (e, o0, on)


def test_estacionario_igual_ao_limite_do_transitorio():
    p = _merge()
    q = 2.0 * (0.36 * 1361.0 * 0.3 * F_VIEW_90 * 0.5 + 0.9 * 237.0 * F_VIEW_90)
    pl = Plant(60.0, None, static_drive(1300.0, q, 5500.0))
    ys = pl.steady()
    assert ys is not None
    r = pl.simulate(np.concatenate(([320.0], np.full(pl.N + 1, 300.0))), 40000.0, 100.0, 0.5)
    assert float(np.max(np.abs(r["y"] - ys))) < 1e-5, float(np.max(np.abs(r["y"] - ys)))
    assert abs(describe(pl, ys)["balanco_resid_W"]) < 1e-4, describe(pl, ys)["balanco_resid_W"]


def test_orbita_periodica_e_balanco_de_energia():
    p = _merge()
    drive = orbit_drive(p, 1300.0, 600.0, 5500.0, 1500.0)
    pl = Plant(60.0, p, drive)
    r = pl.simulate(np.concatenate(([320.0], np.full(pl.N + 1, 300.0))), 5 * p["t_orb"], 25.0, 0.5, orbit_breaks(p, 5))
    T = p["t_orb"]
    mx = [float(r["tj"][(r["t"] >= k * T) & (r["t"] <= (k + 1) * T)].max()) for k in range(5)]
    assert abs(mx[4] - mx[3]) < 1e-3, mx
    assert r["rel_resid"] < 1e-7, r["rel_resid"]
    assert float(r["y"].min()) >= p["t_space"] and bool(np.all(np.isfinite(r["y"])))


def test_monotonia_da_area_de_projeto():
    a1, a2, a3 = (design_area(x)["a_pl"] for x in (20000.0, 30000.0, 38000.0))
    assert a1 < a2 < a3, (a1, a2, a3)
    assert design_area(30000.0, {"coating": "Z93_BOL"})["a_pl"] < design_area(30000.0, {"coating": "AgFEP_BOL"})["a_pl"]


def test_projeto_respeita_limite_e_limites_fisicos():
    d = design_area(30000.0)
    assert abs(d["tj_C"] - 80.0) < 0.01 and d["t_tip_C"] <= d["t_root_C"] <= d["tj_C"] and 0.0 < d["eta_fin"] <= 1.0 and d["bi_face"] < 1e-3, d
    assert abs(d["balanco_resid_W"]) < 1e-3 * d["q_rej_W"], d


def test_dominio():
    ok = {"coating": "Z93_EOL30"}
    for bad in (lambda: Plant(-1.0, None, static_drive(1000.0, 0.0, 0.0)), lambda: Plant(10.0, {"n_nodes": 1}, static_drive(1000.0, 0.0, 0.0)),
                lambda: Plant(10.0, {"coating": "x"}, static_drive(1000.0, 0.0, 0.0)), lambda: Plant(10.0, {"t_ecl": 4000.0}, static_drive(1000.0, 0.0, 0.0)),
                lambda: Plant(10.0, {"foo": 1}, static_drive(1000.0, 0.0, 0.0)), lambda: design_area(-5.0), lambda: bus_heat(-1.0, 10.0),
                lambda: Plant(10.0, ok, None), lambda: Plant(10.0, ok, static_drive(1000.0, 0.0, 0.0)).simulate(np.full(42, 300.0), -1.0, 1.0)):
        try:
            bad()
        except ValueError:
            continue
        raise AssertionError("deveria recusar entrada invalida")


def test_reprodutibilidade():
    p = _merge()
    pl = Plant(60.0, p, orbit_drive(p, 1300.0, 600.0, 5500.0, 1500.0))
    y0 = np.concatenate(([320.0], np.full(pl.N + 1, 300.0)))
    a = pl.simulate(y0, 3000.0, 50.0, 0.5, orbit_breaks(p, 1))
    b = pl.simulate(y0, 3000.0, 50.0, 0.5, orbit_breaks(p, 1))
    assert bool(np.array_equal(a["tj"], b["tj"])) and bool(np.array_equal(a["y"], b["y"]))


def _cross_check():
    try:
        import modelo_ia as ia
        import modelo_barramento as mb
    except ImportError:
        return "SKIPPED (modelo_ia.py / modelo_barramento.py ausentes)"
    for f in (400.0, 900.0, 1410.0):
        for tj in (60.0, 80.0):
            assert abs(ia.p_mod(f, tj) - ia_p(f, tj)[0]) < 1e-9
    r = mb.bus_flow(50000.0, 40000.0, eta_pol=0.9)
    assert abs(r["q_node"] - bus_heat(50000.0, 40000.0, eta_pol=0.9)["q_node"]) < 1e-6
    return "PASS"


def test_conferencia_cruzada():
    assert _cross_check().split()[0] in ("PASS", "SKIPPED")


TESTS = [test_regressao_2b_2c_contra_saidas_locais_do_usuario, test_modulo_calibrado_em_400_W, test_lumped_resfriamento_radiativo_analitico_ordem_temporal,
         test_aleta_linear_analitica_ordem_espacial, test_aleta_radiante_nao_linear_vs_integral_primeira, test_estacionario_igual_ao_limite_do_transitorio,
         test_orbita_periodica_e_balanco_de_energia, test_monotonia_da_area_de_projeto, test_projeto_respeita_limite_e_limites_fisicos, test_dominio, test_reprodutibilidade]


def run_selftests():
    res = {}
    for f in TESTS:
        try:
            f()
            res[f.__name__] = "PASS"
        except AssertionError as e:
            res[f.__name__] = "FAIL: " + str(e)[:400]
        except Exception as e:
            res[f.__name__] = "ERROR " + type(e).__name__ + ": " + str(e)[:400]
    try:
        res["test_conferencia_cruzada"] = _cross_check()
    except Exception as e:
        res["test_conferencia_cruzada"] = "FAIL: " + type(e).__name__ + " " + str(e)[:300]
    return res


# ================= relatorios =================
def _r(o, nd=7):
    if isinstance(o, (float, np.floating)):
        o = float(o)
        return float("%.*g" % (nd, o)) if math.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, dict):
        return {k: _r(v, nd) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_r(v, nd) for v in o]
    if isinstance(o, np.ndarray):
        return _r(o.tolist(), nd)
    return o


def report_closure():
    out = {}
    for nome, ep in (("eta_pol_0.90_antes", 0.90), ("eta_pol_1.00_decisao", 1.00)):
        a, s = analytic_orbit(T_ECL, T_ORB, 18000.0, eta_pol=ep), size_bess(T_ECL, T_ORB, 18000.0, eta_pol=ep)
        s105 = size_bess(2023.3335273230243, T_ORB, 18000.0, eta_pol=ep)
        out[nome] = {"L_ecl_bus_18kW_W": a["L_ecl_bus_W"], "TI_sol_sustentavel_W": a["P_TI_sun_sust_W"], "TI_max_sol_sem_carregar_W": (0.91 * 50000.0 - 1000.0) * 0.9655 * ep,
                     "BESS_LTAN18": {k: s[k] for k in ("n_par", "n_cells", "mass_cells_kg", "E_inst_BOL_Wh", "E_req_EOL_Wh", "dod_eol")},
                     "BESS_LTAN10.5_18kW": {k: s105[k] for k in ("n_par", "n_cells", "mass_cells_kg")},
                     "q_node_TI40kW_PVon_W": bus_heat(50000.0, 40000.0, eta_pol=ep)["q_node"],
                     "q_gen_TI34kW_W": 34000.0 / (0.9655 * ep) + Q_OISL, "q_gen_TI40kW_W": 40000.0 / (0.9655 * ep) + Q_OISL}
    return out


def report_loads():
    rows = []
    for ti in (18000.0, 25000.0, 30000.0, 34000.0, 38000.0, 40000.0):
        s, e = bus_heat(50000.0, ti), bus_heat(0.0, ti)
        rows.append([ti, s["q_dab"], s["q_buck"], s["q_bd"], s["q_node"], e["q_bd"], e["q_node"]])
    return {"cargas_termicas_eta_pol1_[P_TI | Q_DAB | Q_Buck_sol | Q_BESS_sol | Q_no_sol | Q_BESS_ecl | Q_no_ecl]_W": rows}


def report_design():
    p_sun = analytic_orbit(T_ECL, T_ORB, 18000.0)["P_TI_sun_sust_W"]
    rows = []
    for ti in (18000.0, 25000.0, 30000.0, p_sun, 40000.0):
        d = design_area(ti)
        rows.append({k: d[k] for k in ("p_ti_W", "f_MHz", "a_pl_m2", "a_emit_2faces_m2", "t_root_C", "t_tip_C", "eta_fin", "q_rej_W", "q_rej_Wm2_planta", "t_eq_C", "R_mod_eq_K_W", "R_rad_eq_K_W_por_modulo", "bi_face", "balanco_resid_W")})
    return {"P_TI_sol_sustentavel_eta_pol1": p_sun, "projeto_Tj80_[coating Z93_EOL30]": rows}


def report_sens():
    p_sun = analytic_orbit(T_ECL, T_ORB, 18000.0)["P_TI_sun_sust_W"]
    base = design_area(p_sun)["a_pl"]
    rows = []
    for nome, par in (("Z93_BOL", {"coating": "Z93_BOL"}), ("AgFEP_BOL", {"coating": "AgFEP_BOL"}), ("LEGADO a0.14/e0.88", {"coating": "LEGADO"}),
                      ("r_int 0.04", {"r_int": 0.04}), ("r_int 0.16", {"r_int": 0.16}), ("k 100", {"k": 100.0}), ("k 400", {"k": 400.0}), ("delta 1mm", {"delta": 1e-3}),
                      ("delta 4mm", {"delta": 4e-3}), ("L 0.05", {"L": 0.05}), ("L 0.10", {"L": 0.10}), ("Q_IR 230", {"q_ir": 230.0}), ("Q_IR 244", {"q_ir": 244.0}),
                      ("albedo 0.28", {"albedo": 0.28}), ("albedo 0.32", {"albedo": 0.32}), ("beta 20 (albedo maior)", {"beta_deg": 20.0}), ("n_nodes 20", {"n_nodes": 20}),
                      ("kappa_dab (eta_dab 0.95)", {"kappa_dab": 1.0 / 0.95 - 1.0})):
        a = design_area(p_sun, par)["a_pl"]
        rows.append([nome, a, 100.0 * (a / base - 1.0)])
    return {"base_a_pl_m2": base, "sensibilidade_OAT_[caso | a_pl_m2 | dA_%]": rows}


def _orbit_run(a, fs, fe, qs, qe, n_orb, dt, y0=None, par=None, theta=0.5):
    p = _merge(par)
    pl = Plant(a, p, orbit_drive(p, fs, fe, qs, qe))
    y0 = np.concatenate(([320.0], np.full(pl.N + 1, 300.0))) if y0 is None else y0
    return pl, pl.simulate(y0, n_orb * p["t_orb"], dt, theta, orbit_breaks(p, n_orb))


def _orbit_stats(r, T, n_orb):
    m = r["t"] >= (n_orb - 1) * T
    return {"tj_max_C": r["tj"][m].max() - T0K, "tj_min_C": r["tj"][m].min() - T0K, "t_root_max_C": r["t0"][m].max() - T0K, "t_root_min_C": r["t0"][m].min() - T0K,
            "t_tip_max_C": r["tN"][m].max() - T0K, "tj_max_por_orbita_C": [float(r["tj"][(r["t"] >= k * T) & (r["t"] <= (k + 1) * T)].max() - T0K) for k in range(n_orb)],
            "rel_resid_energia": r["rel_resid"], "max_step_resid_J": r["max_step_resid_J"], "passos": r["steps"], "newton_its": r["newton_its"]}


def report_orbit(n_orb=6, dt=25.0):
    T = T_ORB
    p_sun = analytic_orbit(T_ECL, T_ORB, 18000.0)["P_TI_sun_sust_W"]
    des = design_area(p_sun)
    a, n = des["a_pl"], PAR["n_mod"]
    f_sun, f_ecl = ia_f_for_power(p_sun / n, 80.0)[0], ia_f_for_power(18000.0 / n, 80.0)[0]
    hs, he = bus_heat(50000.0, p_sun), bus_heat(0.0, 18000.0)
    qs, qe = hs["heat_aux"] + hs["q_buck"] + hs["q_bd"], he["heat_aux"] + he["q_buck"] + he["q_bd"]
    out = {"a_pl_projeto_m2": a, "f_sol_MHz": f_sun, "f_ecl_MHz": f_ecl, "q_rest_sol_W": qs, "q_rest_ecl_W": qe}
    pl, r = _orbit_run(a, f_sun, f_ecl, qs, qe, n_orb, dt)
    out["S1_energeticamente_viavel_TIsol_sust_TIecl_18kW"] = _orbit_stats(r, T, n_orb)
    h4, h4e = bus_heat(50000.0, 40000.0), bus_heat(0.0, 40000.0)
    pl0, r0 = _orbit_run(a, IA["f_nom"], IA["f_nom"], h4["heat_aux"] + h4["q_buck"] + h4["q_bd"], h4e["heat_aux"] + h4e["q_buck"] + h4e["q_bd"], n_orb, dt)
    out["S0_pico_40kW_continuo_(sem_limite_energetico)_mesma_area"] = _orbit_stats(r0, T, n_orb)
    out["dominio_PINN"] = {"x_m": [0.0, PAR["L"]], "t_s": [0.0, T], "T_K_min_max_S1": [float(r["y"].min()), float(max(r["t0"].max(), r["tj"].max()))]}
    return out


def report_conv():
    out = {}
    dts = (100.0, 50.0, 25.0, 12.5)
    cn, be = _lumped_errors(dts, 0.5), _lumped_errors(dts, 1.0)
    out["temporal_lumped_nao_linear_[dt | erro_CN_K | erro_BE_K]"] = [[d, c, b] for d, c, b in zip(dts, cn, be)]
    out["ordem_temporal_CN_BE"] = [_orders(cn), _orders(be)]
    ns = (10, 20, 40, 80, 160)
    el, en = _linear_fin_errors(ns), _nonlinear_fin_errors(ns)
    out["espacial_linear_[N | erro_raiz_K | erro_ponta_K]"] = [[n, x[0], x[1]] for n, x in zip(ns, el)]
    out["espacial_nao_linear_vs_integral_primeira_[N | erro_raiz_K | erro_ponta_K]"] = [[n, x[0], x[1]] for n, x in zip(ns, en)]
    out["ordem_espacial_linear"] = [_orders([x[0] for x in el]), _orders([x[1] for x in el])]
    out["ordem_espacial_nao_linear"] = [_orders([x[0] for x in en]), _orders([x[1] for x in en])]
    p_sun = analytic_orbit(T_ECL, T_ORB, 18000.0)["P_TI_sun_sust_W"]
    a = design_area(p_sun)["a_pl"]
    n = PAR["n_mod"]
    hs, he = bus_heat(50000.0, p_sun), bus_heat(0.0, 18000.0)
    args = (ia_f_for_power(p_sun / n, 80.0)[0], ia_f_for_power(18000.0 / n, 80.0)[0], hs["heat_aux"] + hs["q_buck"] + hs["q_bd"], he["heat_aux"] + he["q_buck"] + he["q_bd"])
    ref = _orbit_run(a, *args, 1, 3.125)[1]["edge_y"]
    rows = []
    for d in (50.0, 25.0, 12.5, 6.25):
        e = _orbit_run(a, *args, 1, d)[1]["edge_y"]
        rows.append([d, max(float(np.max(np.abs(x - y))) for x, y in zip(e, ref))])
    out["temporal_orbita_acoplada_[dt | max|dT|_vs_dt=3.125_nos_6_pontos_de_quebra_K]"] = rows
    out["ordem_temporal_orbita"] = _orders([x[1] for x in rows])
    return out


SECTION = "all"     # sandbox: "tests" | "closure" | "loads" | "design" | "sens" | "orbit" | "conv" ; uso normal: "all"
VIA_ERRO = False       # sandbox: devolve o resultado pela mensagem de erro (unica via que ficou visivel para o assistente); uso normal: False

def _brief(sec, pl):
    """Resumo curto (so sandbox) para caber na via de saida visivel."""
    g = lambda x: "%.4g" % x
    L = lambda v: "/".join("%.2f" % x for x in v)
    pre = lambda p: next(v for k, v in pl.items() if k.startswith(p))
    if sec == "closure":
        s = ""
        for k in ("eta_pol_0.90_antes", "eta_pol_1.00_decisao"):
            d = pl[k]
            b = d["BESS_LTAN18"]
            s += "%s|TIsol=%s TImax=%s Lecl=%s n=%d m=%s E=%s q40=%s n105=%d;" % (k[8:12], g(d["TI_sol_sustentavel_W"]), g(d["TI_max_sol_sem_carregar_W"]), g(d["L_ecl_bus_18kW_W"]),
                                                                                 b["n_cells"], g(b["mass_cells_kg"]), g(b["E_inst_BOL_Wh"]), g(d["q_node_TI40kW_PVon_W"]), d["BESS_LTAN10.5_18kW"]["n_cells"])
        return s
    if sec == "design":
        rows = pre("projeto_Tj80")
        s = "Psol=%s;" % g(pl["P_TI_sol_sustentavel_eta_pol1"])
        for i in (0, 3, 4):
            r = rows[i]
            s += "ti=%s f=%s A=%s Tr=%s Tt=%s eta=%s q=%s Rm=%s Rr=%s;" % (g(r["p_ti_W"]), g(r["f_MHz"]), g(r["a_pl_m2"]), g(r["t_root_C"]), g(r["t_tip_C"]), g(r["eta_fin"]),
                                                                          g(r["q_rej_Wm2_planta"]), g(r["R_mod_eq_K_W"]), g(r["R_rad_eq_K_W_por_modulo"]))
        return s
    if sec == "orbit":
        s = "A=%s fs=%s fe=%s;" % (g(pl["a_pl_projeto_m2"]), g(pl["f_sol_MHz"]), g(pl["f_ecl_MHz"]))
        for tag in ("S1", "S0"):
            d = pre(tag)
            s += "%s Tjmax=%s Tjmin=%s Trmax=%s Trmin=%s Ttip=%s res=%s orb=%s;" % (tag, g(d["tj_max_C"]), g(d["tj_min_C"]), g(d["t_root_max_C"]), g(d["t_root_min_C"]), g(d["t_tip_max_C"]),
                                                                                  g(d["rel_resid_energia"]), "/".join(g(x) for x in d["tj_max_por_orbita_C"]))
        return s
    if sec == "conv":
        ol, on = pre("espacial_linear"), pre("espacial_nao_linear")
        return "CN=%s BE=%s linR=%s linT=%s nlR=%s nlT=%s orb=%s e=%s tl=%s lin=%s nl=%s" % (
            L(pl["ordem_temporal_CN_BE"][0]), L(pl["ordem_temporal_CN_BE"][1]), L(pl["ordem_espacial_linear"][0]), L(pl["ordem_espacial_linear"][1]),
            L(pl["ordem_espacial_nao_linear"][0]), L(pl["ordem_espacial_nao_linear"][1]), L(pl["ordem_temporal_orbita"]),
            "/".join(g(r[1]) for r in pre("temporal_orbita")), "/".join(g(x) for x in pre("temporal_lumped")[-1]), "/".join(g(x) for x in ol[-1]), "/".join(g(x) for x in on[-1]))
    if sec == "sens":
        rows = pre("sensibilidade")
        return "base=%s;" % g(pl["base_a_pl_m2"]) + ";".join("%s=%+.1f" % (r[0][:8].replace(" ", ""), r[2]) for r in rows)
    return "?"


if __name__ == "__main__":
    secs = {"tests": run_selftests, "closure": report_closure, "loads": report_loads, "design": report_design, "sens": report_sens, "orbit": report_orbit, "conv": report_conv}
    names = list(secs) if SECTION == "all" else [SECTION]
    payload = {n: secs[n]() for n in names}
    fails = {k: v for k, v in payload.get("tests", {}).items() if not (v == "PASS" or v.startswith("SKIPPED"))}
    txt = json.dumps(_r(payload), ensure_ascii=False, separators=(",", ":"))
    if fails:
        raise RuntimeError("TESTES FALHARAM: " + json.dumps(fails, ensure_ascii=False))
    if VIA_ERRO:
        if SECTION == "tests":
            t = payload["tests"]
            raise RuntimeError("PLACAR %d/%d PASS; nao-PASS=%s" % (sum(1 for v in t.values() if v == "PASS"), len(t), [k + ":" + v[:60] for k, v in t.items() if v != "PASS"]))
        raise RuntimeError("B:" + _brief(SECTION, payload[SECTION]))
    print(txt)
    result = {"json": txt}

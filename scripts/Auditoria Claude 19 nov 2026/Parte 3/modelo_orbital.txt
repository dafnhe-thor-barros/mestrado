"""modelo_orbital.py -- Marco 1 (versao final, escopo enxuto). Reconstruido do zero.

Etiquetas: [FONTE] valor em fonte identificada | [FONTE-RESUMO] valor visto apenas em resumo de busca, texto integral NAO lido por mim
| [CALCULADO] | [HIPOTESE DECLARADA] | [PENDENTE] | [NAO VERIFICADO].
Escopo: Sol (USNO), SSO circular (J2 secular), LTAN->RAAN->beta, eclipse (cone + penumbra + Terra oblata WGS84), fator de visao da Terra,
fluxos ambientais no radiador. NAO modelados (limitacoes declaradas): J3-J6, arrasto, terceiro corpo, escurecimento de limbo,
J2 de curto periodo no eclipse operacional, TT/UT1. Nao validado experimentalmente (sem TVAC/orbita).
Testes: test_modelo_orbital.py (pytest).  Relatorio de relevancia: marco1_relatorio.py.
"""
import math
import json
from datetime import datetime, timezone, timedelta
import numpy as np

# ---------------------------------------------------------------- constantes
WGS84_A_M = 6378137.0            # [FONTE-RESUMO] NGA.STND.0036_1.0.0_WGS84 (2014-07-08), Tab. 3.1: semi-eixo maior [m]
WGS84_INV_F = 298.257223563      # [FONTE-RESUMO] idem, 1/f [-]
WGS84_GM = 3.986004418e14        # [FONTE-RESUMO] idem, GM [m3/s2]
WGS84_J2 = 1.082629821313e-3     # [NAO VERIFICADO] valor lembrado do WGS84; candidato alternativo 1.08262999e-3 (dif. relativa ~1.6e-7, i_SSO muda < 1e-5 deg)
WGS84_B_M = WGS84_A_M * (1.0 - 1.0 / WGS84_INV_F)   # [CALCULADO] semi-eixo menor [m]
R_SUN_M = 6.957e8                # [FONTE-RESUMO] IAU 2015 Res. B3, raio solar nominal [m]
AU_M = 149597870700.0            # [FONTE-RESUMO] IAU 2012 Res. B2 (exato) [m]
S_NOMINAL = 1361.0               # [FONTE-RESUMO] IAU 2015 B3 / Kopp & Lean 2011: 1360.8 +- 0.5 W/m2
SUN_MEAN_LON_RATE_DEG_DAY = 0.98564736   # [NAO VERIFICADO] coef. de memoria (USNO q = 280.459 + 0.98564736 D); resumo de busca deu 0.9856474 (Almanaque), coerente; pagina USNO nao vista
SUN_RATE_RAD_S = math.radians(SUN_MEAN_LON_RATE_DEG_DAY) / 86400.0   # [CALCULADO] taxa do Sol medio = alvo da precessao SSO
ALBEDO = 0.30                    # [FONTE-RESUMO] NASA SP-8067 (1971), valor medio anual 0.30 +- 0.02 (texto nao visto por mim)
Q_IR = 237.0                     # [FONTE-RESUMO] idem, emissao media anual 237 +- 7 W/m2 (legado usa 234, sem fonte)
EPS_RAD = 0.88                   # [PENDENTE] emissividade do radiador (legado, sem fonte)
ALPHA_RAD = 0.14                 # [PENDENTE] absortividade solar EOL (legado, sem fonte)
THR_ECL_S = 11.28                # [HIPOTESE DECLARADA] limiar de eclipse aprovado pelo usuario [s]
THR_Q_WM2 = 2.8                  # [HIPOTESE DECLARADA] limiar de fluxo absorvido aprovado pelo usuario [W/m2]


def rad(x):
    return math.radians(x)


def deg(x):
    return math.degrees(x)


# ---------------------------------------------------------------- tempo e Sol
def jd_from_datetime(dt):
    """JD = 2440587.5 + t_unix/86400 (UTC, sem segundos intercalares tratados). [CALCULADO]"""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return 2440587.5 + dt.timestamp() / 86400.0


def sun_usno(jd):
    """Algoritmo USNO de baixa precisao ('Computing Approximate Solar Coordinates', ~1 arcmin) [NAO VERIFICADO: constantes de memoria; coerentes com o resumo de busca (Almanaque: 280.460, 0.9856474, 357.528, 0.9856003) a ~1e-3 deg].
    D = JD-2451545; g = 357.529+0.98560028 D; q = 280.459+0.98564736 D; L = q+1.915 sin g+0.020 sin 2g;
    R = 1.00014-0.01671 cos g-0.00014 cos 2g [UA]; e = 23.439-0.00000036 D; EqT = q/15 - RA [h]."""
    D = jd - 2451545.0
    g = rad((357.529 + 0.98560028 * D) % 360.0)
    q = (280.459 + 0.98564736 * D) % 360.0
    L = rad((q + 1.915 * math.sin(g) + 0.020 * math.sin(2.0 * g)) % 360.0)
    R = 1.00014 - 0.01671 * math.cos(g) - 0.00014 * math.cos(2.0 * g)
    e = rad(23.439 - 0.00000036 * D)
    ra = math.atan2(math.cos(e) * math.sin(L), math.cos(L))
    dec = math.asin(math.sin(e) * math.sin(L))
    eqt_h = (q / 15.0 - deg(ra) / 15.0 + 12.0) % 24.0 - 12.0
    vec = np.array([math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)])
    return {"vec": vec, "ra_deg": deg(ra) % 360.0, "dec_deg": deg(dec), "R_au": R, "q_deg": q, "eqt_min": eqt_h * 60.0}


# ---------------------------------------------------------------- orbita
def kepler_period_s(alt_m):
    """T = 2 pi sqrt(a^3/GM), a = a_eq + h, orbita circular [CALCULADO]."""
    a = WGS84_A_M + alt_m
    return 2.0 * math.pi * math.sqrt(a ** 3 / WGS84_GM)


def secular_raan_rad_s(alt_m, inc_deg, ecc=0.0):
    """dOmega/dt = -(3/2) J2 (Re/p)^2 n cos i, p = a(1-e^2).
    Deducao em marco1_documento_tecnico.md (Lagrange + media orbital). [NAO VERIFICADO em fonte lida; verificado por propagador numerico]"""
    a = WGS84_A_M + alt_m
    n = math.sqrt(WGS84_GM / a ** 3)
    p = a * (1.0 - ecc ** 2)
    return -1.5 * WGS84_J2 * (WGS84_A_M / p) ** 2 * n * math.cos(rad(inc_deg))


def secular_raan_deg_day(alt_m, inc_deg, ecc=0.0):
    return deg(secular_raan_rad_s(alt_m, inc_deg, ecc)) * 86400.0


def sso_inclination_deg(alt_m, ecc=0.0):
    """Condicao heliossincrona: dOmega/dt = taxa do Sol medio  =>  cos i = -Omega_sun / (1.5 J2 (Re/p)^2 n). [CALCULADO]"""
    a = WGS84_A_M + alt_m
    n = math.sqrt(WGS84_GM / a ** 3)
    p = a * (1.0 - ecc ** 2)
    c = -SUN_RATE_RAD_S / (1.5 * WGS84_J2 * (WGS84_A_M / p) ** 2 * n)
    return deg(math.acos(c))


def orbit_frame(raan, inc):
    """Base (P, Q, h): P = direcao do no ascendente, Q = P rotacionado de 90 deg no plano, h = P x Q (normal). [rad]"""
    cO, sO, ci, si = math.cos(raan), math.sin(raan), math.cos(inc), math.sin(inc)
    P = np.array([cO, sO, 0.0])
    Q = np.array([-sO * ci, cO * ci, si])
    return P, Q, np.cross(P, Q)


def orbit_geometry(dt, ltan_h, alt_m, inc_deg=None, raan_deg=None, convention="mean", const_sun_dist=False,
                   shift_s=0.0, rsun_m=R_SUN_M, au_m=AU_M):
    """LTAN -> RAAN -> beta.  Omega = alpha_ref + 15 deg (LTAN - 12 h).
    convention='mean' (padrao, [HIPOTESE DECLARADA] por decisao do usuario): alpha_ref = q (longitude media do Sol, USNO).
    convention='apparent': alpha_ref = ascensao reta do Sol verdadeiro. beta = asin(h . s)."""
    s = sun_usno(jd_from_datetime(dt) + shift_s / 86400.0)
    inc = sso_inclination_deg(alt_m) if inc_deg is None else inc_deg
    ref = s["q_deg"] if convention == "mean" else s["ra_deg"]
    raan = (ref + 15.0 * (ltan_h - 12.0)) % 360.0 if raan_deg is None else raan_deg
    _, _, h = orbit_frame(rad(raan), rad(inc))
    beta = deg(math.asin(max(-1.0, min(1.0, float(h @ s["vec"])))))
    R_au = 1.0 if const_sun_dist else s["R_au"]
    rs = math.asin(rsun_m / (R_au * au_m))
    return {"beta_deg": beta, "raan_deg": raan, "inc_deg": inc, "sun": s["vec"], "rs_rad": rs,
            "T_s": kepler_period_s(alt_m), "q_deg": s["q_deg"], "eqt_min": s["eqt_min"], "R_au": R_au}


# ---------------------------------------------------------------- eclipse
def sunlit_fraction(g, rs, re, legacy=False, clip=True):
    """Fracao do disco solar visivel. g = separacao angular Sol-centro da Terra, rs/re = raios angulares (rad); re pode ser array.
    Area da interseccao de dois discos: A = rs^2 acos(a1) + re^2 acos(a2) - 1/2 sqrt((-d+rs+re)(d+rs-re)(d-rs+re)(d+rs+re)).
    legacy=True reproduz o 3o fator duplicado do legado (d-re+rs). Aproximacao de lente planar (discos pequenos)."""
    g = np.atleast_1d(np.asarray(g, float))
    re = np.broadcast_to(np.asarray(re, float), g.shape)
    out = np.ones_like(g)
    out[g <= re - rs] = 0.0
    pen = (g > re - rs) & (g < re + rs)
    if pen.any():
        s = g[pen]
        e = re[pen]
        a1 = np.clip((s * s + rs * rs - e * e) / (2.0 * s * rs), -1.0, 1.0)
        a2 = np.clip((s * s + e * e - rs * rs) / (2.0 * s * e), -1.0, 1.0)
        f3 = (s - e + rs) if legacy else (s - rs + e)
        rd = np.maximum(0.0, (-s + rs + e) * (s + rs - e) * f3 * (s + e + rs))
        ov = rs * rs * np.arccos(a1) + e * e * np.arccos(a2) - 0.5 * np.sqrt(rd)
        out[pen] = 1.0 - ov / (math.pi * rs * rs)
    return np.clip(out, 0.0, 1.0) if clip else out


def eclipse_cyl_beta_s(beta_deg, alt_m):
    """Sombra cilindrica, Terra esferica de raio a: t = (T/pi) acos( sqrt(1-(a/r)^2) / cos(beta) ). [CALCULADO]"""
    r = WGS84_A_M + alt_m
    x = math.sqrt(1.0 - (WGS84_A_M / r) ** 2) / math.cos(rad(beta_deg))
    return kepler_period_s(alt_m) / math.pi * math.acos(x) if x < 1.0 else 0.0


def eclipse_cone_beta_s(beta_deg, alt_m, rs, n=20000, legacy=False):
    """Eclipse equivalente por orbita, Terra esferica, cone com penumbra: integral de (1-f) dt; cos(gamma) = -cos(beta) cos(theta)."""
    r = WGS84_A_M + alt_m
    T = kepler_period_s(alt_m)
    re = math.asin(WGS84_A_M / r)
    th = (np.arange(n) + 0.5) * 2.0 * math.pi / n
    g = np.arccos(np.clip(-math.cos(rad(beta_deg)) * np.cos(th), -1.0, 1.0))
    return float((1.0 - sunlit_fraction(g, rs, re, legacy)).sum() * T / n)


def shadow_basis(s, a=WGS84_A_M, b=WGS84_B_M):
    """Silhueta do elipsoide projetada ao longo de s (raios paralelos): elipse de semi-eixos a (eixo e2 = s x e1) e
    B = sqrt(a^2 sz^2 + b^2 (1-sz^2)) (eixo e1 = z projetado), sz = s.z (funcao suporte do elipsoide). [CALCULADO]"""
    sz = float(s[2])
    e1 = np.array([0.0, 0.0, 1.0]) - sz * s
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(s, e1)
    return e1, e2, math.sqrt(a * a * sz * sz + b * b * (1.0 - sz * sz))


def eclipse_equivalent_s(raan_deg, inc_deg, sun_vec, alt_m, rs, a=WGS84_A_M, b=WGS84_B_M, n=20000,
                         mode="cone", legacy_lens=False):
    """Eclipse operacional por orbita [s]. Orbita circular de raio a_eq+h. Terra oblata via RAIO EFETIVO de sombra:
    para cada ponto, R_eff = 1/sqrt((c2/a)^2 + (c1/B)^2) (c1,c2 = componentes da direcao radial de p_perp na base da silhueta);
    r_e,eff = asin(R_eff/r) entra na formula de penumbra esferica. mode='cyl' (Sol pontual) e exato para o elipsoide;
    mode='cone' e [HIPOTESE DECLARADA] (penumbra do elipsoide aproximada pela da esfera de raio efetivo). b=a recupera a esfera."""
    r = WGS84_A_M + alt_m
    T = kepler_period_s(alt_m)
    P, Q, _ = orbit_frame(rad(raan_deg), rad(inc_deg))
    u = (np.arange(n) + 0.5) * 2.0 * math.pi / n
    p = r * (np.cos(u)[:, None] * P + np.sin(u)[:, None] * Q)
    s = np.asarray(sun_vec, float)
    s = s / np.linalg.norm(s)
    along = p @ s
    perp = p - along[:, None] * s
    rho = np.linalg.norm(perp, axis=1)
    uh = perp / np.maximum(rho, 1e-9)[:, None]
    e1, e2, B = shadow_basis(s, a, b)
    Reff = 1.0 / np.sqrt((uh @ e2 / a) ** 2 + (uh @ e1 / B) ** 2)
    if mode == "cyl":
        frac = 1.0 - ((along < 0.0) & (rho < Reff)).astype(float)
    else:
        re = np.arcsin(np.minimum(1.0, Reff / r))
        g = np.arccos(np.clip(-along / r, -1.0, 1.0))
        frac = sunlit_fraction(g, rs, re, legacy_lens)
    return float((1.0 - frac).sum() * T / n)


def eclipse_ray_s(raan_deg, inc_deg, sun_vec, alt_m, a=WGS84_A_M, b=WGS84_B_M, n=100000):
    """Teste INDEPENDENTE (raio x elipsoide, Sol pontual): z' = z a/b transforma o elipsoide em esfera de raio a; eclipse se o raio
    p + t s (t>0) tem raizes reais com B<0 (interseccao a frente)."""
    r = WGS84_A_M + alt_m
    T = kepler_period_s(alt_m)
    P, Q, _ = orbit_frame(rad(raan_deg), rad(inc_deg))
    u = (np.arange(n) + 0.5) * 2.0 * math.pi / n
    pos = r * (np.cos(u)[:, None] * P + np.sin(u)[:, None] * Q)
    k = a / b
    p = pos.copy()
    p[:, 2] *= k
    d = np.array(sun_vec, float)
    d[2] *= k
    A = float(d @ d)
    Bq = 2.0 * (p @ d)
    C = (p * p).sum(1) - a ** 2
    return float(((Bq < 0.0) & (Bq * Bq - 4.0 * A * C > 0.0)).sum() * T / n)


def scan_year(ltan_h, alt_m, year=2026, n=20000, convention="mean", oblate=True, mode="cone", step_days=1):
    """Varredura anual (00:00 UTC): beta e eclipse por dia."""
    d0 = datetime(year, 1, 1, tzinfo=timezone.utc)
    b = WGS84_B_M if oblate else WGS84_A_M
    rows = []
    for k in range(0, 365, step_days):
        g = orbit_geometry(d0 + timedelta(days=k), ltan_h, alt_m, convention=convention)
        t = eclipse_equivalent_s(g["raan_deg"], g["inc_deg"], g["sun"], alt_m, g["rs_rad"], b=b, n=n, mode=mode)
        rows.append((k, g["beta_deg"], t))
    arr = np.array(rows)
    return {"day": arr[:, 0], "beta_deg": arr[:, 1], "ecl_s": arr[:, 2]}


# ---------------------------------------------------------------- fator de visao
def view_factor_closed(alt_m, phi_deg, legacy=False):
    """F(h, phi) da Terra (esfera lambertiana uniforme) para elemento plano com normal a angulo phi do nadir.
    Zona 1 (phi <= 90-theta_m): cos(phi)/H^2. Zona 2: forma fechada com o 2o arccos de argumento -sqrt(H^2-1) cos(phi)/sin(phi)
    (legacy=True usa o sinal +, que e o defeito do legado). Fonte da forma fechada: [PENDENTE]; verificada contra a quadratura."""
    H = (WGS84_A_M + alt_m) / WGS84_A_M
    thm = math.asin(1.0 / H)
    root = math.sqrt(H * H - 1.0)
    phi = rad(phi_deg)

    def cl(v, lo, hi):
        return max(lo, min(hi, v))

    if phi <= math.pi / 2 - thm:
        return cl(math.cos(phi) / H ** 2, 0.0, 1.0)
    if phi <= math.pi / 2 + thm:
        sp = math.sin(phi)
        a1 = cl(root / (H * sp), -1.0, 1.0)
        a2 = cl((1.0 if legacy else -1.0) * root * math.cos(phi) / sp, -1.0, 1.0)
        rd = max(0.0, 1.0 - H * H * math.cos(phi) ** 2)
        f = (0.5 - math.asin(a1) / math.pi + math.cos(phi) / (math.pi * H * H) * math.acos(a2)
             - root * math.sqrt(rd) / (math.pi * H * H))
        return cl(f, 0.0, 1.0)
    return 0.0


def view_factor_quad(alt_m, phi_deg, na=1200, nz=1200):
    """F = (1/pi) Int_calota max(0, cos psi) dOmega (definicao; quadratura por pontos medios). [CALCULADO]"""
    H = (WGS84_A_M + alt_m) / WGS84_A_M
    thm = math.asin(1.0 / H)
    phi = rad(phi_deg)
    al = (np.arange(na) + 0.5) * thm / na
    ze = (np.arange(nz) + 0.5) * 2.0 * math.pi / nz
    A, Z = np.meshgrid(al, ze, indexing="ij")
    c = np.sin(A) * np.cos(Z) * math.sin(phi) + np.cos(A) * math.cos(phi)
    return float((np.maximum(c, 0.0) * np.sin(A)).sum() * (thm / na) * (2.0 * math.pi / nz) / math.pi)


# ---------------------------------------------------------------- radiador
def radiator_env(theta, beta_deg, alt_m, law="s_cross_r", legacy_vf=False, albedo=ALBEDO, q_ir=Q_IR,
                 eps=EPS_RAD, alpha=ALPHA_RAD, G=S_NOMINAL):
    """Fluxos absorvidos por UMA face de normal n [W/m2]. Referencial do legado: r = (cos th, sin th, 0), s = (cos b, 0, sin b), h = z.
    law='s_cross_r' (opcao 1, PADRAO aprovado): n = s x r / |s x r|  => phi = 90 deg, sem Sol direto.
    law='inertial_fixed' (opcao 2): n = componente de h perpendicular a s, fixa no inercial.  law='legacy': n = -s.
    Albedo avaliado no ponto subsatelite: q_alb = alpha G a F max(0, r.s)  [HIPOTESE DECLARADA, herdada do legado: 1a aproximacao].
    Dupla face: a face oposta tem normal -n (campo F_back)."""
    b = rad(beta_deg)
    r = np.array([math.cos(theta), math.sin(theta), 0.0])
    s = np.array([math.cos(b), 0.0, math.sin(b)])
    h = np.array([0.0, 0.0, 1.0])
    if law == "s_cross_r":
        n = np.cross(s, r)
        nn = float(np.linalg.norm(n))
        n = h if nn < 1e-9 else n / nn
    elif law == "inertial_fixed":
        n = h - float(h @ s) * s
        n = n / np.linalg.norm(n)
    else:
        n = -s
    phi = deg(math.acos(max(-1.0, min(1.0, float(-n @ r)))))
    F = view_factor_closed(alt_m, phi, legacy=legacy_vf)
    Fb = view_factor_closed(alt_m, 180.0 - phi, legacy=legacy_vf)
    mu = max(0.0, float(r @ s))
    return {"phi_deg": phi, "F": F, "F_back": Fb, "q_ir": eps * q_ir * F, "q_alb": alpha * G * albedo * F * mu,
            "q_sol": alpha * G * max(0.0, float(n @ s)), "n": n}


# ---------------------------------------------------------------- propagador J2 (apenas verificacao)
def _acc_j2(t, y):
    r = y[:3]
    rn = float(np.linalg.norm(r))
    k = 1.5 * WGS84_J2 * (WGS84_A_M / rn) ** 2
    z2 = (r[2] / rn) ** 2
    f = -WGS84_GM / rn ** 3
    return np.concatenate([y[3:], [f * r[0] * (1 - k * (5 * z2 - 1)), f * r[1] * (1 - k * (5 * z2 - 1)),
                                    f * r[2] * (1 - k * (5 * z2 - 3))]])


def propagate_j2(raan_deg, inc_deg, alt_m, t_end_s, step_s):
    """Integra a aceleracao J2 (forma classica, [NAO VERIFICADO em fonte]) com DOP853 (scipy). Elementos iniciais OSCULADORES."""
    from scipy.integrate import solve_ivp
    r0 = WGS84_A_M + alt_m
    P, Q, _ = orbit_frame(rad(raan_deg), rad(inc_deg))
    y0 = np.concatenate([r0 * P, math.sqrt(WGS84_GM / r0) * Q])
    ts = np.arange(0.0, t_end_s + step_s, step_s)
    sol = solve_ivp(_acc_j2, (0.0, ts[-1]), y0, method="DOP853", t_eval=ts, rtol=1e-10, atol=1e-4)
    return ts, sol.y[:3].T, sol.y[3:].T


def j2_drift_deg_day(ts, R, V):
    hv = np.cross(R, V)
    om = np.unwrap(np.arctan2(hv[:, 0], -hv[:, 1]))
    return deg(np.polyfit(ts, om, 1)[0]) * 86400.0


# ---------------------------------------------------------------- resumo de projeto
def design_summary(dt, ltan_h, alt_m, n=40000):
    g = orbit_geometry(dt, ltan_h, alt_m)
    args = (g["raan_deg"], g["inc_deg"], g["sun"], alt_m, g["rs_rad"])
    out = {"data_utc": dt.isoformat(), "ltan_h": ltan_h, "alt_km": alt_m / 1e3, "beta_deg": g["beta_deg"], "inc_deg": g["inc_deg"],
           "raan_deg": g["raan_deg"], "T_s": g["T_s"],
           "ecl_elipsoide_cone_s": eclipse_equivalent_s(*args, b=WGS84_B_M, n=n),
           "ecl_esfera_cone_s": eclipse_equivalent_s(*args, b=WGS84_A_M, n=n),
           "ecl_esfera_cilindro_s": eclipse_cyl_beta_s(g["beta_deg"], alt_m)}
    e = radiator_env(0.0, g["beta_deg"], alt_m)
    out["radiador_opcao1"] = {"phi_deg": e["phi_deg"], "F": e["F"], "q_ir_Wm2": e["q_ir"], "q_sol_Wm2": e["q_sol"]}
    return out


if __name__ == "__main__":
    _out = design_summary(datetime(2026, 6, 1, tzinfo=timezone.utc), 10.5, 700e3, n=20000)
    print(json.dumps(_out, indent=1, default=float))
    result = {"json": json.dumps(_out, default=float)}

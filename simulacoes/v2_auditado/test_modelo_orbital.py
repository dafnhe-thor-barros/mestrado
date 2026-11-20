"""Testes pytest do modelo_orbital (Marco 1). Rodar:  pytest -q test_modelo_orbital.py
Verificacao = codigo implementa a equacao (casos analiticos, convergencia, implementacoes independentes). Nao e validacao experimental."""
import math
from datetime import datetime, timezone, timedelta
import numpy as np
import pytest
import modelo_orbital as mo

ALT = 700e3


def _utc(y, m, d, h=0):
    return datetime(y, m, d, h, tzinfo=timezone.utc)


# ---------- tempo, Sol, constantes ----------
def test_jd_j2000():
    assert abs(mo.jd_from_datetime(_utc(2000, 1, 1, 12)) - 2451545.0) < 1e-9


def test_b_wgs84_derivado_e_literatura():
    assert abs(mo.WGS84_B_M - mo.WGS84_A_M * (1 - 1 / mo.WGS84_INV_F)) < 1e-9
    assert abs(mo.WGS84_B_M - 6356752.3142) < 1e-3   # valor de literatura lembrado [NAO VERIFICADO]


def test_periodo_kepler():
    a = mo.WGS84_A_M + ALT
    n = math.sqrt(mo.WGS84_GM / a ** 3)
    assert abs(mo.kepler_period_s(ALT) - 2 * math.pi / n) < 1e-9
    assert 5925.0 < mo.kepler_period_s(ALT) < 5928.0


def test_sol_equinocio_solsticio_equacao_do_tempo():
    d19 = mo.sun_usno(mo.jd_from_datetime(_utc(2026, 3, 19, 12)))["dec_deg"]
    d22 = mo.sun_usno(mo.jd_from_datetime(_utc(2026, 3, 22, 12)))["dec_deg"]
    ds = mo.sun_usno(mo.jd_from_datetime(_utc(2026, 6, 21, 12)))["dec_deg"]
    eq = [mo.sun_usno(mo.jd_from_datetime(_utc(2026, 1, 1) + timedelta(days=k)))["eqt_min"] for k in range(365)]
    assert d19 < 0.0 < d22
    assert 23.35 < ds < 23.50
    assert 16.0 < max(eq) < 17.0 and -15.0 < min(eq) < -13.5


def test_distancia_solar_perielio_afelio():
    r_jan = mo.sun_usno(mo.jd_from_datetime(_utc(2026, 1, 4)))["R_au"]
    r_jul = mo.sun_usno(mo.jd_from_datetime(_utc(2026, 7, 6)))["R_au"]
    assert r_jan < 0.985 and r_jul > 1.015


# ---------- SSO por J2 ----------
def test_sso_inclinacao_consistente_com_taxa_do_sol():
    i = mo.sso_inclination_deg(ALT)
    assert 98.0 < i < 98.4                                   # decisao do usuario: i ~ 98.19 deg
    assert abs(mo.secular_raan_rad_s(ALT, i) - mo.SUN_RATE_RAD_S) < 1e-15


def test_sso_inclinacao_cresce_com_altitude():
    assert mo.sso_inclination_deg(500e3) < mo.sso_inclination_deg(700e3) < mo.sso_inclination_deg(1000e3)


def test_j2_taxa_secular_vs_propagador():
    pytest.importorskip("scipy")
    g = mo.orbit_geometry(_utc(2026, 3, 1), 18.0, ALT)
    ts, R, V = mo.propagate_j2(g["raan_deg"], g["inc_deg"], ALT, 86400.0, 20.0)
    num = mo.j2_drift_deg_day(ts, R, V)
    ana = mo.secular_raan_deg_day(ALT, g["inc_deg"])
    # tolerancia 2%: elementos iniciais osculadores (semi-eixo osculador != medio) -> diferenca esperada ~0.5% [hipotese]
    assert abs(num - ana) / abs(ana) < 0.02, (num, ana)


# ---------- LTAN, beta ----------
def test_beta_por_ltan_em_2026():
    d0 = _utc(2026, 1, 1)
    res = {}
    for L in (18.0, 6.0, 10.5):
        b = [abs(mo.orbit_geometry(d0 + timedelta(days=k), L, ALT)["beta_deg"]) for k in range(365)]
        res[L] = (min(b), max(b))
    assert 57.0 < res[18.0][0] < 60.0 and res[18.0][1] > 85.0
    assert 57.0 < res[6.0][0] < 60.0 and res[6.0][1] > 85.0
    assert res[10.5][0] > 10.0 and res[10.5][1] < 30.0


def test_ltan_convencao_media_vs_aparente():
    d = _utc(2026, 2, 11)    # proximo do minimo da equacao do tempo
    gm = mo.orbit_geometry(d, 18.0, ALT, convention="mean")
    ga = mo.orbit_geometry(d, 18.0, ALT, convention="apparent")
    dra = abs((gm["raan_deg"] - ga["raan_deg"] + 180.0) % 360.0 - 180.0)
    assert abs(dra - abs(gm["eqt_min"]) * 0.25) < 0.05      # 1 min de tempo = 0.25 deg


def test_beta_normal_do_plano_ortogonal():
    P, Q, h = mo.orbit_frame(math.radians(70.0), math.radians(98.19))
    assert abs(P @ Q) < 1e-12 and abs(np.linalg.norm(h) - 1.0) < 1e-12


# ---------- eclipse ----------
def test_lente_vs_integracao_esferica_e_defeito_do_legado():
    rs = math.asin(mo.R_SUN_M / mo.AU_M)
    re = math.asin(mo.WGS84_A_M / (mo.WGS84_A_M + ALT))
    xs = np.linspace(-rs, rs, 301)
    X, Y = np.meshgrid(xs, xs)
    r = np.hypot(X, Y)
    m = r <= rs
    ph = np.arctan2(Y, X)
    vx, vz = np.sin(r) * np.cos(ph), np.cos(r)
    worst = 0.0
    for k in (-0.9, -0.5, 0.0, 0.5, 0.9):
        g = re + k * rs
        ex = float(((np.arccos(np.clip(vx * math.sin(g) + vz * math.cos(g), -1, 1)) > re) & m).sum() / m.sum())
        worst = max(worst, abs(float(mo.sunlit_fraction(g, rs, re)[0]) - ex))
    assert worst < 4e-3
    assert float(mo.sunlit_fraction(re, rs, re, legacy=True, clip=False)[0]) < -10.0


def test_cone_tende_ao_cilindro_quando_rsol_tende_a_zero():
    assert abs(mo.eclipse_cone_beta_s(30.0, ALT, 1e-7, n=200000) - mo.eclipse_cyl_beta_s(30.0, ALT)) < 0.5


def test_operacional_esfera_cilindrica_igual_analitico():
    g = mo.orbit_geometry(_utc(2026, 6, 1), 10.5, ALT)
    t_op = mo.eclipse_equivalent_s(g["raan_deg"], g["inc_deg"], g["sun"], ALT, 1e-9, b=mo.WGS84_A_M, n=200000, mode="cyl")
    assert t_op > 1000.0
    assert abs(t_op - mo.eclipse_cyl_beta_s(g["beta_deg"], ALT)) < 0.3


@pytest.mark.parametrize("ltan,data", [(10.5, (2026, 6, 1)), (14.0, (2026, 9, 1)), (18.0, (2026, 11, 20))])
def test_elipsoide_operacional_vs_raycast_independente(ltan, data):
    g = mo.orbit_geometry(_utc(*data), ltan, ALT)
    t_op = mo.eclipse_equivalent_s(g["raan_deg"], g["inc_deg"], g["sun"], ALT, 1e-9, n=200000, mode="cyl")
    t_ray = mo.eclipse_ray_s(g["raan_deg"], g["inc_deg"], g["sun"], ALT, n=200000)
    assert abs(t_op - t_ray) < 0.3, (t_op, t_ray)


def test_convergencia_em_n():
    rs = math.asin(mo.R_SUN_M / mo.AU_M)
    ref = mo.eclipse_cone_beta_s(63.8, ALT, rs, n=640000)
    err = {n: abs(mo.eclipse_cone_beta_s(63.8, ALT, rs, n=n) - ref) for n in (5000, 10000, 20000, 40000)}
    assert err[40000] < 0.5
    assert err[40000] <= err[5000] + 1e-9


def test_eclipse_monotono_em_beta_e_limites():
    rs = math.asin(mo.R_SUN_M / mo.AU_M)
    v = [mo.eclipse_cone_beta_s(b, ALT, rs, n=40000) for b in np.arange(0.0, 64.6, 0.5)]
    assert all(v[i + 1] <= v[i] + 0.05 for i in range(len(v) - 1))
    assert v[0] < mo.kepler_period_s(ALT) / 2.0 and min(v) >= 0.0


def test_fracao_solar_em_0_1():
    rs = math.asin(mo.R_SUN_M / mo.AU_M)
    g = np.linspace(0.0, 3.1, 4000)
    f = mo.sunlit_fraction(g, rs, math.asin(mo.WGS84_A_M / (mo.WGS84_A_M + ALT)))
    assert f.min() >= 0.0 and f.max() <= 1.0


# ---------- fator de visao e radiador ----------
def test_fator_de_visao_vs_quadratura():
    worst = max(abs(mo.view_factor_closed(ALT, p) - mo.view_factor_quad(ALT, p)) for p in (10, 26, 30, 45, 60, 90, 120, 150, 154))
    assert worst < 5e-4, worst
    assert abs(mo.view_factor_closed(ALT, 30, legacy=True) - mo.view_factor_quad(ALT, 30)) > 0.3


def test_fator_de_visao_limites_e_continuidade():
    H = (mo.WGS84_A_M + ALT) / mo.WGS84_A_M
    thm = math.degrees(math.asin(1 / H))
    assert abs(mo.view_factor_closed(ALT, 0.0) - 1 / H ** 2) < 1e-12
    assert mo.view_factor_closed(ALT, 170.0) == 0.0
    assert abs(mo.view_factor_closed(ALT, 90 - thm - 1e-6) - mo.view_factor_closed(ALT, 90 - thm + 1e-6)) < 1e-4
    assert abs(mo.view_factor_closed(ALT, 90 + thm - 1e-6) - mo.view_factor_closed(ALT, 90 + thm + 1e-6)) < 1e-4
    assert all(0.0 <= mo.view_factor_closed(ALT, p) <= 1.0 for p in np.linspace(0, 180, 361))


def test_radiador_opcao1_phi90_sem_sol_direto_faces_iguais():
    for th in np.linspace(0.0, 2 * math.pi, 37):
        e = mo.radiator_env(th, 58.4, ALT, law="s_cross_r")
        assert abs(e["phi_deg"] - 90.0) < 1e-6 and abs(e["q_sol"]) < 1e-9
        assert abs(e["F"] - e["F_back"]) < 1e-9


def test_reprodutibilidade():
    rs = math.asin(mo.R_SUN_M / mo.AU_M)
    assert mo.eclipse_cone_beta_s(60.0, ALT, rs) == mo.eclipse_cone_beta_s(60.0, ALT, rs)

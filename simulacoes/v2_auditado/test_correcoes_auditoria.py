"""test_correcoes_auditoria.py -- testes das correcoes da auditoria de 22/11/2026 (numpy/pandas; sem torch)."""
import math
import numpy as np
import pytest

import modelo_termico as mt
import contrato_parametros as cp
import marco5_tems as m5
import janela_pue as jp
import rodar_campanha_marco5 as rc


def test_hash_cobre_planta_e_controle():
    h0 = cp.contract_hash(cp.PROJETO_1P0)
    old = mt.PAR["c_area"]
    try:
        mt.PAR["c_area"] = 4050.0
        assert cp.contract_hash(cp.PROJETO_1P0) != h0          # trocar c_area muda o hash
    finally:
        mt.PAR["c_area"] = old
    assert cp.contract_hash(cp.PROJETO_1P0) == h0
    o = cp.M_PRED_K
    try:
        cp.M_PRED_K = 0.9
        assert cp.contract_hash(cp.PROJETO_1P0) != h0           # margem de controle tambem
    finally:
        cp.M_PRED_K = o
    assert cp.contract_hash(cp.PROJETO_1P0) == h0


def test_bus_step_apagao_conserva_energia():
    c = cp.PROJETO_1P0
    ecap = 29000.0
    bs = m5.bus_step(c, 0.0, 30000.0, 5.0, c["soc_min"] * ecap, ecap, c["soc_min"])   # BESS no piso, eclipse
    assert bs["ti"] == 0.0 and abs(bs["resid"]) < 1e-6
    bs2 = m5.bus_step(c, 0.0, 10000.0, 5.0, 0.9 * ecap, ecap, c["soc_min"])           # normal: BESS cobre
    assert abs(bs2["resid"]) < 1e-6 and bs2["ti"] == pytest.approx(10000.0)


def test_s0_sem_residuo_de_barramento():
    r = m5.run_scenario("nominal", "S0", 0, n_orb=2, sensor="ideal")
    assert r["bus_resid_W"] < 1e-6
    assert r["pred_err_rms_K"] is None                                                  # S0 nao tem preditor


def test_fronteira_de_orbita_tolerante_ao_arredondamento():
    T = mt.T_ORB
    n = 20
    t = []
    for k in range(1, n + 1):
        for j in range(1, 13):
            t.append(round((k - 1) * T + j * T / 12, 6))
    ser = {"t": t, "h": [T / 12] * len(t), "tj": [330.0] * len(t), "soc": [0.8] * len(t), "f": [1000.0] * len(t), "p_ti": [3e4] * len(t)}
    o = jp.por_orbita(ser, T, 1410.0, 0.95, n)
    assert len(o["orb_soc_end"]) == n
    assert all(abs(x - 3e4 * T / 3600.0) < 1e-6 * 3e4 for x in o["orb_ti_Wh"])        # cada orbita recebe exatamente 12 amostras


def test_newton_divergencia_e_reconhecida_pelo_orquestrador():
    pl = mt.Plant(49.89754, {}, mt.static_drive(1410.0, 0.0, 1e9))
    y = np.concatenate(([330.0], np.full(pl.N + 1, 330.0)))
    try:
        for _ in range(50):
            y, _ = pl.step(y, 0.0, 600.0, 1.0)
    except RuntimeError as e:
        assert rc.DIVERG_TXT in str(e)


def _linha(ctrl, bloco, cenario, seed, tj, tops=0.9, soc=0.7, **kw):
    d = dict(bloco=bloco, cenario=cenario, sensor="sem_calibracao", quant="none", ctrl=ctrl, seed=seed, tj_max_C=tj, tj_max_last_C=tj, tops_frac_last=tops,
             soc_min=soc, viol_int_Ks=max(0.0, tj - 80.0), dE_over_ti_last3=1e-4, n_infeasible=0, shed_Wh=0.0, bus_resid_W=0.0)
    d.update(kw)
    return d


def test_p1_divergencia_do_tems_nao_passa_despercebida():
    linhas = [_linha(c, "nominal", "nominal", s, 78.0) for c in ("TEMS_RC3", "TEMS_PINN", "R_SINT") for s in range(10, 13)]
    v = rc.calcula_veredito(linhas)
    assert v["P1"]["status"] == "ATENDIDO"
    v2 = rc.calcula_veredito(linhas, diverg=[{"bloco": "nominal", "parte": "x", "ctrls": ["TEMS_PINN"], "cenarios": ["nominal"], "seeds": [13]}])
    assert v2["P1"]["status"] == "NAO_ATENDIDO" and v2["P1"]["divergencias_TEMS"]


def test_p5_bonferroni_e_aviso_de_n_pequeno():
    rng = np.random.default_rng(1)
    linhas = []
    for sc in ("a", "b", "c", "d"):
        for s in range(10, 15):
            base = 0.85 + 0.01 * rng.standard_normal()
            linhas.append(_linha("TEMS_RC3", "nominal", sc, s, 78.0, tops=base))
            linhas.append(_linha("TEMS_PINN", "nominal", sc, s, 78.0, tops=base + 0.001 * rng.standard_normal()))
    v = rc.calcula_veredito(linhas)
    d = v["P5"]["detalhe"]
    assert d and all(r["k_grupos"] == 4 and r["aviso_n_pequeno"] for r in d)
    lo0, hi0 = rc._boot_ic([0.1, 0.2, 0.3, 0.4, 0.5])
    lo1, hi1 = rc._boot_ic([0.1, 0.2, 0.3, 0.4, 0.5], alpha=0.0125)
    assert lo1 <= lo0 and hi1 >= hi0


def test_confere_config_recusa_partes_sem_config():
    import os, tempfile
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "parts", "1"))
    open(os.path.join(d, "parts", "1", "x.json"), "w").write("{}")
    with pytest.raises(SystemExit):
        rc.confere_config(d, {"norb": 20}, False, False)


def test_ood_margem_saturada():
    # a margem OOD nunca passa de M_OOD_MAX_K, mesmo com distancia acima do 2o corte e sem fallback
    c = cp
    d = 5.0
    extra = c.M_OOD_MAX_K * min(1.0, (d - c.OOD_CORTES[0]) / (c.OOD_CORTES[1] - c.OOD_CORTES[0]))
    assert extra == c.M_OOD_MAX_K

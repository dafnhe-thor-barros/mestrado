"""test_marco5.py -- testes do Marco 5 (contrato + malha fechada TEMS, DVFS-only).
Uso:  python -m pytest test_marco5.py -q        (ou: python test_marco5.py)
Os 12 primeiros testes usam so numpy/scipy (+ pytest). Os 2 ultimos (API da PINN e hash dos checkpoints) exigem torch + modelo_pinn.py / os .pt e sao PULADOS se ausentes.
Ordem de grandeza: ~30 s no total (a malha fechada com TEMS e o que pesa).
"""
import os

import numpy as np
import pytest

from modelo_termico import bus_heat, ia_p, _merge
from contrato_parametros import (PROJETO_1P0, LEGADO_2B_0P90, CHECKPOINTS_SHA256, contract_hash, verify_checkpoint, sensor_cfg, sha256_file)
from marco5_tems import (A_PL_PROJETO, RC3, SolverPred, Sensors, RC3 as _RC3, CONTROLADORES, battery_for, sun_cap, bus_step, env_q_abs, plant_par,
                         q_rest_model, make_pred, run_closed_loop, scenario_cfg)


import json
from marco5_tems import run_scenario, valida_sintonia, carrega_sintonia, chave_sintonia, ultima_orbita
from sintonia_reativo import escolhe


def _short(ctrl, scen="nominal", seed=0, sensor="ideal", t_end=None, n_orb=1, m_total=None):
    cfg = scenario_cfg(scen, seed, PROJETO_1P0, n_orb, sensor)
    if t_end:
        cfg["t_end"] = float(t_end)
    kind, pk = CONTROLADORES[ctrl]
    cfg["ctrl"], cfg["m_total"] = kind, m_total
    if pk:
        cfg["pred"] = make_pred(pk, PROJETO_1P0, cfg.get("plant_par") if pk == "SOLVER" else None, A_PL_PROJETO)
    return run_closed_loop(cfg)


def _sem_tempo(r):
    return {k: v for k, v in r.items() if not k.startswith("decision_ms")}


def test_contrato_reproduz_regressoes_confirmadas():
    h = bus_heat(50000.0, 40000.0, eta_pol=0.90)
    assert abs(h["q_node"] - 51483.17) < 0.02, h["q_node"]
    s, _ = battery_for(LEGADO_2B_0P90)
    assert s["n_cells"] == 2268 and abs(s["mass_cells_kg"] - 260.82) < 0.01, s
    assert abs(sun_cap(LEGADO_2B_0P90) - 33904.08) < 0.05, sun_cap(LEGADO_2B_0P90)
    s1, _ = battery_for(PROJETO_1P0)
    assert s1["n_par"] == 19 and s1["n_cells"] == 2052, s1
    assert abs(sun_cap(PROJETO_1P0) - 38176.18) < 0.05, sun_cap(PROJETO_1P0)


def test_hash_do_contrato_muda_com_parametro():
    h1 = contract_hash(PROJETO_1P0)
    assert h1 == contract_hash(dict(PROJETO_1P0)) and contract_hash(dict(PROJETO_1P0, eta_buck=0.90)) != h1


def test_checkpoint_com_hash_divergente_aborta(tmp_path):
    p = tmp_path / "v2c.pt"
    p.write_bytes(b"nao e o checkpoint")
    with pytest.raises(RuntimeError):
        verify_checkpoint(str(p))
    with pytest.raises(ValueError):
        verify_checkpoint(str(tmp_path / "xx.pt"))


def test_bus_step_fecha_balanco_e_respeita_piso():
    c = PROJETO_1P0
    _, ecap = battery_for(c)
    rng = np.random.default_rng(3)
    for _ in range(300):
        pv = float(rng.choice([0.0, 50000.0]))
        ti, E = float(rng.uniform(0.0, 42000.0)), float(rng.uniform(0.72, 0.95)) * ecap
        r = bus_step(c, pv, ti, 5.0, E, ecap, 0.70)
        assert abs(r["resid"]) < 1e-6 and r["ti"] <= ti + 1e-9 and r["curtail"] >= -1e-9, r
        assert 0.70 * ecap - 1e-6 <= r["E"] <= 0.95 * ecap + 1e-6, r


def test_sensor_vies_inseguro_e_quantizacao():
    s = Sensors(sensor_cfg("sem_calibracao", sigma_tj=0.0, sigma_rtd=0.0), 0)
    tj, tr = s.read(350.0, 320.0)
    assert abs(tj - 348.5) < 1e-12 and abs(tr - 319.75) < 1e-12
    q = Sensors(sensor_cfg("ideal", quant_mode="floor"), 0)
    assert abs(q.read(273.15 + 79.9, 300.0)[0] - (273.15 + 79.0)) < 1e-9


def test_tj_cresce_com_f_nos_preditores():
    c = PROJETO_1P0
    _, ecap = battery_for(c)
    p = _merge(plant_par(c))
    for pr in (RC3(plant_par(c), A_PL_PROJETO), SolverPred(plant_par(c), A_PL_PROJETO)):
        st = pr.init_state(330.0, 310.0, 7.5, None)
        fim = []
        for f in (500.0, 900.0, 1300.0, 1410.0):
            qr = q_rest_model(c, ecap, True, c["n_mod"] * ia_p(f, 57.0)[0])
            fim.append(float(pr.window(st, f, env_q_abs(1000.0, p), qr, 300.0, 15.0)["tj"][-1]))   # T_j no FIM da janela
        assert all(fim[i] < fim[i + 1] for i in range(3)), (pr.name, fim)


def test_rc3_concorda_com_solver_em_janela():
    c = PROJETO_1P0
    _, ecap = battery_for(c)
    p = _merge(plant_par(c))
    sol, eu, ra = SolverPred(plant_par(c), A_PL_PROJETO), RC3(plant_par(c), A_PL_PROJETO), RC3(plant_par(c), A_PL_PROJETO, "radau")
    worst = 0.0
    for f in (600.0, 1400.0):
        qr = q_rest_model(c, ecap, True, c["n_mod"] * ia_p(f, 60.0)[0])
        ref = sol.window(sol.init_state(340.0, 315.0, 7.5, None), f, env_q_abs(1000.0, p), qr, 300.0, 5.0)["tj"][-1]
        for pr in (eu, ra):
            d = abs(pr.window(pr.init_state(340.0, 315.0, 7.5, None), f, env_q_abs(1000.0, p), qr, 300.0, 18.75)["tj"][-1] - ref)
            worst = max(worst, float(d))
    assert worst < 2.0, worst      # limite FOLGADO de sanidade (nao e meta de exatidao)


def test_oraculo_sem_margem_respeita_80C_antes_do_eclipse():
    r = _short("TEMS_SOLVER", t_end=2400.0, m_total=0.0)
    assert r["tj_max_C"] <= 80.0 + 0.5, r["tj_max_C"]       # tolerancia de discretizacao: DIAGNOSTICO, nao calibrada


def test_energia_e_faixas_na_malha_fechada():
    r = _short("R")
    assert r["f_min_MHz"] >= 300.0 - 1e-6 and r["f_max_MHz"] <= 1410.0 + 1e-6, r
    assert r["soc_min"] >= 0.70 - 1e-6 and r["energy_resid_rel"] < 1e-7 and r["bus_resid_W"] < 1e-6, r
    assert _short("S0")["t_soc_below_min_s"] > 0.0      # sem supervisor o piso de SoC e violado: o supervisor faz diferenca


def test_determinismo_por_semente():
    a, b = _short("R", "burst", 3, "calibrado", t_end=1200.0), _short("R", "burst", 3, "calibrado", t_end=1200.0)
    c = _short("R", "burst", 4, "calibrado", t_end=1200.0)
    assert _sem_tempo(a) == _sem_tempo(b), "mesma semente deveria reproduzir o resultado (exceto latencia)"
    assert _sem_tempo(a) != _sem_tempo(c), "semente diferente deveria mudar o resultado"


def test_barreira_ood_dispara_fallback():
    class OOD(_RC3):
        def domain_distance(self, case):
            return 0.5
    c = PROJETO_1P0
    cfg = scenario_cfg("nominal", 0, c, 1)
    cfg.update({"t_end": 600.0, "ctrl": "TEMS", "pred": OOD(plant_par(c), A_PL_PROJETO), "fallback": RC3(plant_par(c), A_PL_PROJETO)})
    assert run_closed_loop(cfg)["n_fallback"] > 0
    assert _short("TEMS_RC3", t_end=600.0, m_total=0.0)["n_fallback"] == 0


def test_margem_maior_nao_aumenta_tj():
    a, b = _short("TEMS_RC3", t_end=1500.0, m_total=0.0), _short("TEMS_RC3", t_end=1500.0, m_total=3.0)
    assert b["tj_max_C"] <= a["tj_max_C"] + 1e-6 and b["tops_frac_mean"] <= a["tops_frac_mean"] + 1e-6, (a["tj_max_C"], b["tj_max_C"])


# ---- exigem torch / arquivos do repositorio: PULADOS se ausentes
def test_pinn_api_presente():
    """Aponta divergencia entre a API ASSUMIDA por PINNPred e o modelo_pinn.py real (o assistente nao leu esse arquivo)."""
    pytest.importorskip("torch")
    mp = pytest.importorskip("modelo_pinn")
    faltam = [n for n in ("load", "W_S", "predict_state", "estimate_dtip", "RNG_", "DJR") if not hasattr(mp, n)]
    assert not faltam, "modelo_pinn.py nao tem: %s -> ajustar PINNPred em marco5_tems.py" % faltam


def test_checkpoints_do_repositorio_batem_com_o_contrato():
    pasta = os.environ.get("CKPT_DIR", os.path.join("simulacoes", "v2_auditado"))
    achou = False
    for nome, h in CHECKPOINTS_SHA256.items():
        p = os.path.join(pasta, nome)
        if os.path.exists(p):
            achou = True
            assert sha256_file(p) == h.upper(), nome
    if not achou:
        pytest.skip("nenhum .pt em %s (defina CKPT_DIR)" % pasta)


# ---- testes novos (R_SINT, sintonia, quantizacao, pue_last3, dt_ctrl, series) -- ESCRITOS SEM TER SIDO EXECUTADOS
def test_valida_e_carrega_sintonia(tmp_path):
    t = {"contrato_sha256": contract_hash(PROJETO_1P0), "seeds_cal": [100, 101],
         "sensores": {"sem_calibracao|none": {"viavel": True, "kp": 20.0, "ki": 0.5, "margem": 2.0}, "calibrado|trunc": {"viavel": False}}}
    ok = {"r_kp": 20.0, "r_ki": 0.5, "r_margin": 2.0}
    assert valida_sintonia(t, "p", [0, 1], "sem_calibracao", "none") == ok
    for args in (([0, 100], "sem_calibracao", "none"), ([0], "calibrado", "trunc"), ([0], "ideal", "none")):
        with pytest.raises(RuntimeError):
            valida_sintonia(t, "p", *args)
    with pytest.raises(RuntimeError):
        valida_sintonia(dict(t, contrato_sha256="0" * 64), "p", [0], "sem_calibracao", "none")
    p = tmp_path / "s.json"
    p.write_text(json.dumps(t), encoding="utf-8")
    assert carrega_sintonia(str(p), [0, 1], "sem_calibracao", "none") == ok
    assert chave_sintonia("calibrado", "trunc") == "calibrado|trunc"


def test_quantizacao_round_e_trunc():
    k = 273.15 + 79.6
    r = Sensors(sensor_cfg("ideal", quant_mode="round"), 0).read(k, 300.0)[0]
    f = Sensors(sensor_cfg("ideal", quant_mode="trunc"), 0).read(k, 300.0)[0]
    assert abs(r - (273.15 + 80.0)) < 1e-9 and abs(f - (273.15 + 79.0)) < 1e-9


def test_pue_last3_exige_janela_com_aquecimento():
    a, b = run_scenario("nominal", "S1", 0, n_orb=4), run_scenario("nominal", "S1", 0, n_orb=3)
    assert b["pue_last3"] is None and b["dE_over_ti_last3"] is None
    assert a["pue_last3"] > 1.0 and a["dE_over_ti_last3"] < 1e-3


def test_dt_ctrl_300_e_series_da_ultima_orbita():
    c = dict(PROJETO_1P0, nome="PROJETO_1P0_dt300", dt_ctrl=300.0)
    cfg = scenario_cfg("nominal", 0, c, 1)
    cfg.update({"t_end": 1500.0, "ctrl": "R", "series": True})
    r = run_closed_loop(cfg)
    assert r["n_dec"] == 5 and r["n_infeasible"] == 0 and r["infeasible_t_s"] == []
    u = ultima_orbita(r["series"], 1500.0)
    assert sorted(u) == ["f_MHz", "soc", "t_s", "tj_C", "tr_C"] and len({len(v) for v in u.values()}) == 1
    assert abs(u["t_s"][-1] - 1500.0) < 1e-6


def test_r_sint_com_os_mesmos_ganhos_reproduz_r():
    ov = {"r_kp": 20.0, "r_ki": 0.5, "r_margin": 2.0, "t_end": 1200.0}
    a = run_scenario("nominal", "R", 0, n_orb=1, sensor="calibrado", cfg_over=ov)
    b = run_scenario("nominal", "R_SINT", 0, n_orb=1, sensor="calibrado", cfg_over=ov)
    assert _sem_tempo(a) == _sem_tempo(b)


def test_sintonia_escolhe_maior_vazao_viavel_e_avisa_borda():
    L = [dict(kp=20.0, ki=0.5, margem=m, viavel=v, tops_last_medio=tp, tj_max_pior_C=tj)
         for m, v, tp, tj in ((0.0, False, 0.90, 80.2), (2.0, True, 0.87, 79.0), (3.0, True, 0.85, 78.0))]
    b, av = escolhe(L)
    assert b["margem"] == 2.0 and av == []
    b2, av2 = escolhe([dict(x, viavel=(x["margem"] == 3.0)) for x in L])
    assert b2["margem"] == 3.0 and len(av2) == 1
    assert escolhe([dict(x, viavel=False) for x in L])[0] is None


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))

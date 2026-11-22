"""test_janela_pue.py -- testes puros (sem simulacao, sem torch, sem pandas) de janela_pue.py e rodar_campanha_longa.py. ESCRITO SEM TER SIDO EXECUTADO.
Uso: python -m pytest test_janela_pue.py -q -rs
"""
import json
import os

import numpy as np

import janela_pue as jp
import rodar_campanha_longa as rl


def _serie(n=6, T=100.0, soc_fim=None, p_ti=1000.0, f=705.0):
    soc_fim = soc_fim or [0.5 + 0.05 * i for i in range(n)]
    t = np.arange(1.0, n * T + 1.0, 1.0)
    xs = np.concatenate(([0.0], np.arange(1, n + 1) * T))
    ys = np.concatenate(([0.5], soc_fim))
    return {"t": list(t), "h": [1.0] * len(t), "tj": [350.0] * len(t), "soc": list(np.interp(t, xs, ys)), "f": [f] * len(t), "p_ti": [p_ti] * len(t)}


def test_por_orbita_e_enriquece():
    ser = _serie(n=6, soc_fim=[0.6, 0.7, 0.8, 0.8, 0.8, 0.8])
    r = {"series": ser, "ecap_Wh": 1000.0}
    jp.enriquece(r, 100.0, 1410.0, 0.5, 6, 273.15)
    assert len(r["orb_soc_end"]) == 6 and abs(r["orb_soc_end"][0] - 0.6) < 1e-6
    assert abs(r["orb_tops"][0] - 0.5) < 1e-6 and abs(r["orb_tj_max_C"][0] - 76.85) < 1e-6
    assert abs(r["orb_ti_Wh"][0] - 1000.0 * 100.0 / 3600.0) < 1e-3
    assert "series" in r and r["orb_T_s"] == 100.0


def test_rel_ponta():
    soc, ti = [0.5, 0.6, 0.7, 0.8, 0.8], [1000.0] * 5
    assert abs(jp.rel_ponta(soc, ti, 1000.0, 3, 4) - 0.1) < 1e-12
    assert jp.rel_ponta(soc, ti, 1000.0, 3, 3) is None
    assert jp.rel_ponta(soc, [0.0] * 5, 1000.0, 3, 5) is None
    assert jp.rel_ponta([0.8] * 12, [1000.0] * 12, 1000.0, 9, 12) == 0.0


def test_avalia_run_final_e_fracao():
    r = {"orb_soc_end": [0.5, 0.6, 0.7, 0.8] + [0.8] * 8, "orb_ti_Wh": [1000.0] * 12, "ecap_Wh": 1000.0}
    a = jp.avalia_run(r, 3)
    assert a["N"] == 12 and a["n_tot"] == 9 and a["final_ok"] is True
    assert a["n_pass"] == 6 and a["max_rel"] > 1e-3          # n = 4..6 ainda veem a rampa
    curta = {"orb_soc_end": [0.8] * 5, "orb_ti_Wh": [1000.0] * 5, "ecap_Wh": 1000.0}
    assert jp.avalia_run(curta, 9)["final"] is None and jp.avalia_run({}, 9) is None


def _grava_parts(out, linhas):
    d = os.path.join(out, "parts", "1")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "x.json"), "w", encoding="utf-8") as fh:
        json.dump({"resultados": linhas}, fh)


def _reg(ctrl, plano, N=20):
    soc = [0.84] * N if plano else [0.80 + 0.04 * (i % 2) for i in range(N)]   # |dSoC|=0,04 -> rel ~2e-3 > 1e-3 (falha)
    return {"cenario": "nominal", "sensor": "sem_calibracao", "ctrl": ctrl, "seed": 10, "mtotal": None, "dt_ctrl": 60.0,
            "orb_soc_end": soc, "orb_ti_Wh": [57000.0] * N, "ecap_Wh": 26593.9}


def test_secao_pue(tmp_path):
    out = str(tmp_path)
    _grava_parts(out, [_reg("TEMS_RC3", True), _reg("TEMS_PINN", False), {"cenario": "nominal", "ctrl": "S0", "seed": 10}])
    s = rl.secao_pue(out, 9)
    assert "JANELA DE 9 ORBITAS" in s
    j = json.load(open(os.path.join(out, "m5_pue_janela.json")))
    assert j["por_ctrl"]["TEMS_RC3"]["final_ok"] == 1 and j["por_ctrl"]["TEMS_PINN"]["final_ok"] == 0
    assert j["por_ctrl"]["S0"]["sem_orb"] == 1
    assert any(g.startswith("nominal|sem_calibracao|none|nan|60.0|TEMS_RC3") for g in j["grupos"])
    assert "TEMS_RC3,1,0,0,1.0000" in s


def test_extrai_janela():
    assert rl.extrai_janela(["--norb", "20", "--janela", "9", "--out", "x"]) == (["--norb", "20", "--out", "x"], 9)
    assert rl.extrai_janela(["--norb", "20"]) == (["--norb", "20"], 9)


def test_roda_um_pula_existentes(tmp_path):
    p = str(tmp_path / "a.json")
    open(p, "w").write("{}")
    assert rl.roda_um(p, [], "res", str(tmp_path))[1] == "pulado"
    q = str(tmp_path / "b.json")
    open(q + ".divergiu", "w").write("{}")
    assert rl.roda_um(q, [], "res", str(tmp_path))[1] == "pulado_divergiu"


def test_confere_config_inclui_worker_e_janela(tmp_path):
    out = str(tmp_path)
    cfg = {"norb": 20, "norb_sint": 6, "seeds": [10, 11, 12], "burst_per_orbit": 12, "quick": False, "ckpt": "v2c.pt"}
    rl.JANELA[0] = 9
    rl.confere_config(out, cfg, False, False)
    g = json.load(open(os.path.join(out, "campanha_config.json")))
    assert g["worker"] == "marco5_orbitas.py" and g["janela_pue"] == 9 and g["norb"] == 20
    rl.JANELA[0] = 4
    try:
        rl.confere_config(out, cfg, False, False)
        assert False, "deveria recusar janela diferente da gravada"
    except SystemExit:
        pass
    rl.JANELA[0] = 9

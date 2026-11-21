"""test_resumo_v8.py -- teste do --resumo (v8). ESCRITO SEM TER SIDO EXECUTADO. Sem torch, sem simulacao.
Garante que linhas com mtotal ausente (nominal) entram no AGREGADO junto com as do pareto.
Uso: python -m pytest test_resumo_v8.py -q   (exige rodar_campanha_marco5.py = v8 com o nome canonico)
"""
import json
import os

import rodar_campanha_marco5 as rc


def _r(cen, ctrl, seed, tj, tops, mt=None, com_quant=True):
    d = {"cenario": cen, "sensor": "sem_calibracao", "ctrl": ctrl, "seed": seed, "mtotal": mt, "dt_ctrl": 60.0,
         "tj_max_C": tj, "tj_max_last_C": tj, "viol_int_Ks": 0.0, "t_above_80_s": 0.0, "tops_frac_last": tops, "soc_min": 0.7,
         "n_infeasible": 0, "n_fallback": 0, "dE_over_ti_last3": 5e-4, "pred_err_rms_last_K": 0.5, "decision_ms_mean": 10.0, "infeasible_diag": []}
    if com_quant:
        d["quant"] = "none"
    return d


def _grava(out, bloco, nome, linhas):
    p = os.path.join(out, "parts", str(bloco))
    os.makedirs(p, exist_ok=True)
    with open(os.path.join(p, nome), "w", encoding="utf-8") as fh:
        json.dump({"contrato_sha256": "X", "ckpt_sha256": "Y", "rtune_sha256": "Z", "resultados": linhas}, fh)


def test_agregado_nao_descarta_linhas_sem_mtotal(tmp_path):
    out = str(tmp_path)
    _grava(out, 1, "nom.json", [_r("nominal", "R_SINT", 10, 79.5, 0.88, com_quant=False), _r("nominal", "TEMS_RC3", 10, 78.4, 0.878)])
    _grava(out, 5, "par.json", [_r("nominal", "TEMS_RC3", 10, 79.0, 0.879, mt=0.5), _r("nominal", "TEMS_RC3", 11, 79.1, 0.880, mt=0.5)])
    s = rc.gerar_resumo(out, imprimir=False)
    assert "# cobertura: 4 de 4 execucoes" in s, s[:1500]
    assert "nominal,nominal,sem_calibracao,none,nan,60.0,R_SINT" in s
    assert "nominal,nominal,sem_calibracao,none,nan,60.0,TEMS_RC3" in s
    assert "pareto,nominal,sem_calibracao,none,0.5,60.0,TEMS_RC3" in s


def test_texto_chave():
    assert rc._texto_chave(None, "nan") == "nan"
    assert rc._texto_chave(float("nan"), "nan") == "nan"
    assert rc._texto_chave(60.0, "x") == "60.0"
    assert rc._texto_chave("trunc", "none") == "trunc"

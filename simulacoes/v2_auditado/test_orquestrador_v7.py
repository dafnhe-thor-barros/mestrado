"""test_orquestrador_v7.py -- testes do orquestrador v7 (puros, sem torch, sem rodar simulacao). ESCRITO SEM TER SIDO EXECUTADO.
Importa rodar_campanha_marco5 (salve o v7 com esse nome) e sintonia_reativo (so para conferir a grade; exige contrato_parametros e marco5_tems).
Uso: python -m pytest test_orquestrador_v7.py -q
"""
import json
import os

import pytest

import rodar_campanha_marco5 as rc


def _cfg(seeds, chunk=5):
    return {"seeds": list(seeds), "seeds_burst": list(seeds), "norb": 6, "ckpt": "v2c.pt", "rtune": "x.json", "chunk": chunk, "burst_per_orbit": 12,
            "sensores": ["calibrado", "sem_calibracao"], "quants": ["trunc", "round"], "mt": [0, 0.5, 1, 2, 3], "pert": rc.PERT}


def test_valida_seeds():
    assert rc.valida_seeds([10, 11, 12, 13, 14]) == [10, 11, 12, 13, 14]
    for ruim in ([1, 2], [1, 1, 2], [0, 1, 100], [101, 5, 6]):
        with pytest.raises(SystemExit):
            rc.valida_seeds(ruim)


def test_lista_de_argumentos():
    a = ["--scenarios", "x", "y", "--seeds", "3", "4", "--sensor", "s"]
    assert rc._lista(a, "--scenarios") == ["x", "y"] and rc._lista(a, "--seeds") == ["3", "4"] and rc._lista(a, "--zzz") == []


def test_monta_s1_por_semente_e_sem_semente_de_calibracao():
    B = rc.monta(_cfg([10, 11, 12, 13, 14], chunk=5))
    s1 = [p for p in B[4] if "pert_S1_" in p[0]]
    assert len(s1) == 6 * 5, len(s1)                                     # 6 cenarios x 5 sementes, um processo por semente
    for nome, args, tipo in s1:
        assert len(rc._lista(args, "--seeds")) == 1 and rc._lista(args, "--ctrls") == ["S1"]
    outros = [p for p in B[4] if "pert_S1_" not in p[0]]
    assert len(outros) == 6 * 3                                          # R_SINT, RC3, PINN: 1 parte (chunk 5) por cenario
    assert len([p for p in B[6] if "ood_S1_" in p[0]]) == 5
    for b in B.values():
        for nome, args, tipo in b:
            sd = {int(s) for s in rc._lista(args, "--seeds")}
            assert not (sd & set(rc.SEEDS_CAL)), nome
    assert any("ref_det" in p[0] for p in B[1])


def test_roda_um_pula_partes_existentes_e_divergidas(tmp_path):
    p = str(tmp_path / "a.json")
    open(p, "w").write("{}")
    assert rc.roda_um(p, [], "res", str(tmp_path))[1] == "pulado"
    q = str(tmp_path / "b.json")
    open(q + ".divergiu", "w").write("{}")
    assert rc.roda_um(q, [], "res", str(tmp_path))[1] == "pulado_divergiu"


def test_confere_config_nova_e_retomada(tmp_path):
    out = str(tmp_path)
    cfg = {"norb": 6, "norb_sint": 6, "seeds": [10, 11, 12], "burst_per_orbit": 12, "quick": False, "ckpt": "v2c.pt"}
    rc.confere_config(out, cfg, nova=True, ignora=False)                 # pasta vazia: grava a config
    assert json.load(open(os.path.join(out, "campanha_config.json")))["norb"] == 6
    rc.confere_config(out, cfg, nova=False, ignora=False)                # mesma config: ok
    with pytest.raises(SystemExit):
        rc.confere_config(out, dict(cfg, norb=3), nova=False, ignora=False)
    rc.confere_config(out, dict(cfg, norb=3), nova=False, ignora=True)   # ignorar a recusa
    os.makedirs(os.path.join(out, "parts", "1"))
    open(os.path.join(out, "parts", "1", "x.json"), "w").write("{}")
    with pytest.raises(SystemExit):
        rc.confere_config(out, cfg, nova=True, ignora=False)             # --nova com partes existentes


def test_estado_da_clausula_estrita_x_tolerancia():
    assert rc._estado_clausula(False, True, ["EMPATE", "EMPATE"]) == "NAO_ACIONADA"
    assert rc._estado_clausula(True, True, ["EMPATE", "EMPATE"]) == "ACIONADA"
    assert rc._estado_clausula(True, True, ["EMPATE", "TEMS_MAIOR"]) == "PARCIAL"
    assert rc._estado_clausula(True, True, []) == "INDETERMINADA"


def _lin(bloco, cen, ctrl, sensor, seed, tj, tops, tjl=None, mt=None, quant="none", due=5e-4, soc=0.7):
    return {"bloco": bloco, "cenario": cen, "ctrl": ctrl, "sensor": sensor, "seed": seed, "tj_max_C": tj, "tj_max_last_C": tj if tjl is None else tjl, "tops_frac_last": tops,
            "soc_min": soc, "mtotal": mt, "quant": quant, "dt_ctrl": 60.0, "n_infeasible": 0, "infeasible_t_s": [], "dE_over_ti_last3": due, "decision_ms_mean": 10.0}


def _sint(rsint_pior=79.5, pinn_due=-5e-3):
    L = []
    delta = [0.001, -0.001, 0.002, -0.002, 0.0]          # R_SINT em torno do ponto interpolado do Pareto (0,88611): diferencas com media 0 => EMPATE
    for sd in range(5):
        tj_r = rsint_pior if sd == 0 else 79.5
        L.append(_lin("nominal", "nominal", "R_SINT", "sem_calibracao", sd, tj_r, 0.88611 + delta[sd], tjl=79.5))
        L.append(_lin("nominal", "nominal", "TEMS_RC3", "sem_calibracao", sd, 78.2, 0.880 + 0.001 * sd))
        L.append(_lin("nominal", "nominal", "TEMS_PINN", "sem_calibracao", sd, 78.2, 0.8805 + 0.001 * sd, due=pinn_due))
        for m, tj, tp in ((0.0, 80.2, 0.890), (1.0, 79.3, 0.885), (2.0, 78.4, 0.880)):
            L.append(_lin("pareto", "nominal", "TEMS_RC3", "sem_calibracao", sd, tj - 0.1, tp, tjl=tj, mt=m))
        L.append(_lin("pert", "r_int_alto", "S1", "sem_calibracao", sd, 70.0, 0.85))           # S1 nao viola por T_j
        L.append(_lin("pert", "c_mod_mais30", "S1", "sem_calibracao", sd, 70.0, 0.85))
        for ct in ("TEMS_RC3", "TEMS_PINN", "R_SINT"):
            L.append(_lin("pert", "combinada", ct, "sem_calibracao", sd, 82.0, 0.74))             # TEMS viola onde o S1 divergiu
            L.append(_lin("pert", "r_int_alto", ct, "sem_calibracao", sd, 78.0, 0.88))
            L.append(_lin("pert", "c_mod_mais30", ct, "sem_calibracao", sd, 78.0, 0.88))
    return L


def test_p2_conta_divergencia_do_s1_e_mostra_as_duas_leituras():
    div = [{"bloco": "pert", "parte": "pert_S1_combinada_s%d-%d.json" % (s, s), "cenarios": ["combinada"], "ctrls": ["S1"], "seeds": [s]} for s in range(5)]
    v = rc.calcula_veredito(_sint(), diverg=div)
    assert v["P2"]["cenarios_em_que_S1_viola"] == ["combinada"] and v["P2"]["status"] == "NAO_ATENDIDO"
    assert v["P2"]["S1_viola_por_T_j"] == [] and v["P2"]["leitura_sem_divergencia"]["status"] == "NAO_APLICAVEL"
    assert v["P2"]["S1_divergiu_por_cenario"] == {"combinada": [0, 1, 2, 3, 4]}
    assert rc.calcula_veredito(_sint())["P2"]["status"] == "NAO_APLICAVEL"        # sem divergencia e sem violacao do S1


def test_clausula_estrita_nao_aciona_com_80_12_mas_sensibilidade_sim():
    v = rc.calcula_veredito(_sint(rsint_pior=80.12))
    assert v["P1"]["R_SINT_sem_violacao"] is False and v["P1"]["R_SINT_pior_tj_C"] == pytest.approx(80.12)
    assert v["P4"]["detalhe"][0]["classe"] == "EMPATE"
    assert v["clausula_de_falsificacao"]["estado"] == "NAO_ACIONADA"
    assert v["clausula_de_falsificacao"]["sensibilidade_tol_0p25K_posthoc"]["estado"] == "ACIONADA"
    assert rc.calcula_veredito(_sint(rsint_pior=79.5))["clausula_de_falsificacao"]["estado"] == "ACIONADA"


def test_p5_ressalva_de_periodicidade_e_pue_com_modulo():
    v = rc.calcula_veredito(_sint(pinn_due=-5e-3))                      # dE negativo grande: so reprova com o valor absoluto
    g = [r for r in v["P5"]["detalhe"] if r["bloco"] == "nominal"][0]
    assert g["pue_ok_PINN"] == 0.0 and g["pue_ok_RC3"] == 1.0 and g["ressalva_periodicidade"] is True
    assert v["PUE_periodico_fracao_por_ctrl(|dE|<%g)" % rc.LIM_PUE]["TEMS_PINN"] == pytest.approx(0.75)      # 5 execucoes nominais reprovam, 15 de pert (due padrao 5e-4) passam
    assert v["PUE_periodico_fracao_por_ctrl_SEM_abs(comparavel_com_piloto)"]["TEMS_PINN"] == 1.0


def test_grade_estendida_da_sintonia():
    sr = pytest.importorskip("sintonia_reativo")
    assert sr.GRADE["kp"] == (2.0, 5.0, 20.0, 60.0, 120.0) and sr.GRADE["ki"] == (0.1, 0.5, 1.5, 3.0) and sr.GRADE["margem"] == (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
    assert min(sr.GRADE["margem"]) == 0.0           # sem margens negativas (decisao do usuario)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))

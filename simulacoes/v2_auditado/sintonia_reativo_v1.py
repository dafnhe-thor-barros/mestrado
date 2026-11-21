"""sintonia_reativo.py -- Marco 5: sintonia JUSTA do controlador reativo (R) em sementes DISJUNTAS das de teste.
Criterio (declarado ANTES de ver resultados; [HIPOTESE DECLARADA]): para cada (sensor, quantizacao), escolher (kp, ki, margem) que MAXIMIZA a vazao media
da ultima orbita (tops_frac_last, media sobre cenarios x sementes de calibracao) SUJEITO a: nenhum run de calibracao com T_j > 80 C (viol_int_Ks == 0, estrito)
e nenhum run abaixo do piso de SoC. O R recebe a mesma exigencia de seguranca do TEMS (zero tempo acima de 80 C) e a margem e escolhida POR sensor.
Se nenhum ponto da grade for viavel, a entrada fica viavel=false e o R_SINT se recusa a rodar naquele sensor (nao existe sintonia "quase boa").
Uso:  python sintonia_reativo.py --sensors ideal calibrado sem_calibracao --quants none --seeds-cal 100 101 102 --norb 3 --jobs 4 --out sintonia_R.json
      python sintonia_reativo.py --quick ...   (grade reduzida, so para teste de fumaca)
Sementes de calibracao DEVEM ser disjuntas das de teste (marco5_tems.carrega_sintonia recusa sobreposicao).
"""
import argparse
import itertools
import json

from contrato_parametros import PROJETO_1P0, SENSOR_MARGEM_K, contract_hash
from marco5_tems import run_scenario, chave_sintonia

SEEDS_CAL_DEFAULT = (100, 101, 102)
CENARIOS_CAL = ("nominal", "burst")
GRADE = {"kp": (5.0, 20.0, 60.0), "ki": (0.1, 0.5, 1.5), "margem": (0.0, 0.5, 1.0, 1.5, 2.0, 3.0)}
GRADE_QUICK = {"kp": (20.0,), "ki": (0.5,), "margem": (0.0, 1.0, 2.0, 3.0)}
CAMPOS = ("viol_int_Ks", "t_above_80_s", "t_soc_below_min_s", "soc_min", "tops_frac_last", "tj_max_C", "tj_max_last_C", "dE_bess_last_Wh", "energy_resid_rel")


def _um_run(arg):
    sc, seed, sensor, quant, norb, kp, ki, margem = arg
    qkw = {} if quant == "none" else {"quant_mode": quant}
    r = run_scenario(sc, "R", seed, n_orb=norb, sensor=sensor, cfg_over={"r_kp": kp, "r_ki": ki, "r_margin": margem}, **qkw)
    return {k: r[k] for k in CAMPOS}


def avalia_grade(sensor, quant, grade, seeds, norb, ex=None):
    """Roda TODOS os (ponto da grade x cenario x semente) de uma vez (melhor balanceamento entre processos)."""
    pontos = list(itertools.product(grade["kp"], grade["ki"], grade["margem"]))
    tarefas = [(sc, sd, sensor, quant, norb, kp, ki, m) for (kp, ki, m) in pontos for sc in CENARIOS_CAL for sd in seeds]
    rs = list(ex.map(_um_run, tarefas)) if ex is not None else [_um_run(t) for t in tarefas]
    n = len(CENARIOS_CAL) * len(seeds)
    linhas = []
    for i, (kp, ki, m) in enumerate(pontos):
        r = rs[i * n:(i + 1) * n]
        ok = all(x["viol_int_Ks"] <= 0.0 and x["t_soc_below_min_s"] <= 0.0 for x in r)
        linhas.append({"kp": kp, "ki": ki, "margem": m, "viavel": bool(ok), "n_runs": len(r),
                       "tops_last_medio": sum(x["tops_frac_last"] for x in r) / len(r),
                       "tj_max_pior_C": max(x["tj_max_C"] for x in r), "viol_int_pior_Ks": max(x["viol_int_Ks"] for x in r),
                       "t_above_80_pior_s": max(x["t_above_80_s"] for x in r), "soc_min_pior": min(x["soc_min"] for x in r),
                       "dE_last_pior_abs_Wh": max(abs(x["dE_bess_last_Wh"]) for x in r)})
    return linhas


def escolhe(linhas):
    """Maior vazao entre os viaveis; empate -> menor T_j maximo. Retorna (melhor | None, avisos)."""
    viaveis = [x for x in linhas if x["viavel"]]
    if not viaveis:
        return None, ["nenhum ponto da grade e viavel: ampliar a grade (margem maior / ganhos menores) antes de comparar com o TEMS"]
    best = max(viaveis, key=lambda x: (round(x["tops_last_medio"], 6), -x["tj_max_pior_C"]))
    av = []
    for nome in ("kp", "ki", "margem"):
        vals = sorted(set(x[nome] for x in linhas))
        if len(vals) > 1 and best[nome] in (vals[0], vals[-1]):
            av.append("melhor %s=%g esta na BORDA da grade %s: o otimo pode estar fora" % (nome, best[nome], (vals[0], vals[-1])))
    return best, av


def sintoniza(sensores, quants, seeds_cal, norb, grade, jobs=1):
    out = {"contrato_sha256": contract_hash(PROJETO_1P0), "seeds_cal": [int(s) for s in seeds_cal], "cenarios_cal": list(CENARIOS_CAL), "norb": norb,
           "grade": {k: list(v) for k, v in grade.items()},
           "criterio": "max tops_frac_last medio s.t. viol_int_Ks == 0 e SoC >= piso em todos os runs de calibracao", "sensores": {}}
    ex = None
    if jobs > 1:
        from concurrent.futures import ProcessPoolExecutor
        ex = ProcessPoolExecutor(max_workers=jobs)
    try:
        for s in sensores:
            for q in quants:
                linhas = avalia_grade(s, q, grade, seeds_cal, norb, ex)
                best, av = escolhe(linhas)
                ent = {"sensor": s, "quant": q, "viavel": best is not None, "avisos": av, "grade_completa": linhas, "margem_sensor_TEMS_K": SENSOR_MARGEM_K[s]}
                if best:
                    ent.update({k: best[k] for k in ("kp", "ki", "margem", "tops_last_medio", "tj_max_pior_C")})
                out["sensores"][chave_sintonia(s, q)] = ent
    finally:
        if ex is not None:
            ex.shutdown()
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description="Sintonia justa do reativo (sementes de calibracao disjuntas das de teste)")
    ap.add_argument("--sensors", nargs="+", default=["ideal", "calibrado", "sem_calibracao"], choices=sorted(SENSOR_MARGEM_K))
    ap.add_argument("--quants", nargs="+", default=["none"], choices=["none", "floor", "trunc", "round"])
    ap.add_argument("--seeds-cal", nargs="+", type=int, default=list(SEEDS_CAL_DEFAULT))
    ap.add_argument("--norb", type=int, default=3)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="sintonia_R.json")
    a = ap.parse_args(argv)
    res = sintoniza(a.sensors, a.quants, a.seeds_cal, a.norb, GRADE_QUICK if a.quick else GRADE, a.jobs)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        json.dump(res, fh, ensure_ascii=False)
    for k, e in res["sensores"].items():
        print(k, ("VIAVEL kp=%g ki=%g margem=%g tops=%.4f Tjmax_pior=%.2f" % (e["kp"], e["ki"], e["margem"], e["tops_last_medio"], e["tj_max_pior_C"])) if e["viavel"] else "INVIAVEL", e["avisos"])
    return res


if __name__ == "__main__":
    main()

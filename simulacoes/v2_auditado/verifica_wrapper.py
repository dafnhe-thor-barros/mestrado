"""verifica_wrapper.py -- confere que marco5_orbitas.py (wrapper) reproduz o marco5_tems.py original nos campos ja existentes.
ESCRITO SEM TER SIDO EXECUTADO. Roda 2 simulacoes curtas (TEMS_RC3, nominal, sem_calibracao, semente 10, 4 orbitas) e compara.
Ignora so: chaves que comecam com 'orb_' (novas) e decision_ms_* (tempo de parede, varia). Qualquer outra diferenca = FALHA.
Uso (dentro de simulacoes/v2_auditado, venv ativo): python verifica_wrapper.py
"""
import json
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
ARGS = ["--scenarios", "nominal", "--ctrls", "TEMS_RC3", "--seeds", "10", "--sensor", "sem_calibracao", "--norb", "4", "--ckpt", "v2c.pt"]
IGNORA_PREFIXO = ("orb_",)
IGNORA_EXATO = ("decision_ms_mean", "decision_ms_p99")


def roda(script, saida):
    if os.path.exists(saida):
        os.remove(saida)
    r = subprocess.run([sys.executable, script] + ARGS + ["--out", saida], capture_output=True, text=True, cwd=AQUI)
    if r.returncode != 0 or not os.path.exists(saida):
        raise SystemExit("FALHOU %s (codigo %s):\n%s" % (script, r.returncode, (r.stderr or r.stdout)[-1500:]))
    return json.load(open(saida, encoding="utf-8"))


def main():
    a = roda("marco5_tems.py", os.path.join(AQUI, "_vw_orig.json"))["resultados"][0]
    b = roda("marco5_orbitas.py", os.path.join(AQUI, "_vw_wrap.json"))["resultados"][0]
    dif = []
    for k in sorted(set(a) | set(b)):
        if k.startswith(IGNORA_PREFIXO) or k in IGNORA_EXATO:
            continue
        if a.get(k) != b.get(k):
            dif.append(k)
    novas = sorted(k for k in b if k.startswith("orb_"))
    print("chaves novas no wrapper:", novas)
    print("n orbitas nas listas:", {k: len(b[k]) for k in novas if isinstance(b[k], list)})
    if dif:
        print("DIFERENCAS em campos originais:", dif)
        for k in dif[:10]:
            print(" ", k, "orig=", str(a.get(k))[:120], "| wrapper=", str(b.get(k))[:120])
        raise SystemExit(1)
    print("EQUIVALENTE: todos os campos originais (exceto decision_ms_*) sao identicos.")
    for p in ("_vw_orig.json", "_vw_wrap.json"):
        try:
            os.remove(os.path.join(AQUI, p))
        except OSError:
            pass


if __name__ == "__main__":
    main()

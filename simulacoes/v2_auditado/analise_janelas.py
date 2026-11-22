"""analise_janelas.py -- compara janelas candidatas para o criterio de PUE periodico a partir dos conv_*.json de teste_convergencia_orbitas.py.
SO LEITURA: nao simula e nao importa marco5_tems. ESCRITO SEM TER SIDO EXECUTADO.

Para cada arquivo, janela k e variante, e cada numero de orbitas n em que a simulacao poderia terminar:
  ponta: |SoC_fim(n) - SoC_fim(n-k)| * ecap / soma(E_TI das orbitas n-k+1..n)   (mesma definicao de marco5_tems.summarize.janela, L559-571; exige n >= k+1)
  media: |media(SoC_fim das k ultimas orbitas) - media(SoC_fim das k orbitas anteriores)| * ecap / soma(E_TI das k ultimas)   (exige n >= 2k)  [variante NOVA, so deste script]
Limiar: 1e-3 (LIM_PUE). Informa quantos n passam, o maior valor e os n que falham. As janelas e a variante 'media' sao escolha deste script (pos-hoc).
Uso:
  python analise_janelas.py conv_curto --nmin 5
  python analise_janelas.py conv_curto --janelas 3 4 5 8 9 --nmin 9
  python analise_janelas.py --autoteste
"""
import argparse
import glob
import json
import os
import sys

LIM = 1e-3


def rel_janela(soc, ti, ecap, k, n, var):
    """soc[i], ti[i]: SoC no fim e E_TI [Wh] da orbita i+1. Retorna o valor relativo ou None se n for pequeno demais."""
    if var == "ponta":
        if n < k + 1:
            return None
        a, b = soc[n - 1], soc[n - k - 1]
    elif var == "media":
        if n < 2 * k:
            return None
        a = sum(soc[n - k:n]) / k
        b = sum(soc[n - 2 * k:n - k]) / k
    else:
        raise ValueError(var)
    tt = sum(ti[n - k:n])
    return abs(a - b) * ecap / tt if tt > 0 else None


def analisa(path, janelas, nmin):
    d = json.load(open(path, encoding="utf-8"))
    if "por_orbita" not in d:
        print("%s: sem por_orbita (%s)" % (os.path.basename(path), d.get("erro", "?")))
        return
    rows, ecap = d["por_orbita"], d["ecap_Wh"]
    soc, ti = [r["soc_end"] for r in rows], [r["ti_Wh"] for r in rows]
    N = len(rows)
    print("\n== %s | soc0=%.2f norb=%d | n a partir de %d ==" % (os.path.basename(path), d["soc0"], N, nmin))
    print("var    k  n_validos  passam  maior_valor  n_que_falham")
    for var in ("ponta", "media"):
        for k in janelas:
            vals = []
            for n in range(max(nmin, 1), N + 1):
                v = rel_janela(soc, ti, ecap, k, n, var)
                if v is not None:
                    vals.append((n, v))
            if not vals:
                print("%-6s %2d  (nenhum n valido)" % (var, k))
                continue
            falha = [n for n, v in vals if v >= LIM]
            print("%-6s %2d  %9d  %6d  %11.2e  %s" % (var, k, len(vals), len(vals) - len(falha), max(v for _, v in vals), falha))


def autoteste():
    soc, ti = [0.5, 0.6, 0.7, 0.8, 0.8], [1000.0] * 5
    assert abs(rel_janela(soc, ti, 1000.0, 3, 4, "ponta") - 0.1) < 1e-12
    assert rel_janela(soc, ti, 1000.0, 3, 3, "ponta") is None
    assert abs(rel_janela(soc, ti, 1000.0, 2, 4, "media") - 0.1) < 1e-12
    assert rel_janela(soc, ti, 1000.0, 2, 3, "media") is None
    assert rel_janela([0.8] * 6, [1000.0] * 6, 1000.0, 3, 6, "ponta") == 0.0
    print("autoteste OK")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("pasta", nargs="?", default=None)
    ap.add_argument("--janelas", nargs="+", type=int, default=[1, 2, 3, 4, 5, 8, 9])
    ap.add_argument("--nmin", type=int, default=5, help="menor numero de orbitas simuladas considerado como fim possivel")
    ap.add_argument("--autoteste", action="store_true")
    a = ap.parse_args(argv)
    if a.autoteste:
        autoteste()
        return 0
    if not a.pasta:
        raise SystemExit("informe a pasta com os conv_*.json")
    arqs = sorted(glob.glob(os.path.join(a.pasta, "conv_*.json")))
    if not arqs:
        raise SystemExit("nenhum conv_*.json em %s" % a.pasta)
    for p in arqs:
        analisa(p, a.janelas, a.nmin)
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""rodar_campanha_longa.py -- orquestrador para campanhas com 20+ orbitas e janela de PUE periodico de k orbitas (padrao 9).
ESCRITO SEM TER SIDO EXECUTADO. NAO edita rodar_campanha_marco5.py (v8): importa o modulo e troca tres funcoes (roda_um, confere_config, gerar_resumo).

Diferencas em relacao ao v8:
  1) os processos filhos rodam marco5_orbitas.py (wrapper) em vez de marco5_tems.py; a latencia (--latencia) continua no marco5_tems.py original;
  2) campanha_config.json grava tambem worker e janela_pue (o cache continua nao olhando o numero de orbitas: use pasta NOVA para cada --norb);
  3) o --resumo acrescenta a secao 'PUE PERIODICO, JANELA DE k ORBITAS' e grava m5_pue_janela.json; a fracao original (dE_over_ti_last3, janela 3) continua no resumo.
  4) opcao nova --janela K (padrao 9), tratada aqui e removida antes de chamar o main do v8.
Uso:
  python rodar_campanha_longa.py --dry-run --blocos 1 2 3 4 5 6 7 --norb 20 --seeds 10 11 12 13 14 --chunk 1 --out m5_n20
  python rodar_campanha_longa.py --blocos 1 2 3 4 5 6 7 --norb 20 --norb-sint 6 --seeds 10 11 12 13 14 --jobs 5 --chunk 1 --janela 9 --out m5_n20
  python rodar_campanha_longa.py --resumo --out m5_n20
"""
import glob
import json
import os
import subprocess
import sys
import time

import rodar_campanha_marco5 as rc
import janela_pue as jp

WORKER = "marco5_orbitas.py"
JANELA = [9]


def roda_um(parte, args, tipo, cwd):
    """Copia de rc.roda_um (v8) trocando o script filho por WORKER."""
    if os.path.exists(parte):
        return parte, "pulado", 0.0
    if os.path.exists(parte + ".divergiu"):
        return parte, "pulado_divergiu", 0.0
    t0 = time.time()
    tmp = parte + ".tmp"
    cmd = list(args) + (["--save_traces", tmp, "--out", parte + ".main.tmp"] if tipo == "trace" else ["--out", tmp])
    r = subprocess.run([sys.executable, WORKER] + cmd, capture_output=True, text=True, env=rc.ENV, cwd=cwd)
    err = r.stderr or ""
    with open(parte + ".err", "w", encoding="utf-8", newline="\n") as fh:
        fh.write(err[-6000:])
    if r.returncode == 0 and os.path.exists(tmp):
        os.replace(tmp, parte)
        if os.path.exists(parte + ".main.tmp"):
            os.remove(parte + ".main.tmp")
        return parte, "ok", time.time() - t0
    if rc.DIVERG_TXT in err or rc.DIVERG_TXT in (r.stdout or ""):
        meta = {"parte": os.path.basename(parte), "motivo": rc.DIVERG_TXT, "cenarios": rc._lista(args, "--scenarios"), "ctrls": rc._lista(args, "--ctrls"),
                "seeds": [int(s) for s in rc._lista(args, "--seeds")], "sensor": (rc._lista(args, "--sensor") or [None])[0], "args": [str(a) for a in args]}
        with open(parte + ".divergiu", "w", encoding="utf-8", newline="\n") as fh:
            json.dump(meta, fh, ensure_ascii=False)
        return parte, "DIVERGIU", time.time() - t0
    return parte, "ERRO: " + ((r.stderr or r.stdout or "")[-700:]), time.time() - t0


_orig_confere = rc.confere_config


def confere_config(out, cfg_id, nova, ignora):
    return _orig_confere(out, dict(cfg_id, worker=WORKER, janela_pue=JANELA[0]), nova, ignora)


def ler_linhas(out):
    """Registros de execucao: parts/1..6 se existirem; senao m5_campanha_*.json."""
    linhas = []
    arqs = sorted(glob.glob(os.path.join(out, "parts", "[1-6]", "*.json")))
    if not arqs:
        arqs = sorted(glob.glob(os.path.join(out, "m5_campanha_*.json")))
    for f in arqs:
        d = json.load(open(f, encoding="utf-8"))
        linhas += d.get("resultados", [])
    return linhas


def ler_janela(out):
    p = os.path.join(out, "campanha_config.json")
    if os.path.exists(p):
        try:
            return int(json.load(open(p, encoding="utf-8")).get("janela_pue", JANELA[0]))
        except Exception:
            pass
    return JANELA[0]


def secao_pue(out, k):
    linhas = ler_linhas(out)
    res = jp.resume(linhas, k)
    with open(os.path.join(out, "m5_pue_janela.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(res, fh, ensure_ascii=False)
    L = ["##### PUE PERIODICO, JANELA DE %d ORBITAS (EMENDA POS-HOC; |dE|/E_TI < %g; exige n >= %d orbitas) #####" % (k, jp.LIM_PUE, k + 1),
         "execucoes lidas: %d. 'frac_passa_no_n_final' = fracao das execucoes que passam com a simulacao terminando na ultima orbita;" % len(linhas),
         "'frac_media_dos_n_que_passam' = media, por execucao, da fracao dos n = %d..N em que a janela passaria. O criterio antigo (janela 3) segue acima." % (k + 1),
         "--- por controlador", jp.CABECALHO_CSV]
    for c, e in sorted(res["por_ctrl"].items()):
        L.append(jp.linha_csv(c, e))
    L += ["--- por grupo (cenario|sensor|quant|mtotal|dt_ctrl|ctrl)", jp.CABECALHO_CSV]
    for g, e in sorted(res["grupos"].items()):
        L.append(jp.linha_csv(g, e))
    return "\n".join(L)


_orig_resumo = rc.gerar_resumo


def gerar_resumo(out, imprimir=True):
    s = _orig_resumo(out, False)
    try:
        s = s + "\n\n" + secao_pue(out, ler_janela(out))
    except Exception as e:
        s = s + "\n\n##### PUE PERIODICO, JANELA K #####\nERRO nesta secao: %s: %s" % (type(e).__name__, str(e)[:300])
    with open(os.path.join(out, "m5_resumo_para_colar.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(s)
    if imprimir:
        print("##### INICIO #####\n" + s + "\n##### FIM (%d caracteres) #####" % len(s))
    return s


def extrai_janela(argv):
    argv = list(argv)
    k = 9
    if "--janela" in argv:
        i = argv.index("--janela")
        if i + 1 >= len(argv):
            raise SystemExit("--janela exige um inteiro")
        k = int(argv[i + 1])
        del argv[i:i + 2]
    if k < 1:
        raise SystemExit("--janela >= 1")
    return argv, k


rc.roda_um = roda_um
rc.confere_config = confere_config
rc.gerar_resumo = gerar_resumo


def main(argv=None):
    argv, k = extrai_janela(sys.argv[1:] if argv is None else argv)
    JANELA[0] = k
    return rc.main(argv)


if __name__ == "__main__":
    sys.exit(main())

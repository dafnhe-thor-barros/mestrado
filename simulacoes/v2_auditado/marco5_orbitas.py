"""marco5_orbitas.py -- WRAPPER do marco5_tems.py (nao edita o arquivo congelado). Mesmo CLI; so acrescenta ao resultado de cada execucao
as listas por orbita (orb_soc_end, orb_ti_Wh, orb_tops, orb_tj_max_C) para calcular janelas de PUE periodico depois (janela_pue.py).

Como funciona: troca marco5_tems.run_scenario por uma versao que forca series=True, calcula as listas por orbita, e remove a serie completa do
resultado se o chamador nao a pediu (--save_traces continua funcionando). A simulacao em si e a mesma chamada ao run_scenario original.
VERIFIQUE no servidor (comando no roteiro): os campos originais (tj_max_C, tops_frac_last, ...) devem ser IDENTICOS aos do marco5_tems.py puro.
Uso: python marco5_orbitas.py <os mesmos argumentos do marco5_tems.py>
"""
import os
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

import marco5_tems as m
import janela_pue as jp

_orig = m.run_scenario


def run_scenario_orb(*args, **kw):
    quer_serie = bool(kw.get("series", False))
    kw["series"] = True
    r = _orig(*args, **kw)
    c = kw.get("c") or m.PROJETO_1P0
    over = kw.get("cfg_over") or {}
    soc0 = float(over.get("soc0", c["soc_max"]))
    jp.enriquece(r, m.T_ORB, c["f_max"], soc0, int(kw.get("n_orb", 3)), m.T0K)
    if not quer_serie:
        r.pop("series", None)
    return r


m.run_scenario = run_scenario_orb

if __name__ == "__main__":
    m.main()

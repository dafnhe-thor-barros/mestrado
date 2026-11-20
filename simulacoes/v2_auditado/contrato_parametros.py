"""contrato_parametros.py -- Marco 5: contrato de parametros (DVFS-only).
Fonte unica dos parametros da malha fechada. O hash SHA256 do contrato e gravado em todo JSON de resultado;
mudou parametro => mudou hash => resultados antigos nao valem.
Etiquetas: [USUARIO] confirmado localmente | [HIPOTESE DECLARADA] sem fonte (sensibilidade obrigatoria) | [CALCULADO] | [LEGADO] valor do legado, sem fonte.
Autocontido (stdlib). Nao importa torch.
"""
import hashlib
import json

A_PL_PROJETO = 49.89754          # m2, area de projeto do Marco 3 (38.176,18 W TI, Tj = 80 C) [USUARIO]

CONTRATO_BASE = {
    "eta_buck": 0.91, "eta_dab": 0.9655, "eta_bd": 0.97, "eta_ch": 0.995, "eta_dis": 0.97,   # [HIPOTESE DECLARADA] (Marco 2b/2c)
    "p_aux": 1000.0, "q_aux": 903.2, "q_oisl": 343.2, "p_pv": 50000.0,                          # [LEGADO 4_1] / decisao do usuario
    "fade": 0.20, "n_ser": 108, "e_cell": 16.2, "soc_min": 0.70, "soc_max": 0.95,               # BESS (Marco 2c)
    "n_mod": 100, "tj_lim_c": 80.0, "f_min": 300.0, "f_max": 1410.0, "tops_mod": 624.0,         # Marco 2d / limite de controle
    "w_s": 300.0, "dt_ctrl": 60.0, "dt_sim": 5.0,                                               # horizonte W e periodos [HIPOTESE DECLARADA]
    "bat_losses_in_qrest": True,                                                                 # default aprovado: perdas eta_ch/eta_dis em q_rest [HIPOTESE DECLARADA]
}
LEGADO_2B_0P90 = dict(CONTRATO_BASE, nome="LEGADO_2B_0P90", eta_pol=0.90)   # reproduz as regressoes confirmadas localmente
PROJETO_1P0 = dict(CONTRATO_BASE, nome="PROJETO_1P0", eta_pol=1.0)           # decisao do usuario: eta_POL,ext = 1.0

M_PRED_K = 0.72                  # [CALCULADO] pior subestimacao do pico na auditoria LHS (300 casos); valor RELATADO, nao recalculado aqui
M_OOD_MAX_K = 3.0                # [HIPOTESE DECLARADA] alargamento linear da margem entre os cortes de OOD
OOD_CORTES = (0.15, 0.35)        # [HIPOTESE DECLARADA]
DTIP_PRIOR, DTIP_RANGE = 7.5, (0.0, 15.0)
SENSOR_MARGEM_K = {"ideal": 0.0, "calibrado": 0.5, "sem_calibracao": 1.5}   # vies de T_j: [HIPOTESE DECLARADA]

CHECKPOINTS_SHA256 = {   # commit a0f34c6; conferir com simulacoes/v2_auditado/CHECKSUMS.txt (commit b9339ee)
    "v2a.pt": "92B770AA99FCDDE99AE464D00376567D9425162D18D0DB66B96DC64FABC02513",
    "v2b.pt": "C5E8714ACBA3F0002AA898CF50A5D808DB60542E01BE9A6197305639EF7C8529",
    "v2c.pt": "05E32E16C6C0C7FB7F7486E418AEF4BCE7370C5F3D865726D00B4C2E7C29FF98"}   # v2c = congelado PROVISORIAMENTE (aguarda Colab)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest().upper()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest().upper()


def contract_hash(c):
    return sha256_bytes(json.dumps(c, sort_keys=True, default=str).encode("utf-8"))


def verify_checkpoint(path, nome=None):
    """Aborta se o arquivo nao for exatamente o checkpoint travado no contrato."""
    nome = nome or str(path).replace("\\", "/").split("/")[-1]
    exp = CHECKPOINTS_SHA256.get(nome)
    if exp is None:
        raise ValueError("checkpoint sem hash no contrato: %s" % nome)
    got = sha256_file(path)
    if got != exp.upper():
        raise RuntimeError("hash divergente para %s: %s != %s" % (nome, got, exp))
    return got


def sensor_cfg(nome, **over):
    """Sensores: T_j (diodo/NVML) e T_raiz (Pt100 4 fios). Sinal -1 = leitura abaixo do real (lado inseguro)."""
    base = {"ideal": dict(sigma_tj=0.0, sigma_rtd=0.0, bias_tj=0.0, bias_rtd=0.0),
            "calibrado": dict(sigma_tj=0.5, sigma_rtd=0.05, bias_tj=0.5, bias_rtd=0.25),
            "sem_calibracao": dict(sigma_tj=0.5, sigma_rtd=0.05, bias_tj=1.5, bias_rtd=0.25)}[nome]
    d = dict(base, bias_sign_tj=-1.0, bias_sign_rtd=-1.0, quant_mode=None, quant_step=1.0, nome=nome)
    d.update(over)
    return d

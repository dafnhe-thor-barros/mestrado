"""contrato_parametros.py -- Marco 5: contrato de parametros (DVFS-only).
Fonte unica dos parametros da malha fechada. O hash SHA256 do contrato e gravado em todo JSON de resultado;
mudou parametro => mudou hash => resultados antigos nao valem.
Etiquetas: [USUARIO] confirmado localmente | [HIPOTESE DECLARADA] sem fonte (sensibilidade obrigatoria) | [CALCULADO] | [LEGADO] valor do legado, sem fonte.
Autocontido (stdlib + modelo_termico.py, para a impressao digital da planta). Nao importa torch.
AUDITORIA 22/11/2026: o hash agora cobre TAMBEM (a) a planta termica (PAR, COATINGS, IA, BUS, BESS e constantes de modelo_termico, ex.: c_area) e
(b) as constantes de controle (M_PRED_K, M_OOD_MAX_K, OOD_CORTES, DTIP_*, SENSOR_MARGEM_K, hashes de checkpoint). Antes, trocar c_area (4050 -> 4860)
NAO mudava o hash e checkpoints/sintonias antigos continuavam aceitos.
"""
import hashlib
import json
import os
import warnings

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

M_PRED_K = 1.09                  # [CALCULADO] pior subestimacao do pico da PINN v2c.pt na auditoria LHS300 com c_area=4860 (-1.0865 C -> 1.09 K)
M_PRED_K_RC3 = 0.33              # [CALCULADO] pior subestimacao do pico do RC3 na auditoria LHS300 com c_area=4860 (-0.3246 C -> 0.33 K)
M_OOD_MAX_K = 3.0                # [HIPOTESE DECLARADA] alargamento linear da margem entre os cortes de OOD
OOD_CORTES = (0.15, 0.35)        # [CALCULADO/ALINHADO] cortes em norma L2 (mesma metrica de auditoria_marco4_robusta.box_distance_v2)
DTIP_PRIOR, DTIP_RANGE = 7.5, (0.0, 15.0)
SENSOR_MARGEM_K = {"ideal": 0.0, "calibrado": 0.5, "sem_calibracao": 1.5}   # vies de T_j: [HIPOTESE DECLARADA]

CHECKPOINTS_SHA256 = {   # commit a0f34c6; conferir com simulacoes/v2_auditado/CHECKSUMS.txt (commit b9339ee)
    "v2a.pt": "92B770AA99FCDDE99AE464D00376567D9425162D18D0DB66B96DC64FABC02513",
    "v2b.pt": "C5E8714ACBA3F0002AA898CF50A5D808DB60542E01BE9A6197305639EF7C8529",
    "v2c.pt": "934DB51B4156BBBA410505F9C89056724ABB786DE9E234408DE8331048767C60"}   # v2c retreinado com c_area = 4860 (planta_sha256 gravado)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest().upper()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest().upper()


def plant_fingerprint():
    """SHA256 dos parametros FISICOS da planta (modelo_termico): PAR, COATINGS, IA, BUS, BESS e constantes orbitais/ambientais."""
    import modelo_termico as mt
    d = {k: getattr(mt, k) for k in ("PAR", "COATINGS", "IA", "BUS", "BESS", "AUX_HEAT_FRAC", "Q_OISL", "T_ORB", "T_ECL", "BETA_DEG", "F_VIEW_90", "G_SUN", "ALBEDO", "Q_IR", "SIGMA")}
    return sha256_bytes(json.dumps(d, sort_keys=True, default=str).encode("utf-8"))


def control_fingerprint():
    d = {"M_PRED_K": M_PRED_K, "M_PRED_K_RC3": M_PRED_K_RC3, "M_OOD_MAX_K": M_OOD_MAX_K, "OOD_CORTES": OOD_CORTES, "DTIP_PRIOR": DTIP_PRIOR, "DTIP_RANGE": DTIP_RANGE,
         "SENSOR_MARGEM_K": SENSOR_MARGEM_K, "CHECKPOINTS_SHA256": CHECKPOINTS_SHA256}
    return sha256_bytes(json.dumps(d, sort_keys=True, default=str).encode("utf-8"))


def contract_hash(c):
    """Hash do contrato + planta + constantes de controle (ver docstring do modulo)."""
    return sha256_bytes(json.dumps({"contrato": c, "planta": plant_fingerprint(), "controle": control_fingerprint()}, sort_keys=True, default=str).encode("utf-8"))


# Impressao digital da planta COM A QUAL cada checkpoint foi treinado. None = nao registrada (treino anterior a auditoria: c_area = 4050 presumido).
# Depois de retreinar com c_area = 4860, registre aqui o valor impresso por modelo_pinn.py (campo "planta_sha256" do .pt) e o novo SHA256 do arquivo.
CKPT_PLANTA_SHA256 = {"v2a.pt": None, "v2b.pt": None, "v2c.pt": "3A4B013FEC539B08D056A5F32F0663900E3A87D628E28E7C3BAF022D4F40AB5D"}


def verify_checkpoint(path, nome=None):
    """Aborta se o arquivo nao for exatamente o checkpoint travado no contrato. Alem disso confere se foi treinado na planta ATUAL:
    divergencia => UserWarning; com M5_STRICT=1 (definido pelo orquestrador em --nova) => RuntimeError."""
    nome = nome or str(path).replace("\\", "/").split("/")[-1]
    exp = CHECKPOINTS_SHA256.get(nome)
    if exp is None:
        raise ValueError("checkpoint sem hash no contrato: %s" % nome)
    got = sha256_file(path)
    if got != exp.upper():
        raise RuntimeError("hash divergente para %s: %s != %s" % (nome, got, exp))
    if CKPT_PLANTA_SHA256.get(nome) != plant_fingerprint():
        msg = "checkpoint %s nao registra a planta atual (treinado com outra c_area/PAR?): retreine e atualize CKPT_PLANTA_SHA256" % nome
        if os.environ.get("M5_STRICT") == "1":
            raise RuntimeError(msg)
        warnings.warn(msg, UserWarning)
    return got


def sensor_cfg(nome, **over):
    """Sensores: T_j (diodo/NVML) e T_raiz (Pt100 4 fios). Sinal -1 = leitura abaixo do real (lado inseguro)."""
    base = {"ideal": dict(sigma_tj=0.0, sigma_rtd=0.0, bias_tj=0.0, bias_rtd=0.0),
            "calibrado": dict(sigma_tj=0.5, sigma_rtd=0.05, bias_tj=0.5, bias_rtd=0.25),
            "sem_calibracao": dict(sigma_tj=0.5, sigma_rtd=0.05, bias_tj=1.5, bias_rtd=0.25)}[nome]
    d = dict(base, bias_sign_tj=-1.0, bias_sign_rtd=-1.0, quant_mode=None, quant_step=1.0, nome=nome)
    d.update(over)
    return d

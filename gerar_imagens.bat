#!/bin/bash
set -e
echo "1. A executar simulacoes no Octave..."
cd simulacoes
octave-cli --eval "simulacao_3_1_barramento_hvdc"
octave-cli --eval "simulacao_3_2_gerador_fotovoltaico_am0"
octave-cli --eval "simulacao_3_3_bateria_dp_aukf"
octave-cli --eval "simulacao_3_4_conversores_gan_middlebrook"
octave-cli --eval "simulacao_4_1_fluxos_orbitais"
octave-cli --eval "simulacao_4_6_fadiga_bga_protecao"

echo "2. A compor figuras finais em Python..."
python gerar_figuras_compostas_capitulo3.py
python gerar_figuras_compostas_capitulo4.py

echo "Processo concluido: todas as figuras foram atualizadas em figuras/."
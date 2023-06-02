#!/bin/bash
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../config.sh
voc=$1
parameters="${path}/data/stage3_processed/Renewal_aug_params/Parameters_${voc}.csv"
covid="${path}/data/stage3/epi_owid.csv"
waiting="15"
lineages="${path}/data/stage3_processed/sequences/${voc}_sequences_fraction"
seeds="${path}/data/stage3_processed/seeds" 
output="${path}/results/Introductions/${voc}/${voc}"

#this script generates 'Introductions' files containing the number of daily arrivals in the target country
#from the source country who are infected with the variant of concern. It does this by comparing flight fluxes,
#estimated prevalence, share of lineage sequences in the source country
python3 Work_with_data.py -p ${parameters} -o "${output}" -c ${covid} -l ${lineages} -a ${seeds} -w ${waiting} 

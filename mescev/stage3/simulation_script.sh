#!/bash/bin
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../../config.sh
voc=$1
outfile="${path}/results/Introductions/${voc}/${voc}" #outfile prefix
arrival="${path}/data/stage3_processed/seeds" #arrivals from source (modify in the future)
parameters="${path}/data/stage3_processed/Renewal_aug_params/Parameters_${voc}" #File with parameters
folder_results="${path}/results/Projections/Renewal_base/${voc}/" #Store results
output_name="${voc}_epidemics"
line_start=1  #line of parameters file to start from
line_end=0    #ending line (0 = read until the end)

#Runs the renewal process
Rscript Renewal_process.R "${parameters}.csv" "${outfile}" "${folder_results}" ${line_start} ${line_end} "${output_name}"
#Rscript Renewal_process_augmented.R  "${parameters}.csv" "${outfile}" "${folder_results}" ${line_start} ${line_end} "${output_name}"
#python Consolidate_output.py -s "${folder_result${path}/epidemic_results_${date_start}_${date_end}" -m 1 4 -p ../Data/Parameters.txt -l 1 2

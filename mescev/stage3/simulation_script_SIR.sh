#!/bash/bin
cd "$(dirname "${BASH_SOURCE[0]}")"
source ../config.sh 
voc=$1
parameters="${path}/data/processed/SIR_params_Phylo/Parameters_${voc}.csv" #File with parameters
covid_file="${path}/Data/epi_owid.csv"
lineages="${path}/Data/sequences_fractions/${voc}_sequences_fraction"
introductions="${path}/Results/Introductions/${voc}/${voc}" #outfile prefix
folder_results="${path}/Results/Projections/SIR/${voc}/${voc}_epidemics" #Store results
line_start=1  #line of parameters file to start from
line_end=0    #ending line (0 = read until the end)

#Runs the renewal process
python multipop_SIR.py -p ${parameters} -c ${covid_file} -l ${lineages} -i ${introductions} -o ${folder_results}
#python Consolidate_output.py -s "${folder_result${path}/epidemic_results_${date_start}_${date_end}" -m 1 4 -p ../Data/Parameters.txt -l 1 2

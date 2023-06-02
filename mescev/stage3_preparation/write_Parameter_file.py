# write Parameter file for simulation_script.sh

import math
import pandas as pd
import sys
import os
from datetime import datetime, date, timedelta
from pathlib import Path


# import numpy as np
# import matplotlib.pyplot as plt

# example use: python3 write_parameter_file.py 'Alpha' 20 '../data/stage3/' '../data/stage3_processed/' '../data/stage3_processed/Renewal_aug_params/'

def find_undercount(country_target, date_start, data_folder):
    if (date_start>="2020-01-01") & (date_start<"2020-07-01"):
        file_name = '2020-1st'
    elif (date_start>="2020-07-01") & (date_start<"2021-01-01"):
        file_name = '2020-2nd'
    elif (date_start>="2021-01-01") & (date_start<"2021-07-01"):
        file_name = '2021-1st'
    elif (date_start>="2021-07-01") & (date_start<"2022-01-01"):
        file_name = '2021-2nd'
    elif (date_start>="2022-01-01") & (date_start<"2022-07-01"):
        file_name = '2022-1st'
    else:
        print('date not in range')

    df_ihme = pd.read_csv( data_folder + 'ihme_cases_aggregated_by_country/ihme_cases_aggregated_by_country_' + file_name + '.csv', sep = ',', engine = 'python')
    lower = df_ihme.loc[df_ihme['ISO_A3']==country_target, 'under.lower'].values[0]
    upper = df_ihme.loc[df_ihme['ISO_A3']==country_target, 'under.upper'].values[0]
    undercounting_multipliers = ["{:.2f}".format(lower), "{:.2f}".format(upper)]
    return undercounting_multipliers

def get_date(decimal_date):
    day_of_the_year, year = math.modf(decimal_date)
    if int(year)%4==0:
        a = datetime(int(year), 1, 1) + timedelta(round(day_of_the_year*365.242))
    else:
        a = datetime(int(year), 1, 1) + timedelta(round(day_of_the_year*365.242)-1)
    return a.strftime('%Y-%m-%d')

def get_Rt_from_growth_rate(gr):
    media = 5.2 # gen time gamma dist
    std = 2.14
    shape = media**2/std**2
    rate = media/std**2
    R1 = (((gr/365.242) + rate)**shape)/rate**shape
    R2 = media*(gr/365.242) +1
    print(shape)
    print(rate)
    return R1


def main(argv):
    variant_name = argv[0] # example: 'Alpha'
    sample=int(argv[1]) # example: 20
    data_folder = argv[2] # example: './stage3/'
    processed_folder = argv[3] # example: './stage3_processed/'
    output_folder = argv[4] # example: './stage3_processed/Renewal_aug_params/'

    # some variable that depends on variants
    if variant_name == 'Alpha':
        variant_code = 'B.1.1.7' # check
        lineage = variant_name + ' ' + variant_code
        country_source = 'GBR'
        country_source_long = 'England'
        Ve = '0' # Vaccine efficacy
        Ne = '0' # Natural infection efficacy
    elif variant_name == 'Delta':
        variant_code = 'B.1.617.2'
        lineage = variant_name + ' ' + variant_code
        country_source = 'IND'
        country_source_long = 'India'
        Ve = '0.5'
        Ne = '1'
    elif variant_name == 'BA1-Omicron':
        variant_code = 'B.1.1.529.1' # Or 21K
        lineage = 'Omicron BA.1'
        country_source = 'ZAF'
        country_source_long = 'South Africa'
        Ve = '0.5'
        Ne = '0.56'

    elif variant_name == 'BA2-Omicron':
        variant_code = 'B.1.1.529.2' # Or 21L # Earliest Date: 2021-10-22 from covariants.org
        lineage = 'Omicron BA.2'
        country_source = 'ZAF'
        country_source_long = 'South Africa'
        Ve = '0.5' # no differences with BA1
        Ne = '0.56'

    elif variant_name == 'BA5-Omicron':
        variant_code = 'B.1.1.529.5' # or 22B  # Earliest Date: 2021-11-15 from covariants.org
        country_source = 'ZAF'
        lineage = 'Omicron BA.5'
        country_source_long = 'South Africa'
        Ve = "{:.2f}".format(0.5/4.2)
        Ne = "{:.2f}".format(0.56/4.2)
    else:
        print('Error: variant not found')

    # load "our world in data" data
    file_name = 'epi_owid'
    df_owid = pd.read_csv(data_folder + file_name + '.csv', sep = ',', engine = 'python')

    # load Philo data
    file_name='BEAST_epi_parameters' # PHYLO_finalResults
    phylo_file = Path(processed_folder).parent / f'stage1/{file_name}.csv'   # = mescev/data/stage1/BEAST_epi_parameters.csv
    df_phylo = pd.read_csv(phylo_file, sep = ',', engine = 'python')

    # find start date and growth date
    decimal_date = df_phylo.loc[(df_phylo['Parameter']=='tMRCA')&(df_phylo['Lineage']==lineage)&(df_phylo['Sample']==sample),'Mean'].values[0]
    growth_rate_mean = df_phylo.loc[(df_phylo['Parameter']=='growthRate')&(df_phylo['Lineage']==lineage)&(df_phylo['Sample']==sample),'Mean'].values[0]
    growth_rate_min = df_phylo.loc[(df_phylo['Parameter']=='growthRate')&(df_phylo['Lineage']==lineage)&(df_phylo['Sample']==sample),'LowHPD'].values[0]
    growth_rate_max = df_phylo.loc[(df_phylo['Parameter']=='growthRate')&(df_phylo['Lineage']==lineage)&(df_phylo['Sample']==sample),'UpHPD'].values[0]
    date_start = get_date(decimal_date)
    Rt_mean=get_Rt_from_growth_rate(growth_rate_mean)
    Rt_min=get_Rt_from_growth_rate(growth_rate_min)
    Rt_max=get_Rt_from_growth_rate(growth_rate_max)
    Rt_values = [Rt_min, Rt_mean, Rt_max]

    # load list of countries
    list_countries = [x[:3] for x in os.listdir(processed_folder +'sequences/'+ variant_name + '_sequences_fraction/')]
    list_countries.sort()

    # write
    undercounting_multipliers_source = find_undercount(country_source, date_start,data_folder)

    line = []

    for country_target in list_countries:
        print(country_target)
        if country_target!=country_source:

            # date_end: date_start_target + N days
            #--------------------------------------
            try:

                undercounting_multipliers_target = find_undercount(country_target, date_start,data_folder)

                # end sim when sequences are found in the target country
                file_name_data = country_target + '_sequences_fraction'
                df_seq_all = pd.read_csv(processed_folder +'sequences/'+ variant_name + '_sequences_fraction/'  + file_name_data + '.csv', sep = ' ', engine = 'python')
                df_seq_all.columns=['week', 'sequences_fraction', 'total_sequences']

                date_end = df_seq_all.loc[df_seq_all['sequences_fraction']>0, 'week'].tolist()[2] # third point
                print('date_end=' + date_end )


                # Evaluate immunity: if t1 is the time of first reported case in country then :
                # Removed = undercounting_multiplier_target * Cum Cases[t1-6 months] + Ve *(FullyVax[t1])
                #-----
                date_start_estimation_recov = pd.to_datetime(date_start) - timedelta(days = 180)
                date_start_estimation_recov = date_start_estimation_recov.strftime('%Y-%m-%d')
                Cum_cases = df_owid.loc[(df_owid['iso_code']==country_target) & (df_owid['date']>=date_start_estimation_recov) & (df_owid['date']<=date_start), 'new_cases'].sum()

                PopFullyVax =  df_owid.loc[(df_owid['iso_code']==country_target) & (df_owid['date']==date_start),'people_fully_vaccinated'].values[0]
                # if 'people_fully_vaccinated' is nan, put 0
                if math.isnan(PopFullyVax):
                    PopFullyVax=0
                pop = df_owid.loc[(df_owid['iso_code']==country_target) & (df_owid['date']==date_start),'population'].values[0]


                # Rt
                #-----
                for Rt in Rt_values:
                    Rt = round(Rt, 2)

                    for undercounting_multiplier_source in undercounting_multipliers_source:
                        for undercounting_multiplier_target in undercounting_multipliers_target:
                            line.append([variant_code, country_source , country_target, date_start, date_end, undercounting_multiplier_source, undercounting_multiplier_target, Ne, Ve, str(Cum_cases), str(PopFullyVax), str(pop), str(Rt)])

            except IndexError:
                print(country_target)

    df_par = pd.DataFrame(line, columns=['variant','source_country','target_country','date_start','date_end','undercounting_multiplier_source','undercounting_multiplier_target','Ne', 'Ve','Cum_cases', 'PopFullyVax', 'pop', 'Rt'])

    # Save
    do_save=True
    if do_save:
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        df_par.to_csv(output_folder + 'Parameters_' + variant_name + '.csv', index=None, sep=',')



"""
RUN MAIN
"""
if __name__== "__main__":
    main(sys.argv[1:])

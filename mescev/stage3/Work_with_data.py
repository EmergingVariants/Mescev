#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 19:14:19 2022

@author: fra
"""

#Simple routine that takes a date and two countries as inputs and generates number of infected
#seeds per day from country one to country two for n days starting from that
#period

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import argparse
import os

#arguments="-p ../Data/Parameters_v1.txt -o ../Data/Parameters_v1.txt -w 15 -c ../Data/Data_ENG_IT/epi_ovid.csv -l ../Data/Data_ENG_IT/lineages_by_ltla_and_week.tsv -a ../Data/Data_ENG_IT/seeds"
'''
parameters = "-p /home/fra/Projects/EmergingVariants/Omicron/Params/Renewal_aug_params/Parameters_BA5-Omicron.csv "
covid ="-c /home/fra/Projects/EmergingVariants/Omicron/Data/epi_ovid.csv "
waiting="-w 15 "
lineages="-l /home/fra/Projects/EmergingVariants/Omicron/Data/sequences_fraction/BA2-Omicron_sequences_fraction "
seeds="-a /home/fra/Projects/EmergingVariants/Omicron/Data/seeds "
output="-o /home/fra/Projects/EmergingVariants/Omicron/Results/Introductions/Data"
arguments = parameters+covid+waiting+lineages+seeds+output
runfile('Work_with_data.py', args =arguments )
'''

variant_dictionary={'Alpha':'B.1.1.7', 'Beta':'B.1.351', 'Delta':'B.1.617.2', 'Gamma':'P.1', 'Omicron':'B.1.1.529'}
def get_keys_from_value(d, val):
    return [k for k, v in d.items() if v == val]

if __name__=="__main__":

    # Process the input argument
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--parameters_file", type = str, required = True, \
       help = "Parameters filename")

    parser.add_argument("-o", "--outfilename", type = str, required = False, \
       help = "Output filename prefix (may include a path)", default = 'output')

    parser.add_argument("-w", "--window_size", type = int, required = False, \
       help = "moving windows of n days for cumulative incidence (default=15)", default = 15)

    parser.add_argument("-c", "--coviddata", type = str, required = True, \
       help = "covid data file formatted as in owid data")

    parser.add_argument("-l", "--lineages", type = str, required = True, \
       help = "lineage share through time")

    parser.add_argument("-a", "--arrival", type = str, required = True, \
       help = "arrival times folder (with subfolders for months)")


    parser.usage = parser.format_help()
    args = parser.parse_args()


    #ALL THIS STUFF NEEDS TO BE CHANGED PROGRAMMATICALLY IN THE NEAR FUTURE
    paramter_file = args.parameters_file
    covid_data = args.coviddata      #OWID epi-related data
    lineages = args.lineages         #Lineage files
    arrival_folder = args.arrival    #Arriving seeds from src to country
    outfilename = args.outfilename   #Output file
    windows_size = args.window_size #For cumulative incidence


    parameters=pd.read_csv(paramter_file,sep=',')

    #We want to minimize the number of files to be created at this step, so
    #we look at the unique combinations of parameters in terms of data

    countries_source_unique = np.unique(parameters.source_country)
    countries_target_unique = np.unique(parameters.target_country)


    owid_file = pd.read_csv(covid_data, sep=',')
    owid_file.date = pd.to_datetime(owid_file.date)
    for country_source in countries_source_unique:
        for country_target in countries_target_unique:
            print(country_target)
            starting_date_unique = np.unique(parameters.query('source_country == @country_source and target_country==@country_target').date_start)
            ending_date_unique = np.unique(parameters.query('source_country == @country_source and target_country==@country_target').date_end)
            variants_unique = np.unique(parameters.query('source_country == @country_source and target_country==@country_target').variant)
            incidence_multiplier_in_unique = np.unique(parameters.query('source_country == @country_source and target_country==@country_target').undercounting_multiplier_source)
            incidence_multiplier_out_unique = np.unique(parameters.query('source_country == @country_source and target_country==@country_target').undercounting_multiplier_target)
            params=parameters.query('source_country==@country_source and target_country==@country_target')

            for starting_date in starting_date_unique:
                for ending_date in ending_date_unique:
                    for variant in variants_unique:
                        for incidence_multiplier_in in incidence_multiplier_in_unique:
                            for incidence_multiplier_out in incidence_multiplier_out_unique:

                                #How many days we want from the starting date (no more than a threshold - 90 days)
                                #end_date=min(datetime.fromisoformat(ending_date),datetime.fromisoformat(starting_date)+timedelta(days=90))
                                #if (os.path.exists("{}_{}_{}_{}_{}_{}_{}_{}.csv".format(outfilename,country_source,country_target,starting_date,end_date,variant,incidence_multiplier_in,incidence_multiplier_out)
   # )):
                                   # continue
                                end_date = datetime.fromisoformat(ending_date)
                                #end_date=datetime.fromisoformat(starting_date) + timedelta(days=90)
                                if end_date < datetime.fromisoformat(starting_date) + timedelta(days=28):
                                    end_date = datetime.fromisoformat(starting_date) + timedelta(days=28)
                                number_of_days = (end_date-datetime.fromisoformat(starting_date)).days
                                #Consider data from the starting date onwards
                                target_epi_data_in = owid_file.query('iso_code == @country_source and date > @starting_date' )
                                population = target_epi_data_in.population.to_list()[0]

                                #For population out take past six months
                                sixmonthsearly_begin=datetime.fromisoformat(starting_date)+timedelta(days=-180)
                                target_epi_data_out = owid_file.query('iso_code == @country_target and date >=@sixmonthsearly_begin and date <@end_date' )

                                target_epi_data_out_active_recovered=target_epi_data_out['total_cases']-target_epi_data_out['total_cases'].shift(180)
                                target_epi_data_out_active_recovered = np.array(target_epi_data_out_active_recovered[180:].to_list())
                                target_epi_data_out_active_recovered[target_epi_data_out_active_recovered<0]=0
                                if(len(target_epi_data_out_active_recovered)<number_of_days):
                                    for i in range( number_of_days - len(target_epi_data_out_active_recovered)):
                                        target_epi_data_out_active_recovered = np.insert(target_epi_data_out_active_recovered,0,0)

                                target_epi_data_out_vax=owid_file.query('iso_code == @country_target and date >=@sixmonthsearly_begin and date <@end_date' ).fillna(0)
                                target_epi_data_out_vax=target_epi_data_out_vax['people_fully_vaccinated']-target_epi_data_out_vax['people_fully_vaccinated'].shift(180)
                                target_epi_data_out_vax = np.array(target_epi_data_out_vax[180:].to_list())
                                target_epi_data_out_vax[target_epi_data_out_vax<0]=0
                                if(len(target_epi_data_out_vax)<number_of_days):
                                    for i in range( number_of_days - len(target_epi_data_out_vax)):
                                        target_epi_data_out_vax = np.insert(target_epi_data_out_vax,0,0)

                                target_epi_data_out = owid_file.query('iso_code == @country_target and date >@starting_date' )
                                population_TRG = target_epi_data_out.population.to_list()[0]
                                target_epi_data_out_active_recovered = target_epi_data_out_active_recovered/population_TRG*incidence_multiplier_out
                                target_epi_data_out_active_recovered[target_epi_data_out_active_recovered>=1]=1
                                target_epi_data_out_vax = target_epi_data_out_vax/population_TRG

                                #Incidence (cut at n days)
                                Incidence = target_epi_data_in['new_cases'].fillna(0).to_list()[:number_of_days]


                                #Use it only for ALPHA_INCIDENCE, WAIT FOR SEB for the rest!!

                                lineages_through_time = pd.read_csv(lineages+"/{}_sequences_fraction.csv".format(country_source), sep=' ')
                                lineages_through_time['WeekEndDate'] =  pd.to_datetime(lineages_through_time.week)

                                #These are weekly updated absolute numbers of cases by variant, we need to get
                                #to the day to day basis

                                min_date= lineages_through_time.loc[(lineages_through_time.WeekEndDate
                                                  -pd.to_datetime(starting_date)).abs().idxmin(),'WeekEndDate']

                                max_date = min_date + timedelta(days=number_of_days+7)
                                #Select only rows after the date of interest
                                lineages_through_time_after_start = lineages_through_time.query('WeekEndDate>@min_date and WeekEndDate<@max_date').copy()

                                '''
                                lineages_through_time_after_start['variant'] = [True if variant  in x else False for x in lineages_through_time_after_start['Lineage']]
                                number_of_sequences = lineages_through_time_after_start.groupby(['WeekEndDate'])['Count'].sum()

                                #number_of_sequences_that_are_alpha
                                number_of_sequences_variant =  lineages_through_time_after_start[ lineages_through_time_after_start['variant']==True]

                                number_of_sequences_variant =  number_of_sequences_variant.groupby(['WeekEndDate'])['Count'].sum()



                                proportion = number_of_sequences_variant/number_of_sequences
                                np.nan_to_num(proportion, copy=False, nan=0.0)

                                '''
                                proportion = lineages_through_time_after_start.sequences_fraction

                                #Now we need to multiply daily incidence by number of sequences
                                #put weeks into days
                                proportion_in_days = np.repeat(np.array(proportion.to_list()),7)[:number_of_days]

                                incidence_by_variant = proportion_in_days*Incidence

                                #Consider cumulative sum over 15 days

                                cumulative_incidence = np.cumsum(incidence_by_variant)*incidence_multiplier_in

                                cumulative_incidence[windows_size:] -= cumulative_incidence[:len(cumulative_incidence)-windows_size]


                                dates=pd.date_range(start=pd.to_datetime(starting_date),  end=pd.to_datetime(max_date))[:number_of_days]


                                months = dates.strftime('%Y-%m')
                                #Take daily passengers from Seb files
                                daily_passengers = np.zeros(len(dates))

                                try:
                                    for index,month in enumerate(months):
                                        #variant_iso = get_keys_from_value(variant_dictionary,variant)[0]
                                        file_departures = pd.read_csv("{}/{}/daily_passengers_from_{}_month_{}.csv".format(arrival_folder,month, country_source,month), sep=' ')

                                        daily_passengers[index] = file_departures.query('iso_a3_target == @country_target').daily_seeds.to_list()[0]
                                except:
                                    print("no month:{}\n".format(month))
                                multiplier_in=incidence_multiplier_in
                                multiplier_out=incidence_multiplier_out

                                out_name="{}_{}_{}_{}_{}_{}_{}.csv".format(outfilename,country_source,country_target,starting_date,ending_date,multiplier_in,multiplier_out)
                                dataframe_to_save= {'days':dates,'cumulative_incidence_SRC':cumulative_incidence*daily_passengers/population,
                                                    'cumulative_incidence_TRG':target_epi_data_out_active_recovered,'cumulative_incidence_vax': target_epi_data_out_vax}
                                df = pd.DataFrame(dataframe_to_save)
                                df.to_csv(out_name, sep=',', index=True)







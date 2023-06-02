#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 31 14:33:41 2022

@author: fra
"""

import numpy as np
from scipy import integrate
import pandas as pd
from datetime import datetime, timedelta
import argparse


# import os
# import matplotlib.pyplot as plt


parameters = "-p /home/fra/Projects/EmergingVariants/Omicron/Params/SIR_params/Parameters_Omicron.csv "
covid ="-c /home/fra/Projects/EmergingVariants/Omicron/Data/epi_ovid.csv "
lineages="-l /home/fra/Projects/EmergingVariants/Omicron/Data/sequences_fractions/Omicron_sequences_fraction "
introductions="-i /home/fra/Projects/EmergingVariants/Omicron/Results/Introductions/Omicron/Omicron "
output="-o /home/fra/Projects/EmergingVariants/Omicron/Results/Projections/SIR/Omicron/Omicron_sir"
arguments = parameters+covid+introductions+lineages+output
#runfile('multipop_SIR.py', args =arguments )


if __name__=="__main__":


    parser = argparse.ArgumentParser()

    parser.add_argument("-p", "--parameters_file", type = str, required = True, \
       help = "Parameters filename")

    parser.add_argument("-c", "--coviddata", type = str, required = True, \
       help = "covid data file formatted as in owid data")

    parser.add_argument("-l", "--lineages", type = str, required = True, \
       help = "lineage share through time")
    parser.add_argument("-i", "--introductions", type = str, required = True, \
       help = "variant introductions from source")

    parser.add_argument("-o", "--outfilename", type = str, required = True, \
       help = "arrival times folder (with subfolders for months)")

    parser.usage = parser.format_help()
    args = parser.parse_args()


    #ALL THIS STUFF NEEDS TO BE CHANGED PROGRAMMATICALLY IN THE NEAR FUTURE
    paramter_file = args.parameters_file
    covid_data = args.coviddata      #OWID epi-related data
    introductions = args.introductions    #Arriving seeds from src to country
    outfilename = args.outfilename   #Output file
    lineages=args.lineages
    parameters=pd.read_csv(paramter_file,sep=',')
    cov_data=pd.read_csv(covid_data)
    #arrivals = pd.read_csv("/home/fra/Projects/EmergingVariants/Omicron/Results/Introductions/Alpha/Alpha_GBR_ITA_2020-09-20_2020-12-19_B.1.1.7_4_4.csv")
    def sir_two_variants(v,t,parameters):
            t_0,t_1,beta_0, beta_1, gamma, alpha,N=parameters
            vdot = np.zeros(len(v))
            lambda_0 = beta_0*v[1]
            lambda_1 = beta_1*v[2]
            vdot[0]  = - v[0]*(lambda_0 + lambda_1)
            #print(lambda_0)
            vdot[1]  = lambda_0*v[0] - gamma*v[1]
            vdot[2]  = lambda_1*v[0] - gamma*v[2] + (1-alpha)*lambda_1*v[3]
            vdot[3] = gamma*v[1] -(1-alpha)*lambda_1*v[3]
            vdot[4] = gamma*v[2]
            return vdot


    def sir_odeint(model, t_f, incond_0, arrivals, parameters,nt = 500):
            incond = incond_0
            sol=np.zeros((int(t_f),5))
            #sol=np.zeros((int(t_f*nt),5))


            for t in np.arange(0,t_f):
                incond[2] = incond[2]+arrivals[t]

                t_day = np.linspace(0,1,nt)
                y= integrate.odeint(model, incond, t_day, args=(parameters,))
                incond=y[-1]
                #sol[int(t*nt):int((t+1)*nt)] = y
                sol[t]=y[-1]
            return sol

    unique_targets = np.unique(parameters.target_country)
    for target in unique_targets:
        new_parameters = parameters[parameters.target_country==target].copy().reset_index(drop=True)
        for index,row in new_parameters.iterrows():

            try:
                Ve = new_parameters.Ve[index]     #Vaccine efficacy
                Ne = new_parameters.Ne[index]    #Natural immunity efficacy
                population = new_parameters['pop'][index]   #Population in target country
                source = new_parameters.source_country[index]
                target = new_parameters.target_country[index]
                multiplier_in = new_parameters.undercounting_multiplier_source[index]
                multiplier_out = new_parameters.undercounting_multiplier_target[index]
                date_start = new_parameters.date_start[index]
                date_end = new_parameters.date_end[index]
                variant = new_parameters.variant[index]

                gamma = 1/parameters.Tg[index]
                seeds=pd.read_csv(lineages+"/"+target+"_sequences_fraction.csv", sep=" ")
                seeds.week=pd.to_datetime(seeds.week)
                seeds = seeds[(seeds.week)>=datetime.fromisoformat(date_start)+timedelta(days=-7)]
                seeds = seeds[seeds.week<=datetime.fromisoformat(date_end)+timedelta(days=7)]

                expanded_seeds=np.repeat(seeds.sequences_fraction,7)
                df_ovid = cov_data.loc[(cov_data['iso_code']==target) & (cov_data['date']>=date_start) &  (cov_data['date']<=date_end)].reset_index(drop=True)

                infected_at_0 = df_ovid.new_cases[0]*multiplier_out
                cum_cases = df_ovid.total_cases[0]*multiplier_out
                beta_0 = population*gamma*df_ovid.reproduction_rate[0]/(population-df_ovid.total_cases[0]*multiplier_out)
                beta_var = new_parameters.beta_alpha[index]

                introduction_name="{}_{}_{}_{}_{}_{}_{}.csv".format(introductions,
                    source,target,date_start, date_end,multiplier_in,multiplier_out)
                #print(introduction_name)

                introductions_read = pd.read_csv(introduction_name)
                #print(introductions_read)
                # beta_0 = 1.3/6.6
                # beta_1 = 1.5*beta_0
                # gamma=1/6.6
                # alpha=0.9
                # parameters = [t_0,t_1,t_2,beta_0, beta_1, gamma, alpha,N]


                number_of_days = (datetime.fromisoformat(date_end)-datetime.fromisoformat(date_start)).days
                parameters_for_equation = [0,number_of_days,beta_0, beta_var*beta_0, gamma, 1,1]
                incond=[(population-cum_cases-infected_at_0)/population,infected_at_0/population,0,cum_cases/population,0]
                arrivals = np.array(introductions_read.cumulative_incidence_SRC.to_list())/population

                y=sir_odeint(sir_two_variants, number_of_days, incond, arrivals, parameters_for_equation,nt = 200)
                Incidence = np.diff(y[:,2],prepend=0)

                datelist = np.array(pd.date_range(datetime.fromisoformat(date_start), periods=number_of_days).tolist())

                epi_dictionary = {'date':datelist, 'I':y[:,2],'Incidence':Incidence*population, 'imported':arrivals*population}
                epidemics = pd.DataFrame.from_dict(epi_dictionary)

                name = "{}_{}_{}_{}_{}_{}_undercounting_SRC_{}_undercounting_TRG_{}_beta_var_{}_Tg_{}_ne_{}_ve_{}.csv".format(outfilename,source,target,date_start,date_end,variant,multiplier_in,multiplier_out,beta_var,parameters.Tg[index],Ne,Ve)
                epidemics.to_csv(name,index=False )
            except:
                print("country of {} does not find some file (check introductions)\n".format(target))
                1/0
                continue

            '''
            plt.plot(y[:,2])
            try:

                #plt.scatter(np.arange(len(y)+1),df_ovid.new_cases*multiplier_out/population)
                cases_to_plot=np.array(df_ovid.new_cases.to_list())[:len(y)+1]
                plt.scatter(np.arange(len(y)+1),expanded_seeds[:len(y)+1]*cases_to_plot*multiplier_out/population)
                plt.title(target)
                plt.xlabel('days')
                plt.ylabel('Incidence')
            except:
                continue
            '''






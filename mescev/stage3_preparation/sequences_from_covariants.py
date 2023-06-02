# -*- coding: utf-8 -*-
"""
Functions to read covariants data json and merge Delta sequences.
Read JSON from Hodcroft's cluster tables

fields in covariants datasets (json read as dicts):
```country_keys = ['week', 'total_sequences', 'cluster_sequences', 'unsmoothed_cluster_sequences', 'unsmoothed_total_sequences']```
"""

import numpy as np
import pandas as pd
import json
import os
import datetime


# from scipy import stats
# import sys
# import networkx as nx
# import scipy as sc
# import matplotlib.pyplot as plt

def generate_directory(directory, variant):
    '''
    Auxiliary function to generate directory where to store data
    '''
    if not os.path.exists(directory):
        os.mkdir(directory)
        print(f"Generated directory to save variant data {variant}. Now saving data ..")
    else:
        print(f"Directory '{directory}' already exists. Now saving data ..")


def generate_mapping(countries, countries_codes):
    '''
    Auxiliary function between covariants and who names used for countries
    '''
    mapping = {}
    for c in countries:

        # correction as in country_keys of covariants data it was stored as USA
        if c == "USA":
            c = "United States"
        value = countries_codes[countries_codes['Country']==c]['Alpha-3 code'].values[0]

        # slicing the correct value
        iso_a3_code = str(value[2:5])
        mapping[c] = iso_a3_code

    return mapping


def generate_ext_lists(weeks, sequences_fractions, total_sequences):
    '''
    (documentation to be improved)
    Function to interpolate covariants data for fractions and weekdates
    from 2-week sampling to 1-week sampling. Interpolation of sequences fraction
    is done by taking average between two dates:
    eg: additional point mid_s_f = (sequences_fractions[i] + sequences_fractions[i+1])/2.0
    ---
    Parameters:
    - weeks [list]
    - sequences_fractions [list]
    - total_sequences [list]
    Returns:
    - weeks_ext [list], sequences_fractions_ext [list], total_sequences_ext [list]
    '''

    week_delta = datetime.timedelta(days=7)
    weeks_ext, sequences_fractions_ext, total_sequences_ext = [],[],[]

    # Weeks addition
    for w in weeks[:-1]:

        weeks_ext.append(w)
        w_datetime = datetime.datetime.strptime(w, "%Y-%m-%d")
        w_datetime_delta = w_datetime + week_delta
        weeks_ext.append(w_datetime_delta.strftime("%Y-%m-%d"))

    # Add last value missing
    weeks_ext.append(weeks[-1])

    # Sequences fraction interpolation
    for i in range(len(sequences_fractions)-1):
        s_f = sequences_fractions[i]
        mid_s_f = (sequences_fractions[i] + sequences_fractions[i+1])/2.0

        sequences_fractions_ext.append(s_f)
        sequences_fractions_ext.append(round(mid_s_f,3))

    # Add last value missing
    sequences_fractions_ext.append(sequences_fractions[-1])

    # Total sequences not interpolated
    for i in range(len(total_sequences)-1):
        total_sequences_ext.append(total_sequences[i])
        total_sequences_ext.append(total_sequences[i])

    # Add last value missing
    total_sequences_ext.append(total_sequences[-1])

    return weeks_ext, sequences_fractions_ext, total_sequences_ext



def save_country_sequences(data, country, mapping, path):

    cluster_fraction = []
    for i in range(len(data[country]['cluster_sequences'])):

        if data[country]['total_sequences'][i] != 0.0:
            fraction = np.round(data[country]['cluster_sequences'][i]/float(data[country]['total_sequences'][i]),3)
            cluster_fraction.append(fraction)
        else:
            cluster_fraction.append(0.0)

    # Extend fields with 1-week interpolation
    weeks_ext, sequences_fractions_ext, total_sequences_ext = generate_ext_lists(
                                                                weeks = list(data[country]['week']),
                                                                sequences_fractions = list(np.array(cluster_fraction)),
                                                                total_sequences = list(data[country]['total_sequences']))

    country_dict = {'week':weeks_ext, 'sequences_fraction':sequences_fractions_ext, 'total_sequences': total_sequences_ext}
    df = pd.DataFrame.from_dict(country_dict)

    # correction as in country_keys of covariants data it was stored as USA
    if country == "USA": country = "United States"
    iso_a3_code = mapping[country]

    # Save
    df.to_csv(path +f'{iso_a3_code}_sequences_fraction.csv', sep=' ', index = None)


####
# ---

variants_df = pd.read_csv('../data/stage3/variants.csv', sep = ',')
variants = dict(zip(variants_df['name'], variants_df['file_covariants']))

# WHO name mapping from covariants codes (omicron currently)
omicron_names_mapping_df = pd.read_csv('../data/stage3/omicron_names_mapping.csv', sep = ',')
omicron_names_mapping = dict(zip(omicron_names_mapping_df['covariants_name'], omicron_names_mapping_df['who_name']))

# Generate sequences fractions
for variant_key in list(variants.keys()):

    file = '../data/stage3/covariants_org/'+variants[variant_key]

    # rename
    if variant_key in list(omicron_names_mapping.keys()):
        file_name = omicron_names_mapping[variant_key]
    else:
        file_name = variant_key

    # path where directory is generated and sequences are stored
    path = f"../data/stage3_processed/sequences/{file_name}_sequences_fraction/"

    # Generate directory to save variant
    generate_directory(path, file_name)

    # open file as dict
    with open(file) as data_file:
        data = json.load(data_file)

    # set of countries with sequencing data in covariants tables
    countries = list(data.keys())

    # countries codes
    countries_codes = pd.read_csv('../data/stage3/countries_codes_and_coordinates.csv', sep = ',')

    # remove certain countries for which I don't have iso-a3-coding (to be solved)
    countries = [c for c in countries if c not in ['Curacao','Bonaire' ,'Sint Maarten','Botswana']]
    mapping = generate_mapping(countries, countries_codes)

    # store sequencing for each country
    for country in countries:
        save_country_sequences(data, country, mapping, path)


# ---
# For Delta aggregate sequencing for different subvariants

delta_keys = ['21A-Delta','21I-Delta','21J-Delta']
path_delta = f"../data/stage3_processed/sequences/Delta_sequences_fraction/"

# Generate directory to save variant
generate_directory(path_delta, 'Delta')

for country in countries:

    if country == "USA": country = "United States"
    iso_a3_code = mapping[country]

    dataframes, fractions = [], []

    for key in delta_keys:
        file =  f"../data/stage3_processed/sequences/{key}_sequences_fraction/{iso_a3_code}_sequences_fraction.csv"
        df = pd.read_csv(file, sep = ' ')

        dataframes.append(df)
        fractions.append(np.array(df['sequences_fraction'].values))

    df_ref = dataframes[0]

    fractions_total = np.zeros(len(fractions[0]))
    for f in fractions:
        fractions_total = fractions_total + f

    # round decimals
    fractions_total = np.around(fractions_total, decimals=4)

    combined_dict = {'week':list(df_ref['week'].values),
                'sequences_fraction':fractions_total,
                'total_sequences':list(df_ref['total_sequences'].values)}

    df_delta = pd.DataFrame.from_dict(combined_dict)

    # Save
    df_delta.to_csv(path_delta +f'{iso_a3_code}_sequences_fraction.csv', sep=' ', index = None)












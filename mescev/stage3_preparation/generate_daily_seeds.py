# -*- coding: utf-8 -*-
"""
Functions to obtain daily seeds from import risk
to be fed to epidemic models
"""
import pandas as pd
import os

# import the IR functions
from importrisk_load_precomputed import precomputed_import_risk_cntr_agg

# get as input the list of months (as a csv) and the set of sources
# and run the code to generate seeds and store them


def generate_directory(directory, month):
    '''
    auxiliary function
    '''

    if not os.path.exists(directory):
        os.mkdir(directory)
        print(f"Generated directory to save seeds for month {month}")
    else:
        print(f"Directory for month {month} already exists. Not generating")


def generate_seeds_from_SOURCE_in_MONTH(SOURCE, month, isoA3_to_isoA2, isoA2_to_isoA3,
                                    PATH):
    '''
    (documentation to be improved)
    function to generate number of seeds from precomputed IR and save in PATH
    ---
    parameters:
    - SOURCE
    - month
    - isoA3_to_isoA2 mapping dictionary
    - isoA2_to_isoA3 mapping dictionary
    - PATH to store seeds as .csv files
    '''

    # get precomputed IR
    out = precomputed_import_risk_cntr_agg(month, isoA3_to_isoA2[SOURCE])

    # first element in tuple is df of import risks where in this df indexes labels are ISOA2 of countries
    df = out[0]

    # second element is the monthly outflux
    monthly_outflux = out[1]
    daily_outflux = monthly_outflux/30.0

    # generate columns to be saved and compute seeds
    daily_seeds = [ round(IR*daily_outflux,5) for IR in list(df.values) ]
    iso_a3_targets = [ isoA2_to_isoA3[c] for c in list(df.index) ]

    # generate dict
    dict_df = {'iso_a3_target':iso_a3_targets,'daily_seeds':daily_seeds}
    df_to_save = pd.DataFrame.from_dict(dict_df)

    # save
    df_to_save.to_csv(PATH + f'{month}/daily_passengers_from_{SOURCE}_month_{month}.csv',index=False, sep = ' ')



# Main

if __name__ == "__main__":

    # Load parameters
    MONTHS_df = pd.read_csv('../data/stage3/months.csv', sep =',',index_col=False, names=['month'])
    MONTHS = list(MONTHS_df['month'])

    SOURCES_df = pd.read_csv('../data/stage3/seeds_sources.csv', sep = ',',index_col=False, names=['iso_a3'])
    SOURCES = list(SOURCES_df['iso_a3'])

    for month in MONTHS:
        generate_directory(f'../data/stage3_processed/seeds/{month}', month)

    # Dictionaries between countries ISOA* codes
    country_info = pd.read_csv("../data/stage3/countryInfo.csv")
    isoA2_to_isoA3, isoA3_to_isoA2 = {}, {}

    for index, row in country_info.iterrows():
        isoA2_to_isoA3[row['iso_alpha2']] = row['iso_alpha3']
        isoA3_to_isoA2[row['iso_alpha3']] = row['iso_alpha2']

    # Generate seeds
    for SOURCE in SOURCES:
        for month in MONTHS:
            generate_seeds_from_SOURCE_in_MONTH(SOURCE,month,isoA3_to_isoA2, isoA2_to_isoA3,\
                                                PATH= '../data/stage3_processed/seeds/')


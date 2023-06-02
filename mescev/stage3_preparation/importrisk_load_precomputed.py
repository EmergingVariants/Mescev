# -*- coding: utf-8 -*-
"""
Efficient computation of import probabilities given a flux network.
"""

import pandas as pd
from pathlib import Path


# local variables
d_data_ir = Path(__file__).absolute().parent.parent.parent / 'data/stage2/import_risks'


def round_date_2_month(date, fraction=2/3):
    lim = 30 * fraction
    d = pd.to_datetime(date)
    return f'{d.year}-{d.month + int(d.day > lim)}-01'


def get_file_names_precomputed(date, paras_risk=None,
                               top_traffic=0.99,
                               d_data=d_data_ir):
    '''
    INPUT:
        date
        paras_risk dict
            dictionary of optional parameters for ImportRiskAnalyzer
        d_data pathlib.Path object
            directory to where the precomputed data is located
    '''
    if paras_risk is None:
        paras_risk = {'exit_prob_weighted_by': 'geo_dist'}
    f_id = f'top_traffic_{top_traffic}__' + '__'.join([f'{k}_{v}' for k, v in paras_risk.items()])
    d = pd.to_datetime(round_date_2_month(date))
    dd = f'{d.year}-{d.month}'
    assert d_data.exists(), f'folder {d_data} does not exist'
    f_name_ir      = d_data / (f_id + f'__import_risk_{dd}.csv')
    f_name_outflux = d_data / (f_id + f'__outflux_{dd}.csv')
    return f_name_ir, f_name_outflux, d



def precomputed_import_risk_cntr_agg_all(date):
    '''
    load import risk that was computed previously via   "import_risk_all_cntr_aggregate"
    (takes between 10-20 minutes to compute 1 month)
    INPUT:
        date str
            e.g.: '2022-01-01' or '2022-01-22'
    OUTPUT:
        df_ir pd.DataFrame
            columns = source_iso2
            index   = target_iso2
        df_outflux pd.
    '''
    # standard parameters used in computation
    top_traffic = 0.99
    paras_risk = {'exit_prob_weighted_by': 'geo_dist'}
    f_name_ir, f_name_outflux, d = get_file_names_precomputed(date, paras_risk=paras_risk,
                                                          top_traffic=top_traffic,
                                                          d_data=d_data_ir)
    # load import risk
    f_name = f_name_ir
    assert f_name.exists(), f'{f_name} does not exists! Computation needed via "import_risk_all_cntr_aggregate"'
    df_ir = pd.read_csv(f_name)
    df_ir = df_ir.rename(columns={df_ir.columns[0]: (ci:= 'target_iso2')})\
                 .set_index(ci)
    # load outflux matrix
    f_name = f_name_outflux
    assert f_name.exists(), f'{f_name} does not exist! this should not happen since outflux is always computed with import risk'
    df_outflux = pd.read_csv(f_name)
    df_outflux = df_outflux.drop(columns=[df_outflux.columns[0]])\
                           .set_index('source_iso2')
    return df_ir, df_outflux


def precomputed_import_risk_cntr_agg(date, cntr):
    df_ir, df_outflux = precomputed_import_risk_cntr_agg_all(date)
    assert cntr in df_ir.columns, f'{cntr} is wrong iso2 code for country (or not in network)'
    return df_ir[cntr], df_outflux['outflux'][cntr]

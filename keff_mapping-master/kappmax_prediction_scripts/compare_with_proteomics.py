
# coding: utf-8

import json
import re
import pandas as pd
import numpy as np

from os.path import dirname, abspath
from six import string_types

import cobra

from kappmax_prediction_scripts.proteomics_helper_functions import \
    get_all_protein_localizations, get_metabolic_reactions_from_gene, \
    get_membrane_transport_genes, get_rna_modfication_genes, \
    get_sum_of_metabolic_fluxes_of_gene

import kappmax_prediction_scripts

scripts_dir = dirname(abspath(kappmax_prediction_scripts.__file__))
home_dir = dirname(scripts_dir)
resource_dir = home_dir + '/data/'
ijo = cobra.io.load_json_model('%s/iJO1366.json' % scripts_dir)

with open('%s/gene_name_to_bnum.json' % resource_dir, 'r') as f:
    bnum_to_gene = json.load(f)


def get_sim_raw_df(model, sim_filename):

    with open(sim_filename, 'r') as f:
        me_sim = json.load(f)

    series_dict = dict()
    model.solution = cobra.core.Solution(me_sim['biomass_dilution'],
                                         x_dict=me_sim, status='optimal')
    for key, value in model.get_translation_flux().items():
        gene = key.replace('translation_', '')
        series_dict[gene] = value

    sim_series = pd.Series(series_dict)

    return sim_series


def get_old_sim_raw_df(model, media, old_me_sims_filename):
    media = media.map_media_to_old_me_df[media]

    old_sim_df = pd.read_pickle(old_me_sims_filename)
    sim_df_filtered = old_sim_df[media]
    series_dict = {}
    for gene in sim_df_filtered.index:
        protein_id = 'protein_' + gene
        if protein_id not in model.metabolites:
            continue
        series_dict[gene] = sim_df_filtered[gene]

    old_sim_series = pd.Series(series_dict)
    return old_sim_series


def get_proteomics_data_raw_df(model, media, proteomics_data_path):
    # Load proteomics data and filter genes not modeled in ME-model
    data_df = pd.read_excel(proteomics_data_path)
    data_df = data_df.set_index('Gene').rename(index=bnum_to_gene)
    data_series = data_df[media].copy()
    cog_column = ['Annotated functional COG group (description)']

    data_dict = {}
    cog_dict = {}
    genes = [i.id.replace('RNA_', '') for i in
             model.metabolites.query(re.compile('RNA_b[0-9]'))]
    for gene in genes:
        protein_id = 'protein_' + gene
        # Some RNAs do not code proteins, skip these
        if protein_id not in model.metabolites:
            continue
        if gene in data_series.index:
            # Some genes have two entries. Take the average of these
            data_dict[gene] = data_series[gene].mean()

            # Handle genes with duplicate entries
            cog_value = data_df.loc[gene, cog_column]
            if isinstance(cog_value, pd.DataFrame):
                cog_value = cog_value.values[0]
            if isinstance(cog_value, pd.Series):
                cog_value = cog_value.values[0]
            if isinstance(cog_value, np.ndarray):
                cog_value = cog_value[0]
            if type(cog_value) == list:
                cog_value = cog_value[0]
            if not isinstance(cog_value, string_types) and \
                    type(cog_value) != float:
                raise UserWarning('Cog is bad', cog_value, type(cog_value))
            cog_dict[gene] = cog_value
    filtered_data_series = pd.Series(data_dict)
    cog_series = pd.Series(cog_dict)
    return filtered_data_series, cog_series


def transform_df_to_mass_or_mol_fraction(df, model, comparison_column_list,
                                         mass_fraction=True):
    for comparison in comparison_column_list:
        for gene in df.index:

            gene_obj = model.metabolites.get_by_id('protein_' + gene)
            mass = gene_obj.formula_weight / 1000. if mass_fraction else 1.

            df.loc[gene, comparison] *= mass

    df[comparison_column_list] = \
        df[comparison_column_list] / df[comparison_column_list].sum()


def filter_dataframe(df, filter_nonmetabolic=False, filter_lpp=False,
                     filter_membrane=False, filter_modifications=False,
                     filter_transport=False, low_metabolic_flux_filter=0):
    if filter_nonmetabolic:
        df = df.loc[df['Metabolic'] is True]

    if filter_lpp:
        df.drop('b1677', inplace=True)

    if filter_membrane:
        df = df.loc[df['Membrane_complex_associated'] is False]

    if filter_modifications:
        df = df.loc[df['Modification_gene'] is False]

    if filter_transport:
        df = df.loc[df['Transport_gene'] is False]

    if low_metabolic_flux_filter:
        df = df.loc[df['Metabolic flux'] > low_metabolic_flux_filter]

    return df


def add_dataframe_annotations(model, df):
    mod_genes = get_rna_modfication_genes(model)
    mem_genes = get_membrane_transport_genes(model)

    for gene in df.index:
        # Skip genes not in model
        if 'protein_' + str(gene) not in model.metabolites:
            continue

        df.loc[gene, 'Metabolic'] = True if gene in ijo.genes else False
        df.loc[gene, 'Modification_gene'] = True if gene in mod_genes else False
        df.loc[gene, 'Transport_gene'] = True if gene in mem_genes else False
        rxns = model.metabolites.get_by_id('protein_' + gene).metabolic_reactions
        df.loc[gene, 'Metabolic_reactions'] = str([i.stoichiometric_data.id for i in rxns])
        df.loc[gene, 'keffs'] = str([i.keff for i in rxns])
        df.loc[gene, 'Reaction_string'] = str([i.reaction for i in rxns])
        membrane_complex = False
        for protein in get_all_protein_localizations(model, gene):
            for r in protein.reactions:
                if 'formation_' not in r.id:
                    continue
                cplx = r.complex
                if cplx.compartment != 'c':
                    membrane_complex = True

        df.loc[gene, 'Membrane_complex_associated'] = membrane_complex
        df.loc[gene, 'Metabolic flux'] = \
            get_sum_of_metabolic_fluxes_of_gene(model, gene)


def return_dataframe(model, media, proteomics_data_path, simulation_path,
                     old_simulation_path=None):
    compare_dict = dict()
    data_df, cog_df = get_proteomics_data_raw_df(model, media,
                                                 proteomics_data_path)

    simulated_df = get_sim_raw_df(model, simulation_path)
    if old_simulation_path:
        old_simulated_df = get_old_sim_raw_df(model, media,
                                              old_simulation_path)

    compare_dict['Measured'] = data_df.to_dict()
    if old_simulation_path:
        compare_dict['Simulated_old_model'] = old_simulated_df.to_dict()

    compare_dict['Simulated'] = simulated_df.to_dict()
    compare_dict['COG'] = cog_df.to_dict()
    compare_df = pd.DataFrame(compare_dict)

    add_dataframe_annotations(model, compare_df)
    return compare_df


# coding: utf-8

# In[1]:

import cobrame
import pickle
import json
import numpy as np
import pandas as pd
from os.path import exists, dirname, abspath
from collections import defaultdict

import cobra

here = dirname(abspath(__file__))
parent = '/'.join(here.split('/')[:-1])
map_df = pd.read_csv('%s/data/david_m_id_to_me_id.csv' % parent,
                     index_col=0).fillna(0)
bigg_map_df = \
    pd.read_csv('%s/data/bigg_model_changes.csv' % parent)[['old_reaction',
                                                            'new_reaction']]


def data_frame_to_dict(df, m_to_me=True):
    out_dict = defaultdict(list)
    for i, v in df.iterrows():
        if m_to_me:
            out_dict[i].append(v['me_id'])
        else:
            out_dict[v['me_id']].append(i)
    return out_dict


david_id_to_me = data_frame_to_dict(map_df)
me_to_david_id = data_frame_to_dict(map_df, m_to_me=False)

bigg_id_to_me = bigg_map_df.set_index('new_reaction')

ijo = cobra.io.load_json_model('%s/iJO1366.json' % here)
iml = cobra.io.load_json_model('%s/iML1515.json' % here)
ijo_bigg = cobra.io.load_json_model('%s/data/iJO1366_bigg.json' % parent)


def append_me_keff_to_m_model_mapping(df, me_rxn_list, m_rxn):
    for r in me_rxn_list:
        df.loc[r, 'm_reaction'] = m_rxn


def handle_iron_sulfur_keffs(me, keff_series, out_df):

    def return_reactions(met_id, skip=list()):
        met = me.metabolites.get_by_id(met_id)
        rs = []
        for r in met.reactions:
            if met in r.reactants and r.id not in skip and not isinstance(r,
                                                                          cobrame.tRNAChargingReaction):
                rs.append(r.id)
        return rs

    def flatten_reactions(intermediate_reactions):
        return set(
            [item for sublist in intermediate_reactions for item in sublist])

    def add_keffs_to_carriers(carrier_complex_list, keff):
        for c in carrier_complex_list:
            me.process_data.get_by_id(c + '_carrier_activity').keff = keff

    # SCYSDS
    iml_rxn = 'SCYSDS'
    carrier_complexes = ['CPLX0-246_CPLX0-1342_mod_pydx5p']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # ICYSDS
    iml_rxn = 'ICYSDS'
    carrier_complexes = ['IscS_mod_2:pydx5p']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S2FE2SS (set carrier reaction keffs = S2FE2SS keff)
    # TODO consider maybe scaling this by number of reactions
    iml_rxn = 'S2FE2SS'
    carrier_complexes = ['CPLX0-1341', 'CPLX0-1341_mod_1:fe2',
                         'CPLX0-1341_mod_2:fe2', 'CPLX0-1341_mod_1:2fe1s',
                         'CPLX0-246_CPLX0-1342_mod_pydx5p_mod_1:SH']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])

    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S2FE2SR default to this keff over the above
    iml_rxn = 'S2FE2SR'
    carrier_complexes = ['CPLX0-1341_mod_1:2fe1s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I2FE2SS
    iml_rxn = 'I2FE2SS'
    # TODO EG11653-MONOMER is involved in Fe loading but not in iML1515.
    # Should this stay or be set to 65?
    carrier_complexes = ['IscU', 'EG11653-MONOMER_mod_1:fe2',
                         'IscU_mod_1:fe2', 'IscU_mod_2:fe2',
                         'IscS_mod_2:pydx5p_mod_1:SH', 'IscU_mod_1:2fe1s',
                         'IscU_mod_1:2fe1s', 'EG11653-MONOMER']
    # Reactions in skip use IscS_mod_2:pydx5p_mod_1:SH
    skip = ['BTS6_FWD_BIOTIN-SYN-CPLX_mod_4fe4s_mod_2fe2s',
            'THZPSN31_FWD_THIH-MONOMER_THIF-MONOMER_THII-MONOMER_THIS-MONOMER',
            'MOADSUx1_FWD_CPLX_dummy']
    flattened_reactions = \
        flatten_reactions([return_reactions(c, skip) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I2FE2SR default to this keff over the above
    iml_rxn = 'I2FE2SR'
    carrier_complexes = ['IscU_mod_1:2fe1s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c, skip) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S2FE2ST or I2FE2ST (keffs = average)
    set_keff = (keff_series['S2FE2ST'] + keff_series['I2FE2ST']) / 2.
    for mod in ['mod_2fe2s_c_G6712-MONOMER', 'mod_2fe2s_c']:
        me.process_data.get_by_id(mod).keff = set_keff
        append_me_keff_to_m_model_mapping(out_df, [mod],
                                          'Average of S2FE2ST and I2FE2ST')

    # I2FE2ST
    iml_rxn = 'I2FE2ST'
    carrier_complexes = ['IscA_tetra_mod_1:2fe2s',
                         'IscA_tetra',
                         'CPLX0-7824',
                         'CPLX0-7824_mod_1:2fe2s',
                         'IscU_mod_1:2fe2s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S2FE2ST
    iml_rxn = 'S2FE2ST'
    carrier_complexes = ['CPLX0-1341_mod_1:2fe2s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I2FE2ST or I4FE4ST (keffs = average).
    # This complex catalyzes the transfer of IscU iron sulfur clusters
    # to chaperones
    set_keff = (keff_series['I2FE2ST'] + keff_series['I4FE4ST']) / 2.
    carrier_complexes = ['EG12130-MONOMER_EG12131-MONOMER']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    for r in flattened_reactions:
        me.reactions.get_by_id(r).keff = set_keff

    # S2FE2SS2
    iml_rxn = 'S2FE2SS2'
    carrier_complexes = ['CPLX0-1341_mod_1:2fe2s',
                         'CPLX0-1341_mod_1:2fe2s_mod_1:fe2',
                         'CPLX0-1341_mod_1:2fe2s_mod_2:fe2',
                         'CPLX0-1341_mod_1:2fe2s_mod_1:2fe1s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S4FE4SR
    iml_rxn = 'S4FE4SR'
    carrier_complexes = ['CPLX0-1341_mod_2:2fe2s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I2FE2SS2
    iml_rxn = 'I2FE2SS2'
    carrier_complexes = ['IscU_mod_1:2fe2s', 'IscU_mod_1:2fe2s_mod_1:fe2',
                         'IscU_mod_1:2fe2s_mod_2:fe2',
                         'IscU_mod_1:2fe2s_mod_1:2fe1s']
    skip = ['2Fe2S_to_SufA_by_IscU_FWD_EG12130-MONOMER_EG12131-MONOMER',
            '2Fe2S_to_ErpA_by_IscU_FWD_EG12130-MONOMER_EG12131-MONOMER',
            '2Fe2S_to_IscA_by_IscU_FWD_EG12130-MONOMER_EG12131-MONOMER']
    flattened_reactions = \
        flatten_reactions([return_reactions(c, skip) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I4FE4SR
    iml_rxn = 'I4FE4SR'
    carrier_complexes = ['IscU_mod_2:2fe2s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # S4FE4ST or I4FE4ST (keffs = average)
    # generic transfer complexes are ['CPLX0-7617_mod_1:4fe4s',
    # 'CPLX0-7824_mod_1:4fe4s', 'IscA_tetra_mod_1:4fe4s']
    set_keff = (keff_series['S4FE4ST'] + keff_series['I4FE4ST']) / 2.
    for mod in ['mod_4fe4s_c', 'mod_3fe4s_c']:
        me.process_data.get_by_id(mod).keff = set_keff
        append_me_keff_to_m_model_mapping(out_df, [mod],
                                          'Average of S4FE4ST and I4FE4ST')

    # S4FE4ST
    iml_rxn = 'S4FE4ST'
    carrier_complexes = ['CPLX0-1341_mod_1:4fe4s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])

    # I4FE4ST
    iml_rxn = 'I4FE4ST'
    carrier_complexes = ['IscA_tetra_mod_1:4fe4s',
                         'CPLX0-7824_mod_1:4fe4s',
                         'CPLX0-7617',
                         'IscU_mod_1:4fe4s']
    flattened_reactions = \
        flatten_reactions([return_reactions(c) for c in carrier_complexes])
    append_me_keff_to_m_model_mapping(out_df, flattened_reactions, iml_rxn)
    add_keffs_to_carriers(carrier_complexes, keff_series[iml_rxn])


def map_m_id_to_me(me, orig_r, keff):
    rs = david_id_to_me.get(orig_r, [orig_r])
    for r in rs:
        if not r:
            return

        elif r == 'LIPOS':
            me.process_data.get_by_id(
                'CPLX0-782_mod_2:4fe4s_carrier_activity').keff = keff
            me.process_data.get_by_id(
                'EG50003-MONOMER_mod_pan4p_mod_oc_carrier_activity').keff = keff

        elif r == 'LIPAMPL':
            me.process_data.get_by_id('mod_lipo_c').keff = keff

        elif r == 'LIPOCT':
            me.process_data.get_by_id('mod_lipo_c_alt').keff = keff

        else:

            r_trunc = r.replace('_b', '').replace('_f', '')
            try:
                data = me.process_data.get_by_id(r_trunc)
            except:
                print(me.reactions.query(r))
                print(r, 'not in ijo')
                return

            for rxn in data.parent_reactions:
                if r.endswith('_b') and '_REV_' in rxn.id:
                    rxn.keff = keff
                    rxn.update()

                elif r.endswith('_f') and '_FWD_' in rxn.id:
                    rxn.keff = keff
                    rxn.update()

                elif '_FWD_' in rxn.id and not r.endswith('_b') \
                        and not r.endswith('_f'):
                    rxn.keff = keff
                    rxn.update()

                elif '_FWD_' in rxn.id and r.endswith('_b'):
                    pass

                elif '_REV_' in rxn.id and r.endswith('_f'):
                    pass

                # These reactions changed reversibility between iML and iJO
                elif r in ['5DGLCNR', 'ACOAD1f', 'LCARR', 'MOX', 'QUINDH',
                           'PPKr', 'PPK2r', 'ILEt2rpp', 'INDOLEt2rpp',
                           'LCTStpp', 'LEUt2rpp', 'PIt2rpp', 'SERt2rpp',
                           'THRt2rpp', 'VALt2rpp']:
                    rxn.keff = keff
                    rxn.update()
                elif '_REV_' in rxn.id and not r.endswith('_b'):
                    pass

                else:
                    print(rxn, r, rxn.reaction)
                    raise UserWarning('Error')


def get_keffs_from_m_reactions(enzyme, reactions, keff_series, m_keffs,
                               m_rxns):
    """
    For subreactions (carriers) that participate in many different reactions,
    use the average keff of all the m-model reactions that they are involved
    in.

    """
    import sympy
    for r in reactions:
        stoich = r._metabolites[enzyme]
        if isinstance(stoich, sympy.Basic) and stoich.subs(cobrame.mu, 0) > 0:
            continue
        elif stoich > 0:
            continue
        if not isinstance(r, cobrame.MetabolicReaction):

            continue
        if 'FWD' in r.id:
            stoich_id = r.stoichiometric_data.id
            m_rxn = [stoich_id] if stoich_id in ijo_bigg.reactions else None
            if not m_rxn:
                m_rxn = me_to_david_id.get(stoich_id, None)
            if not m_rxn:
                m_rxn = me_to_david_id.get(stoich_id + '_f', None)

            m_keffs.append(keff_series[m_rxn[0]])
            m_rxns.append(m_rxn[0])

        elif 'REV' in r.id:
            m_rxn = me_to_david_id.get(r.stoichiometric_data.id + '_b', None)
            m_keffs.append(keff_series[m_rxn[0]])
            m_rxns.append(m_rxn[0])

    return m_keffs, m_rxns


# Skip this. It's only involved in a metabolic process not included in the
# M-model
skip_process_data = ['CPLX0-782_mod_1:2fe2s_mod_1:4fe4s_carrier_activity']


def update_all_keffs(me, keff_series, objective_rxn='ATPM',
                     transporters=list()):
    out_df = pd.DataFrame()

    me.objective = objective_rxn

    for r in me.reactions:
        if hasattr(r, 'keff'):
            r.keff = 65.
    for r in me.process_data:
        if hasattr(r, 'keff'):
            r.keff = 65.

    # If has carrier considered then set keff of dummy for rxn very high
    for r in me.metabolites.CPLX_dummy.reactions:
        if isinstance(r, cobrame.ComplexFormation):
            continue
        if len(r.stoichiometric_data.subreactions) > 0:
            r.keff = 6000000.
            r.update()

    handle_iron_sulfur_keffs(me, keff_series, out_df)

    for r, keff in keff_series.items():
        # Ignore all membrane proteins in model
        if r.startswith('DM_') or 'BIOMASS' in r:
            continue

        # Keff series is for iML reactions. The ME models is only for iJO
        # reactions. Skip reactions only in iML
        if r not in ijo_bigg.reactions \
                and r.replace('_f', '') not in ijo_bigg.reactions \
                and r.replace('_b', '') not in ijo_bigg.reactions:
            continue

        if r in transporters:
            print(r)
            continue

        # Iron sulfur formation uses carriers as catalysts.
        #  These are highly interchangable. Used average keff
        if r.startswith('I2FE2') or r.startswith('I4FE4') or 'SCYSDS' in r:
            continue

        elif r.startswith('S2FE2') or r.startswith('S4FE4') or 'ICYSDS' in r:
            continue

        map_m_id_to_me(me, r, keff)

    for d in me.process_data:
        m_keffs = []
        m_rxns = []
        if isinstance(d, cobrame.TranslocationData) or d.id in \
                skip_process_data:
            continue
        if hasattr(d, 'keff') and d.enzyme and d.keff == 65. and '_at_' not in d.id:
            enzymes = [d.enzyme] if type(d.enzyme) == str else d.enzyme
            for enzyme_id in enzymes:
                enzyme = me.metabolites.get_by_id(enzyme_id)
                if isinstance(enzyme, cobrame.GenericComponent):
                    for e_id in me.process_data.get_by_id(enzyme.id).component_list:
                        e = me.metabolites.get_by_id(e_id)
                        if e.compartment.lower() != 'c':
                            continue
                        reactions = e.reactions

                        m_keffs, m_rxns = \
                            get_keffs_from_m_reactions(e, reactions,
                                                       keff_series, m_keffs,
                                                       m_rxns)
                else:
                    if enzyme.compartment.lower() != 'c':
                        continue
                    reactions = enzyme.reactions

                    m_keffs, m_rxns = \
                        get_keffs_from_m_reactions(enzyme, reactions,
                                                   keff_series, m_keffs,
                                                   m_rxns)
            if len(m_keffs) > 0:
                d.keff = np.array(m_keffs).mean()
                append_me_keff_to_m_model_mapping(out_df, [d.id], 'Average of ' + str(m_rxns))

            else:
                print(d.id, ' no metabolic processes')

    # Set keff of ACP carrier subreactions = to the mean of all keffs its
    # involved in
    vals = []
    for d in me.process_data.query('EG50003-MONOMER'):
        if '_carrier_activity' in d.id:
            vals.append(d.keff)
    for d in me.process_data.query('EG50003-MONOMER'):
        if '_carrier_activity' in d.id:
            d.keff = np.array(vals).mean()

    # TODO evaluate how necessary the code below is
    #for r in me.process_data.ATPS4rpp.parent_reactions:
    #    r.keff = 232.
    #    r.update()
    #    print(r.id, r.keff)

    #me.update()
    #for r in me.reactions:
    #    if r.upper_bound == 0 and 'PPKr' not in r.id and 'FHL' not in r.id:
    #        print(r.id)
    #        r.upper_bound = 1000.

    #me.process_data.MOX.lower_bound = 0
    #for r in me.process_data.MOX.parent_reactions:
    #    r.update()

    #for r in me.process_data.PDH.parent_reactions:
    #    r.keff = 3000.
    #    r.update()
    #    print(r.id, r.keff)
    me.update()
    print('DONE')

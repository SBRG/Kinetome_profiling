LB_media = {
    "EX_ni2_e": -1000,
    "EX_dcyt_e": -1000,
    "EX_hg2_e": -1000,
    "EX_ins_e": -1000,
    "EX_cd2_e": -1000,
    "EX_so4_e": -1000,
    "EX_uri_e": -1000,
    "EX_tungs_e": -1000,
    "EX_glu__L_e": -1000,
    "EX_slnt_e": -1000,
    "EX_trp__L_e": -1000,
    "EX_dad__2_e": -1000,
    "EX_mobd_e": -1000,
    "EX_val__L_e": -1000,
    "EX_cobalt2_e": -1000,
    "EX_gln__L_e": -1000,
    "EX_co2_e": -1000,
    "EX_k_e": -1000,
    "EX_cu2_e": -1000,
    "EX_sel_e": -1000,
    "EX_na1_e": -1000,
    "EX_cl_e": -1000,
    "EX_fe3_e": -1000,
    "EX_arg__L_e": -1000,
    "EX_pnto__R_e": -1000,
    "EX_lys__L_e": -1000,
    "EX_ala__L_e": -1000,
    "EX_gal_e": -1000,
    "EX_cbl1_e": -1000,
    "EX_ser__L_e": -1000,
    "EX_adn_e": -1000,
    "EX_thr__L_e": -1000,
    "EX_pi_e": -1000,
    "EX_thymd_e": -1000,
    "EX_mn2_e": -1000,
    "EX_phe__L_e": -1000,
    "EX_leu__L_e": -1000,
    "EX_ura_e": -1000,
    "EX_h_e": -100,
    "EX_h2o_e": -100,
    "EX_aso3_e": -1000,
    "EX_hxan_e": -1000,
    "EX_glc__D_e": -1000,
    "EX_nac_e": -1000,
    "EX_his__L_e": -1000,
    "EX_o2_e": -1000,
    "EX_pro__L_e": -1000,
    "EX_mg2_e": -1000,
    "EX_asp__L_e": -1000,
    "EX_gly_e": -1000,
    "EX_cys__L_e": -1000,
    "EX_fe2_e": -1000,
    "EX_ca2_e": -1000,
    "EX_tyr__L_e": -1000,
    "EX_zn2_e": -1000,
    "EX_fru_e": -1000,
    "EX_met__L_e": -1000,
    "EX_ile__L_e": -1000
}

aas = {"EX_glyc_e": -1000,
       "EX_asp__L_e": -1000,
       "EX_gly_e": -1000,
       "EX_cys__L_e": -1000,
       "EX_met__L_e": -1000,
       "EX_ile__L_e": -1000,
       "EX_tyr__L_e": -1000,
       "EX_pro__L_e": -1000,
       "EX_his__L_e": -1000,
       "EX_phe__L_e": -1000,
       "EX_leu__L_e": -1000,
       "EX_ser__L_e": -1000,
       "EX_arg__L_e": -1000,
       "EX_lys__L_e": -1000,
       "EX_ala__L_e": -1000,
       "EX_gln__L_e": -1000,
       "EX_glu__L_e": -1000,
       "EX_trp__L_e": -1000,
       "EX_val__L_e": -1000,
       "EX_thr__L_e": -1000,
       "EX_asn__L_e": -1000
       }

# Mapping of Aerbersold media conditions to exchange reaction
media_dict = {'Glucose': 'EX_glc__D_e', 'Acetate': 'EX_ac_e',
              'Pyruvate': 'EX_pyr_e', 'Glycerol': 'EX_glyc_e',
              'Fumarate': 'EX_fum_e', 'Succinate': 'EX_succ_e',
              'LB': '', 'Glucosamine': 'EX_gam_e',
              'Mannose': 'EX_man_e', 'Xylose': 'EX_xyl__D_e',
              'Fructose': 'EX_fru_e', 'Glycerol + AA': '',
              'Galactose': 'EX_gal_e', 'Gluconate': 'EX_glcn_e'}

map_media_to_old_me_df = {
    'Glucose': 'base', 'Acetate': 'Acetate', 'Fumarate': 'Fumarate',
    'Glycerol': 'Glycerol', 'Pyruvate': 'Pyruvate', 'Succinate': 'Succinate'
}


def set_media(model, name, value=-1000):
    model.reactions.EX_glc__D_e.lower_bound = 0

    reactions_changed = []
    if name in model.reactions:
        model.reactions.get_by_id(name).lower_bound = value
        reactions_changed.append(name)
    elif name == 'Glycerol + AA':
        for r, v in aas.items():
            model_rxn = model.reactions.get_by_id(r)
            if model_rxn.lower_bound == 0:
                model_rxn.lower_bound = v
                reactions_changed.append(r)
    elif name == 'LB':
        for r, v in LB_media.items():
            model_rxn = model.reactions.get_by_id(r)
            if model_rxn.lower_bound == 0:
                model_rxn.lower_bound = v
                reactions_changed.append(r)
    elif name in media_dict:
        model.reactions.get_by_id(media_dict[name]).lower_bound = value
        reactions_changed.append(media_dict[name])
    else:
        raise UserWarning('Media (s) not valid' % name)

    return reactions_changed

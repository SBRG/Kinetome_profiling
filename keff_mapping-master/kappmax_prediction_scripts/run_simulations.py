from __future__ import absolute_import

from os.path import exists
import pickle

import numpy as np
import json

from kappmax_prediction_scripts.media import set_media
from cobrame.solve.algorithms import binary_search, compile_expressions


def solve_at_sup(model, media, sweep_simulation_filename, sup=1000,
                 checkexist=True, solver='qminos', compiled_expressions=None):

    if checkexist and exists(sweep_simulation_filename):
        print('%s already exists, skipping' % sweep_simulation_filename)
        return

    set_media(model, media, value=-sup)

    sol = solve_model(model, solver=solver,
                      compiled_expressions=compiled_expressions)

    if not sol:
        print(sweep_simulation_filename, ' is infeasible')
        return

    with open(sweep_simulation_filename, 'w') as f:
        json.dump(sol.x_dict, f)


def solve_model(me, precision=1e-2, solver='qminos',
                compiled_expressions=None):

    if solver == 'qminos':
        from qminospy.me1 import ME_NLP1
        me_nlp = ME_NLP1(me, growth_key='mu')
        me_nlp.compiled_expressions = compiled_expressions
        x, status, hs = me_nlp.solvelp(.01)

        if status != 'optimal':
            print(status)
            return None
        print(me.solution.x_dict['biomass_dilution'])

        muopt, hs, xopt, cache = \
            me_nlp.bisectmu(precision=precision, basis=hs)
        me.solution.f = me.solution.x_dict['biomass_dilution']

    else:

        binary_search(me, min_mu=.01, max_mu=1.5,  mu_accuracy=precision,
                      solver=solver, compiled_expressions=compiled_expressions,
                      FeasibilityTol=1e-9)

    return me.solution


def maximize_growth_rate(model, media, simulation_filename,
                         checkexist=True, precision=1e-2, solver='qminos',
                         compiled_expressions=None):

    if checkexist and exists(simulation_filename):
        print('%s already exists, skipping' % simulation_filename)
        return

    # change media, function sets glucose uptake to 0
    reactions_changed = set_media(model, media)

    sol = solve_model(model, precision, solver,
                      compiled_expressions=compiled_expressions)

    if not sol:
        print(simulation_filename, ' is infeasible')
        return

    with open(simulation_filename, 'w') as f:
        json.dump(sol.x_dict, f)

    with open(simulation_filename.replace('.json', '_met_flux.json'), 'w') as f:
        json.dump(model.get_metabolic_flux(solution=sol), f)

    for r in reactions_changed:
        model.reactions.get_by_id(r).lower_bound = 0

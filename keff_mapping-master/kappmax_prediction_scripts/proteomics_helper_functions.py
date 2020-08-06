import cobrame

# possible protein compartments
compartments = ['', '_Inner_Membrane', '_Outer_Membrane', '_Periplasm',
                '_lipoprotein_Inner_Membrane',
                '_lipoprotein_Outer_Membrane']

def get_all_protein_localizations(model, gene):
    """Get all possible locations for a given gene in model without using
    query. Query is slow"""

    proteins = set()

    for compartment in compartments:
        protein_id = 'protein_' + gene + compartment
        if protein_id in model.metabolites:
            proteins.add(model.metabolites.get_by_id(protein_id))

    return proteins


def get_metabolic_reactions_from_gene(model, gene):
    def get_gene_obj(gene_id):
        genes = list(get_all_protein_localizations(model, gene_id))
        if len(genes) == 2:
            # second gene corresponds to compartment (inner or outer membrane)
            genes.remove(model.metabolites.get_by_id('protein_' + gene_id))
        elif len(genes) == 1:
            pass
        else:
            # this corresponds to lipoproteins
            #  assume lipoprotein is last in list
            return genes[-1]
        return genes[0]

    def get_complexes(gene_object):
        list_of_complexes = []
        for reaction in gene_object.reactions:
            if hasattr(reaction, 'complex'):
                list_of_complexes.append(reaction.complex)
        return list_of_complexes

    def get_metabolic_reactions(list_of_complexes):
        metabolic_reactions = []
        for complexes in list_of_complexes:
            metabolic_reactions.extend(complexes.metabolic_reactions)
        return metabolic_reactions

    gene_obj = get_gene_obj(gene)
    complex_list = get_complexes(gene_obj)
    return set(get_metabolic_reactions(complex_list))


def get_sum_of_metabolic_fluxes_of_gene(model, gene):
    met_rxns = get_metabolic_reactions_from_gene(model, gene)

    if not met_rxns:
        return

    total = 0
    for r in met_rxns:
        total += model.solution.x_dict[r.id]

    return total


def get_proteins_in_membrane_complexes(model):
    membrane_proteins = set()
    for gene in model.metabolites.query('protein_'):
        for r in gene.reactions:
            if not r.id.startswith('formation_'):
                continue
            if r.complex.compartment in ['p', 'e', 'im', 'om', 'mc', 'm']:
                membrane_proteins.add(gene.id.split('_')[1])
    return membrane_proteins


# Identify genes that should be skipped
def get_membrane_transport_genes(model):
    skip_genes = []
    for r in model.reactions:
        if isinstance(r, cobrame.MetabolicReaction):
            if r.subsystem == 'Transport, Inner Membrane':
                if r.complex_data:
                    skip_genes.extend([i.split('_')[1] for i in
                                       r.complex_data.stoichiometry.keys()])
            if r.subsystem == 'Inorganic Ion Transport and Metabolism' and \
                    'pp' in r.id:
                if r.complex_data:
                    skip_genes.extend([i.split('_')[1] for i in
                                       r.complex_data.stoichiometry.keys()])
            if r.subsystem == 'Transport, Outer Membrane':
                if r.complex_data:
                    skip_genes.extend([i.split('_')[1] for i in
                                       r.complex_data.stoichiometry.keys()])

    skip_genes = set(skip_genes)
    return skip_genes


def get_rna_modfication_genes(model):

    # Get process data corresponding to each rna modification
    rna_modification_data = set()
    for d in model.tRNA_data:
        for mod in d.subreactions:
            if 'at_' in mod:
                rna_modification_data.add(mod)

    for mod in model.process_data.ribosome.subreactions:
        if 'at_' in mod:
            rna_modification_data.add(mod)

    # Get genes associated with each modification
    mod_genes = set()
    for mod in rna_modification_data:
        data = model.process_data.get_by_id(mod)
        if not data.enzyme:
            continue
        enzymes = [data.enzyme] if type(data.enzyme) == str else data.enzyme

        for enzyme in enzymes:
            # include carrier activity in analysis
            if enzyme + '_carrier_activity' in model.process_data:
                continue

            # generic versions of some enzymes
            data = model.process_data.get_by_id(enzyme)
            if isinstance(data, cobrame.GenericData):
                for enz in data.component_list:
                    enzyme_data = model.process_data.get_by_id(enz)
                    mod_genes.update(set([i.replace('protein_', '') for i in
                                          enzyme_data.stoichiometry.keys()]))
            else:
                enzyme_data = data
                mod_genes.update(set([i.replace('protein_', '') for i in
                                      enzyme_data.stoichiometry.keys()]))

    return mod_genes

"""
Created on Mon Dec  4 11:52:50 2023

Functions to analyse sequence properties of different AMPs and correlate with phenotype at single-cell level

@author: fragasso
"""
import numpy as np
import pandas as pd


aa_weights = {     # amino acid molecular weighths in [g/mol]
    "A": 89.09,    # Alanine
    "R": 174.20,   # Arginine
    "N": 132.12,   # Asparagine
    "D": 133.10,   # Aspartic Acid
    "C": 121.16,   # Cysteine
    "Q": 146.15,   # Glutamine
    "E": 147.13,   # Glutamic Acid
    "G": 75.07,    # Glycine
    "H": 155.16,   # Histidine
    "I": 131.17,   # Isoleucine
    "L": 131.17,   # Leucine
    "K": 146.19,   # Lysine
    "F": 165.19,   # Phenylalanine
    "M": 149.21,   # Methionine
    "P": 115.13,   # Proline
    "S": 105.09,   # Serine
    "T": 119.12,   # Threonine
    "W": 204.23,   # Tryptophan
    "Y": 181.19,   # Tyrosine
    "V": 117.15    # Valine
}

aa_list = list(aa_weights.keys())

aa_charges = {
    "A": 0,  # Alanine
    "R": +1, # Arginine (basic)
    "N": 0,  # Asparagine
    "D": -1, # Aspartic Acid (acidic)
    "C": 0,  # Cysteine
    "Q": 0,  # Glutamine
    "E": -1, # Glutamic Acid (acidic)
    "G": 0,  # Glycine
    "H": 0, # Histidine (basic, but can be neutral depending on the environment)
    "I": 0,  # Isoleucine
    "L": 0,  # Leucine
    "K": +1, # Lysine (basic)
    "M": 0,  # Methionine
    "F": 0,  # Phenylalanine
    "P": 0,  # Proline
    "S": 0,  # Serine
    "T": 0,  # Threonine
    "W": 0,  # Tryptophan
    "Y": 0,  # Tyrosine
    "V": 0   # Valine
}

# Hydropathy values based on the Kyte-Doolittle scale
hydropathy_values = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

# classification based on https://www.pnas.org/doi/10.1073/pnas.2019053118#sec-3
# 'Learning the molecular grammar of protein condensates from sequence determinants and embeddings', by Saar et al., 2021
polar_amino_acids = {'S', 'T', 'N', 'Q', 'D', 'E', 'K', 'R', 'H', 'P', 'G', 'Y'}
hydrophobic_amino_acids = {'A','I', 'L', 'M', 'F', 'V', 'C','W'}
aromatic_amino_acids = {'W', 'Y', 'F'}
cationic_residues = {'K', 'R', 'H'}
anionic_residues = {'D', 'E'}


Mw_H2O = 18.015 # [g/mol] molecular weight of water . Each peptide bond formation causes release of a water molecule.
amidation_weight = 0.98 #[g/mol]  weight loss when amidation of the carboxylic C-terminal group takes place
disulfide_weight = 2.02 #[g/mol] weight loss due to loss of two hydrogen atom when disulfide bridge forms


def get_gravy(sequence):
    total_hydropathy = sum(hydropathy_values[aa] for aa in sequence if aa in hydropathy_values)
    return total_hydropathy / len(sequence)

def get_polar_fraction(sequence):
    polar_count = sum(1 for aa in sequence if aa in polar_amino_acids)
    return polar_count / len(sequence)

def get_hydrophobic_fraction(sequence):
    hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic_amino_acids)
    return hydrophobic_count / len(sequence)

def get_aromatic_fraction(sequence):
    aromatic_count = sum(1 for aa in sequence if aa in aromatic_amino_acids)
    return aromatic_count / len(sequence)

def get_cationic_fraction(sequence,amidation):
    cationic_count = sum(1 for aa in sequence if aa in cationic_residues) + amidation
    return cationic_count / len(sequence)

def get_anionic_fraction(sequence):
    anionic_count = sum(1 for aa in sequence if aa in anionic_residues)
    return anionic_count / len(sequence)

def get_proline_fraction(sequence):
    proline_count = sum(1 for aa in sequence if aa == 'P')
    return proline_count / len(sequence)

def get_valine_fraction(sequence):
    proline_count = sum(1 for aa in sequence if aa == 'V')
    return proline_count / len(sequence)


def get_arginine_fraction(sequence):
    arginine_count = sum(1 for aa in sequence if aa == 'R')
    return arginine_count / len(sequence)

def get_lysine_fraction(sequence):
    lysine_count = sum(1 for aa in sequence if aa == 'K')
    return lysine_count / len(sequence)

def amino_acid_composition(sequence):
    amino_acids = aa_list

    total_count = len(sequence)

    # Initialize composition dictionary with zeros
    composition = {aa: 0 for aa in amino_acids}

    # Count occurrences of each amino acid in the sequence
    for aa in sequence:
        if aa in composition:
            composition[aa] += 1

    for aa in composition:
       composition[aa] = (composition[aa] / total_count) * 100

    composition_array = np.array([composition[aa] for aa in amino_acids])

    return composition, composition_array


def get_protein_df(label,aa_seq,amidated,disulfide_bridges, source = []):
    df = pd.DataFrame()
    N = []
    weights =[]
    net_charge= []
    gravy = []
    polar = []
    hydrop = []
    aromatic = []
    cationic = []
    anionic = []
    proline = []
    valine = []
    arginine = []
    lysine = []
    composition = []
    comp_array = []
    for k in range(len(label)):
        seq = aa_seq[k]
        amidation = amidated[k]
        weight = sum(aa_weights[p] for p in seq) - Mw_H2O*(len(seq)-1)-amidation*amidation_weight-disulfide_bridges[k]*disulfide_weight
        charge = sum(aa_charges[p] for p in seq) + amidation
        N.append(len(seq))
        weights.append(weight)
        net_charge.append(charge)
        gravy.append(get_gravy(seq))
        polar.append(get_polar_fraction(seq))
        hydrop.append(get_hydrophobic_fraction(seq))
        aromatic.append(get_aromatic_fraction(seq))
        cationic.append(get_cationic_fraction(seq,amidation))
        anionic.append(get_anionic_fraction(seq))
        proline.append(get_proline_fraction(seq))
        valine.append(get_valine_fraction(seq))
        arginine.append(get_arginine_fraction(seq))
        lysine.append(get_lysine_fraction(seq))
        composition.append(amino_acid_composition(seq)[0])
        comp_array.append(amino_acid_composition(seq)[1])
    df['Name'] = label
    if len(source)>0:
        df['Source'] = source
    df['Sequence'] = aa_seq
    df['Length'] = N
    df['Molecular weight (Da)'] = weights
    df['Net charge'] = net_charge
    df['Gravy index']=gravy
    df['Polar fraction']=polar
    df['Hydrophobic fraction']=hydrop
    df['Aromatic fraction']=aromatic
    df['Cationic fraction']=cationic
    df['anionic_fraction']=anionic
    df['Proline fraction']=proline
    df['valine_fraction']=valine
    df['Arginine fraction']=arginine
    df['Lysine fraction']=lysine
    df['aa_composition']=composition
    df['composition_array']=comp_array


    df['NH2'] = amidated
    df['Cys-Cys (#)'] = disulfide_bridges

    return df

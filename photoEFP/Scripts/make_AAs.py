# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:13:19 2024

@author: jackl

This script creates .inp files for amino acids or molecules located in the EFP region.
Non-amino acid molecules are treated with no virtual bonds (i.e. no broken bonds).

Reads:
    structure file with only EFP-region atoms and residues (.g96)
    full structure file will all atoms and residues (.g96)
    User-created text file with QM atoms and QM-MM covalent bond definitions (.txt)
    Topology file from GROMACS molecular dynamics (.top or itp)
            #WARNING# If you have your topology separated into several .itp files instead of one master .top file,
            this script will have problems because atom numbers do not match between the structure and topology.
            For general use, this is not a problem. If you have special QM-MM covalent "bridge-bonds" within 
            an .itp file (with different atom numbers than the full structure file), then this script will fail.
    
"""

import numpy as np
import sys

# ---------------------------
# Input Files and Command Line Arguments
# ---------------------------
# Get input filenames from command line
# Sample execution: python make_AAs.py efp_pair53004.g96 confout_pair53004.g96 user_defined.txt topol.top

efp_g96 = sys.argv[1]     # EFP file ("efp_pair53004.g96")
full_g96 = sys.argv[2]    # Full configuration file ("confout_pair53004.g96")
settings_file = sys.argv[3]   # User settings file
topol_file = sys.argv[4]      # Topology file ("topol.top")

# Read the EFP file lines
with open(efp_g96, 'r') as f:
    efp_lines = f.readlines()

# Read the full configuration file lines
with open(full_g96, 'r') as f:
    full_lines = f.readlines()

# ---------------------------
# Global Data and Parameters
# ---------------------------
# Dictionary of amino acid charges. Assumes any N-terminal is +1 and any C-terminal is -1.
AA_charge = {
    'ASP': '-1', 'GLU': '-1', 'HIP': '1', 
    'LYS': '1', 'ARG': '1', 'HISH': '1'
}
# Amino acids with non-zero charge
spec_AAs = ['ASP', 'GLU', 'HIP', 'LYS', 'ARG', 'HISH']

# Map full amino acid names to one-letter codes for filenames
amino_acid_dict = {
    'ALA': 'a', 'ARG': 'r', 'ASN': 'n', 'ASP': 'd', 'CYS': 'c', 
    'GLN': 'q', 'GLU': 'e', 'GLY': 'g', 'HIS': 'h', 'ILE': 'i', 
    'LEU': 'l', 'LYS': 'k', 'MET': 'm', 'PHE': 'f', 'PRO': 'p', 
    'SER': 's', 'THR': 't', 'TRP': 'w', 'TYR': 'y', 'VAL': 'v', 
    'HIP': 'hp', 'HID': 'hd', 'HIE': 'he', 'NVAL': 'v', 'CGLN': 'q'
}

# List of known amino acids
known_amino_acids = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
    'TYR', 'VAL', 'HIP', 'HID', 'HIE', 'HISE', 'HISD', 'HISH'
]

# Atomic symbols and associated (string) numbers
at_sym = {
    'H': '1.0', 'C': '6.0', 'N': '7.0', 'O': '8.0', 'MG': '12.0', 
    'P': '15.0', 'S': '16.0', 'FE': '26.0', 'NA': '11.0', 'CL': '17.0'
}

# ---------------------------
# Function Definitions
# ---------------------------
def qm_atoms(set_file):
    """
    Grab the QM_atoms from the user-defined settings file.
    
    Returns:
        QMs: List of QM atom index numbers
        pol_rem: List of atom indices for polarization removal
    """
    qm = 0
    QMs = []
    pol_rem = []
    border = 0
    i = 3
    with open(set_file, 'r') as input_file:
        inp_lines = input_file.readlines()
    for line in inp_lines:
        # After the border, only one line is expected per atom.
        if border == 1:
            if 'boundary' in line:
                i = 0
            elif i < 1:
                if int(line.split()[3]) in QMs:
                    i += 1
                else:
                    print('Atom number: ' + line.split()[3] + ' not found in QM atoms. Check bridge atom order.')
                    break
            elif i == 1:
                QMs.append(int(line.split()[3]))
                pol_rem.append(int(line.split()[3]))
        # Look for QM_atoms section.
        elif qm == 1:
            if 'QM-MM' in line:
                border = 1
                qm = 0
            if len(line.split()) > 3:
                QMs.append(int(line.split()[3]))
        elif 'QM_atoms' in line:
            qm = 1
    return QMs, pol_rem


def cut_frag(head, tail):
    """
    Calculate coordinates for a virutal hydrogen  using a desired bond distance.
    
    Parameters:
        head: Line containing the head atom information.
        tail: Line containing the tail atom information.
    
    Returns:
        A list with the virtual coordinate [x, y, z] for the H-atom.
    """
    desired_dist = 1.07886  # Desired bond distance
    # Multiply by 10 to scale the coordinates.
    xh, yh, zh = [float(head.split()[i]) * 10 for i in range(4, 7)]
    xt, yt, zt = [float(tail.split()[i]) * 10 for i in range(4, 7)]
    # Calculate the magnitude of the distance vector.
    dist_mag = np.sqrt((xh - xt)**2 + (yh - yt)**2 + (zh - zt)**2)
    # Calculate the new coordinate at the desired distance.
    h_t = [
        ((xt - xh) * desired_dist / dist_mag) + xh,
        ((yt - yh) * desired_dist / dist_mag) + yh,
        ((zt - zh) * desired_dist / dist_mag) + zh
    ]
    return h_t


def make_inp(fragment, QMs, POLs):
    """
    Generate an input (.inp) file for a fragment.
    
    Parameters:
        fragment: List of lines containing atoms in the fragment.
        QMs: List of QM atom names to be marked for removal.
        POLs: List of atom names to have polarization points removed.
    """
    lines_out = []
    num_virtuals = 0
    charge = 0

    # Count virtual atoms (labeled as 'H000')
    for atom in fragment:
        if 'H000' in atom:
            num_virtuals += 1

    # If all atoms are either QM or virtual atoms, do not write an input.
    if len(QMs) + num_virtuals == len(fragment):
        return

    # Determine the charge based on the residue type.
    res_info = fragment[4].split()[1]
    if res_info in spec_AAs:
        charge = AA_charge[res_info]
    elif len(res_info) == 4:
        # Assume N-terminal or C-terminal if the residue code starts with N or C.
        if res_info[0] == 'N':
            charge = 1
        elif res_info[0] == 'C':
            charge = -1

    # Construct the filename based on residue information. 
    #Use 1 letter identifier for Amino acids, otherwise filename will match residue name
    if res_info in known_amino_acids:
        filename = amino_acid_dict[res_info] + '_' + fragment[4].split()[0] + '_' + fragment[0].split()[3] + '.inp'
    else:
        filename = res_info.lower() + '_' + fragment[4].split()[0] + '_' + fragment[0].split()[3] + '.inp'

    # Build header sections for the input file.
    # You may want to adjust these parameters.
    lines_out.append(
        f" $contrl units=angs local=boys runtyp=makefp \n"
        f"       mult=1 icharg={charge} coord=cart icut=11 $end\n"
        f" $system timlim=99999 mwords=200 $end\n"
        f" $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06 $end\n"
        f" $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n"
        f" $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n"
        f" $MAKEFP POL=.t. DISP=.f. CHTR=.f. EXREP=.f. $end\n"
        f" $data\n {filename.split('.')[0]}\n C1\n"
    )

    # Process each atom line.
    for atom in fragment:
        if 'H000' in atom:
            # Virtual atoms are written as is.
            lines_out.append(atom)
        else:
            # Format the atom line.
            # First column: atom identifier from the fragment.
            col1 = f" {atom.split()[2]}".ljust(6)
            # Determine atomic number based on element symbol.
            # Two-letter atomnames will need to be added here; this is a lazy solution.
            if atom.split()[2][0] == 'M':
                col2 = at_sym['MG']
            else:
                col2 = at_sym[atom.split()[2][0]]
            # Get coordinates, scaling by 10. Then format
            x, y, z = [float(atom.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(17)
            col4 = f"{y:.8f}".rjust(18)
            col5 = f"{z:.8f}".rjust(18)
            lines_out.append(f"{col1}{col2}{col3}{col4}{col5}\n")

    # Append comments regarding atoms to be removed later.
    lines_out.append(" $end \n!comment atoms to be erased:")
    for atomname in QMs:
        lines_out.append(" " + str(atomname))
    lines_out.append("\n")
    # Append comments regarding atoms that will have polarization removed.
    lines_out.append("!polarization points to remove:")
    for atomname in POLs:
        lines_out.append(" " + str(atomname))

    # Write the constructed lines to the file.
    with open(filename, 'w') as outfile:
        outfile.writelines(lines_out)


def QM_MM_covalent(MM, QMs, topol_file):
    """
    For a given MM atom that is covalently bound to a QM atom, find other MM atoms that are covalently bonded
    to the given atom but are not in the QM region.
    
    This function parses the topology file for bonds and returns a list of atom IDs.
        WARNING this script assumes matching atom IDs between the full structure .g96 and the topology .top files;
        -the script cannot see how your topology is defined.
    
    Parameters:
        MM: The MM atom index.
        QMs: List of QM atom indices.
        topol_file: Path to the topology file.
    
    Returns:
        bond_IDs: List of MM atom indices covalently bonded to the given MM atom.
    """
    bond_IDs = []
    bonds_section = False
    with open(topol_file, 'r') as f:
        for line in f:
            # Identify the bonds section.
            if '[ bonds ]' in line:
                bonds_section = True
                continue
            elif bonds_section and len(line.split()) < 2:
                # End of bonds section.
                return bond_IDs
            elif bonds_section and line.strip().startswith(';'):
                # Skip comment lines.
                continue
            elif bonds_section:
                atom1 = int(line.split()[0])
                atom2 = int(line.split()[1])
                # If one atom is given MM and the other is not in the QM region, add the non‐QM atom.
                if atom1 == MM and atom2 not in QMs:
                    bond_IDs.append(atom2)
                elif atom2 == MM and atom1 not in QMs:
                    bond_IDs.append(atom1)
    return bond_IDs

# ---------------------------
# Main Script Processing
# ---------------------------

# Determine QM region atoms and polarization removal atoms from the settings file.
qm_IDs, MM_remove = qm_atoms(settings_file)

# For each MM atom covalently bound to the QM region, find the other MM atoms bound
# to it. These atoms will have polarization points removed.
pol_remove = []
for atom in MM_remove:
    found_pol_remove = QM_MM_covalent(atom, qm_IDs, topol_file)
    for bound_atom in found_pol_remove:
        pol_remove.append(bound_atom)

# ---------------------------
# Process EFP Region Residues
# ---------------------------
# Get a list of unique residue IDs from the EFP file.
#   Note that atom IDs do not match between EFP and full structures. However, residue IDs do match
efp_resis = []
prevres = '0'
start = False
for line in efp_lines:
    if start:
        if len(line.split()) < 3:
            break
        #only count residues once
        if line.split()[0] != prevres:
            efp_resis.append(line.split()[0])
            prevres = line.split()[0]
    if 'POSITION' in line:
        start = True

# ---------------------------
# Process Fragments from Full Configuration File
# ---------------------------
i = 0              # Index for current EFP residue
prev_co = []       # List to store previous "C" or "O" atom lines; these will be 1 residue number "off" from .g96 convention
frag = []          # List to append all fragment lines
CAs = []           # List to accumulate "CA" (alpha carbon) lines
start = False

#"C", "O", and "CA" atoms are significant. We redefine residues using these atom names

for line in full_lines:
    if 'END' in line:
        if start:
            break
    #elif 'SOL' in line:
    elif 'SOL' in line or 'QSL' in line:        #Sol is standard water, QSL is my QM region water, this is not standardized
        continue
    elif start:
        # Collect CA atoms.
        if line.split()[2] == 'CA':
            CAs.append(line)
        # Collect "O" atoms.
        elif line.split()[2] == 'O':
            prev_co.append(line)
        # Process "C" atoms. This is the "end" of a residue as we would like to define it.
        elif line.split()[2] == 'C':
            prev_co.append(line)
            # If the current line belongs to the current EFP residue.
            if line.split()[0] == efp_resis[i]:
                # standard amino acid fragments (non-terminals, non-cofactors) need two virtual hydrogens.
                if line.split()[1] in known_amino_acids:
                    vH1 = cut_frag(frag[0], CAs[-2])      #Cut between current C and previous CA
                    vH2 = cut_frag(CAs[-1], line)         #Cut between current CA and next C
                    frag.append(f" H000 1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
                    frag.append(f" H000 1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
                # If the residue identifier (field 1) has length 4 (a non-standard amino acid),
                # determine if it is N-terminal or C-terminal.
                elif len(line.split()[1]) == 4:
                    if line.split()[1][1:4] in known_amino_acids:
                        if line.split()[1][0] == 'N':  # N-terminal residue
                            vH2 = cut_frag(CAs[-1], line)
                            frag.append(f" H000    1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
                        elif line.split()[1][0] == 'C':  # C-terminal residue
                            vH1 = cut_frag(frag[0], CAs[-2])
                            frag.append(f" H000    1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
                # If there is enough fragment data, determine which atoms to remove.
                if len(frag) > 1:
                    qm_names = []
                    pol_names = []
                    for atom in frag:
                        if 'H000' in atom:               #All virtual hydrogens will be removed, no need to specify in comment
                            continue
                        elif int(atom.split()[3]) in qm_IDs:
                            qm_names.append(atom.split()[2])
                        elif int(atom.split()[3]) in pol_remove:
                            pol_names.append(atom.split()[2])
                    make_inp(frag, qm_names, pol_names)
                    frag = []  # Reset fragment lines
                    i += 1  # Move to the next EFP residue.
        # Check if the fragment is complete (with more than one line)
        # Check that this frag is the current EFP residue (and end of that residue)
        # Check that this frag is not an amino acid (is a cofactor, lipid, etc) [redundant]
        #         No virtual hydrogens are added here
        if len(frag) > 1 and frag[-1].split()[0] == efp_resis[i]:
            if (frag[0].split()[0] != line.split()[0]) and frag[0].split()[1] not in known_amino_acids:
                if frag[0].split()[1][1:] not in known_amino_acids:
                    qm_names = []
                    pol_names = []
                    for atom in frag:
                        if int(atom.split()[3]) in qm_IDs:
                            qm_names.append(atom.split()[2])
                        elif int(atom.split()[3]) in pol_remove:
                            pol_names.append(atom.split()[2])
                    make_inp(frag, qm_names, pol_names)
                    frag = []
                    i += 1
        # If the current line belongs to the current EFP residue, add it to the fragment.
        if line.split()[0] == efp_resis[i]:
            if len(prev_co) > 1 and line.split()[1] in known_amino_acids:
                # Add the previous two "C" or "O" atoms to the fragment if the current fragment is an amino acid
                frag.append(prev_co[-2])
                frag.append(prev_co[-1])
                prev_co = []
            frag.append(line)
    elif 'POSITION' in line:
        start = True
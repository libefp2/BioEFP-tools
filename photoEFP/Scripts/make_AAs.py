# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 15:13:19 2024

@author: jackl & lyuda

This script takes as input full system g96 file and efp-shell g96 file, as well as 
the file with g96-like lines defining the QM region and QM-EFP boundaries, and system protein topology file.
Additionally, if some of the cofactors are not included in the protein topology file, 
the path to directory with their itp files needs to be adjusted (PATH_TO_AMBER parameter).

Sample execution: 
    python make_AAs.py shell_bchl361-79002.g96 bchl361-79002.g96 user_defined.txt topol.top 50000 wt

- The script prepares makefp input files for amino acids located in the EFP region, taking care of 
proper capping with hydrogens.
- Non-standard ligands are assumed to be a single fragment; however, fragmentation can be introduced 
if needed by adding information into 'define_ligand_cut' fucntion. 
- The script also creates MM region fragment (prot.efp)
- The script prepares text files with sections for Q-Chem QM/EFP calculations - 
QM region section, fragment section (including water molecules and prot.efp fragment)

Check carefully and adjust Global Data and Parameters section in case you have any non-standard 
amino acids, co-factors, etc. 

"""

import numpy as np
import sys
import os
import argparse

# ---------------------------
# Input Files and Command Line Arguments
# ---------------------------
# Get input filenames from command line

parser = argparse.ArgumentParser(
    description='Prepare GAMESS input files for QM/EFP calculations.'
)
parser.add_argument('efp_g96',      help='EFP shell file in g96 format')
parser.add_argument('full_g96',     help='Full system structure file in g96 format')
parser.add_argument('settings_file',help='User made file defining QM region and QM-EFP boundaries')
parser.add_argument('topol_file',   help='Topology file (topol.top)')
parser.add_argument('timestamp',    help='Timestamp to differentiate output files (50000)')
parser.add_argument('mutant',       help='Mutant label (wt)')
args = parser.parse_args()

efp_g96       = args.efp_g96
full_g96      = args.full_g96
settings_file = args.settings_file
topol_file    = args.topol_file
timestamp     = args.timestamp
mutant        = args.mutant

# ---------------------------
# Global Data and Parameters
# ---------------------------

# directory with itp topolgies of non-standard cofactors
PATH_TO_AMBER = '/depot/lslipche/data/PS1/WT/amber03.ff/'

# g96 atom names with more than 1 character that need to be renamed for Q-Chem calculations
# onl relevant for atoms in QM region
ATOM_EXCEPTIONS = {'mg1' : 'Mg', 'MG' : 'Mg'}

# charge and multiplicity of QM region in QM/EFP calculation
qm_charge = 0
qm_multiplicity = 1

nm2Bohr = 18.897161646321  # nm -> Bohr

# RESNAMEs in g96 that should be completely ignored for either efp, mm or qm regions
 # (XXX is a link atom, generally)
skip_resnames=['XXX']

#All non-amino acid residue names with separate topologies should be in this list
special_cofactors = ['ECH', '45D', 'EQ3', 'C7Z', 'CLA', 'PQN', 'BCR', 'QLA', 'LHG', 'LMG', 'SQD', 'LMT']

# 'SOL' returns the charge of oxygen in the water model (TIP3P). H charge is taken to be -1/2 * oxygen charge.
# Change these for different water models.
WATER_AND_IONS = {
    'SOL': -0.834,
    'CL': -1.0,
    'NA': 1.0
}

ions = ['NA', 'CL']

# used water names
waters = ['SOL', 'QSL']

# List of known amino acids
# these aminoacids follow standard amino acid fragmentation and 
# found in main topogly file
known_amino_acids = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 
    'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 
    'TYR', 'VAL', 'HIP', 'HID', 'HIE', 'HISE', 'HISD', 'HISH', 'ACYS'
    ]

# Atomic symbols and associated nuclear charges in string format
at_sym = {
    'H': '1.0', 'C': '6.0', 'N': '7.0', 'O': '8.0', 'MG': '12.0', 
    'P': '15.0', 'S': '16.0', 'FE': '26.0', 'NA': '11.0', 'CL': '17.0'
}


# ligands whose efp fragments should be splitted, (typically becase of size considerations)
# make sure to include rules for those cuts in define_ligand_cut() function
cut_ligands = ['CLA','LHG','LMG','BCL']   

# residues that require special EFP treatment (different makefp setup etc). Develop as needed
fancy_ligands = ['SF4']

# non-zero charges in ligands that will be split into several fragments
cut_ligand_charges = {'LHG':['tail', -1]}  # LGH tail charge -1

# non-zero charges in non-fragmented ligands 
ligand_charges = {}

# Dictionary of amino acid charges. Assumes any N-terminal is +1 and any C-terminal is -1.
AA_charge = {
    'ASP': -1, 'GLU': -1, 'HIP': 1, 
    'LYS': 1, 'ARG': 1, 'HISH': 1
}

# Amino acids with non-zero charge
charged_AAs = ['ASP', 'GLU', 'HIP', 'LYS', 'ARG', 'HISH']

# Map full amino acid names to one-letter codes for efp filenames
amino_acid_dict = {
    'ALA': 'a', 'ARG': 'r', 'ASN': 'n', 'ASP': 'd', 'CYS': 'c', 
    'GLN': 'q', 'GLU': 'e', 'GLY': 'g', 'HIS': 'h', 'ILE': 'i', 
    'LEU': 'l', 'LYS': 'k', 'MET': 'm', 'PHE': 'f', 'PRO': 'p', 
    'SER': 's', 'THR': 't', 'TRP': 'w', 'TYR': 'y', 'VAL': 'v', 
    'HIP': 'hp', 'HID': 'hd', 'HIE': 'he', 'NVAL': 'v', 'CGLN': 'q',
    'HISE':'he', 'HISD':'hd','HISH':'hp','ACYS':'acys'
}

# distance in Angstroms for adding capping Hydrogens
CH_distance = 1.07886
SH_distance = 1.34

# ---------------------------
# Function Definitions
# ---------------------------

# Makefp header. Works for most cases. 
def makefp_header(charge, comment):
    return (
        f" $contrl units=angs local=boys runtyp=makefp \n"
        f"       mult=1 icharg={charge} coord=cart icut=11 $end\n"
        f" $system timlim=99999 mwords=200 $end\n"
        f" $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06 $end\n"
        f" $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n"
        f" $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n"
        f" $MAKEFP POL=.t. DISP=.f. CHTR=.f. EXREP=.f. $end\n"
        f" $data\n {comment}\n C1\n"
    )

# Makefp header. Works for SF4 ligand. 
def makefp_special_header(comment):
    return (
        f" $contrl units=angs runtyp=makefp maxit=200 local=boys pp=sbkjc\n"
        f"      mult=1 icharg=-2 coord=cart icut=11 $end\n"
        f" $local maxloc=600 cvgloc=1.0E-6 $end\n"
        f" $system timlim=99999 memddi=0 mwords=200 $end\n"
        f" $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-05 fdiff=.f.\n"
        f"   shift=.t.  $end\n"
        f" $basis gbasis=sbkjc ndfunc=1 $end\n"
        f" $cphf ncphf=200 $end\n"
        f" $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n"
        f" $MAKEFP  POL=.t. DISP=.f. CHTR=.f.  EXREP=.f. $end\n"
        f" $data\n"
        f" sf4_1_{comment}\n"
        f" C1\n"
    )


def define_ligand_cut(RESNAME):
    Rings = []
    tailside = ''
    headside = ''

    # Settings for fragment splitting. 
    #CLA. the bond between atoms C5 and C6 is where the fragments split.
    if RESNAME=='CLA':
        tailside = 'CAA'             # Atom name defining the tail side boundary
        headside = 'C2A'             # Atom name defining the head side boundary
        Rings = [
        'mg1','h11','nc1','MG','CHA','CHB','H1','CHC','H2','CHD','H3','NA','C1A','C2A','H4',
        'C3A','H5','C4A','CMA','H6','H7','H8','NB','C1B','C2B','C3B','C4B',
        'CMB','H13','H14','H15','CAB','H16','CBB','H17','H18','NC','C1C',
        'C2C','C3C','C4C','CMC','H19','H20','H21','CAC','H22','H23','CBC',
        'H24','H25','H26','ND','C1D','C2D','C3D','C4D','CMD','H27','H28',
        'H29','CAD','OBD','CBD','H30','CGD','O1D','O2D','CED','H31','H32','H33'
        ]

    # BCL
    elif RESNAME=='BCL':
        tailside = 'CAA'             # Atom name defining the tail side boundary
        headside = 'C2A'             # Atom name defining the head side boundary
        Rings = ['MG','CHA','CHB','HB','CHC','HC','CHD','HD','NA','C1A','C2A','H2A','C3A','H3A','C4A',
       'CMA','HMA1','HMA2','HMA3','NB','C1B','C2B','C3B','C4B','CMB','HMB1','HMB2','HMB3',
       'CAB','OBB','CBB','HBB1','HBB2','HBB3','NC','C1C','C2C','H2C','C3C','H3C','C4C',
       'CMC','HMC1','HMC2','HMC3','CAC','HAC1','HAC2','CBC','HBC1','HBC2','HBC3','ND',
       'C1D','C2D','C3D','C4D','CMD','HMD1','HMD2','HMD3','CAD','OBD','CBD','HBD','CGD',
       'O1D','O2D','CED','HED1','HED2','HED3']

    
    #LMG
    elif RESNAME=='LMG':
        tailside = 'C11'             # Atom name defining the tail side boundary
        headside = 'C12'             # Atom name defining the head side boundary
        Rings = ['C12','H19','H20','C13','H21','H22','C14','C15','C16','C17','C18','C19','C20',
              'C21','C22','C23','C24','C25','C26','C27','H23','H24','H25','H26','H27','H28',
              'H29','H30','H31','H32','H33','H34','H35','H36','H37','H38','H39','H40','H41',
              'H42','H43','H44','H45','H46','H47','H48','H49','H50','H51'
        ]
    #LHG
    elif RESNAME=='LHG':
        tailside = 'C8'             # Atom name defining the tail side boundary
        headside = 'C9'             # Atom name defining the head side boundary
        Rings = ['C9','C10','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24',
             'H15','H16','H17','H18','H21','H22','H23','H24','H25','H26','H27','H28','H29','H30',
             'H31','H32','H33','H34','H35','H36','H37','H38','H39','H40','H41','H42','H43','H44',
             'H45'
        ]
    else:
        print('Wrong resname given to define_ligand_cut. Exit')
        sys.exit()
    return headside, tailside, Rings


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
                if len(line.split()) > 3:
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


def cut_frag(head, tail, distance):
    """
    Calculate coordinates for a virutal hydrogen  using a desired bond distance.
    
    Parameters:
        head: Line containing the head atom information.
        tail: Line containing the tail atom information.
    
    Returns:
        A list with the virtual coordinate [x, y, z] for the H-atom.
    """
    desired_dist = distance   # Desired bond distance
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

def make_standard_inp(fragment, res_type):
    """
    Generate standard makefp input 

    Parameters
    ----------
    fragment : list of atom lines in g96 format
    res_type : string with residue type

    Returns
    -------
    None.

    """
    if len(fragment) < 3:
        print(f'Fragment {fragment[0]} has less than 3 atoms. Standard makefp input will not work well. Take care!')
        sys.exit()
    
    charge = 0
    filename = ''

    qm_atoms = []
    mm_atoms = []
    H_caps = 0
    for atom in fragment:
        if 'H000' in atom:
            H_caps += 1
            continue
        if int(atom.split()[3]) in qm_IDs:
            qm_atoms.append(atom.split()[2])
        if int(atom.split()[3]) in MM_remove:
            mm_atoms.append(atom.split()[2])
        
    if H_caps + len(qm_atoms) == len(fragment):
        print(f'Fragment {fragment[0]} is completely in QM region. We will not build its makefp file. Continue.')
        return
    
    # Determine the charge based on the residue type.
    resname = fragment[0].split()[1]
    if res_type == 'ligand':    # normal ligand 
        if resname in ligand_charges.keys():
            charge = ligand_charges[resname]
    elif res_type == 'head' or res_type == 'tail':    # cut ligand 
        if resname in cut_ligand_charges.keys():
            if cut_ligand_charges[resname][0] == res_type:
                charge = cut_ligand_charges[resname][1]
    # taking care of AA charge. Account for charged AAs, terminal AAs and combinations
    elif res_type in ['aa', 'naa', 'caa']:   
        if resname[0] == 'N' and resname[1:] in known_amino_acids :  # N terminal AA
            charge += 1
        if resname[0] == 'C' and resname[1:] in known_amino_acids :  # C terminal AA
            charge -= 1
        if resname in charged_AAs :
            charge += AA_charge[resname]
        if resname[1:] in charged_AAs:
            charge += AA_charge[resname[1:]]
    else:
        print(f'Do not know how to determine charge for residue {resname}! Assume charge 0. Check!')


    # Construct the filename based on residue information. 
    #Use 1 letter identifier for Amino acids, otherwise filename will match residue name
    if resname in known_amino_acids:
        filename = amino_acid_dict[resname] + '_' + fragment[2].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
    # adding c or n to terminal AAs
    elif resname[1:] in known_amino_acids:
        filename = resname[0].lower() + amino_acid_dict[resname[1:]] + '_' + fragment[2].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
    elif res_type == 'ligand':
        filename = resname.lower() + '_' + fragment[0].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
    elif res_type == 'head':
        filename = resname.lower() + 'h_' + fragment[0].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
    elif res_type == 'tail':
        filename = resname.lower() + 't_' + fragment[0].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
    else :
        filename = resname.lower() + '_' + fragment[0].split()[0] + '_' + fragment[0].split()[3] + '_' + timestamp + '_'+mutant
        print(f'Do not know how to determine makefp filename for residue {resname}! Filename {filename} given.')

    
    # add fragname to frag_list for future use
    frag_list.append(filename)
    
    # add fragment to dictionary with fragment makefp file naame and first three coordiante lines for easy creating (QM/)EFP inputs
    frag_list_info.update(make_frag_info(filename, fragment, qm_atoms))

    filename=filename+'.inp'

    lines_out = []
    # Build header sections for the input file.
    lines_out.append(makefp_header(charge, filename.split('.')[0]))

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
            if 'MG' in atom.split()[2].upper():
                col2 = at_sym['MG']
            elif 'FE' in atom.split()[2].upper():
                col2 = at_sym['FE']
            else:
                col2 = at_sym[atom.split()[2][0].upper()]
            # Get coordinates, scaling by 10. Then format
            x, y, z = [float(atom.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(17)
            col4 = f"{y:.8f}".rjust(18)
            col5 = f"{z:.8f}".rjust(18)
            lines_out.append(f"{col1}{col2}{col3}{col4}{col5}\n")

    # Append comments regarding atoms to be removed later.
    lines_out.append(" $end \n!comment atoms to be erased:")
    for atomname in qm_atoms:
        lines_out.append(" " + str(atomname))
    lines_out.append("\n")
    # Append comments regarding atoms that will have polarization removed.
    lines_out.append("!polarization points to remove:")

    # Write the constructed lines to the file.
    with open(filename, 'w') as outfile:
        outfile.writelines(lines_out)
    

# compares the distance between two atoms given in g96 format
def pos_match(line1, line2):
    match = False
    part1 = line1.split()
    part2 = line2.split()
    # 4,5,6 columns
    pos1 = np.array([float(part1[4]), float(part1[5]), float(part1[6])])
    pos2 = np.array([float(part2[4]), float(part2[5]), float(part2[6])])
    dist = np.linalg.norm(pos2-pos1)
    # do not use shift for the first match
    if dist < 0.01 :
        match = True
    return match
    
# compares the resname and resid between two atoms given in g96 format
def residue_match(line1, line2):
    match = False
    part1 = line1.split()
    part2 = line2.split()
    # 0 and 1 columns
    if part1[0] == part2[0] and part1[1] == part2[1]:
        match = True
    return match

def read_ligand_topol(resname):
    #print(f'Read topology of ligand {resname}')
    start=0
    out_charges=[]
    with open(PATH_TO_AMBER +resname+'.itp','r') as itp:
        itp_lines=itp.readlines()
    for line in itp_lines:
        if start==1 and len(line.split()) >= 8:
            if ';' not in line.split()[0]:
                out_charges.append(float(line.split()[6]))
            #out_dict[line.split()[0]]=line.split()[6]
        if '[ atoms ]' in line:
            start=1
        if '[ bonds ]' in line:
            break
    return out_charges

def read_ion_water_charges(resname):
    if resname == 'SOL':
        o_charge = WATER_AND_IONS['SOL']
        charges = [o_charge, -o_charge/2, -o_charge/2]
    elif resname in WATER_AND_IONS:
        charges = [WATER_AND_IONS[resname]]
    else:
        print(f'Unknown resname {resname} in function read_ion_water_charges. Make new charge entry in WATER_AND_IONS!!')
        sys.exit()
    return charges
        
def residue_type(resname):
    if not isinstance(resname, str):
        print(f"Resname {resname} is not a string. Error.")
        sys.exit()
    res_type = 'none'
    if resname in skip_resnames:
        res_type = 'skip'
        return res_type
    if [aa for aa in known_amino_acids if aa in resname]:   # this line should take care for normal and terminal AAs
        if resname[0] == 'N' and resname[1:] in known_amino_acids :  # N terminal AA
            res_type = 'naa'
        elif resname[0] == 'C' and resname[1:] in known_amino_acids :  # C terminal AA
            res_type = 'caa'
        elif resname in known_amino_acids :
            res_type = 'aa'
        else: 
            res_type = 'ligand'
            print(f'Residue {resname} is assigned as ligand. Check whether this is corrcet.')
    elif resname in cut_ligands:
        res_type = 'cut_ligand'
    elif resname in fancy_ligands:
        res_type = 'fancy'
    elif resname in waters:
        res_type = 'water'
    elif resname in ions:
        res_type = 'ion'
    else: 
        res_type = 'ligand'
    
    return res_type
    
def update_charges(line1, line2, charge):
    charge1 = float(line1.split()[-1]) + charge/2
    charge2 = float(line2.split()[-1]) + charge/2
    #dcharge = (charge1 + charge2)/2
    part1 = line1.rsplit(' ', 1)
    if len(part1) > 1:
        new_line1 = part1[0] + ' ' + str(charge1)  # assume first atom is C and its charge is larger (positive) than charge of O
    else:
        print(f'Trouble splitting {line1}')
        sys.exit()
    part2 = line2.rsplit(' ', 1)
    if len(part2) > 1:
        new_line2 = part2[0] + ' ' + str(charge2) 
    else:
        print(f'Trouble splitting {line2}')
        sys.exit()
    return new_line1, new_line2

def update_charge(line1, charge):
    charge1 = float(line1.split()[-1]) + charge
    part1 = line1.rsplit(' ', 1)
    if len(part1) > 1:
        new_line1 = part1[0] + ' ' + str(charge1)  # assume first atom is C and its charge is larger (positive) than charge of O
    else:
        print(f'Trouble splitting {line1}')
        sys.exit()
    return new_line1

def print_qchem_qm_lines(qm_info_file, g96_lines):
    """
    Process the QM input file lines to generate QM coordinates.
    
    The function scans through the QM file until it finds a line with 'boundary'.
    For lines in the QM_atoms section (after 'QM_atoms' is encountered),
    it converts coordinate values and formats them.
    
    Parameters:
        qm_lines (list of str): Lines from the user deifned text file.
        g96_lines (list of str): Lines from the cleaned g96 list (with charges and no headers.
    
    Returns:
        A list of formatted coordinate lines.
    """

    qm_lines = []
    with open(qm_info_file, 'r') as file:
        qm_lines = file.readlines()

    bridge_IDs=[]
    qm_IDs=[]
    bridge_pairs = []
    start='none'
    for line in qm_lines:
        if 'QM_atoms' in line:
            start = 'QM'
            continue
        # After 'QM_atoms' is encountered, process lines with sufficient columns.
        if start=='QM' and (len(line.split()) > 4):
            # Format the atom label:
            qm_IDs.append(int(line.split()[3]))
        # When we reach the boundary marker, finish the QM section, find bonds to cap.
        if 'QM-MM boundary boundary' in line:
            continue
        if 'boundary' in line:
            start='boundary'
            continue
        if start=='boundary':
            bridge_IDs.append(line.split()[3])
            #print(bridge_IDs)
            if(len(bridge_IDs)%2==0):
                bridge_pairs.append([int(bridge_IDs[-2]), int(bridge_IDs[-1])])
                start='none'
    
    #print(bridge_pairs)
    # assume g96 atomid is continious and unique!!!
    g96_qm_lines = []
    for atom in qm_IDs:
        if atom == int(g96_lines[atom-1].split()[3]):
            g96_qm_lines.append(g96_lines[atom-1])
        else:
            print(f'{atom} does not match g96_lines[atom-1]')

    # format QM atoms for Q-Chem
    # converting from nm to Angstroms
    outlines = []         
    for line in g96_qm_lines:
        # Format the atom label:
        parts = line.split()
        # Use full label if the atom is in exceptions, otherwise take the first letter.
        if parts[2] in ATOM_EXCEPTIONS.keys():
            atom_name = ATOM_EXCEPTIONS[parts[2]]
        else:
            atom_name = parts[2][0].upper()
        outlines.append(f'{atom_name}    {float(parts[4])*10:.8f}   {float(parts[5])*10:.8f}   {float(parts[6])*10:.8f}')

    
    for pair in bridge_pairs:
        qm_line = ''
        mm_line = ''
        if pair[0] == int(g96_lines[pair[0]-1].split()[3]):
            qm_line = g96_lines[pair[0]-1]
        else:
            print(f'{pair[0]} does not match g96_lines[pair[0]-1]')
        if pair[1] == int(g96_lines[pair[1]-1].split()[3]):
            mm_line = g96_lines[pair[1]-1]
        else:
            print(f'{pair[1]} does not match g96_lines[pair[1]-1]')
        #print('Bridge')
        #print(qm_line)
        #print(mm_line)
        virt_coords=cut_frag(qm_line,mm_line,CH_distance)
        outlines.append(f'H    {virt_coords[0]:.8f}   {virt_coords[1]:.8f}   {virt_coords[2]:.8f}')

    return outlines

def make_frag_info(filename, frag, qm_atoms):
    lines = []
    lines.append(filename)
    counter = 0
    for i in range(len(frag)):
        if frag[i].split()[2] not in qm_atoms:
            atom_name = 'A0' + str(i+1) + frag[i].split()[2]
            x,y,z = [float(frag[i].split()[j])*10.0 for j in range(4, 7)]
            line = f'{atom_name}     {x:.8f}    {y:.8f}    {z:.8f}'
            lines.append(line)
            counter += 1
        if counter == 3:
            break
    return {filename:lines}

# ---------------------------
# Main Script Processing
# ---------------------------

# Read the EFP file lines
with open(efp_g96, 'r') as f:
    efp_lines = f.readlines()

# Read the full configuration file lines
with open(full_g96, 'r') as f:
    full_lines = f.readlines()


# Determine QM region atoms and MM removal atoms from the settings file.
qm_IDs, MM_remove = qm_atoms(settings_file)

# storage for EFP fragment names and first three coordinates
frag_list_info = {}

# ---------------------------
# Process EFP Region Residues
# ---------------------------
# Get a list of unique residue IDs first lines from the EFP file.
#   Note that atom IDs do not match between EFP and full structures. However, residue IDs do match

prev_resid = ''
prev_resname = ''
efp_residues = []
start = False
for line in efp_lines:
    if start:
        if 'END' in line:
            break
        
        if line.split()[0] != prev_resid or line.split()[1] != prev_resname:
            if residue_type(line.split()[1]) != 'skip':
                efp_residues.append(line)
            prev_resid = line.split()[0]
            prev_resname = line.split()[1]
               
    if 'POSITION' in line:
        start = True

##################
######## Processing full g96 and topology files to make a new list containing both atomids and charges
##################
# Read the topology file 
topol_lines = []
atomid_charge_dict = {}
with open(topol_file, 'r') as f:
    start = 0
    while f:
        line = f.readline()
        if '[ atoms ]' in line:
            start = 1
            continue
        if '[ bonds ]' in line:
            break
        if start == 1 and len(line.split()) >= 8 :
            if ';' not in line.split()[0]:
                topol_lines.append(line)
                atomid_charge_dict.update({int(line.split()[5]) : float(line.split()[6])})
    
# making g96_charges with combined info of (stripped from top and bottom) g96 and MM charges
g96_charges = []
start = 0
counter = 0
charges = []
resid_old = -1

for line in full_lines:
    if 'POSITION' in line:   #start processing
        start = 1
        continue
    if start == 1:
        if 'END' in line:
            break
        # working with protein topology
        if int(line.split()[3]) <= len(atomid_charge_dict):  # process protein topolgy lines. Assume that g96 atom ids are in order
            atomid = int(line.split()[3])
            charge = atomid_charge_dict[atomid]
            g96_charges.append(line.rstrip() + '  ' + str(charge))
        # working with cofactor topologies
        else:
            if line.split()[1] not in special_cofactors and line.split()[1] not in WATER_AND_IONS:   # Error
                print(f'Error! cofactor {line.split()[1]} is unknown')
                print(line)
                sys.exit()
            elif resid_old != line.split()[0]:
                if line.split()[1] in special_cofactors:
                    charges = read_ligand_topol(line.split()[1])
                elif line.split()[1] in WATER_AND_IONS:
                    charges = read_ion_water_charges(line.split()[1])
                counter = 0
                resid_old = line.split()[0]
            g96_charges.append(line.rstrip() + '  ' + str(charges[counter]))
            counter += 1

###################
            
# working with g96_charges list now                
mm_lines = []  # not EFP atoms
efp_lines = []  # atoms that will be in EFP region
efp_counter = 0

resid_old = ''
resname_old = ''
efp_residue = 'none' # 'none', 'aa', 'naa', 'caa', 'fancy','ligand', 'cut_ligand', 'water', 'ion'
CO = []
CA = []
efp_residue_old = 'none'

# for cut fragments
headside = ''
tailside = ''
ring = []
head = []
tail = []
head_cut = ''
tail_cut = ''

frag = []
efp_water_list = []  # list to store EFP waters
efp_ion_list = []  # list to store EFP ions if any
frag_list = []

accumulated_charge = 0.0
CO_charge = 0.0
CO_charge_next = 0.0
CA_index = -1
CA_index_next  = -1
CA_charge_index_list = []

fancy_list = []

for j in range(len(g96_charges)):
    line = g96_charges[j]
    
    # check whether resid changed. If it did, do all work on (possibly) fragment that just ended 
    if resid_old != line.split()[0] or resname_old != line.split()[1]:
        # take care of previous residue
        if efp_residue == 'aa':
            # take care of CO things 
                        
            vH2 = cut_frag(CA[-1], frag[-2], CH_distance)      #Cut between current CA and next C 
            vH1 = cut_frag(CO[0], CA[-2], CH_distance)         #Cut between current C and previous CA
            frag.pop()
            frag.pop()
            
            CO_charge_next = float(g96_charges[j-1].split()[-1]) + float(g96_charges[j-2].split()[-1])
            #line1, line2 = update_charges(g96_charges[j-2], g96_charges[j-1], 0.0)
            mm_lines.append(g96_charges[j-2])
            mm_lines.append(g96_charges[j-1])
            
            accumulated_charge += float(g96_charges[j-1].split()[-1]) + float(g96_charges[j-2].split()[-1])
            
            for co in CO:
                frag.append(co)
            frag.append(f" H000 1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
            frag.append(f" H000 1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
            
        elif efp_residue == 'naa':  # N-terminal residue
            vH2 = cut_frag(CA[-1], frag[-2], CH_distance)
            frag.pop()
            frag.pop()
            
            CO_charge_next = float(g96_charges[j-1].split()[-1]) + float(g96_charges[j-2].split()[-1])
            #line1, line2 = update_charges(g96_charges[j-2], g96_charges[j-1], 0.0)
            mm_lines.append(g96_charges[j-2])
            mm_lines.append(g96_charges[j-1])
            #mm_lines.append(g96_charges[j-2])
            #mm_lines.append(g96_charges[j-1])
            accumulated_charge += float(g96_charges[j-1].split()[-1]) + float(g96_charges[j-2].split()[-1])
            
            frag.append(f" H000    1.0      {vH2[0]:.8f}       {vH2[1]:.8f}       {vH2[2]:.8f}\n")
            
        elif efp_residue == 'caa':  # C-terminal residue
            vH1 = cut_frag(CO[0], CA[-2], CH_distance)
            for co in CO:
                frag.append(co)
            frag.append(f" H000    1.0      {vH1[0]:.8f}       {vH1[1]:.8f}       {vH1[2]:.8f}\n")
            
        elif efp_residue == 'cut_ligand': 
            # separate in head and tail
            head = []
            tail = []
            tail_cut = ''
            head_cut = ''
                       
            for atom in frag:
                # Get the current atom name, add line to the correct fragment
                atomname = atom.split()[2]
            
                # Append the line into the desired fragment based on the Rings list.
                if atomname in ring:
                    head.append(atom)
                else:
                    tail.append(atom)
                # Save the line for later use if it matches the headside or tailside.
                if atomname == tailside:
                    tail_cut = atom
                elif atomname == headside:
                    head_cut = atom

            if not head_cut:
                print(f'head_cut atom is not found in special ligand {efp_residues[efp_counter-1]}')
                sys.exit()
            if not tail_cut:
                print(f'tail_cut atom is not found in special ligand {efp_residues[efp_counter-1]}')
                sys.exit()
            head_coord = cut_frag(head_cut, tail_cut, CH_distance)
            tail_coord = cut_frag(tail_cut, head_cut, CH_distance)
            # Append the virtual hydrogen lines to the corresponding fragment lists.
            head.append(
                " H000 1.0" +
                f"{head_coord[0]:.8f}".rjust(17) +
                f"{head_coord[1]:.8f}".rjust(18) +
                f"{head_coord[2]:.8f}".rjust(18) + "\n"
            )
            tail.append(
                " H000 1.0" +
                f"{tail_coord[0]:.8f}".rjust(17) +
                f"{tail_coord[1]:.8f}".rjust(18) +
                f"{tail_coord[2]:.8f}".rjust(18) + "\n"
            )
        
        # won't make makefp inputs for water and ions, but count and save this info for later
        elif efp_residue == 'water':
            QM_water = False
            for atom in frag:
                if int(atom.split()[3]) in qm_IDs:
                    print(f'Found QM water in EFP list {atom}, continue...')
                    QM_water = True
                    break
            if not QM_water:      # append in EFP water list only EFP waters which fragments need to be in EFP region
                efp_water_list.append(frag)

        elif efp_residue == 'ion':
            if not int(frag[0].split()[3]) in qm_IDs:   # only EFP ions in the list
                efp_ion_list.append(frag)
        
        # now send prepared fragment lists to build makefp inputs
        if efp_residue in ['aa', 'naa', 'caa', 'ligand']:    # normal case
            make_standard_inp(frag, efp_residue)
        elif efp_residue == 'cut_ligand':
            make_standard_inp(head, 'head')
            make_standard_inp(tail, 'tail')
        elif efp_residue == 'fancy':
            print(f'Fancy frag {efp_residues[efp_counter-1]}')
            fancy_list.append([frag[0].split()[1], frag])
        if efp_residue != 'none':
            efp_lines.append(frag)
            
########    

        # now do work to prepare a new fragment if any    
        efp_residue_old = efp_residue
        efp_residue = 'none'
        resid_old = line.split()[0]
        resname_old = line.split()[1]
        frag = []
        CO_charge = 0.0
        
        res_type = residue_type(line.split()[1])
        if res_type != 'skip':            
            if efp_counter < len(efp_residues):  # check that we have any efp fragments left to search for
                if pos_match(line, efp_residues[efp_counter]):  # found matching EFP fragment
                    
                    efp_residue = res_type   
                    
                    # aminoacid
                    if efp_residue in ['aa','caa','naa']:
                        CO = []
                    if efp_residue in ['aa', 'caa']:
                        CO.append(g96_charges[j-2])
                        CO.append(g96_charges[j-1])
                        CO_charge = float(g96_charges[j-2].split()[-1]) + float(g96_charges[j-1].split()[-1])
                        # pop from mm_charges
                        mm_lines.pop()
                        mm_lines.pop()
                        accumulated_charge -= float(g96_charges[j-1].split()[-1]) + float(g96_charges[j-2].split()[-1])
                        # add charge to previous MM residue
                        if efp_residue_old == 'none' : 
                            CA_charge_index_list.append([CA_index, CO_charge])
                    
                    # cut ligand
                    if efp_residue == 'cut_ligand':
                        headside, tailside, ring = define_ligand_cut(line.split()[1])
                        
                    # no special treatment for other types
                    efp_counter += 1
                         
#  done with new-residue condition   
#  take care of saving lines into fragments
  
    if efp_residue != 'none' and efp_residue != 'skip':
        frag.append(line)
    else:
        if efp_residue != 'skip':
            mm_lines.append(line)
        
    # maybe add another check
    if line.split()[2] == 'CA':
        CA.append(line)
        
        # this part takes care of (integer)izing charge on classical AAs that were near EFP AAs. The charge of CO groups left from EFP AA fragments 
        # is added/subtracted to the Ca of the nearest classical AA 
        if efp_residue == 'none':   # prepare for the case if the next residue will be EFP amino acid
            CA_index = len(mm_lines) - 1
        if efp_residue_old in ['aa', 'naa'] and efp_residue == 'none' :  # this is the case when the previous residue was EFP AA, but the present one is MM AA
            CA_index_next = len(mm_lines) - 1
            CA_charge_index_list.append([CA_index_next, -CO_charge_next])
            CO_charge_next = 0
                

###
# Updating charges in mm_lines to make total charge of MM region integer and opposite of total charge of (EFP + QM) region

dcharge = 0.0
for i in CA_charge_index_list:
    dcharge += i[1]
print(f'Accumulated delta charge due to CO groups {dcharge}. We will fix it next.')
    
j = 0
for i in range(len(mm_lines)):
    for line in CA_charge_index_list:
        if i == line[0]:
            mm_lines[i] = update_charge(mm_lines[i], line[1])
        
# take care of 'fancy' residues now
# this is very system specific and non-standard
for ligand in fancy_list:
    if ligand[0] == 'SF4':        
        print('Updating SF4 and ACYS fragments.')
        # part 1. Update SF4 fragment by adding ACYS residue groups (CB, HB1, HB2, SG) and H link atoms
        acys_lines=[]
        sf4_lines=[]
        resids=['574','556','565','583']
        atoms=['CB','HB1','HB2','SG']
        bridges=[]
        at_start = -1

        for line in g96_charges:
            if line.split()[1]=='ACYS' and line.split()[0] in resids:
                if line.split()[2] in atoms:
                    acys_lines.append(line)
                elif line.split()[2]=='CA':
                    bridges.append(line)
            elif line.split()[1]=='SF4' and line.split()[0]=='1':
                if(line.split()[2]=='FE1'):
                    at_start=int(line.split()[3])
                sf4_lines.append(line)

        virt_Hs=[]
        i=0
        for line in acys_lines:
            if 'CB' in line:
                head=line
                tail=bridges[i]
                virt_coords=cut_frag(head,tail, CH_distance)
                virt_Hs.append(virt_coords)
                i+=1

        outlines=[]
        outlines.append(makefp_special_header(at_start))

        for line in sf4_lines:
            col1 = f" {line.split()[2]}".ljust(6)
            if line.split()[2][0].upper() == 'F':
                col2 = at_sym['FE'].ljust(4)
            else:
                col2 = at_sym[line.split()[2][0]].ljust(4)
            x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(16)
            col4 = f"{y:.8f}".rjust(15)
            col5 = f"{z:.8f}".rjust(15)
            outlines.append(f"{col1}{col2}{col3}{col4}{col5}\n")
        for line in acys_lines:
            col1 = f" {line.split()[2]}".ljust(6)
            col2 = at_sym[line.split()[2][0]].ljust(4)
            x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(16)
            col4 = f"{y:.8f}".rjust(15)
            col5 = f"{z:.8f}".rjust(15)
            outlines.append(f"{col1}{col2}{col3}{col4}{col5}\n")
        for line in virt_Hs:
            col1=' H000 '
            col2 = '1.0 '
            col3 = f"{line[0]:.8f}".rjust(16)
            col4 = f"{line[1]:.8f}".rjust(15)
            col5 = f"{line[2]:.8f}".rjust(15)
            outlines.append(f"{col1}{col2}{col3}{col4}{col5}\n")
            
        outlines.append(' $end\n')
        outlines.append('!comment atoms to be erased:\n')
        outlines.append('!polarization points to remove:\n')

        filename ='sf4_1_'+str(at_start)+'_'+timestamp+'_'+mutant
        frag_list.append(filename)
        
        # adding to frag_list_info
        lines = []
        lines.append(filename)
        counter = 0
        start = 0
        for line in outlines:
            if 'FE1' in line:
                start = 1
            if start == 1:
                if 'end' in line or 'END' in line:
                    start = 0
                    break
                else:
                    atom_name = 'A0' + str(counter+1) + line.split()[0]
                    x,y,z = [float(line.split()[j]) for j in range(2, 5)]
                    ll = f'{atom_name}     {x:.8f}    {y:.8f}    {z:.8f}'
                    lines.append(ll)
                    counter += 1
            if counter == 3:
                break
        frag_list_info.update({filename:lines})

        filename=filename+'.inp'

        with open(filename,'w') as out:
            for line in outlines:
                out.write(line)
        
        # Part 2. Updating ACYS fragments. Removing cysteine residue part (CB, HB1, HB2, SG) that were added to SF4 fragment. 
        # Adding capping hydrogen between CA and CB

        for filename in os.listdir('.'):
            if filename.startswith('acys') and filename.endswith('.inp'):
                outlines=[]
                head = ''
                tail = ''
                with open(filename,'r') as acys:
                    acys_lines=acys.readlines()
                start = 0
                for line in acys_lines:
                    if line.strip() == 'C1':
                        start = 1
                    if start == 1:
                        if 'CA' in line:
                            head = line
                        if 'CB' in line:
                            tail = line                            
                        if line.split()[0] in ['CB', 'HB1', 'HB2', 'SG']:
                            continue
                        if 'end' in line or 'END' in line:
                            break
                    outlines.append(line)
                    
                desired_dist = CH_distance  
                xh, yh, zh = [float(head.split()[i]) for i in range(2, 5)] 
                xt, yt, zt = [float(tail.split()[i]) for i in range(2, 5)] 
                dist_mag = np.sqrt((xh - xt)**2 + (yh - yt)**2 + (zh - zt)**2)
                virt = [
                    ((xt - xh) * desired_dist / dist_mag) + xh,
                    ((yt - yh) * desired_dist / dist_mag) + yh,
                    ((zt - zh) * desired_dist / dist_mag) + zh
                ]

                outlines.append(f" H000 1.0      {virt[0]:.8f}       {virt[1]:.8f}       {virt[2]:.8f}\n")
                outlines.append(' $end\n')
                outlines.append('!comment atoms to be erased:\n')
                outlines.append('!polarization points to remove:\n')
                                
                outfile=filename #DANGEROUS. Rewriting inp file for acys
                with open(outfile,'w') as out:
                    for line in outlines:
                        out.write(line)

    else:
        print(f'Do not know what to do about fancy fragment {ligand}. Take care!')
    

#############
# ANALYSIS

print('\n The main work is done. Some analysis now. Check it carefully!') 
 
total_charge = 0
for line in g96_charges:
    total_charge += float(line.split()[-1])
print(f'Total charge {total_charge}')

mm_charge = 0
for line in mm_lines:
    mm_charge += float(line.split()[-1])
print(f'MM charge {mm_charge}')


### checking charge of EFP region based on prepared makefp input files 
lines = []
for fname in os.listdir('.'):
    if fname.endswith('.inp'):
        with open(fname, 'r') as f:
            for line in f:
                if 'icharg' in line:
                    lines.append(fname + ':' + line)
                    break

efp_charge = 0
charged_list = []
for line in lines:
    parts = line.split()[2]
    res_charge = int(parts.split('=')[1])
    if res_charge != 0:
        charged_list.append(line)
    efp_charge += res_charge


print(f'Total charge of EFP region based on prepared makefp inputs: {efp_charge}')


print('\n Prepare files needed for building Q-Chem or libefp inputs... ')

with open('frag_list.txt','w') as file:
    for elem in frag_list:
        file.write(elem+'\n')

# $prot.efp file
with open('prot.efp', 'w') as file:
    print(' Writing prot.efp file.')
    file.write(' $prot\n')
    file.write(' Generated by make_AAs.py\n')
    file.write('  COORDINATES (BOHR)\n')
    # convert coordinates to Bohr units
    # make name as resid_resname_atom_name_atomid
    for line in mm_lines:
        parts = line.split()
        if len(parts) != 8:
            print(f'problem with mm_lines line {line}! Exit.')
            sys.exit()
        name = parts[0] + '_' + parts[1] + '_' + parts[2] + '_' + parts[3]
        file.write(f'{name}   {float(parts[4])*nm2Bohr:6.10}   {float(parts[5])*nm2Bohr:6.10}   {float(parts[6])*nm2Bohr:6.10}   0.0001  0.0\n')
    file.write(' STOP\n')
    # charges 
    file.write(' MONOPOLES\n')
    for line in mm_lines:
        parts = line.split()
        if len(parts) != 8:
            print(f'problem with mm_lines line {line}! Exit.')
            sys.exit()
        name = parts[0] + '_' + parts[1] + '_' + parts[2] + '_' + parts[3]
        file.write(f'{name}   {float(parts[7]):2.6}   0.0\n')
    file.write(' STOP\n')
    file.write(' $end\n')
    
# writing prot fragment info to file
with open('prot.txt','w') as file:
    print(' Writing prot.txt file.')
    file.write('prot\n')
    # convert coordinates to Angstrom units
    # make name as resid_resname_atom_name_atomid
    counter = 0
    for line in mm_lines:
        if counter == 3:
            break
        parts = line.split()
        if len(parts) != 8:
            print(f'problem with mm_lines line {line}! Exit.')
            sys.exit()
        name = parts[0] + '_' + parts[1] + '_' + parts[2] + '_' + parts[3]
        file.write(f'{name}   {float(parts[4])*10:.8f}   {float(parts[5])*10:.8f}   {float(parts[6])*10:.8f}\n')
        counter += 1
    
    
# QM region for Q-Chem input
# define charge and multiplicity of the QM region
with open('qm_region_qchem.txt', 'w') as file:
    print(' Writing QM region file.')
    qm_lines = print_qchem_qm_lines(settings_file, g96_charges)
    for line in qm_lines:
        file.write(line + '\n')
        
# prepare EFP water fragments
# Assume that efp water file uses TIP3P atom names, A01OW, A02HW1, A03HW2
# This is Q-Chem format
# adjust as needed
with open('water_efp_region.txt', 'w') as file:
    print(' Writing EFP water fragments into file.')
    for water in efp_water_list:
        ox,oy,oz = [float(water[0].split()[i])*10.0 for i in range(4, 7)] 
        h1x,h1y,h1z = [float(water[1].split()[i])*10.0 for i in range(4, 7)] 
        h2x,h2y,h2z = [float(water[2].split()[i])*10.0 for i in range(4, 7)] 
        file.write('water\n')
        file.write(f'A01OW    {ox:.8f}   {oy:.8f}   {oz:.8f}\n')
        file.write(f'A02HW1   {h1x:.8f}   {h1y:.8f}   {h1z:.8f}\n')
        file.write(f'A03HW2   {h2x:.8f}   {h2y:.8f}   {h2z:.8f}\n')
        
# prepare EFP fragment output, in Q-Chem format
with open('efp_region.txt', 'w') as file:
    print(' Writing EFP fragments into file.')
    for frag in frag_list_info.values():
        for line in frag:
            file.write(line + '\n')
        
    
    
    

    
    


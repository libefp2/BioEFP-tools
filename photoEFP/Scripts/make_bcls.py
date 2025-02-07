# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 14:42:41 2024

@author: jackl

This script reads in a .g96 file of the EFP region and a full configuration .g96 file,
and then creates .inp files for bacteriochlorophyll a molecules. In this process, the
head and tail groups are treated as separate fragments. 

This is a "hard-coded" script, that is, this script is very specific to chlorophylls 
or bacteriochlorophylls (specifically chloropyll a and bacteriochlorophyll a). 
If you intend to separate a molecule into multiplefragment .efp files, you will 
need to edit this script.
"""

import numpy as np
import sys

# ---------------------------
# Input files and cutting definitions
# ---------------------------
# Command-line arguments:
# Sample execution: python make_bcls.py efp_opt_83855.g96 optimized_83855.g96
g96_file = sys.argv[1]        # EFP region file (e.g., "efp_opt_83855.g96")
full_g96_file = sys.argv[2]   # Full configuration file (e.g., "optimized_83855.g96")

# Settings for fragment splitting. the bond between atoms C5 and C6 is where the fragments split.
tailside = 'C6'             # Atom name defining the tail side boundary
headside = 'C5'             # Atom name defining the head side boundary

# List of QM residue IDs. If you want all CLA/BCL residues to be made into fragments, make this empty
site = ['752', '754', '753', '790', '792', '793']  
oldres = '751'              # Initial residue number for tracking changes

# ---------------------------
# Function Definitions
# ---------------------------
def cut_frag(head, tail):
    """    
    This function uses the coordinates from two atoms (head and tail) to compute
    the positions for adding virtual hydrogens to cap the now-separated fragments
    (along the vector that is the C5-C6 bond).
    
    Parameters:
        head: A string line with the head atom lines.
        tail: A string line with the tail atom lines.
    
    Returns:
        h_t: coordinates [x, y, z] for the virtual hydrogen on the head side.
        t_h: coordinates [x, y, z] for the virtual hydrogen on the tail side.
    """
    desired_dist = 1.07886  # Desired bond distance
    # Scale coordinates by 10 (nm -> angstrom)
    xh, yh, zh = [float(head.split()[i]) * 10 for i in range(4, 7)]
    xt, yt, zt = [float(tail.split()[i]) * 10 for i in range(4, 7)]
    # Calculate the magnitude of the distance vector between head and tail.
    dist_mag = np.sqrt((xh - xt)**2 + (yh - yt)**2 + (zh - zt)**2)
    # Compute the new coordinates along the head-to-tail vector.
    h_t = [
        ((xt - xh) * desired_dist / dist_mag) + xh,
        ((yt - yh) * desired_dist / dist_mag) + yh,
        ((zt - zh) * desired_dist / dist_mag) + zh
    ]
    # Compute the new coordinates along the tail-to-head vector.
    t_h = [
        ((xh - xt) * desired_dist / dist_mag) + xt,
        ((yh - yt) * desired_dist / dist_mag) + yt,
        ((zh - zt) * desired_dist / dist_mag) + zt
    ]
    return h_t, t_h

def make_inp(fragment):
    """
    Generate an input (.inp) file for a fragment of bacteriochlorophyll a.
    The filename is based on the fragment's residue number and atom ID.
    The file includes various control sections and a coordinate list.
    
    Parameters:
        fragment: Lines of atoms comprising the fragment.
    """
    txt = []
    # Extract residue number and atom ID from the first line of the fragment.
    resnum = fragment[0].split()[0]
    atomid = fragment[0].split()[3]
    
    # Choose a filename prefix based on whether the first fragment line contains 'MG'
    if 'MG' in fragment[0]:
        namestart = 'clah_'
    else:
        namestart = 'clat_'
    filename = namestart + resnum + '_' + atomid + '.inp'
    
    # Construct the header for the input file.
    # You may want to customize these paramaters.
    header = (
        " $contrl units=angs local=boys runtyp=makefp \n"
        "       mult=1 icharg=0 coord=cart icut=11 $end\n"
        " $system timlim=99999   mwords=4000 $end\n"
        " $scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06  $end\n"
        " $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n"
        " $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n"
        " $MAKEFP  POL=.t. DISP=.f. CHTR=.f.  EXREP=.f. $end\n"
        " $data\n"
        " " + filename.split('.')[0] + "\n"
        " C1\n"
    )
    txt.append(header)
    
    # Append each atom line.
    for atom in fragment:
        if 'H000' in atom:
            # Virtual hydrogen atoms are taken as is.
            txt.append(atom)
        else:
            # Format the lines.
            col1 = f" {atom.split()[2]}".ljust(6)
            # Determine the atomic number using the at_sym dictionary.
            #    WARNING: Atoms with >1 letter symbols must be added if you have them.
            if atom.split()[2][0] == 'M':
                col2 = at_sym['MG']
            else:
                col2 = at_sym[atom.split()[2][0]]
            # nm -> angstrom.
            x, y, z = [float(atom.split()[i]) * 10 for i in range(4, 7)]
            col3 = f"{x:.8f}".rjust(17)
            col4 = f"{y:.8f}".rjust(18)
            col5 = f"{z:.8f}".rjust(18)
            txt.append(f"{col1}{col2}{col3}{col4}{col5}\n")
    
    # Append comments regarding atoms to be removed.
    txt.append(" $end \n!comment atoms to be erased:")
    # The QM region might include the entire chlorophyll or just the headring.
    # In case of only headring, a tail fragment must still be created, and
    # C6 is an MM atom covalently bound to C5, a QM atom.
    if fragment[0].split()[0] in site:
        txt.append(" C6")
    txt.append("\n")
    txt.append("!polarization points to remove:")
    if fragment[0].split()[0] in site:
        #H42, H43, C7 are bound to C6, polarizability must be romoved
        txt.append(" H42 H43 C7")
    
    # Write the output to a file.
    with open(filename, 'w') as outfile:
        outfile.writelines(txt)

# ---------------------------
# Global Dictionaries and Lists
# ---------------------------
# Dictionary to format atomic symbols to numbers (for output).
at_sym = {
    'H': '1.0', 'C': '6.0', 'N': '7.0', 'O': '8.0', 'MG': '12.0',
    'P': '15.0', 'S': '16.0', 'FE': '26.0', 'NA': '11.0', 'CL': '17.0'
}

# List of headring atoms. List includes names for YB's BCL and Reppert group's CLA.
Rings = [
    'MG','CHA','CHB','CHC','CHD','NA','C1A','C2A','C3A','C4A','CMA','NB',
    'C1B','C2B','C3B','C4B','CMB','CAB','CBB','NC','C1C','C2C','C3C','C4C',
    'CMC','CAC','CBC','ND','C1D','C2D','C3D','C4D','CMD','CAD','OBD','CBD',
    'CGD','O1D','O2D','CED','CAA','CBA','CGA','O1A','O2A','C1','C2','C3','C4',
    'C5','H1','H2','H3','H4','H5','H6','H7','H8','H9','H10','H11','H12','H13',
    'H14','H15','H16','H17','H18','H19','H20','H21','H22','H23','H24','H25',
    'H26','H27','H28','H29','H30','H31','H32','H33','H34','H35','H36','H37',
    'H38','H39','H40','H41'
]

# ---------------------------
# File Reading and Fragment Identification
# ---------------------------
# Read the EFP region file.
with open(g96_file, 'r') as efp:
    efp_lines = efp.readlines()

# Read the full configuration file.
with open(full_g96_file, 'r') as g96:
    full_lines = g96.readlines()

# Identify residue numbers (from lines containing 'CLA') in the EFP file that are not in the 'site' list.
frag_cla_resnums = []
for line in efp_lines:
    if 'CLA' in line:
        res_num = line.split()[0]
        if res_num in frag_cla_resnums:
            continue
        elif res_num not in site:
            frag_cla_resnums.append(res_num)

# ---------------------------
# Fragment Processing for Head and Tail Groups
# ---------------------------
# 'tailside' and 'headside' are the atoms where bonds are to be broken to create
# separate fragment files.
curr_head = []   # Append lines for the "head" fragment
curr_tail = []   # append lines for the "tail" fragment
head_cut = ''    # Temporary storage for the head atom used for virtual bond creation

# Process each line in the full configuration file.
for line in full_lines:
    if 'CLA' in line:
        # Process only for fragments with residue numbers in our list.
        if line.split()[0] in frag_cla_resnums and (oldres != line.split()[0]):
            # If no head_cut has been found, update oldres and continue.
            if len(head_cut) < 3:
                oldres = line.split()[0]
                continue
            # Once both a head and tail line have been recorded, compute virtual hydrogen positions.
            head_coord, tail_coord = cut_frag(head_cut, tail_cut)
            # Append the virtual hydrogen lines to the corresponding fragment lists.
            curr_head.append(
                f" H000 1.0" +
                f"{head_coord[0]:.8f}".rjust(17) +
                f"{head_coord[1]:.8f}".rjust(18) +
                f"{head_coord[2]:.8f}".rjust(18) + "\n"
            )
            curr_tail.append(
                f" H000 1.0" +
                f"{tail_coord[0]:.8f}".rjust(17) +
                f"{tail_coord[1]:.8f}".rjust(18) +
                f"{tail_coord[2]:.8f}".rjust(18) + "\n"
            )
            # Reset head_cut for the next fragment.
            head_cut = ''
            # Depending on whether the residue is in the 'site' list, write both fragments or only the tail.
            if oldres not in site:
                make_inp(curr_head)
                make_inp(curr_tail)
            else:
                make_inp(curr_tail)
            # Reset the current fragment accumulators.
            curr_head = []
            curr_tail = []
            oldres = line.split()[0]
        
        # Get the current atom name.
        atomname = line.split()[2]
        # Append the line into the desired fragment based on the Rings list.
        if line.split()[0] in frag_cla_resnums:
            if atomname in Rings:
                curr_head.append(line)
            else:
                curr_tail.append(line)
            # Save the line for later use if it matches the headside or tailside.
            if atomname == tailside:
                tail_cut = line
            elif atomname == headside:
                head_cut = line

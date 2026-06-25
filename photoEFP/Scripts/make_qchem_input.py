# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:33:41 2024

@author: lyuda

Sample execution:
    python make_qchem_input.py qchem_file_name
    
The script will search for the following files in the current directory:
    qm_region_qchem.txt 
    efp_region.txt 
    water_efp_region.txt 
    prot.txt

"""

import sys
import os
import numpy as np

if len(sys.argv) != 2:
    print('Provide qchem file name')
    sys.exit()
    
qchem_file_name = sys.argv[1]

## Check charge, multiplicity, header, n_states (number of excited states)

charge = 0
multiplicity = 1
n_states = 5


# adjust as needed!
def build_header(n_states):
    """
    Return a list of header lines (as strings) that will be prepended to the output.
    Hard coded, you may want to adjust these parameters.    
    """
    header_lines = [
'$rem\n',
'JOBTYPE SP\n',
'METHOD wB97x\n',
'BASIS 6-31G*\n',
'purecart = 1111\n',
'MAX_SCF_CYCLES=100\n',
'SCF_ALGORITHM = diis\n',
f'cis_n_roots {n_states}\n',
'cis_singlets true\n',
'cis_triplets false\n',
'max_cis_cycles 100\n',
'CIS_CONVERGENCE = 7\n',
'rpa 2\n',
'SYM_IGNORE TRUE\n',
'mem_total 250000  ! 32 proc on lslipche, 2 GB per processor\n',
'mem_static 1000\n',
'gui 2\n',
'EFP_COORD_XYZ          TRUE\n',
'EFP                    TRUE\n',
'EFP_FRAGMENTS_ONLY     FALSE\n',
'EFP_DISP               FALSE\n',
'EFP_EXREP              FALSE\n',
'EFP_QM_DISP            FALSE\n',
'EFP_QM_EXREP           FALSE\n',
'EFP_POL TRUE\n',
'EFP_QM_POL TRUE\n',
'efp_pairwise = 0\n',
'efp_print = 2\n',
'efp_pol_field_update = 3\n',
'EFP_POL_DAMP=1\n',
'EFP_POL_DAMP_TT_VALUE=300\n',
'   MAKE_CUBE_FILES  true   ! triggers writing of cube files\n',
'   PLOTS            true\n',
'$end\n',
'\n']
    return header_lines


def get_plot_section(qm_outlines, n_states):
    """
    Prepares $plots section basen on QM coordinates

    Parameters
    ----------
    qm_outlines : $molecule section

    Returns
    -------
    $plots section lines

    """
    outlines = []
    xmin = ymin = zmin = 100000
    xmax = ymax = zmax = -100000
    for line in qm_outlines:
        line = line.rsplit()
        if len(line) == 4:
            #print(line)
            if float(line[1]) < xmin: xmin = float(line[1])
            if float(line[1]) > xmax: xmax = float(line[1])
            if float(line[2]) < ymin: ymin = float(line[2])
            if float(line[2]) > ymax: ymax = float(line[2])
            if float(line[3]) < zmin: zmin = float(line[3])
            if float(line[3]) > zmax: zmax = float(line[3])
    
    outlines.append('$plots\n')
    outlines.extend(f' grid_range ({int(xmin) - 2},{int(xmax) + 2}) ({int(ymin) - 2},{int(ymax) + 2}) ({int(zmin) - 2},{int(zmax) + 2})\n')
    gridx = (int(xmax)-int(xmin)+4)*3
    gridy = (int(ymax)-int(ymin)+4)*3
    gridz = (int(zmax)-int(zmin)+4)*3
    outlines.extend(f' grid_points   {gridx}  {gridy}  {gridz}\n')
    outlines.extend(f' attachment_detachment_density  1-{n_states}\n')
    outlines.extend('$end\n\n')
    
    return outlines
            

#####################################
# main

#if len(sys.argv) != 4:
#    print('Please provide a full structure .g96 file, an EFP structure .g96 file, and a user_defined text file with QM atoms and QM-MM boundary atoms')
#    sys.exit()

qm_lines = []
efp_lines = []
water_lines = []
prot_lines = []

with open('qm_region_qchem.txt', 'r') as qm_file:
    try:  
        qm_lines = qm_file.readlines()
    except:
        print('Could not find qm_region_qchem.txt file. Exit.')
        sys.error()

with open('efp_region.txt', 'r') as efp_file:
    try:  
        efp_lines = efp_file.readlines()
    except:
        print('Could not find efp_region.txt file. Take care! Continue...')

with open('water_efp_region.txt', 'r') as water_file:
    try:  
        water_lines = water_file.readlines()
    except:
        print('Could not find water_efp_region.txt file. Take care! Continue...')

with open('prot.txt', 'r') as prot_file:
    try:  
        prot_lines = prot_file.readlines()
    except:
        print('Could not find water_efp_region.txt file. Take care! Continue...')


# Build the output list by starting with the header.
header_lines = build_header(n_states)
outlines = []
for line in header_lines:
    outlines.append(line)

# add qm section
outlines.append('$molecule\n')
outlines.append(f'{charge} {multiplicity}\n')
outlines.extend(qm_lines)
outlines.append('$end\n\n')


######    comment if not needed! ########
# prepare $plot section
plot_lines = get_plot_section(qm_lines, n_states)
outlines.extend(plot_lines)

# adding efp section and filling it with fragments, waters, prot fragment
outlines.append('\n$efp_fragments\n')
if efp_lines:
    outlines.extend(efp_lines)
if water_lines:
    outlines.extend(water_lines)
if prot_lines:
    outlines.extend(prot_lines)
outlines.append('$end\n\n')


'''
#THESE LINES BELOW ARE FOR EFP PAIRWISE!! When 2 jobs are run from 1 file
outlines.append('\n')
outlines.append('@@@\n')
outlines.append('\n')
# then repeat all input sections starting with modified header
'''

# Write the final output to a file.
#outname=g96_filename.replace('.g96','.in')
outname = qchem_file_name
print(f'Writing Q-Chem QM/EFP file {qchem_file_name}')
with open(outname, 'w') as f:
    for line in outlines:
        f.write(line)
    #f.write('$end')

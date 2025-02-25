# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:33:41 2024

@author: jackl

Sample execution:
    python cut_qm.py pair_83855.g96 efp_pair_83855.g96 qm_file.inp

This script reads in three files:
  - A structure file in .g96 format (e.g. pair_83855.g96)
  - An EFP structure file (e.g. efp_pair_83855.g96)
  - A user_defined text file (user_defined.txt) with QM atoms and QM-MM boundary atoms
  
It then builds an output file containing:
  - A header section (with parameters)
  - A modified coordinate section (with fragments marked for removal)
  - Information for water molecules and a conversion of coordinates from another file.
"""

import sys
import os

# List of atom symbols to treat as exceptions (printed as full two-letter symbol 
#    rather than first character)
ATOM_EXCEPTIONS = ['MG']

def build_header(efp_order):
    """
    Return a list of header lines (as strings) that will be prepended to the output.
    Hard coded, you may want to adjust these parameters.
    """
    header_lines = [
        '$rem\n',
        'JOBTYPE SP\n',
        'METHOD PBE0\n',
        'BASIS 6-31G*\n',
        'purecart = 1111\n',
        'xc_grid 000096000302\n',
        'SCF_ALGORITHM = diis\n',
        'cis_n_roots 2\n',
        'cis_singlets true\n',
        'cis_triplets false\n',
        'max_cis_cycles 100\n',
        'CIS_CONVERGENCE = 7\n',
        'rpa 2\n',
        'SYM_IGNORE TRUE\n',
        'mem_total 20000\n',
        'mem_static 5000\n',
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
        'efp_pairwise = 1\n',
        'efp_order = '+efp_order+' \n',
        '$end\n\n',
        '$molecule\n',
        '0 1\n'
    ]
    return header_lines

def read_files(g96_filename, efp_filename, qm_filename):
    """
    Read the contents of the three input files.
    
    Parameters:
        g96_filename (str): Name of the structure (.g96) file.
        efp_filename (str): Name of the reference EFP file.
        qm_filename (str): Name of the QM input file.
    
    Returns:
        Tuple of lists: (g96_lines, efp_lines, qm_lines)
    """
    with open(g96_filename, 'r') as f:
        g96_lines = f.readlines()
    with open(efp_filename, 'r') as f:
        efp_lines = f.readlines()
    with open(qm_filename, 'r') as f:
        qm_lines = f.readlines()
    return g96_lines, efp_lines, qm_lines

def build_efp_atom_lists():
    """
    Build two data structures from the current directory:
      - A dictionary mapping the starting atom number (as string) to the filename (without extension).
      - A list of atom numbers (as strings) for EFP atoms. For each found file, we add the starting atom
        and the two subsequent atom numbers.
    Returns:
        (efp_dict, efp_atoms)
    """
    efp_atoms = []
    efp_dict = {}
    for filename in os.listdir('./'):
        # Skip water and classical region files.
        if filename in ('water.efp', 'prot.efp'):
            continue
        if filename.startswith('cla'):
            continue
        if filename.endswith('.efp'):
            # Remove the extension and extract the starting atom number.
            fragname = filename.split('.')[0]
            efp_atom_start = int(fragname.split('_')[2])
            key = str(efp_atom_start)
            efp_dict[key] = fragname
            # For each fragment file, assume three atoms: starting number, +1 and +2.
            efp_atoms.extend([str(efp_atom_start), str(efp_atom_start + 1), str(efp_atom_start + 2)])
    return efp_dict, efp_atoms


def process_qm_file_lines(qm_lines):
    """
    Process the QM input file lines to generate QM coordinates.
    
    The function scans through the QM file until it finds a line with 'boundary'.
    For lines in the QM_atoms section (after 'QM_atoms' is encountered),
    it converts coordinate values and formats them.
    
    Parameters:
        qm_lines (list of str): Lines from the user deifned text file.
    
    Returns:
        A list of formatted coordinate lines.
    """
    outlines = []
    start = False
    for line in qm_lines:
        # When we reach the boundary marker, finish the QM section.
        if 'boundary' in line:
            break
        # After 'QM_atoms' is encountered, process lines with sufficient columns.
        elif start and (len(line.split()) > 4):
            # Format the atom label:
            tokens = line.split()
            # Use full label if the atom is in exceptions, otherwise take the first letter.
            if tokens[2] in ATOM_EXCEPTIONS:
                col1 = tokens[2].rjust(6)
            else:
                col1 = (tokens[2][0] + ' ').rjust(6)
            # (nm -> angstrom)
            x, y, z = [float(tokens[i]) * 10 for i in range(4, 7)]
            col2 = f"{x:.8f}".rjust(15)
            col3 = f"{y:.8f}".rjust(15)
            col4 = f"{z:.8f}".rjust(15)
            outlines.append(f"{col1}{col2}{col3}{col4}\n")
        elif 'QM_atoms' in line:
            start = True
    outlines.append('$end\n')
    outlines.append(' $efp_fragments\n')
    return outlines

def process_structure_coords(g96_lines, efp_atoms, efp_dict):
    """
    Process the structure (.g96) file lines to add EFP fragment coordinates.
    
    The function scans for the 'POSITION' marker and then grabs the first three coordinates.
    
    Parameters:
        g96_lines (list of str): Lines from the structure file.
        efp_atoms (list of str): List of EFP atom numbers (as strings) from the EFP files.
        efp_dict (dict): Mapping from starting atom number (as string) to fragment name.
    
    Returns:
        A list of formatted EFP coordinates.
    """
    outlines = []
    start = False
    atomcounter = 3  # Initialize counter (used to control how many atoms to output)
    for line in g96_lines:
        if start:
            if 'END' in line:
                break
            # For the first few atoms (if atomcounter < 3) add coordinate lines.
            elif atomcounter < 3:
                atomcounter += 1
                col1 = ('A0' + str(atomcounter) + line.split()[2]).ljust(8)
                x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
                col2 = f"{x:.8f}".rjust(13)
                col3 = f"{y:.8f}".rjust(13)
                col4 = f"{z:.8f}".rjust(13)
                outlines.append(f"{col1}{col2}{col3}{col4}\n")
            # If the fourth token of the line (index 3) is one of the EFP atoms...
            elif line.split()[3] in efp_atoms:
                # Retrieve the corresponding fragment name.
                frag_name = efp_dict[line.split()[3]]
                outlines.append(frag_name + '\n')
                col1 = ('A01' + line.split()[2]).ljust(8)
                x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
                col2 = f"{x:.8f}".rjust(13)
                col3 = f"{y:.8f}".rjust(13)
                col4 = f"{z:.8f}".rjust(13)
                outlines.append(f"{col1}{col2}{col3}{col4}\n")
                atomcounter = 1
        if 'POSITION' in line:
            start = True
    return outlines

def process_water_coords(efp_lines):
    """
    Process water coordinates from the EFP structure file.
    
    Looks for lines with 'SOL' and 'OW' to identify water oxygen,
    then appends corresponding water hydrogen coordinates.
    -only EFP waters are desired
    
    Parameters:
        efp_lines (list of str): Lines from the EFP (truncated) structure file.
    
    Returns:
        A list of formatted water coordinate lines.
    """
    outlines = []
    found_water = 0
    for line in efp_lines:
        if(len(line.split())<2):
            continue
        #if 'SOL   OW' in line:
        elif(line.split()[1]=='SOL' and line.split()[2]=='OW'):
            found_water = 2
            #Every water EFP fragment is called "water"; residue numbers are lost
            outlines.append('water\n')
            col1 = 'A01O1'.ljust(8)
            x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
            col2 = f"{x:.8f}".rjust(13)
            col3 = f"{y:.8f}".rjust(13)
            col4 = f"{z:.8f}".rjust(13)
            outlines.append(f"{col1}{col2}{col3}{col4}\n")
        elif found_water > 0:
            # For subsequent lines, assign hydrogen labels.
            col1 = ('A01H' + str(3 - found_water)).ljust(8)
            x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
            col2 = f"{x:.8f}".rjust(13)
            col3 = f"{y:.8f}".rjust(13)
            col4 = f"{z:.8f}".rjust(13)
            outlines.append(f"{col1}{col2}{col3}{col4}\n")
            found_water -= 1
    return outlines

def process_prot(prot):
    """
    Process lines from the classical efp file (prot.efp) to grab coordinates.
    
    For each coordinate line, convert bohr -> angstrom
    
    Parameters:
        prot (str): Path to the prot.efp file.
    
    Returns:
        A list of formatted converted coordinate lines.
    """
    outlines = []
    start = -1
    with open(prot, 'r') as f:
        test_lines = f.readlines()
    for line in test_lines:
        # Once 'COORDINATES' is encountered, set start counter to 3.
        if start > 0:
            col1 = line.split()[0].ljust(8)
            # Convert coordinates from Bohr to Angstroms
            x, y, z = [float(line.split()[i]) * 0.529177249 for i in range(1, 4)]
            col2 = f"{x:.8f}".rjust(13)
            col3 = f"{y:.8f}".rjust(13)
            col4 = f"{z:.8f}".rjust(13)
            outlines.append(f"{col1}{col2}{col3}{col4}\n")
            start -= 1
        elif start == 0:
            break
        elif 'COORDINATES' in line:
            start = 3
    return outlines

def main(g96_filename,efp_filename,qm_filename):
    # Read input files.
    g96_lines, efp_lines, qm_lines = read_files(g96_filename, efp_filename, qm_filename)
    
    # Build the list/dictionary of EFP atom numbers from all .efp files in the current directory.
    efp_dict, efp_atoms = build_efp_atom_lists()

    # Build the output list by starting with the header.
    header_lines = build_header('1')
    outlines = []
    for line in header_lines:
        outlines.append(line)

    # Process QM input file lines to create coordinate lines for the fragment.
    qm_outlines = process_qm_file_lines(qm_lines)
    outlines.extend(qm_outlines)
    
    # Process structure coordinates from the g96 file.
    g96_coords = process_structure_coords(g96_lines, efp_atoms, efp_dict)
    outlines.extend(g96_coords)
    
    # Process water molecule coordinates from the structure file.
    water_coords = process_water_coords(efp_lines)
    outlines.extend(water_coords)
    
    # Process classical fragment file 'prot.efp' and append converted coordinate lines.
    test_mm_lines = process_prot('prot.efp')
    outlines.extend(test_mm_lines)
    
    # Append updated header lines (with efp_order updated) at the end.
    updated_header = build_header('2')
    for line in updated_header:
        outlines.append(line)
    
    # Write the final output to a file.
    with open('test_file', 'w') as f:
        for line in outlines:
            f.write(line)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3])
    #g96_file, efp_g96_file, user_defined.txt

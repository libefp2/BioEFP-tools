# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:33:41 2024

@author: jackl

Sample execution:
    python make_final.py efp_pair_83855.g96 pair_83855.g96 user_defined.txt

This script reads in three files:
  - An EFP structure file (e.g. efp_pair_83855.g96)
  - A structure file in .g96 format (e.g. pair_83855.g96)
  - A user_defined text file (user_defined.txt) with QM atoms and QM-MM boundary atoms
  
It then builds an output file containing:
  - A header section (with parameters)
  - A modified coordinate section (with fragments marked for removal)
  - Information for water molecules and a conversion of coordinates from another file.
"""

import sys
import os
import numpy as np

# List of atom symbols to treat as exceptions (printed as full two-letter symbol 
#    rather than first character)
ATOM_EXCEPTIONS = ['MG']
RESNAME='bcl'

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
        (efp_dict, efp_atoms, efp_nums)
    """
    efp_atoms = {}
    efp_dict = {}
    efp_nums={}
    for filename in os.listdir('./'):
        # Skip water and classical region files.
        if filename in ('water.efp', 'prot.efp'):
            continue
        '''
        if filename.startswith(RESNAME):
            #Logic to include CLAs or BCLs into final
            
            continue
        '''
        if filename.endswith('.efp'):
            # Remove the extension and extract the starting atom number.
            fragname = filename.split('.')[0]
            efp_atom_start = int(fragname.split('_')[2])
            key = str(efp_atom_start)
            efp_dict[key] = fragname
            names=[]
            nums=[]
            with open(filename,'r') as fragment:
                for line in fragment.readlines():
                    if(line[0]=='A'):
                        atomnum=int(line[1:3])
                        nums.append(atomnum+efp_atom_start-1)
                        names.append(line.split()[0])
                    if(len(names)==3):
                        efp_atoms[key]=names
                        efp_nums[key]=nums
                        break
            # For each fragment file, assume three atoms: starting number, +1, and +2.
            #efp_atoms.extend([str(efp_atom_start), str(efp_atom_start + 1), str(efp_atom_start + 2)])
    return efp_dict, efp_atoms, efp_nums

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
    return h_t

def get_qm_lines(user_lines,g96):
    """
    Process the QM input file lines to generate QM coordinates.
    
    The function scans through the QM file until it finds a line with 'boundary'.
    For lines in the QM_atoms section (after 'QM_atoms' is encountered),
    it converts coordinate values and formats them.
    
    Parameters:
        qm_lines (list of str): Lines from the user deifned text file.
        g96_lines (list of str): Lines from the full .g96 file.
    
    Returns:
        A list of formatted coordinate lines.
    """
    bridge_IDs=[]
    bridge_pairs={}
    qm_IDs=[]
    qm_bridge=[]
    start=0
    for line in user_lines:
        # When we reach the boundary marker, finish the QM section, find bonds to cap.
        if 'boundary' in line:
            start=2
        elif start==2:
            bridge_IDs.append(line.split()[3])
            #print(bridge_IDs)
            if(len(bridge_IDs)%2==0):
                bridge_pairs[bridge_IDs[-2]]=bridge_IDs[-1]
                qm_bridge.append(bridge_IDs[-2])
                #bridge_IDs=[]
                start=0
        # After 'QM_atoms' is encountered, process lines with sufficient columns.
        elif start==1 and (len(line.split()) > 4):
            # Format the atom label:
            qm_IDs.append(line.split()[3])
        elif 'QM_atoms' in line:
            start = 1
    outlines=[]
    bridge_atoms=[]
    #bridge_lines=[]
    #print(bridge_IDs)
    for line in g96:
        if(len(line.split())<4):
            continue
        elif(line.split()[3] in qm_IDs):
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
        if(line.split()[3] in bridge_IDs):
            bridge_atoms.append(line)
            '''
            if(len(bridge_atoms)==2):
                virt_coords=cut_frag(bridge_atoms[0],bridge_atoms[1])
                col1='H '.rjust(6)
                col2 = f"{virt_coords[0]:.8f}".rjust(15)
                col3 = f"{virt_coords[1]:.8f}".rjust(15)
                col4 = f"{virt_coords[2]:.8f}".rjust(15)
                bridge_lines.append(f"{col1}{col2}{col3}{col4}\n")
                bridge_atoms=[]
                start=0
            '''
    #print(qm_bridge)
    #print(bridge_atoms)
    for qm in qm_bridge:
        for line in bridge_atoms:
            if line.split()[3]==qm:
                head=line
            elif line.split()[3]==bridge_pairs[qm]:
                tail=line
        virt_coords=cut_frag(head,tail)
        col1='H '.rjust(6)
        col2 = f"{virt_coords[0]:.8f}".rjust(15)
        col3 = f"{virt_coords[1]:.8f}".rjust(15)
        col4 = f"{virt_coords[2]:.8f}".rjust(15)
        outlines.append(f"{col1}{col2}{col3}{col4}\n")
    #for virt_H in bridge_lines:
    #    outlines.append(virt_H)
    outlines.append('$end\n')
    outlines.append(' $efp_fragments\n')
    return outlines

def process_structure_coords(g96_lines, efp_atoms, efp_dict, efp_nums):
    """
    Process the structure (.g96) file lines to add EFP fragment coordinates.
    
    The function scans for the 'POSITION' marker and then grabs the first three coordinates.
    
    Parameters:
        g96_lines (list of str): Lines from the structure file.
        efp_atoms (dict): Mapping from filename atom ID (string) to the .efp file atomnames.
        efp_dict (dict): Mapping from filename atom ID (string) to fragment full name.
        efp_nums (dict): Mapping from filename atom ID (string) to .g96 IDs for the 3 reference atoms.
    
    Returns:
        A list of formatted EFP coordinates.
    """
    outlines = []
    start = False
    search= False
    #atomcounter = 3  # Initialize counter (used to control how many atoms to output)
    for line in g96_lines:
        if start:
            # Search until "END"
            if 'END' in line:
                break
            #If current atom index is found in fragment file atom, get fragment name, atomname, atom IDs
            #    -note, generally, the first atom ID will match file name, HOWEVER, in cases where cut_qm.py 
            #     has removed the first atom, this will be offset
            elif line.split()[3] in efp_dict:
                #print(line)
                fragname=efp_dict[line.split()[3]]
                atomnames=efp_atoms[line.split()[3]]
                atomIDs=efp_nums[line.split()[3]]
                outlines.append(fragname + '\n')
                #print(atomIDs[0],line.split()[3])
                j=0
                search=True
            # Fragment indicator is found, now find specific atoms and coordinates
            if search:
                if int(line.split()[3]) == atomIDs[j]:
                    col1 = (atomnames[j]).ljust(8)
                    x, y, z = [float(line.split()[i]) * 10 for i in range(4, 7)]
                    col2 = f"{x:.8f}".rjust(13)
                    col3 = f"{y:.8f}".rjust(13)
                    col4 = f"{z:.8f}".rjust(13)
                    outlines.append(f"{col1}{col2}{col3}{col4}\n")
                    j+=1
                    if(j==3):
                        search = False
        if 'POSITION' in line:
            start = True
    return outlines

def water_dist(x,y,z,x2,y2,z2):
    dist=np.sqrt((x-x2)**2+(y-y2)**2+(z-z2)**2)
    return dist

def process_water_coords(efp_lines,user_lines):
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
    qm_sol=[]
    for line in user_lines:
        parts=line.split()
        if 'QM-MM' in line:
            break
        elif(len(parts)>6):
            if line.split()[1]=='SOL' and line.split()[2]=='OW':
                qm_sol.append([float(parts[4]),float(parts[5]),float(parts[6])])
    outlines = []
    found_water = 0
    for line in efp_lines:
        if(len(line.split())<2):
            continue
        #if 'SOL   OW' in line, and atom ID is not a QM atom:
        elif(line.split()[1]=='SOL' and line.split()[2]=='OW') and line.split()[3]:
            nearest_dist=100.0
            for atom in qm_sol:
                dist=water_dist(atom[0],atom[1],atom[2],float(line.split()[4]),float(line.split()[5]),float(line.split()[6]))
                if(nearest_dist>dist):
                    nearest_dist=dist
            if(nearest_dist>0.001):
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
            col1 = ('A0'+str(4-found_water)+'H' + str(3 - found_water)).ljust(8)
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
    efp_dict, efp_atoms, efp_nums = build_efp_atom_lists()

    # Build the output list by starting with the header.
    header_lines = build_header('1')
    outlines = []
    for line in header_lines:
        outlines.append(line)

    # Process QM input file lines to create coordinate lines for the fragment.
    qm_outlines = get_qm_lines(qm_lines,g96_lines)
    outlines.extend(qm_outlines)
    
    # Process structure coordinates from the g96 file.
    g96_coords = process_structure_coords(g96_lines, efp_atoms, efp_dict, efp_nums)
    outlines.extend(g96_coords)
    
    # Process water molecule coordinates from the structure file.
    water_coords = process_water_coords(efp_lines,qm_lines)
    outlines.extend(water_coords)
    
    # Process classical fragment file 'prot.efp' and append converted coordinate lines.
    test_mm_lines = process_prot('prot.efp')
    outlines.extend('prot\n')
    outlines.extend(test_mm_lines)
    
    outlines.append('$end\n')
    outlines.append('\n')
    outlines.append('@@@\n')
    outlines.append('\n')
    
    # Append updated header lines (with efp_order updated) at the end.
    updated_header = build_header('2')
    for line in updated_header:
        outlines.append(line)
    
    # Process QM input file lines to create coordinate lines for the fragment.
    outlines.extend(qm_outlines)
    
    # Process structure coordinates from the g96 file.
    outlines.extend(g96_coords)
    
    # Process water molecule coordinates from the structure file.
    outlines.extend(water_coords)
    
    # Process classical fragment file 'prot.efp' and append converted coordinate lines.
    outlines.extend('prot\n')
    outlines.extend(test_mm_lines)
    
    # Write the final output to a file.
    with open('test_file', 'w') as f:
        for line in outlines:
            f.write(line)
        f.write('$end')

if __name__ == "__main__":
    main(sys.argv[2],sys.argv[1],sys.argv[3])
    # g96_file, efp_g96_file, user_defined.txt

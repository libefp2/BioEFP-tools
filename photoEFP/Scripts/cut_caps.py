# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 07:42:51 2024

@author: jackl & lyuda

Sample execution:
    cut_caps.py makefp.inp makefp.efp mode
                       -- mode is 'check' (default) or 'fix'

The script prepares efp prameter file for efp or qmefp calculations. 
1. It checks efp file for errors and potential problems 
  - script exits with error if the efp file is not complete ($end statement is missing), 
  fragment charge is not integer, or partial charges exceed |3|)
  - warning is printed if large (>20) polarizability values are found
  
!!!!  Check/adjust charge_cutoff and polab_cutoff settings if needed !!!!
  
  Nothing more happens in mode 'check'
  
2. In 'fix' mode, the script continues by cutting off capping hydrogens and other atoms 
    following comments in inp file and cleans corresponding multipoles and polarizability 
    points. The residual charge of the fragment is corrected. 
    
    Additioanlly, POLAB section is added if large polarizablities are found.
    The resulting efp potential is saved in cut_{makefp.efp} file.
"""

import sys
import numpy as np
import itertools
import os



#########################################
### USER SETTINGS

# charge_cutoff is the value of a partial charge when the potential is 
# considered bad without hope (probably meaning that makefp job did not 
# run properly)
charge_cutoff = 4.0

# polab_cutoff is a used to assign polab values based on the values of polarizabilities
# { polarizabiloity_value : polab_value}
# provide as many entries as wanted
polab_cutoff = {20.0 : 0.2, 50.0 : 0.1, 100.0 : 0.05}

#########################################

if len(sys.argv) < 3:
    print('Provide 1) input file, 2) efp file, and optionally, mode of operation (check or fix, default: check). Exit')
    sys.exit()

inp = sys.argv[1]
efp = sys.argv[2]
if len(sys.argv) > 3:
    mode = sys.argv[3]
else: 
    mode = 'check'



def extract_atomnums_B(line):
    """
    Given a line starting with 'B', extract two atom numbers based on the length
    of the first token.
    """
    token = line.split()[0]
    if len(token) == 4:
        return int(token[2]), int(token[3])
    elif len(token) == 5:
        return int(token[2:4]), int(token[4])
    elif len(token) == 6:
        return int(token[2:4]), int(token[4:6])
    else:
        return None, None
         
def extract_bond_atoms(bond_name):
    """
    Given a line starting with 'B', extract two atom numbers based on the length
    of the first token.
    """
    token = bond_name
    if len(token) == 4:
        return int(token[2]), int(token[3])
    elif len(token) == 5:
        return int(token[2:4]), int(token[4])
    elif len(token) == 6:
        return int(token[2:4]), int(token[4:6])
    else:
        return None, None

def distance(x1, y1, z1, x2, y2, z2):
    dist = np.sqrt((float(x1) - float(x2))**2 +
                   (float(y1) - float(y2))**2 +
                   (float(z1) - float(z2))**2)
    return dist

def get_specials(lines):
    """
    Grab input file lines to obtain the atom IDs that should be removed entirely
    and the atom indexes for which only the polarizable points should be removed.
    
    The function looks for keywords 'erased:' and 'remove:' in the lines.
        -make_AAs_V2.py added these comments.
    
    Reads:
        lines: input file (.inp) lines.
    Returns:
        atomid (list): Atom indexes for complete removal.
    """
    atoms = []  # List for names of atoms to be removed entirely.
    
    # First, extract the atom names from lines that contain 'erased:' and 'remove:'
    for line in lines:
        j = -1
        if 'erased:' in line:
            # Walk backwards from the end of the line until reaching 'erased:'
            while line.split()[j] != 'erased:':
                atoms.append(line.split()[j])
                j -= 1

    # Scan through the file lines to grab the positions of the pesky atoms/polarizabilities.
    i = 0
    atomid = []  # List to record the atom index numbers for full removal.
    start = 0    # Flag to indicate when the coordinates section starts.
    for line in lines:
        if '$end' in line:
            start = 0
            continue
        if start == 1:
            i += 1  # Increment a counter for each coordinate line
            # Check the atom name (first token in the line)
            if line.split()[0] in atoms:
                atomid.append(i)
            # Additionally, if the atom label contains 'H000', mark it for removal
            #     -All virual hydrogens are to be removed.
            elif 'H000' in line:
                atomid.append(i)
        # The coordinates section is assumed to start when a line: 'C1'
        if line.strip() == 'C1':
            start = 1
    return atomid

def get_coords(lines, atoms):
    """
    Grab the coordinate lines while separating those that should be removed (rem_coords) 
    from those to be kept (coords).
    
    Reads:
        lines: Lines from the EFP file.
        atoms: Atom indexes to remove entire atom.
    
    Returns:
        coords, rem_coords; where:
          - coords: List of coordinate lines to keep.
          - rem_coords: List of coordinate lines that are flagged for removal.
    """
    start = False
    coords = []
    rem_coords = []

    # Process each line in the file.
    for line in lines:
        # Exit when "STOP" is encountered.
        if 'STOP' in line:
            return coords, rem_coords
        # Begin processing with the "COORDINATES" section.
        if 'COORDINATES' in line:
            start = True
            continue
        if not start:
            continue

        # If the line starts with 'A', extract the two-digit atom number.
        if line.startswith('A'):
            atomnum = int(line[1:3])
            if atomnum in atoms:
                rem_coords.append(line)
            else:
                coords.append(line)

        # For lines starting with 'B', use the helper to extract atom numbers.
        elif line.startswith('BO'):
            atomnum, atomnum2 = extract_atomnums_B(line)
            if atomnum is None:
                continue  # Skip if extraction failed
            # If either atom number is in the removal list, add to rem_coords.
            if atomnum in atoms or atomnum2 in atoms:
                rem_coords.append(line)
            else:
                coords.append(line)

    return coords, rem_coords

def get_monopoles(efp_lines, keep_coords, remove_coords):
    
    # extract lines containing monopoles
    start = 0
    mono_lines = []
    for line in efp_lines:
        if 'MONOPOLE' in line: 
            start = 1
            continue
        if start == 1:
            if 'STOP' in line:
                start = 0
                break
            else:
                mono_lines.append(line)

    # create master list with all important info
    monopoles = []
    for line in mono_lines:
        parts = line.split()
        if parts[0][0] == 'A': # atom
            monopoles.append(['atom', int(parts[0][1:3]), parts[0], float(parts[1]), float(parts[2])])
        elif parts[0][:2] == 'BO': # bond mid point
            atom1, atom2=extract_atomnums_B(line)
            monopoles.append(['bond', [atom1,atom2], parts[0], float(parts[1]), float(parts[2])])
        else:
            print(f'Could not assign line {line} to atom or bond. Error')
            sys.exit()
    
    keep_names = []
    for coord in keep_coords:
        keep_names.append(coord.split()[0])
    remove_names = []
    for coord in remove_coords:
        remove_names.append(coord.split()[0])

    # analyze names and create arrays connecting names and atom numers
    H_caps = []
    bonds = []
    non_caps = []
    keep_atoms = []
    for name in remove_names:
        if name[-4:] == 'H000': # capping H 
            H_caps.append(int(name[1:3]))  # real atom number
        elif name[:2] == 'BO': # bond mid point
            atom1, atom2=extract_bond_atoms(name)
            bonds.append([atom1, atom2])
        else:
            #print(f'non-capping atom to be removed {name}')
            non_caps.append(int(name[1:3]))
    for name in keep_names:
        if name[0] == 'A': # atom
            keep_atoms.append(int(name[1:3]))
            
    # points with charges already redistributed
    removed_charge_names = []
    
    # combine H_caps and corresponding bonds and compute those total charges
    for cap in H_caps:
        removed_charge_names_tmp = []
        cap_charge = 0.0
        for mon in monopoles:
            if cap == mon[1]:
                #print('cap',cap, mon)
                removed_charge_names_tmp.append(mon[2])
                cap_charge += mon[3] + mon[4]
                break
        # atom to add charge to 
        atom_to = -1
        for bond in bonds:
            if cap in bond:  # found corresponding bond
                for mon in monopoles:
                    if mon[0] == 'bond':
                        if bond == mon[1]:
                            #print('bond', bond, mon)
                            removed_charge_names_tmp.append(mon[2])
                            cap_charge += mon[3] + mon[4]
                            # atom_to is another number in the bond
                            if cap == bond[0]:
                                atom_to = bond[1]
                            else:
                                atom_to = bond[0]
                            #print(f'atom_to {atom_to}')
                            break
        # working with atom_to charge
        if atom_to == -1: # something went wrong
            print(f'did not find matching atom for H_cap {cap}. Error')
            sys.exit()
        else:
            if atom_to in keep_atoms:
                for mon in monopoles:
                    if mon[0] == 'atom':
                        if atom_to == mon[1]:
                            #print(f'adding charge {cap_charge} to {mon}')
                            # adding cap and bond names to removed_charge array only if there is a matching atoms
                            # to add charge to
                            for name in removed_charge_names_tmp:
                                removed_charge_names.append(name)
                            mon[3] += cap_charge
                            break
                        
            # in case the cap is connected to the atom that also will be removed,
            # simply add it to non_caps array and treat there
            elif atom_to in non_caps:
                non_caps.append(cap)
            else: 
                print(f'Could not find where to place a charge of H_cap {mon}. Error')
                sys.exit()
                
    # heavy atoms to be removed
    # two tasks: 1. compute total charge to be shifted
    # 2. find the atom to which to place it
    # start with step 2. 
    # There might be multiple bonds along which to-be-kept and 
    # to-be-removed regions are connected
    
    if non_caps:
        atoms_to = []
        for atom in keep_atoms:  # remaining atoms
            for bond in bonds:  # to be removed bonds
                if atom in bond: # found atoms that is connected to the one to be removed
                    atoms_to.append(atom)
                    #print(f'Charge will be placed to atom {atom_to}')
        if not atoms_to:
            print('did not find matching atom for deleted heavy atoms. Error')
            sys.exit()
        
        # total charge to be shifted
        charge = 0.0
        for name in remove_names:
            if name not in removed_charge_names:  # not worked with this monopoles during H_caps removal
                for mon in monopoles:
                    if name == mon[2]:
                        charge += mon[3] + mon[4]
        #print('charge', charge)
        
        # adding this charge to atoms in atoms_to // split charge equallyif multiple atoms are there
        charge /= len(atoms_to)
        for mon in monopoles:
            if mon[0] == 'atom':
                for atom in atoms_to:
                    if mon[1] == atom:
                        mon[3] += charge
                        #print(f'adding charge {charge} to {mon}')
                             
    # print out resulting info
    outlines = []
    for mon in monopoles:
        if mon[2] in keep_names:
            outlines.append(f'{mon[2]}   {mon[3]:0.10f}   {mon[4]:0.5f}\n')
    return outlines
            
            
def get_dipoles(lines, coords):
    """
    Extract the dipoles section lines that should be removed.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinate lines of atoms to remove.
    
    Returns:
        dipoles (list of str): Lines from the DIPOLES section matching the kept atoms.
    """
    dipoles = []
    keep_names = [atom.split()[0] for atom in coords]
    start = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return dipoles
            # append lines if the atom name is in the list of coordinates.
            if line.split()[0] in keep_names:
                dipoles.append(line)
        if 'DIPOLES' in line:
            start = 1
    return dipoles

def get_quadrupoles(lines, coords):
    """
    Extract the quadrupoles section lines.
    
    Uses counters (k and j) to control how many lines to skip or include
    when the atom name is encountered. Quardupoles are multi-line entries.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinate lines of atoms to remove.
    
    Returns:
        quadrupoles (list of str): Lines from the QUADRUPOLES section.
    """
    quadrupoles = []
    keep_names = [atom.split()[0] for atom in coords]
    start = 0
    k = 0  # Counter to skip one line after a match
    j = 0  # Counter to include one line after a non-match
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return quadrupoles
            elif k == 1:
                k -= 1
            elif j == 1:
                quadrupoles.append(line)
                j -= 1
            # When encountering a cut name, set k to skip the next line; do not append.
            elif line.split()[0] not in keep_names:
                k = 1
            else:
                quadrupoles.append(line)
                j = 1
        if 'QUADRUPOLES' in line:
            start = 1
    return quadrupoles

def get_octupoles(lines, coords):
    """
    Extract the octupoles section lines.
    
    Uses counters (j and k) similar to quadrupoles to control skipping/inclusion.
    -octupoles are also multi-lin entries.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinate lines of atoms to remove.
    
    Returns:
        octupoles (list of str): Lines from the OCTUPOLES section.
    """
    octupoles = []
    keep_names = [atom.split()[0] for atom in coords]
    start = 0
    j = 0
    k = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return octupoles
            elif k > 0:
                k -= 1
            elif j > 0:
                octupoles.append(line)
                j -= 1
            elif line.split()[0] not in keep_names:
                k = 2
            else:
                octupoles.append(line)
                j = 2
        if 'OCTUPOLES' in line:
            start = 1
    return octupoles


def get_polarpts(lines, cut_coords):
    """
    Extract the polarizable points section.
    
    For each polarizable point, compute the distance to each atom in the cut_coords.
    If the minimum distance is greater than cutoff, the point is retained.
    -Polarizable points do not have "atomnames" the same as the multipole parameters.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinates of atoms to remove.
    
    Returns:
        polars (list of str): Lines from the POLARIZABLE POINTS section that pass the filter.
    """
    
    # atoms to be removed 
    atom_cut_coords = []
    for line in cut_coords:
        if line.rsplit()[0][0] == 'A': # atom
            atom_cut_coords.append(line)

    # H000s
    Hs=[]
    for line in cut_coords:
        if 'H000' in line:
            Hs.append(line)
            
    polars = []
    #remove_names = [atom.split()[0] for atom in cut_coords]
    start = 0
    j = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return polars
            elif j > 0:
                polars.append(line)
                j -= 1
            # Check lines starting with "CT" (assumed polarizable point lines)
            elif 'CT' in line:
                j=3
##########################################################################
#########  ATTENTION: CUTOFFS to remove polarizability points  ###########
#########   different cutoofs for H000 and other removed atoms ###########
#########    want to remove nond and LP LMOs around the removed atom #####
#########    might need to adjust for heavier atoms                   ####
##########################################################################
                # severe cut: remove polarizability points on neigboring atoms
                for atom in atom_cut_coords:
                    #j=3
                    # Calculate the distance between the polarizable point and an atom in cut_coords.
                    current_dist = distance(line.split()[1], line.split()[2], line.split()[3],
                                            atom.split()[1], atom.split()[2], atom.split()[3])
                    if current_dist < 1.88973 * 1.6:  # in Bohr
                        j = 0
                        break
                
                # soft cut here: only delete the pol point at the nearest bond mid-point
                for atom in Hs:
                    current_dist = distance(line.split()[1], line.split()[2], line.split()[3],
                                            atom.split()[1], atom.split()[2], atom.split()[3])
                    if current_dist < 1.88973 * 0.8:  
                        j = 0
                        break

                if(j!=0):
                    polars.append(line)                
        if 'POLARIZABLE POINTS' in line:
            start = 1
    return polars

def get_screen(lines, coords, title):
    """
    Extract screening parameter lines from the file.
    
    Only lines with atom names in coords are kept.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        coords (list of str): Coordinates used to filter the screening parameters.
        title (str): A string that identifies the screen section (e.g., 'SCREEN ' or 'SCREEN2').
    
    Returns:
        params (list of str): Filtered screening parameter lines.
    """
    params = []
    keep_names = [atom.split()[0] for atom in coords]
    start = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return params
            if line.split()[0] in keep_names:
                params.append(line)
        if title in line:
            start = 1
    return params

def get_header(lines, input_name):
    """
    Construct the header for the output file by replacing the placeholder '$FRAGNAME'
    with the fragment name (derived from the input file name); otherwise retain the same
    heading lines as-is from the .efp file.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        input_name (str): The input file name.
    
    Returns:
        header (list of str): The modified header lines.
    """
    header = []
    fragname = input_name.split('.')[0]
    for line in lines:
        header.append(line.replace('$FRAGNAME', '$'+fragname))
        # Stop adding header lines once the COORDINATES section is encountered.
        if 'COORDINATES (BOHR)' in line:
            return header

def check_completeness(efp_lines):

    start = 0
    charges = []
    polarizabilities = []
    pol_name = ''
    error = 0
    
    # check monopoles
    for line in efp_lines:
        if 'MONOPOLES' in line:
            start = 1
            continue
        if 'STOP' in line:
            start = 0
            continue
        if start == 1:
            parts = line.rsplit()
            charges.append([float(parts[1]), float(parts[2])])
        if 'POLARIZABLE POINTS' in line:
            start = 2
            continue
        if start == 2:
            if 'CT' in line:
                parts = line.rsplit()
                pol_name = parts[0]
                pol = []
                continue
            else:
                parts = line.rsplit()
                if len(parts) == 5:
                    pol.append([ float(parts[0]), float(parts[1]),float(parts[2]),float(parts[3]) ])
                elif len(parts) == 1:
                    pol.append( [float(parts[0])] )
                    pol = list(itertools.chain.from_iterable(pol))
                    polarizabilities.append([pol_name, pol])
                    pol_name = ''
        if '$END' in line or '$end' in line:
            start = -1
    
    # analysis and info
    if start != -1:
        #print(f'ERROR: {efp} is incomplete!')
        error = 1
        
    return error, charges, polarizabilities

def check_charges(charges, max_charge):
    error = 0
    charges_np = np.array(charges)
    chg = 0.0
    total_charge = 0.0
    error_charge = 0.0
    for c in charges_np:
        ch = c[0] + c[1]
        if abs(ch) > chg:
            chg = abs(ch)
        total_charge += ch
    if abs(total_charge - round(total_charge)) > 0.001:
        #print(f'WARNING: {efp} has total non-integer charge {total_charge}')
        error_charge = abs(total_charge - round(total_charge))
        error = 1
    if chg > max_charge:
        #print(f'WARNING: {efp} has large partial charge {chg}')
        error_charge = chg
        error = 2
    return error, error_charge

def check_pol(polarizabilities, polab_cutoff):
    error = 0
    max_pol = -1
    polab = -1
    max_pol_final = -1
    
    if not polarizabilities:
        error = 1
    else:    
        for pol in polarizabilities:
            for p in pol[1]:
                if abs(p) > max_pol:
                    max_pol = abs(p)
            
        polab_cutoff_ordered = dict(sorted(polab_cutoff.items()))
    
        for key in polab_cutoff_ordered:
            if max_pol >= key:
                polab = polab_cutoff_ordered[key]
                max_pol_final = max_pol
                
    return error, max_pol_final, polab



# Read input files
with open(inp, 'r') as inp_:
    inp_lines = inp_.readlines()
with open(efp, 'r') as efp_:
    efp_lines = efp_.readlines()
    

# check efp potential for completeness and sanity
# Exit with error if potential is incomplete (does not have $end finisher)
# or if partial charges are very large

polab = -1 
max_pol = -1
error_complete, charges, polarizabilities = check_completeness(efp_lines)
error_elec, max_charge = check_charges(charges, charge_cutoff)
error_pol, max_pol, polab = check_pol(polarizabilities, polab_cutoff)

# conclusions based on the original assessment of the efp potential
if error_complete:
    error_inp = 'err_' + inp
    os.system(f'cp {inp} {error_inp}') 
    print(f'{efp} is incomplete. Error file {error_inp} is created.')
    sys.exit()
    
if error_elec: 
    error_inp = 'err_' + inp
    os.system(f'cp {inp} {error_inp}') 
    if error_elec == 1:
        print(f'{efp} has non-integer charge of {max_charge}. Error file {error_inp} is created.')
    elif error_elec == 2:
        print(f'{efp} has large partial charge of {max_charge}. Error file {error_inp} is created.')
    sys.exit()
    
if error_pol:
    error_inp = 'err_' + inp
    os.system(f'cp {inp} {error_inp}') 
    print(f'WARNING! {efp} has no polarizability section. Error file {error_inp} is created.')
    # do not exit here in case no pol section is intended
    
if max_pol > 0:
    print(f'WARNING! {efp} has large polarizability point of {max_pol}. POLAB value of {polab} will be added.')
    
# exit if not 'fix' mode
if mode != 'fix':
    sys.exit()
    
# changing the efp potential in 'fix' mode
else:                
    # Construct header by replacing '$FRAGNAME' with the actual fragment name.
    header = get_header(efp_lines, inp)
    
    # Get atom indexes for complete removal 
    # Polarizable point removals will happen based on removed atoms
    rem_atoms = get_specials(inp_lines)
    
    # Extract coordinate lines:
    # keep_coords: coordinate lines to be kept,
    # rem_coords: coordinate lines flagged for removal.
    keep_coords, rem_coords = get_coords(efp_lines, rem_atoms)
    
    # Extract sections for monopoles, dipoles, quadrupoles, octupoles, and polarizable points.
    keep_monop = get_monopoles(efp_lines, keep_coords, rem_coords)
    
    # run check_charges again on a new charge section
    charges_updated = []
    for ll in keep_monop:
        charges_updated.append([float(ll.split()[1]), float(ll.split()[2])])
    error_elec, max_charge = check_charges(charges_updated, charge_cutoff)
    if error_elec: 
        error_inp = 'err_' + inp
        os.system(f'cp {inp} {error_inp}') 
        if error_elec == 1:
            print(f'WARNING! Modified {efp} file (cut_{efp}) has non-integer charge of {max_charge}. Check it!')
        elif error_elec == 2:
            print(f'WARNING! Modified {efp} file (cut_{efp}) has large partial charge of {max_charge}. Check it!')
    
    keep_dip = get_dipoles(efp_lines, keep_coords)
    keep_quadrup = get_quadrupoles(efp_lines, keep_coords)
    keep_octup = get_octupoles(efp_lines, keep_coords)
    
    # Hidden here distance cutoffs to eliminate bond LMOs around atoms
    keep_pols = get_polarpts(efp_lines, rem_coords)
    
    # Extract screening parameters from two different screen sections.
    keep_screen = get_screen(efp_lines, keep_coords, 'SCREEN ')
    keep_screen2 = get_screen(efp_lines, keep_coords, 'SCREEN2')
    
    # Write the final output file with the appropriate sections.
    with open('cut_' + efp, 'w') as outfile:
        # Write header lines.
        for outline in header:
            outfile.write(outline)
        # Write coordinates.
        for outline in keep_coords:
            outfile.write(outline)
        # Write monopoles section.
        outfile.write(' STOP\nMONOPOLES\n')
        for outline in keep_monop:
            outfile.write(outline)
        # Write dipoles section.
        outfile.write(' STOP\nDIPOLES\n')
        for outline in keep_dip:
            outfile.write(outline)
        # Write quadrupoles section.
        outfile.write(' STOP\nQUADRUPOLES\n')
        for outline in keep_quadrup:
            outfile.write(outline)
        # Write octupoles section.
        outfile.write(' STOP\nOCTUPOLES\n')
        for outline in keep_octup:
            outfile.write(outline)
        # Write polarizable points section.
        outfile.write(' STOP\nPOLARIZABLE POINTS\n')
        for outline in keep_pols:
            outfile.write(outline)
        # Write screening sections.
        outfile.write(' STOP\nSCREEN2      (FROM VDWSCL=   0.700)\n')
        for outline in keep_screen2:
            outfile.write(outline)
        outfile.write('STOP\nSCREEN       (FROM VDWSCL=   0.700)\n')
        for outline in keep_screen:
            outfile.write(outline)
        outfile.write('STOP\n')
        if polab != -1:
            outfile.write(f'POLAB {polab:.2f}\nSTOP\n')
        outfile.write(' $END\n')
    

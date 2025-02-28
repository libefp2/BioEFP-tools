# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 07:42:51 2024

@author: jackl

Sample execution:
    cut_qm.py ala_33_473.inp a0001.efp

This file reads in an input file and a reference EFP file that should be a good geometric match (low RMSD),
AND the .inp file has listed atoms and polarizable points that should be deleted.
Atoms that will be added to the QM region should be removed from the fragment file.
Additionally, a bridge between a QM-atom and an MM-atom that are covalently bound should have the QM atom removed,
and some information of the MM-atom should be removed. The coordinates, monopole, and screen parameter sections
for the MM-atom are retained, while the dipoles, quadrupoles, octupoles, and polarizable point sections 
should be removed.
"""

import sys
import numpy as np

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
        polid (list): Atom indexes for polarizable points removal.
    """
    atoms = []  # List for names of atoms to be removed entirely.
    pols = []   # List for names of atoms whose polarizable points are to be removed.
    
    # First, extract the atom names from lines that contain 'erased:' and 'remove:'
    for line in lines:
        j = -1
        if 'erased:' in line:
            # Walk backwards from the end of the line until reaching 'erased:'
            while line.split()[j] != 'erased:':
                atoms.append(line.split()[j])
                j -= 1
        elif 'remove:' in line:
            # Walk backwards from the end of the line until reaching 'remove:'
            while line.split()[j] != 'remove:':
                pols.append(line.split()[j])
                j -= 1

    # Scan through the file lines to grab the positions of the pesky atoms/polarizabilities.
    i = 0
    atomid = []  # List to record the atom index numbers for full removal.
    polid = []   # List to record the atom index numbers for polarizable point removal.
    start = 0    # Flag to indicate when the coordinates section starts.
    for line in lines:
        if '$end' in line:
            start = 0
        elif start:
            i += 1  # Increment a counter for each coordinate line
            # Check the atom name (first token in the line)
            if line.split()[0] in atoms:
                atomid.append(i)
            elif line.split()[0] in pols:
                polid.append(i)
            # Additionally, if the atom label contains 'H000', mark it for removal
            #     -All virual hydrogens are to be removed.
            elif 'H000' in line:
                atomid.append(i)
        # The coordinates section is assumed to start when a line: 'C1'
        elif 'C1' == line.split()[0]:
            start = 1
    return atomid, polid

def get_coords(lines, atoms, pols):
    """
    Grab the coordinate lines while separating those that should be removed (rem_coords) 
    from those to be kept (coords).
    
    Reads:
        lines: Lines from the EFP file.
        atoms: Atom indexes to remove entire atom.
        pols: Atom indexes for which only polarizable points should be removed.
    
    Returns:
        coords, rem_coords; where:
          - coords: List of coordinate lines to keep.
          - rem_coords: List of coordinate lines that are flagged for removal.
    """
    start = False
    coords = []
    rem_coords = []
    
    # Helper function to extract two atom numbers from bond midpoints, bond identifier starts with 'B'
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
            elif atomnum in pols:
                coords.append(line)
                rem_coords.append(line)
            else:
                coords.append(line)

        # For lines starting with 'B', use the helper to extract atom numbers.
        elif line.startswith('B'):
            atomnum, atomnum2 = extract_atomnums_B(line)
            if atomnum is None:
                continue  # Skip if extraction failed
            # If either atom number is in the removal list, add to rem_coords.
            if atomnum in atoms or atomnum2 in atoms:
                rem_coords.append(line)
            # If either atom number is in the polarizable removal list, add to both lists.
            elif atomnum in pols or atomnum2 in pols:
                coords.append(line)
                rem_coords.append(line)
            else:
                coords.append(line)

    return coords, rem_coords

def get_monopoles(lines, coords):
    """
    Extract the monopole section for the atoms given coordinates.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        coords (list of str): List of coordinate lines to keep.
    
    Returns:
        monopoles (list of str): Lines from the MONOPOLES section corresponding to the kept atoms.
    """
    monopoles = []
    keep_names = [atom.split()[0] for atom in coords]
    start = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return monopoles
            if line.split()[0] in keep_names:
                monopoles.append(line)
        if 'MONOPOLES' in line:
            start = 1

def get_dipoles(lines, cut_coords):
    """
    Extract the dipoles section lines that should be removed.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinate lines of atoms to remove.
    
    Returns:
        dipoles (list of str): Lines from the DIPOLES section not matching the kept atoms.
    """
    dipoles = []
    cut_names = [atom.split()[0] for atom in cut_coords]
    start = 0
    for line in lines:
        if start == 1:
            if 'STOP' in line:
                return dipoles
            # Skip lines if the atom name is in the list of coordinates to cut.
            if line.split()[0] in cut_names:
                continue
            else:
                dipoles.append(line)
        if 'DIPOLES' in line:
            start = 1

def get_quadrupoles(lines, cut_coords):
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
    cut_names = [atom.split()[0] for atom in cut_coords]
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
            elif line.split()[0] in cut_names:
                k = 1
            else:
                quadrupoles.append(line)
                j = 1
        if 'QUADRUPOLES' in line:
            start = 1

def get_octupoles(lines, cut_coords):
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
    cut_names = [atom.split()[0] for atom in cut_coords]
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
            elif line.split()[0] in cut_names:
                k = 2
            else:
                octupoles.append(line)
                j = 2
        if 'OCTUPOLES' in line:
            start = 1

def get_polarpts(lines, cut_coords):
    """
    Extract the polarizable points section.
    
    For each polarizable point, compute the distance to each atom in the cut_coords.
    If the minimum distance is greater than 3, the point is retained.
    -Polarizable points do not have "atomnames" the same as the multipole parameters.
    
    Reads:
        lines (list of str): Lines from the EFP file.
        cut_coords (list of str): Coordinates of atoms to remove.
    
    Returns:
        polars (list of str): Lines from the POLARIZABLE POINTS section that pass the filter.
    """
    polars = []
    #remove_names = [atom.split()[0] for atom in cut_coords]
    start = 0
    j = 0
    for line in lines:
        mindist = 20.0  # Initialize minimum distance with a high value.
        if start == 1:
            if 'STOP' in line:
                return polars
            elif j > 0:
                polars.append(line)
                j -= 1
            # Check lines starting with "CT" (assumed polarizable point lines)
            elif line[0:2] == 'CT':
                for atom in cut_coords:
                    # Calculate the distance between the polarizable point and an atom in cut_coords.
                    current_dist = distance(line.split()[1], line.split()[2], line.split()[3],
                                            atom.split()[1], atom.split()[2], atom.split()[3])
                    if current_dist < mindist:
                        mindist = current_dist
                # If the minimum distance is greater than 3, keep this polarizable point.
                if mindist > 3:
                    polars.append(line)
                    j = 3  # Use a counter to skip the next 3 lines.
        if 'POLARIZABLE POINTS' in line:
            start = 1

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
        header.append(line.replace('$FRAGNAME', fragname))
        # Stop adding header lines once the COORDINATES section is encountered.
        if 'COORDINATES (BOHR)' in line:
            return header

def main(inp, efp):
    """
    Main processing function.
    
    1. Read the input (.inp) file and the reference EFP file.
    2. Build the header by replacing '$FRAGNAME' with the fragment name.
    3. Parse the input file to determine which atoms/polarizable points are to be removed.
    4. Extract coordinate lines (keeping and removal sets).
    5. Extract monopoles, dipoles, quadrupoles, octupoles, and polarizable points.
    6. Extract screening parameter sections.
    7. Write the output file 'cut_<efp>' containing the header and the filtered sections.
    """
    # Read input files
    with open(inp, 'r') as inp_:
        inp_lines = inp_.readlines()
    with open(efp, 'r') as efp_:
        efp_lines = efp_.readlines()

    # Construct header by replacing '$FRAGNAME' with the actual fragment name.
    header = get_header(efp_lines, inp)
    
    # Get atom indexes for complete removal and for polarizable point removal.
    rem_atoms, rem_pols = get_specials(inp_lines)
    
    # Extract coordinate lines:
    # keep_coords: coordinate lines to be kept,
    # rem_coords: coordinate lines flagged for removal.
    keep_coords, rem_coords = get_coords(efp_lines, rem_atoms, rem_pols)
    
    # Extract sections for monopoles, dipoles, quadrupoles, octupoles, and polarizable points.
    keep_monop = get_monopoles(efp_lines, keep_coords)
    keep_dip = get_dipoles(efp_lines, rem_coords)
    keep_quadrup = get_quadrupoles(efp_lines, rem_coords)
    keep_octup = get_octupoles(efp_lines, rem_coords)
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
        outfile.write(' STOP\nQUADRPOLES\n')
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
        for outline in keep_screen:
            outfile.write(outline)
        outfile.write('STOP\nSCREEN       (FROM VDWSCL=   0.700)\n')
        for outline in keep_screen2:
            outfile.write(outline)
        outfile.write('STOP\n $END\n')
    
if __name__ == "__main__":
    # Execute the main function using command-line arguments.
    # Example execution: python cut_qm.py ala_33_473.inp a0001.efp
    main(sys.argv[1], sys.argv[2])

# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 2024

@author: jackl
"""

import sys
import numpy as np
import os




# Directory containing the library files. CHANGE#######################################
base_directory = '/depot/lslipche/data/yb_boss/flexible_efp/efpdb/'

ang_cutoff = 0.20                        # Minimum RMSD allowed for a "good" match (in Angstroms)
bohr_cutoff = ang_cutoff * 1.8897259886  # Convert Angstrom cutoff to Bohr

#Add to these if there are more 
amino_acid_dict = {
    'a': 'ala', 'r': 'arg', 'n': 'asn', 'd': 'asp', 'c': 'cys',
    'q': 'gln', 'e': 'glu', 'g': 'gly', 'h': 'hip', 'i': 'ile',
    'l': 'leu', 'k': 'lys', 'm': 'met', 'f': 'phe', 'p': 'pro',
    's': 'ser', 't': 'thr', 'w': 'trp', 'y': 'tyr', 'v': 'val',
    'hp': 'hip', 'hd': 'hid', 'he': 'hie'
}

# "M" below is Mg. If atoms are present in your system but not defined here
#     you will need to add them here
atom_weights = {
    'H': 1.0, 'O': 15.9949100, 'C': 12.0000000,
    'N': 12.0030700, 'S': 32.0650000, 'M': 23.3040000
}


# Get input filename from command-line argument
inp = sys.argv[1]    
#              For testing, you might use a fixed input file:
# inp = 'g_952.inp'  # e.g., a_22_304.inp. Input file with fragment coordinates

# Read input file
with open(inp, "r") as orig_file:
    input_lines = orig_file.readlines()

# If good case found, output file name will be the same but .efp instead of .inp
# (Note: This file is only used if a match is found.)
outname = inp.replace('.inp', '.efp')




def get_RMSD(coords1,coords2,mass):
    """
    Compute the weighted root-mean-square deviation (RMSD) between two sets of coordinates.
    
    Parameters:
        coords1, coords2: 2D lists or arrays containing coordinate data.
        weights: List of weights for each coordinate point.
    
    Returns:
        Weighted RMSD.
    """
    tot=0.0
    length=len(coords1)
    dim=len(coords1[0])
    totmass=0.0
    
    for i in range(length):
        #Sum distances for each atom, weighted by atomic mass
        tot+=mass[i]*sum([(coords1[i][j]-coords2[i][j])**2.0 for j in range(dim)])
        totmass+=mass[i]
        
    # Normalized RMSD calculation.
    rmsd=np.sqrt(tot)/np.sqrt(totmass)
    return rmsd

def kabsch_algorithm(coords, coords2):
    """
    Compute the optimal rotation matrix using the Kabsch algorithm to align coords onto coords2
    
    Parameters:
        coords: Coordinates to be rotated 
        coords2: Reference coordinates 
    
    Returns:
        rotation_matrix: Optimal rotation matrix.
        com1: Center-of-mass (mean) of the first coordinates
        com2: Center-of-mass (mean) of the second coordinates
    """
    # Compute centers of mass
    com1 = np.mean(coords, axis=0)
    com2 = np.mean(coords2, axis=0)
    
    # Center the coordinates
    coords1_centered = coords - com1
    coords2_centered = coords2 - com2
    
    # Compute the covariance matrix
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    
    # Compute the SVD of the covariance matrix
    u, _, vh = np.linalg.svd(covariance_matrix)
    rotation_matrix = np.dot(u, vh)
    
    #Check for unwanted reflections/make sure rotation is "right-handed"
    if np.linalg.det(rotation_matrix) < 0:
        rotation_matrix[:, -1] *= -1
        
    return rotation_matrix, com1, com2

def apply_transform(coords, rotation_matrix, com1, com2):
    """
    Reads:
        coords: Coordinates to transform
        rotation_matrix: Rotation matrix computed from the Kabsch algorithm
        com1: Center-of-mass of the original coordinates
        com2: Center-of-mass of the target coordinates
    
    Returns:
        Transformed coordinates
    """
    # Translate coordinates to the origin
    translated_coords = coords - com1
    # Rotate the translated coordinates.
    rotated_coords = np.dot(translated_coords, rotation_matrix)
    # Translate to the target center of mass
    transformed_coords = rotated_coords + com2
    return transformed_coords


amino_acid_dict = {
    'a': 'ala', 'r': 'arg', 'n': 'asn', 'd': 'asp', 'c': 'cys',
    'q': 'gln', 'e': 'glu', 'g': 'gly', 'h': 'hip', 'i': 'ile',
    'l': 'leu', 'k': 'lys', 'm': 'met', 'f': 'phe', 'p': 'pro',
    's': 'ser', 't': 'thr', 'w': 'trp', 'y': 'tyr', 'v': 'val',
    'hp': 'hip', 'hd': 'hid', 'he': 'hie'
}

atom_weights = {
    'H': 1.0, 'O': 15.9949100, 'C': 12.0000000,
    'N': 12.0030700, 'S': 32.0650000, 'M': 23.3040000
}




# Extract the residue code from the input filename (assumes filename starts with a residue code)
try:
    res = amino_acid_dict[inp.split('_')[0]]
except KeyError:
    print('No library folder for: ' + inp)
    sys.exit()

directory = os.path.join(base_directory, res)

# Define the file extension for matching EFP files.
file_extension = '.efp'

# --------------------------
# Read coordinates from the input file
# --------------------------

# (1 Angstrom = 0.529177249 Bohr)
ANGSTROM_TO_BOHR = 0.529177249

orig_coords = []  # list of [x, y, z] coordinates in Bohr
weights = []      # list of atomic weights corresponding to each coordinate

reading_coords = False
for line in input_lines:
    if '$end' in line:
        reading_coords = False
    elif reading_coords:
        parts = line.split()
        # Convert coordinates from Angstroms to Bohr
        x = float(parts[2]) / ANGSTROM_TO_BOHR
        y = float(parts[3]) / ANGSTROM_TO_BOHR
        z = float(parts[4]) / ANGSTROM_TO_BOHR
        
        orig_coords.append([x, y, z])
        # Use the first character of the atom type to lookup the weight
        try:
            weights.append(atom_weights[parts[0][0]])
        except:
            print('Atom '+parts[0]+' is not defined. Define it and its atomic mass in `atom_weights`.')
    elif 'C1' in line:
        # Start reading coordinates after encountering 'C1'
        reading_coords = True

# --------------------------
# Loop through EFP fragment files in the specified directory
# --------------------------

min_rmsd = 20.0  # Initialize with a large value; will be updated if a better match is found.
match = None     # Will store the file path of the best match

# Loop over files with the specified extension in the directory.
for filename in os.listdir(directory):
    if not filename.endswith(file_extension):
        continue  # Skip files that do not match the extension

    full_path = os.path.join(directory, filename)

    # Read the candidate file.
    with open(full_path, "r") as candidate_file:
        candidate_lines = candidate_file.readlines()

    # Grab coordinates from the fragment file
    # The coordinates are assumed to be in the section following the header "COORDINATES (BOHR)" 
    # and ending when a line containing 'BO21' is encountered, this is the standard name for the first "bond midpoint"
    candidate_coords = []
    reading_candidate_coords = False
    for line in candidate_lines:
        if 'BO21' in line:
            reading_candidate_coords = False
        elif reading_candidate_coords:
            parts = line.split()
            x = float(parts[1])
            y = float(parts[2])
            z = float(parts[3])
            candidate_coords.append([x, y, z])
        elif 'COORDINATES (BOHR)' in line:
            reading_candidate_coords = True

    # Convert lists to numpy arrays
    candidate_coords = np.array(candidate_coords)
    fragment_coords = np.array(orig_coords)

    # Use the Kabsch algorithm to compute the optimal rotation to impost structure 1 onto structure 2
    rotation_matrix, com_candidate, com_fragment = kabsch_algorithm(candidate_coords, fragment_coords)

    # Apply the transformation
    aligned_candidate_coords = apply_transform(candidate_coords, rotation_matrix, com_candidate, com_fragment)

    # Compute the weighted RMSD between the transformed coordinates and the input
    current_rmsd = get_RMSD(aligned_candidate_coords, fragment_coords, weights)

    # Update best match if current RMSD is lower.
    if current_rmsd < min_rmsd:
        min_rmsd = current_rmsd
        match = full_path

# --------------------------
# Transform library .efp file to input coordinates if good match found.
# --------------------------
if min_rmsd < bohr_cutoff:
    # If a good match is found, run the next script with the matching file and the input file.
    os.system("python step4.Flexible_V5.py " + match + " " + inp + " -d -xr")
else:
    # If no match is found, print a message (or launch GAMESS).
    print('No match, run GAMESS for: ' + inp)
    # Uncomment the following line if you want to run GAMESS:
    # os.system('gms_slurm -p 20 -q standby -v 2023 ' + inp)

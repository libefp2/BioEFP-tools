# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 12:28:23 2024

@author: jackl

This script reads in EFP region file (.g96), a full structure file (.g96), and a topology file (.top or .itp)
to extract MM (molecular mechanics) coordinates, charges, and screening parameters.The extracted information 
is then written to an output file ("prot.efp").
"""

import sys

# Global dictionaries and lists

# 'SOL' returns the charge of oxygen in the water model (TIP3P). H charge is taken to be -1/2 * oxygen charge.
# Change these for different water models.
water_and_ions = {
    'SOL': -0.834,
    'CL': -1.0,
    'NA': 1.0
}

#All non-amino acid residue names should be in this list
known_cofactors = ['ECH', '45D', 'EQ3', 'C7Z', 'CLA', 'PQN', 'BCR', 'QLA', 'LHG', 'LMG', 'SQD', 'LMT']


def get_EFPs(efp_lines):
    """
    Grab the EFP file lines to extract unique residue IDs and names.
    
    Parameters:
        efp_lines (list of str): Lines from an EFP .g96 file.
        
    Returns:
        efp_resis (list of [resid, resname]): A list of residue information.
    """
    efp_resis = []
    start = False
    prev_res = None
    for line in efp_lines:
        if 'END' in line:
            if len(efp_resis) > 1:
                return efp_resis
        elif start:
            # When in the POSITION block, add a new residue if different from previous.
            parts = line.split()
            if parts and (parts[0] != prev_res):
                resid = parts[0]
                resname = parts[1]
                efp_resis.append([resid, resname])
                prev_res = resid
        if 'POSITION' in line:
            start = True
    return efp_resis


def get_MM_coords(mm_indexes, g96_lines):
    """
    Extract MM atom coordinates from a full configuration .g96 file.
    
    For each residue (given by mm_indexes), this function extracts coordinates 
    for the MM atoms. A conversion factor (18.897161646321) converts nm -> Bohr.
    
    Parameters:
        mm_indexes (list of [resid, resname]): Residues to process.
        g96_lines (list of str): Lines from the full configuration .g96 file.
        
    Returns:
        MMs (list of str): Formatted coordinate lines.
    """
    MMs = []
    temp_MMs = []
    prev_MM_flag = 0
    index = 0
    prev_resid = None
    conversion = 18.897161646321  # nm -> Bohr

    for line in g96_lines:
        parts = line.split()
        if len(parts) < 4 or line[0] != ' ':
            continue

        # If residue changes and a MM atom was processed, advance the index.
        if parts[0] != prev_resid and prev_MM_flag == 1:
            index += 1

        # Process lines with atom type "C" or "O"; we do not know yet if these are needed.
        if parts[2] in ('C', 'O'):
            # Build a formatted atom label: first letter of atom type + atom ID (from field 3)
            col1 = parts[2][0] + parts[3]
            # Pad the label to length 7.
            col1 = col1.ljust(7)
            # Convert coordinates with the conversion factor.
            '''
            col2 = float(parts[4]) * conversion
            col3 = float(parts[5]) * conversion
            col4 = float(parts[6]) * conversion
            temp_MMs.append(col1 + '%20.12f' % col2 + '%20.12f' % col3 + '%20.12f' % col4 + '\n')
            '''
            x, y, z = [float(parts[i]) * conversion for i in range(4, 7)]
            col2 = f"{x:.12f}".rjust(20)
            col3 = f"{y:.12f}".rjust(20)
            col4 = f"{z:.12f}".rjust(20)
            temp_MMs.append(f"{col1}{col2}{col3}{col4}\n")

        else:
            # For non-"C"/"O" lines, check if the residue matches the expected mm_indexes information.
            if parts[0] == mm_indexes[index][0] and parts[1] == mm_indexes[index][1]:
                if len(temp_MMs) > 1:
                    # Append the last two entries from the temporary list, only the most recent "C" and "O"
                    MMs.append(temp_MMs[-2])
                    MMs.append(temp_MMs[-1])
                    temp_MMs = []
                # Append current atom
                col1 = parts[2][0] + parts[3]
                col1 = col1.ljust(7)
                x, y, z = [float(parts[i]) * conversion for i in range(4, 7)]
                col2 = f"{x:.12f}".rjust(20)
                col3 = f"{y:.12f}".rjust(20)
                col4 = f"{z:.12f}".rjust(20)
                MMs.append(f"{col1}{col2}{col3}{col4}\n")
                prev_MM_flag = 1
            else:
                prev_MM_flag = 0
        prev_resid = parts[0]
    return MMs


def get_MM_charges(mm_indexes, topol_lines):
    """
    Extract MM atom charges from the topology file.
    Charges for water ("SOL") and ions are assigned from the water_and_ions dictionary.
    For each residue (from mm_indexes) the corresponding charge information is appended.
    
    Parameters:
        mm_indexes (list of [resid, resname]): Residues to process.
        topol_lines (list of str): Lines from the topology file.
        
    Returns:
        MMs (list of str): Formatted charge lines.
    """
    MMs = []
    temp_MMs = []
    prev_MM_flag = 0
    index = 0
    prev_resid = None

    for line in topol_lines:
        # When reaching the [ bonds ] section, finish processing by adding charges for remaining mm_indexes.
        if '[ bonds ]' in line:
            # Parse last processed charge to determine atom id for subsequent water/ion entries.
            atomid_str = MMs[-1].split()[0]
            # Example: if atom label is like "O12", extract "12"
            atomid = int(atomid_str[1:])  
            while index < len(mm_indexes):
                # No need to "read" water charges. TIP3P is desired model.
                if mm_indexes[index][1] == 'SOL':
                    atomid += 1
                    col1 = ('O' + str(atomid)).ljust(7)
                    col2 = water_and_ions['SOL']
                    MMs.append(col1 + '%21.10f\n' % col2)
                    # For water, assign half the oxygen charge, changed to positive
                    col2 = water_and_ions['SOL'] / (-2)
                    atomid += 1
                    col1 = ('H' + str(atomid)).ljust(7)
                    MMs.append(col1 + '%21.10f\n' % col2)
                    atomid += 1
                    col1 = ('H' + str(atomid)).ljust(7)
                    MMs.append(col1 + '%21.10f\n' % col2)
                else:
                    atomid += 1
                    col1 = (mm_indexes[index][1][0] + str(atomid)).ljust(7)
                    # Ion charge found by atom name
                    col2 = water_and_ions[mm_indexes[index][1]]
                    MMs.append(col1 + '%21.10f\n' % col2)
                index += 1
            return MMs

        # Skip non-data lines.
        if line[0] != ' ':
            continue

        parts = line.split()
        # If residue changes and previous MM was processed, move to the next index.
        if parts[2] != prev_resid and prev_MM_flag == 1:
            index += 1
            prev_MM_flag = 0

        # Process lines with atom types "C" or "O"; not known yet which C and O needed.
        if parts[1] in ('C', 'O'):
            col1 = parts[1][0] + parts[0]
            col1 = col1.ljust(7)
            col2 = float(parts[6])
            temp_MMs.append(col1 + '%21.10f\n' % col2)
        else:
            # Check if the current residue matches mm_indexes.
            if parts[2] == mm_indexes[index][0] and parts[3] == mm_indexes[index][1]:
                if len(temp_MMs) > 1:
                    MMs.append(temp_MMs[-2])
                    MMs.append(temp_MMs[-1])
                    temp_MMs = []
                col1 = parts[4][0] + parts[0]
                col1 = col1.ljust(7)
                col2 = float(parts[6])
                MMs.append(col1 + '%21.10f\n' % col2)
                prev_MM_flag = 1
            else:
                prev_MM_flag = 0
        prev_resid = parts[2]
    return MMs


def charges_from_spec_topol(topol_lines, last_atom):
    """
    Generate charge lines from a specialized topology file.(in the case that several .itp files are used)
    Each line is formatted with an atom name, the charge (from the topology), and a fixed zero multipole.
    
    Parameters:
        topol_lines (list of str): Lines from a specialized topology file.
        last_atom (int): The atom counter from which to continue numbering.
        
    Returns:
        outlines (list of str): Formatted charge lines.
    """
    i = last_atom
    col3 = '        0.0000000000'
    outlines = []
    for line in topol_lines:
        # Bonds section beginning means atom charges section is done
        if '[ bonds ]' in line:
            return outlines
        # Atom lines have first character as space, ignore any other lines
        if line[0] != ' ':
            continue
        # Atom ID in topology is not "correct." Instead of using this, take last_atom to be the 
        #    'continuation' point.
        else:
            i += 1
            col2 = f"{float(line.split()[6]):14.10f}"
            atomname = line.split()[4][0] + str(i)
            col1 = f"{atomname}".ljust(14)
            outlines.append(f"{col1}{col2}{col3}\n")
    return outlines


def get_screen(charges):
    """
    Generate screening lines for a list of charge lines.
    For each atom in charges, a screening line is produced with fixed screening parameters.
    Screening terms are not "read" or generated. They are all the same.
    
    Parameters:
        charges (list of str): List of formatted charge lines.
        
    Returns:
        screens (list of str): Formatted screening lines.
    """
    screens = []
    #every atom that has charges listed also needs screen paramters.
    for atom in charges:
        col1 = atom.split()[0]
        col1 = col1.ljust(7)
        # Screening parameters are fixed (e.g., vdW scaling and cutoff)
        col2 = '   1.0000000000  10.0000000000\n'
        screens.append(col1 + col2)
    return screens


def main(efp_g96, full_g96, topol_file):
    """
    Main routine:
      - Reads the EFP region structure file, full structure file, and topology file.
      - Extracts residue information and determines which residues belong to MM.
      - Extracts MM coordinates, charges, and screening parameters.
      - Writes the results to "prot.efp".
    """
    with open(efp_g96, 'r') as inp:
        shell_lines = inp.readlines()
    with open(full_g96, 'r') as efp:
        full_lines = efp.readlines()
    with open(topol_file, 'r') as top:
        topol_lines = top.readlines()

    # Get residue information from the EFP and full structure files.
    efp_residues = get_EFPs(shell_lines)
    all_residues = get_EFPs(full_lines)

    # Separate residues into those for MM processing and those considered as cofactors.
    mm_residues = []
    separate_topol = []
    for res in all_residues:
        # Skip residues with name "XXX"; these are link atoms
        if res[1] == 'XXX':
            continue
        # Skip residues that are in the EFP region; not needed for classical region
        elif res in efp_residues:
            continue
        # If residue name is not in known cofactors, add to MM residues.
        # known_cofactors are residues that have separate topology and will not be found 
        #     in the standard topology file. This script expeects to find an .itp file
        #     for every known_cofactor that will be used instead of the master topology.
        elif res[1] not in known_cofactors:
            mm_residues.append(res)
        else:
            separate_topol.append(res)

    # Extract MM coordinates and charges.
    MM_coords = get_MM_coords(mm_residues, full_lines)
    MM_charge = get_MM_charges(mm_residues, topol_lines)

    # For residues in separate_topol, obtain additional charges from their specific topology.
    for residue in separate_topol:
        # Use the last atom from MM_charge to set the numbering
        last_atom_str = MM_charge[-1].split()[0]
        # Remove extra spaces and extract the numeric part (assuming format like "X<number>")
        last_ID = int(last_atom_str.strip()[1:])
        # This example expects .itp contained in a folder named "amber03.ff"
        # Change this as needed!
        topol_filename = 'amber03.ff/' + residue[1] + '.itp'
        with open(topol_filename, 'r') as toplines_file:
            spec_topol_lines = toplines_file.readlines()
        temp_charges = charges_from_spec_topol(spec_topol_lines, last_ID)
        for atom_charge in temp_charges:
            MM_charge.append(atom_charge)

    MM_screen2 = get_screen(MM_charge)

    # Write the output file with coordinates, charges, and screening information.
    # Headings and sections are written explicitly here.
    with open('prot.efp', 'w') as outfile:
        outfile.write(' $PROT\n')
        outfile.write('TITLE\n')
        outfile.write('  COORDINATES (BOHR)\n')
        for outline in MM_coords:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write('MONOPOLES\n')
        for outline in MM_charge:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write('SCREEN2      (FROM VDWSCL=   0.700)\n')
        for outline in MM_screen2:
            outfile.write(outline)
        outfile.write('STOP\n')
        outfile.write(' $END')


if __name__ == "__main__":
    # Example usage: python script.py efp_pair53004.g96 confout_pair53004.g96 edit_topol.itp
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    # For testing:
    # main('efp_pair53004.g96', 'confout_pair53004.g96', 'edit_topol.itp')

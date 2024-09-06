# ------------------------------------------------------------------------
# Name: bioMAKEFP.py
# Authors: Andres S. Urbina
# Email: asurbinab@gmail.com
# Created Date: 2024-09-05
# Last Modified Date: 2024-09-05
# Version: 1.0.0
# Description: Generation of GAMESS MAKEFP input files associated to
# amino acids, ligands, and water molecules existing in the solvation
# shell of interest.
# ------------------------------------------------------------------------

"""
Change Log:
Version 1.0.0 : 2024-09-05
    - Initial creation of the script with basic functionality. (ASU)
Version 1.1.0 : YYYY-MM-DD
    - Added new feature X.
    - Fixed bug Y.
"""
import time
import configparser
import sys
import math
import subprocess
import re

# Define the version
__version__ = "1.0.0"

def main(ifile1,ifile2,ifile3):
    start_time = time.time()

    print("Running script version:", __version__)

    ligands = []

    # Open and read the ligands file
    with open('ligands', 'r') as file:
        ligands = [line.strip() for line in file]

    taas = []
    with open('taas', 'r') as file:
        taas = [line.strip() for line in file]

    config = configparser.ConfigParser()
    config.read('settings')
    ligs = config.get('Settings', 'ligands')
    sf = config.get('Settings', 'sf')

    regex = re.compile(r'\d+')
    s = regex.search(ifile2).group(0)

    with open(ifile1,'r') as f1, open(ifile2,'r') as f2, open(ifile3,'r') as f3:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        lines3 = f3.readlines()

    lastlines1 = len(lines1) - 4
    lastlines2 = len(lines2) - 4
    lenlines3 = len(lines3)

    maplist_resnum, maplist_resnum_lig, maplist_resname, maplist_xcoord, maplist_ycoord, maplist_zcoord = [], [], [], [], [], []
    maplist_resnum_ca, maplist_xcoord_ca, maplist_resnum_c, maplist_xcoord_c = [], [], [], []
    maplist_linenum, maplist_no_superfrag = [], []

    # Process lines2 in a single pass to fill various maps

    for i in range(4,lastlines2,1):
        line2 = lines2[i].split()
        if ligs == 'yes':
            if line2[0] not in maplist_resnum_lig and line2[1] in ligands:
                maplist_resnum_lig.append(line2[0])
        if line2[0] not in maplist_resnum and line2[1] not in taas and line2[1] not in ligands:
            maplist_resnum.append(line2[0])
            maplist_resname.append(line2[1])
            maplist_xcoord.append(line2[4])
            maplist_ycoord.append(line2[5])
            maplist_zcoord.append(line2[6])
        if line2[2] == 'CA' and line2[1] not in taas and line2[1] not in ligands:
            maplist_resnum_ca.append(line2[0])
            maplist_xcoord_ca.append(line2[4])
        if line2[2] == 'C' and line2[1] not in taas and line2[1] not in ligands:
            maplist_resnum_c.append(line2[0])
            maplist_xcoord_c.append(line2[4])

    if ligs == 'yes':
        lenmaplistlig = len(maplist_resnum_lig)
    lenmaplist = len(maplist_resnum)
    lenmaplist_resnum_ca = len(maplist_resnum_ca)

    # Process lines1 in a single pass to fill line numbers

    for i in range(7,lastlines1,1):
        line1 = lines1[i].split()
        j=i+1
        if line1[4] in maplist_xcoord and line1[5] in maplist_ycoord and line1[6] in maplist_zcoord:
            maplist_linenum.append(j)

#================================================
#EFP fragments: Geometries of capped molecules
#existing in the solvation shell in MAKEFP 
#input file format
#================================================

    fragment_data = {}

    #GAMESS input file header
    if ligs == 'yes':
        for k in range(0,lenmaplistlig,1):
            fragfile = 'f_%s_%s.inp' %(maplist_resnum_lig[k],s)
            header = (f' $contrl units=angs local=boys runtyp=makefp\n mult=1 '
                           f'icharg={-1 if maplist_resname[k] in {"ASP", "GLU"} else (1 if maplist_resname[k] in {"ARG", "LYS"} else 0)} '
                           f'coord=cart icut=11 $end\n $system timlim=99999 mwords=500 $end\n '
                           f'$scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06 $end\n $basis gbasis=n31 ngauss=6\n ndfunc=1 $end\n '
                           f'$local maxloc=1000 $end\n $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n '
                           f'$MAKEFP POL=.t. DISP=.f. CHTR=.f. EXREP=.f. DISP7=.f. $end\n $data\n fragment\nC1\n')
            fragment_data[fragfile] = [header]

    for k in range(0,lenmaplist,1):
        fragfile = 'f_%s_%s.inp' %(maplist_resnum[k],s)
        header = (f' $contrl units=angs local=boys runtyp=makefp\n mult=1 '
                       f'icharg={-1 if maplist_resname[k] in {"ASP", "GLU"} else (1 if maplist_resname[k] in {"ARG", "LYS"} else 0)} '
                       f'coord=cart icut=11 $end\n $system timlim=99999 mwords=500 $end\n '
                       f'$scf soscf=.f. dirscf=.t. diis=.t. CONV=1.0d-06 $end\n $basis gbasis=n31 ngauss=6\n ndfunc=1 $end\n '
                       f'$local maxloc=1000 $end\n $DAMP IFTTYP(1)=2,0 IFTFIX(1)=1,1 thrsh=500.0 $end\n '
                       f'$MAKEFP POL=.t. DISP=.f. CHTR=.f. EXREP=.f. DISP7=.f. $end\n $data\n fragment\nC1\n')
        fragment_data[fragfile] = [header]

    #Ligands, do not need capping
    if ligs == 'yes':
        for i in range(4,lastlines2,1):
            line2 = lines2[i].split()
            if line2[0] in maplist_resnum_lig:
                x = float(line2[4])*10
                y = float(line2[5])*10
                z = float(line2[6])*10
                fragfile = 'f_%s_%s.inp' %(line2[0],s)
                if line2[2][0] == 'C':
                    if line2[2][1] == 'L':
                        anum = 17.0
                        fragment_data[fragfile].append(f' Cl     {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                    else:
                        anum = 6.0
                        fragment_data[fragfile].append(f' C      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'H':
                    anum = 1.0
                    fragment_data[fragfile].append(f' H      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'O':
                    anum = 8.0
                    fragment_data[fragfile].append(f' O      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'N':
                    anum = 7.0
                    fragment_data[fragfile].append(f' N      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'S':
                    anum = 16.0
                    fragment_data[fragfile].append(f' S      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'F':
                    anum = 9.0
                    fragment_data[fragfile].append(f' F      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'B':
                    if line2[2][1] == 'R':
                        anum = 35.0
                        fragment_data[fragfile].append(f' Br    {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                elif line2[2][0] == 'M':
                    if line2[2][1] == 'G':
                        anum = 12.0
                        fragment_data[fragfile].append(f' Mg    {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                maplist_no_superfrag.append(line2[4])

    #Amino acids that need capping and water molecules
    for i in range(4,lastlines2,1):
        line2 = lines2[i].split()
        if line2[0] in maplist_resnum and line2[2] not in {'C', 'O'} or line2[1] in {'SOL', 'TP3', 'QSL'}:
            x = float(line2[4])*10
            y = float(line2[5])*10
            z = float(line2[6])*10
            fragfile = 'f_%s_%s.inp' %(line2[0],s)
            if line2[2][0] == 'C':
                if line2[2][1] == 'L':
                    anum = 17.0
                    fragment_data[fragfile].append(f' Cl     {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
                else:
                    anum = 6.0
                    fragment_data[fragfile].append(f' C      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'H':
                anum = 1.0
                fragment_data[fragfile].append(f' H      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'O':
                anum = 8.0
                fragment_data[fragfile].append(f' O      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'N':
                anum = 7.0
                fragment_data[fragfile].append(f' N      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'S':
                anum = 16.0
                fragment_data[fragfile].append(f' S      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'F':
                anum = 9.0
                fragment_data[fragfile].append(f' F      {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'B':
                if line2[2][1] == 'R':
                    anum = 35.0
                    fragment_data[fragfile].append(f' Br    {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            elif line2[2][0] == 'M':
                if line2[2][1] == 'G':
                    anum = 12.0
                    fragment_data[fragfile].append(f' Mg    {anum}   {x:.8f}   {y:.8f}   {z:.8f}\n')
            maplist_no_superfrag.append(line2[4])

    #Capping from above

    for k in range(0,lenmaplist_resnum_ca,1):
        fragfile = 'f_%s_%s.inp' %(maplist_resnum_ca[k],s)
        for j in range(maplist_linenum[k],7,-1):
            line1 = lines1[j].split()
            if line1[2] == 'O' and line1[1] not in {'SOL', 'TP3', 'QSL'} and line1[1] not in ligands:
                x1 = float(line1[4])*10
                y1 = float(line1[5])*10
                z1 = float(line1[6])*10
                maplist_no_superfrag.append(line1[4])
                break
        for j in range(maplist_linenum[k],7,-1):
            line1 = lines1[j].split()
            if line1[2] == 'C':
                x2 = float(line1[4])*10
                y2 = float(line1[5])*10
                z2 = float(line1[6])*10
                maplist_no_superfrag.append(line1[4])
                break
        for j in range(maplist_linenum[k],7,-1):
            line1 = lines1[j].split()
            if line1[2] == 'CA':
                x3 = float(line1[4])*10
                y3 = float(line1[5])*10
                z3 = float(line1[6])*10
                break

        actual_length = math.sqrt(((x2-x3)**2)+((y2-y3)**2)+((z2-z3)**2))
        desired_length = 1.07886
        x3new = ((x3-x2)*desired_length/actual_length)+x2
        y3new = ((y3-y2)*desired_length/actual_length)+y2
        z3new = ((z3-z2)*desired_length/actual_length)+z2
        fragment_data[fragfile].append(f' C      6.0   {x2:.8f}   {y2:.8f}   {z2:.8f}\n')
        fragment_data[fragfile].append(f' O      8.0   {x1:.8f}   {y1:.8f}   {z1:.8f}\n')
        fragment_data[fragfile].append(f' H000   1.0   {x3new:.8f}   {y3new:.8f}   {z3new:.8f}\n')

    #Capping from below

    for k in range(0,lenmaplist_resnum_ca,1):
        fragfile = 'f_%s_%s.inp' %(maplist_resnum_ca[k],s)
        for i in range(4,lastlines2,1):
            line2 = lines2[i].split()
            if line2[4] == maplist_xcoord_ca[k] and line2[0] == maplist_resnum[k]:
                x1 = float(line2[4])*10
                y1 = float(line2[5])*10
                z1 = float(line2[6])*10
                break
        for i in range(4,lastlines2,1):
            line2 = lines2[i].split()
            if line2[4] == maplist_xcoord_c[k] and line2[0] == maplist_resnum[k]:
                x2 = float(line2[4])*10
                y2 = float(line2[5])*10
                z2 = float(line2[6])*10
                break
        actual_length = math.sqrt(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))
        desired_length = 1.07886
        x1new = ((x2-x1)*desired_length/actual_length)+x1
        y1new = ((y2-y1)*desired_length/actual_length)+y1
        z1new = ((z2-z1)*desired_length/actual_length)+z1
        fragment_data[fragfile].append(f' H000   1.0   {x1new:.8f}   {y1new:.8f}   {z1new:.8f}\n')

    subprocess.Popen(['mkdir','efp_%s' %s],stdin=subprocess.PIPE).wait()

    # Write all GAMESS input files at once
    for fragfile, content in fragment_data.items():
        with open(fragfile, 'w') as f:
            f.writelines(content)
            f.write(' $end')
        subprocess.Popen(['mv', fragfile, f'efp_{s}'], stdin=subprocess.PIPE).wait()

#=================================================
# Superfragment (MM point charges in the atoms
# that are not in the solvation shell)
#=================================================

    if sf == 'yes':

        superfragfile = 'sf_%s.efp' %s
        sfganum = -1
        coordinates = []
        monopoles = []
        screen2 = []
  
        maplist_anum_sf = []
        maplist_aname_sf = []
        maplist_aname2_sf = []
        maplist_charges = []

        # Prepare coordinates (in Bohr), monopoles, and screen 2 parameters of the superfragment

        for i in range(1, lenlines3,1):
            line3 = lines3[i].split()
            if line3[5] not in maplist_anum_sf:
                maplist_anum_sf.append(line3[5])
                maplist_charges.append(line3[6])

        for i in range(7, lastlines1):
            line1 = lines1[i].split()
    
            if line1[4] not in maplist_no_superfrag and line1[2] not in {'LA'}:
                x1 = float(line1[4]) * 10 * 1.8897259886
                y1 = float(line1[5]) * 10 * 1.8897259886
                z1 = float(line1[6]) * 10 * 1.8897259886
                sfganum += 1
                coordinates.append(f' O{sfganum}   {x1:.8f}   {y1:.8f}   {z1:.8f}   0.00000001   0.000000005\n')
                screen2.append(f' O{sfganum}          1.00000000      10.0000000\n')
   
                if line1[3] in maplist_anum_sf:
                    index = maplist_anum_sf.index(line1[3])
                    monopoles.append(f' O{sfganum}   {float(maplist_charges[index]):.8f}  {0.000:.8f}\n')

                elif line1[1] in {'SOL', 'TP3', 'QSL'} and line1[2][0] in {'O'}:
                    monopoles.append(f' O{sfganum}   {-0.834:.8f}  {0.000:.8f}\n')

                elif line1[1] in {'SOL', 'TP3', 'QSL'} and line1[2][0] in {'H'}:
                    monopoles.append(f' O{sfganum}   {0.417:.8f}  {0.000:.8f}\n')

                elif line1[1] in {'CL', 'Cl-'}:
                    monopoles.append(f' O{sfganum}   {-1.000:.8f}  {0.000:.8f}\n')

        # Write superfragfile
        with open(superfragfile, 'a') as f:
            f.write('          RUNTYP=MAKEFP EFFECTIVE FRAGMENT POTENTIAL DATA FOLLOWS...\n          FRAGNAMEEFP\n $FRAGNAME\nEFP DATA FOR FRAGNAME\n COORDINATES (BOHR)\n')
            f.writelines(coordinates)
            f.write(' STOP\n MONOPOLES\n')
            f.writelines(monopoles)
            f.write(' STOP\n SCREEN2      (FROM VDWSCL=   0.700)\n')
            f.writelines(screen2)
            f.write(' STOP\n $END')
    
        subprocess.Popen(['mv', superfragfile, 'efp_%s' % s], stdin=subprocess.PIPE).wait()

    end_time = time.time()  # End timing
    print(f"Execution time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])

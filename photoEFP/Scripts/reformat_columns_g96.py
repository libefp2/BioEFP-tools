# Sample execution: python format.py bchl361-79002.g96
# Output will be: formed_bchl361-79002.g96

#Occasionally, gmx traj commands will yield .g96 files that are incorrectly read by VMD to create 'sheets' 
#of solvent water molecules due to a column spacing inconsistency.
#This script reads a .g96 file and outputs a new file with columns forced to be uniform throughout. 

import sys

#file='bchl361-79002.g96'
file=sys.argv[1]

with open(file, "r") as orig:
    orgs = orig.readlines()

start=0
i=0
temp=[]
for line in orgs:
    if(len(line.split())>4):
        col0=line.split()[0]
        while(len(col0)<5):
              col0=' '+col0
        col0=col0+' '
        col1=line.split()[1]
        while(len(col1)<6):
              col1=col1+' '
        col2=line.split()[2]
        while(len(col2)<6):
              col2=col2+' '
        col3=line.split()[3]
        while(len(col3)<6):
              col3=' '+col3
        col4=line.split()[4]
        while(len(col4)<15):
              col4=' '+col4
        col5=line.split()[5]
        while(len(col5)<15):
              col5=' '+col5
        col6=line.split()[6]
        while(len(col6)<15):
              col6=' '+col6
        temp.append(col0+col1+col2+col3+col4+col5+col6+'\n')
    else:
        temp.append(line)
        

with open('formed_'+file, 'w') as out:
    for line in temp:
        out.write(line)
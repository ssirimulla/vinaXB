#! /usr/bin/env python

# pass two conformations of the same ligand in pdbqt format
# ex: python rmsd.py native.pdbqt conformation.pdbqt
# atoms must be in the same order to get accurate results
# result is in Angstroms

import math
import fileinput

dict_atom = { }
rmsd_pieces = []
first = True
readFirst = False
n1 = 0
n2 = 0

def rmsd_piece(x1,x2,y1,y2,z1,z2):
    return (x2-x1)**2 +(y2-y1)**2 +(z2-z1)**2

for line in fileinput.input():
    if fileinput.isfirstline() and first and not readFirst:
        filename = fileinput.filename()         #filename of the first file
        readFirst = True
    elif fileinput.filename() != filename and readFirst:
        first = False
        filename = fileinput.filename()         #filename of the second file (what is printed at termination)
    record_type=line[0:6]
    if first and record_type == 'HETATM':
        newLine = ' '.join(line.split()).split(' ')
        dict_atom[n1] = {'number':newLine[1],'name':newLine[2],'x':line[30:38],'y':line[38:46],'z':line[46:54]}
        n1 += 1
    elif record_type == 'HETATM':
        newLine = ' '.join(line.split()).split(' ')
        #check if the order of the atoms are consistent
        if dict_atom[n2]['number'] == newLine[1] and dict_atom[n2]['name'] == newLine[2]:
            rmsd_pieces.append(rmsd_piece(float(dict_atom[n2]['x']),float(line[30:38]),float(dict_atom[n2]['y']),float(line[38:46]),float(dict_atom[n2]['z']),float(line[46:54])))
        else:
            print "ERROR: order of atoms is inconsistent"
        n2 += 1

if n1 == n2:
    print filename,"RMSD:",math.sqrt(sum(rmsd_pieces)/n1), "A"
else:
    print "ERROR: number of atoms is inconsistent"
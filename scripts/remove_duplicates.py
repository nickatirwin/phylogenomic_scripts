# remove_duplicates.py
#
# Remove duplicated sequences in a fasta file (i.e., sequences
# with the same headers).
#
# Outputs *.unique.fasta
#
# Usage: python remove_duplicates.py [fasta]
# Example: python remove_duplicates.py mysequences.fasta
#

import sys
import subprocess

n = 0

file = open(sys.argv[1], 'r')
lines = file.readlines()
file.close()
output = open(sys.argv[1] + '.sl', 'w')

for i in lines:
	if i.startswith('>') == True:
		if n == 0:
			output.write(i)
			n += 1
		else:
			output.write('\n' + i)
			n += 1
	else:
		output.write(i.strip('\n'))
		n += 1
output.write('\n')
output.close()

out = open(sys.argv[1].split('.fa')[0]+'.unique.fasta','w')
infile = open(sys.argv[1]+'.sl','r').readlines()

n = 0
records = []
while n < len(infile):
    if infile[n] in records:
        print(infile[n].strip() + ' deleted')
        n = n + 2
    else:
        out.write(infile[n])
        out.write(infile[n+1])
        records.append(infile[n])
        n = n + 2

out.close()

subprocess.call('rm ' + sys.argv[1] + '.sl', shell = True)   

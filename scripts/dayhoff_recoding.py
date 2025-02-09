# dayhoff_recoding.py
#
# Recode an amino acid alignment file with Dayhoff groups
#
# usage: python dayhoff_recoding.py [protein alignment file] [# of categories: 4 or 6]

import sys

dayhoff6 = {
'C':'0',
'G':'1',
'S':'1',
'T':'1',
'A':'1',
'P':'1',
'D':'2',
'E':'2',
'N':'2',
'Q':'2',
'R':'3',
'H':'3',
'K':'3',
'L':'4',
'V':'4',
'M':'4',
'I':'4',
'Y':'5',
'F':'5',
'W':'5',
'X':'-',
}

dayhoff4 = {
'C':'C',
'G':'A',
'S':'A',
'T':'A',
'A':'A',
'P':'A',
'D':'G',
'E':'G',
'N':'A',
'Q':'G',
'R':'G',
'H':'C',
'K':'G',
'L':'T',
'V':'T',
'M':'T',
'I':'T',
'Y':'C',
'F':'T',
'W':'C',
'X':'-',
}


if sys.argv[2] == '6':
    dayhoff = dayhoff6
else:
    dayhoff = dayhoff4

alignment = open(sys.argv[1],'r').read().split('>')
out = open(sys.argv[1]+'.dayhoff'+sys.argv[2],'w')

for seq in alignment[1:]:
    out.write('>'+seq.split('\n')[0].strip()+'\n')
    for c in seq.split('\n',1)[1]:
        if c.upper() in dayhoff:
            out.write(dayhoff[c])
        else:
            out.write(c)
out.close()



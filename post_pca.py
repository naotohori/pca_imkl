#!/usr/bin/env python 

import sys

if len(sys.argv) != 4 :
    print ''
    print ' Usage : SCRIPT [input pca file] [# of eigen values] [prefix of output files]'
    print ''
    sys.exit(2)

f_in = open(sys.argv[1], 'r')
n_vec = int(sys.argv[2].strip())
prefix = sys.argv[3]

f_value = open(prefix+'.value', 'w')
fs_vec = []
for i in xrange(n_vec) :
    fs_vec.append(open(prefix+'_%i.vec'%(i+1,) , 'w'))

i = 0
f_out = f_value
for line in f_in :
    if line.find('#eigenvector') != -1 :
        i += 1
        if i > n_vec :
            break
        f_out = fs_vec[i-1] 
    f_out.write(line)

f_in.close()
f_value.close()
for f in fs_vec :
    f.close()



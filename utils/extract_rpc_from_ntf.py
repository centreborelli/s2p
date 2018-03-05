#!/usr/bin/python
# tools for extracting RPC from NTIF/NTF file
# Copyright (C) 2016, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>

import sys, subprocess
if len(sys.argv)<2:
      print(" extracts RPC from NTIF/NTF file"%sys.argv[0])
      print("    usage: %s file.NTF > file.RPC"%sys.argv[0])
      exit()
x=subprocess.Popen(["gdalinfo", sys.argv[1]], stdout=subprocess.PIPE).communicate()[0]
x=x.splitlines()
#x=sys.stdin.readlines()
for l in x:
      if ('SAMP_' not in l) and ('LINE_' not in l) and ('HEIGHT_' not in l) and ('LAT_' not in l) and ('LONG_' not in l) and ('MAX_' not in l) and ('MIN_' not in l):
              continue
      y = l.strip().replace('=',': ')
      if 'COEFF' in y:
              z = y.split(' ')
              t=1
              for j in z[1:]:
                      print('%s_%d: %s'%(z[0][:-1],t,j))
                      t+=1
      else:
              print(y)

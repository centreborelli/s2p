#!/usr/bin/env python
import piio

d = piio.read('testimg.tif')
print d.shape
print d[:,:,0] 
piio.write('testimg2.png',d)


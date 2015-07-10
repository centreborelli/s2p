# Primitive python wrapper for iio
# Copyright 2013, Gabriele Facciolo <facciolo@cmla.ens-cachan.fr>


import os,sys,ctypes

#if sys.platform.startswith('win'):
#   lib_ext = '.dll'
#elif sys.platform == 'darwin':
#   lib_ext = '.dylib'
#else:
#   lib_ext = '.so'
#gcc -std=c99 iio.c -shared -o iio.dylib -lpng -ltiff -ljpeg


lib_ext = '.so'
here  = os.path.dirname(__file__)
libiiofile= os.path.join(here, 'libiio'+lib_ext)

### HACK TO BUILD libiio ON THE FLY
try:
   os.stat(libiiofile)
except OSError:
   print 'BUILDING PIIO...'
   os.system('cd %s; python setup.py build'%here)

libiio   = ctypes.CDLL(libiiofile)
del libiiofile, here, lib_ext




def read(filename):
   '''
   IIO: numpyarray = read(filename)
   '''
   from numpy import array, zeros, ctypeslib
   from ctypes import c_int, c_float, c_void_p, POINTER, cast, byref

   iioread = libiio.iio_read_image_float_vec
   
   w=c_int()
   h=c_int()
   nch=c_int()
   
   iioread.restype = c_void_p  # it's like this
   tptr = iioread(str(filename),byref(w),byref(h),byref(nch))
   c_float_p = POINTER(c_float)       # define a new type of pointer
   ptr = cast(tptr, c_float_p)
   #print w,h,nch
   
   #nasty read data into array using buffer copy
   #http://stackoverflow.com/questions/4355524/getting-data-from-ctypes-array-into-numpy
   #http://docs.scipy.org/doc/numpy/reference/generated/numpy.frombuffer.html
   
   # this numpy array uses the memory provided by the c library, which will be freed
   data_tmp = ctypeslib.as_array( ptr, (h.value,w.value,nch.value) )
   # so we copy it to the definitive array before the free
   data = data_tmp.copy()
   
   # free the memory
   iiofreemem = libiio.freemem
   iiofreemem(ptr)
   return data


def write(filename,data):
   '''
   IIO: write(filename,numpyarray)
   '''
   from ctypes import c_char_p, c_int, c_float
   from numpy.ctypeslib import ndpointer

   iiosave = libiio.iio_save_image_float_vec

   h  =data.shape[0]
   w  =data.shape[1]
   nch=1
   if (len(data.shape)>2):
      nch=data.shape[2]

   iiosave.restype = None
   iiosave.argtypes = [c_char_p, ndpointer(c_float),c_int,c_int,c_int]
   iiosave(str(filename), data.astype('float32'), w, h, nch)


#d = piio.read('kk.tif')
#print d.shape
#print d[:,:,0] 
#piio.write('kk2.tif',d)

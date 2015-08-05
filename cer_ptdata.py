import ctypes
import numpy as np

lib = '/f/CER/cat/linux64/lib/cer_ptdata_linux.so'
a = ctypes.cdll.LoadLibrary(lib)
#method = 'idl_cer_ichord_'
method = a.idl_cer_ichord_
schord = 'T03'
byte_array = [ord(i) for i in schord.upper()]
int_array_tmp = ctypes.c_byte * len(schord)
int_array_c = int_array_tmp(*byte_array)
ilen_c = ctypes.c_int32(len(schord))
shot_c = ctypes.c_int32(160409)
ichord_c = ctypes.c_int32()
lib = '/f/CER/cat/linux64/lib/libcerdata.so'
ctypes.c_char_p('T03')

##b = ctypes.cdll.LoadLibrary(lib)
#ctypes.c
#method = b.cer_config_cvt_string_to_ichord_

#method(ctypes.byref(shot_c), ctypes.byref(ilen_c), ctypes.byref(int_array_c), ctypes.byref(ichord_c))

#byte_array = byte(strupcase(schord))
#ilen = n_elements(byte_array)

#dum = call_external (image, executable, long(ishot), long(ilen),    $
#                      byte_array, ichord ) 

# ;     argv(1)  is DIII-D shot number
# ;     argv(2)  is number of characters in user's string
# ;     argv(3)  is array of ascii positions (integers) for characters in
# ;                 user's string
# ; Output:
# ;     argv(4) is ichord number - 0 is returned if there is
# ;                 an error in conversion.

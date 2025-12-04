import sys
from ctypes import c_double, CDLL

# Function to call C function to compute the angular distance between 2 objects
def calculate_angular_distance_caller(ra1, dec1, ra2, dec2):
    return angdist_lib.calculateAngularDistance(
        float(ra1),
        float(dec1),
        float(ra2),
        float(dec2)
    )

lib_path = 'Practices/Practice_4/calculate_angular_distance.so'

try:
    angdist_lib = CDLL(lib_path)
except OSError:
    print('Shared library not found on platform %s' % sys.platform)
    raise

# Specify argument and return types for the C function
angdist_lib.calculateAngularDistance.argtypes = [c_double, c_double, c_double, c_double]
angdist_lib.calculateAngularDistance.restype = c_double
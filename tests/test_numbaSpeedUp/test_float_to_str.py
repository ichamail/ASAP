from numba import njit
import src.numbaSpeedUp

@njit
def print_numba_float(x):
    s = "This is the number: " + str(x)
    print(s)

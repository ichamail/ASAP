import numpy as np
from numba import types, njit
from numba.cpython.unicode import _empty_string, _set_code_point, PY_UNICODE_1BYTE_KIND
from numba.extending import overload_method

DIGITS_START = 48
DIGITS_END = 58
DASH = 45
DOT = 46
PLUS = 43
E_CHAR = 101

@njit(cache=True)
def get_n_digits(x):
    l1,l2 = 0,-1
    _x = x
    while _x > 0:
        _x = _x // 10
        l1 += 1

    _x = x % 10     
    while _x > 1e-10:
        _x = (_x * 10) % 10
        l2 += 1
        if(l2 >= 16): break
    return l1, l2

@njit
def float_to_str(x):
    if(x == np.inf):
        return 'inf'
    elif(x == -np.inf):
        return '-inf'

    isneg = int(x < 0.0)
    x = np.abs(x)
    
    if(x != 0.0):
        # There is probably a more efficient way to do this
        e = np.floor(np.log10(x))
        if(10**e - x > 0): e -= 1
    else:
        e = 0
    
    is_exp, is_neg_exp = e >= 16, e <= -16

    exp_chars = 0
    if (is_exp or is_neg_exp):
        exp_chars = 4
        if(e >= 100 or e <= -100): exp_chars = 5


    if(is_exp):
        offset_x = np.around(x * (10.0**-(e)),15)
        l1, l2 = get_n_digits(offset_x)
    elif(is_neg_exp):
        offset_x = np.around(x * (10**-(e)),15)
        l1, l2 = get_n_digits(offset_x)
    else:
        offset_x = x
        l1,l2 = get_n_digits(x)
        l2 = max(1,l2) # Will have at least .0 
    
    use_dec = l2 > 0

    # print("<<", e, offset_x, l2)

    l = l1+l2+use_dec
    length = l+isneg+exp_chars
    s = _empty_string(PY_UNICODE_1BYTE_KIND,length)
    if(isneg): _set_code_point(s,0,DASH)

    _x = offset_x
    for i in range(l1):
        digit = int(_x % 10)
        _set_code_point(s,(isneg+l1)-i-1,digit + DIGITS_START)
        _x = _x // 10

    if(use_dec):
        _set_code_point(s,l1+isneg,DOT)

    _x = offset_x % 10
    for i in range(l2):
        _x = (_x * 10) % 10
        digit = int(_x)  
        
        _set_code_point(s,(isneg+l1)+i+use_dec,digit + DIGITS_START)

    if(is_exp or is_neg_exp):
        i = (isneg+l1+use_dec+l2)
        _set_code_point(s,i,E_CHAR)        
        if(is_exp):
            _set_code_point(s,i+1,PLUS)     
        if(is_neg_exp):
            _set_code_point(s,i+1,DASH)           

        i = length-1
        exp = np.abs(e)
        while(exp > 0):
            digit = exp % 10
            _set_code_point(s,i,digit + DIGITS_START)
            exp = exp // 10
            i -= 1

    return s

# Monkey Patch type conversions so that float(str) and str(float)
# work like in normal python

@overload_method(types.Float, '__str__')
def overload_float_to_str(x):
    def impl(x):
        return float_to_str(x)
    return impl


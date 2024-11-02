import numpy as np
from scipy.stats import beta
from numpy.linalg import inv
from ..myMath import Vector


def LeastSquares(A, b):
    transposed = A.T
    inverted = inv(transposed @ A)
    x = (inverted @ transposed) @ b
    return x

def light_vector(magnitude, alpha, beta):
    
    """
    alpha: angle [degs] of light vector with x-axis in x-z plane
    beta: angle [degs] of light vector with x-axis in x-y plane
    """
    alpha = np.deg2rad(alpha)
    beta = np.deg2rad(beta)
    x_component = magnitude * np.cos(alpha) * np.cos(beta)
    y_component = magnitude * np.cos(alpha) * np.sin(beta)
    z_component = - magnitude * np.sin(alpha)
    light_vector = Vector(x_component, y_component, z_component)
    
    return light_vector

def cubic_function(x):
    
    """
    cubic function used for airfoil interpolations
    
    from
    May-Fun Liou, Hyoungjin Kim, ByungJoon Lee, and Meng-Sing Liou
    "Aerodynamic Design of the Hybrid Wing Body Propulsion-Airframe Integration"
    """
    return x**2 * (3 - 2*x)

def interpolation(rootValue, tipValue, spanPercentage, type:str="linear"):
    
    """
    linear interpolation:
    x = {(y_tip - y)*x_root + (y - y_root)*x_tip}/(y_tip - y_root)
    x = {y_tip*x_root + y*(x_tip - x_root)}/y_tip
    x = x_root + y/y_tip * (x_tip - x_root)
    """
    
    if type == "linear":
        section_value = rootValue + (tipValue - rootValue) * spanPercentage
    elif type == "cubic":
        section_value = (
            rootValue + (tipValue 
                        - rootValue) * cubic_function(spanPercentage)
            )
    
    return section_value

def DenserAtBoundaries(start, end, num_points, alpha):
    '''
    Beta distribution
    
    Cumulative distribution function of beta distribution
    
    alpha exists in (-oo, +oo)
    when alpha is 1 => evenly spaced
    when alpha < 1 denser at boundaries
    when alpha > 1 denser at midpoint
    '''
    x = np.linspace(0, 1, num_points)
    a = b = 2-alpha
    return start + beta.cdf(x, a, b) * (end-start)

def DenserAtWingTips(start, end, num_points, factor=5):

    # produces same results as flow5 inv_exp spacing when factor=5

    x = np.linspace(0, 1, num_points)
    if start < 0  and end <=0:
        a = 1 + 0.1*factor
        b = 1
    elif start>= 0 and end>0 :
        a = 1
        b = 1 + 0.1*factor

    else:
        raise Exception("start, end >=0 or start, end <=0")


    return start + beta.cdf(x, a, b) * (end-start)

def DenserAtWingRoot(start, end, num_points, factor=5):
    x = np.linspace(0, 1, num_points)
    if start < 0 :
        a = 1
        b = 1 + 0.1*factor
    else:
        a = 1 + 0.1*factor
        b = 1

    return start + beta.cdf(x, a, b) * (end-start)  

def logspace(start, end, num_points):
    x = np.linspace(1, np.e, num_points)
    x = np.log(x)
    if  0 <= start < end:
        x = x*(end-start)
        x = x + start
    elif start < end <= 0:
        x = x * (start-end)
        x = x + end
        x = np.flip(x)
    return x

def cosspace(start, end, num_points):
    mean = (start+end)/2
    amp = (end-start)/2
    return mean + amp * np.cos(np.linspace(np.pi, 0, num_points))

    

if __name__ == "__main__":
    A = np.array([[0, 1],
                [1, 1],
                [2, 1],
                [3, 1]])
    b = np.array([-1, 0.2, 0.9, 2.1])

    x = np.linalg.lstsq(A, b, rcond=None)[0]

    print(x)
    
    b = np.array([[-1],
                  [0.2],
                  [0.9],
                  [2.1]])
    
    x = LeastSquares(A, b)
    print(x)
    

from src.myMath import Vector
import numpy as np

def test_Vector():
    
    e_x = Vector(1, 0, 0)
    e_y = Vector(0, 1, 0)
    e_z = Vector(0, 0, 1)
    
    print("ex = ", e_x)
    print("ey = ", e_y)
    print("ez = ", e_z)
    
    print("2 * e_x = ", 2*e_x)
    print("ex + ey + ez = ", e_x + e_y + e_z)
    print("ex dot ey = ", e_x.dot(e_y))
    print("ex dot ex = ", e_x.dot(e_x))
    print("ex cross ey = ", e_x.cross(e_y))
    print("ey cross ez = ", e_y.cross(e_z))
    print("ez cross ex = ", e_z.cross(e_x))
    
    pass

import numpy as np
from numba.experimental import jitclass
from numba import float64

spec = [('x', float64), ('y', float64), ('z', float64)]

@jitclass(spec)
class Vector:
    
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def __str__(self):
        
        """Human-readable string representation of the vector."""
        # return '({:g})ex + ({:g})ey + ({:g})ez'.format(self.x, self.y, self.z)
        
        # using numba and overloaded str() from ulilities
        string = "(" + str(round(self.x, 4)) + ")ex + ("
        string += str(round(self.y, 4)) + ")ey + ("
        string += str(round(self.z, 4)) + ")ez"
        
        return string
    
    def norm(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)
          
    def changeBasis(self, Matrix:np.ndarray[float, np.dtype[float]]):
        
        tranformedVector = Vector(
          Matrix[0][0] * self.x + Matrix[0][1] * self.y + Matrix[0][2] * self.z,
          Matrix[1][0] * self.x + Matrix[1][1] * self.y + Matrix[1][2] * self.z,
          Matrix[2][0] * self.x + Matrix[2][1] * self.y + Matrix[2][2] * self.z
        )
        
        return tranformedVector
    
    # @classmethod  # class methods are not yet supported in numba
    # def addition(cls, vec1, vec2):
    #     return Vector(vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z)
    
    # @classmethod  # class methods are not yet supported in numba
    # def dot_product(cls, vec1, vec2):
    #     return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z
    
    # @classmethod  # class methods are not yet supported in numba
    # def cross_product(cls, vec1, vec2):
        
    #     crossProductVector = Vector(
    #         vec1.y*vec2.z - vec1.z*vec2.y,
    #         - (vec1.x*vec2.z - vec1.z*vec2.x),
    #         vec1.x*vec2.y - vec1.y*vec2.x  
    #     )
        
    #     return crossProductVector
    
    def __add__(self, vector):
        
        return Vector(self.x + vector.x, self.y + vector.y, self.z + vector.z)
    
    def __iadd__(self, vector):
        self.x += vector.x
        self.y += vector.y
        self.z += vector.z
        return self
    
    def __sub__(self, vector):
        
        return Vector(self.x - vector.x, self.y - vector.y, self.z - vector.z)
    
    def __isub__(self, vector):
        self.x -= vector.x
        self.y -= vector.y
        self.z -= vector.z
        return self
    
    def __mul__(self, scalar):
        
        return Vector(scalar * self.x, scalar * self.y, scalar * self.z)
    
    def __rmul__(self, scalar):
        
        return Vector(scalar * self.x, scalar * self.y, scalar * self.z)
    
    def __rmatmul__(self, Matrix:np.ndarray[float, np.dtype[float]]):
        
        tranformedVector = Vector(
          Matrix[0][0] * self.x + Matrix[0][1] * self.y + Matrix[0][2] * self.z,
          Matrix[1][0] * self.x + Matrix[1][1] * self.y + Matrix[1][2] * self.z,
          Matrix[2][0] * self.x + Matrix[2][1] * self.y + Matrix[2][2] * self.z
        )
        
        return tranformedVector
        
    def __truediv__(self, scalar):
        
        return Vector(self.x/scalar, self.y/scalar, self.z/scalar)
    
    def __pos__(self):
        return self
       
    def __neg__(self):
        
        return Vector(-self.x, -self.y, -self.z)
    
    def __eq__(self, vector) -> bool:
        
        return (
            self.x == vector.x and self.y == vector.y and self.z == vector.z
        )
    
    def __ne__(self, vector) -> bool:
        return not self == vector
     
    def dot(self, vector):
        return self.x * vector.x + self.y * vector.y + self.z * vector.z
    
    def cross(self, vector):
        
        crossProductVector = Vector(
            self.y * vector.z - self.z * vector.y,
            - (self.x * vector.z - self.z * vector.x),
            self.x * vector.y - self.y * vector.x
        )
        
        return crossProductVector

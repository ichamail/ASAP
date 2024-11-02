import numpy as np
from .icosphere import icosphere
from .cube import Cube


class Sphere:
    
    def __init__(self, center:tuple[float], radius:float) -> None:
        self.center = center
        self.radius = radius
        pass
    
    def meshUVSurface(self, numOfMeridians:int, numOfParallels:int):
        
        def vertexIndex(meridian, parallel):
            return (
                (parallel*numOfMeridians + meridian%numOfMeridians)
                + 1
            )
        
        l = 0
        def addFace(*vertexIds):
            nonlocal l
            for i in range(len(vertexIds)):
                face[l][i] = vertexIds[i]
            l  = l + 1
        
        vertex = np.zeros(
            shape=(numOfParallels * numOfMeridians + 2, 3),
            dtype=float
        )
        face = -np.ones(
            shape=(numOfMeridians*(numOfParallels+1), 4),
            dtype=int
        )
        
        vertex[0] = 0, 0, self.radius
        vertex[-1] = 0, 0, -self.radius
        
        k=0
        for i in range(numOfParallels):
            
            thetaParallel = np.pi * (i+1)/(numOfParallels+1)
            sinTheta = np.sin(thetaParallel)
            cosTheta = np.cos(thetaParallel)
            
            for j in range(numOfMeridians):
                
                k = k+1
                
                phiMeridian = 2 * np.pi * j/numOfMeridians
                sinPhi = np.sin(phiMeridian)
                cosPhi = np.cos(phiMeridian)
                
                vertex[k][0] = self.center[0] + self.radius * sinTheta * cosPhi
                vertex[k][1] = self.center[1] + self.radius * sinTheta * sinPhi
                vertex[k][2] = self.center[2] + self.radius * cosTheta
                    
        for j in range(numOfMeridians):
            addFace(
                0,
                vertexIndex(meridian=j, parallel=0),
                vertexIndex(meridian=j+1, parallel=0)
            )
            
        for i in range(numOfParallels-1): 
            for j in range(numOfMeridians):
                addFace(
                    vertexIndex(meridian=j, parallel=i),
                    vertexIndex(meridian=j, parallel=i+1),
                    vertexIndex(meridian=j+1, parallel=i+1),
                    vertexIndex(meridian=j+1, parallel=i)
                )
        
        for j in range(numOfMeridians):
            addFace(
                vertexIndex(meridian=j, parallel=numOfParallels-1),
                vertex.shape[0]-1,
                vertexIndex(meridian=j+1, parallel=numOfParallels-1)
            )
        
        return vertex, face
    
    def meshIcoSurface(self, nu = 1, nr_verts:int|None = None):
        
        vertex, face = icosphere(nu, nr_verts)
        
        vertex = self.radius * vertex +  self.center
        
        return vertex, face
        
    def meshSpherifiedCube(self, numOfVerticesPerEdge):
        vertex, face = Cube(2, center=(0, 0, 0)).meshSurface(
            numOfVerticesPerEdge
        )
        for i in range(len(vertex)):
            x, y, z = vertex[i]
            x2, y2, z2 = x**2, y**2, z**2
            vertex[i][0] = x * np.sqrt(1 - (y2 + z2)/2 + (y2 * z2)/3)
            vertex[i][1] = y * np.sqrt(1 - (x2 + z2)/2 + (x2 * z2)/3)
            vertex[i][2] = z * np.sqrt(1 - (x2 + y2)/2 + (x2 * y2)/3)
            
            vertex[i][0] = vertex[i][0] * self.radius + self.center[0]
            vertex[i][1] = vertex[i][1] * self.radius + self.center[1]
            vertex[i][2] = vertex[i][2] * self.radius + self.center[2]
        pass
        
        
        
        return vertex, face

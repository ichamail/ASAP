import numpy as np


class CubeVertex:
    def __init__(self, x:float, y:float, z:float, id:int) -> None:
        self.x = x
        self.y = y
        self.z = z
        self.id = id
        pass


class CubeEdge:
    def __init__(self, vertexStart:CubeVertex, vertexStop:CubeVertex, id:int) -> None:
        self.id = id
        self.vertex = np.array([vertexStart, vertexStop])
        pass
    
    def addVertices(self, numOfVertices):
                
        x = np.linspace(
            start=self.vertex[0].x,
            stop=self.vertex[-1].x,
            num=numOfVertices+2,
        )
        y = np.linspace(
            start=self.vertex[0].y,
            stop=self.vertex[-1].y,
            num=numOfVertices+2,
        )
        z = np.linspace(
            start=self.vertex[0].z,
            stop=self.vertex[-1].z,
            num=numOfVertices+2,
        )
        
        self.vertex = np.array(
            [
                self.vertex[0],
                *[
                    CubeVertex(
                        x=x[i+1], y=y[i+1], z=z[i+1],
                        id=8+self.id*numOfVertices + i
                    )
                    for i in range(numOfVertices)
                ],
                self.vertex[-1]
            ]
        )
        pass
    
        
class CubeFace:
    def __init__(
        self,
        leftEdge:CubeEdge, rightEdge:CubeEdge, upperEdge:CubeEdge, lowerEdge:CubeEdge, id:int
    ) -> None:
        self.id = id
        
        self.leftEdge = leftEdge
        self.rightEdge = rightEdge
        self.upperEdge = upperEdge
        self.lowerEdge = lowerEdge
                        
        self.setVertex()
        pass
    
    def setVertex(self):
        self.numOfVerticesPerEdge = len(self.upperEdge.vertex)
        
        self.vertex = np.zeros(
            shape=(self.numOfVerticesPerEdge, self.numOfVerticesPerEdge),
            dtype=CubeVertex
        )                
        if (
            self.upperEdge.vertex[0].id == self.leftEdge.vertex[0].id
            and
            self.upperEdge.vertex[-1].id == self.rightEdge.vertex[0].id
            and
            self.lowerEdge.vertex[0].id == self.leftEdge.vertex[-1].id
            and
            self.lowerEdge.vertex[-1].id == self.rightEdge.vertex[-1].id
        ):
            
            for i in range(self.numOfVerticesPerEdge):
                self.vertex[0][i] = self.upperEdge.vertex[i]
                self.vertex[-1][i] = self.lowerEdge.vertex[i]
                self.vertex[i][0] = self.leftEdge.vertex[i]
                self.vertex[i][-1] = self.rightEdge.vertex[i]
        elif (
            self.upperEdge.vertex[-1].id == self.rightEdge.vertex[0].id
            and
            self.rightEdge.vertex[-1].id == self.lowerEdge.vertex[0].id
            and
            self.lowerEdge.vertex[-1].id == self.leftEdge.vertex[0].id
            and
            self.leftEdge.vertex[-1].id == self.upperEdge.vertex[0].id
        ):
            for i in range(self.numOfVerticesPerEdge):
                self.vertex[0][i] = self.upperEdge.vertex[i]
                self.vertex[i][-1] = self.rightEdge.vertex[i] 
                self.vertex[-1][i] = self.lowerEdge.vertex[-1-i]
                self.vertex[i][0] = self.leftEdge.vertex[-1-i]
        else:
            for i in range(self.numOfVerticesPerEdge):
                self.vertex[0][i] = self.upperEdge.vertex[-1-i]
                self.vertex[i][-1] = self.rightEdge.vertex[-1-i] 
                self.vertex[-1][i] = self.lowerEdge.vertex[i]
                self.vertex[i][0] = self.leftEdge.vertex[i]
        pass
                           
    def addVertices(self):
        
        self.setVertex()
        
        di = (
            self.vertex[1][0].x - self.vertex[0][0].x,
            self.vertex[1][0].y - self.vertex[0][0].y,
            self.vertex[1][0].z - self.vertex[0][0].z
        )
        
        dj = (
            self.vertex[0][1].x - self.vertex[0][0].x,
            self.vertex[0][1].y - self.vertex[0][0].y,
            self.vertex[0][1].z - self.vertex[0][0].z
        )
        
        id = (
            8 + 12*(self.numOfVerticesPerEdge-2) 
            + self.id * (
                (self.numOfVerticesPerEdge-2) * (self.numOfVerticesPerEdge-2)
            )
        )
        
        for i in range(1, self.numOfVerticesPerEdge-1):
            
            for j in range(1, self.numOfVerticesPerEdge-1):
                
                dx = i * di[0] + j * dj[0]
                dy = i * di[1] + j * dj[1]
                dz = i * di[2] + j * dj[2]
                
                x = self.vertex[0][0].x + dx
                y = self.vertex[0][0].y + dy
                z = self.vertex[0][0].z + dz
                
                self.vertex[i][j] = CubeVertex(x, y, z, id)
                
                id = id+1
        pass
        

class Cube:
    
    numOfCubeEdges = 12
    numOfCubeFaces = 6
    numOfCubeVertices = 8
    
    def __init__(self, edgeLength:float, center:tuple[float]=(0, 0, 0)) -> None:
        
        self.edgeLength = edgeLength
        
        xPlus, xMinus = center[0] + edgeLength/2, center[0] - edgeLength/2
        yPlus, yMinus = center[1] + edgeLength/2, center[1] - edgeLength/2
        zPlus, zMinus = center[2] + edgeLength/2, center[2] - edgeLength/2
        
        self.vertex = np.array(
            [
                CubeVertex(x=xMinus, y=yPlus, z=zMinus, id=0),
                CubeVertex(x=xPlus, y=yPlus, z=zMinus, id=1),
                CubeVertex(x=xPlus, y=yMinus, z=zMinus, id=2),
                CubeVertex(x=xMinus, y=yMinus, z=zMinus, id=3),
                CubeVertex(x=xMinus, y=yPlus, z=zPlus, id=4),
                CubeVertex(x=xPlus, y=yPlus, z=zPlus, id=5),
                CubeVertex(x=xPlus, y=yMinus, z=zPlus, id=6),
                CubeVertex(x=xMinus, y=yMinus, z=zPlus, id=7),
                
            ]
        )
        
        self.edge = np.array(
            [
                *[
                    CubeEdge(
                        vertexStart=self.vertex[i],
                        vertexStop=self.vertex[(i+1)%4],
                        id=i
                    )
                    for i in range(4)
                ],
                *[
                    CubeEdge(
                        vertexStart=self.vertex[i + 4],
                        vertexStop=self.vertex[(i+1)%4 + 4],
                        id=i+4
                    )
                    for i in range(4)
                ],
                *[
                    CubeEdge(
                        vertexStart=self.vertex[i%4],
                        vertexStop=self.vertex[i%4+4],
                        id=i%4 + 8
                    )
                    for i in range(4)
                ]
            ]
        )
        
        self.face = np.array(
            [
                *[
                    CubeFace(
                        leftEdge=self.edge[i+8],
                        rightEdge=self.edge[(i+1)%4 + 8],
                        upperEdge=self.edge[i],
                        lowerEdge=self.edge[i+4],
                        id=i
                    )
                    for i in range(4)  
                ],
                CubeFace(
                    leftEdge=self.edge[3],
                    rightEdge=self.edge[1],
                    upperEdge=self.edge[2],
                    lowerEdge=self.edge[0],
                    id=4
                ),
                CubeFace(
                    leftEdge=self.edge[7],
                    rightEdge=self.edge[5],
                    upperEdge=self.edge[4],
                    lowerEdge=self.edge[6],
                    id=5
                )
            ]
        )
        pass
          
    def meshSurface(self, numOfVerticesPerEdge):
           
        for i in range(self.numOfCubeEdges):
            self.edge[i].addVertices(numOfVerticesPerEdge-2)        
            
        for i in range(self.numOfCubeFaces):
            
            self.face[i].addVertices()

        numOfVertices = (
                self.numOfCubeVertices 
                + self.numOfCubeEdges * (numOfVerticesPerEdge - 2)
                + self.numOfCubeFaces * (numOfVerticesPerEdge - 2)**2
        )
        numOfFaces = self.numOfCubeFaces * (numOfVerticesPerEdge - 1)**2
        
        vertex = np.zeros(shape=(numOfVertices, 3), dtype=float)
        face = - np.ones(shape=(numOfFaces, 4), dtype=int)
        
        for k in range(self.numOfCubeFaces):
            for i in range(self.face[k].numOfVerticesPerEdge):
                for j in range(self.face[k].numOfVerticesPerEdge):
                    id = self.face[k].vertex[i][j].id
                    vertex[id][0] = self.face[k].vertex[i][j].x
                    vertex[id][1] = self.face[k].vertex[i][j].y
                    vertex[id][2] = self.face[k].vertex[i][j].z
                    
        l = 0
        for k in range(self.numOfCubeFaces):
            for i in range(self.face[k].numOfVerticesPerEdge-1):
                for j in range(self.face[k].numOfVerticesPerEdge-1):
                    face[l][0] = self.face[k].vertex[i][j].id
                    face[l][1] = self.face[k].vertex[i+1][j].id
                    face[l][2] = self.face[k].vertex[i+1][j+1].id
                    face[l][3] = self.face[k].vertex[i][j+1].id
                    l = l + 1
                    
        return vertex, face
                    

if __name__=="__main__":
    from matplotlib import pyplot as plt
    from mesh_class import Mesh
    
    
    Mesh(
        *Cube(2, center=(0, 0, 0)).meshSurface(numOfVerticesPerEdge=5)
    ).display()
    
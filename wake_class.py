from vector_class import Vector
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LightSource as LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plot_functions import set_axes_equal
from panel_class import WakeQuadPanel, WakeTriPanel


"""
F:X,Y,Z (e_X, e_Y, e_Z): inertial frame of reference
f:x,y,z (e_x, e_y, e_z): translating frame of reference
(e_X = e_x , e_Y = e_y, e_Z = e_z)

f':x',y',z' (e_x', e_y', e_z'): body-fixed frame of reference

    [e_x * e_x', e_y * e_x', e_z * e_x'
A =  e_x * e_y', e_y * e_y', e_z * e_y'
     e_x * e_z', e_y * e_z', e_z * e_z']
"""


class WakeLine:
    
    def __init__(self, ro:Vector, length:float, numOfVertices:int) -> None:
        
                
        self.ro = ro
        self.Vo = Vector(0, 0, 0)
        self.A = np.identity(3, dtype=float)
                
        self.vertex = np.column_stack(
            (
                np.linspace(0, length, numOfVertices),
                np.zeros(numOfVertices, dtype=float),
                np.zeros(numOfVertices, dtype=float)
            )
        )
        
        self.numOfVertices = numOfVertices
        
        pass
        
    def setRefFrameOrientation(self, theta_x, theta_y, theta_z):
        
        theta_x = np.deg2rad(theta_x)
        theta_y = np.deg2rad(theta_y)
        theta_z = np.deg2rad(theta_z)
        
        Ax = np.array([[1, 0, 0],
                       [0, np.cos(theta_x), np.sin(theta_x)],                [0, - np.sin(theta_x), np.cos(theta_x)]])
        
        Ay = np.array([[np.cos(theta_y), 0, -np.sin(theta_y)],
                       [0, 1, 0],
                       [np.sin(theta_y), 0, np.cos(theta_y)]])
        
        Az = np.array([[np.cos(theta_z), np.sin(theta_z), 0],
                       [-np.sin(theta_z), np.cos(theta_z), 0],
                       [0, 0, 1]])
        
        self.A = Az @ Ay @ Ax
        
    def getVertexPositionVector(self, index:int):
        return self.ro + Vector(*self.vertex[index]).changeBasis(self.A.T)
        
    def getVertex(self, index:int):
        positionVector = self.getVertexPositionVector(index)
        return positionVector.x, positionVector.y, positionVector.z
        
    def addVertex(self, x, y, z):
        
        positionVector = (Vector(x, y, z) - self.ro).changeBasis(self.A)
        
        self.vertex = np.vstack(
            (
                self.vertex, 
                [positionVector.x, positionVector.y, positionVector.z]
            )
        )
        
        self.numOfVertices = self.numOfVertices + 1

    def addTrailingVertex(self, x, y, z):
    
        positionVector = (Vector(x, y, z) - self.ro).changeBasis(self.A)
        
        self.vertex = np.vstack(
            ( 
                [positionVector.x, positionVector.y, positionVector.z],
                self.vertex
            )
        )
        
        self.numOfVertices = self.numOfVertices + 1

    def moveWakeFixedFrame(self, dt):
        self.ro = self.ro + dt*self.Vo
        
    def moveVertex(self, vertexIndex:int, dr:Vector):
        dr = dr.changeBasis(self.A)
        self.vertex[vertexIndex][0] = self.vertex[vertexIndex][0] + dr.x
        self.vertex[vertexIndex][1] = self.vertex[vertexIndex][1] + dr.y
        self.vertex[vertexIndex][2] = self.vertex[vertexIndex][2] + dr.z
    

class WakeRow:
    
    def __init__(self, wakeLines:list[WakeLine], faceType:str="Quads") -> None:
        
        self.wakeLine = np.array(wakeLines)
        self.faceType = faceType
        
    @property
    def numOfFaces(self):
        if self.faceType=="Quads":
            return self.wakeLine[0].numOfVertices - 1
        elif self.faceType=="Trias":
            return (self.wakeLine[0].numOfVertices - 1) * 2              
       
    def getFace(self, index:int):
        
        if self.faceType=="Quads":
            
            return self.getQuadFace(index)
        
        elif self.faceType=="Trias":
            print("Not implementedye")
            return -1   
        
    def getQuadFace(self, index:int):
        
        if 0 <= index < self.numOfFaces:
            
            return np.array(
                [
                    self.wakeLine[0].getVertex(index),
                    self.wakeLine[1].getVertex(index),
                    self.wakeLine[1].getVertex(index+1),
                    self.wakeLine[0].getVertex(index+1)
                ]
            )
            
        elif index == -1:
            
            return np.array(
                [
                    self.wakeLine[0].getVertex(self.numOfFaces - 1),
                    self.wakeLine[1].getVertex(self.numOfFaces - 1),
                    self.wakeLine[1].getVertex(self.numOfFaces),
                    self.wakeLine[0].getVertex(self.numOfFaces)
                ]
            )
            
        else:
            
            print("index out of bounds")


class WakePanelRow(WakeRow):
    
    def __init__(self, wakeLines: list[WakeLine], faceType: str = "Quads") -> None:
        super().__init__(wakeLines, faceType)
        
        self.setPanels()
        
        pass
    
    def setPanels(self):
        if self.faceType == "Quads":
            self.panel = np.array(
                [
                    WakeQuadPanel(vertices=self.getFace(i), CCW=True)
                    for i in range(self.numOfFaces)
                ]
            )
        elif self.faceType == "Trias":
            self.panel = np.array(
                [
                    WakeTriPanel(vertices=self.getFace(i), CCW=True)
                    for i in range(self.numOfFaces)
                ]
            )
        
        pass
    
    def addTrailingPanel(self):
        
        if self.faceType == "Quads":
            
            self.panel = np.hstack(
                (
                   WakeQuadPanel(vertices=self.getFace(0), CCW=True),
                   self.panel
                )
            )
            
        elif self.faceType == "Trias":
            
            self.panel = np.hstack(
                (
                    [
                        WakeTriPanel(vertices=self.getFace(0), CCW=True),
                        WakeTriPanel(vertices=self.getFace(1), CCW=True)
                    ],
                    self.panel
                )
            )
    
    def updatePanelsPosition(self):
        for i in range(self.numOfFaces):
            self.panel[i].set_vertices(vertices=self.getFace(i), CCW=True)
                    
class Wake:
    
    def __init__(
        self, trailingEdgeVertex, length, numOfWakeFaces, faceType="Quads"
    ) -> None:
        
        self.setWakeLines(trailingEdgeVertex, length, numOfWakeFaces)
        
        self.numOfWakeLines = len(self.wakeLine)
        
        self.numOfWakeRows = self.numOfWakeLines - 1
        
        self.setWakeRows(faceType)
        
        V_inf = Vector(0, 0 , 0)
        
        pass
        
    @property
    def numOfWakeVertices(self):
        # number of wake vertices per wake line
        return self.wakeLine[0].numOfVertices
               
    @property
    def numOfWakeFaces(self):
        # number of wake faces per wake row
        return self.wakeRow[0].numOfFaces
    
    @property
    def faceType(self):
        return self.wakeRow[0].faceType
    
    @property
    def Vinf(self):
        return self.wakeLine[0].Vo
    
    @Vinf.setter
    def Vinf(self, Vinf):
        for i in range(self.numOfWakeLines):
            self.wakeLine[i].Vo = Vector(Vinf.x, Vinf.y, Vinf.z)
    
    def setWakeLines(self, trailingEdgeVertex, length, numOfWakeFaces):
        
        self.wakeLine = np.array(
            [
                WakeLine(
                    ro=Vector(*trailingEdgeVertex[i]),
                    length=length, numOfVertices = numOfWakeFaces + 1
                )
                for i in range(len(trailingEdgeVertex))
            ]
        )
    
    def setWakeRows(self, faceType):
                    
        self.wakeRow = np.array(
            [
                WakeRow([self.wakeLine[i], self.wakeLine[i+1]], faceType)
                for i in range(self.numOfWakeRows)
            ]
        )
         
    def getFace(self, wakeRowIndex:int, faceIndex:int):
        
        if self.faceType == "Quads":
            
            return self.getQuadFace(wakeRowIndex, faceIndex)

    def getQuadFace(self, wakeRowIndex:int, faceIndex:int):
        
        if (
            0 <= wakeRowIndex < self.numOfWakeRows
            and 0 <= faceIndex < self.numOfWakeFaces
        ):
            # just for readability
            wakeLineIndex = wakeRowIndex 
            VertexIndex = faceIndex
            
            return [
                (wakeLineIndex, VertexIndex),
                (wakeLineIndex+1, VertexIndex),
                (wakeLineIndex+1, VertexIndex+1),
                (wakeLineIndex, VertexIndex+1)
            ]
    
    def getFaceVertices(self, wakeRowIndex:int, faceIndex:int):
        return self.wakeRow[wakeRowIndex].getFace(faceIndex)

    def getFacesVertices(self):
        
        # # straight forward way but slower
        # return [
        #     self.getFaceVertices(wakeRowIndex=j, faceIndex=i)
        #     for j in range(self.numOfWakeRows)
        #     for i in range(self.numOfWakeFaces)
        # ]
        
        vertex = np.zeros(
            shape=(self.numOfWakeLines, self.numOfWakeVertices),
            dtype=tuple
        )

        for i in range(self.numOfWakeLines):
            for j in range(self.numOfWakeVertices):
                vertex[i][j] = self.wakeLine[i].getVertex(j)
        
        return [
            np.array(
                [
                    vertex[i][j]
                    for i, j in self.getFace(wakeRowIndex=I, faceIndex=J)
                ]
            )
            for I in range(self.numOfWakeRows)
            for J in range(self.numOfWakeFaces)
        ]

    def plot(self, elevation=30, azimuth=-60):
        
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        
        # unit vectors of body-fixed frame of reference f'
        e_x = Vector(1, 0, 0)
        e_y = Vector(0, 1, 0)
        e_z = Vector(0, 0, 1)
        
        # plot wake surface in body-fixed frame of reference f'
        ax.add_collection(
            Poly3DCollection(
                self.getFacesVertices(),
                facecolors = "steelblue",
                edgecolors = "black",
                alpha=0.5
            )
        )
        
        # plot body-fixed frame of reference f'
        ax.quiver(
            0, 0, 0, e_x.x, e_x.y, e_x.z,
            color="slateblue", label="$e_{x'}$"
        )
        ax.quiver(
            0, 0, 0, e_y.x, e_y.y, e_y.z,
            color="darkseagreen", label="$e_{y'}$"
        )
        ax.quiver(
            0, 0, 0, e_z.x, e_z.y, e_z.z,
            color="darkorange", label="$e_{z'}$"
        )
        
        
        ax.view_init(elevation, azimuth)
        ax.set_xlabel("$x'$")
        ax.set_ylabel("$y'$")
        ax.set_zlabel("$z'$")
        
        
        ax.set_xlim3d(
            min(
                [
                    0,
                    self.wakeLine[0].getVertex(0)[0],
                    self.wakeLine[-1].getVertex(0)[0]
                ]
            ),
            max(
                [
                    0,
                    self.wakeLine[0].getVertex(-1)[0],
                    self.wakeLine[-1].getVertex(-1)[0]
                ]
            )
        )
        
        ax.set_ylim3d(
            min(
                [
                    0,
                    self.wakeLine[0].getVertex(0)[1],
                    self.wakeLine[-1].getVertex(0)[1]
                ]
            ),
            max(
                [
                    0,
                    self.wakeLine[0].getVertex(-1)[1],
                    self.wakeLine[-1].getVertex(-1)[1]
                ]
            )
        )
        
        ax.set_zlim3d(
            -1+min(
                [
                    0,
                    self.wakeLine[0].getVertex(0)[2],
                    self.wakeLine[-1].getVertex(0)[2]
                ]
            ),
            1+max(
                [
                    0,
                    self.wakeLine[0].getVertex(-1)[2],
                    self.wakeLine[-1].getVertex(-1)[2]
                ]
            )
        )
        
        set_axes_equal(ax)
        
        return ax, fig
    
    def display(self, elevation=30, azimuth=-60):
        ax, fig = self.plot(elevation, azimuth)
        plt.show()

    def shed(self, trailingEdgeVertex):
        for i in range(self.numOfWakeLines):
            self.wakeLine[i].addTrailingVertex(*trailingEdgeVertex[i])     

    def moveWakeFixedFrames(self, dt):
        for i in range(self.numOfWakeLines):
            self.wakeLine[i].moveWakeFixedFrame(dt)


class PanelWake(Wake):
    
    def setWakeRows(self, faceType):
        
        self.wakeRow = np.array(
            [
                WakePanelRow([self.wakeLine[i], self.wakeLine[i+1]], faceType)
                for i in range(self.numOfWakeRows)
            ]
        )
        
        pass
    
    def shed(self, trailingEdgeVertex):
        
        super().shed(trailingEdgeVertex)
        
        for i in range(self.numOfWakeRows):
                
            self.wakeRow[i].addTrailingPanel()
            
        pass
    
if __name__=="__main__":
    
    wake = Wake(
        trailingEdgeVertex=np.array(
            [
                [2, i, 1] for i in range(-2, 2, 1)
            ]
        ),
        length=10,
        numOfWakeFaces=5,
        faceType="Quads",
    )
    
    wake.display(elevation=30, azimuth=-60)
    
    pass

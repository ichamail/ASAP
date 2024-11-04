import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LightSource as LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import stl
from src.utilities import set_axes_equal
from .panel_class import SurfaceQuadPanel, SurfaceTriPanel


"""
F:X,Y,Z (e_X, e_Y, e_Z): inertial frame of reference
f:x,y,z (e_x, e_y, e_z): translating frame of reference
(e_X = e_x , e_Y = e_y, e_Z = e_z)

f':x',y',z' (e_x', e_y', e_z'): body-fixed frame of reference

    [e_x * e_x', e_y * e_x', e_z * e_x'
A =  e_x * e_y', e_y * e_y', e_z * e_y'
     e_x * e_z', e_y * e_z', e_z * e_z']
"""


class Mesh:
    
    def __init__(
        self, vertex:np.ndarray[(any, 3), np.dtype:float],
        face:np.ndarray[np.dtype:int]
    ) -> None:
        
                
        self.vertex = np.copy(vertex)
        self.face = np.copy(face)
        
        self.numOfVertices = len(self.vertex)
        self.numOfFaces = len(self.face)
        
        self.adjacencyMatrix:np.ndarray[np.dtype:int]
        self.setAdjacencyMatrix()
        
        pass
    
    def getFace(self, face_id:int) -> list[int]:
        return [vertex_id for vertex_id in self.face[face_id] if vertex_id!=-1]
    
    def getFaceVertices(self, face_id:int) -> list[np.ndarray[float]]:
        return np.array(
            [ self.vertex[vertex_id] for vertex_id in self.getFace(face_id) ]
        )
    
    def getFacesVertices(self) -> list[np.ndarray[float]]:
        
        return [
            self.getFaceVertices(face_id=i) 
            for i in range(self.numOfFaces)
        ]
        
    @staticmethod
    def doIntersect(face_i, face_j) -> bool:
        
        # intersections = 0
        # for vertex_id in face_i:
        #     if vertex_id in face_j:
        #         intersections = intersections + 1
                
        # if intersections > 1 :
        #     return True
        # else:
        #     return False
        
        return sum(vertex_id in face_j for vertex_id in face_i)>1  

    def setAdjacencyMatrix(self) -> None:
        
        self.adjacencyMatrix = np.zeros((self.numOfFaces, self.numOfFaces))
        
        for face_i_id in range(self.numOfFaces):
            for face_j_id in range(self.numOfFaces):
                if (
                    face_i_id != face_j_id 
                    and self.doIntersect(
                        self.getFace(face_i_id), self.getFace(face_j_id)
                    )
                ):
                    self.adjacencyMatrix[face_i_id][face_j_id] = 1
                                                
        # # should be faster but when meassured it isn't               
        # for face_i_id in range(self.numOfFaces):
        #     for face_j_id in range(self.numOfFaces):
        #         if self.adjacency_matrix[face_i_id][face_j_id] == 0:
        #             if (
        #                 face_i_id != face_j_id
        #                 and self.doIntersect(
        #                     self.getFace(face_i_id), self.getFace(face_j_id)
        #                 )
        #             ):
        #                 self.adjacency_matrix[face_i_id][face_j_id] = 1
        #                 self.adjacency_matrix[face_j_id][face_i_id] = 1
        
        pass
    
    def setVSAeroAdjacencyMatrix(
        self, numOfChordWiseFaces:int, numOfSpanWiseFaces:int
    ):
        numOfChordWiseFaces = 2 * numOfChordWiseFaces
        numOfSpanWiseFaces = 2 * numOfSpanWiseFaces
                
        self.adjacencyMatrix = np.zeros((self.numOfFaces, self.numOfFaces))
        
        leftTipStartingIndex = numOfChordWiseFaces*numOfSpanWiseFaces
        rightTipStartingIndex = leftTipStartingIndex + (self.numOfFaces - leftTipStartingIndex)//2
        
        # left wing tip
        for face_i_id in range(leftTipStartingIndex, rightTipStartingIndex):
            for face_j_id in range(leftTipStartingIndex, rightTipStartingIndex):
                if (
                    face_i_id != face_j_id 
                    and self.doIntersect(
                        self.getFace(face_i_id), self.getFace(face_j_id)
                    )
                ):
                    self.adjacencyMatrix[face_i_id][face_j_id] = 1

        # right wingtip
        for face_i_id in range(rightTipStartingIndex, self.numOfFaces):
            for face_j_id in range(rightTipStartingIndex, self.numOfFaces):
                if (
                    face_i_id != face_j_id 
                    and self.doIntersect(
                        self.getFace(face_i_id), self.getFace(face_j_id)
                    )
                ):
                    self.adjacencyMatrix[face_i_id][face_j_id] = 1
        
        
        def faceId(chordWiseIndex, spanWiseIndex):
            return spanWiseIndex + chordWiseIndex*numOfSpanWiseFaces
        
        
        for i in range(1, numOfChordWiseFaces-1):
            for j in range(1, numOfSpanWiseFaces-1):
                
                id = faceId(chordWiseIndex=i, spanWiseIndex=j)
                
                self.adjacencyMatrix[id][
                    faceId(chordWiseIndex=i+1, spanWiseIndex=j)
                ] = 1
                
                self.adjacencyMatrix[id][
                    faceId(chordWiseIndex=i-1, spanWiseIndex=j)
                ] = 1
                
                self.adjacencyMatrix[id][
                    faceId(chordWiseIndex=i, spanWiseIndex=j+1)
                ] = 1
                
                self.adjacencyMatrix[id][
                    faceId(chordWiseIndex=i, spanWiseIndex=j-1)
                ] = 1
        
        for i in range(1, numOfChordWiseFaces-1):
            
            j = 0
            
            id = faceId(chordWiseIndex=i, spanWiseIndex=j)
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i+1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i-1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j+1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j+2)
            ] = 1   
            
            
            j=numOfSpanWiseFaces-1
            
            id = faceId(chordWiseIndex=i, spanWiseIndex=j)
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i+1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i-1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j-1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j-2)
            ] = 1
            
        for j in range(1, numOfSpanWiseFaces-1):
            
            i=0
            
            id = faceId(chordWiseIndex=i, spanWiseIndex=j)
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j+1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j-1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i+1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i+2, spanWiseIndex=j)
            ] = 1
            
            i=numOfChordWiseFaces-1
            
            id = faceId(chordWiseIndex=i, spanWiseIndex=j)
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j+1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i, spanWiseIndex=j-1)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i-1, spanWiseIndex=j)
            ] = 1
            
            self.adjacencyMatrix[id][
                faceId(chordWiseIndex=i-2, spanWiseIndex=j)
            ] = 1
                
        i, j = 0, 0
        id = faceId(chordWiseIndex=i, spanWiseIndex=j)
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i+1, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i+2, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j+1)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j+2)]=1
        
        i, j = numOfChordWiseFaces-1, 0
        id = faceId(chordWiseIndex=i, spanWiseIndex=j)
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i-1, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i-2, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j+1)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j+2)]=1
        
        i, j = 0, numOfSpanWiseFaces-1
        id = faceId(chordWiseIndex=i, spanWiseIndex=j)
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i+1, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i+2, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j-1)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j-2)]=1
        
        i, j = numOfChordWiseFaces -1, numOfSpanWiseFaces-1
        id = faceId(chordWiseIndex=i, spanWiseIndex=j)
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i-1, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i-2, spanWiseIndex=j)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j-1)]=1
        self.adjacencyMatrix[id][faceId(chordWiseIndex=i, spanWiseIndex=j-2)]=1
        
        pass
        
    def getFacesAdjacentFaces(self, faceIndex):
        return [
            index for index in range(self.numOfFaces)
            if self.adjacencyMatrix[faceIndex][index] == 1
        ]
    
    def __copy__(self):
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj
    
    def copy(self):
        return self.__copy__()

    def plot(self, elevation=30, azimuth=-60):
        
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.add_collection(
            Poly3DCollection(
                verts=self.getFacesVertices(),
                shade=True,
                lightsource=LightSource(azdeg=315, altdeg=45),
                facecolors = "gray",
                edgecolors = "black"
            )
        )
        
        ax.view_init(elevation, azimuth)
        ax.set_xlabel("$x'$")
        ax.set_ylabel("$y'$")
        ax.set_zlabel("$z'$")
        
        ax.set_xlim3d(self.vertex[:, 0].min(), self.vertex[:, 0].max())
        ax.set_ylim3d(self.vertex[:, 1].min(), self.vertex[:, 1].max())
        ax.set_zlim3d(self.vertex[:, 2].min(), self.vertex[:, 2].max())
        
        set_axes_equal(ax)
        
        return ax, fig
    
    def display(self, elevation=30, azimuth=-60):
        ax, fig = self.plot(elevation, azimuth)
        plt.show()

    @classmethod
    def generateFromSTLfile(cls, fileName:str, filePath="src/geometry/STL_models/"):
        
        fileName = filePath + fileName + ".stl"
        
        solid_stl_file = open(fileName, "r")
        solid_stl_obj = stl.read_ascii_file(solid_stl_file)

        vertex_list = []
        face_list = []

        for facet in solid_stl_obj.facets:
            for vertex in facet.vertices:
                is_vertex_in_vertex_list = False
                for vertex_else in vertex_list:
                    dx = abs(vertex[0] - vertex_else[0])
                    dy = abs(vertex[1] - vertex_else[1])
                    dz = abs(vertex[2] - vertex_else[2])
                    if dx < 10**(-6) and dy < 10**(-6) and dz < 10**(-6):
                        is_vertex_in_vertex_list = True
                        break
                if not is_vertex_in_vertex_list:
                    vertex_list.append([vertex[0], vertex[1], vertex[2]])
     
        for facet in solid_stl_obj.facets:
            face = []
            for vertex in facet.vertices:                
                for vetrex_id, vertex_else in enumerate(vertex_list):
                    dx = abs(vertex[0] - vertex_else[0])
                    dy = abs(vertex[1] - vertex_else[1])
                    dz = abs(vertex[2] - vertex_else[2])
                    if dx < 10**(-6) and dy < 10**(-6) and dz < 10**(-6):
                        break
                    
                face.append(vetrex_id)
                
            face_list.append(face)
        
        return cls(np.array(vertex_list), np.array(face_list))

class PanelMesh(Mesh):
    
    def __init__(self, vertex: np.ndarray, face: np.ndarray) -> None:
        super().__init__(vertex, face)
        
        self.panel:np.ndarray[SurfaceTriPanel|SurfaceQuadPanel] = np.array(
            [
                SurfaceTriPanel(
                    vertices=self.getFaceVertices(face_id=id), CCW=True
                )
                if len(self.getFace(face_id=id)) == 3 else
                SurfaceQuadPanel(
                    vertices=self.getFaceVertices(face_id=id), CCW=True
                )
                for id in range(self.numOfFaces)
            ]
        )

    def getPanelsAdjacentPanels(self, panelIndex):
        return [
            self.panel[id] for id in self.getFacesAdjacentFaces(panelIndex)
        ]

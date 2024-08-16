import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LightSource as LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plot_functions import set_axes_equal
from panel_class import SurfaceQuadPanel, SurfaceTriPanel, WakeQuadPanel, WakeTriPanel


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
        
        self.adjacency_matrix:np.ndarray[np.dtype:int]
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
        
        self.adjacency_matrix = np.zeros((self.numOfFaces, self.numOfFaces))
        
        for face_i_id in range(self.numOfFaces):
            for face_j_id in range(self.numOfFaces):
                if (
                    face_i_id != face_j_id 
                    and self.doIntersect(
                        self.getFace(face_i_id), self.getFace(face_j_id)
                    )
                ):
                    self.adjacency_matrix[face_i_id][face_j_id] = 1
                                                
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


class PanelMesh(Mesh):
    
    def __init__(self, vertex: np.ndarray, face: np.ndarray) -> None:
        super().__init__(vertex, face)
        
        self.panel = np.array(
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


if __name__=="__main__":
    from airfoil_class import Airfoil
    from wing_class import Wing
    
    node, face = Wing(
        root_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
        tip_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
        halfSpan=0.5
    ).meshSurface(
        numOfChordWiseFaces=5,
        numOfSpanWiseFaces=5,
        faceType="quadrilateral",
        chordWiseSpacing="cosine",
        spanWiseSpacing="uniform",
        mesh_MainSurface=True,
        mesh_WingTips=True,
    )
    
    wingMesh = Mesh(node, face)
    
    wingMesh.display(elevation=30, azimuth=-60)
    
    pass

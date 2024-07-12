from vector_class import Vector
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

    def setAdjacencyMatrix(self):
        
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
        #                 and self.do_intersect(
        #                     self.get_face(face_i_id), self.get_face(face_j_id)
        #                 )
        #             ):
        #                 self.adjacency_matrix[face_i_id][face_j_id] = 1
        #                 self.adjacency_matrix[face_j_id][face_i_id] = 1

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
                self.getFacesVertices(),
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
    
    def dispaly(self, elevation=30, azimuth=-60):
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

    
class RigidBody:
    
    def __init__(self, surfaceMesh:Mesh, name:str="rigid body") -> None:
        
        self.name = name
        self._surfaceMesh = surfaceMesh
        
        self.ro = Vector(0, 0, 0) # position vector of mesh origin
        self.A = np.identity(3, dtype=float) # orientation matrix A, A^T = R
        
        # velocity vector of body-fixed frame origin
        self.Vo = Vector(0, 0, 0)
        
        # body-fixed frame's angular velocity vector
        self.omega = Vector(0, 0, 0)
        
        self.trailingEdge = None
        self.sheddingFaces = None
        
        pass
    
    @property
    def surface(self):
        return self._surfaceMesh
    
    def getFaceVertices(self, face_id, bodyFixedFrame=True):
        
        if bodyFixedFrame:
            
            return self.surface.getFaceVertices(face_id)
        
        else:
            
            vertex = self.surface.getFaceVertices(face_id)
                        
            for i in range(len(vertex)):
                r = (
                    self.ro
                    + Vector(*vertex[i]).changeBasis(self.A.T)
                )
                vertex[i] = r.x, r.y, r.z
            
            return vertex
    
    def getFacesVertices(self, bodyFixedFrame=True):
        
        if bodyFixedFrame:
            
            return self.surface.getFacesVertices()
            
        else:
            
            vertex = np.zeros_like(self.surface.vertex)
            for id in range(self.surface.numOfVertices):
                r = (
                    self.ro 
                    + Vector(*self.surface.vertex[id]).changeBasis(self.A.T)
                )
                vertex[id] = r.x, r.y, r.z
            
            return [
                np.array(
                    [
                        vertex[vertex_id] 
                        for vertex_id in self.surface.getFace(face_id)
                    ]
                )
                for face_id in range(self.surface.numOfFaces)
            ]
            
    def set_BodyFixedFrame_origin(self, xo, yo, zo):
        self.ro = Vector(xo, yo, zo)
        
    def set_BodyFixedFrame_orientation(self, theta_x, theta_y, theta_z):
        """
            A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
            A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) 
                        
            where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
            
            if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
            initial orientation of the body is such that f' coincides with f
            
            if A(t=0) =/= I then
            A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx) 
                        
            where Δθx, Δθy, and Δθz correspond to finite rotations
            
            https://en.wikipedia.org/wiki/Infinitesimal_rotation_matrix
            https://en.wikipedia.org/wiki/Davenport_chained_rotations
        """
        
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
    
    def set_BodyFixedFrame_origin_velocity(self, Vo_x, Vo_y, Vo_z):
        self.Vo = Vector(Vo_x, Vo_y, Vo_z)
        
    def set_BodyFixedFrame_angular_velocity(self, omega_x, omega_y, omega_z):
        self.omega = Vector(
            np.deg2rad(omega_x), np.deg2rad(omega_y), np.deg2rad(omega_z)
        )
        
    def move_BodyFixedFrame(self, dt):
        """
            A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
            A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) =>
                        
            where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
            
            if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
            initial orientation of the body is such that f' coincides with f
            
            if A(t=0) =/= I then
            A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx)
            
            where Δθx, Δθy, and Δθz correspond to finite rotations
            
            https://en.wikipedia.org/wiki/Infinitesimal_rotation_matrix
            https://en.wikipedia.org/wiki/Davenport_chained_rotations
        """
        
        self.ro = self.ro + self.Vo*dt
        dtheta =  self.omega*dt
        
        Ax = np.array([[1, 0, 0],
                       [0, np.cos(dtheta.x), np.sin(dtheta.x)],                [0, - np.sin(dtheta.x), np.cos(dtheta.x)]])
        
        Ay = np.array([[np.cos(dtheta.y), 0, -np.sin(dtheta.y)],
                       [0, 1, 0],
                       [np.sin(dtheta.y), 0, np.cos(dtheta.y)]])
        
        Az = np.array([[np.cos(dtheta.z), np.sin(dtheta.z), 0],
                       [-np.sin(dtheta.z), np.cos(dtheta.z), 0],
                       [0, 0, 1]])
        
        self.A = Az @ Ay @ Ax @ self.A

    def set_trailingEdge(self, trailingEdgeVerticesIDs:list[int]):
        self.trailingEdge = np.array(trailingEdgeVerticesIDs)
       
    def locateSheddingFaces(self):
        
        self.sheddingFaces = np.zeros(
            shape=( len(self.trailingEdge) - 1, 2 ), dtype=int
        )

        for i in range(self.sheddingFaces.shape[0]):
            j = 0 
            for face_id in range(self.surface.numOfFaces):
                if sum(
                    [
                        vertex_id in [self.trailingEdge[i], self.trailingEdge[i+1]] 
                        for vertex_id in self.surface.getFace(face_id)
                    ]
                ) == 2:
                    
                    self.sheddingFaces[i][j] = face_id
                    j = j + 1
    
    def setSheddingFace(
        self, upperSheddingFaces_idList, lowerSheddingFaces_idList
        ):
        self.sheddingFaces = np.column_stack(
            (upperSheddingFaces_idList, lowerSheddingFaces_idList)
        )
    
    def plot_withRespectToBodyFixedFrame(self, elevation=30, azimuth=-60):
        
        # unit vectors of inertial frame of reference F  expressed in the body-fixed frame of reference f'
        ro = -self.ro.changeBasis(self.A)  # ro: r_oo' -> r_o'o = -roo'
        e_X = Vector(1, 0, 0).changeBasis(self.A)
        e_Y = Vector(0, 1, 0).changeBasis(self.A)
        e_Z = Vector(0, 0, 1).changeBasis(self.A)
        
                    
        # unit vectors of body-fixed frame of reference f'
        e_x = Vector(1, 0, 0)
        e_y = Vector(0, 1, 0)
        e_z = Vector(0, 0, 1)
        
        # plot surface in body-fixed frame of reference f'
        ax, fig = self.surface.plot(elevation, azimuth)
        
        # plot inertial frame of reference F in the body-fixed frame of reference f'
        ax.quiver(
            ro.x, ro.y, ro.z, e_X.x, e_X.y, e_X.z,
            color='b', label="$e_x$"
        )
        
        ax.quiver(
            ro.x, ro.y, ro.z, e_Y.x, e_Y.y, e_Y.z,
            color='g', label="$e_y$"
        )
        ax.quiver(
            ro.x, ro.y, ro.z, e_Z.x, e_Z.y, e_Z.z,
            color='r', label="$e_z$"
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
        
        ax.set_xlim3d(
            min([ax.get_xlim()[0], ro.x]),
            max([ax.get_xlim()[1], ro.x])
        )
        ax.set_ylim3d(
            min([ax.get_ylim()[0], ro.y]),
            max([ax.get_ylim()[1], ro.y])
        )
        ax.set_zlim3d(
            min([ax.get_zlim()[0], ro.z]),
            max([ax.get_zlim()[1], ro.z])
        )
        
        set_axes_equal(ax)
        
        return ax, fig
        
    def plot_withRespectToInertialFrame(self, elevation=30, azimuth=-60):
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        ax.view_init(elevation, azimuth)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_zlabel("$z$")
        
        
        # unit vectors of Inertial frame of reference F 
        e_X = Vector(1, 0, 0)
        e_Y = Vector(0, 1, 0)
        e_Z = Vector(0, 0, 1)
        
                    
        # unit vectors of body-fixed frame of reference f' expressed in the inertial frame of reference F
        e_x = Vector(1, 0, 0).changeBasis(self.A.T)
        e_y = Vector(0, 1, 0).changeBasis(self.A.T)
        e_z = Vector(0, 0, 1).changeBasis(self.A.T)
        
        # surface vertices epressed in the inertial frame of reference F
        vertex = np.zeros_like(self.surface.vertex)
        for id in range(self.surface.numOfVertices):
            r = (
                self.ro 
                + Vector(*self.surface.vertex[id]).changeBasis(self.A.T)
            )
            vertex[id] = r.x, r.y, r.z
        
        
        # plot surface in the inertial frame of reference F           
        ax.add_collection(
            Poly3DCollection(
                self.getFacesVertices(bodyFixedFrame=False),
                shade=True,
                lightsource=LightSource(azdeg=315, altdeg=45),
                facecolors = "gray",
                edgecolors = "black"
            )
        )
        
        # plot inertial frame of reference F
        e_X = ax.quiver(
            0, 0, 0, e_X.x, e_X.y, e_X.z,
            color='b', label="$e_x$"
        )
        
        e_Y = ax.quiver(
            0, 0, 0, e_Y.x, e_Y.y, e_Y.z,
            color='g', label="$e_y$"
        )
        e_Z = ax.quiver(
            0, 0, 0, e_Z.x, e_Z.y, e_Z.z,
            color='r', label="$e_z$"
        )
        
        # plot body-fixed frame of reference f' in the inertial frame of reference F
        e_x = ax.quiver(
            self.ro.x, self.ro.y, self.ro.z, e_x.x, e_x.y, e_x.z,
            color="slateblue", label="$e_{x'}$"
        )
        e_y = ax.quiver(
            self.ro.x, self.ro.y, self.ro.z, e_y.x, e_y.y, e_y.z,
            color="darkseagreen", label="$e_{y'}$"
        )
        e_z = ax.quiver(
            self.ro.x, self.ro.y, self.ro.z, e_z.x, e_z.y, e_z.z,
            color="darkorange", label="$e_{z'}$"
        )
        
        ax.set_xlim3d(min([0, *vertex[:, 0]]), max([0, *vertex[:, 0]]))
        ax.set_ylim3d(min([0, *vertex[:, 1]]), max([0, *vertex[:, 1]]))
        ax.set_zlim3d(min([0, *vertex[:, 2]]), max([0, *vertex[:, 2]]))
                    
        set_axes_equal(ax)
        
        return ax, fig
    
    def plot(self, elevation=30, azimuth=-60, bodyFixedFrame=False):
        if bodyFixedFrame:
            return self.plot_withRespectToBodyFixedFrame(elevation, azimuth)    
        else:
            return self.plot_withRespectToInertialFrame(elevation, azimuth)
    
    def display(self, elevation=30, azimuth=-60, bodyFixedFrame=False):
        self.plot(elevation, azimuth, bodyFixedFrame)
        plt.show()
                    
    
class WakeLine:
    
    def __init__(self, ro:Vector, length:float, numOfVertices:int) -> None:
        
                
        self.ro = ro
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
        else:
            self.panel = np.array(
                [
                    WakeTriPanel(vertices=self.getFace(i), CCW=True)
                    for i in range(self.numOfFaces)
                ]
            )
                    
class Wake:
    
    def __init__(
        self, rigidBody:RigidBody, length, numOfWakeFaces, faceType="Quads"
    ) -> None:
        
        
        self.wakeLine = np.array(
            [
                WakeLine(
                    ro=Vector(*rigidBody.surface.vertex[id]),
                    length=length, numOfVertices = numOfWakeFaces + 1
                )
                for id in rigidBody.trailingEdge
            ]
        )
        
        self.numOfWakeLines = len(self.wakeLine)
        self.numOfWakeRows = self.numOfWakeLines - 1
        
        self.wakeRow = np.array(
            [
                WakeRow([self.wakeLine[i], self.wakeLine[i+1]], faceType)
                for i in range(self.numOfWakeRows)
            ]
        )
        
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
    
    def dispaly(self, elevation=30, azimuth=-60):
        ax, fig = self.plot(elevation, azimuth)
        plt.show()


class FreeWake(Wake):

    def __init__(
        self, rigidBody:RigidBody, length, numOfWakeFaces, faceType="Quads"
    ) -> None:       
        
        self.wakeLine = np.array(
            [
                WakeLine(
                    ro=rigiBody.r_o + Vector(*rigidBody.surface.vertex[id]).changeBasis(rigiBody.A.T),
                    length=length, numOfVertices = numOfWakeFaces + 1
                )
                for id in rigidBody.trailingEdge
            ]
        )
        
        self.numOfWakeLines = len(self.wakeLine)
        self.numOfWakeRows = self.numOfWakeLines - 1
        
        self.wakeRow = np.array(
            [
                WakeRow([self.wakeLine[i], self.wakeLine[i+1]], faceType)
                for i in range(self.numOfWakeRows)
            ]
        )
        
    def plot(self, elevation=30, azimuth=-60):
        
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        
        # unit vectors of inertial frame of reference f'
        e_x = Vector(1, 0, 0)
        e_y = Vector(0, 1, 0)
        e_z = Vector(0, 0, 1)
        
        # plot wake surface in inertial frame of reference f'
        ax.add_collection(
            Poly3DCollection(
                self.getFacesVertices(),
                facecolors = "steelblue",
                edgecolors = "black",
                alpha=0.5
            )
        )
        
        # plot inertial frame of reference F
        e_x = ax.quiver(
            0, 0, 0, e_x.x, e_x.y, e_x.z,
            color='b', label="$e_x$"
        )
        
        e_y = ax.quiver(
            0, 0, 0, e_y.x, e_y.y, e_y.z,
            color='g', label="$e_y$"
        )
        
        e_z = ax.quiver(
            0, 0, 0, e_z.x, e_z.y, e_z.z,
            color='r', label="$e_z$"
        )
        
        
        ax.view_init(elevation, azimuth)
        ax.set_xlabel("$x$")
        ax.set_ylabel("$y$")
        ax.set_zlabel("$z$")
        
        
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

       
    
if __name__=="__main__":
    from airfoil_class import Airfoil
    from wing_class import Wing
    
    node, face = Wing(
        root_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
        tip_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
        halfSpan=0.5
    ).meshSurface(
        numOfChordWiseFaces=5,
        numOfSpanWiseFaces=2,
        faceType="quadrilateral",
        chordWiseSpacing="cosine",
        spanWiseSpacing="uniform",
        mesh_MainSurface=True,
        mesh_WingTips=True,
    )
    
    wingMesh = Mesh(node, face)
    
    wingMesh.dispaly(elevation=30, azimuth=-60)
    
    rigidBody = RigidBody(surfaceMesh=wingMesh)
    rigidBody.set_BodyFixedFrame_origin(2, 1, 1)
    
    rigidBody.display(bodyFixedFrame=False)
    
    rigidBody.set_trailingEdge(
        trailingEdgeVerticesIDs=[id for id in range(5)]
    )
    
    wake = Wake(
        rigidBody=rigidBody, length=10, numOfWakeFaces=5, faceType="Quads"
    )
    
    
    
    wake.dispaly()

    pass

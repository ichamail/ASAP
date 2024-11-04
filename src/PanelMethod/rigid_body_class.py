import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LightSource as LightSource
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from src.utilities import set_axes_equal
from src.myMath import Vector
from .mesh_class import Mesh, PanelMesh
from .wake_class import Wake, PanelWake


"""
F:X,Y,Z (e_X, e_Y, e_Z): inertial frame of reference
f:x,y,z (e_x, e_y, e_z): translating frame of reference
(e_X = e_x , e_Y = e_y, e_Z = e_z)

f':x',y',z' (e_x', e_y', e_z'): body-fixed frame of reference

    [e_x * e_x', e_y * e_x', e_z * e_x'
A =  e_x * e_y', e_y * e_y', e_z * e_y'
     e_x * e_z', e_y * e_z', e_z * e_z']
"""


class RigidBody:
    
    def __init__(self, surfaceMesh:Mesh|PanelMesh, name:str="rigid body") -> None:
        
        self.name = name
        self._surfaceMesh = surfaceMesh
        
        self.ro = Vector(0, 0, 0) # position vector of mesh origin
        self.A = np.identity(3, dtype=float) # orientation matrix A, A^T = R
        
        # velocity vector of body-fixed frame origin
        self.Vo = Vector(0, 0, 0)
        
        # body-fixed frame's angular velocity vector
        self.omega = Vector(0, 0, 0)
        
        pass
    
    @property
    def surface(self):
        return self._surfaceMesh
    
    def getVertex(self, vertex_id, bodyFixedFrame=True):
        if bodyFixedFrame:
            return self.surface.vertex[vertex_id]
        else:
            r = (
                self.ro 
                + Vector(*self.getVertex(vertex_id)).changeBasis(self.A.T)
            )
            return r.x, r.y, r.z
            
    def getFaceVertices(self, face_id, bodyFixedFrame=True):
        
        if bodyFixedFrame:
            
            return self.surface.getFaceVertices(face_id)
        
        else:
                        
            return np.array(
                [
                    self.getVertex(vertex_id, bodyFixedFrame)
                    for vertex_id in self.surface.getFace(face_id)
                ]
            )
    
    def getFacesVertices(self, bodyFixedFrame=True):
        
        if bodyFixedFrame:
            
            return self.surface.getFacesVertices()
            
        else:
            
            # # straight forward way but slower
            # return [
            #     self.getFaceVertices(face_id, bodyFixedFrame)
            #     for face_id in range(self.surface.numOfFaces)
            # ]
            
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
        ax.legend()
        
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
                verts=self.getFacesVertices(bodyFixedFrame=False),
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
        ax.legend()
        
        return ax, fig
    
    def plot(self, elevation=30, azimuth=-60, bodyFixedFrame=False):
        if bodyFixedFrame:
            return self.plot_withRespectToBodyFixedFrame(elevation, azimuth)    
        else:
            return self.plot_withRespectToInertialFrame(elevation, azimuth)
    
    def display(self, elevation=30, azimuth=-60, bodyFixedFrame=False):
        ax, fig = self.plot(elevation, azimuth, bodyFixedFrame)
        plt.show()


class RigidAerodynamicBody(RigidBody):
    
    def __init__(self, surfaceMesh:Mesh|PanelMesh, name: str = "rigid body") -> None:
        super().__init__(surfaceMesh, name)
        self.trailingEdge = None
        self.sheddingFaces = None
        self.wake = None
        self.isWakeFree = False
    
    def setTrailingEdge(self, trailingEdgeVerticesIDs:list[int]):
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
    
    def setSheddingFaces(
        self, upperSheddingFacesIDs:list[int], lowerSheddingFacesIDs:list[int]
        ):
        self.sheddingFaces = np.column_stack(
            (upperSheddingFacesIDs, lowerSheddingFacesIDs)
        )
    
    def setWake(
        self, length:float, numOfWakeFaces:int,
        faceType:str="Quads", isWakeFree:bool=False):
                
        if isWakeFree:
            self.isWakeFree = True
            bodyFixedFrame = False
        else:
            self.isWakeFree = False
            bodyFixedFrame = True
            
        if isinstance(self.surface, PanelMesh):
            Wake = PanelWake
            
        self.wake = Wake(
            trailingEdgeVertex = np.array(
                [
                    self.getVertex(vertex_id, bodyFixedFrame)
                    for vertex_id in self.trailingEdge
                ] 
            ),
            length = length,
            numOfWakeFaces = numOfWakeFaces,
            faceType = faceType       
        )
                   
        pass
        
    def getWakeVertex(self, wakeLineIndex, vertexIndex, bodyFixedFrame=True):
        
        if (bodyFixedFrame and not self.isWakeFree) or (not bodyFixedFrame and self.isWakeFree):
            
            return self.wake.wakeLine[wakeLineIndex].getVertex(vertexIndex)
        
        elif not bodyFixedFrame and not self.isWakeFree:
            
            r = (
                self.ro 
                + self.wake.wakeLine[wakeLineIndex].getVertexPositionVector(
                    vertexIndex).changeBasis(self.A.T)
            )
            return r.x, r.y, r.z
        
        elif bodyFixedFrame and self.isWakeFree:
            r = (
                self.wake.wakeLine[wakeLineIndex].getVertexPositionVector(
                    vertexIndex
                )
                - self.ro
            ).changeBasis(self.A)
            
            return r.x, r.y, r.z
    
    def getWakeFaceVertices(self, wakeRowIndex, faceIndex, bodyFixedFrame=True):
        
        return np.array(
            [
                self.getWakeVertex(wakeLineIndex, vertexIndex, bodyFixedFrame)
                for wakeLineIndex, vertexIndex in self.wake.getFace(
                    wakeRowIndex, faceIndex
                )
            ]
        )
    
    def getWakeFacesVertices(self, bodyFixedFrame=True):
        
        # # straight forward way but slower
        # return [
        #     self.getWakeFaceVertices(wakeRowIndex=j, faceIndex=i, bodyFixedFrame=bodyFixedFrame)
        #     for j in range(self.wake.numOfWakeRows)
        #     for i in range(self.wake.numOfWakeFaces)
        # ]
        
        vertex = np.zeros(
            shape=(self.wake.numOfWakeLines, self.wake.numOfVerticesPerWakeLine),
            dtype=tuple
        )

        for wakeLineIndex in range(self.wake.numOfWakeLines):
            for vertexIndex in range(self.wake.numOfVerticesPerWakeLine):
                vertex[wakeLineIndex][vertexIndex] = self.getWakeVertex(
                    wakeLineIndex, vertexIndex, bodyFixedFrame
                )
              
        return [
            np.array(
                [
                    vertex[wakeLineIndex][vertexIndex]
                    for wakeLineIndex, vertexIndex in self.wake.getFace(
                        wakeRowIndex, faceIndex
                    )
                ]
            )
            for wakeRowIndex in range(self.wake.numOfWakeRows)
            for faceIndex in range(self.wake.numOfFacesPerWakeRow)
        ]
             
    def plot(
        self, elevation=30, azimuth=-60, bodyFixedFrame=False, plotWake=False
    ):
        
        if plotWake:
            
            ax, fig =  super().plot(elevation, azimuth, bodyFixedFrame)
            
            # plot wake surface in body-fixed frame of reference f'
            ax.add_collection(
                Poly3DCollection(
                    verts=self.getWakeFacesVertices(bodyFixedFrame),
                    facecolors = "steelblue",
                    edgecolors = "black",
                    alpha=0.5
                )
            )
            
            verts = np.array(
                [
                    self.getWakeVertex(
                        wakeLineIndex, self.wake.numOfVerticesPerWakeLine-1, bodyFixedFrame
                    )
                    for wakeLineIndex in range(self.wake.numOfWakeLines)
                ]
            )
            
            x_limits = ax.get_xlim3d()
            y_limits = ax.get_ylim3d()           
            z_limits = ax.get_zlim3d()
            
            ax.set_xlim3d(
                min([x_limits[0], *verts[:, 0]]),
                max([x_limits[1], *verts[:, 0]])
            )
            ax.set_ylim3d(
                min([y_limits[0], *verts[:, 1]]),
                max([y_limits[1], *verts[:, 1]])
            )
            ax.set_zlim3d(
                min([z_limits[0], *verts[:, 2]]),
                max([z_limits[1], *verts[:, 2]])
            )
            set_axes_equal(ax)
            
            return ax, fig
        
        else:
            
            return super().plot(elevation, azimuth, bodyFixedFrame)
    
    def display(
        self, elevation=30, azimuth=-60, bodyFixedFrame=False, displayWake=True
    ):
        ax, fig = self.plot(
            elevation, azimuth, bodyFixedFrame, plotWake=displayWake
        )
        plt.show()

    def set_BodyFixedFrame_origin(self, xo, yo, zo):
        super().set_BodyFixedFrame_origin(xo, yo, zo)
        
        if self.wake != None and self.isWakeFree:
                            
            self.wake.setWakeLinesRefFrameOrigin(
                trailingEdgeVertex=np.array(
                    [
                        self.getVertex(
                            vertex_id, bodyFixedFrame=False
                        )
                        for vertex_id in self.trailingEdge
                    ]
                )
            )
    
    def set_BodyFixedFrame_orientation(self, theta_x, theta_y, theta_z):
        
        super().set_BodyFixedFrame_orientation(theta_x, theta_y, theta_z)
        
        if self.wake != None:
            
            if self.isWakeFree:
                
                self.wake.setWakeLinesRefFrameOrientation(
                    np.identity(3)
                )
                
                self.wake.setWakeLinesRefFrameOrigin(
                    trailingEdgeVertex=np.array(
                        [
                            self.getVertex(
                                vertex_id, bodyFixedFrame=False
                            )
                            for vertex_id in self.trailingEdge
                        ]
                    )
                )
                
            else:
                
                self.wake.setWakeLinesRefFrameOrientation(self.A.T)
    
    def move_BodyFixedFrame(self, dt):
        super().move_BodyFixedFrame(dt)
        self.ravelWake(dt)
        
    def ravelWake(self, dt):
        if self.isWakeFree:
            self.ravelFreeWake(dt)
        else:
            self.ravelSurfaceFixedWake(dt)
    
    def ravelFreeWake(self, dt):
        
        # self.wake.moveWakeFixedFrames(dt) # this method updates panels' positions
        
        if self.wake.Vinf != Vector(0, 0, 0):
            self.wake.moveWakeLinesRefFrames(dt) # this method updates panels' positions
            
        self.wake.shed(
            trailingEdgeVertex=np.array(
                [
                    self.getVertex(vertex_id=id, bodyFixedFrame=False)
                    for id in self.trailingEdge
                ]
            )
        )
        
    def ravelSurfaceFixedWake(self, dt):
        
        for wakeLineIndex in range(self.wake.numOfWakeLines):
            for vertexIndex in range(self.wake.numOfVerticesPerWakeLine):
                                
                # r = Vector(
                #     *self.getWakeVertex(
                #         wakeLineIndex, vertexIndex, bodyFixedFrame=False
                #     )
                # )
                # V =  self.wake.Vinf - (
                #     self.Vo + self.omega.cross(r - self.ro) 
                # )
                
                # faster - less vector operations
                # r = Vector(
                #     *self.getWakeVertex(
                #         wakeLineIndex, vertexIndex, bodyFixedFrame=True
                #     )
                # )
                r = self.wake.wakeLine[wakeLineIndex].getVertexPositionVector(
                    vertexIndex
                )
                
                V =  self.wake.Vinf - (
                    self.Vo + self.omega.cross(r.changeBasis(self.A.T)) 
                )
                                              
                dr = dt * V.changeBasis(self.A)
                                                    
                self.wake.wakeLine[wakeLineIndex].moveVertex(vertexIndex,dr)
        
        self.wake.updatePanelsPosition() # free wake doesn't need update
                
        self.wake.shed(
            trailingEdgeVertex=np.array(
                [
                    self.getVertex(vertex_id=id, bodyFixedFrame=True)
                    for id in self.trailingEdge
                ]
            )
        )
        
        pass
    

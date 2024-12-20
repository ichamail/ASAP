import numpy as np
from matplotlib import pyplot as plt, cm
from .rigid_body_class import RigidBody, RigidAerodynamicBody
from src.myMath import Vector
from src.numbaSpeedUp import jit_computeSurfaceInfluence


class BoundaryElementMethod:
    
    def __init__(self, rigidBody:RigidBody|RigidAerodynamicBody) -> None:
        self.rigidBody = rigidBody
        
        self.Bij = np.zeros(
            shape=(self.surface.numOfFaces, self.surface.numOfFaces)
        )
        self.Cij = np.zeros_like(self.Bij)
        
        self.computeSurfaceInfluence()
        
        self.Vinf = Vector(0, 0, 0)
        
        pass
        
    def setVinf(
        self, angleOfAttack:float, sideSlipAngle:float, magnitude:float
    ):
        self.rigidBody.set_BodyFixedFrame_orientation(0, 0, 0)
        self.rigidBody.set_BodyFixedFrame_angular_velocity(0, 0, 0)
        self.rigidBody.set_BodyFixedFrame_origin_velocity(0, 0, 0)
        
        alpha, beta = np.deg2rad(angleOfAttack), np.deg2rad(sideSlipAngle)
        
        self.Vinf = magnitude * Vector(
            np.cos(alpha) * np.cos(beta),
            np.cos(alpha) * (- np.sin(beta)),
            np.sin(alpha)
        )
        
        pass
    
    def setVfs(
        self, angleOfAttack:float, sideSlipAngle:float, magnitude:float
    ):
        self.rigidBody.set_BodyFixedFrame_orientation(
            theta_x=0, theta_y=angleOfAttack, theta_z=sideSlipAngle
        )
        self.rigidBody.set_BodyFixedFrame_angular_velocity(0, 0, 0)
        
        self.rigidBody.set_BodyFixedFrame_origin_velocity(
            Vo_x=-magnitude, Vo_y=0, Vo_z=0
        )
        
        self.Vinf = Vector(0, 0, 0)
        
        pass
    
    def Vfs_atBodyFixedFrameOrigin(self):
        return (self.Vinf - self.rigidBody.Vo).changeBasis(self.rigidBody.A)
    
    def Vfs_at(self, r_rigidBodyPoint:Vector):
        
        Vfs_atPointOfRigidBody = (
            self.Vinf
            - (
                self.rigidBody.Vo
                + self.rigidBody.omega.cross(
                    r_rigidBodyPoint.changeBasis(self.rigidBody.A.T)
                )
            )
        ).changeBasis(self.rigidBody.A)
        
        return Vfs_atPointOfRigidBody
        
    @property
    def surface(self):
        return self.rigidBody.surface
    
    @property
    def surfacePanel(self):
        return self.surface.panel
    
    @property
    def Aij(self):
        return self.Cij
    
    @property
    def RHS(self):
        return -self.Bij @ np.array(
            [panel.sigma for panel in self.surfacePanel]
        )
    
    def computeSurfaceInfluence(self):
        
        for i in range(self.surface.numOfFaces):
            
            for j in range(self.surface.numOfFaces):
                
                self.Bij[i][j], self.Cij[i][j] = (
                    self.surface.panel[j].unitStrength_inducedVelocityPotential(
                        self.surface.panel[i].r_cp
                    )
                )
                if i == j:
                    print("i=", i, ", j=", j, " : ", self.Bij[i][j],self.Cij[i][j])
        pass
    
    def computeSurfaceInfluence(self):
        return jit_computeSurfaceInfluence(self)
    
    # overwrite method with numba jit compilation
    def computeSourceStrengths(self):
        if self.rigidBody.omega.norm() != 0:
            for panel in self.surfacePanel:
                panel.sigma = - panel.e_n.dot(self.Vfs_at(panel.r_centroid))
        else:
            Vfs = self.Vfs_atBodyFixedFrameOrigin()
            for panel in self.surfacePanel:
                panel.sigma = - panel.e_n.dot(Vfs)
        pass
    
    def solveLinearSystem(self):
        mu = np.linalg.solve(self.Aij, self.RHS)
        for i in range(self.surface.numOfFaces):
            self.surface.panel[i].mu = mu[i]
        pass
    
    def advanceSolution(self):
        self.computeSourceStrengths()
        self.solveLinearSystem()
    
    def computeSurfaceVelocity(self):
        
        if self.rigidBody.omega.norm() != 0:
            
            for i in range(self.surface.numOfFaces):
                self.surface.panel[i].setSurfaceVelocity(
                    self.surface.getPanelsAdjacentPanels(panelIndex=i),
                    self.Vfs_at(self.surface.panel[i].r_centroid)
                )        
        else:
            
            Vfs = self.Vfs_atBodyFixedFrameOrigin()
            for i in range(self.surface.numOfFaces):
                self.surface.panel[i].setSurfaceVelocity(
                    self.surface.getPanelsAdjacentPanels(panelIndex=i), Vfs
                )
        pass
    
    def computeSurfacePressure(self):
        if self.rigidBody.omega.norm() != 0:
            for panel in self.surfacePanel:
                Vfs = self.Vfs_at(panel.r_centroid)
                panel.Cp = 1 - ( panel.V.norm()/Vfs.norm() )**2
        else:
            Vfs = self.Vfs_atBodyFixedFrameOrigin()
            for panel in self.surfacePanel:
                panel.Cp = 1 - ( panel.V.norm()/Vfs.norm() )**2
        pass
    
    def solve(self):
        self.advanceSolution()
        self.computeSurfaceVelocity()
        self.computeSurfacePressure()
        pass

    def inducedVelocity(self, r_p:Vector) -> Vector:
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: induced velocity vector of point P with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        Vinduced = Vector(0, 0, 0)
        for panel in self.surfacePanel:
            Vinduced = Vinduced + panel.inducedVelocity(r_p)
        
        return Vinduced
    
    def inducedVelocityPotential(self, r_p:Vector) -> float:
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            float: induced velocity potential at P
        """
        
        phi_induced = 0
        for panel in self.surfacePanel:
            phi_induced = phi_induced + panel.inducedVelocityPotential(r_p)
        
        return phi_induced
    
    def velocity(self, r_p:Vector) -> Vector:
        """_summary_

        Args:
            r_p (Vector): position vector of point P with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: velocity vector of point P with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        if self.rigidBody.omega.norm() != 0:
            Vfs = self.Vfs_at(r_p)                
        else:
            Vfs = self.Vfs_atBodyFixedFrameOrigin()               
              
        return Vfs + self.inducedVelocity(r_p)
    
    def velocityPotential(self, r_p:Vector) -> float:
        """_summary_

        Args:
            r_p (Vector): position vector of point with respect to body-fixed frame of reference. Coordinates of r_p are expressed, also, with respect to the basis of body-fixed of reference

        Returns:
            Vector: velocity potential with respect to body-fixed frame of reference. Coordinates of velocity vector are expressed, also, with respect to the basis of body-fixed of reference
        """
        
        if self.rigidBody.omega.norm() != 0:
            print("error: unknown phi_infty when body is rotating")
            Vfs = self.Vfs_at(r_p)
            phi_fs = r_p.x * Vfs.x + r_p.y * Vfs.y                
        else:
            Vfs = self.Vfs_atBodyFixedFrameOrigin()
            phi_fs = r_p.x * Vfs.x + r_p.y * Vfs.y                
              
        return phi_fs + self.inducedVelocityPotential(r_p)
    
    def plotSurfacePressureCoefficientContour(
        self, elevation=30, azimuth=-60, bodyFixedFrame=False
    ):
        ax, fig = self.rigidBody.plot(elevation, azimuth, bodyFixedFrame)
        Cp = [panel.Cp for panel in self.surfacePanel]
        CpMin, CpMax = min(Cp), max(Cp)
        NormCp = [(float(Cp_i)-CpMin)/(CpMax-CpMin) for Cp_i in Cp]
        facecolor = plt.cm.coolwarm(NormCp)
        ax.collections[0].set_facecolor(facecolor)
        
        m = cm.ScalarMappable(cmap=cm.coolwarm)
        m.set_array([CpMin, CpMax])
        m.set_clim(vmin=CpMin,vmax=CpMax)
        Cbar = fig.colorbar(m, ax=ax)
        Cbar.set_ticks(np.linspace(CpMin, CpMax, 6))
        Cbar.set_ticklabels(
            [str(round(x,2)) for x in np.linspace(CpMin, CpMax, 6)]
        )
        Cbar.set_label("Cp", rotation=0)
        
        return ax, fig
    
    def displaySurfacePressureCoefficientContour(
        self, elevation=30, azimuth=-60, bodyFixedFrame=False
    ):
        self.plotSurfacePressureCoefficientContour(
            elevation, azimuth, bodyFixedFrame
        )
        
        plt.show()


class PanelMethod(BoundaryElementMethod):
    
    def __init__(self, rigidBody: RigidAerodynamicBody) -> None:
        super().__init__(rigidBody)
        
        self.steadyState = None
        
        if self.wake:
            
            self.Cij = np.pad(self.Cij, ((0, 0), (0, self.wake.numOfFaces)))
    
    @property
    def wake(self):
        return self.rigidBody.wake
        
    def computeWakeInfluence(self):
               
        self.Cij = np.pad(
            self.Cij[:, 0:self.surface.numOfFaces],
            ((0, 0), (0, self.wake.numOfFaces))
        )
        
        index = self.wakeFaceIndex
        
        if not self.rigidBody.isWakeFree:
            
            for i in range(self.surface.numOfFaces):
                
                for j in range(self.wake.numOfWakeRows):
                    
                    for k in range(self.wake.numOfFacesPerWakeRow):
                        
                        self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=k)] = \
                            self.wake.wakeRow[j].panel[k].unitStrength_inducedVelocityPotential(
                                self.surface.panel[i].r_cp
                            )
        
        else:
            
            for i in range(self.surface.numOfFaces):
                
                for j in range(self.wake.numOfWakeRows):
                    
                    for k in range(self.wake.numOfFacesPerWakeRow):
                        
                        self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=k)] = \
                            self.wake.wakeRow[j].panel[k].unitStrength_inducedVelocityPotential(
                                self.rigidBody.ro 
                                + self.surface.panel[i].r_cp.changeBasis(
                                    self.rigidBody.A.T
                                )
                            )
                        
        pass
    
    @property
    def Aij(self):
        
        Aij = np.array(self.Cij[:, :self.surface.numOfFaces])
        index = self.wakeFaceIndex
        
        if self.steadyState :
            
            for i in range(self.surface.numOfFaces):
                
                for j in range(self.wake.numOfWakeRows):
                    
                    upperSheddingFaceIndex, lowerSheddingFaceIndex = \
                        self.rigidBody.sheddingFaces[j]
                    
                    for k in range(self.wake.numOfFacesPerWakeRow):
                        
                        Aij[i][upperSheddingFaceIndex] += \
                            self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=k)]
                            
                        Aij[i][lowerSheddingFaceIndex] -= \
                            self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=k)]
                        
        elif not self.steadyState:
            
            for i in range(self.surface.numOfFaces):
                
                for j in range(self.wake.numOfWakeRows):
                    
                    upperSheddingFaceIndex, lowerSheddingFaceIndex = \
                        self.rigidBody.sheddingFaces[j]
                    
                    Aij[i][upperSheddingFaceIndex] += \
                        self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=0)]
                        
                    Aij[i][upperSheddingFaceIndex] -= \
                        self.Cij[i][index(wakeRowIndex=j, wakefaceIndex=0)]
        
        return Aij                 
    
    @property
    def RHS(self):
        if self.steadyState:
            
            return super().RHS
        
        else:            
            sigma_columnVector = np.array(
                [panel.sigma for panel in self.surfacePanel]
            )
            
            mu_columnVector = np.array(
                [
                    self.wake.wakeRow[i].panel[j].mu 
                    for i in range(self.wake.numOfWakeRows) 
                    for j in range(1, self.wake.numOfFacesPerWakeRow)
                ]
            )
                       
            indexList = [
                    self.wakeFaceIndex(wakeRowIndex=i, wakefaceIndex=j) 
                    for i in range(self.wake.numOfWakeRows) 
                    for j in range(1, self.wake.numOfFacesPerWakeRow)
                ]
              
            RHS = - self.Bij @ sigma_columnVector \
                - self.Cij[:, indexList] @ mu_columnVector
            
            return RHS 
    
    def wakeFaceIndex(self, wakeRowIndex, wakefaceIndex):
        return (
            self.surface.numOfFaces 
            + wakeRowIndex * self.wake.numOfFacesPerWakeRow + wakefaceIndex
        )
    
    def setTrailingEdge(self, trailingEdgeVerticesIDs:list[int]):
        self.rigidBody.setTrailingEdge(
            trailingEdgeVerticesIDs
        )
    
    def locateSheddingFaces(self):
        self.rigidBody.locateSheddingFaces()
    
    def setSheddingFaces(
        self, upperSheddingFacesIDs:list[int], lowerSheddingFacesIDs:list[int]
    ):
        self.rigidBody.setSheddingFaces(
            upperSheddingFacesIDs, lowerSheddingFacesIDs
        )
    
    def setWake(
        self, length:float, numOfWakeFaces:int,
        faceType:str="Quads", isWakeFree:bool=False
    ):
        self.rigidBody.setWake(
            length, numOfWakeFaces, faceType, isWakeFree
        )  
            
    def setVfs(
        self, angleOfAttack: float, sideSlipAngle: float, magnitude: float
    ):
        super().setVfs(angleOfAttack, sideSlipAngle, magnitude)
                
        self.rigidBody.wake.Vinf = self.Vinf
    
    def setVinf(
        self, angleOfAttack: float, sideSlipAngle: float, magnitude: float
    ):
        
        super().setVinf(angleOfAttack, sideSlipAngle, magnitude)
        
        self.rigidBody.wake.Vinf = self.Vinf
        
        e_1 = self.Vinf/self.Vinf.norm()
        if e_1.x != 0:
            e_2 = Vector(0, 1, 0).cross(e_1)
            
        elif e_1.y != 0:
            e_2 = Vector(0, 0, 1).cross(e_1)
            
        elif e_1.z != 0:
            e_2 = Vector(1, 0, 0).cross(e_1)
            
        e_3 = e_2.cross(e_1)
        
        e_x, e_y, e_z = Vector(1, 0, 0), Vector(0, 1, 0), Vector(0, 0, 1)
        self.rigidBody.wake.setWakeLinesRefFrameOrientation(
            np.array(
                [
                    [e_x.dot(e_1), e_y.dot(e_1), e_z.dot(e_1)],
                    [e_x.dot(e_2), e_y.dot(e_2), e_z.dot(e_2)],
                    [e_x.dot(e_3), e_y.dot(e_3), e_z.dot(e_3)]
                ]
            )
        )        
        pass
        
    def inducedVelocity(self, r_p: Vector) -> Vector:
        
        surfaceInducedVelocity = super().inducedVelocity(r_p)      
        
        if self.rigidBody.isWakeFree:
            r_p = self.rigidBody.ro + r_p.changeBasis(self.rigidBody.A.T)
        
        wakeInducedVelocity = Vector(0, 0, 0)
        
        for i in range(self.wake.numOfWakeRows):
                
            for j in range(self.wake.numOfFacesPerWakeRow):
                
                wakeInducedVelocity = wakeInducedVelocity + \
                    self.wake.wakeRow[i].panel[j].inducedVelocity(r_p)
        
        
        if self.rigidBody.isWakeFree:
            wakeInducedVelocity = wakeInducedVelocity.changeBasis(
                self.rigidBody.A
        )
        
        return surfaceInducedVelocity + wakeInducedVelocity
     
    def inducedVelocityPotential(self, r_p: Vector) -> float:
        
        phi_induced =  super().inducedVelocityPotential(r_p)
        
        if self.rigidBody.isWakeFree:
            
            r_p = self.rigidBody.ro + r_p.changeBasis(self.rigidBody.A.T)
            
        for i in range(self.wake.numOfWakeRows):
                
            for j in range(self.wake.numOfFacesPerWakeRow):
                
                phi_induced = phi_induced + \
                    self.wake.wakeRow[i].panel[j].inducedVelocityPotential(r_p)

    def solveLinearSystem(self):
        
        self.computeWakeInfluence()
        
        super().solveLinearSystem()
        
        if self.steadyState:
            
            for wakeRowIndex in range(self.wake.numOfWakeRows):
                
                upperSheddingFaceIndex, lowerSheddingFaceIndex = \
                        self.rigidBody.sheddingFaces[wakeRowIndex]
                
                mu = self.surfacePanel[upperSheddingFaceIndex].mu \
                    - self.surfacePanel[lowerSheddingFaceIndex].mu
                        
                for wakeFaceIndex in range(self.wake.numOfFacesPerWakeRow):
                    
                    self.wake.wakeRow[wakeRowIndex].panel[wakeFaceIndex].mu = mu
        
        else:
                    
            for wakeRowIndex in range(self.wake.numOfWakeRows):
                
                upperSheddingFaceIndex, lowerSheddingFaceIndex = \
                        self.rigidBody.sheddingFaces[wakeRowIndex]
                
                self.wake.wakeRow[wakeRowIndex].panel[0].mu = (
                    self.surfacePanel[upperSheddingFaceIndex].mu
                    - self.surfacePanel[lowerSheddingFaceIndex].mu
                )
                     
        pass       
    
    def solve(self, steadyState=True, iters=0):
        
        if steadyState:
            
            if iters==0:
                
                self.steadyState = steadyState
                
                return super().solve()
            
            else:

                return self.solveIteratively(iters)
            
        else:
            
            V = self.Vfs_at(
                r_rigidBodyPoint=Vector(
                    *self.rigidBody.getVertex(
                        vertex_id=self.rigidBody.trailingEdge[
                            self.wake.numOfWakeLines//2
                        ],
                        bodyFixedFrame=True
                    )
                )
            ).norm()
                        
            length = max([panel.charLength for panel in self.surfacePanel])
            
            dt = length/V 
                    
            self.solveUnsteady(dt, iters)
            
        pass
       
    def solveIteratively(self, iters:int = 10):
        
        self.steadyState = True
        
        if self.rigidBody.isWakeFree:
            iterWake = self.iterFreeWake
        else:
            iterWake = self.iterSufaceFixedWake
        
        self.advanceSolution()
        
        for i in range(iters):
            
            iterWake()
            self.rigidBody.display(bodyFixedFrame=True)
            self.solveLinearSystem()
        
        self.computeSurfaceVelocity()
        self.computeSurfacePressure()
        
        pass
            
    def iterSufaceFixedWake(self):
        
        for i in range(self.rigidBody.wake.numOfWakeLines):
            
            ds = np.array(
                [
                    (
                        self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j+1)
                        - self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                    ).norm()
                    for j in range(
                        self.rigidBody.wake.numOfVerticesPerWakeLine-1
                    )
                ]
            )
                                    
            for j in range(self.rigidBody.wake.numOfVerticesPerWakeLine-1):
                
                Vs = (
                    self.velocity(
                        r_p=self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                    )
                    +
                    self.velocity(
                        r_p=self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j+1)
                    )
                )/2
                
                e_s = Vs/Vs.norm()
                dr = ds[j] * e_s
                
                self.rigidBody.wake.wakeLine[i].changeVertexPosition(
                    vertexIndex=j+1,
                    r=self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)+dr
                )
                
        self.rigidBody.wake.updatePanelsPosition()
        
        pass
    
    def iterFreeWake(self):
        
        for i in range(self.rigidBody.wake.numOfWakeLines):
            
            ds = np.array(
                [
                    (
                        self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j+1)
                        - self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                    ).norm()
                    for j in range(
                        self.rigidBody.wake.numOfVerticesPerWakeLine-1
                    )
                ] 
            )
            
                                
            for j in range(self.rigidBody.wake.numOfVerticesPerWakeLine-1):
                
                Vs = (
                        (
                            self.velocity(
                                r_p=(
                                    self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                                    - self.rigidBody.ro
                                ).changeBasis(self.rigidBody.A)
                            )
                            +
                            self.velocity(
                                r_p=(
                                    self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j+1)
                                    - self.rigidBody.ro
                                ).changeBasis(self.rigidBody.A)
                            )
                        )/2
                ).changeBasis(self.rigidBody.A.T)
                
                e_s = Vs/Vs.norm()
                dr = ds[j] * e_s
                
                self.rigidBody.wake.wakeLine[i].changeVertexPosition(
                    vertexIndex=j+1,
                    r=self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)+dr
                )
                
        self.rigidBody.wake.updatePanelsPosition()
        
        pass
    
    def solveUnsteady(self, dt:float, iters:int):
        
        self.steadyState = False
        
        if self.rigidBody.isWakeFree:
            rollWake = self.rollFreeWake
        else:
            rollWake = self.rollSurfaceFixedWake
        
        for i in range(iters):
            
            self.rigidBody.move_BodyFixedFrame(dt)
            
            self.rigidBody.display()
            
            self.advanceSolution()
            
            rollWake(dt)
            
            self.rigidBody.display()
        
        self.computeSurfaceVelocity()
        self.computeSurfacePressure()
            
        
        pass
    
    def rollSurfaceFixedWake(self, dt:float):
        
        for i in range(self.rigidBody.wake.numOfWakeLines):
            
            for j in range(1, self.rigidBody.wake.numOfVerticesPerWakeLine):
                                
                self.rigidBody.wake.wakeLine[i].moveVertex(
                    vertexIndex=j,
                    dr = dt * self.inducedVelocity(
                        r_p=self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                    )
                )
        
        self.rigidBody.wake.updatePanelsPosition()
        
        pass
    
    def rollFreeWake(self, dt:float):
        
        for i in range(self.rigidBody.wake.numOfWakeLines):
            
            for j in range(1, self.rigidBody.wake.numOfVerticesPerWakeLine):
                
                self.rigidBody.wake.wakeLine[i].moveVertex(
                    vertexIndex=j,
                    dr = dt * self.inducedVelocity(
                        r_p=(
                            self.rigidBody.wake.wakeLine[i].getVertexPositionVector(j)
                            - self.rigidBody.ro
                        ).changeBasis(self.rigidBody.A)
                    ).changeBasis(self.rigidBody.A.T)
                )
        
        self.rigidBody.wake.updatePanelsPosition()
        
        pass
        

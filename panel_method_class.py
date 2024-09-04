import numpy as np
from rigid_body_class import RigidBody
from vector_class import Vector


class BoundaryElementMethod:
    
    def __init__(self, rigidBody:RigidBody) -> None:
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
        for panel_i in self.surfacePanel:
            for panel_j in self.surfacePanel:
                (
                    self.Bij[panel_i.id][panel_j.id],
                    self.Cij[panel_i.id][panel_j.id]
                ) = (
                    panel_j.unitStrength_inducedVelocityPotential(panel_i.r_cp)
                )
        pass
    
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
        for panel in self.surfacePanel:
            panel.mu = mu[panel.id]
        pass
    
    def advanceSolution(self):
        self.computeSourceStrengths()
        self.solveLinearSystem()
    
    def computeSurfaceVelocity(self):
        
        if self.rigidBody.omega.norm() != 0:
            
            for panel in self.surfacePanel:
                panel.setSurfaceVelocity(
                    self.surface.getPanelsAdjacentPanels(panel.id),
                    self.Vfs_at(panel.r_centroid)
                )        
        else:
            
            Vfs = self.Vfs_atBodyFixedFrameOrigin()
            for panel in self.surfacePanel:
                panel.setSurfaceVelocity(
                    self.surface.getPanelsAdjacentPanels(panel.id), Vfs
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
        for panel in self.surface.panel:
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
    
    
if __name__=="__main__":
    from sphere_class import Sphere
    from mesh_class import PanelMesh
    from matplotlib import pyplot as plt
    
    bem = BoundaryElementMethod(
        rigidBody=RigidBody(
            PanelMesh(
                *Sphere(
                    center=(0, 0, 0),
                    radius=1
                ).meshIcoSurface(3)
            ),
            name="Unit Sphere"
        )
    )
    
    bem.setVfs(
        angleOfAttack=0,
        sideSlipAngle=0,
        magnitude=1
    )
    
    # bem.setVinf(
    #     angleOfAttack=0,
    #     sideSlipAngle=0,
    #     magnitude=1
    # )
    
    bem.solve()
    
    ax, fig = bem.rigidBody.plot(bodyFixedFrame=False)
    
    Cp = [panel.Cp for panel in bem.surfacePanel]
    Cp_norm = [(float(Cp_i)-min(Cp))/(max(Cp)-min(Cp)) for Cp_i in Cp]
    facecolor = plt.cm.coolwarm(Cp_norm)
    ax.collections[0].set_facecolor(facecolor)
    
    plt.show()
    
    pass

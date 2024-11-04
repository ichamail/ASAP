import numpy as np
from matplotlib import pyplot as plt
from src.myMath import Vector
from src.utilities import set_axes_equal, is_inside_polygon, LeastSquares
        
class Panel:
    
    collocationPoint_offset:float = 10**(-8)
    
    def __init__(self, vertices:np.ndarray, CCW:bool=True):
                
        self.CCW:bool
        self.numOfVertices:int
        self.r:np.ndarray[np.dtype[Vector]]
        self.r_centroid:Vector # position vector of panel's centroid
        self.e_n:Vector # position vector of panel's normal unit vector
        self.e_l:Vector # position vector of panel's longitudinal unit vector
        self.e_m:Vector  # position vector of panel's transverse unit vector
        self.A = np.identity(3)
        
             
        self.charLength:float # maximum diagonal length or maximum edge length
        self.area:float # area of the panel
        
        self.r_cp:Vector # position vector of panel's control point/collocation point
        self.V:Vector # Velocity Vector at collocation point
        self.Cp:float # Pressure coefficient at collocation point
                
        self.set_vertices(vertices, CCW)
              
    def set_centroid(self):
        self.r_centroid = Vector(0, 0, 0)
        for i in range(self.numOfVertices):
            self.r_centroid = self.r_centroid + self.r[i]
        
        self.r_centroid = self.r_centroid/self.numOfVertices
    
    def set_e_n_and_area(self):
        
        normal = Vector(0, 0, 0)
        for i in range(self.numOfVertices):
            j = (i+1)%self.numOfVertices
            r_i = self.r[i] - self.r_centroid
            r_j = self.r[j] - self.r_centroid
            normal = normal + r_i.cross(r_j)
        
        self.area = normal.norm()/2
        self.e_n = normal/normal.norm()
        if not self.CCW:
            self.e_n = - self.e_n
        
    def set_e_l(self):
        dr = self.r[0] - self.r[1]
        self.e_l = dr/dr.norm()
    
    def set_e_m(self):
        self.e_m = self.e_n.cross(self.e_l)
        
    def set_A(self):
        self.A = np.array(
            [[self.e_l.x, self.e_l.y, self.e_l.z],
             [self.e_m.x, self.e_m.y, self.e_m.z],
             [self.e_n.x, self.e_n.y, self.e_n.z]]
        )
    
    def set_char_length(self):
        
        self.charLength = np.sqrt(self.area)
    
    def set_collocationPoint(self):
        # (r_cp - r_centroid).dot(e_n) --> 0^-
        self.r_cp = self.r_centroid - self.collocationPoint_offset * self.e_n
    
    def set_atributes(self):
        self.set_centroid()
        self.set_e_n_and_area()
        self.set_e_l()
        self.set_e_m()
        self.set_A()
        self.set_char_length()
        self.set_collocationPoint()
         
    def set_vertices(self, vertices:np.ndarray, CCW=True):
        self.CCW = CCW
        self.numOfVertices = len(vertices) # panel's number of vertices
        
        self.r = np.array(
            [Vector(*vertices[i]) for i in range(self.numOfVertices)]
        )
        
        self.set_atributes()

    def plot(self, panelFixedFrame=False):
        
        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        
        
        if panelFixedFrame:
            
            ax.set_title("Panel in Panel-fixed frame of reference $f'_j$")
            
            # Body-fixed frame of reference f'
            ro = -self.r_centroid
            e_x = Vector(1, 0, 0).changeBasis(self.A.T)
            e_y = Vector(0, 1, 0).changeBasis(self.A.T)
            e_z = Vector(0, 0, 1).changeBasis(self.A.T)
            
            ax.quiver(
                ro.x, ro.y, ro.z, e_x.x, e_x.y, e_x.z, label="$e_{x'}$",
                color="slateblue", zorder=3
            )
            ax.quiver(
                ro.x, ro.y, ro.z, e_y.x, e_y.y, e_y.z, label="$e_{y'}$",
                color="darkseagreen", zorder=3
            )
            ax.quiver(
                ro.x, ro.y, ro.z, e_z.x, e_z.y, e_z.z, label="$e_{z'}$",
                color="darkorange", zorder=3
            )
            
                       
            # Panel-fixed frame of reference f'_j
            e_l = self.e_l.changeBasis(self.A)
            e_m = self.e_m.changeBasis(self.A)
            e_n = self.e_n.changeBasis(self.A)
            
            ax.quiver(
                0, 0, 0, e_l.x, e_l.y, e_l.z, label="$e_{l}$",
                color="b", zorder=2
            )
            ax.quiver(
                0, 0, 0, e_m.x, e_m.y, e_m.z, label="$e_{m}$",
                color="g", zorder=2
            )
            ax.quiver(
                0, 0, 0, e_n.x, e_n.y, e_n.z, label="$e_{n}$",
                color="r", zorder=2
            )
            
            
            # Panel's Edges
            r = [
                (self.r[i] - self.r_centroid).changeBasis(self.A)
                for i in range(self.numOfVertices)
            ]
            
            ax.plot3D(
                [r[i%self.numOfVertices].x for i in range(self.numOfVertices+1)],
                [r[i%self.numOfVertices].y for i in range(self.numOfVertices+1)],
                [r[i%self.numOfVertices].z for i in range(self.numOfVertices+1)],
                "k", linewidth=2, zorder=1, label="panel"
            )
            
        else:
            
            ax.set_title("Panel in body-fixed frame of reference $f'$")
            
            # Body-fixed frame of reference f'
            e_x = Vector(1, 0, 0)
            e_y = Vector(0, 1, 0)
            e_z = Vector(0, 0, 1)
            
            ax.quiver(
                0, 0, 0, e_x.x, e_x.y, e_x.z, label="$e_{x'}$",
                color="slateblue", zorder=3
            )
            ax.quiver(
                0, 0, 0, e_y.x, e_y.y, e_y.z, label="$e_{y'}$",
                color="darkseagreen", zorder=3
            )
            ax.quiver(
                0, 0, 0, e_z.x, e_z.y, e_z.z, label="$e_{z'}$",
                color="darkorange", zorder=3
            )
            
                        
            # Panel-fixed frame of reference f'_j
                        
            ax.quiver(
                self.r_centroid.x, self.r_centroid.y, self.r_centroid.z,
                self.e_l.x, self.e_l.y, self.e_l.z, label="$e_{l}$",
                color="b", zorder=2
            )
            ax.quiver(
                self.r_centroid.x, self.r_centroid.y, self.r_centroid.z,
                self.e_m.x, self.e_m.y, self.e_m.z, label="$e_{m}$",
                color="g", zorder=2
            )
            ax.quiver(
                self.r_centroid.x, self.r_centroid.y, self.r_centroid.z,
                self.e_n.x, self.e_n.y, self.e_n.z, label="$e_{n}$",
                color="r", zorder=2
            )
            
            
            # Panel's Edges
            ax.plot3D(
                [
                    self.r[i%self.numOfVertices].x 
                    for i in range(self.numOfVertices+1)
                ],
                [
                    self.r[i%self.numOfVertices].y 
                    for i in range(self.numOfVertices+1)
                ],
                [
                    self.r[i%self.numOfVertices].z 
                    for i in range(self.numOfVertices+1)
                ],
                "k", linewidth=2, zorder=1, label="panel"
            )
        
        set_axes_equal(ax)
        
        ax.set_xlabel('X')
        ax.set_ylabel('Y')

        ax.legend()
            
        return ax, fig
    
    def display(self, panelFixedFrame=False):
        
        ax, fig = self.plot(panelFixedFrame)
        
        plt.show()
        
        pass
    
    def display_point_P(self, r_p:Vector, panelFixedFrame=False):
        
        if panelFixedFrame:
            r_p = (r_p - self.r_centroid).changeBasis(self.A) 
        
        ax, fig = self.plot(panelFixedFrame)
        ax.scatter(r_p.x, r_p.y, r_p.z, c="r", edgecolors="k")
        
        ax.set_xlim3d(
            min([ax.get_xlim()[0], r_p.x]),
            max([ax.get_xlim()[1], r_p.x])
        )
        ax.set_ylim3d(
            min([ax.get_ylim()[0], r_p.y]),
            max([ax.get_ylim()[1], r_p.y])
        )
        ax.set_zlim3d(
            min([ax.get_zlim()[0], r_p.z]),
            max([ax.get_zlim()[1], r_p.z])
        )
        
        set_axes_equal(ax)
        plt.show()
        
        pass

   
class QuadPanel(Panel):
    
    def __init__(self, vertices:np.ndarray, CCW:bool=True):
        
        # VSAero atributes
        self.SMP:float
        self.SMQ:float
        self.T:Vector
        
        super().__init__(vertices, CCW)
       
    def set_e_n_and_area(self):
        
        normal = (self.r[2] - self.r[0]).cross(self.r[3] - self.r[1])
        self.area = normal.norm()/2
        self.e_n = normal/normal.norm()
        if not self.CCW:
            self.e_n = - self.e_n
       
    def set_e_m(self):
               
        dr = (self.r[2] + self.r[3])/2 - self.r_centroid
        self.SMQ = dr.norm()
        self.e_m = dr/self.SMQ
        
    def set_e_l(self):
        self.e_l = self.e_m.cross(self.e_n)
        
    def set_T(self):
        self.T = (self.r[2] + self.r[1])/2 - self.r_centroid
        self.SMP = self.T.norm()
           
    def set_char_length(self):
        self.charLength = np.max(
            [(self.r[3] - self.r[1]).norm(), (self.r[0] - self.r[2]).norm()]
        )
      
    def set_atributes(self):
        self.set_centroid()
        self.set_e_n_and_area()
        self.set_e_m()
        self.set_e_l()
        self.set_A()
        self.set_char_length()
        self.set_T()
        self.set_collocationPoint()

        
class TriPanel(Panel):
    
    def __init__(self, vertices:np.ndarray, CCW:bool=True):
        
        super().__init__(vertices, CCW)
        
    def set_e_n_and_area(self):
        
        normal = (self.r[0] - self.r[1]).cross(self.r[0] - self.r[2])
        self.area = normal.norm()/2
        self.e_n = normal/normal.norm()
        if not self.CCW:
            self.e_n = - self.e_n
    

class Source(Panel):
    
    farFieldFactor:float = 10
    
    def __init__(self, vertices: np.ndarray, CCW: bool = True):
        super().__init__(vertices, CCW)
        
        self.sigma:float = 0
    
    def unitStrength_inducedVelocityPotential(self, r_p:Vector) -> float:
        
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
                
        if r_p.norm() > self.farFieldFactor * self.charLength:
            
            phi = self.area/r_p.norm()
        
        else:
            
            phi = 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                d_ab = r_ab.norm()
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
                                                
                if d_ap + d_bp - d_ab != 0:
                    phi += (
                        (r_ap.x * r_ab.y - r_ap.y * r_ab.x)/d_ab
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )
                          
                
                if r_ab.x != 0 and r_p.z != 0:
                    # paper of Lothar birk  (correct)
                    # Katz & Plotkin formula is wrong
                    
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    
                    phi += - r_p.z * (
                        np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                                - r_ap.x * r_ap.y
                            )
                            /
                            (
                                r_p.z*d_ap
                            )
                        )
                        - np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                                     
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering

        return  - 1/(4 * np.pi) * phi

    def inducedVelocityPotential(self, r_p:Vector) -> float:
        return self.sigma * self.unitStrength_inducedVelocityPotential(r_p)

    def inducedVelocity(self, r_p:Vector) -> Vector:
        # slower but for a rectangular panel of unit strength it produces v_z = 0.5 on panel's vertices and edges
        
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
        if r_p.norm() > self.farFieldFactor * self.charLength:
            v_p = self.area/(r_p.norm()**3) * Vector(r_p.x, r_p.y, r_p.z)
        
        elif r_p.z == 0:
            
            v_x, v_y, v_z = 0, 0, 0
                        
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                d_ab = r_ab.norm()
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
                
                
                if d_ap + d_bp - d_ab != 0:
                    
                    v_x += (
                        r_ab.y/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )
                    
                    v_y += (
                        - r_ab.x/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )

            if is_inside_polygon(
                    points=
                    [(r[i].x, r[i].y) for i in range(self.numOfVertices)],
                    p=(r_p.x, r_p.y)   
                ):
                   
                v_z =  2*np.pi
            
            v_p = Vector(v_x, v_y, v_z)

        else:
            
            v_x, v_y, v_z = 0, 0, 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                d_ab = r_ab.norm()
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
                
                
                if d_ap + d_bp - d_ab != 0:
                    
                    v_x += (
                        r_ab.y/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )
                    
                    v_y += (
                        - r_ab.x/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )

                if r_ab.x != 0:
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    v_z += (
                        np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                                - r_ap.x * r_ap.y
                            )
                            /
                            (
                                r_p.z * d_ap
                            )
                        )
                        - np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                                                    
            if self.CCW: v_z = - v_z # Hess and Smith integrals are calculated with clock wise ordering
            
            v_p = Vector(v_x, v_y, v_z)
            
        return self.sigma/(4 * np.pi) * v_p.changeBasis(self.A.T)
    
    def inducedVelocity(self, r_p:Vector) -> Vector:
        # faster but for a rectangular panel of unit strength it produces v_z = 0 on panel's vertices and edges
        
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
        if r_p.norm() > self.farFieldFactor * self.charLength:
            v_p = self.area/(r_p.norm()**3) * Vector(r_p.x, r_p.y, r_p.z)
        
        else:
            
            v_x, v_y, v_z = 0, 0, 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                d_ab = r_ab.norm()
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
                
                
                if d_ap + d_bp - d_ab != 0:
                    
                    v_x += (
                        r_ab.y/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )
                    
                    v_y += (
                        - r_ab.x/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )

                if r_ab.x != 0 and r_p.z !=0:
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    v_z += (
                        np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                                - r_ap.x * r_ap.y
                            )
                            /
                            (
                                r_p.z * d_ap
                            )
                        )
                        - np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                    
                                             
            if self.CCW: v_z = - v_z # Hess and Smith integrals are calculated with clock wise ordering
            
            v_p = Vector(v_x, v_y, v_z)
            
        return self.sigma/(4 * np.pi) * v_p.changeBasis(self.A.T)
    
    def inducedVelocity(self, r_p:Vector) -> Vector:
        # faster but for a rectangular panel of unit strength it produces v_z = 0.25 on panel edges and v_z = 0.125 on panel vertices
                
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
        if r_p.norm() > self.farFieldFactor * self.charLength:
            v_p = self.area/(r_p.norm()**3) * Vector(r_p.x, r_p.y, r_p.z)
        
        else:
            
            v_x, v_y, v_z = 0, 0, 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                d_ab = r_ab.norm()
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
                
                
                if d_ap + d_bp - d_ab != 0:
                    
                    v_x += (
                        r_ab.y/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )
                    
                    v_y += (
                        - r_ab.x/d_ab 
                        * np.log((d_ap + d_bp + d_ab)/(d_ap + d_bp - d_ab))
                    )

                if r_ab.x != 0:
                    
                    v_z += (
                        arctan(
                            r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                            - r_ap.x * r_ap.y
                            ,
                            r_p.z * d_ap
                        )
                        - arctan(
                            r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                            - r_bp.x * r_bp.y
                            ,
                            r_p.z * d_bp
                        )
                    )
                    
                                             
            if self.CCW: v_z = - v_z # Hess and Smith integrals are calculated with clock wise ordering
            
            v_p = Vector(v_x, v_y, v_z)
            
        return self.sigma/(4 * np.pi) * v_p.changeBasis(self.A.T)

       
class Doublet(Panel):
    
    farFieldFactor:float = 10
    vortexCoreRadius:float = 10**(-6)
    
    def __init__(self, vertices: np.ndarray, CCW: bool = True):
        super().__init__(vertices, CCW)
        
        self.mu:float = 0
        
    def unitStrength_inducedVelocityPotential(self, r_p:Vector) -> float:
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
                
        if r_p.norm() > self.farFieldFactor * self.charLength:
            
            phi = self.area * r_p.z/(r_p.norm()**3)
        
        elif r_p.z == 0:
            
            if is_inside_polygon(
                    points=
                    [(r[i].x, r[i].y) for i in range(self.numOfVertices)],
                    p=(r_p.x, r_p.y)   
                ):
                   
                # phi =  2 * np.pi   # r_p.z --> 0^+
                phi = - 2 * np.pi  # r_p.z --> 0^-
                
            else:
                
                phi = 0

        else:
            
            phi = 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
            
                if r_ab.x != 0:
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    phi += (
                        np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                                - r_ap.x * r_ap.y
                            )
                            /
                            (
                                r_p.z * d_ap
                            )
                        )
                        - np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return  1/(4*np.pi) * phi
    
    def unitStrength_inducedVelocityPotential(self, r_p:Vector) -> float:
        # faster but produces phi = 0 on panel vertices and edges
        
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
                
        if r_p.norm() > self.farFieldFactor * self.charLength:
            
            phi = self.area * r_p.z/(r_p.norm()**3)
        
        elif r_p.norm() == 0:
            
            # phi =  2 * np.pi   # r_p.z --> 0^+
            phi = - 2 * np.pi  # r_p.z --> 0^-

        else:
            
            phi = 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
            
                if r_ab.x != 0 and r_p.z !=0:
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    phi += (
                        np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                                - r_ap.x * r_ap.y
                            )
                            /
                            (
                                r_p.z * d_ap
                            )
                        )
                        - np.arctan(
                            (
                                r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return  1/(4*np.pi) * phi
    
    def unitStrength_inducedVelocityPotential(self, r_p:Vector) -> float:
                
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
                
        if r_p.norm() > self.farFieldFactor * self.charLength:
            
            phi = self.area * r_p.z/(r_p.norm()**3)
        
        else:
            
            phi = 0
            
            for i in range(self.numOfVertices):
                                         
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()
            
                if r_ab.x != 0:
                    # the inverse tangents in equation are evaluated in 
                    # (-π/2 < θ < π/2) so we should use arctan and not arctan2
                    # θ = arctan2, (-π < θ =< π)
                    phi += (
                        arctan(
                            r_ab.y/r_ab.x * (r_ap.x**2 + r_p.z**2)
                            - r_ap.x * r_ap.y
                            ,
                            r_p.z * d_ap
                        )
                        - arctan(
                            r_ab.y/r_ab.x * (r_bp.x**2 + r_p.z**2)
                            - r_bp.x * r_bp.y
                            ,
                            r_p.z * d_bp
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return  1/(4*np.pi) * phi
    
    def inducedVelocityPotential(self, r_p:Vector) -> float:
        return self.mu * self.unitStrength_inducedVelocityPotential(r_p)

    def inducedVelocity(self, r_p:Vector) -> Vector:
        
        r_p = (r_p - self.r_centroid).changeBasis(self.A)
        r = np.array([
            (self.r[i] - self.r_centroid).changeBasis(self.A)
            for i in range(self.numOfVertices)
        ])
        
        if r_p.norm() > self.farFieldFactor * self.charLength:
            
            v_p = - self.area/(r_p.norm()**5) * Vector(
               - 3 * r_p.x * r_p.z,
               - 3 * r_p.y * r_p.z,
               r_p.x**2 + r_p.y**2 - 2 * r_p.z**2
            )
                        
        else:
            v_x, v_y, v_z = 0, 0, 0
            for i in range(self.numOfVertices):
                a, b = i, (i+1)%self.numOfVertices
                r_ab = r[b] - r[a]
                r_ap, r_bp = r_p - r[a], r_p - r[b]
                d_ap, d_bp = r_ap.norm(), r_bp.norm()

                denominator = (
                    d_ap * d_bp
                    * (d_ap*d_bp + r_ap.x * r_bp.x + r_ap.y * r_bp.y + r_p.z**2)
                )
                
                # if denominator != 0:
                epsilon = self.vortexCoreRadius * self.charLength
                if not ( - epsilon < denominator < epsilon ):
                    
                    term = (d_ap + d_bp)/denominator
                    
                    v_x +=  - r_p.z * r_ab.y * term
                    
                    v_y += r_p.z * r_ab.x * term
                    
                    v_z += (r_bp.x * r_ap.y - r_ap.x * r_bp.y) * term
            
            v_p = Vector(v_x, v_y, v_z)
            
            if self.CCW: v_p = - v_p # Hess and Smith integrals are calculated with clock wise ordering        
            
        
        return - self.mu/(4 * np.pi) * v_p.changeBasis(self.A.T)


class SurfacePanel(Source, Doublet):
    
    def unitStrength_inducedVelocityPotential(self, r_p: Vector) -> float:
                
        return (
            Source.unitStrength_inducedVelocityPotential(self, r_p),
            Doublet.unitStrength_inducedVelocityPotential(self, r_p)
        )
        
    def inducedVelocityPotential(self, r_p: Vector) -> float:
        
        phi = self.unitStrength_inducedVelocityPotential(r_p)
                
        return self.sigma * phi[0] + self.mu * phi[1]
    
    def inducedVelocity(self, r_p: Vector) -> Vector:
        return (
            Source.inducedVelocity(self, r_p) 
            + Doublet.inducedVelocity(self, r_p)
        )

    def setSurfaceVelocity(self, adjacentPanels, Vfs):
        
        """
        V = Vx*ex + Vy*ey + Vz*ez or u*i + V*j  + w*k (body-fixed frame of ref)
        
        V = Vl*l + Vm*m + Vn*n (panel's local frame of reference)
        
        sigma = (e_n * nabla)(φ - φ_i) = (e_n * nabla)(φ - φ_infty) =>
        sigma = e_n * (V - V_infty) = e_n * (v + V_infty - V_infty) =>
        sigma = e_n * v = vn
        
        μ = φ - φ_i = φ - φ_infty => nabla μ = nabla (φ - φ_infty) = V - V_infty
        nabla μ = v + V_infty - V_infty => nabla μ = v
        
        (r_ij * nabla)μ = μ_j - μ_i
        
        (r_ij * nabla)μ = Δl_ij*dμ/dl + Δm_ij*dμ/dm + Δn_ij*dμ/dn 
                        =~ Δl_ij*dμ/dl + Δm_ij*dμ/dm
        
                
        (r_ij * nabla)μ = Δl_ij(vl) + Δm_ij(vm) = μ_j - μ_i
        
        [[Δl_i1 , Δm_i1]                [[μ_1 - μ_i]
        [Δl_i2 , Δm_i2]      [[vl]       [μ_2 - μ_i]
        [Δl_i3 , Δm_i3]  =    [vm]]  =   [μ_3 - μ_i]
            ....                            ....
        [Δl_iN , Δm_iN]]                 [μ_4 - μ_i]]
        
        least squares method ---> vl, vm
        """
        
        n = len(adjacentPanels)
        A = np.zeros((n, 2))
        b = np.zeros((n, 1))
        
        for j in range(n):
            panel_j = adjacentPanels[j]
            # Rodrigues' rotation formula
            phi = np.arccos(self.e_n.dot(panel_j.e_n))
            isConvex = False
            if self.e_n.dot(panel_j.r_centroid - self.r_centroid) <= 0:
                isConvex = True
            
            index = []
            for ii in range(self.numOfVertices):
                for jj in range(panel_j.numOfVertices):
                    if self.r[ii] == panel_j.r[jj]:
                        index.append(ii)
                        break
            
            if index[1]-index[0]==1:
                k = self.r[index[1]] - self.r[index[0]]
            else:
                k = self.r[index[0]] - self.r[index[1]]
                
            if not self.CCW:
                k = - k
            
            if isConvex:
                k = -k/k.norm()
            else:
                k = k/k.norm()
                
            v = panel_j.r_centroid - self.r[index[0]]
            # Rodrigues' rotation formula
            v_rot = v * np.cos(phi) + k.cross(v)*np.sin(phi) + k * (k.dot(v)) * (1 - np.cos(phi))
            
            r_centroid_rot = v_rot + self.r[index[0]]
            r_ij = (r_centroid_rot - self.r_centroid).changeBasis(self.A)
            A[j][0], A[j][1], = r_ij.x, r_ij.y
            b[j][0] = panel_j.mu - self.mu
          
        # nabla_mu = LeastSquares(A, b)
        nabla_mu, _, _, _ = np.linalg.lstsq(A, b)
        
        vl, vm, vn = nabla_mu[0][0], nabla_mu[1][0], self.sigma
        
        self.V = Vector(vl, vm, vn).changeBasis(self.A.T) + Vfs
        
        pass

    # def setSurfaceVelocity(self, adjacentPanels, Vfs):
        
    #     """
    #     V = Vx*ex + Vy*ey + Vz*ez or u*i + V*j  + w*k (body-fixed frame of ref)
        
    #     V = Vl*l + Vm*m + Vn*n (panel's local frame of reference)
        
    #     sigma = (e_n * nabla)(φ - φ_i) = (e_n * nabla)(φ - φ_infty) =>
    #     sigma = e_n * (V - V_infty) = e_n * (v + V_infty - V_infty) =>
    #     sigma = e_n * v = vn
        
    #     μ = φ - φ_i = φ - φ_infty => nabla μ = nabla (φ - φ_infty) = V - V_infty
    #     nabla μ = v + V_infty - V_infty => nabla μ = v
        
    #     (r_ij * nabla)μ = μ_j - μ_i
        
    #     (r_ij * nabla)μ = Δl_ij*dμ/dl + Δm_ij*dμ/dm + Δn_ij*dμ/dn 
    #                     =~ Δl_ij*dμ/dl + Δm_ij*dμ/dm
        
                
    #     (r_ij * nabla)μ = Δl_ij(vl) + Δm_ij(vm) = μ_j - μ_i
        
    #     [[Δl_i1 , Δm_i1]                [[μ_1 - μ_i]
    #     [Δl_i2 , Δm_i2]      [[vl]       [μ_2 - μ_i]
    #     [Δl_i3 , Δm_i3]  =    [vm]]  =   [μ_3 - μ_i]
    #         ....                            ....
    #     [Δl_iN , Δm_iN]]                 [μ_4 - μ_i]]
        
    #     least squares method ---> vl, vm
    #     """
        
    #     n = len(adjacentPanels)
    #     A = np.zeros((n, 2))
    #     b = np.zeros((n, 1))
        
    #     for j in range(n):
                        
    #         r_ij = (
    #             adjacentPanels[j].r_centroid - self.r_centroid
    #         ).changeBasis(self.A)
            
    #         A[j][0], A[j][1], = r_ij.x, r_ij.y
            
    #         b[j][0] = adjacentPanels[j].mu - self.mu
        
    #     # nabla_mu = LeastSquares(A, b)
    #     nabla_mu, _, _, _ = np.linalg.lstsq(A, b)
        
    #     vl, vm, vn = nabla_mu[0][0], nabla_mu[1][0], self.sigma
        
    #     self.V = Vector(vl, vm, vn).changeBasis(self.A.T) + Vfs
        
    #     pass


class WakePanel(Doublet):
    pass


class SurfaceTriPanel(TriPanel, SurfacePanel):
    pass


class SurfaceQuadPanel(QuadPanel, SurfacePanel):
    
    def setSurfaceVelocity(self, adjacentPanels, Vfs):
                
        if len(adjacentPanels) == 4:
            
            areAllAdjacentPanelsQuad = True
            for panel in adjacentPanels:
                if not isinstance(panel, SurfaceQuadPanel):
                    areAllAdjacentPanelsQuad = False
                    break
            
            if areAllAdjacentPanelsQuad:
                self.setVSAeroSurfaceVelocity(adjacentPanels, Vfs)
            else:
                return super().setSurfaceVelocity(adjacentPanels, Vfs)
            
        else:
            return super().setSurfaceVelocity(adjacentPanels, Vfs)
    
    def setVSAeroSurfaceVelocity(self, adjacentPanels, Vfs):
        """
        QuadPanels are a special case of Panels than can be used for
        structured surface meshes.
        
        get_surfaceVelocity() use surface numerical differentiation with 
        central, forward and rearward 2nd order finite difference schemes suitable for structured surface meshes.
        
        This function computes the surface velocity, following the notation of NASA Contractor Report 4023 "Program VSAERO theory Document,
        A Computer Program for Calculating Nonlinear Aerodynamic Characteristics
        of Arbitrary Configurations, Brian Maskew"
        
        check pages 48-50 and 23-25
        
        Also
        sigma = (e_n * nabla)(φ - φ_i) = (e_n * nabla)(φ - φ_infty) =>
        sigma = e_n * (V - V_infty) = e_n * (v + V_infty - V_infty) =>
        sigma = e_n * v = vn
        
        μ = φ - φ_i = φ - φ_infty => nabla μ = nabla (φ - φ_infty) = V - V_infty
        nabla μ = v + V_infty - V_infty => nabla μ = v
        """
       
        e_T = self.T/self.T.norm() 
        errorBound = 10**(-3) * self.charLength 
                  
        for j in range(len(adjacentPanels)):
            
            r_kj = adjacentPanels[j].r_centroid - self.r_centroid
        
            if abs(e_T.dot(r_kj)) > abs(self.e_m.dot(r_kj)):
                
                if e_T.dot(r_kj) > 0:
                    
                    if r_kj.norm() <= (
                            self.SMP + adjacentPanels[j].SMP + errorBound
                        ):
                        N2 = adjacentPanels[j]    
                    else:
                        N4 = adjacentPanels[j]
                        
                else:
                    
                    if r_kj.norm() <= (
                        self.SMP + adjacentPanels[j].SMP + errorBound
                    ):
                        N4 = adjacentPanels[j]
                    else:
                        N2 = adjacentPanels[j]   
                                
            else:
                
                if self.e_m.dot(r_kj) > 0:
                    
                    if r_kj.norm() <= (
                        self.SMQ + adjacentPanels[j].SMQ + errorBound
                    ):
                        N3 = adjacentPanels[j]
                    else:
                        N1 = adjacentPanels[j]
                        
                else:
                    
                    if r_kj.norm() <= (
                        self.SMQ + adjacentPanels[j].SMQ + errorBound
                    ):
                        N1 = adjacentPanels[j]
                    else:
                        N3 = adjacentPanels[j]
                
                pass
        
        second_order_central_finite_difference_scheme = False
        second_order_forward_finite_difference_scheme = False
        second_order_backward_finite_difference_scheme = False
        
        if (
            (N1.r_centroid - self.r_centroid).dot(N3.r_centroid - self.r_centroid) < 0
        ):
           second_order_central_finite_difference_scheme = True  
        else:
            if self.e_m.dot(N3.r_centroid - self.r_centroid) > 0:
                second_order_forward_finite_difference_scheme = True
            else:
                second_order_backward_finite_difference_scheme = True
         
        DELQ = 0
        if second_order_central_finite_difference_scheme:
            
            SMQ_k, SMQ_n1, SMQ_n3 = self.SMQ, N1.SMQ, N3.SMQ
            SA = - (SMQ_k + SMQ_n1)
            SB = SMQ_k + SMQ_n3
            DA = (N1.mu - self.mu)/SA
            DB = (N3.mu - self.mu)/SB
            DELQ = (DA * SB - DB * SA)/(SB - SA)
            
            # panel_j_minus1 = N1
            # panel_j_plus1 = N3
            
            # x1 = 0
            # x0 = x1 - self.SMQ - panel_j_minus1.SMQ
            # x2 = x1 + self.SMQ + panel_j_plus1.SMQ
            # mu0 = panel_j_minus1.mu
            # mu1 = self.mu
            # mu2 = panel_j_plus1.mu
            
            # DELQ = mu0 * (x1 - x2)/(x0 - x1)/(x0 - x2) \
            #         + mu1 * (2*x1 - x0 - x2)/(x1 - x0)/(x1 - x2) \
            #         + mu2 * (x1 - x0)/(x2 - x0)/(x2 - x1)
            
        elif second_order_forward_finite_difference_scheme:
            
            panel_j_plus1 = N3
            panel_j_plus2 = N1
            
            x0 = 0
            x1 = x0 + self.SMQ + panel_j_plus1.SMQ
            x2 = x1 + panel_j_plus1.SMQ + panel_j_plus2.SMQ
            
            mu0 = self.mu
            mu1 = panel_j_plus1.mu
            mu2 = panel_j_plus2.mu
            
            DELQ = mu0 * (2*x0 - x1 - x2)/(x0 - x1)/(x0 - x2) \
                    + mu1 * (x0 - x2)/(x1 - x0)/(x1 - x2) \
                    + mu2 * (x0 - x1)/(x2 - x0)/(x2 - x1)
        
        elif second_order_backward_finite_difference_scheme:
            
            panel_j_minus1 = N1
            panel_j_minus2 = N3
            
            x2 = 0
            x1 = x2 - self.SMQ - panel_j_minus1.SMQ
            x0 = x1  - panel_j_minus1.SMQ - panel_j_minus2.SMQ
            
            mu0 = panel_j_minus2.mu
            mu1 = panel_j_minus1.mu
            mu2 = self.mu
            
            DELQ = mu0 * (x2 - x1)/(x0 - x1)/(x0 - x2) \
                    + mu1 * (x2 - x0)/(x1 - x0)/(x1 - x2) \
                    + mu2 * (2*x2 - x0 - x1)/(x2 - x0)/(x2 - x1)
        
        
        second_order_central_finite_difference_scheme = False
        second_order_forward_finite_difference_scheme = False
        second_order_backward_finite_difference_scheme = False    
        
        if (
            (N4.r_centroid - self.r_centroid).dot(N2.r_centroid - self.r_centroid) < 0
        ):
           second_order_central_finite_difference_scheme = True
        else:
            if self.e_l.dot(N2.r_centroid - self.r_centroid) > 0:
                second_order_forward_finite_difference_scheme = True
            else:
                second_order_backward_finite_difference_scheme = True
        
                
        DELP = 0
        if second_order_central_finite_difference_scheme:
            
            SMP_k, SMP_n2, SMP_n4 = self.SMP, N2.SMP, N4.SMP
            SA = - (SMP_k + SMP_n4)
            SB = SMP_k + SMP_n2
            DA = (N4.mu - self.mu)/SA
            DB = (N2.mu - self.mu)/SB
            DELP = (DA * SB - DB * SA)/(SB - SA)
            
            # panel_i_minus1 = N4
            # panel_i_plus1 = N2
            
            # x1 = 0
            # x0 = x1 - self.SMP - panel_i_minus1.SMP
            # x2 = x1 + self.SMP + panel_i_plus1.SMP
            
            # mu0 = panel_i_minus1.mu
            # mu1 = self.mu
            # mu2 = panel_i_plus1.mu
            
            # DELP = mu0 * (x1 - x2)/(x0 - x1)/(x0 - x2) \
            #         + mu1 * (2*x1 - x0 - x2)/(x1 - x0)/(x1 - x2) \
            #         + mu2 * (x1 - x0)/(x2 - x0)/(x2 - x1)
                
        elif second_order_forward_finite_difference_scheme:
            
            panel_i_plus1 = N2
            panel_i_plus2 = N4
            
            x0 = 0
            x1 = x0 + self.SMP + panel_i_plus1.SMP
            x2 = x1 + panel_i_plus1.SMP + panel_i_plus2.SMP
            
            mu0 = self.mu
            mu1 = panel_i_plus1.mu
            mu2 = panel_i_plus2.mu
            
            DELP = mu0 * (2*x0 - x1 - x2)/(x0 - x1)/(x0 - x2) \
                    + mu1 * (x0 - x2)/(x1 - x0)/(x1 - x2) \
                    + mu2 * (x0 - x1)/(x2 - x0)/(x2 - x1)
            
        elif second_order_backward_finite_difference_scheme:
            
            panel_i_minus1 = N4
            panel_i_minus2 = N2
            
            x2 = 0
            x1 = x2 - self.SMP - panel_i_minus1.SMP
            x0 = x1 - panel_i_minus1.SMP - panel_i_minus2.SMP
            
            mu0 = panel_i_minus2.mu
            mu1 = panel_i_minus1.mu
            mu2 = self.mu
            
            DELP = mu0 * (x2 - x1)/(x0 - x1)/(x0 - x2) \
                    + mu1 * (x2 - x0)/(x1 - x0)/(x1 - x2) \
                    + mu2 * (2*x2 - x0 - x1)/(x2 - x0)/(x2 - x1)

        
        T = self.T
        T = T.changeBasis(self.A)
        TM = T.y
        TL = T.x
        
        VL = (self.SMP * DELP - TM * DELQ)/TL 
        VM = DELQ
        VN = self.sigma
                
        self.V = Vector(VL, VM, VN).changeBasis(self.A.T) + Vfs
        
    pass


class WakeTriPanel(TriPanel, WakePanel):
    pass


class WakeQuadPanel(QuadPanel, WakePanel):
    pass


def arctan(y:float, x:float) -> float:
    
    """
    θ = arctan(y/x), -π/2 <= θ <= π/2
    """
    
    arctan = np.arctan2(y,x)
    if x<0:
        if y<0:
            arctan += np.pi 
        else:
            arctan -= np.pi
        
            
    return arctan

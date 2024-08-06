import numpy as np
from vector_class import Vector
from typing import Any
from plot_functions import set_axes_equal
from matplotlib import pyplot as plt
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches
from is_inside_polygon import is_inside_polygon
        
class Panel:
    
    collocationPoint_offset:float = 10**(-15)
    
    def __init__(self, vertices:np.ndarray, CCW:bool=True):
        
        self.id:int = -1  # panel's identity
        
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
                                r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
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
                                r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
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
                                r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
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
                            r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
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
                                r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return - 1/(4*np.pi) * phi
    
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
                                r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
                                - r_bp.x * r_bp.y
                            )
                            /
                            (
                                r_p.z * d_bp
                            )
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return - 1/(4*np.pi) * phi
    
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
                            r_ab.y/r_ab.x * (r_bp.x**2 - r_p.z**2)
                            - r_bp.x * r_bp.y
                            ,
                            r_p.z * d_bp
                        )
                    )
                    
            if self.CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
            
        
        return - 1/(4*np.pi) * phi
    
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
                
                if denominator != 0:
                    term = (d_ap + d_bp)/denominator
                    
                    v_x +=  - r_p.z * r_ab.y * term
                    
                    v_y += r_p.z * r_ab.x * term
                    
                    v_z += (r_bp.x * r_ap.y - r_ap.x * r_bp.y) * term
            
            v_p = Vector(v_x, v_y, v_z)
            
            if self.CCW: v_p = - v_p # Hess and Smith integrals are calculated with clock wise ordering        
            
        
        return self.mu/(4 * np.pi) * v_p.changeBasis(self.A.T)


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

    def get_surfaceVelocity(self, adjacentPanels):
        pass


class WakePanel(Doublet):
    pass


class SurfaceTriPanel(TriPanel, SurfacePanel):
    pass


class SurfaceQuadPanel(QuadPanel, SurfacePanel):
    
    def get_surfaceVelocity(self, adjacentPanels):
        """
        QuadPanels are a special case of Panels than can be used for
        structured surface meshes.
        
        get_surfaceVelocity() use surface numerical differentiation with 
        central, forward and rearward 2nd order finite difference schemes suitable for structured surface meshes.
        """
    
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

def testHexagonPanel():
       
    # Polygon panel
    def hexagon(center=(2, 2, 2), r=1):
        x0, y0, z0 = center
        numS = 6 # number of sides
        # theta0 = (360/(numB-1))/2
        theta0 = 0
        theta = np.linspace(0, 360, numS+1)
        theta = theta + theta0
        theta = theta*(np.pi/180)
        x = x0 + r* np.cos(theta)
        y = y0 + r* np.sin(theta)
        return np.array([[x[i], y[i], z0] for i in range(numS)])
    
    panel = Panel(hexagon())
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(4, 4, 4), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(4, 4, 4), panelFixedFrame=True
    )
    
def testTriPanel():
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2]
    ])
    
    panel = TriPanel(vertices, CCW=True)
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=True
    )
    
def testQuadPanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = QuadPanel(vertices, CCW=True)
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=True
    )
    
    pass

def testSourcePanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = Source(vertices, CCW=True)
    panel.sigma = 1
    
    
    # Velocity Potential - Scalar field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        100
    )
        
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        100
    )
    
    X, Z = np.meshgrid(x, z, indexing="ij")
    phi = np.zeros_like(X)
    
    nx, nz = X.shape
    for i in range(nx):
        for j in range(nz):
            phi[i][j] = panel.unitStrength_inducedVelocityPotential(
                r_p=Vector(X[i][j], panel.r_centroid.y, Z[i][j])
            )
    CS = plt.contour(X, Z, phi)
    plt.clabel(CS, inline = 1, fmt ='% 1.2f', fontsize = 8)
    plt.show()
    
    
    # # Velocity Potential - Scalar field and Velocity Vector Field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        10
    )
    y = np.linspace(
        panel.r_centroid.y - 1 * panel.charLength,
        panel.r_centroid.y + 1 * panel.charLength,
        10
    )
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        10
    )
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    phi = np.zeros_like(X)
    v_x, v_y, v_z = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    nx, ny, nz = X.shape
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz): 
                r_p = Vector(X[i][j][k], Y[i][j][k], Z[i][j][k])
                phi[i][j][k] = panel.inducedVelocityPotential(r_p) 
                v = panel.inducedVelocity(r_p)
                v_x[i][j][k], v_y[i][j][k], v_z[i][j][k] = v.x, v.y, v.z
    
    ax, fig = panel.plot(panelFixedFrame=False)
    scatter = ax.scatter(X, Y, Z, c=phi, cmap="jet")
    fig.colorbar(scatter, ax=ax)
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    
    
    ax, fig = panel.plot(panelFixedFrame=False)
    ax.quiver(X, Y, Z, v_x, v_y, v_z, normalize=True, color='m')
    
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    pass

def testDoubletPanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = Doublet(vertices, CCW=True)
    panel.mu = 1
    
    # Velocity Potential - Scalar field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        100
    )
        
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        100
    )
    
    X, Z = np.meshgrid(x, z, indexing="ij")
    phi = np.zeros_like(X)
    
    nx, nz = X.shape
    for i in range(nx):
        for j in range(nz):
            phi[i][j] = panel.unitStrength_inducedVelocityPotential(
                r_p=Vector(X[i][j], panel.r_centroid.y, Z[i][j])
            )
    CS = plt.contour(X, Z, phi)
    plt.clabel(CS, inline = 1, fmt ='% 1.2f', fontsize = 8)
    plt.show()
    
    
    # # Velocity Potential - Scalar field and Velocity Vector Field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        10
    )
    y = np.linspace(
        panel.r_centroid.y - 1 * panel.charLength,
        panel.r_centroid.y + 1 * panel.charLength,
        10
    )
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        10
    )
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    phi = np.zeros_like(X)
    v_x, v_y, v_z = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    nx, ny, nz = X.shape
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz): 
                r_p = Vector(X[i][j][k], Y[i][j][k], Z[i][j][k])
                phi[i][j][k] = panel.inducedVelocityPotential(r_p) 
                v = panel.inducedVelocity(r_p)
                v_x[i][j][k], v_y[i][j][k], v_z[i][j][k] = v.x, v.y, v.z
    
    ax, fig = panel.plot(panelFixedFrame=False)
    scatter = ax.scatter(X, Y, Z, c=phi, cmap="jet")
    fig.colorbar(scatter, ax=ax)
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    
    
    ax, fig = panel.plot(panelFixedFrame=False)
    ax.quiver(X, Y, Z, v_x, v_y, v_z, normalize=True, color='m')
    
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()

def testSingularity(panel:Source|Doublet, r_p:Vector):
    
    panel.display_point_P(r_p, panelFixedFrame=False)   
    phi = panel.unitStrength_inducedVelocityPotential(r_p)
    V = panel.inducedVelocity(r_p)
    
    print("phi = ", phi)
    print("V = ", V)
 
def testSurfacePanel(r_p):
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = SurfacePanel(vertices, CCW=True)
    
    panel.sigma, panel.mu = 1, 1
    
    panel.display_point_P(r_p)
    
    phi_src, phi_dblt = panel.unitStrength_inducedVelocityPotential(r_p)
    v = panel.inducedVelocity(r_p)
    
    print("phi_src = ", phi_src)
    print("phi_dblt = ", phi_dblt)
    print("v_src + v_dblt = ", v)

   
    
if __name__=='__main__':
    
    # testHexagonPanel()
    
    # testTriPanel()
    
    # testQuadPanel()
    
    testSourcePanel()
    
    testDoubletPanel()
    
    # testSurfacePanel(r_p=Vector(2, 2, 2))
    
from numba import njit, typed
import numpy as np
from src.myMath import Vector


def giveArguments_for_unitStrength_inducedVelocityPotential(r_p, panel):
    
    return r_p, panel.r_centroid, typed.List(panel.r), panel.numOfVertices, panel.A, panel.charLength, panel.area, panel.farFieldFactor, panel.CCW

def giveArguments_for_src_inducedVelocity(r_p, panel):
    return r_p, panel.r_centroid, typed.List(panel.r), panel.numOfVertices, panel.A, panel.charLength, panel.area, panel.farFieldFactor, panel.CCW, panel.sigma

def giveArguments_for_dblt_inducedVelocity(r_p, panel):
    return r_p, panel.r_centroid, typed.List(panel.r), panel.numOfVertices, panel.A, panel.charLength, panel.area, panel.farFieldFactor, panel.CCW, panel.mu, panel.vortexCoreRadius


@njit
def src_unitStrength_inducedVelocityPotential(
    r_p, r_centroid, r, numOfVertices, A, charLength, area, farFieldFactor, CCW
):
    
    r = typed.List(r)
    r_p = (r_p - r_centroid).changeBasis(A)
    
    for i in range(numOfVertices):
        r[i] = (r[i] - r_centroid).changeBasis(A)
    
            
    if r_p.norm() > farFieldFactor * charLength:
        
        phi = area/r_p.norm()
    
    else:
        
        phi = 0
        
        for i in range(numOfVertices):
                                        
            a, b = i, (i+1)%numOfVertices
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
                                    
        if CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering

    return  - 1/(4 * np.pi) * phi
    

@njit
def dblt_unitStrength_inducedVelocityPotential(
    r_p, r_centroid, r, numOfVertices, A, charLength, area, farFieldFactor, CCW
):
    r_p = (r_p - r_centroid).changeBasis(A)
    
    for i in range(numOfVertices):
        r[i] = (r[i] - r_centroid).changeBasis(A)
    
            
    if r_p.norm() > farFieldFactor * charLength:
        
        phi = area * r_p.z/(r_p.norm()**3)
    
    else:
        
        phi = 0
        
        for i in range(numOfVertices):
                                        
            a, b = i, (i+1)%numOfVertices
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
                
        if CCW: phi = - phi # Hess and Smith integrals are calculated with clock wise ordering
        
    
    return  1/(4*np.pi) * phi


@njit
def src_inducedVelocity(r_p, r_centroid, r, numOfVertices, A, charLength, area, farFieldFactor, CCW, sigma):
        
    r_p = (r_p - r_centroid).changeBasis(A)
    
    for i in range(numOfVertices):
        r[i] = (r[i] - r_centroid).changeBasis(A)
    
    if r_p.norm() > farFieldFactor * charLength:
        
        v_p =  Vector(r_p.x, r_p.y, r_p.z) * area/(r_p.norm()**3)
    
    else:
        
        v_x, v_y, v_z = 0, 0, 0
        
        for i in range(numOfVertices):
                                        
            a, b = i, (i+1)%numOfVertices
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
                
                                            
        if CCW: v_z = - v_z # Hess and Smith integrals are calculated with clock wise ordering
        
        v_p = Vector(v_x, v_y, v_z)
        
    return  v_p.changeBasis(A.T) * sigma/(4 * np.pi)


@njit
def dblt_inducedVelocity(r_p, r_centroid, r, numOfVertices, A, charLength, area, farFieldFactor, CCW, mu, vortexCoreRadius):
    r_p = (r_p - r_centroid).changeBasis(A)
    
    for i in range(numOfVertices):
        r[i] = (r[i] - r_centroid).changeBasis(A)
    
    if r_p.norm() > farFieldFactor * charLength:
        
        v_p =  Vector(
            - 3 * r_p.x * r_p.z,
            - 3 * r_p.y * r_p.z,
            r_p.x**2 + r_p.y**2 - 2 * r_p.z**2
        ) * (- area/(r_p.norm()**5))
                    
    else:
        v_x, v_y, v_z = 0, 0, 0
        for i in range(numOfVertices):
            a, b = i, (i+1)%numOfVertices
            r_ab = r[b] - r[a]
            r_ap, r_bp = r_p - r[a], r_p - r[b]
            d_ap, d_bp = r_ap.norm(), r_bp.norm()

            denominator = (
                d_ap * d_bp
                * (d_ap*d_bp + r_ap.x * r_bp.x + r_ap.y * r_bp.y + r_p.z**2)
            )
            
            # if denominator != 0:
            epsilon = vortexCoreRadius * charLength
            if not ( - epsilon < denominator < epsilon ):
                
                term = (d_ap + d_bp)/denominator
                
                v_x +=  - r_p.z * r_ab.y * term
                
                v_y += r_p.z * r_ab.x * term
                
                v_z += (r_bp.x * r_ap.y - r_ap.x * r_bp.y) * term
        
        v_p = Vector(v_x, v_y, v_z)
        
        if CCW: v_p = - v_p # Hess and Smith integrals are calculated with clock wise ordering        
        
    
    return  v_p.changeBasis(A.T) * (- mu/(4 * np.pi))


@njit
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

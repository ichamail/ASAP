import numpy as np
from numba import njit, prange, typed
from .numbaPanelFunctions import src_dblt_unitStrength_inducedVelocityPotential


def jit_computeSurfaceInfluence(bemObj):
        
    return computeSurfaceInfluence(
        numOfFaces=bemObj.surface.numOfFaces,
        Bij=bemObj.Bij,
        Cij=bemObj.Cij,
        r=typed.List(
            [
                typed.List(bemObj.surface.panel[i].r)
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        r_centroid=typed.List(
            [
                bemObj.surface.panel[i].r_centroid 
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        r_cp=typed.List(
            [
                bemObj.surface.panel[i].r_cp
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        numOfVertices=np.array(
            [
                bemObj.surface.panel[i].numOfVertices
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        A=np.array(
            [
                bemObj.surface.panel[i].A
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        charLength=np.array(
            [
                bemObj.surface.panel[i].charLength
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        area=np.array(
            [
                bemObj.surface.panel[i].area
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        CCW=np.array(
            [
                bemObj.surface.panel[i].CCW
                for i in range(bemObj.surface.numOfFaces)
            ]
        ),
        farFieldFactor = bemObj.surface.panel[0].farFieldFactor
    )
      
@njit(parallel=True)
def computeSurfaceInfluence(
    numOfFaces, Bij, Cij, r, r_centroid, r_cp, numOfVertices, A,
    charLength, area, CCW, farFieldFactor
):
     
    for i in prange(numOfFaces):

        for j in prange(numOfFaces):
                         
            Bij[i][j],Cij[i][j]=src_dblt_unitStrength_inducedVelocityPotential(
                r_p=r_cp[i],
                r_centroid=r_centroid[j],
                r=typed.List(r[j]), # r should be passed as a copy
                numOfVertices=numOfVertices[j],
                A=A[j],
                charLength=charLength[j],
                area=area[j],
                CCW=CCW[j],
                farFieldFactor=farFieldFactor   
            )                
    pass



##################### TO DO #######################################

def jit_computeWakeInfluence():
    pass

@njit(parallel=True)
def computeWakeInfluence():
    pass


def jit_surfaceInducedVelocity():
    pass

@njit(parallel=True)
def surfaceInducedVelocity():
    pass


def jit_wakeinducedVelocity():
    pass

@njit(parallel=True)
def wakeinducedVelocity():
    pass


def jit_iterSufaceFixedWake():
    pass

@njit(parallel=True)
def iterSufaceFixedWake():
    pass


def jit_iterFreeWake():
    pass

@njit(parallel=True)
def iterFreeWake():
    pass


def jit_rollSurfaceFixedWake():
    pass

@njit(parallel=True)
def rollSurfaceFixedWake():
    pass


def jit_rollFreeWake():
    pass

@njit(parallel=True)
def rollFreeWake():
    pass

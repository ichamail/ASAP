from matplotlib import pyplot as plt
from src.PanelMethod import Mesh
from src.geometry import Cube, Sphere

def plotCube():
    Mesh(
        *Cube(2, center=(0, 0, 0)).meshSurface(numOfVerticesPerEdge=5)
    ).display()
    
def plotSphere():
    
    Mesh(
        *Sphere(
            center=(0, 0, 0), radius=1
        ).meshUVSurface(numOfMeridians=20, numOfParallels=18)
    ).display(elevation=30, azimuth=-60)
    
    Mesh(
        *Sphere(
            center=(0, 0, 0), radius=1
        ).meshIcoSurface(nu=5)
    ).display()
    
    Mesh(
        *Sphere(
            center=(0, 0, 0), radius=1
        ).meshSpherifiedCube(numOfVerticesPerEdge=10)
    ).display()

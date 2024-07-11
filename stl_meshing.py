import stl
import numpy as np
from vector_class import Vector
from airfoil_class import Airfoil
from wing_class import Wing
from mesh_class import PanelMesh, PanelAeroMesh
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from Algorithms import light_vector
from plot_functions import set_axes_equal


def plot_mesh(mesh, shells2highlight_id_list=[], nodes2highlight_id_list=[],
              elevation=-150, azimuth=-120):
    
    shells = []
    vert_coords = []
    
    for panel in mesh.panels:
        shell=[]
        for r_vertex in panel.r_vertex:
            shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
            vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
        shells.append(shell)
    
    
    light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
    face_normals = [panel.n for panel in mesh.panels]
    dot_prods = [-light_vec * face_normal for face_normal in face_normals]
    min = np.min(dot_prods)
    max = np.max(dot_prods)
    target_min = 0.2 # darker gray
    target_max = 0.6 # lighter gray
    shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                + target_min
                for dot_prod in dot_prods]
    facecolor = plt.cm.gray(shading)
    
    # change shell color to red
    for id in shells2highlight_id_list:
        facecolor[id] = [1, 0, 0, 1]
    
    ax = plt.axes(projection='3d')
    poly3 = Poly3DCollection(shells, facecolor=facecolor)
    ax.add_collection(poly3)
    ax.view_init(elevation, azimuth)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    x = [mesh.nodes[id][0] for id in nodes2highlight_id_list]
    y = [mesh.nodes[id][1] for id in nodes2highlight_id_list]
    z = [mesh.nodes[id][2] for id in nodes2highlight_id_list]
    ax.scatter(x, y, z, c="r")

    vert_coords = np.array(vert_coords)
    x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
    ax.set_xlim3d(x.min(), x.max())
    ax.set_ylim3d(y.min(), y.max())
    ax.set_zlim3d(z.min(), z.max())
    set_axes_equal(ax)    
    plt.show()
    

# create wing object   
root_airfoil = Airfoil(name="naca0012_new", chordLength=1)
tip_airfoil = Airfoil(name="naca0012_new", chordLength=1)
wing = Wing(root_airfoil, tip_airfoil, halfSpan=1, sweepAngle=0, dihedralAngle=0)

# generate wing mesh
num_x_bodyShells = 20
num_y_Shells = 25

nodes, shells = wing.meshSurface(
    numOfChordWiseFaces=num_x_bodyShells, numOfSpanWiseFaces=num_y_Shells,
    spanWiseSpacing="uniform", faceType="triangular",
    mesh_MainSurface=True, mesh_WingTips=True, mesh_wake=False 
)

wing_mesh = PanelMesh(nodes, shells)

wing_mesh.plot_mesh_inertial_frame(elevation=-150, azimuth=-120)

facet_list = [
    stl.Facet(
        normal=(panel.n.x, panel.n.y, panel.n.z),
        vertices=(
            (panel.r_vertex[0].x, panel.r_vertex[0].y, panel.r_vertex[0].z),
            (panel.r_vertex[1].x, panel.r_vertex[1].y, panel.r_vertex[1].z),
            (panel.r_vertex[2].x, panel.r_vertex[2].y, panel.r_vertex[2].z)
        )
    )
    for panel in wing_mesh.panels
]
wing_solid = stl.Solid(name="wing", facets=facet_list)

stl_mesh = open("coord_seligFmt/wing.stl", "w")

wing_solid.write_ascii(stl_mesh)





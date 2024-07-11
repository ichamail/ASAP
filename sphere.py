from matplotlib import pyplot as plt
import numpy as np


def sphere(radius, num_longitude, num_latitude,
                     mesh_shell_type='triangular'):
    
    # sqrt = np.sqrt(node_num)
    # num_longitude, num_latitude = np.floor(sqrt), np.ceil(sqrt)
    # num_longitude, num_latitude = np.ceil(sqrt), np.floor(sqrt)
    
    num_theta = int(num_longitude)
    num_phi = int(num_latitude)
    
    thetas = np.linspace(0, 2*np.pi, num_theta)
    phis = np.linspace(0, np.pi, num_phi)
    r = radius
    
    thetas = np.delete(thetas, -1) # 0 = 2π
    num_theta = num_theta-1  # 0 = 2π
    
    nodes = []
    shells = []
    
    for theta in thetas:
        for phi in phis:
            x = r * np.sin(phi) * np.cos(theta)
            y = r * np.sin(phi) * np.sin(theta) 
            z = r * np.cos(phi)
            node = (x, y, z)
            nodes.append(node)
        
    if mesh_shell_type == 'quadrilateral':
        # quadrilateral shells
        
        for i in range(num_theta):
            for j in range(num_phi-1):
                
                if j == 0:
                    # shells.append([i*num_phi + j,
                    #             i*num_phi + j+1,
                    #             ((i+1)%num_theta )* num_phi + j+1])
                    
                    shells.append([j,
                                i*num_phi + j+1,
                                ((i+1)%num_theta )* num_phi + j+1])
                    
                elif j == num_phi-2:            
                    # shells.append([i*num_phi + j,
                    #                i*num_phi + j+1,
                    #             ((i+1)%num_theta )* num_phi + j])
                    
                    shells.append([i*num_phi + j,
                                   j+1,
                                ((i+1)%num_theta )* num_phi + j])
                    
                    # alternatively
                    # shells.append([i*num_phi + j,
                    #             ((i+1)%num_theta )* num_phi + j+1,
                    #             ((i+1)%num_theta )* num_phi + j])
                    
                    
                else:            
                    shells.append([i*num_phi + j,
                                i*num_phi + j+1,
                                ((i+1)%num_theta )* num_phi + j+1,
                                ((i+1)%num_theta )* num_phi + j])
                
    elif mesh_shell_type == 'triangular':
        ## triangular panels
        for i in range(num_theta):
            for j in range(num_phi-1):
                
                if j == 0:            
                    # shells.append([i*num_phi + j,
                    #                i*num_phi + j+1,
                    #                ((i+1)%num_theta )* num_phi + j+1])
                    
                    shells.append([j,
                                   i*num_phi + j+1,
                                   ((i+1)%num_theta )* num_phi + j+1])
                  
                elif j == num_phi-2:
                    
                    # shells.append([i*num_phi + j,
                    #                i*num_phi + j+1,
                    #             ((i+1)%num_theta )* num_phi + j])
                    
                    shells.append([i*num_phi + j,
                                   j+1,
                                   ((i+1)%num_theta )* num_phi + j])
                    
                    # alternatively            
                    # shells.append([i*num_phi + j,
                    #                ((i+1)%num_theta)*num_phi + j+1,
                    #                ((i+1)%num_theta)*num_phi + j])
                    
                else:
                    shells.append([i*num_phi + j,
                                   i*num_phi + j+1,
                                   ((i+1)%num_theta)*num_phi + j+1])
                    
                    shells.append([((i+1)%num_theta)*num_phi + j+1,
                                   ((i+1)%num_theta)*num_phi + j,
                                   i*num_phi + j])
    
    return nodes, shells

if __name__=='__main__':
    from mesh_class import Mesh
    
    radius = 1
    num_longitude, num_latitude = 5, 6
    nodes, shells = sphere(radius, num_longitude, num_latitude,
                                     mesh_shell_type='triangular')
    
    neighbors = []                    
    for i, shell_i in enumerate(shells):
        neighbors.append([])
        for j, shell_j in enumerate(shells):
            if i != j and Mesh.do_intersect(shell_i, shell_j):
                neighbors[-1].append(j)
    
    print("node id, (x, y, x)")
    for i, node in enumerate(nodes):
        print(i, node)
        
    print("shell id, [nodes' ids]")
    for i, shell in enumerate(shells):
        print(i, shell)
    
    print("shell id, node coords")
    for shell_id, shell in enumerate(shells):
        for node_id in shell:
            print(shell_id, ": ",
                  nodes[node_id][0],nodes[node_id][1],nodes[node_id][2])        
    
    print("shell id, [neighbors' ids]")
    for shell_id, neighbor in enumerate(neighbors):
        print(shell_id, neighbor)
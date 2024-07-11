from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import axes3d
from Algorithms import light_vector
from plot_functions import set_axes_equal, move_view
import numpy as np
from vector_class import Vector
from panel_class import Panel, TriPanel, QuadPanel
from copy import deepcopy
import stl




class Mesh:
    
    def __init__(self, nodes:list, shells:list,
                 ):
        self.nodes = nodes
        self.shells = shells
        self.node_num = len(nodes)
        self.shell_num = len(shells)

        self.shell_neighbours = self.locate_shells_adjacency()
        
        
        
        ### unsteady features ###
                
        # position vector of body-fixed frame origin
        self.ro = Vector((0, 0, 0)) 
        
        # velocity vector of body-fixed frame origin
        self.Vo = Vector((0, 0, 0))
        
        # orientation of body-fixed frame : (yaw-angle, pitch-angle, roll-angle)
             
        self.theta = Vector((0, 0, 0))
        
        # Rotation Matrix
        self.R = np.array([[1, 0, 0],
                           [0, 1, 0],
                           [0, 0, 1]])
                
        # body-fixed frame's angular velocity vector
        self.omega = Vector((0, 0, 0))

    @classmethod
    def generate_from_stl_file(cls, fileName:str, filePath="STL_models/"):
        
        fileName = filePath + fileName + ".stl"
        
        solid_stl_file = open(fileName, "r")
        solid_stl_obj = stl.read_ascii_file(solid_stl_file)

        nodes = []
        shells = []

        for facet in solid_stl_obj.facets:
            for vertex in facet.vertices:
                node = (vertex[0], vertex[1], vertex[2])
                is_node_in_nodes_list = False
                for node_else in nodes:
                    dx = abs(node[0] - node_else[0])
                    dy = abs(node[1] - node_else[1])
                    dz = abs(node[2] - node_else[2])
                    if dx < 10**(-6) and dy < 10**(-6) and dz < 10**(-6):
                        is_node_in_nodes_list = True
                        break
                if not is_node_in_nodes_list:
                    nodes.append(node)
     
        for facet in solid_stl_obj.facets:
            shell = []
            for vertex in facet.vertices:
                node = (vertex[0], vertex[1], vertex[2])
                
                for node_id, node_else in enumerate(nodes):
                    dx = abs(node[0] - node_else[0])
                    dy = abs(node[1] - node_else[1])
                    dz = abs(node[2] - node_else[2])
                    if dx < 10**(-6) and dy < 10**(-6) and dz < 10**(-6):
                        break
                    
                shell.append(node_id)
                
            shells.append(shell)
        
        return cls(nodes, shells)
        
    
    @staticmethod
    def do_intersect(shell_i, shell_j):
        # intersections = 0
        # for node_id in shell_i:
        #     if node_id in shell_j:
        #         intersections = intersections + 1
                
        # if intersections > 1 :
        #     return True
        # else:
        #     return False
        
        return sum(node_id in shell_j for node_id in shell_i)>1    
    
    def locate_shells_adjacency(self):
        shells = self.shells
        neighbours=[]
        for i, shell_i in enumerate(shells):
            neighbours.append([])
            for j, shell_j in enumerate(shells):
                if i != j and self.do_intersect(shell_i, shell_j):
                    neighbours[-1].append(j)
        
        return neighbours

    def eliminate_adjacency(self, id_list1, id_list2):
        
        for id in id_list1:
                        
            self.shell_neighbours[id] = [
                id for id in self.shell_neighbours[id] if id not in id_list2
            ]
                    
        for id in id_list2:
            
            self.shell_neighbours[id] = [
                id for id in self.shell_neighbours[id] if id not in id_list1
            ]
    
    def add_extra_neighbours(self):
         
        old_shell_neighbours = {}
        for id_i in range(self.shell_num):
            old_shell_neighbours[id_i] = self.shell_neighbours[id_i].copy()
            for id_j in self.shell_neighbours[id_i]:
                old_shell_neighbours[id_j] = self.shell_neighbours[id_j].copy()
        
        for id_i in range(self.shell_num):
            for id_j in old_shell_neighbours[id_i]:
                for id_k in old_shell_neighbours[id_j]: 
                    if id_k!=id_i and id_k not in self.shell_neighbours[id_i]:
                        self.shell_neighbours[id_i].append(id_k)
    
    def plot_shells(self, elevation=30, azimuth=-60):
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.view_init(elevation, azimuth)
        for shell in self.shells:
            x, y, z = [], [], []
            for id in shell:
                [x_i, y_i, z_i] = self.nodes[id]
                x.append(x_i)
                y.append(y_i)
                z.append(z_i)
            
            id = shell[0]
            [x_i, y_i, z_i] = self.nodes[id]
            x.append(x_i)
            y.append(y_i)
            z.append(z_i)    
            ax.plot3D(x, y, z, color='k')    
        
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()
    
    
    
    ### unsteady features ###
    
    def set_body_fixed_frame_origin(self, xo, yo, zo):
        self.ro = Vector((xo, yo, zo))
    
    def set_body_fixed_frame_orientation(self, roll, pitch, yaw):
        
        self.theta = Vector((roll, pitch, yaw))
        
        Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0],
                       [np.sin(yaw), np.cos(yaw), 0],
                       [0, 0, 1]])
        
        Ry = np.array([[np.cos(pitch), 0, np.sin(pitch)],
                       [0, 1, 0],
                       [-np.sin(pitch), 0, np.cos(pitch)]])
        
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(roll), -np.sin(roll)],
                       [0, np.sin(roll), np.cos(roll)]])
        
        """
            A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
            A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) =>
            A(t+Δt)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx) A(t)]^T =>
            R(t+Δt) = A(t)^T Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t+Δt) = R(t) Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
            
            if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
            initial orientation of the body is such that f' coincides with f
            
            if A(t=0) =/= I then
            A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx) =>
            A(t=0)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx)]^T =>
            R(t=0) = Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
            R(t=0) = Rx(Δθx) Ry(Δθy) Rz(Δθz)
            
            where Δθx, Δθy, and Δθz correspond to finite rotations
            
            for more information check master thesis documentation and
            https://en.wikipedia.org/wiki/Davenport_chained_rotations
        """

        self.R = Rx @ Ry @ Rz      
        
    def set_origin_velocity(self, Vo:Vector):
        self.Vo = Vo
    
    def set_angular_velocity(self, omega:Vector):
        self.omega = omega
    
    def move_body(self, dt):
        
        """
        We chose to work in the body-fixed frame of reference (and not in the inertial frame of reference). Any movement of the body in space isn't observed by someone that stands on the bodyfixed frame. Hence, every coordinate or position vector of the rigid body mesh doesn't change, unless it can move relative to the body-fixed frame (e.g. wake nodes, wake panels).
        
        This method change only the position vector of body-fixed frame's origin and the orintetation of the body-fixed frame.
        Note that origin's position vector and body-fixed frame's orientation are observed from the inertial frame of reference.
        
        
        
        
        A(t+Δt) = Az(t+Δt) Ay(t+Δt) Ax(t+Δt) A(t) =>
        A(t+Δt) = Az(Δθz) Ay(Δθy) Ax(Δθx) A(t) =>
        A(t+Δt)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx) A(t)]^T =>
        R(t+Δt) = A(t)^T Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
        R(t+Δt) = R(t) Rx(Δθx) Ry(Δθy) Rz(Δθz)
        
        where Δθx, Δθy, and Δθz correspond to infinitesimal rotations (Δθx, Δθy,Δθz --> 0)
        
        if A(t=0) = I => R(t=0) = A(t=0)^T = I^T = I
        initial orientation of the body is such that f' coincides with f
        
        if A(t=0) =/= I then
        A(t=0) = Az(Δθx) Ay(Δθx) Ax(Δθx) =>
        A(t=0)^T = [Az(Δθz) Ay(Δθy) Ax(Δθx)]^T =>
        R(t=0) = Ax(Δθx)^T Ay(Δθy)^T Az(Δθz)^T =>
        R(t=0) = Rx(Δθx) Ry(Δθy) Rz(Δθz)
        
        where Δθx, Δθy, and Δθz correspond to finite rotations
        
        for more information check master thesis documentation and
        https://en.wikipedia.org/wiki/Davenport_chained_rotations        
        """
        
        self.ro = self.ro + self.Vo*dt
        dtheta =  self.omega*dt
        self.theta = self.theta + dtheta


        Rz = np.array([[np.cos(dtheta.z), -np.sin(dtheta.z), 0],
                       [np.sin(dtheta.z), np.cos(dtheta.z), 0],
                       [0, 0, 1]])

        Ry = np.array([[np.cos(dtheta.y), 0, np.sin(dtheta.y)],
                       [0, 1, 0],
                       [-np.sin(dtheta.y), 0, np.cos(dtheta.y)]])

        Rx = np.array([[1, 0, 0],
                       [0, np.cos(dtheta.x), -np.sin(dtheta.x)],
                       [0, np.sin(dtheta.x), np.cos(dtheta.x)]])

        self.R = self.R @ Rx @ Ry @ Rz 
    
    def move_node(self, node_id, v_rel, dt):
        
        """
        v_rel = velocity relative to the inertial frame of reference. (e.g. v_rel = wind speed => v_rel + v_body = v_freestream)    
                
        Ιδιότητα Πρόσθεσης:
        ω_f"/f = ω_f"/f' + ω_f'/f
        
        1)
        f"=f  =>  ω_f"/f' = ω_f/f' και ω_f"/f = 0 
        => 0 = ω_f/f' + ω_f'/f
        => ω_f/f' = - ω_f'/f
        
        2)
        F:αδρανειακό συστημα αναφοράς, F:OXYZ
        f:μεταφερόμενο σύστημα αναφοράς, f:O'xyz
        f:προσδεδεμένο σύστημα αναφοράς στο στερεό, f':O'x'y'z'
        
        Από ιδιότητα πρόσθεσης: 
        ω_f'/F = ω_f'/f + ω_f/F
        (ω_f/F=0) => ω_f'/F = ω_f'/f
        
        από 1) ω_f/f' = - ω_f'/f => ω_f'/F = - ω_F/f'
        (ω_F/F = ω_F/f' + ω_f'/F => 0 = ω_F/f' + ω_f'/F => ω_f'/F = - ω_F/f')
        
        3)
        Μεταφορά παραγώγου διανύσματος από το σύστημα αναφοράς f' στο σύστημα ανφοράς f:
        (r_dot)f = (r_dot)f' + (ω)f'/f X r
        
        
        r = ro' + r' => r_dot = ro'_dot + r'_dot =>
        (r_dot)F = (ro'_dot)F + (r'_dot)F => ...
        (r'_dot)f' = (r_dot)F - [ (ro'_dot)F + (ω)f'/F X r']
        r':position vector of node measured from bodyfixed frame of reference f'
        r: position vector of node measured from inertial frame of reference F
        ω: angular velocity of body-fixed frame observed from inetial frame of reference F (or from frame of refernce f)
         
        (Δες και χειρόγραφες σημειώσεις για περισσότερες λεπτομέρειες)        
        """
        node = self.nodes[node_id]
        
        # r': position vector of node measured from bodyfixed frame of reference f'
        r = Vector(node)  
        
        # (r'_dot)f' 
        # r':position vector of node measured from bodyfixed frame of reference f'
        v = v_rel - (self.Vo + Vector.cross_product(self.omega, r))
        
        dr = v*dt
        dr = dr.transformation(self.R.T)
        r = r+dr
        
        self.nodes[node_id] = (r.x, r.y, r.z)

    def move_nodes(self, node_id_list, v_rel, dt):
        
        for node_id in node_id_list:
            self.move_node(node_id, v_rel, dt)
    
    def convect_node(self, node_id, velocity:Vector, dt):
        
        node = self.nodes[node_id]
        r = Vector(node)
        dr = velocity * dt
        r = r + dr
        self.nodes[node_id] = (r.x, r.y, r.z)
    
    def convect_nodes(self, node_id_list, velocity_list, dt):
        for i in range(len(node_id_list)):
            node_id = node_id_list[i]
            velocity = velocity_list[i]
            self.convect_node(node_id, velocity, dt)
        
        
class AeroMesh(Mesh):
    
    def __init__(self, nodes: list, shells: list, nodes_ids:dict):
        super().__init__(nodes, shells)
        
        self.nodes_ids = nodes_ids
        self.shells_ids = {}
        self.TrailingEdge = {}
        self.wake_sheddingShells = {}
        
        self.set_shells_ids()
        self.find_TrailingEdge()
        self.set_wake_sheddingShells()
        
        self.free_TrailingEdge()
        # self.free_LeadingEdge()
        self.eliminate_main_surface_wing_tips_adjacency()
        self.eliminate_body_wake_adjacency()
    
    def find_TrailingEdge(self):
        te = self.nodes_ids["trailing edge"]
        ss = self.nodes_ids["suction side"]
        
        
        ps = self.nodes_ids["pressure side"]
        
        SS_TE = [shell_id for shell_id, shell in enumerate(self.shells)
                 if sum(node_id in te for node_id in shell)>1 
                 and sum(node_id in ss for node_id in shell) >2
                ]
        PS_TE = [shell_id for shell_id, shell in enumerate(self.shells)
                 if sum(node_id in te for node_id in shell)>1 
                 and sum(node_id in ps for node_id in shell) > 2
                ]
                
        self.TrailingEdge = {"suction side": SS_TE, "pressure side": PS_TE}

    def free_TrailingEdge(self):
        """
        if trailing edge shells (or panels) share trailing edge nodes then
        the method "find_shell_neighbours" will assume that suction side trailing shells and pressure side trailing shells are neighbours.
        In panel methods we assume that traling edge is a free edge.
        "free_TrailingEdge" method will remove false neighbour ids from the attribute shell_neighbour
        """
        
        list1 = [id for id in self.TrailingEdge["suction side"]]
        list2 = [id for id in self.TrailingEdge["pressure side"]]
        self.eliminate_adjacency(id_list1=list1, id_list2=list2)

    def free_LeadingEdge(self):
        num_shells_le = len(self.TrailingEdge["suction side"])
        LE_SS = self.shells_ids["suction side"][-num_shells_le:]
        LE_PS = self.shells_ids["pressure side"][0:num_shells_le]
        self.eliminate_adjacency(LE_SS, LE_PS)
        
    def find_suction_side(self):
        suction_side_nodes_ids = self.nodes_ids["suction side"]
        suction_side_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in suction_side_nodes_ids for node_id in shell)>2
        ]
        
        return suction_side_shells_ids
    
    def find_pressure_side(self):
        pressure_side_nodes_ids = self.nodes_ids["pressure side"]
        pressure_side_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in pressure_side_nodes_ids for node_id in shell)>2
        ]
        
        return pressure_side_shells_ids
    
    def find_right_wing_tip(self):
        right_wing_tip_nodes_ids = self.nodes_ids["right wing tip"]
        right_wing_tip_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in right_wing_tip_nodes_ids for node_id in shell)>2
        ]
        return right_wing_tip_shells_ids
    
    def find_left_wing_tip(self):
        left_wing_tip_nodes_ids =  self.nodes_ids["left wing tip"]
        left_wing_tip_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in left_wing_tip_nodes_ids for node_id in shell)>2
        ]
        return left_wing_tip_shells_ids
    
    def find_wake(self):
        wake_nodes_ids = self.nodes_ids["wake"]
        wake_shells_ids = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in wake_nodes_ids for node_id in shell)>2
        ]
        
        return wake_shells_ids
    
    def set_shells_ids(self):
        
        suction_side_shells_ids = self.find_suction_side()
        pressure_side_shells_ids = self.find_pressure_side()
        main_surface_shells_ids = suction_side_shells_ids\
            +pressure_side_shells_ids
        
        right_wing_tip_shells_ids = self.find_right_wing_tip()
        left_wing_tip_shells_ids = self.find_left_wing_tip()
        wing_tips_shells_ids = right_wing_tip_shells_ids \
            + left_wing_tip_shells_ids
        
        body_shells_ids = main_surface_shells_ids + wing_tips_shells_ids
        
        wake_shells_ids = self.find_wake()  
        
        self.shells_ids = {
            "body": body_shells_ids,
            "main surface": main_surface_shells_ids,
            "suction side": suction_side_shells_ids,
            "pressure side": pressure_side_shells_ids,
            "wing tips": wing_tips_shells_ids,
            "right tip": right_wing_tip_shells_ids,
            "left tip" : left_wing_tip_shells_ids,
            "wake": wake_shells_ids
        }
    
    def init_wake_sheddingShells(self):
        # if not empty dict
        if self.TrailingEdge:

            for j in range(len(self.TrailingEdge["suction side"])):

                self.wake_sheddingShells[self.TrailingEdge["suction side"][j]] = []

                self.wake_sheddingShells[self.TrailingEdge["pressure side"][j]] = []

    def set_wake_sheddingShells(self):
        
        self.init_wake_sheddingShells()
        
        for wake_line_id in range(len(self.nodes_ids["wake lines"])-1):
            
            line = self.nodes_ids["wake lines"][wake_line_id]
            next_line = self.nodes_ids["wake lines"][wake_line_id + 1]
            
            for te_shell_id in self.wake_sheddingShells:
                te_shell = self.shells[te_shell_id]
                
                if sum(
                    node_id in [line[0], next_line[0]]
                    for node_id in te_shell
                ) == 2:
                        
                        self.wake_sheddingShells[te_shell_id] = [
                            shell_id for shell_id in self.shells_ids["wake"]
                            if sum(
                                node_id in [*line, *next_line]
                                for node_id in self.shells[shell_id]
                            ) > 2
                        ]        
      
    def add_extra_neighbours(self):
         
        old_shell_neighbours = {}
        for id_i in self.shells_ids["body"]:
            old_shell_neighbours[id_i] = self.shell_neighbours[id_i].copy()
            for id_j in self.shell_neighbours[id_i]:
                old_shell_neighbours[id_j] = self.shell_neighbours[id_j].copy()
        
        for id_i in self.shells_ids["body"]:
            for id_j in old_shell_neighbours[id_i]:
                for id_k in old_shell_neighbours[id_j]: 
                    if id_k!=id_i and id_k not in self.shell_neighbours[id_i]:
                        self.shell_neighbours[id_i].append(id_k)      

    def eliminate_main_surface_wing_tips_adjacency(self):
        
        if "wing tips" in self.shells_ids and "main surface" in self.shells_ids:
            
            self.eliminate_adjacency(
                self.shells_ids["wing tips"], self.shells_ids["main surface"]
            )
    
    def eliminate_body_wake_adjacency(self):
        
        if "main surface" in self.shells_ids and "wake" in self.shells_ids:
            self.eliminate_adjacency(
                self.shells_ids["main surface"], self.shells_ids["wake"]
            )
    
    def give_near_root_shells_id(self):
        
        near_root_nodes_id = []
        for node_id, node in enumerate(self.nodes):
            if abs(node[1]) < 10**(-10):
                near_root_nodes_id.append(node_id)
        
        if self.shells_ids:
            body_shells_id = self.shells_ids["body"]
        else:
            body_shells_id = np.arange(len(self.shells))
        
        near_root_shells_id = []
        for shell_id in body_shells_id:
            for node_id in self.shells[shell_id]:
                if node_id in near_root_nodes_id:
                    near_root_shells_id.append(shell_id)
                    break
        
        return near_root_shells_id
    
    def give_leftSide_near_tip_shells_id(self):

        if self.shells_ids["left tip"]:
            left_wing_tip_shells_ids = self.shells_ids["left tip"]
        else:
            left_wing_tip_shells_ids = self.find_left_wing_tip()


        left_wing_tip_nodes_ids =  self.nodes_ids["left wing tip"]
        # wing tip shells' ids + near tip shells' ids
        left_side_near_tip_shells_id = [
            shell_id for shell_id, shell in enumerate(self.shells)
            if sum(node_id in left_wing_tip_nodes_ids for node_id in shell)>1
        ]

        # remove ids from wing tip shells
        left_side_near_tip_shells_id = [
            id for id in left_side_near_tip_shells_id if id not in left_wing_tip_shells_ids
        ]

        return left_side_near_tip_shells_id
    
    def locate_VSAERO_adjacency(self):
        
        """
        this function locates the neighbours of shells, following the notation of NASA Contractor Report 4023 "Program VSAERO theory Document,
        A Computer Program for Calculating Nonlinear Aerodynamic Characteristics
        of Arbitrary Configurations, Brian Maskew"
        """
        
        span_wise_nodes = len(self.nodes_ids["trailing edge"])
        chrod_wise_nodes = len(self.nodes_ids["main surface"])//span_wise_nodes
        span_wise_shells = span_wise_nodes - 1
        chrod_wise_shells = chrod_wise_nodes
        
        j_max = span_wise_shells - 1
        i_max = chrod_wise_shells - 1
        
        def shell_id(chord_wise_index, span_wise_index):            
            i = chord_wise_index
            j = span_wise_index
            shell_id = i*(j_max+1) + j
            return shell_id
              
              
        for i in range(1, i_max):
            for j in range(1, j_max):
                self.shell_neighbours[shell_id(i, j)] = [
                    shell_id(i, j-1),
                    shell_id(i+1, j),
                    shell_id(i, j+1),
                    shell_id(i-1, j)
                ]
        
        for i in range(1, i_max):
            j = 0
            self.shell_neighbours[shell_id(i, j)] = [
                -1,
                shell_id(i+1, j),
                shell_id(i, j+1),
                shell_id(i-1, j),
                
                shell_id(i, j+2)
            ]
            
            j = j_max
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                shell_id(i+1, j),
                -1,
                shell_id(i-1, j),
                
                shell_id(i, j-2)
            ]
            
        for j in range(1, j_max):
            i = 0
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                shell_id(i+1, j),
                shell_id(i, j+1),
                -1,
                
                shell_id(i+2, j)
            ]
            
            i = i_max
            self.shell_neighbours[shell_id(i, j)] = [
                shell_id(i, j-1),
                -1,
                shell_id(i, j+1),
                shell_id(i-1, j),
                
                shell_id(i-2, j)
            ]
            
        i, j = 0, 0
        self.shell_neighbours[shell_id(i, j)] = [
            -1,
            shell_id(i+1, j),
            shell_id(i, j+1),
            -1,
            
            shell_id(i, j+2),
            shell_id(i+2, j)
        ]
        
        i, j = 0, j_max
        self.shell_neighbours[shell_id(i, j)] = [
            shell_id(i, j-1),
            shell_id(i+1, j),
            -1,
            -1,
            
            shell_id(i, j-2),
            shell_id(i+2, j)
        ]
        
        i, j = i_max, 0
        self.shell_neighbours[shell_id(i, j)] = [
            -1,
            -1,
            shell_id(i, j+1),
            shell_id(i-1, j),
            
            shell_id(i, j+2),
            shell_id(i-2, j)
        ]
        
        i, j = i_max, j_max
        self.shell_neighbours[shell_id(i, j)] = [
            shell_id(i, j-1),
            -1,
            -1,
            shell_id(i-1, j),
            
            shell_id(i-2, j),
            shell_id(i, j-2)
        ]
      
    def give_ChordWiseStrips(self):
        """
        this function returns an array Ny X Nx array of panel ids, where Nx is the number of chordwise panels and Ny the number of spanwise panels.
        
        j-th row corresponds to j-th chordwise strip
        i-th column corresponds to i-th spanwise strip
        
        This function is meaningful only for quadrilateral structured meshes
        """

        span_wise_nodes = len(self.nodes_ids["trailing edge"])
        chrod_wise_nodes = len(self.nodes_ids["main surface"])//span_wise_nodes
        span_wise_shells = span_wise_nodes - 1
        chrod_wise_shells = chrod_wise_nodes

        j_max = span_wise_shells - 1
        i_max = chrod_wise_shells - 1

        def shell_id(chord_wise_index, span_wise_index):            
            i = chord_wise_index
            j = span_wise_index
            shell_id = i*(j_max+1) + j
            return shell_id

        ChordWiseStrips = np.zeros((j_max+1, i_max+1), dtype=int)
        for i in range(i_max + 1):
            for j in range(j_max + 1):
                ChordWiseStrips[j][i] = shell_id(i,j)

        return ChordWiseStrips
    
       
    ### unsteady features  ###
    def initialize_wake_nodes(self):
        new_nodes_ids_list = []
        for id in self.nodes_ids['trailing edge']:
            
            new_node_id = len(self.nodes)
            new_node = self.nodes[id]
            
            new_nodes_ids_list.append(new_node_id)
            self.nodes.append(new_node)
                        
        self.nodes_ids["wake"] = new_nodes_ids_list
                    
        self.nodes_ids["wake lines"] = np.array(
            [[id] for id in self.nodes_ids["wake"]]
        )
             
    def add_wakeNodes(self):     
        new_nodes_ids_list = []
        
        for id in self.nodes_ids['trailing edge']:
            
            new_node_id = self.nodes_ids["wake"][-1] + 1
            new_node = self.nodes[id]  # new_node coords
            
            self.nodes.append(new_node)
            self.nodes_ids["wake"].append(new_node_id)
            
            new_nodes_ids_list.append(new_node_id)
            
        self.nodes_ids["wake lines"] = np.column_stack(
            (new_nodes_ids_list, self.nodes_ids["wake lines"])
        )
        
    def add_wakeShells(self, type="quadrilateral"):
                
        def node_id(chord_wise_index, span_wise_index):
            i = chord_wise_index
            j = span_wise_index
                 
            return self.nodes_ids["wake lines"][j][i]
            
        def add_shell(*node_ids, reverse_order=False):
            # node_id_list should be in counter clock wise order
            
            if reverse_order:
                node_ids = list(node_ids)
                node_ids.reverse()
                
            if len(node_ids) == 4:
                if type == "quadrilateral":
                    self.shells.append(list(node_ids))
                elif type == "triangular":
                    index = node_ids
                    self.shells.append([index[0], index[1], index[2]])
                    self.shells.append([index[2], index[3], index[0]])
                    
            elif len(node_ids) == 3:
                self.shells.append(list(node_ids))
                
        ny = len(self.nodes_ids["trailing edge"])
               
        if type=="quadrilateral":
            
            for j in range(ny-1):
                
                add_shell(
                    node_id(1, j),
                    node_id(0, j),
                    node_id(0, j+1),
                    node_id(1, j+1)
                )
                
                if self.shells_ids["wake"] == []:
                    id = self.shells_ids["body"][-1] + 1
                else:
                    id = self.shells_ids["wake"][-1] + 1
                
                # Προσοχή θα ανανεωθεί και το λεξικό self.panels_ids
                self.shells_ids["wake"].append(id)
                
                value_list_id = id
                
                key_id = self.TrailingEdge["pressure side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                
                key_id = self.TrailingEdge["suction side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                
        elif type=="triangular":
            
            # right side
            for j in range((ny-1)//2):
                
                add_shell(
                    node_id(1, j+1),
                    node_id(1, j),
                    node_id(0, j),
                    node_id(0, j+1)
                )
                
                if self.shells_ids["wake"] == []:
                    id = self.shells_ids["body"][-1] + 1
                else:
                    id = self.shells_ids["wake"][-1] + 1
                
                # Προσοχή θα ανανεωθεί και το λεξικό self.panels_ids
                self.shells_ids["wake"].append(id)
                self.shells_ids["wake"].append(id+1)
                
                value_list_id = id
                
                key_id = self.TrailingEdge["pressure side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                self.wake_sheddingShells[key_id].append(value_list_id + 1)
                
                key_id = self.TrailingEdge["suction side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                self.wake_sheddingShells[key_id].append(value_list_id+1)
                
            # left side
            for j in range((ny-1)//2, ny-1):
                
                add_shell(
                    node_id(1, j),
                    node_id(0, j),
                    node_id(0, j+1),
                    node_id(1, j+1)
                )
                                
                if self.shells_ids["wake"] == []:
                    id = self.shells_ids["body"][-1] + 1
                else:
                    id = self.shells_ids["wake"][-1] + 1
                
                # Προσοχή θα ανανεωθεί και το λεξικό self.panels_ids
                self.shells_ids["wake"].append(id)
                self.shells_ids["wake"].append(id+1)
                
                value_list_id = id
                
                key_id = self.TrailingEdge["pressure side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                self.wake_sheddingShells[key_id].append(value_list_id + 1)
                
                key_id = self.TrailingEdge["suction side"][j]
                self.wake_sheddingShells[key_id].append(value_list_id)
                self.wake_sheddingShells[key_id].append(value_list_id+1)
    
    def ravel_wake(self, v_rel, dt, type="quadrilateral"):
        
        """
        Performs the same operation as shed_wake() but without using a sheding factor. With this method wake panels aren't streched
        """
        
        if self.nodes_ids["wake"] == []:
            self.initialize_wake_nodes()
        
        nodes_to_shed =  self.nodes_ids["wake lines"].flatten()
        self.move_nodes(nodes_to_shed, v_rel, dt)
        self.add_wakeNodes()
        self.add_wakeShells(type)
                     
    def shed_wake(self, v_rel, dt, wake_shed_factor=1, type="quadrilateral"):
                
        if self.nodes_ids["wake"] == []:
            self.initialize_wake_nodes()
            
            nodes_to_shed =  self.nodes_ids["wake lines"].flatten()
            self.move_nodes(nodes_to_shed, v_rel, dt*wake_shed_factor)
            
            self.add_wakeNodes()
            self.add_wakeShells(type)
            
        else:
            
            nodes_to_shed =  self.nodes_ids["wake lines"][:, 0].flatten()
            self.move_nodes(nodes_to_shed, v_rel, dt*wake_shed_factor)
            
            nodes_to_shed =  self.nodes_ids["wake lines"][:, 1:].flatten()
            self.move_nodes(nodes_to_shed, v_rel, dt)
            
            self.add_wakeNodes()
            self.add_wakeShells(type)
    
    def nodes_to_convect(self):
        nodes_to_convect = self.nodes_ids["wake lines"][:, 1:].flatten()
        node_id_list = list(nodes_to_convect)
        return node_id_list
            
    def convect_wake(self, velocity_list, dt):
        node_id_list = self.nodes_to_convect()
        self.convect_nodes(node_id_list, velocity_list, dt)
        
                          
class PanelMesh(Mesh):
    def __init__(self, nodes:list, shells:list):
        super().__init__(nodes, shells)
        self.panels = None
        self.panels_num = None
        self.panel_neighbours = self.shell_neighbours
        self.CreatePanels()         
        
    def CreatePanels(self):
        panels = []
        for shell_id, shell in enumerate(self.shells):
            vertex = []
            for node_id in shell:
                node = self.nodes[node_id]
                vertex.append(Vector(node))

            if len(vertex) == 3:
                panels.append(TriPanel(vertex[0], vertex[1], vertex[2]))
            elif len(vertex) == 4:
                panels.append(QuadPanel(vertex[0], vertex[1],
                                        vertex[2], vertex[3]))
            
            panels[-1].id = shell_id
                
        self.panels = panels
        self.panels_num = len(panels)

    def give_neighbours(self, panel):
        
        neighbours_id_list = self.panel_neighbours[panel.id]
        
        # neighbours_list = []
        # for id in neighbours_id_list:
        #     neighbours_list.append(self.panels[id])
         
        neighbours_list = [self.panels[id] for id in neighbours_id_list]
        
        return neighbours_list
    
    def plot_panels(self, elevation=30, azimuth=-60):
        ax = plt.axes(projection='3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.view_init(elevation, azimuth)
        
        for panel in self.panels:
            
            r_vertex = panel.r_vertex
            
            # plot panels
            if panel.num_vertices == 3:
                x = [r_vertex[0].x, r_vertex[1].x, r_vertex[2].x, r_vertex[0].x]
                y = [r_vertex[0].y, r_vertex[1].y, r_vertex[2].y, r_vertex[0].y]
                z = [r_vertex[0].z, r_vertex[1].z, r_vertex[2].z, r_vertex[0].z]
                ax.plot3D(x, y, z, color='k')
                
            elif panel.num_vertices == 4:
                
                x = [r_vertex[0].x, r_vertex[1].x, r_vertex[2].x, r_vertex[3].x,
                    r_vertex[0].x]
                y = [r_vertex[0].y, r_vertex[1].y, r_vertex[2].y, r_vertex[3].y,
                    r_vertex[0].y]
                z = [r_vertex[0].z, r_vertex[1].z, r_vertex[2].z, r_vertex[3].z,
                    r_vertex[0].z]
                ax.plot3D(x, y, z, color='k') 
                
            # plot normal vectors
            r_cp = panel.r_cp
            n = panel.n
            l = panel.l
            m = panel.m
            scale = 0.1 * panel.char_length
            n = n * scale
            l = l * scale
            m = m * scale
            ax.scatter(r_cp.x, r_cp.y, r_cp.z, color='k', s=5)
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, n.x, n.y, n.z, color='r')
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, l.x, l.y, l.z, color='b')
            ax.quiver(r_cp.x, r_cp.y, r_cp.z, m.x, m.y, m.z, color='g')
                
                   
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()

    def plot_mesh(self, elevation=30, azimuth=-60):
        shells = []
        vert_coords = []
        
        for panel in self.panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
                vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
            shells.append(shell)
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n for panel in self.panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        poly3 = Poly3DCollection(shells, facecolor=facecolor)
        ax.add_collection(poly3)
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()
    
    
    ### unsteady features  ###
      
    def move_panel(self, panel_id, v_rel, dt):
        
        panel = self.panels[panel_id]
        panel.move(v_rel, self.Vo, self.omega, self.R, dt)
        self.panels[panel_id] = panel
    
    def move_panels(self, panel_id_list, v_rel, dt):
        
        for id in panel_id_list:
            self.move_panel(id, v_rel, dt)      
           
    def plot_mesh_inertial_frame(self, elevation=30, azimuth=-60):
        shells = []
        vert_coords = []
                
        for panel in self.panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                r = self.ro + r_vertex.transformation(self.R)
                shell.append((r.x, r.y, r.z))
                vert_coords.append([r.x, r.y, r.z])
            shells.append(shell)
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n.transformation(self.R) for panel in self.panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        poly3 = Poly3DCollection(shells, facecolor=facecolor)
        ax.add_collection(poly3)
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        
        # Body-fixed frame of reference f'
        ex = Vector((1, 0, 0))  # ex'
        ey = Vector((0, 1, 0))  # ey'
        ez = Vector((0, 0, 1))  # ez'
        
        ex = ex.changeBasis(self.R)
        ey = ey.changeBasis(self.R)
        ez = ez.changeBasis(self.R)
        
        ro = self.ro
        ax.quiver(ro.x, ro.y, ro.z, ex.x, ex.y, ex.z, color='b', label="ex'")
        ax.quiver(ro.x, ro.y, ro.z, ey.x, ey.y, ey.z, color='g', label="ey'")
        ax.quiver(ro.x, ro.y, ro.z, ez.x, ez.y, ez.z, color='r', label="ez'")
        
        
        # Inertial frame of reference
        ex = Vector((1, 0, 0))  # eX
        ey = Vector((0, 1, 0))  # eY
        ez = Vector((0, 0, 1))  # eZ
        
        ax.quiver(0, 0, 0, ex.x, ex.y, ex.z, color='b', label='ex')
        ax.quiver(0, 0, 0, ey.x, ey.y, ey.z, color='g', label='ey')
        ax.quiver(0, 0, 0, ez.x, ez.y, ez.z, color='r', label='ez')
        
        vert_coords.append([ex.x, ex.y, ex.z])
        vert_coords.append([ey.x, ey.y, ey.z])
        vert_coords.append([ez.x, ez.y, ez.z])
       
       
        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()
    
    def plot_mesh_bodyfixed_frame(self, elevation=30, azimuth=-60):
        shells = []
        vert_coords = []
                
        for panel in self.panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
                vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
            shells.append(shell)
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n for panel in self.panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        poly3 = Poly3DCollection(shells, facecolor=facecolor)
        ax.add_collection(poly3)
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        
        # Body-fixed frame of reference f'
        ex = Vector((1, 0, 0))  # ex'
        ey = Vector((0, 1, 0))  # ey'
        ez = Vector((0, 0, 1))  # ez'
        
        ax.quiver(0, 0, 0, ex.x, ex.y, ex.z, color='b', label="ex'")
        ax.quiver(0, 0, 0, ey.x, ey.y, ey.z, color='g', label="ey'")
        ax.quiver(0, 0, 0, ez.x, ez.y, ez.z, color='r', label="ez'")
        
        
        # Inertial frame of reference F
        ro = -self.ro  # ro: r_oo' -> r_o'o = -roo'  
        ex = Vector((1, 0, 0))  # eX
        ey = Vector((0, 1, 0))  # eY
        ez = Vector((0, 0, 1))  # eZ
        
        ro = ro.changeBasis(self.R.T)
        ex = ex.changeBasis(self.R.T)
        ey = ey.changeBasis(self.R.T)
        ez = ez.changeBasis(self.R.T)
        
        
        ax.quiver(ro.x, ro.y, ro.z, ex.x, ex.y, ex.z, color='b', label='ex')
        ax.quiver(ro.x, ro.y, ro.z, ey.x, ey.y, ey.z, color='g', label='ey')
        ax.quiver(ro.x, ro.y, ro.z, ez.x, ez.y, ez.z, color='r', label='ez')
        
        vert_coords.append([(ro+ex).x, (ro+ex).y, (ro+ex).z])
        vert_coords.append([(ro+ey).x, (ro+ey).y, (ro+ey).z])
        vert_coords.append([(ro+ez).x, (ro+ez).y, (ro+ez).z])
                    
        
        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()        
         

class PanelAeroMesh(AeroMesh, PanelMesh):
    
    def __init__(self, nodes: list, shells: list, nodes_ids: dict):
        super().__init__(nodes, shells, nodes_ids)
        
        self.panels_ids = self.shells_ids
        self.wake_sheddingPanels = self.wake_sheddingShells
    
    def free_TrailingEdge(self):
        super().free_TrailingEdge()
        self.panel_neighbours = self.shell_neighbours
      
    def give_near_root_panels(self):
        near_root_shells_id = super().give_near_root_shells_id()
        near_root_panels = []
        for id in near_root_shells_id:
            near_root_panels.append(self.panels[id])
        
        return near_root_panels
    
    def give_leftSide_near_root_panels(self):
        
        near_root_panels = self.give_near_root_panels()
        
        for panel in near_root_panels:
            if panel.r_cp.y < 0:
                near_root_panels.remove(panel)
        
        return near_root_panels

    def give_leftSide_near_tip_panels(self):

        leftSide_near_tip_panels = [
            self.panels[id] for id in self.give_leftSide_near_tip_shells_id()
        ]
        return leftSide_near_tip_panels
    
    ### unsteady features ###
    
    def add_wakePanels(self, type="quadrilateral"):
                
        num_TrailingEdge_panels = len(self.TrailingEdge["pressure side"])
        id_end = self.shells_ids["wake"][-1]
        if type == "quadrilateral":
            id_start = id_end - num_TrailingEdge_panels + 1
        elif type == "triangular":
            id_start = id_end - 2*num_TrailingEdge_panels + 1
                
        for shell_id in range(id_start, id_end+1):
            vertex = []
            for node_id in self.shells[shell_id]:
                node = self.nodes[node_id]
                vertex.append(Vector(node))
            
            if len(vertex) == 3:
                self.panels.append(TriPanel(vertex[0], vertex[1], vertex[2]))
            elif len(vertex) == 4:
                self.panels.append(QuadPanel(vertex[0], vertex[1],
                                        vertex[2], vertex[3]))
            
            self.panels[-1].id = shell_id
    
    def ravel_wake(self, v_rel, dt, type="quadrilateral"):
        
        """
        Performs the same operation as shed_wake() but without using a sheding factor. With this method wake panels aren't streched
        """
        
        if self.panels_ids["wake"] == []:
            
            super().ravel_wake(v_rel, dt, type)
            self.add_wakePanels(type)
        
        else:
            
            self.move_panels(self.panels_ids["wake"], v_rel, dt)
            super().ravel_wake(v_rel, dt, type) 
            self.add_wakePanels(type)
    
    def shed_wake(self, v_rel, dt, wake_shed_factor=1, type="quadrilateral"):
        
        if self.panels_ids["wake"] == []:
            
            super().shed_wake(v_rel, dt, wake_shed_factor, type)
            self.add_wakePanels(type)
                        
        else:
                       
            super().shed_wake(v_rel, dt, wake_shed_factor, type)
            self.add_wakePanels(type)
            self.update_wake_panel_vertices()
            
    def convect_wake(self, induced_velocity_function, dt):
        
        # create velocity list 
        node_id_list = self.nodes_to_convect()
        velocity_list = []
        body_panels = [self.panels[id] for id in self.panels_ids["body"]]
        wake_panels = [self.panels[id] for id in self.panels_ids["wake"]]
        for node_id in node_id_list:
            r_p = Vector(self.nodes[node_id])
            induced_velocity = induced_velocity_function(r_p, body_panels,
                                                         wake_panels)
            velocity_list.append(induced_velocity)
        
        # convect wake
        super().convect_wake(velocity_list, dt)
        
        # update panel vertices' location
        self.update_wake_panel_vertices()
    
    def update_wake_panel_vertices(self):
        for shell_id in self.shells_ids["wake"]:
            shell = self.shells[shell_id]          
            vertex_list = [Vector(self.nodes[node_id]) for node_id in shell]
            self.panels[shell_id].update_vertices_location(vertex_list)
                      
    def plot_mesh_inertial_frame(self, elevation=30, azimuth=-60,
                                 plot_wake=False):
        body_shells = []
        wake_shells = []
        vert_coords = []
        
        if self.shells_ids:
            body_panels = [self.panels[id] for id in self.shells_ids["body"]]
            if plot_wake:
                wake_panels = [self.panels[id] for id in self.shells_ids["wake"]]
        else:
            body_panels = self.panels
            
        for panel in body_panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                r = self.ro + r_vertex.transformation(self.R)
                shell.append((r.x, r.y, r.z))
                vert_coords.append([r.x, r.y, r.z])
            body_shells.append(shell)
            
        if plot_wake:
            for panel in wake_panels:
                shell=[]
                for r_vertex in panel.r_vertex:
                    r = self.ro + r_vertex.transformation(self.R)
                    shell.append((r.x, r.y, r.z))
                    vert_coords.append([r.x, r.y, r.z])
                wake_shells.append(shell)
            
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n.transformation(self.R) for panel in body_panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        body_collection = Poly3DCollection(body_shells, facecolor=facecolor)
        ax.add_collection(body_collection)
        
        if plot_wake:
            wake_collection = Poly3DCollection(wake_shells, alpha=0.1)
            wake_collection.set_edgecolors('k')
            ax.add_collection(wake_collection)
            
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        
        
        # Body-fixed frame of reference f'
        ex = Vector((1, 0, 0))  # ex'
        ey = Vector((0, 1, 0))  # ey'
        ez = Vector((0, 0, 1))  # ez'
        
        ex = ex.changeBasis(self.R)
        ey = ey.changeBasis(self.R)
        ez = ez.changeBasis(self.R)
        
        ro = self.ro
        ax.quiver(ro.x, ro.y, ro.z, ex.x, ex.y, ex.z, color='b', label="ex'")
        ax.quiver(ro.x, ro.y, ro.z, ey.x, ey.y, ey.z, color='g', label="ey'")
        ax.quiver(ro.x, ro.y, ro.z, ez.x, ez.y, ez.z, color='r', label="ez'")
        
        
        # Inertial frame of reference
        ex = Vector((1, 0, 0))  # eX
        ey = Vector((0, 1, 0))  # eY
        ez = Vector((0, 0, 1))  # eZ
        
        ax.quiver(0, 0, 0, ex.x, ex.y, ex.z, color='b', label='ex')
        ax.quiver(0, 0, 0, ey.x, ey.y, ey.z, color='g', label='ey')
        ax.quiver(0, 0, 0, ez.x, ez.y, ez.z, color='r', label='ez')
        
        vert_coords.append([ex.x, ex.y, ex.z])
        vert_coords.append([ey.x, ey.y, ey.z])
        vert_coords.append([ez.x, ez.y, ez.z])
        
        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()
    
    def plot_mesh_bodyfixed_frame(self, elevation=30, azimuth=-60,
                                  plot_wake=False):
        body_shells = []
        wake_shells = []
        vert_coords = []
        
        if self.shells_ids:
            body_panels = [self.panels[id] for id in self.shells_ids["body"]]
            if plot_wake:
                wake_panels = [self.panels[id] for id in self.shells_ids["wake"]]
        else:
            body_panels = self.panels
        
        for panel in body_panels:
            shell=[]
            for r_vertex in panel.r_vertex:
                shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
                vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
            body_shells.append(shell)
        
        if plot_wake:     
            for panel in wake_panels:
                shell=[]
                for r_vertex in panel.r_vertex:
                    shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
                    vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
                wake_shells.append(shell)
            
        
        
        
        light_vec = light_vector(magnitude=1, alpha=-45, beta=-45)
        face_normals = [panel.n for panel in body_panels]
        dot_prods = [-light_vec * face_normal for face_normal in face_normals]
        min = np.min(dot_prods)
        max = np.max(dot_prods)
        target_min = 0.2 # darker gray
        target_max = 0.6 # lighter gray
        shading = [(dot_prod - min)/(max - min) *(target_max - target_min) 
                   + target_min
                   for dot_prod in dot_prods]
        facecolor = plt.cm.gray(shading)
        
        ax = plt.axes(projection='3d')
        body_collection = Poly3DCollection(body_shells, facecolor=facecolor)
        ax.add_collection(body_collection)
        
        if plot_wake:
            wake_collection = Poly3DCollection(wake_shells, alpha=0.1)
            wake_collection.set_edgecolors('k')
            ax.add_collection(wake_collection)
            
        ax.view_init(elevation, azimuth)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        
        
        # Body-fixed frame of reference f'
        ex = Vector((1, 0, 0))  # ex'
        ey = Vector((0, 1, 0))  # ey'
        ez = Vector((0, 0, 1))  # ez'
        
        ax.quiver(0, 0, 0, ex.x, ex.y, ex.z, color='b', label="ex'")
        ax.quiver(0, 0, 0, ey.x, ey.y, ey.z, color='g', label="ey'")
        ax.quiver(0, 0, 0, ez.x, ez.y, ez.z, color='r', label="ez'")
        
        
        # Inertial frame of reference F
        ro = -self.ro  # ro: r_oo' -> r_o'o = -roo'  
        ex = Vector((1, 0, 0))  # eX
        ey = Vector((0, 1, 0))  # eY
        ez = Vector((0, 0, 1))  # eZ
        
        ro = ro.changeBasis(self.R.T)
        ex = ex.changeBasis(self.R.T)
        ey = ey.changeBasis(self.R.T)
        ez = ez.changeBasis(self.R.T)
        
        ax.quiver(ro.x, ro.y, ro.z, ex.x, ex.y, ex.z, color='b', label='ex')
        ax.quiver(ro.x, ro.y, ro.z, ey.x, ey.y, ey.z, color='g', label='ey')
        ax.quiver(ro.x, ro.y, ro.z, ez.x, ez.y, ez.z, color='r', label='ez')
        
        
        vert_coords.append([(ro+ex).x, (ro+ex).y, (ro+ex).z])
        vert_coords.append([(ro+ey).x, (ro+ey).y, (ro+ey).z])
        vert_coords.append([(ro+ez).x, (ro+ez).y, (ro+ez).z])
               
        
        vert_coords = np.array(vert_coords)
        x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
        ax.set_xlim3d(x.min(), x.max())
        ax.set_ylim3d(y.min(), y.max())
        ax.set_zlim3d(z.min(), z.max())
        set_axes_equal(ax)
        
        move_view_ = lambda event: move_view(event, ax)
        ax.figure.canvas.mpl_connect("key_press_event", move_view_)
        
        plt.show()
    
    def copy(self):
        nodes = deepcopy(self.nodes)
        shells = deepcopy(self.shells)
        nodes_ids = deepcopy(self.nodes_ids)
                
        mesh = PanelAeroMesh(nodes, shells, nodes_ids)
        
        (xo, yo, zo) = self.origin
        mesh.set_body_fixed_frame_origin(xo, yo, zo)
        
        (roll, pitch, yaw) = self.orientation
        mesh.set_body_fixed_frame_orientation(roll, pitch, yaw)
        
        mesh.set_origin_velocity(self.Vo)
        mesh.set_angular_velocity(self.omega)
        
        return mesh



        
if __name__=='__main__':
    from matplotlib import pyplot as plt
    from sphere import sphere
    nodes, shells = sphere(1, 10, 10, mesh_shell_type='quadrilateral')
    sphere_mesh = PanelMesh(nodes, shells)    
    sphere_mesh.plot_panels()
    sphere_mesh.plot_mesh_bodyfixed_frame()
    
    
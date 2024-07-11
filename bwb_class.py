from airfoil_class import Airfoil
from vector_class import Vector
from Algorithms import DenserAtBoundaries, cosspace, interpolation
import numpy as np
import csv


class WingCrossSection:
    
    def __init__(self, r_leadingEdge:Vector,
                 chord:float, twist:float, airfoil:Airfoil):
        
        self.r_leadingEdge = r_leadingEdge
        self.chord = chord
        self.twist = twist
        self.airfoil = airfoil
    
    def give_RotationMatrix(self, wing_cross_section_prev,
                            wing_cross_section_next):
        
        
        # local frame of reference before twist
        ex_local = Vector((1, 0, 0))
        
        if wing_cross_section_next == None:
            
            vec01 = (self.r_leadingEdge
                    - wing_cross_section_prev.r_leadingEdge)
            vec01_yz_project = Vector((0, vec01.y, vec01.z))
            ey_local = vec01_yz_project/vec01_yz_project.norm()
            
            z_scale = 1
            
        elif wing_cross_section_prev == None:
            
            vec12 = (wing_cross_section_next.r_leadingEdge
                    - self.r_leadingEdge)
            
            vec12_yz_project = Vector((0, vec12.y, vec12.z))
            ey_local = vec12_yz_project/vec12_yz_project.norm()
            
            z_scale = 1
            
        else:
            vec01 = (self.r_leadingEdge
                    - wing_cross_section_prev.r_leadingEdge)
            vec01_yz_project = Vector((0, vec01.y, vec01.z))
            vec_prev = vec01_yz_project/vec01_yz_project.norm()
            
            vec12 = (wing_cross_section_next.r_leadingEdge
                    - self.r_leadingEdge)
            
            vec12_yz_project = Vector((0, vec12.y, vec12.z))
            vec_next = vec12_yz_project/vec12_yz_project.norm()
            
            span_vec = (vec_prev+vec_next)/2
            
            ey_local = span_vec/span_vec.norm()
            
            z_scale = np.sqrt(2/(vec_prev * vec_next + 1))           
        
        ez_local = Vector.cross_product(ex_local, ey_local) * z_scale
        
        
        # local frame of reference after twist
        twist = np.deg2rad(-self.twist)
        
        # first method
        
        # Rotation matrix from axis and angle (https://en.wikipedia.org/wiki/Rotation_matrix)
        
        
        c, s = np.cos(twist), np.sin(twist)
        ux, uy, uz = ey_local.x, ey_local.y, ey_local.z
        
        R_twist = [
            [c + ux ** 2 * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s],
            [uy * ux * (1 - c) + uz * s, c + uy ** 2 * (1 - c), uy * uz * (1 - c) - ux * s],
            [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz ** 2 * (1 - c)]
            ]
        
        R_twist = np.array(R_twist)
        
        # Ά τρόπος
        # ex_local = ex_local.transformation(R_twist)
        # # ey_local = ey_local.transformation(R_twist) 
        # ez_local = ez_local.transformation(R_twist)
        
        # R_twisted = np.array([[ex_local.x, ey_local.x, ez_local.x],
        #                       [ex_local.y, ey_local.y, ez_local.y],
        #                       [ex_local.z, ey_local.z, ez_local.z]])
        
        #΄ Β τρόπος
        R_untwisted = np.array([[ex_local.x, ey_local.x, ez_local.x],
                              [ex_local.y, ey_local.y, ez_local.y],
                              [ex_local.z, ey_local.z, ez_local.z]])
        
        R_twisted = R_twist@R_untwisted
        
        
        # Rotate local frame of reference so z axis points downward
        theta = np.deg2rad(-90)
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(theta), -np.sin(theta)],
                       [0, np.sin(theta), np.cos(theta)]])
        
        RotationMatrix = R_twisted @ Rx        
                  
        return RotationMatrix

    def give_coords(self, x_percent, y_percent,
                    wing_cross_section_prev, wing_cross_section_next):
        
        r_le = self.r_leadingEdge
        R = self.give_RotationMatrix(wing_cross_section_prev,
                                     wing_cross_section_next)
        chord = self.chord
        r = Vector((x_percent*chord, y_percent*chord, 0))
        r_p = r_le + r.changeBasis(R)
        
        return (r_p.x, r_p.y, r_p.z)

    @classmethod
    def blend_WingCrossSections(cls, Xsection_0, Xsection_1, blend_fraction,
                                interpolation_type="linear"):
        
        # interpolation_type = "linear" or "cubic"
                
        Xsection_0_perc = 1 - blend_fraction
        Xsection_1_perc = blend_fraction
        
        name = f"{Xsection_0_perc * 100:.0f}% {Xsection_0.airfoil.name},\
            {Xsection_1_perc * 100:.0f}% {Xsection_1.airfoil.name}"
        
        
        x_coords = interpolation(
            Xsection_0.airfoil.x_coords, Xsection_1.airfoil.x_coords, blend_fraction, type=interpolation_type
        )
        
        y_coords = interpolation(
            Xsection_0.airfoil.y_coords, Xsection_1.airfoil.y_coords,
            blend_fraction, type=interpolation_type
        )
        
        # x_coords = Xsection_0_perc * Xsection_0.airfoil.x_coords \
        #     + Xsection_1_perc * Xsection_1.airfoil.x_coords
        
        # y_coords = Xsection_0_perc * Xsection_0.airfoil.y_coords \
        #     + Xsection_1_perc * Xsection_1.airfoil.y_coords
            
        blend_airfoil = Airfoil(name=name, x_coords=x_coords, y_coords=y_coords)
        
        r_leadingEdge = interpolation(
            Xsection_0.r_leadingEdge, Xsection_1.r_leadingEdge,
            blend_fraction, type="linear" # linear always
        )
        
        twist = interpolation(
            Xsection_0.twist, Xsection_1.twist, blend_fraction,
            type=interpolation_type
        )
        
        chord = interpolation(
            Xsection_0.chord, Xsection_1.chord, blend_fraction, type="linear" 
        ) # linear always
        
        blended_Xsection = WingCrossSection(
            r_leadingEdge, chord, twist, blend_airfoil
        )

        return blended_Xsection
     

class BWB:
    
    def __init__(self, name, wingXsection_list):
        self.name = name
        self.wingXsection_list = wingXsection_list
        self.compute_platform_chars()
            
    def compute_wingXsection_coord(self, id, x_perc, y_perc):
        wing_cross_section = self.wingXsection_list[id]
        
        if id == 0:
            wing_cross_section_prev = None
            wing_cross_section_next = self.wingXsection_list[id + 1]
            R = wing_cross_section.give_RotationMatrix(wing_cross_section_prev,
                                                       wing_cross_section_next)
            
        elif id == len(self.wingXsection_list)-1 or id==-1:
            wing_cross_section_prev = self.wingXsection_list[id - 1]
            wing_cross_section_next = None
            R = wing_cross_section.give_RotationMatrix(wing_cross_section_prev,
                                                       wing_cross_section_next)
        else:
            wing_cross_section_prev = self.wingXsection_list[id - 1]
            wing_cross_section_next = self.wingXsection_list[id + 1]
            R = wing_cross_section.give_RotationMatrix(wing_cross_section_prev,
                                                       wing_cross_section_next)
        
        r_le = wing_cross_section.r_leadingEdge
        chord = wing_cross_section.chord
        r = Vector((x_perc*chord, y_perc*chord, 0))
        r_p = r_le + r.changeBasis(R)
        # print(R)
        return (r_p.x, r_p.y, r_p.z)
    
    def mesh_line(self, x_percent_list, y_percent_list, resolution,
                  spacing="uniform"):
        
        """
        meshes the line that connects the i-th node of every wing cross section
        
        returns line nodes = [[x0 y0 z0]
                              [x1 y1 z1]
                              [x2 y2 z2]]
        """
        if len(x_percent_list) != len(self.wingXsection_list):
            print("number of x-coords must \
                  be equal to number of wing cross sections ")

        if len(y_percent_list) != len(self.wingXsection_list):
            print("number of y-coords must \
                  be equal to number of wing cross sections ")
        
        if spacing == "uniform":
            space = np.linspace
        elif spacing == "cosine":
            space = cosspace
        elif spacing == "beta distribution":
            space = lambda start, end, steps: DenserAtBoundaries(start, end,
                                                                 steps, -0.15)
                    
        wingXsectionS_nodes = []
        
        for id, wingXsection in enumerate(self.wingXsection_list):
            
            x_perc = x_percent_list[id]
            y_perc = y_percent_list[id]

            if id == 0:
                wingXsection_prev=None
                wingXsection_next=self.wingXsection_list[id + 1]
            
            elif id == len(self.wingXsection_list)-1 or id==-1:
                wingXsection_prev = self.wingXsection_list[id - 1]
                wingXsection_next = None
            
            else:
                wingXsection_prev = self.wingXsection_list[id - 1]
                wingXsection_next=self.wingXsection_list[id + 1]
                
            
            wingXsection_node = wingXsection.give_coords(
                x_perc, y_perc, wingXsection_prev, wingXsection_next
            )
            
            wingXsectionS_nodes.append(wingXsection_node)
            
        wingSpanwiseSectionS_nodes = []
        for i in range(len(wingXsectionS_nodes) - 1):
            wingSpanwiseSection_nodes = np.stack(
                [ space(wingXsectionS_nodes[i][dim],
                        wingXsectionS_nodes[i+1][dim],
                        resolution) for dim in [0, 1, 2] ], axis=1 )
            
            # avoid last node of i-th spanwise section to be the first 
            # node of the i-th + 1 spanwise section
            if not i== len(wingXsectionS_nodes)-2:
                wingSpanwiseSection_nodes = wingSpanwiseSection_nodes[0:-1, :]
            
            wingSpanwiseSectionS_nodes.append(wingSpanwiseSection_nodes)
        
        line_nodes = np.concatenate(wingSpanwiseSectionS_nodes)
        
        return line_nodes

    def mesh_body(self, ChordWise_resolution, SpanWise_resolution,
                  ChordWise_spacing = "cosine", SpanWise_spacing="uniform", shellType="quads", mesh_main_surface=True, mesh_tips=True,
                  mesh_wake=False, triangular_wake_mesh = False, wake_resolution=10, planar_wake=False, V_fs=Vector((1,0,0)),
                  wake_length_in_chords = 30, standard_mesh_format=True):
        
        for wingXsection in self.wingXsection_list:
            wingXsection.airfoil.repanel(
                ChordWise_resolution + 1, spacing=ChordWise_spacing
            )
        
        x_perc = np.array([wingXsection.airfoil.x_coords
                           for wingXsection
                           in self.wingXsection_list]).T
        
        y_perc = np.array([wingXsection.airfoil.y_coords
                           for wingXsection
                           in self.wingXsection_list]).T
        
        x_perc = x_perc[0:-1] # removedouble node at trailing edge
        y_perc = y_perc[0:-1]
        
        nodes_ofSpanWiseLineS = []
        for line_x_perc, line_y_perc in zip(x_perc, y_perc):
            
            nodes_ofSpanWiseLine = self.mesh_line(line_x_perc, line_y_perc,
                                                  resolution = SpanWise_resolution + 1,
                                                  spacing=SpanWise_spacing)
            
            nodes_ofSpanWiseLineS.append(nodes_ofSpanWiseLine)
        
        nodes = np.concatenate(nodes_ofSpanWiseLineS)       
        
        
        def node_id(chord_wise_index, span_wise_index):
            
            ny = SpanWise_resolution * (len(self.wingXsection_list)-1) + 1
            nx = len(nodes_ofSpanWiseLineS)
            global j_max, i_max
            j_max = ny-1
            i_max = nx  # i_max = nx-1 when double node at trailing edge
            i = chord_wise_index
            j = span_wise_index
            node_id = (j + i*ny)%(nx*ny)
            return node_id
        
         
        def add_shell(*node_ids, reverse_order=False):
            # node_id_list should be in counter clock wise order

            if reverse_order:
                node_ids = list(node_ids)
                node_ids.reverse()
                
            if len(node_ids) == 4:
                if shellType == "quads":
                    shells.append(list(node_ids))
                elif shellType == "trias":
                    index = node_ids
                    shells.append([index[0], index[1], index[2]])
                    shells.append([index[2], index[3], index[0]])
            
            elif len(node_ids) == 3:
                shells.append(list(node_ids))
            
        
        # call node_id() so i_max and j_max can be accessed
        node_id(0, 0)
        
        shells = []
        if mesh_main_surface:
            
            for i in range(i_max):
                for j in range(j_max):
                    add_shell(
                        node_id(i, j),
                        node_id(i+1, j),
                        node_id(i+1, j+1),
                        node_id(i, j+1)
                    )

        # if mesh_tips:
            
        #     # root or right tip
        #     # trailing edge
        #     add_shell(
        #         node_id(0, 0),
        #         node_id(i_max - 1, 0),
        #         node_id(1, 0)
        #     )
        #     # leading edge
        #     add_shell(
        #         node_id(i_max//2 + 1, 0),
        #         node_id(i_max//2, 0),
        #         node_id(i_max//2 - 1, 0)   
        #     )
            
        #     # tip or left tip
        #     # trailing edge
        #     add_shell(
        #         node_id(0, j_max),
        #         node_id(1, j_max),
        #         node_id(i_max-1, j_max)    
        #     )
        #     # leading edge
        #     add_shell(
        #         node_id(i_max//2 - 1, j_max),
        #         node_id(i_max//2, j_max),
        #         node_id(i_max//2 + 1, j_max)   
        #     )
            
        #     for i in range(1, i_max//2-1):
                
        #         # root or right tip
        #         add_shell(
        #             node_id(i, 0),
        #             node_id(i_max-i, 0),
        #             node_id(i_max - i - 1, 0),
        #             node_id(i+1, 0)   
        #         )
                
        #         # tip or left tip
        #         add_shell(
        #             node_id(i, j_max),
        #             node_id(i+1, j_max),
        #             node_id(i_max - i - 1, j_max),
        #             node_id(i_max-i, j_max)    
        #         )
        
        if mesh_tips:
            last_id = len(nodes) -1  # last node id
            id = last_id
            
            for j in [0, j_max]:
                # j=0 root or right tip
                # j-j_max yip or left tip
                
                if j==0:
                    add_face = add_shell
                elif j==j_max:
                    add_face = lambda *node_ids: add_shell(*node_ids,
                                                            reverse_order=True)
                
                # root or right tip   
                id = id + 1
                   
                # trailing edge   
                i = 0      
                node = (nodes[node_id(i+1, j)]
                        + nodes[node_id(i_max-i-1, j)]) / 2
                nodes = np.vstack([nodes, node])
                                
                add_face(
                    node_id(i, j),
                    id,
                    node_id(i+1, j),
                )
                
                add_face(
                    node_id(i, j),
                    node_id(i_max - i - 1, j),
                    id
                )
                
                for i in range(1, i_max//2 - 1):
                    id = id+1
                    
                    node = (nodes[node_id(i+1, j)]
                            + nodes[node_id(i_max-i-1, j)]) / 2
                                
                    nodes = np.vstack([nodes, node])
                                
                    add_face(
                        node_id(i, j),
                        id-1,
                        id,
                        node_id(i+1, j)   
                    )
                    
                    add_face(
                        id-1,
                        node_id(i_max-i, j),
                        node_id(i_max - i - 1, j),
                        id    
                    )
                    
                # leading edge
                i = i+1           
                add_face(
                    node_id(i, j),
                    id,
                    node_id(i+1, j)
                )
                
                add_face(
                    id,
                    node_id(i+2, j),  # node_id(i_max - i, j),
                    node_id(i+1, j),  # node_id(i_max - i - 1, j)     
                )
            
            extra_wing_tip_node_ids = [(last_id+1)+i for i in range(i_max-2)]
        
        if mesh_wake:
            
            x = np.zeros((wake_resolution+1, j_max+1))
            y = np.zeros_like(x)
            z = np.zeros_like(x)
                        
            C_root = max([wingXsection.chord
                              for wingXsection in self.wingXsection_list])
            
            if not planar_wake:
                
                for j in range(j_max+1):
                    
                    components = tuple(nodes[node_id(0, j)] - nodes[node_id(1, j)])
                    vec1 = Vector(components)
                    
                    components = tuple(nodes[node_id(i_max, j)]
                                    - nodes[node_id(i_max - 1, j)])
                    vec2 = Vector(components)
                    
                    bisector = vec1+vec2
                    
                    bisector = bisector/bisector.norm()
                    
                    bisector = bisector * C_root * wake_length_in_chords
                    
                    (x0, y0, z0) = nodes[node_id(0, j)]
                    x[:, j] = np.linspace(
                        x0, x0 + bisector.x, wake_resolution+1
                    )
                    y[:, j] = np.linspace(
                        y0, y0 + bisector.y, wake_resolution+1
                    )
                    z[:, j] = np.linspace(
                        z0, z0 + bisector.z, wake_resolution+1
                    )
                
            else:
                wake_direction_unit_vec = V_fs/V_fs.norm()
                max_chord = max(
                    [Xsection.chord for Xsection in self.wingXsection_list]
                )
                vec=wake_direction_unit_vec * wake_length_in_chords * max_chord 
                
                for j in range(j_max+1):
                    (x0, y0, z0) = nodes[node_id(0, j)]
                    x[:, j] = np.linspace(x0, x0 + vec.x, wake_resolution+1)
                    y[:, j] = np.linspace(y0, y0 + vec.y, wake_resolution+1)
                    z[:, j] = np.linspace(z0, z0 + vec.z, wake_resolution+1)
            
            num_body_nodes = len(nodes)
            for i in range(1, wake_resolution+1):
                for j in range(j_max + 1):
                    node = (x[i][j], y[i][j], z[i][j])
                    nodes = np.vstack([nodes, np.array(node)])
            
            def wake_node_id(chord_wise_index, span_wise_index):
                ny = SpanWise_resolution * (len(self.wingXsection_list)-1) + 1
                nx = wake_resolution + 1

                global i_max_wake
                i_max_wake = nx - 1
                i = chord_wise_index
                j = span_wise_index
                
                if i==0:
                    id = node_id(i, j)
                
                else:
                                   
                    id = j + i*ny + (num_body_nodes - ny)
                
                return id
            
            # call wake_node_id() so i_max_wake can be accessed
            wake_node_id(0, 0)
            
            if shellType == "quads":
                if triangular_wake_mesh:
                    shellType = "trias"
            
            for i in range(i_max_wake):
                for j in range(j_max):
                    
                    # right side
                    if j < j_max//2:
    
                        add_shell(
                            wake_node_id(i+1, j),
                            wake_node_id(i+1, j+1),
                            wake_node_id(i, j+1),
                            wake_node_id(i, j),
                            reverse_order=True
                        )
                    
                    # left side
                    else:
                        add_shell(
                            wake_node_id(i, j),
                            wake_node_id(i+1, j),
                            wake_node_id(i+1, j+1),
                            wake_node_id(i, j+1),
                            reverse_order=True
                        )
        
        
        # store nodes information in nodes_id_dict
        
        if mesh_main_surface:
            
            suction_side_nodes_ids = [node_id(i, j)
                                    for i in range(i_max//2 + 1)
                                    for j in range(j_max + 1)]
            
            pressure_side_nodes_ids = [node_id(i, j)
                                    for i in range(i_max//2, i_max + 1)
                                    for j in range(j_max + 1)]
                
            main_surface_nodes_ids = [node_id(i, j)
                                    for i in range(i_max)
                                    for j in range(j_max+1)]
            
            trailing_edge_nodes_ids = [node_id(0, j) for j in range(j_max+1)]
            
        else:
            suction_side_nodes_ids = []
            pressure_side_nodes_ids = []
            main_surface_nodes_ids = []
            trailing_edge_nodes_ids = []
        
        
        if mesh_tips:
            
            right_tip_nodes_ids = [node_id(i, 0) for i in range(i_max)] \
                + extra_wing_tip_node_ids[0:len(extra_wing_tip_node_ids)//2]

            left_tip_nodes_ids = [node_id(i, j_max) for i in range(i_max)] \
                + extra_wing_tip_node_ids[len(extra_wing_tip_node_ids)//2:]
                
            wing_tips_nodes_ids = right_tip_nodes_ids + left_tip_nodes_ids   
        
        else:
            extra_wing_tip_node_ids = []
            right_tip_nodes_ids = []
            left_tip_nodes_ids = []
            wing_tips_nodes_ids = []
        
        
        if mesh_wake:
            
            wake_nodes_ids = [wake_node_id(i, j)
                             for i in range(i_max_wake + 1)
                             for j in range(j_max + 1)]
            wake_lines = [[wake_node_id(i, j) for i in range(i_max_wake + 1)]
                          for j in range(j_max + 1)]
            wake_lines = np.array(wake_lines)
            
        else:
            wake_nodes_ids = []
            wake_lines = []
        
        
        body_nodes_ids = main_surface_nodes_ids + extra_wing_tip_node_ids 
        
        nodes_ids_dict = {
            "body": body_nodes_ids,
            "main surface" : main_surface_nodes_ids,
            "suction side": suction_side_nodes_ids,
            "pressure side": pressure_side_nodes_ids,
            "wing tips": wing_tips_nodes_ids,
            "right wing tip": right_tip_nodes_ids,
            "left wing tip": left_tip_nodes_ids,
            "trailing edge": trailing_edge_nodes_ids,
            "wake": wake_nodes_ids,
            "wake lines": wake_lines
        }
        
        nodes = [(node[0], node[1], node[2]) for node in nodes]
        
        if standard_mesh_format:
            return nodes, shells
        else:
            return nodes, shells, nodes_ids_dict

    def subdivide_spanwise_sections(self, div_num:int,
                                    interpolation_type="linear"):
        
        # interpolation_type="linear" or "cubic"
        
        new_Xsections = []
        spanwise_fractions = np.linspace(0, 1, div_num+2)
        
        for Xsection_0, Xsection_1 in zip(self.wingXsection_list[0:-1],
                                          self.wingXsection_list[1:]):
            
            new_Xsections.append(Xsection_0)
    
            for s in spanwise_fractions[1:-1]:    
                new_Xsections.append(
                    WingCrossSection.blend_WingCrossSections(
                        Xsection_0, Xsection_1, s, interpolation_type
                    )
                )
        
        new_Xsections.append(Xsection_1)
        
        self.wingXsection_list = new_Xsections

    def compute_platform_chars(self):
        mac = 0
        S = 0
        for i in range(len(self.wingXsection_list)-1):
            Xsection = self.wingXsection_list[i]
            Xsection_next = self.wingXsection_list[i+1]
            dy = abs(Xsection_next.r_leadingEdge.y - Xsection.r_leadingEdge.y)
            chord = (Xsection.chord + Xsection_next.chord)/2
            S = chord * dy + S
            mac = chord**2 * dy + mac
        self.S_planform = S
        self.MAC = mac/S

class Wing(BWB):
           
    def __init__(self,name:str, root_airfoil:Airfoil, tip_airfoil:Airfoil,
                 half_span:float, sweep_angle:float = 0,
                 dihedral_angle:float = 0, twist_angle:float = 0):
        
        
        self.half_span = half_span
        self.sweep_angle = sweep_angle
        self.dihedral_angle = dihedral_angle
        self.twist_angle = twist_angle
        self.root_wingXsection = None
        self.left_tip_wingXsection = None
        self.right_tip_wingXsection = None        
        super().__init__(name, [self.right_tip_wingXsection,
                                self.root_wingXsection,
                                self.left_tip_wingXsection])
        self.set_wingXsection_list(root_airfoil, tip_airfoil)
        
    def set_root_wingXsection(self, root_airfoil):
        chord = root_airfoil.chord
        x_coords = root_airfoil.x_coords/chord
        y_coords = root_airfoil.y_coords/chord
        name = root_airfoil.name
        airfoil = Airfoil(name, chordLength=1,
                          x_coords=x_coords, y_coords=y_coords)
        r_leadingEdge = Vector((0, 0, 0))
        self.root_wingXsection = WingCrossSection(r_leadingEdge = r_leadingEdge,
                                                  chord = chord,
                                                  twist=0,
                                                  airfoil = airfoil)
    
    def set_tip_wingXsections(self, tip_airfoil):
        chord = tip_airfoil.chord
        x_coords = tip_airfoil.x_coords
        y_coords = tip_airfoil.y_coords
        name = tip_airfoil.name
        airfoil = Airfoil(name, chordLength=1,
                          x_coords=x_coords, y_coords=y_coords)
        
        # left wing tip
        x_le, y_le, z_le = 0, self.half_span, 0
        x_le = self.sweep(x_le, self.half_span, np.deg2rad(self.sweep_angle))
        y_le, z_le = self.rotate(y_le, z_le, (0, 0),
                                 np.deg2rad(self.dihedral_angle))
        
        r_leadingEdge = Vector((x_le, y_le, z_le))
        self.left_tip_wingXsection = WingCrossSection(r_leadingEdge,
                                                      chord=chord,
                                                      twist = self.twist_angle,
                                                      airfoil=airfoil)
        
        # right wing tip
        y_le, z_le = -self.half_span, 0
        y_le, z_le = self.rotate(y_le, z_le, (0, 0),
                                 np.deg2rad(-self.dihedral_angle))
        
        r_leadingEdge = Vector((x_le, y_le, z_le))
        self.right_tip_wingXsection = WingCrossSection(r_leadingEdge,
                                                       chord=chord,
                                                       twist = self.twist_angle,
                                                       airfoil=airfoil)
    
    def set_wingXsection_list(self, root_airfoil, tip_airfoil):
        self.set_root_wingXsection(root_airfoil)
        self.set_tip_wingXsections(tip_airfoil)
        self.wingXsection_list = [self.right_tip_wingXsection,
                                  self.root_wingXsection,
                                  self.left_tip_wingXsection]
    
    @staticmethod
    def rotate(x, y, rotate_location, rotate_angle):
        
        # this fucntion rotates a point(x, y) about z-axis
        # to rotate a point (y, z) about x-axis:
        # y, z = rotate(y, z, (yo, zo), angle)
        # to rotate a point (x, z) about y-axis:
        # z, x = rotate(z, x, (zo, xo), angle)
        
        x_rot = rotate_location[0]
        y_rot = rotate_location[1]
        angle = rotate_angle
        x = (
            (x - x_rot) * np.cos(angle)
            + (y - y_rot) * np.sin(angle)
            + x_rot 
        )
        
        y = (
            -(x - x_rot) * np.sin(angle)
            + (y - y_rot) * np.cos(angle)
            + y_rot                 
        )
        
        return x, y   
    
    @staticmethod
    def sweep(x, span_location, sweep_angle):
        x = x + abs(span_location) * np.tan(sweep_angle)
        return x 


def BWB_reader(filePath, fileName, scale = 0.001):
    fileName = filePath + fileName

    with open(fileName) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")              
        data_list = [row for row in csv_reader]

    data_dict = {
        "airfoil name" : [data_list[i][0] for i in range(len(data_list))],
        
        "chord" : [
            scale * float(data_list[i][1]) for i in range(len(data_list))
        ],
        
        "leading edge coords" : [
            (
                scale * float(data_list[i][2]), scale * float(data_list[i][3]),scale * float(data_list[i][4])
            )
            for i in range(len(data_list))
        ],
        
        "twist" : [
            float(data_list[i][5]) for i in range(len(data_list))
        ]
    }
    
    return data_dict
    

if __name__=="__main__":
    from mesh_class import PanelMesh, PanelAeroMesh
    
    
    # bwb = BWB(name = "first try",
    #           wingXsection_list = [
    #               WingCrossSection(
    #               r_leadingEdge=Vector((0, -1, 0)),
    #               chord=1,
    #               twist=0,
    #               airfoil=Airfoil(name="naca0018_new")),
                  
    #               WingCrossSection(
    #               r_leadingEdge=Vector((0, -0.5, 0)),
    #               chord=1,
    #               twist=0,
    #               airfoil=Airfoil(name="naca0018_new")),
                  
    #               WingCrossSection(
    #               r_leadingEdge=Vector((0, 0.5, 0)),
    #               chord=1,
    #               twist=0,
    #               airfoil=Airfoil(name="naca0018_new")),
                  
    #               WingCrossSection(
    #               r_leadingEdge=Vector((0, 1, 0)),
    #               chord=1,
    #               twist=0,
    #               airfoil=Airfoil(name="naca0018_new"))
    #                                     ]
    #           )
    
    # bwb.subdivide_spanwise_sections(3)
    
    # nodes, shells = bwb.mesh_body(5, 1)
    
    # bwb_mesh = PanelMesh(nodes, shells)
    # bwb_mesh.plot_mesh_inertial_frame(elevation=-150,azimuth=-120)
    
    
    # # test wing class and meshing
    
    # wing = Wing(name="random",
    #             root_airfoil=Airfoil(name="naca0012_new", chord_length=1),
    #             tip_airfoil=Airfoil(name="naca0012_new", chord_length=1),
    #             half_span=1, sweep_angle=0, dihedral_angle=0, twist_angle=0)
    
    # nodes, shells, nodes_id = wing.mesh_body(
    #     5, 1, mesh_wake=True, wake_resolution=3, standard_mesh_format=False,
    #     shellType="quads")
    
    
    # wing_mesh = PanelAeroMesh(nodes, shells, nodes_id)
        
    # wing_mesh.plot_mesh_inertial_frame(elevation=-150, azimuth=-120,
    #                                    plot_wake=True)
    
    
    # RX3 concept
    
    # data_dict = BWB_reader(filePath="BWB/BWB concept/" , fileName= "BWB_X_sections_info")

    # # Change Airfoil's class, class atribute
    # Airfoil.filePath = "BWB/BWB concept/Airfoils/"

    # RX3 = BWB(
    #     name="RX3",
    #     wingXsection_list=[
    #         WingCrossSection(
    #             Vector(data_dict["leading edge coords"][id]),
    #             chord=data_dict["chord"][id],
    #             twist=data_dict["twist"][id],
    #             airfoil=Airfoil(
    #                 name=data_dict["airfoil name"][id]
    #             )
    #         )
            
    #         for id in range(len(data_dict["airfoil name"]))
                
    #     ]
    # )
    
    # RX3.subdivide_spanwise_sections(1, interpolation_type="linear")
    
    # nodes, shells, nodes_ids = RX3.mesh_body(
    #     ChordWise_resolution=20,
    #     SpanWise_resolution=1,
    #     SpanWise_spacing="uniform",
    #     shellType="quads",
    #     mesh_main_surface=True,
    #     mesh_tips=True,
    #     mesh_wake=True,
    #     wake_resolution=1,
    #     planar_wake=True,
    #     V_fs=Vector((1, 0, 0)),
    #     standard_mesh_format=False
    # )

    # rx3_mesh = PanelAeroMesh(nodes, shells, nodes_ids)

    # rx3_mesh.plot_mesh_bodyfixed_frame(
    # elevation=-120, azimuth=-150, plot_wake=False
    # )
    
    
    # BWB40_Sweep40deg
    filePath = "BWB/BWB40_Sweep40deg/"
    fileName = "BWB40_Sweep40deg_21_X_sections_info"
    data_dict = BWB_reader(filePath, fileName, scale=1)
    
    # Change Airfoil's class, class atribute
    Airfoil.filePath = "BWB/BWB40_Sweep40deg/Airfoils/"
    
    BWB40 = BWB(
            name="RX3",
            wingXsection_list=[
                WingCrossSection(
                    Vector(data_dict["leading edge coords"][id]),
                    chord=data_dict["chord"][id],
                    twist=data_dict["twist"][id],
                    airfoil=Airfoil(
                        name=data_dict["airfoil name"][id]
                    )
                )
                
                for id in range(len(data_dict["airfoil name"]))        
            ]
    )

    BWB40.subdivide_spanwise_sections(1, interpolation_type="linear")
    
    nodes, shells, nodes_ids = BWB40.mesh_body(
        ChordWise_resolution=20,
        SpanWise_resolution=1,
        SpanWise_spacing="uniform",
        shellType="quads",
        mesh_main_surface=True,
        mesh_tips=True,
        mesh_wake=True,
        wake_resolution=1,
        wake_length_in_chords=10,
        planar_wake=True,
        V_fs=Vector((1, 0, 0)),
        standard_mesh_format=False
    )

    BWB40_mesh = PanelAeroMesh(nodes, shells, nodes_ids)

    BWB40_mesh.plot_mesh_bodyfixed_frame(
        elevation=-120, azimuth=-150, plot_wake=False
    )
    
    BWB40_mesh.plot_mesh_bodyfixed_frame(
        elevation=-120, azimuth=-150, plot_wake=True
    )

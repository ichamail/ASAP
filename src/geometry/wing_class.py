import numpy as np
from scipy.stats import beta
from src.myMath import Vector
from .airfoil_class import Airfoil
from collections import deque


class Wing:
    
    def __init__(self, root_airfoil:Airfoil, tip_airfoil:Airfoil,
                 halfSpan:float, sweepAngle:float = 0, dihedralAngle:float = 0,
                 twistAngle:float = 0):
        
        self.halfSpan = halfSpan
        self.sweepAngle = np.deg2rad(sweepAngle)
        self.dihedralAngle = np.deg2rad(dihedralAngle)
        self.twistAngle = np.deg2rad(twistAngle)
        self.root_airfoil = root_airfoil
        self.tip_airfoil = tip_airfoil
        self.taperRatio = self.tip_airfoil.chord/self.root_airfoil.chord
        self.RefArea = 0
        self.setReferenceArea()
    
    def setReferenceArea(self):
                
        # self.RefArea = 2 * (
        #     (self.root_airfoil.chord + self.tip_airfoil.chord) * self.halfSpan
        # )/2
        
        # precise way
        xTipTwisted, _ = self.twistAirfoilSection(
            x_coords=self.tip_airfoil.coords[:-1, 0],
            z_coords=- self.tip_airfoil.coords[:-1, 1], # z_coords
            chord=self.tip_airfoil.chord,
            twist_angle=self.twistAngle
        )
                
        C_t = xTipTwisted.max() - xTipTwisted.min() 
        half_span = self.halfSpan - self.halfSpan*(1-np.cos(self.dihedralAngle))
        
        C_r = self.root_airfoil.chord
        
        self.RefArea = 2 * ((C_r + C_t) * half_span)/2
        
    def changeChordWiseResolution(self, num_x_points, spacing="cosine"):
        self.root_airfoil.repanel(num_x_points+1, spacing)
        self.tip_airfoil.repanel(num_x_points+1, spacing) 
        
    def meshSurface(
        self, numOfChordWiseFaces:int, numOfSpanWiseFaces:int, faceType:str="quadrilateral", chordWiseSpacing:str="cosine",
        spanWiseSpacing:str="uniform", mesh_MainSurface=True, mesh_WingTips:bool=True):
        
        
        if spanWiseSpacing == "uniform":
            space = np.linspace
        elif spanWiseSpacing == "cosine":
            space = self.cosspace
        elif spanWiseSpacing == "beta distribution":
            space = lambda start, end, steps: self.DenserAtBoundaries(
                start, end,steps, -0.15
            )
        elif spanWiseSpacing == "denser at wingtips":
            space = lambda start, end, steps: self.DenserAtWingTips(
                start, end, steps, factor=5
            )
        elif spanWiseSpacing == "denser at root":
            space = lambda start, end, steps: self.DenserAtWingRoot(
                start, end, steps, factor=5
            )
        elif spanWiseSpacing == "logarithmic":
            space = self.logspace
                
        self.changeChordWiseResolution(
            numOfChordWiseFaces, spacing=chordWiseSpacing
        )
        
        # wing's coordinate system:
        # origin: root airfoils leading edge
        # x-axis: chord wise direction from leading edge to trailing edge
        # y-axis: spanwise direction form root to tip
        # z-axis: from top to bottom
                
        x_root = self.root_airfoil.coords[:-1, 0]/self.root_airfoil.chord
        z_root = self.root_airfoil.coords[:-1, 1]/self.root_airfoil.chord
        x_tip = self.tip_airfoil.coords[:-1, 0]/self.tip_airfoil.chord
        z_tip = self.tip_airfoil.coords[:-1, 1]/self.tip_airfoil.chord
                                
        span_percentage = np.array(
            [
                *space(-self.halfSpan, 0, numOfSpanWiseFaces + 1),
                *space(0, self.halfSpan, numOfSpanWiseFaces + 1)[1:]
            ]
        )/self.halfSpan
        
        nx = len(x_root)
        ny = len(span_percentage)
        
        X = np.zeros((nx, ny))
        Y = np.zeros((nx, ny))
        Z = np.zeros((nx, ny))
        
        for j in range(ny):
                        
            C_y = self.interpolation(
                self.root_airfoil.chord, self.tip_airfoil.chord,
                abs(span_percentage[j])
            )
                        
            x, y, z = self.moveAirfoilSection(
                x_coords = self.interpolation(
                    x_root, x_tip, abs(span_percentage[j])
                ) * C_y
                ,
                y_coords = span_percentage[j] * self.halfSpan * np.ones(nx)
                ,
                z_coords = self.interpolation(
                    z_root, z_tip, abs(span_percentage[j])
                ) * C_y
                ,
                span_location = span_percentage[j] * self.halfSpan
                ,
                chord = C_y
                ,
                twist_angle = self.interpolation(
                    0, self.twistAngle, abs(span_percentage[j])
                )
                ,
                sweep_angle = self.sweepAngle
                ,
                dihedral_angle = self.dihedralAngle
            )          
            
            X[:, j] = x 
            Y[:, j] = y
            Z[:, j] = z 
        
        node = np.array(
            [
                [X[i][j], Y[i][j], Z[i][j]]
                for i in range(nx)
                for j in range(ny)
            ]
        )
        
        
        def node_id(chordWiseIndex, spanWiseIndex):
            return (spanWiseIndex + chordWiseIndex*ny)%(nx*ny)
         
        def addFace(*node_ids, reverse_order=False):
            # node_id_list should be in counter clock wise order
            node_ids = list(node_ids)
            
            if reverse_order:
                node_ids.reverse()
                
            if len(node_ids) == 4:
                
                if faceType=="quadrilateral":
                    
                    face.append(node_ids)
                    
                elif faceType=="triangular":
                    
                    face.append([node_ids[0], node_ids[1], node_ids[2], -1])
                    face.append([node_ids[2], node_ids[3], node_ids[0], -1])
            
            elif len(node_ids) == 3:
                
                face.append(node_ids + [-1])
        
        
        j_max = ny-1
        i_max = nx
        
        face = []
        
        if mesh_MainSurface:
            
            if faceType=="quadrilateral":
                                        
                # pressure and suction sides
                for i in range(i_max):
                    for j in range(j_max):
                        
                        addFace(
                            node_id(i, j),
                            node_id(i, j+1),
                            node_id(i+1, j+1),
                            node_id(i+1, j)
                        )
                                
            elif faceType=="triangular":
                
                # symmetrical meshing 
                for i in range(i_max):
                    
                    for j in range(j_max):
                        
                        # suction side
                        if i < i_max//2:
                            
                            # right side
                            if j < j_max//2:
                                addFace(    
                                    node_id(i, j),
                                    node_id(i, j+1),
                                    node_id(i+1, j+1),
                                    node_id(i+1, j)
                                )
                            
                            # left side    
                            else:
                                addFace(    
                                    node_id(i+1, j),
                                    node_id(i, j),
                                    node_id(i, j+1),
                                    node_id(i+1, j+1)
                                )

                        # pressure side
                        else:
                            
                            # right side
                            if j < j_max//2:
                                addFace(    
                                    node_id(i+1, j),
                                    node_id(i, j),
                                    node_id(i, j+1),
                                    node_id(i+1, j+1)
                                )
                            
                            # left side
                            else:
                                addFace(    
                                    node_id(i, j),
                                    node_id(i, j+1),
                                    node_id(i+1, j+1),
                                    node_id(i+1, j)
                                )
        
        
        if mesh_WingTips:
            
            id = len(node) -1 # last node id
            
            for j in [0, j_max]:
                # j=0 root or left tip
                # j-j_max tip or right tip
                
                if j==0:
                    addWingTipFace = addFace
                elif j==j_max:
                    
                    def addWingTipFace(*node_ids):
                        node_ids = deque(node_ids)
                        node_ids.rotate(-1)
                        return addFace(*node_ids, reverse_order=True)
                
                # root or left tip   
                id = id + 1
                   
                # trailing edge   
                i = 0
                     
                r_node = (
                    Vector(*node[node_id(i+1, j)])
                    + Vector(*node[node_id(i_max-i-1, j)])
                ) / 2
                
                node = np.vstack(
                    (
                        node, [r_node.x, r_node.y, r_node.z] 
                    )
                )
                
                addWingTipFace(
                    node_id(i, j),
                    node_id(i+1, j),
                    id
                )
                
                addWingTipFace(
                    node_id(i, j),
                    id,
                    node_id(i_max - i - 1, j)
                    
                )
                
                for i in range(1, i_max//2 - 1):
                    id = id+1
                    
                    r_node = (
                        Vector(*node[node_id(i+1, j)])
                        + Vector(*node[node_id(i_max-i-1, j)])
                    ) / 2
                    
                    node = np.vstack(
                        (
                            node, [r_node.x, r_node.y, r_node.z] 
                        )
                    )
                                
                    addWingTipFace(
                        node_id(i, j),
                        node_id(i+1, j),
                        id,
                        id-1
                    )
                    
                    addWingTipFace(
                        id-1,
                        id,
                        node_id(i_max - i - 1, j),
                        node_id(i_max-i, j)
                    )
                    
                # leading edge
                i = i+1           
                addWingTipFace(
                    node_id(i, j),
                    node_id(i+1, j),
                    id
                )
                
                addWingTipFace(
                    id,
                    node_id(i+1, j),  # node_id(i_max - i - 1, j),
                    node_id(i+2, j)  # node_id(i_max - i, j) 
                )
          
        return node, np.array(face)
        
    @staticmethod
    def sweepAirfoilSection(x_coords, span_location, sweep_angle):
        x_coords = x_coords + abs(span_location) * np.tan(sweep_angle)
        return x_coords
          
    @staticmethod
    def rotate(x_coords, y_coords, rotate_location, rotate_angle):

        x_c, y_c = rotate_location[0], rotate_location[1]
        angle = rotate_angle
        
        x = (
            (x_coords - x_c) * np.cos(angle)
            + (y_coords - y_c) * np.sin(angle)
            + x_c 
        )
        
        y = (
            -(x_coords - x_c) * np.sin(angle)
            + (y_coords - y_c) * np.cos(angle)
            + y_c                 
        )
        
        return x, y
    
    def twistAirfoilSection(self, x_coords, z_coords, chord, twist_angle):
        
        z_coords, x_coords = self.rotate(z_coords, x_coords,
                                        (0, 0.25*chord), twist_angle)
        
        return x_coords, z_coords
       
    def rotateAirfoilSection(self, y_coords, z_coords, span_location,
                               dihedral_angle):
        
        if span_location < 0:
            dihedral_angle = - dihedral_angle
        
        # y_coords, z_coords = self.rotate(y_coords, z_coords,
        #                                  (0, 0), dihedral_angle)
        
        # or xflr5's style to avoid  airfoil sections near root to intersect
        root_gamma = 0
        tip_gamma = dihedral_angle
        span_percentage = abs(span_location)/self.halfSpan
        
        section_gamma = self.interpolation(root_gamma, tip_gamma, span_percentage)
                
        y_coords, z_coords = self.rotate(y_coords, z_coords,
                                         (span_location, 0), section_gamma)
        
        z_coords = z_coords - span_location * np.tan(tip_gamma)
        
        return y_coords, z_coords
           
    def moveAirfoilSection(self, x_coords, y_coords, z_coords, span_location,
                             chord, twist_angle, sweep_angle, dihedral_angle):
        # angles in rads
                    
        x_coords, z_coords = self.twistAirfoilSection(x_coords,
                                                        z_coords,
                                                        chord, twist_angle)
        
        x_coords = self.sweepAirfoilSection(x_coords, span_location,
                                              sweep_angle)
        
        y_coords, z_coords = self.rotateAirfoilSection(y_coords, z_coords,
                                                         span_location, dihedral_angle)        
        return x_coords, y_coords, z_coords
    
    @staticmethod
    def cubic_function(x):
        
        """
        cubic function used for airfoil interpolations
        
        from
        May-Fun Liou, Hyoungjin Kim, ByungJoon Lee, and Meng-Sing Liou
        "Aerodynamic Design of the Hybrid Wing Body Propulsion-Airframe Integration"
        """
        return x**2 * (3 - 2*x)

    @staticmethod
    def interpolation(rootValue, tipValue, spanPercentage, type:str="linear"):
        
        """
        linear interpolation:
        x = {(y_tip - y)*x_root + (y - y_root)*x_tip}/(y_tip - y_root)
        x = {y_tip*x_root + y*(x_tip - x_root)}/y_tip
        x = x_root + y/y_tip * (x_tip - x_root)
        """
        
        if type == "linear":
            section_value = rootValue + (tipValue - rootValue) * spanPercentage
        elif type == "cubic":
            section_value = (
                rootValue + (tipValue 
                            - rootValue) * Wing.cubic_function(spanPercentage)
                )
        
        return section_value

    @staticmethod
    def DenserAtBoundaries(start, end, num_points, alpha):
        '''
        Beta distribution
        
        Cumulative distribution function of beta distribution
        
        alpha exists in (-oo, +oo)
        when alpha is 1 => evenly spaced
        when alpha < 1 denser at boundaries
        when alpha > 1 denser at midpoint
        '''
        x = np.linspace(0, 1, num_points)
        a = b = 2-alpha
        return start + beta.cdf(x, a, b) * (end-start)

    @staticmethod
    def DenserAtWingTips(start, end, num_points, factor=5):

        # produces same results as flow5 inv_exp spacing when factor=5

        x = np.linspace(0, 1, num_points)
        if start < 0  and end <=0:
            a = 1 + 0.1*factor
            b = 1
        elif start>= 0 and end>0 :
            a = 1
            b = 1 + 0.1*factor

        else:
            raise Exception("start, end >=0 or start, end <=0")


        return start + beta.cdf(x, a, b) * (end-start)

    @staticmethod
    def DenserAtWingRoot(start, end, num_points, factor=5):
        x = np.linspace(0, 1, num_points)
        if start < 0 :
            a = 1
            b = 1 + 0.1*factor
        else:
            a = 1 + 0.1*factor
            b = 1

        return start + beta.cdf(x, a, b) * (end-start)  

    @staticmethod
    def logspace(start, end, num_points):
        x = np.linspace(1, np.e, num_points)
        x = np.log(x)
        if  0 <= start < end:
            x = x*(end-start)
            x = x + start
        elif start < end <= 0:
            x = x * (start-end)
            x = x + end
            x = np.flip(x)
        return x

    @staticmethod
    def cosspace(start, end, num_points):
        mean = (start+end)/2
        amp = (end-start)/2
        return mean + amp * np.cos(np.linspace(np.pi, 0, num_points))


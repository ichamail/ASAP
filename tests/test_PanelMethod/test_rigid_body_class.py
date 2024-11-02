from src.PanelMethod import RigidAerodynamicBody, PanelMesh
from src.geometry import Wing, Airfoil

def test_rigidAerodynamicBodyWake():
    
    numOfChordWiseFaces=5
    numOfSpanWiseFaces=2
    
    rigidAerodynamicBody = RigidAerodynamicBody(
        surfaceMesh=PanelMesh(
            *Wing(
                root_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                tip_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                halfSpan=0.5 
            ).meshSurface(
                numOfChordWiseFaces,
                numOfSpanWiseFaces,
                faceType="quadrilateral",
                chordWiseSpacing="cosine",
                spanWiseSpacing="uniform",
                mesh_MainSurface=True,
                mesh_WingTips=True,
            )
        ),
        name="flying wing"
    )
    
    rigidAerodynamicBody.setTrailingEdge(
        trailingEdgeVerticesIDs=[i for i in range(2*numOfSpanWiseFaces+1)]
    )
  
    rigidAerodynamicBody.setWake(
        length=1,
        numOfWakeFaces=4,
        faceType="Quads",
        isWakeFree=False
    )
    
    rigidAerodynamicBody.display()   
    
    rigidAerodynamicBody.set_BodyFixedFrame_origin(2, 0, 0)
    rigidAerodynamicBody.set_BodyFixedFrame_orientation(0, 45, 0)
    
    rigidAerodynamicBody.display()   

def test_rigidAerodynamicBodyWakeMotion():
    
    numOfChordWiseFaces=5
    numOfSpanWiseFaces=2
    
    rigidAerodynamicBody = RigidAerodynamicBody(
        surfaceMesh=PanelMesh(
            *Wing(
                root_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                tip_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                halfSpan=0.5 
            ).meshSurface(
                numOfChordWiseFaces,
                numOfSpanWiseFaces,
                faceType="quadrilateral",
                chordWiseSpacing="cosine",
                spanWiseSpacing="uniform",
                mesh_MainSurface=True,
                mesh_WingTips=True,
            )
        ),
        name="flying wing"
    )
    
    rigidAerodynamicBody.setTrailingEdge(
        trailingEdgeVerticesIDs=[i for i in range(2*numOfSpanWiseFaces+1)]
    )
  
    rigidAerodynamicBody.setWake(
        length=0,
        numOfWakeFaces=0,
        faceType="Quads",
        isWakeFree=True
    )
    
    rigidAerodynamicBody.display()   
    
    rigidAerodynamicBody.set_BodyFixedFrame_origin(2, 0, 0)
    rigidAerodynamicBody.set_BodyFixedFrame_orientation(0, 45, 0)
    rigidAerodynamicBody.set_BodyFixedFrame_origin_velocity(-1, 0, 0)
    
    rigidAerodynamicBody.display()
        
    rigidAerodynamicBody.move_BodyFixedFrame(dt=0.5)
    
    rigidAerodynamicBody.display()
    
    rigidAerodynamicBody.move_BodyFixedFrame(dt=0.5)
       
    rigidAerodynamicBody.display()
    
def test_rigidAerodynamicBodyWakeMotion2():
    
    numOfChordWiseFaces=5
    numOfSpanWiseFaces=2
    
    rigidAerodynamicBody = RigidAerodynamicBody(
        surfaceMesh=PanelMesh(
            *Wing(
                root_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                tip_airfoil=Airfoil(name="naca0012 sharp", chordLength=0.5),
                halfSpan=0.5 
            ).meshSurface(
                numOfChordWiseFaces,
                numOfSpanWiseFaces,
                faceType="quadrilateral",
                chordWiseSpacing="cosine",
                spanWiseSpacing="uniform",
                mesh_MainSurface=True,
                mesh_WingTips=True,
            )
        ),
        name="flying wing"
    )
    
    rigidAerodynamicBody.setTrailingEdge(
        trailingEdgeVerticesIDs=[i for i in range(2*numOfSpanWiseFaces+1)]
    )
  
    rigidAerodynamicBody.setWake(
        length=1,
        numOfWakeFaces=4,
        faceType="Quads",
        isWakeFree=True
    )
    
    
    rigidAerodynamicBody.display(bodyFixedFrame=False)
    
    rigidAerodynamicBody.set_BodyFixedFrame_origin(xo=1, yo=1, zo=1)
    
    rigidAerodynamicBody.set_BodyFixedFrame_origin_velocity(
        Vo_x=-0.5, Vo_y=0, Vo_z=0
    )
    
    rigidAerodynamicBody.set_BodyFixedFrame_angular_velocity(20, 0, 0)
    
    # rigidAerodynamicBody.wake.Vinf = Vector(0.5, 0, 0)
    
    
    rigidAerodynamicBody.display(bodyFixedFrame=False)
    
    for i in range(90):
        rigidAerodynamicBody.move_BodyFixedFrame(dt=0.05)
    
    rigidAerodynamicBody.display(bodyFixedFrame=False)    

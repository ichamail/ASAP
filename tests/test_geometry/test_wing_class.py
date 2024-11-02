from src.geometry import Wing, Airfoil
from src.PanelMethod import Mesh, RigidBody

def plotWing():
    
    RigidBody(
        surfaceMesh = Mesh(
            *Wing(
                root_airfoil=Airfoil(name="naca0012 sharp"),
                tip_airfoil=Airfoil(name="naca0012 sharp"),
                halfSpan=1,
                sweepAngle=0,
                dihedralAngle=0,
                twistAngle=0
            ).meshSurface(
                numOfChordWiseFaces=10,
                numOfSpanWiseFaces=5,
                faceType="quadrilateral",
                chordWiseSpacing="cosine",
                spanWiseSpacing="uniform",
                mesh_MainSurface=True,
                mesh_WingTips=True,
            )
        ),
        name="wing"
    ).display(bodyFixedFrame=False)

from airfoil_class import Airfoil
from wing_class import Wing
from mesh_class import PanelMesh
from rigid_body_class import RigidAerodynamicBody
from panel_method_class import PanelMethod
from time import perf_counter


numOfChordWiseFaces=5
numOfSpanWiseFaces=5


panelMethod = PanelMethod(
    rigidBody=RigidAerodynamicBody(
        surfaceMesh=PanelMesh(
            *Wing(
                root_airfoil=Airfoil(
                    name="naca0012 sharp",
                    chordLength=1
                ),
                tip_airfoil=Airfoil(
                    name="naca0012 sharp",
                    chordLength=1
                ),
                halfSpan=1,
                sweepAngle=0,
                dihedralAngle=0,
                twistAngle=0
            ).meshSurface(
                numOfChordWiseFaces,
                numOfSpanWiseFaces,
                faceType="quadrilateral",
                chordWiseSpacing="cosine",
                spanWiseSpacing="uniform",
                mesh_MainSurface=True,
                mesh_WingTips=True
            )
        ),
        name="wing"
    )
)

panelMethod.setTrailingEdge(
    trailingEdgeVerticesIDs=[i for i in range(numOfSpanWiseFaces*2 + 1)]
)

panelMethod.locateSheddingFaces()


panelMethod.rigidBody.surface.setVSAeroAdjacencyMatrix(
    numOfChordWiseFaces, numOfSpanWiseFaces
)

panelMethod.setWake(
    length=0,
    numOfWakeFaces=0,
    faceType="Quads",
    isWakeFree=False
)


panelMethod.rigidBody.display(
    elevation=30,
    azimuth=-60,
    bodyFixedFrame=True
)


panelMethod.setVfs(
    angleOfAttack=10,
    sideSlipAngle=0,
    magnitude=1
)

# panelMethod.setVinf(
#     angleOfAttack=10,
#     sideSlipAngle=0,
#     magnitude=1
# )


panelMethod.rigidBody.display(
    elevation=30,
    azimuth=-60,
    bodyFixedFrame=True
)



panelMethod.solve(steadyState=False, iters=5)

panelMethod.rigidBody.display()

panelMethod.displaySurfacePressureCoefficientContour()


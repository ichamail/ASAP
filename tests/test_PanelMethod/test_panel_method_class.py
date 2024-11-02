from src.geometry import Sphere, Airfoil, Wing
from src.PanelMethod import PanelMesh, RigidBody, BoundaryElementMethod, PanelMethod, RigidAerodynamicBody

def testBoundaryElementMethod():
    
    bem = BoundaryElementMethod(
        rigidBody=RigidBody(
            PanelMesh(
                *Sphere(
                    center=(0, 0, 0),
                    radius=1
                ).meshIcoSurface(3)
            ),
            name="Unit Sphere"
        )
    )
        
    bem.setVfs(
        angleOfAttack=0,
        sideSlipAngle=0,
        magnitude=1
    )
    
    # bem.setVinf(
    #     angleOfAttack=0,
    #     sideSlipAngle=0,
    #     magnitude=1
    # )
    
    bem.solve()
    
    bem.displaySurfacePressureCoefficientContour()
      
def testPanelMethod(steadyState=True, iters=0):
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



    panelMethod.solve(steadyState, iters)

    panelMethod.rigidBody.display()

    panelMethod.displaySurfacePressureCoefficientContour()

def test_PanelMethod_SteadyState_rigidWake():
    
    testPanelMethod(steadyState=True, iters=0)

def test_PanelMethod_SteadyState_iterativeWake():
    
    testPanelMethod(steadyState=True, iters=5)
    
def test_PanelMethod_Unsteady():
    
    testPanelMethod(steadyState=False, iters=5)

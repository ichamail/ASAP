from tests.test_geometry import *
from tests.test_PanelMethod import *

def main():
    
    ####################### Test geometry Package ##########################
    
    # plotCircle()
    # plotAirfoil()
        
    # plotCube()
    # plotSphere()
    
    # plotWing()
    
    # test_BWB40_Sweep40deg()
    # test_BWB_spanwise_subdivision()
    
    
    
    ####################### Test PanelMethod Package ##########################
    
    # testQuadPanel()
    # testTriPanel()
    # testHexagonPanel()
    # testSourcePanel()
    # testDoubletPanel()
    
    # test_rigidAerodynamicBodyWake()
    # test_rigidAerodynamicBodyWakeMotion()
    # test_rigidAerodynamicBodyWakeMotion2()
    
    testBoundaryElementMethod()
    # test_PanelMethod_SteadyState_rigidWake()
    # test_PanelMethod_SteadyState_iterativeWake()
    # test_PanelMethod_Unsteady()
    
    # simulateSphere()
    
    # simulateWing()
    
    # simulateBWB()
    
    pass


if __name__=="__main__":
    
    main()

from tests.test_geometry import *
from tests.test_PanelMethod import *
from tests.test_myMath import *
from tests.test_numbaSpeedUp import *

def main():
    
    test()
           
    return 0
 

def test():
    
    ####################### Test myMath Package ##########################
    
    # test_Vector()
    
    
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
    test_PanelMethod_SteadyState_rigidWake()
    # test_PanelMethod_SteadyState_iterativeWake()
    # test_PanelMethod_Unsteady()
    
    # simulateSphere()
    
    # simulateWing()
    
    # simulateBWB()
    
    
    ###################### Test numbaSpeedUp Package #######################
    
    # print_numba_float(1.2)
    # test_numbaPanel()
    # test_numbaPanelMethod()
    
    pass

   
if __name__=="__main__":
    
    main()

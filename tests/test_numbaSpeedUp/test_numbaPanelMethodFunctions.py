from ..test_PanelMethod import testBoundaryElementMethod, test_PanelMethod_SteadyState_rigidWake, test_PanelMethod_SteadyState_iterativeWake, test_PanelMethod_Unsteady


def test_numbaPanelMethod():
    
    testBoundaryElementMethod()
    test_PanelMethod_SteadyState_rigidWake()
    test_PanelMethod_SteadyState_iterativeWake()
    test_PanelMethod_Unsteady()

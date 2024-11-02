from .test_panel_class import testDoubletPanel, testHexagonPanel, testQuadPanel, testSingularity, testSourcePanel, testSurfacePanel, testTriPanel

from .test_rigid_body_class import test_rigidAerodynamicBodyWake, test_rigidAerodynamicBodyWakeMotion, test_rigidAerodynamicBodyWakeMotion2

from .test_panel_method_class import testBoundaryElementMethod, test_PanelMethod_SteadyState_rigidWake, test_PanelMethod_SteadyState_iterativeWake, test_PanelMethod_Unsteady

from .sphere_case import simulateSphere

from .wing_case import simulateWing

from .bwb_case import simulateBWB


__all__= (
    "testDoubletPanel", "testHexagonPanel", "testQuadPanel", "testSingularity", "testSourcePanel", "testSurfacePanel", "testTriPanel",
    "test_rigidAerodynamicBodyWake", "test_rigidAerodynamicBodyWakeMotion","test_rigidAerodynamicBodyWakeMotion2",
    "testBoundaryElementMethod", "test_PanelMethod_SteadyState_rigidWake", "test_PanelMethod_SteadyState_iterativeWake", "test_PanelMethod_Unsteady",
    "simulateSphere",
    "simulateWing",
    "simulateBWB"
)

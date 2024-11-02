from src.geometry import BWB, BWB_reader, Airfoil, WingCrossSection
from src.PanelMethod import PanelMesh, RigidBody
from src.myMath import Vector


def plot_BWB40_Sweep40deg():
    
    # BWB40_Sweep40deg
    data_dict = BWB_reader(
        filePath="data/BWB/BWB40_Sweep40deg/",
        fileName="BWB40_Sweep40deg_21_X_sections_info",
        scale=1
    )

    # Change Airfoil's class, class atribute
    Airfoil.filePath = "data/BWB/BWB40_Sweep40deg/Airfoils/"

    RigidBody(
        surfaceMesh = PanelMesh(
            *BWB(
                name="RX3",
                wingXsection_list=[
                    WingCrossSection(
                        Vector(*data_dict["leading edge coords"][id]),
                        chord=data_dict["chord"][id],
                        twist=data_dict["twist"][id],
                        airfoil=Airfoil(
                            name=data_dict["airfoil name"][id]
                        )
                    )
                    
                    for id in range(len(data_dict["airfoil name"]))        
                ]
            ).mesh_body(
                ChordWise_resolution=20,
                SpanWise_resolution=1,
                SpanWise_spacing="uniform",
                shellType="quads",
                mesh_main_surface=True,
                mesh_tips=True
            )
        ),
        name="RX3"
    ).display(
        bodyFixedFrame=False
    )

def test_BWB_spanwise_subdivision():
    data_dict = BWB_reader(
        filePath="data/BWB/BWB40_Sweep40deg/",
        fileName="BWB40_Sweep40deg_21_X_sections_info",
        scale=1
    )

    # Change Airfoil's class, class atribute
    Airfoil.filePath = "data/BWB/BWB40_Sweep40deg/Airfoils/"
    
    # BWB40_Sweep40deg
    data_dict = BWB_reader(
        filePath="data/BWB/BWB40_Sweep40deg/",
        fileName="BWB40_Sweep40deg_21_X_sections_info",
        scale=1
    )
    
    rx3 = BWB(
        name="RX3",
        wingXsection_list=[
            WingCrossSection(
                Vector(*data_dict["leading edge coords"][id]),
                chord=data_dict["chord"][id],
                twist=data_dict["twist"][id],
                airfoil=Airfoil(
                    name=data_dict["airfoil name"][id]
                )
            )
            
            for id in range(len(data_dict["airfoil name"]))        
        ]
    )
    
    rx3.subdivide_spanwise_sections(1, interpolation_type="linear")
    
    RigidBody(
        *rx3.mesh_body(
            ChordWise_resolution=20,
            SpanWise_resolution=1,
            SpanWise_spacing="uniform",
            shellType="quads",
            mesh_main_surface=True,
            mesh_tips=True
        ),
        name=rx3.name
    ).display(bodyFixedFrame=False)

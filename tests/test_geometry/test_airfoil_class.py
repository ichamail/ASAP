from src.geometry import Airfoil, Circle

def plotAirfoil():
    
    airfoil =  Airfoil(
        name="naca0012 sharp",
        chordLength=2,
    )
    airfoil.repanel(
        numOfPointsPerSide=5+1,
        spacing="denser at leading edge"
    )
    airfoil.plot()
    

def plotCircle():

    Circle(
        name="circle", center=(0, 0), radius=1, num_points=11
    ).plot()

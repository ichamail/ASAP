import numpy as np
from matplotlib import pyplot as plt
from src.utilities import set_axes_equal
from src.myMath import Vector
from src.PanelMethod import Panel, TriPanel, QuadPanel, Source, Doublet, SurfacePanel


def testHexagonPanel():
       
    # Polygon panel
    def hexagon(center=(2, 2, 2), r=1):
        x0, y0, z0 = center
        numS = 6 # number of sides
        # theta0 = (360/(numB-1))/2
        theta0 = 0
        theta = np.linspace(0, 360, numS+1)
        theta = theta + theta0
        theta = theta*(np.pi/180)
        x = x0 + r* np.cos(theta)
        y = y0 + r* np.sin(theta)
        return np.array([[x[i], y[i], z0] for i in range(numS)])
    
    panel = Panel(hexagon())
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(4, 4, 4), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(4, 4, 4), panelFixedFrame=True
    )
    
def testTriPanel():
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2]
    ])
    
    panel = TriPanel(vertices, CCW=True)
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=True
    )
    
def testQuadPanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = QuadPanel(vertices, CCW=True)
    
    panel.display(panelFixedFrame=False)
    panel.display(panelFixedFrame=True)
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=False
    )
    
    panel.display_point_P(
        r_p=Vector(1, 2, 3), panelFixedFrame=True
    )
    
    pass

def testSourcePanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = Source(vertices, CCW=True)
    panel.sigma = 1
    
    
    # Velocity Potential - Scalar field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        100
    )
        
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        100
    )
    
    X, Z = np.meshgrid(x, z, indexing="ij")
    phi = np.zeros_like(X)
    
    nx, nz = X.shape
    for i in range(nx):
        for j in range(nz):
            phi[i][j] = panel.unitStrength_inducedVelocityPotential(
                r_p=Vector(X[i][j], panel.r_centroid.y, Z[i][j])
            )
    CS = plt.contour(X, Z, phi)
    plt.clabel(CS, inline = 1, fmt ='% 1.2f', fontsize = 8)
    plt.show()
    
    
    # # Velocity Potential - Scalar field and Velocity Vector Field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        10
    )
    y = np.linspace(
        panel.r_centroid.y - 1 * panel.charLength,
        panel.r_centroid.y + 1 * panel.charLength,
        10
    )
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        10
    )
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    phi = np.zeros_like(X)
    v_x, v_y, v_z = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    nx, ny, nz = X.shape
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz): 
                r_p = Vector(X[i][j][k], Y[i][j][k], Z[i][j][k])
                phi[i][j][k] = panel.inducedVelocityPotential(r_p) 
                v = panel.inducedVelocity(r_p)
                v_x[i][j][k], v_y[i][j][k], v_z[i][j][k] = v.x, v.y, v.z
    
    ax, fig = panel.plot(panelFixedFrame=False)
    scatter = ax.scatter(X, Y, Z, c=phi, cmap="jet")
    fig.colorbar(scatter, ax=ax)
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    
    
    ax, fig = panel.plot(panelFixedFrame=False)
    ax.quiver(X, Y, Z, v_x, v_y, v_z, normalize=True, color='m')
    
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    pass

def testDoubletPanel():
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = Doublet(vertices, CCW=True)
    panel.mu = 1
    
    # Velocity Potential - Scalar field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        100
    )
        
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        100
    )
    
    X, Z = np.meshgrid(x, z, indexing="ij")
    phi = np.zeros_like(X)
    
    nx, nz = X.shape
    for i in range(nx):
        for j in range(nz):
            phi[i][j] = panel.unitStrength_inducedVelocityPotential(
                r_p=Vector(X[i][j], panel.r_centroid.y, Z[i][j])
            )
    CS = plt.contour(X, Z, phi)
    plt.clabel(CS, inline = 1, fmt ='% 1.2f', fontsize = 8)
    plt.show()
    
    
    # # Velocity Potential - Scalar field and Velocity Vector Field
    x = np.linspace(
        panel.r_centroid.x - 1 * panel.charLength,
        panel.r_centroid.x + 1 * panel.charLength,
        10
    )
    y = np.linspace(
        panel.r_centroid.y - 1 * panel.charLength,
        panel.r_centroid.y + 1 * panel.charLength,
        10
    )
    z = np.linspace(
        panel.r_centroid.z - 1 * panel.charLength,
        panel.r_centroid.z + 1 * panel.charLength,
        10
    )
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    phi = np.zeros_like(X)
    v_x, v_y, v_z = np.zeros_like(X), np.zeros_like(X), np.zeros_like(X)
    nx, ny, nz = X.shape
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz): 
                r_p = Vector(X[i][j][k], Y[i][j][k], Z[i][j][k])
                phi[i][j][k] = panel.inducedVelocityPotential(r_p) 
                v = panel.inducedVelocity(r_p)
                v_x[i][j][k], v_y[i][j][k], v_z[i][j][k] = v.x, v.y, v.z
    
    ax, fig = panel.plot(panelFixedFrame=False)
    scatter = ax.scatter(X, Y, Z, c=phi, cmap="jet")
    fig.colorbar(scatter, ax=ax)
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()
    
    
    
    ax, fig = panel.plot(panelFixedFrame=False)
    ax.quiver(X, Y, Z, v_x, v_y, v_z, normalize=True, color='m')
    
    ax.set_xlim3d(- panel.charLength + x[0], x[-1] + panel.charLength)
    ax.set_ylim3d(- panel.charLength + y[0], y[-1] + panel.charLength)
    ax.set_zlim3d(- panel.charLength + z[0], z[-1] + panel.charLength)
    set_axes_equal(ax)
    plt.show()

def testSingularity(panel:Source|Doublet, r_p:Vector):
    
    panel.display_point_P(r_p, panelFixedFrame=False)   
    phi = panel.unitStrength_inducedVelocityPotential(r_p)
    V = panel.inducedVelocity(r_p)
    
    print("phi = ", phi)
    print("V = ", V)
 
def testSurfacePanel(r_p):
    
    vertices = np.array([
        [1, 1, 2],
        [3, 1, 2],
        [3, 3, 2],
        [1, 3, 2]
    ])
    
    panel = SurfacePanel(vertices, CCW=True)
    
    panel.sigma, panel.mu = 1, 1
    
    panel.display_point_P(r_p)
    
    phi_src, phi_dblt = panel.unitStrength_inducedVelocityPotential(r_p)
    v = panel.inducedVelocity(r_p)
    
    print("phi_src = ", phi_src)
    print("phi_dblt = ", phi_dblt)
    print("v_src + v_dblt = ", v)

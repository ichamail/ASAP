import numpy as np
from matplotlib import pyplot as plt, cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.patches as mpatches



def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

def plot_Cp_SurfaceContours(panel_list, elevation=30, azimuth=-60):
  shells = []
  vert_coords = []
  for panel in panel_list:
      shell=[]
      for r_vertex in panel.r_vertex:
          shell.append((r_vertex.x, r_vertex.y, r_vertex.z))
          vert_coords.append([r_vertex.x, r_vertex.y, r_vertex.z])
      shells.append(shell)
  
  Cp = [panel.Cp for panel in panel_list]
  Cp_norm = [(float(Cp_i)-min(Cp))/(max(Cp)-min(Cp)) for Cp_i in Cp]
  facecolor = plt.cm.coolwarm(Cp_norm)
  
  fig = plt.figure()
  
  ax = plt.axes(projection='3d')
  poly3 = Poly3DCollection(shells, facecolor=facecolor)
  ax.add_collection(poly3)
  ax.view_init(elevation, azimuth)
  ax.set_xlabel('X')
  ax.set_ylabel('Y')
  ax.set_zlabel('Z')

  vert_coords = np.array(vert_coords)
  x, y, z = vert_coords[:, 0], vert_coords[:, 1], vert_coords[:, 2]
  ax.set_xlim3d(x.min(), x.max())
  ax.set_ylim3d(y.min(), y.max())
  ax.set_zlim3d(z.min(), z.max())
  set_axes_equal(ax)
  
  
  m = cm.ScalarMappable(cmap=cm.coolwarm)
  m.set_array([min(Cp),max(Cp)])
  m.set_clim(vmin=min(Cp),vmax=max(Cp))
  Cbar = fig.colorbar(m, ax=ax)
  # Cbar.set_ticks([round(x,2) for x in np.linspace(min(Cp), max(Cp), 6)])
  Cbar.set_ticks(np.linspace(min(Cp), max(Cp), 6))
  Cbar.set_ticklabels([str(round(x,2)) for x in np.linspace(min(Cp), max(Cp), 6)])
  Cbar.set_label("Cp", rotation=0)
  
  plt.show()

def move_view(event, ax):
  ax.autoscale(enable=False, axis='both') 
  koef = 10
  zkoef = (ax.get_zbound()[0] - ax.get_zbound()[1]) / koef
  xkoef = (ax.get_xbound()[0] - ax.get_xbound()[1]) / koef
  ykoef = (ax.get_ybound()[0] - ax.get_ybound()[1]) / koef
  ## Map an motion to keyboard shortcuts
  if event.key == "ctrl+down":
    ax.set_ybound(ax.get_ybound()[0] + xkoef, ax.get_ybound()[1] + xkoef)
  if event.key == "ctrl+up":
    ax.set_ybound(ax.get_ybound()[0] - xkoef, ax.get_ybound()[1] - xkoef)
  if event.key == "ctrl+right":
    ax.set_xbound(ax.get_xbound()[0] + ykoef, ax.get_xbound()[1] + ykoef)
  if event.key == "ctrl+left":
    ax.set_xbound(ax.get_xbound()[0] - ykoef, ax.get_xbound()[1] - ykoef)
      
  if event.key == "alt+down":
    ax.set_zbound(ax.get_zbound()[0] - zkoef, ax.get_zbound()[1] - zkoef)
  if event.key == "alt+up":
    ax.set_zbound(ax.get_zbound()[0] + zkoef, ax.get_zbound()[1] + zkoef)
  
  # zoom option
  if event.key == "+":
    ax.set_xbound(ax.get_xbound()[0]*0.90, ax.get_xbound()[1]*0.90)
    ax.set_ybound(ax.get_ybound()[0]*0.90, ax.get_ybound()[1]*0.90)
    ax.set_zbound(ax.get_zbound()[0]*0.90, ax.get_zbound()[1]*0.90)
  if event.key == "-":
    ax.set_xbound(ax.get_xbound()[0]*1.10, ax.get_xbound()[1]*1.10)
    ax.set_ybound(ax.get_ybound()[0]*1.10, ax.get_ybound()[1]*1.10)
    ax.set_zbound(ax.get_zbound()[0]*1.10, ax.get_zbound()[1]*1.10)
  
  # Rotational movement
  elev=ax.elev
  azim=ax.azim
  if event.key == "8":
    elev+=10
  if event.key == "2":
    elev-=10  
  if event.key == "4":
    azim+=10
  if event.key == "6":
    azim-=10

  ax.view_init(elev= elev, azim = azim)

  # print which ever variable you want 

  ax.figure.canvas.draw()
  
  import matplotlib.patches as mpatches

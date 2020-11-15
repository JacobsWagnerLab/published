import numpy as np
import matplotlib.pyplot as plt

def plot_track(track):
    """ Plot the conformation of a simulated track viewed from different planes
    - track: N by 3 numpy array that represents the coordinates of the track
    """
    # Extract the positions
    x,y,z = track[:,0],track[:,1],track[:,2]
    # Plot XY plane
    plt.figure(figsize=(6,3))
    plt.plot(x,y,'r',linewidth=0.1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()
    # Plot XZ plane
    plt.figure(figsize=(6,3))
    plt.plot(x,z,'g',linewidth=0.1)
    plt.xlabel('x')
    plt.ylabel('z')
    plt.axis('equal')
    plt.show()
    # Plot YZ plane
    plt.figure(figsize=(6,3))
    plt.plot(y,z,'b',linewidth=0.1)
    plt.xlabel('y')
    plt.ylabel('z')
    plt.axis('equal')
    plt.show()

def rotate(v,phi,theta):
    """ Rotate a vector by phi and theta
    """
    Rx = np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])
    Rz = np.array([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1]])
    return np.matmul(Rz,np.matmul(Rx,v))

def check_angle(u,v):
    """ Check the angle between two vectors
    """
    return np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))

def genstep(b,n):
    """ Generate n steps with a step size of b and a random orientation
    """
    steps = []
    step0 = [0,0,b]
    phi0,phi1 = 0,np.pi
    theta0,theta1 = 0,2*np.pi
    thetas = np.random.random(n)*(theta1-theta0)+theta0
    phis = np.arccos(2*(np.random.random(n)*(((np.cos(phi1)+1)/2)-((np.cos(phi0)+1)/2))+((np.cos(phi0)+1)/2))-1)
    for phi,theta in zip(phis,thetas):
        step = rotate(step0,phi,theta)
        steps.append(step)
    return np.array(steps)

def calc_Rg(track):
    """ Calculate the Rg of a track
    - track: N by 3 numpy array that represents the coordinates of the track
    """
    x,y,z = track[:,0],track[:,1],track[:,2]
    r_mean = [np.mean(x),np.mean(y),np.mean(z)]
    rsq = (x - r_mean[0])**2 + (y - r_mean[1])**2 + (z - r_mean[2])**2
    return np.sqrt(np.sum(rsq)/len(rsq))

def is_in_bound(x,y,z,d=1,L=1):
    """ Check if a position (x,y,z) is within the bound of a spherocylinder
    """
    if abs(x) <= L / 2:
        return y ** 2 + z ** 2 <= d ** 2 / 4
    elif abs(x) > L / 2:
        return (abs(x) - L/2) ** 2 + y ** 2 + z ** 2 <= d ** 2 / 4
    else:
        return False

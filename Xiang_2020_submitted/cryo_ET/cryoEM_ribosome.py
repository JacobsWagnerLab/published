import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import time
from scipy.stats import linregress, multivariate_normal
from scipy.spatial import ConvexHull, convex_hull_plot_2d
from mpl_toolkits.mplot3d import Axes3D
from IPython.display import clear_output
import seaborn as sns
from matplotlib.path import Path
import copy
from scipy import ndimage
import pickle

def center_decorator(centroid=(0,0,0)):
    """ Decorator used to center a group of points
    """
    x0,y0,z0 = centroid
    def center_(res):
        res[:,0] = res[:,0] - np.mean(res[:,0]) + x0
        res[:,1] = res[:,1] - np.mean(res[:,1]) + y0
        res[:,2] = res[:,2] - np.mean(res[:,2]) + z0
        return res
    return center_

def rotate(i,j,res,theta=None):
    """ Rotate the points such that they are all in plane
    """
    if theta is None:
        lm = linregress(res[:,i],res[:,j])
        slope = lm.slope
        theta = -slope
    else:
        theta = theta / 180 * np.pi
    R = [[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]
    Rv = np.matmul(R,res[:,[i,j]].T).T
    res[:,i] = Rv[:,0]
    res[:,j] = Rv[:,1]
    return res

def load_position_file(filename):
    """ Load a position file which should contain 3 columns represent X, Y and Z coordinates
    """
    def process(filename):
        res = []
        with open(filename,"r") as f:
            for line in f.readlines():
                res.append(line.split())
        res = np.array(res,dtype='double')
        if res.shape[1] == 4:
            res = res[:,1:]
        elif res.shape[1] == 5:
            res = res[:,1:]
        return res

    def align(res):
        return rotate(1,2,rotate(0,2,res))

    return align(process(filename))

def calc_distance(x0,y0,z0,x1,y1,z1):
    """ Calculate the Euclidean distance
    """
    return np.sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)

def get_contour(pos,ix=None):
    """ Get the convex contour defined by points
    """
    if ix is None:
        points = pos[:,0:2]
    else:
        points = pos[ix, 0:2]
    print('total points %d' % len(points))
    hull = ConvexHull(points)
    contourx = points[hull.vertices,0]
    contourx = np.append(contourx,contourx[0])
    contoury = points[hull.vertices,1]
    contoury = np.append(contoury,contoury[0])
    return contourx, contoury

def make_plot(pos, bound=None, bin_=100, vlim=None, pad=0, blur=0,savename=None):
    """ Make a 2D histogram showing density of points
    """
    # Extract positions
    xs = pos[:,0]
    ys = pos[:,1]
    zs = pos[:,2]

    # Boundary of the histogram
    if bound is None:
        xmin,xmax = np.min(xs), np.max(xs)
        ymin,ymax = np.min(ys), np.max(ys)
    else:
        xmin,xmax,ymin,ymax = bound
        ix = np.logical_and(ys > ymin, ys < ymax)
        xs = xs[ix]
        ys = ys[ix]
        zs = zs[ix]
        if xmin == -1:
            xmin = xs.min()
        if xmax == -1:
            xmax = xs.max()
    rmin,rmax = min(xmin,ymin),max(xmax,ymax)
    rmin -= pad
    rmax += pad
    print("Image size = %d a.u." % (rmax-rmin))
    print("Bin size = %.2f a.u." % ((rmax-rmin)/bin_))

    # Start making a plot
    f = plt.figure(figsize=(2,2))
    ax = f.add_subplot(111)
    # Construct a histogram
    [counts,_,_,_] = ax.hist2d(xs,ys,bins=bin_,density=True,range=((rmin,rmax),(rmin,rmax)))
    plt.close()

    counts[np.isnan(counts)] = 0
    x = np.flipud(counts.T)

    # Using a Gaussian filter to blur the image
    if blur != 0:
        x = ndimage.gaussian_filter(x,blur)
        x[np.isnan(x)] = 0

    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)
    if vlim is None:
        plt.imshow(x,cmap='jet')
    else:
        plt.imshow(x,cmap='jet',vmin=vlim[0],vmax=vlim[1])
    plt.colorbar(ax=ax,fraction=0.046, pad=0.04)
    plt.axis("equal")
    plt.axis('off')
    plt.tight_layout()
    if savename is not None:
        plt.savefig(savename)
    plt.show()

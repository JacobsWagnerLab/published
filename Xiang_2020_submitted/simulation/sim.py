import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt

def plot_track(track):
    x,y,z = track[:,0],track[:,1],track[:,2]
    plt.figure(figsize=(6,3))
    plt.plot(x,y,'r',linewidth=0.1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()
    plt.figure(figsize=(6,3))
    plt.plot(x,z,'g',linewidth=0.1)
    plt.xlabel('x')
    plt.ylabel('z')
    plt.axis('equal')
    plt.show()
    plt.figure(figsize=(6,3))
    plt.plot(y,z,'b',linewidth=0.1)
    plt.xlabel('y')
    plt.ylabel('z')
    plt.axis('equal')
    plt.show()

def rotate(v,phi,theta):
    Rx = np.array([[1,0,0],[0,np.cos(phi),np.sin(phi)],[0,-np.sin(phi),np.cos(phi)]])
    Rz = np.array([[np.cos(theta),np.sin(theta),0],[-np.sin(theta),np.cos(theta),0],[0,0,1]])
    return np.matmul(Rz,np.matmul(Rx,v))

def check_angle(u,v):
    return np.arccos(np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)))

def genstep(b,n):
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

def calc_R(track):
    x,y,z = track[:,0],track[:,1],track[:,2]
    return np.sqrt((x[-1]-x[0])**2+(y[-1]-y[0])**2+(z[-1]-z[0])**2)

def calc_Rg(track):
    x,y,z = track[:,0],track[:,1],track[:,2]
    r_mean = [np.mean(x),np.mean(y),np.mean(z)]
    rsq = (x - r_mean[0])**2 + (y - r_mean[1])**2 + (z - r_mean[2])**2
    return np.sqrt(np.sum(rsq)/len(rsq))

def calc_angle(track):
    angles = []
    for i in range(2,track.shape[0]):
        angles.append(check_angle(track[i,:]-track[i-1,:], track[i-1,:]-track[i-2,:]))
    return np.array(angles)

def is_in_bound(x,y,z,d=1,L=1):
    if abs(x) <= L / 2:
        return y ** 2 + z ** 2 <= d ** 2 / 4
    elif abs(x) > L / 2:
        return (abs(x) - L/2) ** 2 + y ** 2 + z ** 2 <= d ** 2 / 4
    else:
        return False

def _worker(task):
    b,n,rangle,seed,d,L,trial_limit = task
    np.random.seed(seed)
    return _sim_track(b,n,rangle,d,L,np.inf)

def _manager(worker,tasks):
    with mp.Pool(mp.cpu_count()) as pool:
        res = pool.map(worker,tasks)
    return res

def _sim_track(b,n,rangle,d,L,trial_limit):
    track = [np.array([0,0,b])]
    currPos = np.array([0,0,0])
    for _ in range(n-1):
        found = False
        trial = 0
        while not found:
            candidates = genstep(b,1)
            for candidate in candidates:
                x,y,z = currPos + candidate
                if is_in_bound(x,y,z,d,L) and trial <= trial_limit:
                    trial += 1
                    angle = check_angle(track[-1],candidate)
                    if angle >= rangle[0] and angle <= rangle[1]:
                        found = True
                        trial = 0
                        track.append(candidate)
                        currPos = currPos + candidate
                        break
                elif is_in_bound(x,y,z,d,L) and trial > trial_limit:
                    found = True
                    trial = 0
                    track.append(candidate)
                    currPos = currPos + candidate
                    break  
    track[0] = np.array([0,0,0])
    track = np.cumsum(track,axis=0)
    return track

def sim_bound(b,n,c,solvent,trial_limit=np.inf):
    vol = 4.6e6*650e3/6.02e23/1e-12/c
    d = (3./(np.pi*2))**(1/3)*vol**(1/3)
    L = 2.*d
    if solvent.lower() == "poor":
        rangle = (11./18*np.pi,np.pi)
    elif solvent.lower() == "ideal":
        rangle = (0,np.pi)
    elif solvent.lower() == "good":
        rangle = (0,np.pi/2)
    else:
        return
    return _sim_track(b,n,rangle,d,L,trial_limit)

def sim_free(b,n,solvent):
    if solvent.lower() == "poor":
        rangle = (11./18*np.pi,np.pi)
    elif solvent.lower() == "ideal":
        rangle = (0,np.pi)
    elif solvent.lower() == "good":
        rangle = (0,np.pi/2)
    else:
        return

    track = _sim_track(b,n,rangle,np.inf,np.inf,np.inf)
    if solvent.lower() == "poor":
        tmp = track
        minR = calc_R(track)
        if (minR <= 1):
            return track
        else:
            count = 0
            parallel = mp.cpu_count()
            while minR > 1 and count < 10:
                seed = np.random.randint(0,1e6,parallel)
                tasks = [task for task in zip([b]*parallel,[n]*parallel,[rangle]*parallel,seed,[np.inf]*parallel,[np.inf]*parallel,[np.inf]*parallel)]
                tracks = _manager(_worker,tasks)
                for track in tracks:
                    R = calc_R(track)
                    if (R <= 1):
                        return track
                    elif (R < minR):
                        minR = R
                        tmp = track
                count+=1
            return tmp
    else:
        return track

def sim_rg(b,n,solvent,repeat=1,trial_limit=np.inf):
    if solvent.lower() == "poor":
        rangle = (11./18*np.pi,np.pi)
        d = 1.7
        L = 1.7
        trial_limit=200
    elif solvent.lower() == "ideal":
        rangle = (0,np.pi)
        d = np.inf
        L = np.inf
    elif solvent.lower() == "good":
        rangle = (0,np.pi/2)
        d = np.inf
        L = np.inf
    else:
        return
    seed = np.random.randint(0,1e6,repeat)
    tasks = [task for task in zip([b]*repeat,[n]*repeat,[rangle]*repeat,seed,[d]*repeat,[L]*repeat,[trial_limit]*repeat)]
    tracks = _manager(_worker,tasks)
    Rgs = []
    for track in tracks:
        Rgs.append(calc_Rg(track/np.sqrt(2)))
    return np.array(Rgs)

def sim_scaling(ns,rangle,repeat):
    res = []
    for n in ns:
        seed = np.random.randint(0,1e6,repeat)
        tasks = [task for task in zip([1]*repeat,[n]*repeat,[rangle]*repeat,seed,[np.inf]*repeat,[np.inf]*repeat,[np.inf]*repeat)]
        tracks = _manager(_worker,tasks)
        Rgs = []
        for track in tracks:
            Rgs.append(calc_Rg(track))
        res.append([n,np.mean(Rgs)])
    res = np.array(res)
    res[:,1] /= res[0,1]
    return res


__author__ = 'Papai'

import math
import pymc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymc.Matplot import plot
from scipy.stats import gaussian_kde

from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity

def bw_kde(data=[],start=0.01,end=1.0,cv_size=20):
    grid = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(start,end,cv_size)},cv=cv_size)
    grid.fit(data[:, None])
    return grid.best_params_

def bw_kde_cved(data=[],start=0.01,end=1.0,cv_size=20):
    success = False

    while success == False:
        initial_bw = bw_kde(data,start,end,cv_size)['bandwidth']
        if initial_bw in [start,end]:
            end = initial_bw * 10.0
            start = initial_bw / 10.0
        else :
            success = True

    initial_bw = math.floor(initial_bw*1000.0)/1000.0
    start=max(initial_bw-0.05,0.001)
    end=initial_bw+0.05
    actual_bw = bw_kde(data,start=start,end=end,cv_size=cv_size)[ 'bandwidth' ]
    return actual_bw

########################
def get_KDE(measurements=[],start=0.01,end=1.0,cv_size=20,graph=True):
    vmin = measurements.min()
    vmax = measurements.max()

    #grid = GridSearchCV(KernelDensity(),{'bandwidth': np.linspace(start,end,cv_size)},cv=cv_size)
    #grid.fit(measurements[:, None])
    #bw1 = grid.best_params_['bandwidth']
    bw1 = bw_kde_cved(measurements,start,end,cv_size)

    pos = np.linspace(vmin-0.25,vmax+0.25,num=1500)
    kde = gaussian_kde(measurements)
    kde.evaluate(pos)
    bw2 = kde.silverman_factor()
    bw3 = kde.scotts_factor()

    modes = []
    plot_data = {}
    for bw in [bw1,bw2,bw3]:
        kde.set_bandwidth(bw)
        pdf = kde.evaluate(pos)

        df = pd.DataFrame(np.transpose([pos,pdf]),columns=['position','density'])
        max_density = df[ 'position' ][ df[ 'density' ].argmax() ]
        modes.append(max_density)

        print("bandwidth : %s"%bw, "log density:%s"%np.sum(pdf), "max desnity at:%s"%max_density)

        kde.set_bandwidth(bw)
        plot_data['KDE_%s'%bw] = {'position':pos,'density':pdf}

    if graph == True:
        fig, ax = plt.subplots()
        for key in plot_data:
            inner_data = plot_data[key]
            ax.plot(inner_data['position'], inner_data['density'],label='bw={0}'.format(key), linewidth=3, alpha=0.5)
            ax.hist(measurements, 50, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
            ax.set_xlim(vmin-1.0, vmax+1.0)
            #ax.set_xlim(0.2, 0.3)
            ax.legend(loc='upper left')

    return {'modes':modes,'plot_data':plot_data,'bandwidth':[bw1,bw2,bw3]}

############################
def get_Bayes(measurements=[],chunksize=5,Ndp=5,iter=50000,burn=5000):

    sc = pymc.Uniform('sc',0.1,2.0,value=0.24)
    tau = pymc.Uniform('tau',0.0,1.0,value=0.5)

    concinit = 1.0
    conclo = 0.1
    conchi = 10.0
    concentration = pymc.Uniform('concentration', lower=conclo, upper=conchi,value=concinit)

    # The stick-breaking construction: requires Ndp beta draws dependent on the
    # concentration, before the probability mass function is actually constructed.
    #betas = pymc.Beta('betas', alpha=1, beta=concentration, size=Ndp)
    betas = pymc.Beta('betas', alpha=1, beta=1, size=Ndp-1)

    @pymc.deterministic
    def pmf(betas=betas):
        "Construct a probability mass function for the truncated Dirichlet process"
        # prod = lambda x: np.exp(np.sum(np.log(x))) # Slow but more accurate(?)
        prod = np.prod
        value = map(lambda i,u: u * prod(1.0 - betas[:i]), enumerate(betas))
        value.append(1.0 - sum(value[:])) # force value to sum to 1
        return value

    # The cluster assignments: each data point's estimated cluster ID.
    # Remove idinit to allow clusterid to be randomly initialized:
    Ndata = len(measurements)
    idinit = np.zeros(Ndata, dtype=np.int64)
    clusterid = pymc.Categorical('clusterid', p=pmf, size=Ndata, value=idinit)

    @pymc.deterministic(name='clustermean')
    def clustermean(clusterid=clusterid,sc=sc,Ndp=Ndp):
        return sc*np.arange(1,Ndp+1)[clusterid]

    @pymc.deterministic(name='clusterprec')
    def clusterprec(clusterid=clusterid,sc=sc,tau=tau,Ndp=Ndp):
        return 1.0/(sc*sc*tau*tau*(np.arange(1,Ndp+1)[clusterid]))

    y=pymc.Normal('y',mu=clustermean,tau=clusterprec,observed=True,value=measurements)

    ## for predictive poeterior simulation
    @pymc.deterministic(name='y_sim')
    def y_sim(value=[0],sc=sc,tau=tau,clusterid=clusterid,Ndp=Ndp):
        n = np.arange(1,Ndp+1)[np.random.choice(clusterid)]
        return np.random.normal(loc=sc*n,scale=sc*tau*n)

    m = pymc.Model({"scale":sc,"tau":tau,"betas":betas,"clusterid":clusterid,"normal":y,"pred":y_sim})

    sc_samples=[]
    modes=[]
    simulations=[]

    for i in range(0,chunksize):
        mc = pymc.MCMC(m)
        mc.sample(iter=50000,burn=10000)
        plot(mc)

        sc_sample = mc.trace('sc')[:]
        sc_samples.append(sc_sample)

        simulation=mc.trace('y_sim')[:]
        simulations.append(simulation)

        plt.hist(measurements,50, fc='gray', histtype='stepfilled', alpha=0.3, normed=False)
        plt.hist(simulation, 30, fc='blue', histtype='stepfilled', alpha=0.3, normed=True)
        hist, edges = np.histogram(measurements, bins=100, range=[np.min(measurements) - 0.25, np.max(measurements) + 0.25])

        argm = hist.argmax()
        ( edges[argm] + edges[ argm+1 ] ) / 2
        modes.append(( edges[argm] + edges[ argm+1 ] ) / 2)

    if chunksize <= 1:
        gr = np.nan
    else:
        pymc.gelman_rubin(sc_samples)

    dic = {'gelman_rubin':gr,'modes':modes,'simulations':simulations,'sc_samples':sc_samples}
    return dic


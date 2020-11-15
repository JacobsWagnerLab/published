__author__ = 'Papai'

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import SJ_Utils as sjutils
from sklearn import mixture

data = pd.read_excel('X4gmm.xlsx')
keys = data.columns.tolist()

k=1
#for key 3 in data.columns:
key='day6z3'
measurements = data[key]
#measurements = measurements*10000;
measurements = measurements.dropna()

KDE = sjutils.get_KDE(measurements,graph=False)

if KDE['bandwidth'][0] < 0.002:
    mode = KDE['modes'][1]
else:
    mode = KDE['modes'][0]

outliers = measurements[ measurements > 2 * k * mode ]
#assert (0.0+len(outliers))/len(measurements) < 0.05
measurements=measurements[ measurements <= 2 * k * mode ]
#measurements

#plt.hist(measurements.values,bins=40)

g = mixture.GaussianMixture(n_components=k)
g.fit(measurements.values.reshape(-1,1))

vmin = measurements.min()
vmax = measurements.max()

pos = np.linspace(vmin-0.25, vmax+0.25, num=1500)

pdf = np.exp(g.score(pos.reshape(-1,1)))
pdf = np.exp(g.score_samples(pos.reshape(-1,1)))

#To show distribution with the gmm result
fig, ax = plt.subplots()
ax.hist(measurements.values, 15, fc='gray', histtype='stepfilled', alpha=0.3, normed=True)
ax.plot(pos, pdf, label='GMM', linewidth=3, alpha=0.5)
ax.set_xlim(vmin-1.0, vmax+1.0)
#ax.set_xlim(0.2, 0.3)
ax.legend(loc='upper left')
#for m in g.means_:
#    plt.axvline(m[a0])
[g.weights_, g.means_, g.covariances_]
#len(outliers)
plt.show()

# g.means_ contain the information used for SM normalization.


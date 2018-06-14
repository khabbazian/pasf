import time
import numpy as np
import matplotlib.pyplot as plt
import matlab.engine

eng = matlab.engine.start_matlab()

simOpts = eng.struct('noiseVar', 0.16, 'noiseCorrCoeff', 0,
        'sourceEnergy', 6.3)
data = eng.rotating_energy_sources(70, simOpts)


dirname   = 'RES-rn'
bfilename = dirname + '/res'
eng.mkdir( dirname )

myOpts = eng.struct( 'cmethod', 'phase',
	'boolParfor', 0,
	'boolUseSavedData', 0,
	'errorRate', 0.1,
	'bfilename', bfilename)


eng.tic;
Components, Clusters, ClusterInfo, SDFInfo = eng.pasf(data, 2, myOpts, nargout=4)
eng.toc;

Components = np.array(Components._data).reshape(Components.size, order='F')
print( Components.shape )
d1, d2, d3, d4 = Components.shape

for i in range(d3):
    fig = Components[:,:,i,0]
    for c in range(1, d4):
        fig = np.hstack((fig, Components[:,:,i,c]))

    plt.imshow(fig)
    plt.show(block=False)
    time.sleep(0.5)
    plt.close()


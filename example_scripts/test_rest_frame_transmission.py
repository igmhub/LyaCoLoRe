import numpy as np
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
def test_JAF(N_pixels,QSO_RSD=False):
    ### Load data
    path = '/project/projectdirs/desi/mocks/lya_forest/london/v1.1/*/*/transmission-*-*.fits'
    fi = glob.glob(path)
    fi = np.sort(np.random.choice(fi,replace=False,size=N_pixels))
    ### Get T
    for f in fi:
        plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')
        print(f)
        t=fits.open(f)
        if QSO_RSD==False:
            z=t[1].data['Z_noRSD']
        else:
            z=t[1].data['Z']
            lObs=t[2].data
            print(lObs)
            trans=t[3].data
            n_low=0
            n_high=0
        for i in range(trans.shape[1]):
            print(i,z[i])
            if z[i]>(3550./1215.67 - 1):
                if n_high<20:
                    plt.plot(lObs/(1.+z[i]),trans[:,i],color='red',alpha=0.1)
                    n_high += 1
            else:
                if n_low<20:
                    plt.plot(lObs/(1.+z[i]),trans[:,i],color='black',alpha=0.1)
                    n_low += 1
            print(n_low,n_high,'\n')
            if n_low>10 and n_high>10: break
        plt.ylabel('F')
        plt.xlabel('lRF')
        plt.grid()
        plt.xlim(1200.,1250.)
        plt.title(f)
        plt.savefig('test_JAF_output.pdf')
        plt.show()
    return

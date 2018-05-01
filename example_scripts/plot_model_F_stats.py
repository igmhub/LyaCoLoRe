import numpy as np
import matplotlib.pyplot as plt
import process_functions as functions

z_values = np.linspace(0.0,4.0,101)
l_hMpc = 0.25
Om = 0.3
results = []

for z in z_values:
    print('looking at z={:2.2f}'.format(z),end='\r')
    mean_F_model = functions.get_mean_F_model(z)
    sigma_dF_model = functions.get_sigma_dF_P1D(z,l_hMpc=l_hMpc,Om=Om)
    results += [(z,mean_F_model,sigma_dF_model)]
print('\nDone!')

dtype = [('z', '>f4'), ('mean_F', '>f4'), ('sigma_dF', '>f4')]
results = np.array(results,dtype=dtype)

plt.plot(results['z'],results['mean_F'],label='mean_F')
plt.plot(results['z'],results['sigma_dF'],label='sigma_dF')
plt.plot(results['z'],results['sigma_dF']/2,label='sigma_dF/2')
plt.legend()
plt.grid()
plt.savefig('model_F_stats.pdf')
plt.show()

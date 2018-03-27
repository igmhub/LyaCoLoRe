import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from picca import wedgize
import sys

# export.py --data cf_rmax20_n5.fits.gz --out cf_exp_rmax20_n5.fits.gz
# export.py --data cf_rmax80_n20_n2.fits.gz --out cf_exp_rmax80_n20_n2.fits.gz
# export.py --data cf_rmax140_n35_n2.fits.gz --out cf_exp_rmax140_n35_n2.fits.gz
# export.py --data cf_rmax160_n40_n2.fits.gz --out cf_exp_rmax160_n40_n2.fits.gz
# export.py --data cf_rmax100_n25_n10.fits.gz --out cf_exp_rmax100_n25_n10.fits.gz
# export.py --data cf_rmax160_n40_n10.fits.gz --out cf_exp_rmax160_n40_n10.fits.gz

filename='cf_exp.fits.gz'

if len(sys.argv)>1:
    filename=sys.argv[1]
     
wedgize_args = sys.argv[2]
if len(wedgize_args) != 7:
    print("Not enough arguments for wedgizing. 7 are required:")
    print("(mumin,mumax,rtmax,nrt,rpmax,nrp,nr,rmax)")

h = fits.open(filename)

data = h[1].data

rp=data['RP'][:]
rt=data['RT'][:]
z=data['Z'][:]
xi_grid=data['DA'][:]
cov_grid=data['CO'][:]

#b = wedgize.wedge(mumin=0.0, mumax=1.0,rtmax=20.0,nrt=5,rpmax=20.0,nrp=5,nr=20,rmax=40)
#b = wedgize.wedge(mumin=0.0, mumax=1.0,rtmax=40.0,nrt=10,rpmax=40.0,nrp=10,nr=40,rmax=80)
#b = wedgize.wedge(mumin=0.0, mumax=1.0,rtmax=80.0,nrt=20,rpmax=80.0,nrp=20,nr=80,rmax=160)
#b = wedgize.wedge(mumin=0.0, mumax=1.0,rtmax=100.0,nrt=25,rpmax=100.0,nrp=25,nr=100,rmax=200)
#b = wedgize.wedge(mumin=0.0, mumax=1.0,rtmax=140.0,nrt=35,rpmax=140.0,nrp=35,nr=150,rmax=300)

b = wedgize.wedge(mumin=wedgize_args[0], mumax=wedgize_args[1],rtmax=wedgize_args[2],nrt=wedgize_args[3],rpmax=wedgize_args[4],nrp=wedgize_args[5],nr=wedgize_args[6],rmax=wedgize_args[7])

r,xi_wed,cov_wed = b.wedge(xi_grid,cov_grid)

Nr = len(r)
err_wed = np.zeros(Nr)
for i in range(Nr):
    err_wed[i] = np.sqrt(cov_wed[i][i])
    print(i,err_wed[i],cov_wed[i][i])

cut = err_wed>0
plt.errorbar(r[cut],xi_wed[cut]*(r[cut]**2),yerr=err_wed[cut]*(r[cut]**2),fmt='o',label='from ~7 000 density skewers')
mean=0.0071539422
plt.plot(r[cut],mean*mean*(r[cut]**2),
        label='$\\bar \\delta^2 r^2$, $\\bar \\delta = $'+str(mean))
plt.legend(loc=4)
plt.axhline(y=0,color='gray',ls=':')
plt.xlabel('r [Mpc/h]')
plt.ylabel('r^2 xi(r)')
plt.title('density correlation function')
plt.savefig('xi_wedge.pdf')
plt.show()


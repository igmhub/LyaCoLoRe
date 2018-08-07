import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
import process_functions as pf

lya = 1215.67
IVAR_cutoff = 1150.0

def read_file(basedir,quantity,nside,pix):
    pix_100 = int(pix/100)
    dirname = basedir+'/'+str(pix_100)+'/'+str(pix)+'/'
    suffix = str(nside)+'-'+str(pix)+'.fits'

    if quantity == 'flux':
        filename = dirname+'/transmission-'+suffix
        transmission = fits.open(filename)
        data_rows = transmission[3].data.T
        z_qso = transmission[1].data['Z_noRSD']
        transmission.close()
        filename = dirname+'/picca-{}-'.format(quantity)+suffix
        picca = fits.open(filename)
    else:
        # file with delta for picca
        filename = dirname+'/picca-{}-'.format(quantity)+suffix
        #print('picca flux file',filename)
        picca = fits.open(filename)
        data_rows = picca[0].data.T
        z_qso = picca[3].data['Z']

    #ivar_rows = picca[1].data.T
    loglam = picca[2].data
    picca.close()
    ivar_rows = pf.make_IVAR_rows(IVAR_cutoff,z_qso,loglam)

    zs = (10**loglam)/lya - 1.0

    return zs,data_rows,ivar_rows

def get_means(data_rows,ivar_rows):
    print(data_rows.shape,ivar_rows.shape)
    N_skewers=data_rows.shape[0]
    N_cells=data_rows.shape[1]
    # add tiny value to not have 0/0
    if data_rows.shape[0] != ivar_rows.shape[0]:
        print('data and ivar not same shape')
        weights = np.zeros(data_rows.shape[1])
        mean = np.zeros(data_rows.shape[1])
        mean2 = np.zeros(data_rows.shape[1])
    else:
        weights=np.sum(ivar_rows,axis=0)+1.e-10
        mean=np.sum(data_rows*ivar_rows,axis=0)/weights
        mean2=np.sum((data_rows**2)*ivar_rows,axis=0)/weights
    return weights,mean,mean2

# main folder where the processed files are
basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/process_output_G_hZsmooth_4096_32_sr2.0_bm1_biasG18_picos_nside16/'
#basedir = '/global/cscratch1/sd/jfarr/LyaSkewers/CoLoRe_GAUSS/test/'
#basedir = '/Users/James/Projects/test_data/process_output_G_hZ_4096_32_sr2.0_bm1_nside16/'

if len(sys.argv)>1:
    quantity = sys.argv[1]
else:
    quantity = 'flux'
nside = 16
N_pixels = 3072
pixels = np.sort(np.random.choice(list(range(12*nside**2)),size=N_pixels))
#pixels = np.array([1064,1096,1127,1128,1159,1160,1191,1192,1193,1223,1224,1225,1254,1255,1256,1257,1286,1287,1288,1289])
pixels = np.array(list(range(N_pixels)))

# will combine pixel statistics
sum_mean=None
sum_weights=None
sum_var=None

for i,pix in enumerate(pixels):
    print('Looking at pixel {}: ({} out of {})'.format(pix,i+1,pixels.shape[0]))
    # read file for particular pixel
    zs,data_rows,ivar_rows=read_file(basedir,quantity,nside,pix)

    # compute statistics for pixel
    weights,mean,mean2 = get_means(data_rows,ivar_rows)

    # accumulate stats
    if sum_mean is None:
        sum_mean = mean*weights
        sum_mean2 = mean2*weights
        sum_weights = weights
    else:
        sum_mean += mean*weights
        sum_mean2 += mean2*weights
        sum_weights += weights

overall_mean = sum_mean/sum_weights
overall_mean2 = sum_mean2/sum_weights
overall_var = overall_mean2 - overall_mean**2
overall_sigma = np.sqrt(overall_var)

err = np.sqrt(overall_var/sum_weights)

if quantity == 'flux':
    identifier = 'F'
elif quantity == 'gaussian':
    identifier = 'G'
elif quantity == 'density':
    identifier = 'D'
else:
    identifier = ''

plt.figure(figsize=(12, 8), dpi= 80, facecolor='w', edgecolor='k')

max_QSO_z = 3.79
z_max = IVAR_cutoff*(1+max_QSO_z)/lya - 1
plt.xlim(2.0,z_max + 4 - max_QSO_z)

plt.plot(zs,overall_mean,label='mean_{}'.format(identifier))
if quantity != 'flux':
    plt.plot(zs,overall_sigma,label='sigma_{}'.format(identifier))

plt.title('{} stats: mean and sigma, rest frame cut at {}A'.format(quantity,IVAR_cutoff))
plt.xlabel('z')

predicted_data = fits.open('input_files/tune_small_scale_fluctuations.fits')
z_predicted = predicted_data[1].data['z']
if quantity == 'gaussian':
    sigma_G_predicted = predicted_data[1].data['sigma_G']
    plt.plot(z_predicted,np.zeros(z_predicted.shape),label='predicted mean')
    plt.plot(z_predicted,sigma_G_predicted,label='predicted sigma')
if quantity == 'gaussian-RSD':
    sigma_G_predicted = predicted_data[1].data['sigma_G']
    plt.plot(z_predicted,np.zeros(z_predicted.shape),label='predicted mean')
    plt.plot(z_predicted,sigma_G_predicted,label='predicted sigma')
if quantity == 'gaussian-noRSD':
    sigma_G_predicted = predicted_data[1].data['sigma_G']
    plt.plot(z_predicted,np.zeros(z_predicted.shape),label='predicted mean')
    plt.plot(z_predicted,sigma_G_predicted,label='predicted sigma')
elif quantity == 'flux':
    mean_F_predicted = predicted_data[1].data['mean_F']
    sigma_dF_predicted = predicted_data[1].data['sigma_dF']

    sigma_F_predicted = sigma_dF_predicted*mean_F_predicted

    plt.plot(z_predicted,mean_F_predicted,label='predicted mean_F')
    #plt.plot(z_predicted,sigma_F_predicted,label='predicted sigma_F')

    overall_sigma_dF = overall_sigma/overall_mean

    plt.plot(zs,overall_sigma_dF,label='sigma_dF')
    plt.plot(z_predicted,sigma_dF_predicted,label='predicted sigma_dF')

    """
    mean_F_model = predicted_data[1].data['mean_F_needed']
    sigma_dF_model = predicted_data[1].data['sigma_dF_needed']
    plt.plot(z_predicted,mean_F_model,label='model mean')
    plt.plot(z_predicted,sigma_dF_model,label='model sigma')
    """
#plt.ylim(-0.05,0.05)
plt.legend()
plt.grid()
plt.savefig('mean_var_{}_cut{}.pdf'.format(quantity,IVAR_cutoff))
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import scipy.interpolate

curdir = os.path.dirname(__file__)

def set_to_nan(vcorr, value=0):
    # Sets value to NaN
    vcorr[np.where(vcorr == value)] = np.nan
    return vcorr

#%%

for resdir, partitle in zip(['toftpetersen_pars', 'yiu_pars'], ['Old parameters', 'New parameters']):
    # Loads McPhase data. The .xyt file has the "phase diagram" information
    # The .sps file has the induced moment components.
    # For each (H,T) are 25 numbers. The first is the temperature,
    # next are 4 sets of (Sx, Lx, Sy, Ly, Sz, Lz) moments for each site
    xyt_data = np.loadtxt(os.path.join(curdir, resdir, 'results', 'mcphas.xyt'), 
                          converters={n: lambda x: float(x.decode().split('p')[0]) for n in range(12)})
    xx = np.unique(xyt_data[:,0])
    yy = np.unique(xyt_data[:,1])
    shp = (len(xx), len(yy))
    sps_data = np.loadtxt(os.path.join(curdir, resdir, 'results', 'mcphas.sps'), usecols=0)
    sps_data = np.reshape(sps_data, (int(len(sps_data)/13), 13))
    Cx = np.abs(sps_data[:,1] + sps_data[:,10] - sps_data[:, 7] - sps_data[:, 4]) / 2
    Cy = np.abs(sps_data[:,2] + sps_data[:,11] - sps_data[:, 8] - sps_data[:, 5]) / 2
    Cz = np.abs(sps_data[:,3] + sps_data[:,12] - sps_data[:, 9] - sps_data[:, 6]) / 2
    Az = np.abs(sps_data[:,3] - sps_data[:,12] - sps_data[:, 9] + sps_data[:, 6]) / 2
    Ax = np.abs(sps_data[:,1] - sps_data[:,10] - sps_data[:, 7] + sps_data[:, 4]) / 2
    Ay = np.abs(sps_data[:,2] - sps_data[:,11] - sps_data[:, 8] + sps_data[:, 5]) / 2
    Fx = np.abs(sps_data[:,1] + sps_data[:,10] + sps_data[:, 7] + sps_data[:, 4]) / 2
    Fy = np.abs(sps_data[:,2] + sps_data[:,11] + sps_data[:, 8] + sps_data[:, 5]) / 2
    Fz = np.abs(sps_data[:,3] + sps_data[:,12] + sps_data[:, 9] + sps_data[:, 6]) / 2
    allOPs = (['Ax', 'Ay', 'Az', 'Cx', 'Cy', 'Cz', 'Fx', 'Fy', 'Fz'], 
              [set_to_nan(np.reshape(OP, shp)) for OP in [Ax, Ay, Az, Cx, Cy, Cz, Fx, Fy, Fz]])
    # Plots the order parameters
    fg, axs = plt.subplots(3,3, figsize=(10, 10))
    for ii, OP in enumerate(zip(*allOPs)):
        ax = axs[int(ii / 3), ii % 3]
        pcm = ax.pcolormesh(xx, yy, OP[1].T, cmap='viridis', vmin=0, vmax=3)
        ax.set_xlim([0,80])
        ax.set_ylim([0,120])
        if (int(ii / 3) == 2):
            ax.set_xlabel('Temperature (K)')
        if (ii % 3 == 0):
            ax.set_ylabel('Applied Field (T)')
        ax.set_title(f'${OP[0][0]}_{OP[0][1]}$')
        cb = fg.colorbar(pcm, ax=ax)
        if (ii % 3 == 2):
            cb.set_label('Order parameter ($\mu_B$)')
    fg.suptitle(f'{partitle}, spin-only SIA')
    fg.show()
    fg.savefig(os.path.join(curdir, f'ops_{resdir}.pdf'))
    # Plot the "phase diagram" (Cx - Cy)
    # Use a Diverging colormap... e.g.
    divmaps = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
               'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
    divmaps = ['Spectral']#, 'seismic']
    # Plot the "phase diagram" (Cx - Cy)
    dat = np.nansum(np.dstack((allOPs[1][3], -allOPs[1][4])), 2).T
    for mp in divmaps:
        fg, ax = plt.subplots()
        pcm = ax.pcolormesh(xx, yy, dat, cmap=mp, vmin=-3, vmax=3)
        ax.set_xlim([0,80])
        ax.set_ylim([0,120])
        ax.set_xlabel('Temperature (K)')
        ax.set_ylabel('Applied Field (T)')
        cb = fg.colorbar(pcm, ax=ax)
        cb.set_label('Induced Moment (uB)')
        fg.suptitle(f'{partitle}, spin-only SIA')
        fg.show()
        fg.savefig(os.path.join(curdir, f'pd_{resdir}.pdf'))

import h5py
import matplotlib.pyplot as plt
import numpy as np

from pathlib import Path

fnames = sorted(list(Path('data/20240329_C_10').glob('*.h5')))

regions = ['frontal', 'temporal', 'occipital', 'wm']

num_frms = len(fnames)
num_regions = len(regions)

for frm, fname in enumerate(fnames):
    with h5py.File(fname, 'r') as f:
        if frm == 0:
            num_real = f['C_recon_frontal'].shape[0]    
            C_recon = np.zeros((num_regions, num_frms, num_real))
            C_true = np.zeros((num_regions, num_frms))

        for j, region in enumerate(regions):
            C_true[j, frm] = f[f"C_true_{region}"][()]
            C_recon[j, frm, :] = f[f"C_recon_{region}"][:]
            t_min = f['t_min'][:]
        C_P = f['C_P'][:]
        C_TOT = f['C_TOT'][:]
        fbv = f['fbv'][()]

fig, ax = plt.subplots(3, 1, figsize = (8,8), tight_layout=True, sharex=True)
ax[0].plot(t_min, C_P, '.-', label = 'C_P')
ax[0].plot(t_min, C_TOT, '.-', label = 'C_TOT')
ax[1].plot(t_min, C_P/C_TOT, '.-', label = 'C_P/C_TOT')
for j in range(num_regions):
    for i in range(num_real):
        ax[2].plot(t_min, C_recon[j,:,i], color = plt.cm.tab10(j), lw = 0.1)
    ax[2].plot(t_min, C_true[j,:], '.', color = plt.cm.tab10(j), label = regions[j])
ax[0].set_ylabel('C_P and C_TOT')
ax[1].set_ylabel('C_P / C_TOT')
ax[2].set_ylabel('C_PET = (1-fbv) * C_T + fbv* C_TOT')
ax[2].set_xlabel('t (min)')
for axx in ax.ravel():
    axx.legend()
    axx.grid(ls=":")
fig.show()

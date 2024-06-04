from __future__ import annotations

import nibabel as nib
import pymirc.viewer as pv
import array_api_compat.numpy as np
import h5py
import matplotlib.pyplot as plt
import argparse
import array_api_compat.numpy as xp
from pathlib import Path

def cmdlineparse(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--phantom_file", default="tacs_C.h5")
    args = parser.parse_args(args)
    return args
def main(args=None):
    args = cmdlineparse(args)
    tac_file = args.phantom_file
    # %%
    dev = "cpu"
    fbv = 0.05
    #tac_file = "tacs_AD.h5"
    
    # %%
    # load the segmented brain data
    
    brain_seg_nii = nib.as_closest_canonical(nib.load("brainweb/subject54_crisp_v.nii.gz"))
    brain_seg = np.asarray(brain_seg_nii.get_fdata())
    # exclude structures outside the head
    brain_seg[:, 410:, 250:] = 0
    
    # %%
    # load the freesurfer brain parcellation
    aparc_nii = nib.as_closest_canonical(nib.load("brainweb/aparc+aseg_native.nii.gz"))
    aparc = np.asarray(aparc_nii.get_fdata())
    
    head_mask = np.astype(brain_seg > 0, np.int8)
    vessel_mask = np.astype(brain_seg == 8, np.int8)
    brain_mask = np.astype(aparc > 0, np.int8)
    
    parietal_cortex_numbers = np.asarray(
        [2001, 2008, 2022, 2025, 2029, 2031, 2122, 2123, 2124, 2159, 2177]
    )
    frontal_cortex_numbers = np.asarray(
        [
            2003,
            2012,
            2014,
            2017,
            2018,
            2019,
            2020,
            2024,
            2027,
            2028,
            2032,
            2106,
            2107,
            2108,
            2109,
            2110,
            2154,
            2155,
            2156,
        ]
    )
    temporal_cortex_numbers = np.asarray(
        [
            2006,
            2007,
            2009,
            2015,
            2016,
            2030,
            2033,
            2034,
            2035,
            2119,
            2131,
            2132,
            2143,
            2163,
            2164,
            2179,
            2180,
            2181,
        ]
    )
    
    occipital_cortex_numbers = np.asarray(
        [
            2005,
            2011,
            2013,
            2021,
            2114,
            2115,
            2116,
            2142,
            2160,
            2161,
            2162,
            2169,
        ]
    )
    
    frontal_cortex_mask = np.astype(np.isin(aparc, frontal_cortex_numbers), np.int8)
    frontal_cortex_mask += np.astype(np.isin(aparc, frontal_cortex_numbers - 1000), np.int8)
    
    temporal_cortex_mask = np.astype(np.isin(aparc, temporal_cortex_numbers), np.int8)
    temporal_cortex_mask += np.astype(
        np.isin(aparc, temporal_cortex_numbers - 1000), np.int8
    )
    
    parietal_cortex_mask = np.astype(np.isin(aparc, parietal_cortex_numbers), np.int8)
    parietal_cortex_mask += np.astype(
        np.isin(aparc, parietal_cortex_numbers - 1000), np.int8
    )
    
    occipital_cortex_mask = np.astype(np.isin(aparc, occipital_cortex_numbers), np.int8)
    occipital_cortex_mask += np.astype(
        np.isin(aparc, occipital_cortex_numbers - 1000), np.int8
    )
    
    basal_ganglia_mask = np.astype(np.isin(aparc, [10, 11, 12, 49, 50, 51]), np.int8)
    wm_mask = np.astype(
        np.isin(
            aparc,
            [2, 41, 7, 46, 13, 52, 1002, 2002, 1010, 2010, 1023, 2023, 1026, 2026],
        ),
        np.int8,
    )
    cerebellum_mask = np.astype(np.isin(aparc, [8, 47]), np.int8)
    ventricle_mask = np.astype(np.isin(aparc, [4, 43]), np.int8)
    
    
    bg_mask = 1.0 * head_mask
    bg_mask *= (
        (1 - frontal_cortex_mask)
        * (1 - temporal_cortex_mask)
        * (1 - parietal_cortex_mask)
        * (1 - occipital_cortex_mask)
        * (1 - basal_ganglia_mask)
        * (1 - cerebellum_mask)
        * (1 - wm_mask)
        * (1 - vessel_mask)
        * (1 - ventricle_mask)
    )
    
    frontal_cortex_mask *= 1 - vessel_mask
    temporal_cortex_mask *= 1 - vessel_mask
    parietal_cortex_mask *= 1 - vessel_mask
    occipital_cortex_mask *= 1 - vessel_mask
    basal_ganglia_mask *= 1 - vessel_mask
    cerebellum_mask *= 1 - vessel_mask
    wm_mask *= 1 - vessel_mask
    ventricle_mask *= 1 - vessel_mask
    
    
    # %%
    # read the simualted TACs
    
    with h5py.File(tac_file, "r") as f:
        t_min = f["t"][:].ravel() / 60.0  # time (frame end time) in min
        frm_dur_s = f["frm_dur_s"][:].ravel()  # frame duration in s
        C_T = f["C_T"][:].T  # tissue activity concentration - decay corrected
        C_P = f["C_P"][:].ravel()  # plasma activity concentration - decay corrected
        C_TOT = f["C_TOT"][:].ravel()  # total art. concentration - decay corrected
        K = f["K"][:].T
        m_Biexp = f["m_Biexp"][:].ravel()
        Lam = f["Lam"][:].ravel()
        Mu = f["Mu"][:].ravel()
    
    # %%
    
    act_img = xp.zeros((C_TOT.shape[0],) + bg_mask.shape, dtype=xp.float32, device=dev)
    
    for frm in range(C_TOT.shape[0]):
        print(frm)
        act_img[frm, ...] = vessel_mask * C_TOT[frm]
        act_img[frm, ...] += frontal_cortex_mask * (
            (1 - fbv) * C_T[0, frm] + fbv * C_TOT[frm]
        )
        act_img[frm, ...] += temporal_cortex_mask * (
            (1 - fbv) * C_T[1, frm] + fbv * C_TOT[frm]
        )
        act_img[frm, ...] += occipital_cortex_mask * (
            (1 - fbv) * C_T[2, frm] + fbv * C_TOT[frm]
        )
        act_img[frm, ...] += parietal_cortex_mask * (
            (1 - fbv) * C_T[0, frm] + fbv * C_TOT[frm]
        )
        act_img[frm, ...] += wm_mask * ((1 - fbv) * C_T[3, frm] + fbv * C_TOT[frm])
        act_img[frm, ...] += cerebellum_mask * ((1 - fbv) * C_T[0, frm] + fbv * C_TOT[frm])
        act_img[frm, ...] += basal_ganglia_mask * (
            (1 - fbv) * C_T[0, frm] + fbv * C_TOT[frm]
        )
        act_img[frm, ...] += bg_mask * ((1 - fbv) * C_T[0, frm] / 16.0 + fbv * C_TOT[frm])
        act_img[frm, ...] *= 1 - ventricle_mask
    
        assert xp.isclose(
            act_img[frm, ...][frontal_cortex_mask == 1].mean(),
            ((1 - fbv) * C_T[0, frm] + fbv * C_TOT[frm]),
        )
    
        assert xp.isclose(
            act_img[frm, ...][temporal_cortex_mask == 1].mean(),
            ((1 - fbv) * C_T[1, frm] + fbv * C_TOT[frm]),
        )
    
        assert xp.isclose(
            act_img[frm, ...][occipital_cortex_mask == 1].mean(),
            ((1 - fbv) * C_T[2, frm] + fbv * C_TOT[frm]),
        )
    
        assert xp.isclose(
            act_img[frm, ...][wm_mask == 1].mean(),
            ((1 - fbv) * C_T[3, frm] + fbv * C_TOT[frm]),
        )
    
        assert xp.isclose(act_img[frm, ...][vessel_mask == 1].mean(), C_TOT[frm])
    
    np.savez_compressed(
        Path(tac_file).with_suffix(".npz"),
        head_mask=head_mask,
        bg_mask=bg_mask,
        frontal_cortex_mask=frontal_cortex_mask,
        temporal_cortex_mask=temporal_cortex_mask,
        parietal_cortex_mask=parietal_cortex_mask,
        occipital_cortex_mask=occipital_cortex_mask,
        basal_ganglia_mask=basal_ganglia_mask,
        cerebellum_mask=cerebellum_mask,
        wm_mask=wm_mask,
        vessel_mask=vessel_mask,
        ventricle_mask=ventricle_mask,
        act_img=act_img,
        C_P=C_P,
        C_TOT=C_TOT,
        C_T=C_T,
        t_min=t_min,
        fbv=fbv,
        frm_dur_s=frm_dur_s,
        phantom_voxel_size=brain_seg_nii.header["pixdim"][1:4],
        K=K,
        m_Biexp=m_Biexp,
        Lam=Lam,
        Mu=Mu,
    )
    # %%
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    ax[0].plot(t_min, C_P, ".-", label="C_P")
    ax[0].plot(t_min, C_TOT, ".-", label="C_TOT")
    ax[1].plot(t_min, C_P / C_TOT, ".-", label="C_P / C_TOT")
    for i in range(C_T.shape[0]):
        ax[2].plot(t_min, C_T[i], ".-", label=f"C_T_{i}", color=plt.cm.tab10(i))
        ax[2].plot(
            t_min,
            C_T[i] * (1 - fbv) + C_TOT * fbv,
            ".--",
            label=f"C_T_{i}* (1-fbv) + C_TOT*fbv",
            color=plt.cm.tab10(i),
        )
    for axx in ax.ravel():
        axx.grid(ls=":")
        axx.legend()
        axx.set_xlabel("time (min)")
    fig.show()
    
    # %%
    # display the brain segmentation
    
    vi = pv.ThreeAxisViewer([act_img])

if __name__ == '__main__':
    main()

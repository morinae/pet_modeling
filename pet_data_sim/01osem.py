from array_api_strict._array_object import Array
import array_api_compat.cupy as cp
import parallelproj
import array_api_compat.numpy as np
from copy import copy
from pathlib import Path
import h5py

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--frm", type=int, default=0)
parser.add_argument("--num_seeds", type=int, default=2)
parser.add_argument("--phantom_file", default="tacs_C.npz")
parser.add_argument("--scanner_sens", default=1.0, type = float)
args = parser.parse_args()

# %%
dev = cp.cuda.Device(0)
scanner_res_fwhm_mm = 2.0

# number of OSEM iterations
num_iter = 6
num_subsets = 28

# frame number
frm = args.frm

# seed number
seeds = [x for x in range(args.num_seeds)]

# phantom file 
phantom_file = args.phantom_file

# scanner sensitivity
scanner_sens = args.scanner_sens

# %%
# Setup of the forward model :math:`\bar{y}(x) = A x + s`
# --------------------------------------------------------
#
# We setup a linear forward operator :math:`A` consisting of an
# image-based resolution model, a non-TOF PET projector and an attenuation model
#
# .. note::
#     The OSEM implementation below works with all linear operators that
#     subclass :class:`.LinearOperator` (e.g. the high-level projectors).

num_rings = 45
scanner = parallelproj.RegularPolygonPETScannerGeometry(
    cp,
    dev,
    radius=350.0,
    num_sides=28,
    num_lor_endpoints_per_side=16,
    lor_spacing=4.,
    ring_positions=cp.linspace(-125., 125., num_rings),
    symmetry_axis=2,
)

# %%
# setup the LOR descriptor that defines the sinogram

img_shape = (181,217,181)
voxel_size = (1.0, 1.0, 1.0)

lor_desc = parallelproj.RegularPolygonPETLORDescriptor(
    scanner,
    radial_trim=151,
    max_ring_difference=5,
    sinogram_order=parallelproj.SinogramSpatialAxisOrder.RVP,
)

proj = parallelproj.RegularPolygonPETProjector(
    lor_desc, img_shape=img_shape, voxel_size=voxel_size
)

# %% 
# load the images

phantom_data = np.load(phantom_file)
dyn_act_img = phantom_data["act_img"]

x_true = cp.asarray(dyn_act_img[frm,...], dtype=cp.float32, device=dev)

C_true_frontal = float(x_true[phantom_data["frontal_cortex_mask"] == 1].mean())
C_true_temporal = float(x_true[phantom_data["temporal_cortex_mask"] == 1].mean())
C_true_occipital = float(x_true[phantom_data["occipital_cortex_mask"] == 1].mean())
C_true_wm = float(x_true[phantom_data["wm_mask"] == 1].mean())

# true whole blood concentration
C_TOT_true = float(x_true[phantom_data["vessel_mask"] == 1].mean())

# downsample the true image by a factor of 2 in every dimension

x_true = 0.5*(x_true[::2,:,:] + x_true[1::2,:,:])
x_true = 0.5*(x_true[:,::2,:] + x_true[:,1::2,:])
x_true = 0.5*(x_true[:,:,::2] + x_true[:,:,1::2])

x_att = 0.01*cp.asarray(phantom_data["head_mask"], device = dev, dtype = cp.float32)

x_att = 0.5*(x_att[::2,:,:] + x_att[1::2,:,:])
x_att = 0.5*(x_att[:,::2,:] + x_att[:,1::2,:])
x_att = 0.5*(x_att[:,:,::2] + x_att[:,:,1::2])

# calculate the sensitivity factor due to frame length and tracer decay
frm_dur_s = float(phantom_data["frm_dur_s"][frm])
t_frm_min = float(phantom_data["t_min"][frm])

sens_factor = scanner_sens * (0.5**(t_frm_min/109)) * (frm_dur_s/600)

# Attenuation sinogram
# --------------------

# calculate the attenuation sinogram
att_sino = cp.exp(-proj(x_att))

# adjust sensitivty
att_sino *= sens_factor

# %%
# Complete PET forward model setup
# --------------------------------
#
# We combine an image-based resolution model,
# a non-TOF or TOF PET projector and an attenuation model
# into a single linear operator.

# enable TOF - comment if you want to run non-TOF
proj.tof_parameters = parallelproj.TOFParameters(
    num_tofbins=15, tofbin_width=24.5, sigma_tof=24.5
)

# setup the attenuation multiplication operator which is different
# for TOF and non-TOF since the attenuation sinogram is always non-TOF
if proj.tof:
    att_op = parallelproj.TOFNonTOFElementwiseMultiplicationOperator(
        proj.out_shape, att_sino
    )
else:
    att_op = parallelproj.ElementwiseMultiplicationOperator(att_sino)

res_model = parallelproj.GaussianFilterOperator(
    proj.in_shape, sigma=scanner_res_fwhm_mm / (2.35 * proj.voxel_size)
)

# compose all 3 operators into a single linear operator
pet_lin_op = parallelproj.CompositeLinearOperator((att_op, proj, res_model))

# %%
# Simulation of projection data
# -----------------------------
#
# We setup an arbitrary ground truth :math:`x_{true}` and simulate
# noise-free and noisy data :math:`y` by adding Poisson noise.

# simulated noise-free data
noise_free_data = pet_lin_op(x_true)

# generate a contant contamination sinogram
contamination = cp.full(
    noise_free_data.shape,
    0.5 * float(cp.mean(noise_free_data)),
    device=dev,
    dtype=cp.float32,
)

noise_free_data += contamination

# %%
# Splitting of the forward model into subsets :math:`A^k`
# -------------------------------------------------------
#
# Calculate the view numbers and slices for each subset.
# We will use the subset views to setup a sequence of projectors projecting only
# a subset of views. The slices can be used to extract the corresponding subsets
# from full data or corrections sinograms.

subset_views, subset_slices = proj.lor_descriptor.get_distributed_views_and_slices(
    num_subsets, len(proj.out_shape)
)

_, subset_slices_non_tof = proj.lor_descriptor.get_distributed_views_and_slices(
    num_subsets, 3
)

# clear the cached LOR endpoints since we will create many copies of the projector
proj.clear_cached_lor_endpoints()
pet_subset_linop_seq = []

# we setup a sequence of subset forward operators each constisting of
# (1) image-based resolution model
# (2) subset projector
# (3) multiplication with the corresponding subset of the attenuation sinogram
for i in range(num_subsets):
    #print(f"subset {i:02} containing views {subset_views[i]}")

    # make a copy of the full projector and reset the views to project
    subset_proj = copy(proj)
    subset_proj.views = subset_views[i]

    if subset_proj.tof:
        subset_att_op = parallelproj.TOFNonTOFElementwiseMultiplicationOperator(
            subset_proj.out_shape, att_sino[subset_slices_non_tof[i]]
        )
    else:
        subset_att_op = parallelproj.ElementwiseMultiplicationOperator(
            att_sino[subset_slices_non_tof[i]]
        )

    # add the resolution model and multiplication with a subset of the attenuation sinogram
    pet_subset_linop_seq.append(
        parallelproj.CompositeLinearOperator(
            [
                subset_att_op,
                subset_proj,
                res_model,
            ]
        )
    )

pet_subset_linop_seq = parallelproj.LinearOperatorSequence(pet_subset_linop_seq)

# %%
# EM update to minimize :math:`f(x)`
# ----------------------------------
#
# The EM update that can be used in MLEM or OSEM is given by cite:p:`Dempster1977` :cite:p:`Shepp1982` :cite:p:`Lange1984` :cite:p:`Hudson1994`
#
# .. math::
#     x^+ = \frac{x}{A^H 1} A^H \frac{y}{A x + s}
#
# to calculate the minimizer of :math:`f(x)` iteratively.
#
# To monitor the convergence we calculate the relative cost
#
# .. math::
#    \frac{f(x) - f(x^*)}{|f(x^*)|}
#
# and the distance to the optimal point
#
# .. math::
#    \frac{\|x - x^*\|}{\|x^*\|}.
#
#
# We setup a function that calculates a single MLEM/OSEM
# update given the current solution, a linear forward operator,
# data, contamination and the adjoint of ones.


def em_update(
    x_cur: Array,
    data: Array,
    op: parallelproj.LinearOperator,
    s: Array,
    adjoint_ones: Array,
) -> Array:
    """EM update

    Parameters
    ----------
    x_cur : Array
        current solution
    data : Array
        data
    op : parallelproj.LinearOperator
        linear forward operator
    s : Array
        contamination
    adjoint_ones : Array
        adjoint of ones

    Returns
    -------
    Array
    """
    ybar = op(x_cur) + s
    return x_cur * op.adjoint(data / ybar) / adjoint_ones


# %%
# Run the OSEM iterations
# -----------------------
#
# Note that the OSEM iterations are almost the same as the MLEM iterations.
# The only difference is that in every subset update, we pass an operator
# that projects a subset, a subset of the data and a subset of the contamination.
#
# .. math::
#     x^+ = \frac{x}{(A^k)^H 1} (A^k)^H \frac{y^k}{A^k x + s^k}
#
# The "sensitivity" images are also calculated separately for each subset.

subset_adjoint_ones = cp.zeros((num_subsets,) + pet_lin_op.in_shape, dtype=cp.float32, device=dev)
# calculate A_k^H 1 for all subsets k
for i, op in enumerate(pet_subset_linop_seq):
    print(f"calculating sensitivity image {i:02}", end = "\r")
    subset_adjoint_ones[i,...] = op.adjoint(cp.ones(op.out_shape, dtype=cp.float32, device=dev))
 
C_recon_frontal = np.zeros(len(seeds))
C_recon_temporal = np.zeros(len(seeds))
C_recon_occipital = np.zeros(len(seeds))
C_recon_wm = np.zeros(len(seeds))
C_TOT_recon = np.zeros(len(seeds))

C_recon_frontal_std = np.zeros(len(seeds))
C_recon_temporal_std = np.zeros(len(seeds))
C_recon_occipital_std = np.zeros(len(seeds))
C_recon_wm_std = np.zeros(len(seeds))
C_TOT_recon_std = np.zeros(len(seeds))

for iseed, seed in enumerate(seeds):
    # add Poisson noise
    cp.random.seed(seed)
    y = cp.random.poisson(noise_free_data)
    
    # initialize x
    x = cp.ones(pet_lin_op.in_shape, dtype=cp.float32, device=dev)
   
    # OSEM iterations
    print("")
    for i in range(num_iter):
        for k, sl in enumerate(subset_slices):
            print(f"{iseed} {seed} OSEM iteration {(k+1):03} / {(i + 1):03} / {num_iter:03}", end="\r")
            x = em_update(
                x, y[sl], pet_subset_linop_seq[k], contamination[sl], subset_adjoint_ones[k,...]
            )
    
    # %%
    # Calculation of the negative Poisson log-likelihood function of the reconstruction
    # ---------------------------------------------------------------------------------
    
    x_up = x
    x_up = cp.repeat(x_up, 2, axis=0)
    x_up = cp.repeat(x_up, 2, axis=1)
    x_up = cp.repeat(x_up, 2, axis=2)
    
    C_recon_frontal[iseed] = float(x_up[phantom_data["frontal_cortex_mask"] == 1].mean())
    C_recon_temporal[iseed] = float(x_up[phantom_data["temporal_cortex_mask"] == 1].mean())
    C_recon_occipital[iseed] = float(x_up[phantom_data["occipital_cortex_mask"] == 1].mean())
    C_recon_wm[iseed] = float(x_up[phantom_data["wm_mask"] == 1].mean())
    C_TOT_recon[iseed] = float(x_up[phantom_data["vessel_mask"] == 1].mean())

    C_recon_frontal_std[iseed] = float(x_up[phantom_data["frontal_cortex_mask"] == 1].std())
    C_recon_temporal_std[iseed] = float(x_up[phantom_data["temporal_cortex_mask"] == 1].std())
    C_recon_occipital_std[iseed] = float(x_up[phantom_data["occipital_cortex_mask"] == 1].std())
    C_recon_wm_std[iseed] = float(x_up[phantom_data["wm_mask"] == 1].std())
    C_TOT_recon_std[iseed] = float(x_up[phantom_data["vessel_mask"] == 1].std())
    
with h5py.File(f"./data/{Path(phantom_file).stem[5:]}_{scanner_sens:.1f}/{Path(phantom_file).stem}_frm_{frm:03}_s{seeds[0]}_{seeds[-1]}_sens_{scanner_sens:.1f}_recons.h5", "w") as f:
    dset = f.create_dataset("C_recon_frontal", data = C_recon_frontal)
    dset = f.create_dataset("C_recon_temporal", data = C_recon_temporal)
    dset = f.create_dataset("C_recon_occipital", data = C_recon_occipital)
    dset = f.create_dataset("C_recon_wm", data = C_recon_wm)
    dset = f.create_dataset("C_TOT_recon", data = C_TOT_recon)

    dset = f.create_dataset("C_recon_frontal_std", data = C_recon_frontal_std)
    dset = f.create_dataset("C_recon_temporal_std", data = C_recon_temporal_std)
    dset = f.create_dataset("C_recon_occipital_std", data = C_recon_occipital_std)
    dset = f.create_dataset("C_recon_wm_std", data = C_recon_wm_std)
    dset = f.create_dataset("C_TOT_recon_std", data = C_TOT_recon_std)

    dset = f.create_dataset("C_true_frontal", data = C_true_frontal)
    dset = f.create_dataset("C_true_temporal", data = C_true_temporal)
    dset = f.create_dataset("C_true_occipital", data = C_true_occipital)
    dset = f.create_dataset("C_true_wm", data = C_true_wm)
    dset = f.create_dataset("C_TOT_true", data = C_TOT_true)
    for key in ["K", "Mu", "Lam", "m_Biexp", "frm_dur_s", "t_min", "C_P", "C_TOT", "C_T", "fbv"]:
        dset = f.create_dataset(key, data = phantom_data[key])


## %%
#import pymirc.viewer as pv
#from scipy.ndimage import gaussian_filter
#
#x_np = parallelproj.to_numpy_array(x)
#vi = pv.ThreeAxisViewer([x_np, gaussian_filter(x_np,2.0), gaussian_filter(x_np,4.0)])

# gSlider-SMS offline reconstruction with B1 correction

MATLAB code for the offline super-resolution reconstruction of gSlider-SMS diffusion data. The encoding model uses slice profiles simulated at the B1 level of each voxel, so the reconstruction is corrected for B1 inhomogeneity. An example dataset is included.

## Contents

| Path | What it is |
|---|---|
| `Offline_recon_code/gSlider_load.m` | Loads the gSlider magnitude and phase NIfTI files into one complex array |
| `Offline_recon_code/B1map_load.m` | Loads the B1 map and brain mask, resamples B1 onto the gSlider grid and bins it |
| `Offline_recon_code/demo_gSlider_B1corr_v4.m` | Runs the B1-corrected gSlider reconstruction |
| `Offline_recon_code/RF_5x_90thick_1p00t_4p3mm_Verse1_CLv1_B1_Mzsim_scale1p0_MOut.mat` | Simulated slice profiles of the 5 gSlider RF encodings at 11 B1 scales (0.7 to 1.2) |
| `Offline_recon_code/offline_recon_instruction.pptx` | The original step-by-step slides |
| `Offline_recon_code/library/`, `Nifti_Analyze/`, `imagine/` | Helper functions, NIfTI input/output, image viewer |
| `providedfiles_withgSlider/` | Example data, described below |

## Requirements

- MATLAB with the Signal Processing Toolbox (for `hamming`). Tested with R2024b.
- About 32 GB of free RAM. The example data peaks at 31 GB.
- Linux or macOS. `load_nifti.m` decompresses `.nii.gz` files with `zcat` or `gunzip` into `/tmp`.

## Quick start with the example data

In MATLAB, change into `Offline_recon_code` and run:

```matlab
gSlider_load            % select ..._23001_mag.nii, then ..._24001_ph.nii
B1map_load              % select ..._tfl_b1map_..._27001.nii, then bet_..._mask.nii.gz
demo_gSlider_B1corr_v4  % select RF_5x_90thick_..._MOut.mat from this folder
```

Run the three scripts in this order in one session. `gSlider_load.m` starts with `clear`, and each later script uses variables the earlier ones leave in the workspace.

The result, `Img_SuperRes`, stays in the workspace and is not written to disk. For the example data it is 224 × 224 × 130 × 7: the 26 slabs of 5 mm become 130 slices of 1 mm, for each of the 7 diffusion volumes. The example takes about 8 minutes on an 80-core Linux server. The last line of the demo swaps the first two axes, so `Img_SuperRes` is transposed in-plane relative to the input NIfTI files.

## Preparing your own data

1. Convert the DICOMs to NIfTI (`.nii`): the gSlider magnitude, the gSlider phase and the B1 (flip-angle) map, for example with dcm2niix.
2. Make a brain mask from the gSlider magnitude with FSL BET, for example `bet mag.nii mag_brain -f 0.25 -m`, which writes `mag_brain_mask.nii.gz`.
3. Check that:
   - the fourth dimension of the magnitude and phase files has the 5 RF encodings innermost: encodings 1 to 5 of the first diffusion volume, then encodings 1 to 5 of the next one, and so on (`nrf_d = 5` in the demo);
   - the B1 map covers the whole gSlider slab volume. Its position is taken from the NIfTI headers, so its slices do not have to line up with the gSlider slabs. Outside the B1 field of view, values are extended from the nearest edge;
   - the brain mask is on the gSlider grid, which is the case when BET is run on the gSlider magnitude;
   - the phase is either raw Siemens values (0 to 4095) or dcm2niix output (−4096 to 4094, stored as raw values with `scl_slope = 2` and `scl_inter = −4096`). `gSlider_load.m` tells the two apart by whether any value is negative, and scales both to −π to π.

## How the B1 map is used

`B1map_load.m` resamples the B1 map onto the gSlider grid in world coordinates. It uses the NIfTI sform of the B1 map and of the brain mask, and linear interpolation. It then divides by 800 (the Siemens `tfl_b1map` stores the flip angle × 10, for a nominal 80°), applies the mask, and snaps each voxel to the nearest of 0.7:0.05:1.2, the B1 scales of the simulated slice profiles. Voxels outside the mask are set to 1.

This replaces the B1 handling of the original toolbox. The original reordered the B1 slices as if they were stored in interleaved acquisition order (the even-numbered slices, then the odd-numbered ones) and stretched the B1 field of view onto the gSlider one. B1 maps converted with dcm2niix are already in spatial slice order, so the reordering scrambled them. The old code still runs, but its result is overwritten.

Assumptions:

- `flag_b1mask = 0`, the default in `B1map_load.m`.
- The head does not move between the B1 and gSlider scans. The two are aligned from their headers, not registered.

## Changes from the original toolbox

- `B1map_load.m`: B1 resampling in world coordinates, described above.
- `gSlider_load.m`: phase scaling. The original formula, (p − 2048)/2048 · π, assumes raw values. With dcm2niix files it doubled the phase and shifted it by π, so the phase spanned −3π to π.
- `demo_gSlider_B1corr_v4.m`: the last line, which averaged volumes 2 to 31 into `meanDWI`, is commented out. It failed for datasets with fewer than 31 volumes, such as the example.

## Example data

All files are in `providedfiles_withgSlider/`.

| File | Contents |
|---|---|
| `24021520_ep2d_gslider_mgh_1mm_R3MB2_pf68_20240215134355_23001_mag.nii` | gSlider-SMS diffusion magnitude: 224 × 224 × 26 slabs × 35 volumes (7 diffusion volumes × 5 RF encodings), 1 mm in-plane, 5 mm slabs |
| `24021520_ep2d_gslider_mgh_1mm_R3MB2_pf68_20240215134355_24001_ph.nii` | The matching phase |
| `bet_24021520_ep2d_gslider_mgh_1mm_R3MB2_pf68_20240215134355_23001_mask.nii.gz` | FSL BET brain mask on the gSlider grid |
| `24021520_tfl_b1map_20240215134355_27001.nii` | Siemens `tfl_b1map` flip-angle map: 64 × 64 × 26, 3.5 × 3.5 × 5 mm |

## Credits and references

The original toolbox was written by Congyu Liao and Kawin Setsompop at Massachusetts General Hospital. Third-party code keeps its own license: `imagine` (Christian Wuerslin, `imagine/license.txt`), the NIfTI tools (Jimmy Shen, `Nifti_Analyze/license.txt`) and `load_nifti.m` (FreeSurfer, MGH).

- Setsompop K, et al. High-resolution in vivo diffusion imaging of the human brain with generalized slice dithered enhanced resolution simultaneous multislice (gSlider-SMS). *Magn Reson Med.* 2018;79(1):141–151. [PMC5585027](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5585027/)
- Liao C, et al. High-fidelity, high-isotropic-resolution diffusion imaging through gSlider acquisition with B1+ and T1 corrections and integrated ΔB0/Rx shim array. *Magn Reson Med.* 2020;83(1):56–67.

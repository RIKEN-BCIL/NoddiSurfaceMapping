# NoddiSurfaceMapping

NODDI (neurite orientation dispersion and density imaging) of the cerebral
cortex, fitted per vertex on the cortical surface and mapped to CIFTI
grayordinates, for data preprocessed with the HCP Pipelines.

**v2.0.0-beta** — see [RELEASE_NOTES.md](RELEASE_NOTES.md). v2 no longer uses
AMICO. The v1 code (AMICO, volume fit) is at tag `v1.0`.

![noddi](https://user-images.githubusercontent.com/16514166/40781643-a102ca94-6517-11e8-81e1-82e9556199de.png)

## What it does

1. **prep** — samples the diffusion signal once per vertex (cubic) on a
   native-mesh layer surface between white and pial (default a quarter of the
   way from white, `--equi-frac-distance=0.25`), and packs the vertices into a
   NIfTI for the fitter. With `--gradnonlin` the vertices are grouped into bins
   of the effective b-value scale, each with its own b-table.
2. **fit** — NODDI (Watson, parallel diffusivity 1.1 µm²/ms for grey matter)
   per vertex, with cuDIMOT on the GPU, or in MATLAB (the NODDI toolbox, or a Rician NODDI
   of your own via `--rnoddidir`).
3. **map** — the fitted maps are smoothed on the layer surface by its own mean
   vertex spacing and resampled to `32k_fs_LR` CIFTI (MSMSulc or MSMAll).

Outputs: `noddi_ficvf` (NDI), `noddi_odi`, `noddi_fiso`, `noddi_kappa`, and the
quality maps `noddi_snr`, `noddi_gmpv` (grey fraction of the sample) and
`noddi_valid`, under `MNINonLinear/Results/<dwi>/SurfaceFitMapping<suffix>/`.
Nothing is excluded by default; the quality maps are there for the decision to
be made downstream (`-S` flags low-SNR vertices if wanted).

## Requirements

- A subject processed with the HCP Pipelines (PreFreeSurfer, FreeSurfer,
  PostFreeSurfer, DiffusionPreprocessing); `HCPPIPEDIR` and `EnvironmentScript`
  exported
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl) 6.0
- [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)
  1.5 or later (`CARET7DIR`); **2.1.0 or later for `--mapping=hippunfold`**
- python3 with numpy and nibabel (and scipy for `--mapping=hippunfold`)
- For the default fitters: [cuDIMOT](https://users.fmrib.ox.ac.uk/~moisesf/cudimot)
  with the `CorticalNODDI_Watson` model (below) and an NVIDIA GPU
- For `--fitter=matlab`: MATLAB, the
  [NODDI toolbox](http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.NODDImatlab) and
  niftimatlib (`--noddidir=`, `NIFTIMATLIB=`; default `./NODDI/`)
- For `--fitter=matlab-rician`: MATLAB and a Rician NODDI implementation, given
  with `--rnoddidir=<dir>` (or `RNODDIDIR`). The folder must provide
  `noddi_rician_bin(datafile, maskfile, bvalsfile, bvecsfile, outprefix)`
  writing `<outprefix>_ficvf.nii`, `_fiso.nii`, `_kappa.nii`, `_odi.nii`,
  `_fibredirs_{x,y,z}vec.nii` and `_fmin.nii`, as the NODDI toolbox's
  `SaveParamsAsNIfTI` does; it is called once per bin
- For `--mapping=hippunfold`: HippUnfold surfaces in
  `<Study>/<Subject>/T1w/HippUnfold/T1w/<density>/`

### cuDIMOT model

`CorticalNODDI_Watson` is cuDIMOT's own `NODDI_Watson` with one change: the
parallel diffusivity is 1.1 µm²/ms (grey matter; Fukutomi et al., 2018)
instead of 1.7. Nothing else in the model is modified.

Use the maintained cuDIMOT, <https://github.com/SPMIC-UoN/cudimot> (the one
NoddiSurfaceMapping is developed with). Its `NODDI_Watson` includes a 2022 fix
of the kappa boundary between the exact and the approximate Watson spherical
harmonics (0.1 → 0.4 in `WatsonFunctions.h`); older copies of cuDIMOT, such as
the one at git.fmrib.ox.ac.uk, still have 0.1, which can give large errors in
the model prediction at high dispersion.

1. Clone it and set up the build environment as its README describes:
   `FSLDEVDIR` (where `make install` puts the programs), FSL's development
   settings (`source $FSLDIR/etc/fslconf/fsl-devel.sh`) and `CUDA` (the CUDA
   toolkit).
2. Make the model, in the cuDIMOT source tree:

   ```
   cp -r mymodels/NODDI_Watson mymodels/CorticalNODDI_Watson
   cd mymodels/CorticalNODDI_Watson
   mv Pipeline_NODDI_Watson.sh Pipeline_CorticalNODDI_Watson.sh
   mv NODDI_Watson_finish.sh CorticalNODDI_Watson_finish.sh
   sed -i 's/NODDI_Watson/CorticalNODDI_Watson/g' Pipeline_CorticalNODDI_Watson.sh CorticalNODDI_Watson_finish.sh
   sed -i 's/#define Dparallel 0.0017/#define Dparallel 0.0011/' diffusivities.h
   cd ../..
   ```

   (The two scripts are renamed because the Makefile looks for
   `Pipeline_<model>.sh` and `<model>_finish.sh`.)
3. Build and install, then point NoddiSurfaceMapping at the install:

   ```
   modelname=CorticalNODDI_Watson make
   modelname=CorticalNODDI_Watson make install      # into $FSLDEVDIR/bin
   export CUDIMOTDIR=$FSLDEVDIR                      # default /usr/local/CUDIMOT
   ```

NoddiSurfaceMapping needs these in `$CUDIMOTDIR/bin`: `CorticalNODDI_Watson`,
`split_parts_CorticalNODDI_Watson`, `merge_parts_CorticalNODDI_Watson` and
`cart2spherical`. It runs cuDIMOT's fitting itself — the DTI that initialises
the orientation (with `--gradnonlin` when given), the grid search, the two
fractions, then every parameter, with the options and initialisation of
cuDIMOT's `Pipeline_<model>.sh` — one step after another in the same job, so
cuDIMOT's `Pipeline_<model>.sh`, `jobs_wrapper.sh` and `Run_dtifit.sh` are not
used.

For `--mapping=grayordinate+whiteordinate`, build the unmodified `NODDI_Watson`
(1.7 µm²/ms) too, into a second `FSLDEVDIR`, and set `CUDIMOTDIR_WM` to it.

Notes (built and installed this way in 2026-10 from SPMIC-UoN/cudimot master
5f9e4ff, with FSL 6.0.7 and CUDA 12.6):

- CUDA 12.6 accepts host compilers up to gcc 13, while FSL's conda environment
  puts its own gcc 14 first on `PATH`. Put a supported `g++` first before
  `make` (e.g. `export PATH=/usr/bin:$PATH` where the system gcc is 11).
- The Makefile builds for GPU architectures sm_50 to sm_86 (`GPU_CARDs`); add
  newer ones there (e.g. sm_89, sm_90) if your GPUs need them.
- With several fits on one machine, the GPUs in persistence mode
  (`nvidia-smi -pm 1`) and exclusive-process compute mode (`nvidia-smi -c 3`)
  keep runs from sharing a device.

## Installation

```
git clone https://github.com/RIKEN-BCIL/NoddiSurfaceMapping.git
```

No build step; run `NoddiSurfaceMapping.sh` from the clone (it finds
`scripts/` next to itself).

## Usage

```
NoddiSurfaceMapping.sh [options] <StudyFolder> <SubjectID>

# HCP-YA, MSMAll, with the gradient nonlinearity tensor
NoddiSurfaceMapping.sh -M --gradnonlin=<Study>/<Subject>/T1w/Diffusion/grad_dev.nii.gz <Study> <Subject>

# cortex only, prep on a CPU node, fit on a GPU node, map anywhere
NoddiSurfaceMapping.sh --mapping=cortex --runmode=prep <Study> <Subject>
NoddiSurfaceMapping.sh --mapping=cortex --runmode=fit  <Study> <Subject>
NoddiSurfaceMapping.sh --mapping=cortex --runmode=map  <Study> <Subject>

# hippocampus (wb_command >= 2.1.0)
WORKBENCHDIR_HIPP=<dir of wb_command 2.1> NoddiSurfaceMapping.sh --mapping=hippunfold <Study> <Subject>

# v1-style volume fit and ribbon mapping (Fukutomi et al., 2018)
NoddiSurfaceMapping.sh --mapping=volume <Study> <Subject>
```

Run it with no arguments for the full option list (`--mapping`, the sampling
depth, `--fitter`, `--gradnonlin`, `--runmode`, `-a` species, `-S`, ...).

### Environment variables

| Variable | Meaning |
|---|---|
| `HCPPIPEDIR`, `EnvironmentScript` | HCP Pipelines (required) |
| `CARET7DIR` | Workbench, normally set by the HCP environment script |
| `WORKBENCHDIR_HIPP` | a wb_command >= 2.1.0, for `--mapping=hippunfold` |
| `CUDIMOTDIR`, `CUDIMOTDIR_WM` | cuDIMOT installs (grey matter; white matter) |
| `CUDIMOT_FSLDIR` | the FSL cuDIMOT was built against, if not `FSLDIR` |
| `CUDA_VISIBLE_DEVICES` | the GPU(s) to use, as the scheduler sets it (a site-specific SGE example is in the script, commented out) |
| `NSM_MATLAB_QUEUE` | `fsl_sub` queue for the MATLAB bins |
| `NODDIDIR`, `NIFTIMATLIB` | NODDI toolbox and niftimatlib, for `--fitter=matlab` |
| `RNODDIDIR` | the Rician NODDI, for `--fitter=matlab-rician` |

## License

MIT for the files of this repository; see [LICENSE](LICENSE) for the
third-party software it calls.

## References

Fukutomi H, Glasser MF, Zhang H, Autio JA, Coalson TS, Okada T, Togashi K,
Van Essen DC, Hayashi T (2018) Neurite imaging reveals microstructural
variations in human cerebral cortical gray matter. *NeuroImage* 182:488–499.
[doi:10.1016/j.neuroimage.2018.02.017](https://doi.org/10.1016/j.neuroimage.2018.02.017)
([BALSA](https://balsa.wustl.edu/study/show/k77v))

Surface fit and CIFTI mapping (v2): Hayashi et al., in prep.; Yoshida et al., in prep.

Hernandez-Fernandez M, et al. (2019) Using GPUs to accelerate computational
diffusion MRI: from microstructure estimation to tractography and connectomes.
*NeuroImage* 188:598–615. (cuDIMOT)

Zhang H, Schneider T, Wheeler-Kingshott CA, Alexander DC (2012) NODDI:
practical in vivo neurite orientation dispersion and density imaging of the
human brain. *NeuroImage* 61:1000–1016.

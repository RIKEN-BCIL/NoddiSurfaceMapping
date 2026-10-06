# Release notes

## v2.0.0-beta (2026-10)

A rewrite. v1 fitted NODDI in the volume with AMICO and mapped the result to the
surface; v2 maps the diffusion signal to the surface first and fits NODDI per
vertex, so that no voxel straddling white matter or CSF enters the fit. This is
a beta: the method is the one described in Hayashi et al. (in prep.) and Yoshida
et al. (in prep.), but option names and output names may still change before
v2.0.0.

### New
- **Surface fit** (`--mapping=cortex|grayordinate|grayordinate+whiteordinate`,
  default `grayordinate`). One cubic sample of the diffusion signal per vertex
  on a native-mesh layer surface, NODDI fitted per vertex, the fitted maps
  smoothed on that surface by its own mean vertex spacing (smoothing after the
  fit, where it costs no bias), then resampled to 32k_fs_LR CIFTI.
  `grayordinate` adds the subcortical grey (wmparc) fitted in the volume;
  `grayordinate+whiteordinate` adds the white matter, fitted with a parallel
  diffusivity of 1.7 µm²/ms instead of the cortical 1.1.
- **Sampling depth**: `--equi-frac-distance` (default 0.25 of white to pial),
  `--equi-vol-distance`, `--equi-abs-distance` with `--abs-distance-cap`. The
  depth is part of every output name (`_fracdist0.25`, ...).
- **Hippocampus** (`--mapping=hippunfold`): the same per-vertex fit on the
  HippUnfold hippocampal and dentate surfaces, fitted at the finest density and
  downsampled on the surface to the others. Needs Connectome Workbench
  **wb_command 2.1.0 or later** (native HIPPOCAMPUS / HIPPOCAMPUS_DENTATE CIFTI
  structures); set `WORKBENCHDIR_HIPP` if that is not the one in `CARET7DIR`.
- **Fitters** (`--fitter=`): `cudimot-mcmc` (default; cuDIMOT on the GPU,
  MCMC with Rician noise), `cudimot-mle` (Levenberg-Marquardt), `matlab` (the
  NODDI toolbox, CPU) and `matlab-rician` (a Rician NODDI in MATLAB, supplied
  by the user with `--rnoddidir`; the interface is in the README). Each writes
  its own folder and CIFTI suffix.
- **Gradient nonlinearity** (`--gradnonlin=<grad_dev>`, `--gnlbin=`): DTIFit
  gets the tensor, and the NODDI fits are done in bins of the effective b-value
  scale, each bin with its own b-table. Works with every fitter.
- **Quality maps**, always written and carried to CIFTI: `noddi_snr` (b0 SNR
  per vertex), `noddi_gmpv` (grey-matter fraction of the sample) and
  `noddi_valid`. Exclusion is left to the user: `-S <frac>` flags low-SNR
  vertices (default 0, nothing flagged), `-N` keeps them NaN.
- **Stages** (`--runmode=prep,fit,map`), a CPU thread cap (`--nthr`), and a
  scratch directory for the one-off decompression of the data (`--tmpdir`).
- **Species**: `-a` 0 human, 1 macaque (Mac30BS), 1.1 cynomolgus, 1.2 rhesus,
  1.3 Japanese (snow) macaque, 2 marmoset, 3 night monkey.

### Changed
- **AMICO removed.** v1's README asked for AMICO with its parallel diffusivity
  edited to 1.1e-3 mm²/s; v2 does not use AMICO at all.
- v1's route (volume fit, then myelin-style ribbon mapping; Fukutomi et al.,
  2018) is kept as `--mapping=volume`, to reproduce earlier results.
- cuDIMOT's fitting steps (initialising DTI, grid search, the two fractions,
  every parameter) are run by NoddiSurfaceMapping itself, one after another in
  the same job; only cuDIMOT's programs are needed (see the README for building
  the `CorticalNODDI_Watson` model from https://github.com/SPMIC-UoN/cudimot).
- One subject per call: `NoddiSurfaceMapping.sh [options] <StudyFolder> <SubjectID>`
  (v1 took several subject IDs).
- `-s` (mapping only) is replaced by `--runmode=map`.
- DTI b-value thresholds (`-t`) default to `1050,100,50`.
- Outputs of the surface route are under
  `MNINonLinear/Results/<dwi>/SurfaceFitMapping<suffix>/`; the volume route
  writes `RibbonVolumeToSurfaceMapping<suffix>/`.

### Known limitations of the beta
- Developed and tested mainly on human HCP-style 3 T data (HCP-YA, 1.25 mm;
  an HCP-style protocol at 1.7 mm). The non-human species, the white-matter
  ordinates and the hippocampal route have had less testing.
- `--fitter=matlab` runs on the CPU and is much slower than cuDIMOT.

## v1.0 (2018-05)

NODDI with AMICO in the volume, cortical ribbon mapping to the surface
(Fukutomi et al., NeuroImage 2018). Last change on master: 2023-01 (DTI tensor
saved, DTI volumes sorted by label).

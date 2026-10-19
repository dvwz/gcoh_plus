# gcoh+: Global Coherence Toolbox Plus

MATLAB toolbox for global coherence analysis of multichannel EEG. It covers:
- multitaper cross-spectral matrices;
- frequency-domain bootstrap replicates;
- eigendecomposition into coherent modes;
- rank-reduced reconstruction (coherent power);
- angular distances between coherent subspaces.

Author: David W. Zhou

## Requirements

- MATLAB (Statistics and Machine Learning Toolbox not required by the toolbox itself)
- [Chronux](http://chronux.org) on the MATLAB path (`getparams`, `getfgrid`, `dpsschk`, `mtfftc`, `change_row_to_column`)

## Installation

```matlab
addpath('/path/to/gcoh_plus')
addpath(genpath('/path/to/chronux'))
```

## Usage

```matlab
% data: nsamples x nchannels
mtparams   = struct('Fs', 250, 'tapers', [3 5], 'fpass', [0 50], ...
                    'pad', -1, 'movingwin', [3 3]);
gcohparams = struct('version', 'epoch', 'inference', 'fdb', 'Niter', 200);

gcohout = gcoh_plus(data, [], mtparams, gcohparams);

% global coherence: fraction of power in the principal mode
gcoh = squeeze(gcohout.eigenvalues(1,:,:) ./ sum(gcohout.eigenvalues,1));   % Nfreqs x Niter

% rank-1 reconstruction (coherent cross-spectrum) from the first mode
Pxy1 = recon_mtx(gcohout.eigenvectors_l, gcohout.eigenvalues, gcohout.eigenvectors_r, 1);
```

`gcohparams` options:

| Field | Values | Meaning |
|---|---|---|
| `version` | `'window'` (default), `'epoch'` | Per 3-s window, or one estimate per recording (median across windows) |
| `inference` | `'none'` (default), `'fdb'` | No resampling, or frequency-domain bootstrap across tapers and windows |
| `Niter` | default 1 (`none`) or 200 (`fdb`) | Number of bootstrap replicates |
| `step` | `'all'` (default), `'mtx'` | Full analysis, or return the cross-spectra only (`gcohout.cross_spectra`) |

## Contents

| File | Description |
|---|---|
| `gcoh_plus.m` | Top-level function: windows the data, computes cross-spectra, optional bootstrap, and eigendecomposition per frequency |
| `pkmtxsptrm.m` | Multitaper cross-spectrum for all tapers (Nfreqs × Ntapers × Nchans × Nchans) |
| `fdb_mtx.m` | Frequency-domain bootstrap: resamples tapers (and windows, for `'epoch'`), averages across tapers, and takes the median across windows of the real and imaginary parts separately |
| `recon_mtx.m` | Reconstructs the cross-spectral matrix from selected modes (e.g. coherent power = diagonal of the rank-1 reconstruction) |
| `angdist.m` | Pairwise angular distance (principal angle) between eigenvectors |
| `modalproj.m` | Modal projection and network-weighted spectrum of multitaper estimates onto coherent modes |
| `tutorial/` | Live-script tutorial (`tutorial.mlx`) with an example EEG test dataset (`tutorial/data/`) |

## Citation

If you use gcoh+, please cite:

- Zhou, D. W., Conte, M. M., Curley, W. H., Spencer-Salmon, C. A., Chatelle, C., Rosenthal, E. S., Bodien, Y. G., Victor, J. D., Schiff, N. D., Brown, E. N., & Edlow, B. L. (in press). Alpha coherence is a network signature of cognitive recovery from disorders of consciousness. *Cortex*.
- Weiner, V. S., Zhou, D. W., Kahali, P., et al. (2023). Propofol disrupts alpha dynamics in functionally distinct thalamocortical networks during loss of consciousness. *PNAS*, 120(11), e2207831120.

Background on global coherence:
- Mitra, P. (2007). *Observed Brain Dynamics*. Oxford University Press.
- Cimenser, A., Purdon, P. L., Pierce, E. T., et al. (2011). Tracking brain states under general anesthesia by using global coherence analysis. *PNAS*, 108(21), 8832–8837.

The version used in Zhou et al. (*Cortex*) is tagged `cortex-2026` (commit `556860c`).

## License

MIT. See [LICENSE](LICENSE).

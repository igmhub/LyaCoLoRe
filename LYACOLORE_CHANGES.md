# LyaCoLoRe: faster `make_transmission`, and fixes for numpy 2 and scipy 1.14 (branch `gpu-port`)

Three commits and this file on top of `master` of `igmhub/LyaCoLoRe` (`a600bb8`, pull request #103 "lyacolore-dr2-fs"
merged on 2026-10-09). No GPU is involved; nothing changes in the command lines or in the output format. The
transmission files are byte-identical to the ones made without the speed-ups (see below).

## Result

A 1024^3 CoLoRe box (7,446,917 objects in `master.fits`, 3072 pixels), `make_transmission.py` with 128 processes on
one CPU node, the same input and seed, the original code (with the scipy fix below, which the new scipy needs) and
this branch:

| | Original | This branch |
|---|---|---|
| `make_master.py` | 236.8 s | 243.0 s |
| `make_transmission.py` (3072 pixels) | 1607 s | 467 s (3.4 times faster) |
| Summed peak memory of the processes | 199 GB | 198 GB |
| Output | | all 3072 `transmission-*` files identical, compared HDU by HDU and column by column with zero tolerance |

For three pixels in one process the time went from 176 s to 36 s. One run each; the times are single measurements.

## What changed

| Where | Change |
|---|---|
| `py/utils.py` | `make_MOCKID` was a Python loop over every row of the CoLoRe file, repeated for every pixel and every input file. It is now one array expression with the same result. |
| `py/RSD.py` | The RSD weights without thermal effects (`thermal = False`) found the overlapped cells with an `np.where` over all cells for every cell. They now use `searchsorted`. The weight matrices were compared directly (`indptr`, `indices` and `data` identical). The thermal path is not changed. |
| `py/independent.py`, `py/simulation_data.py`, `py/tuning.py`, `py/utils.py` | numpy 2: `np.in1d` became `np.isin` and `np.trapz` became `np.trapezoid`. |
| `py/DLA.py` | `scipy.interpolate.interp2d` was removed in scipy 1.14, so the DLA tables could not be read with a current environment. It is replaced by a small function that does the same linear interpolation on scattered points with `bisplrep` (`kx = ky = 1`, `s = 0`) and `bisplev`, with the sorting and transposition that `interp2d` did. |

## What was checked

- Three pixels (`transmission-16-*.fits.gz`) compared with zero tolerance between the original and the changed
  code, in two Python environments; the RSD weight matrices compared directly; then the whole 1024^3 mock above
  (3072 of 3072 files identical).
- The speed-up commit was checked again on top of the current `master` (24 pixels of the same box, DLAs switched
  off because of the table size above): 25 of 25 files identical with and without it, 24.6 s against 67.7 s.
- The `interp2d` replacement has **not** been compared with the output of a scipy older than 1.14: both sides of
  my comparisons used the replacement, because the original code does not run on the current scipy. A comparison
  needs an environment with an older scipy. On the current `master` the DLA part itself was not run (see above).

## Things noticed

- **Row order of `master.fits`.** `scripts/make_master.py` lists the input files with an unsorted `glob.glob`, so the
  order of the rows of `master.fits` follows the directory listing. Two runs gave the same rows in a different order
  (the first two files swapped). The content is the same and the skewer files are not affected, but the order feeds
  anything that is drawn row by row downstream (for example the quasar catalogue seed), so a `sorted(...)` around
  the `glob` would make it reproducible. The Lya group knows about it and judged it not worth changing so far; I did
  not change it here.
- **`F_LYA` depends on the environment.** The same pixels give different `F_LYA` in DESI 26.9 and in the updated
  `vega_master` (more than 1e-6 in 16% of the pixels, largest difference 0.2), with the original and the changed
  code agreeing exactly inside each environment. I have not found which package causes it.
- **The 2LPT DLA code of the merged pull request.** (1) `scripts/make_transmission.py` on `master` has the guard
  "Adding DLAs from density input skewers is not possible yet!" switched on again (`if (args.skewer_type=='density')
  & args.add_DLAs:`), while the commented-out version is in the fork this work started from, so the 2LPT
  configuration with DLAs stops there. (2) The new DLA code reads `input_files/dla_files/expected_2lpt_flagged_cells_vs_z_cell.txt`,
  a table with one entry per redshift cell of the production grid (4096^3); on a 1024^3 test box the sizes do not
  match (`boolean index did not match indexed array`). So the DLA part cannot be tested on a small box.

## Commits

`2e4e479` numpy 2 (`np.in1d`, `np.trapz`); `241da75` vectorised `make_MOCKID` and RSD weights; `c8b751a`
`interp2d` replacement in `DLA.py`. (The same three changes are `31525fe`, `6f345d9` and `6c79f36` in the fork
`mruizherrerab/LyaCoLoRe`, where they sit on top of Manuel's earlier version of the 2LPT DLA code.)

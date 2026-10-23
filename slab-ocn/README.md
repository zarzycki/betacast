# Slab ocean w/ empirical mixing

Citation: **C. M. Zarzycki**. Tropical cyclone intensity errors associated with lack of two-way ocean coupling in high-resolution global simulations. Journal of Climate, 290 (23):0 8589--8610, 2016. [10.1175/JCLI-D-16-0273.1](http://dx.doi.org/10.1175/JCLI-D-16-0273.1).

## General workflow

1. Generate a DOCN file that contains three variables: SST_cpl, hblt, and qdp.
2. Either copy `docn_comp_mod.F90` to SourceMods or patch in the modifications and build the model.
3. `./xmlchange DOCN_MODE="som"`
4. Copy `user_docn.streams.txt.som` to your case directory.
5. (optional) create/update "qdp" bias correction (see below).

NOTE: You may want to ensure the ice datastream is consistent -- this slab ocean will not try to prognose sea ice (and in fact, can have SSTs below freezing, albeit marginally!)

### Modify slab settings

To modify the slab ocean configuration (see Zarzycki 2016), you need to edit the source code in `docn_comp_mod.F90`. Here, you can either change the turbulent cooling function (X_cool) or any of the parameters (e.g., relaxation time, reference mixed layer depth, etc.).

### Generate SST datafile

`h-qdp.py` will add the hblt (GODAS or Z16) and qdp (Z16 or zero) fields to an SST-only file specified by the user. If the SST is transient (i.e., AMIP) it will just add the vars as a repeating monthly climatology. This is called automatically by Betacast when `do_slab = true`.

```
python h-qdp.py --sst_file sst.nc --output_file sst_slab.nc [--which_h GODAS|Z16] [--which_q Z16|NONE]
```

If starting from an external monthly SST/ice dataset (variables `sst` and `icec`), `external-to-docn.py` will fill missing values and add `date`/`datesec` at the midpoint of each month (see `batch-alcf-docn.sh` for an example workflow through `bcgen`).

```
python external-to-docn.py --infile sst_cice.mon.mean.nc [--output_file out.nc]
```

### Creating bias correction fields

On the slab data stream, one could add a monthly bias correction field under the variable "qdp." This was used in Zarzycki (2016) to improve the correlation of the slab ocean climatology with the fixedSST climatology.

To generate:

1. Run an F simulation (>= 1 year, used 10 years in Z16) with the prescribed SST file (i.e., FIXEDSST), outputting monthly averages.
2. Run an F simulation (>= 1 year, used 10 years in Z16) with the slab ocean turned on (i.e., SLAB) (with qdp set to 0.0).
3. (optional) if being lazy and using linear interp'ed files, clean up the SST field for reasonableness with `ncap2 -O -s where(SST < 271.34 && SST > 0.0001) SST=271.34`
4. Create a monthly climatology of the difference `BIAS = SLAB - FIXEDSST`.
5. Add the monthly climatology to the SOM stream as "qdp" with `add-sst-correction-to-slab.py`, which expects files `diff01.nc` ... `diff12.nc` (variable `SST`) in `--bias_dir`, bilinearly interpolates them to the slab file grid (`xc`, `yc`), and overwrites `qdp` in place (months without a file are set to 0).

```
python add-sst-correction-to-slab.py --slab_file climo_1x1_slab.nc --bias_dir /path/to/bias_files/
```

## Troubleshooting

### CICE convergence errors

When using the slab with CESM3, originally ran into errors:

```
(picard_nonconvergence):picard convergence failed!
...
(picard_solver) picard_solver: Picard solver non-convergence
(icepack_warnings_setabort) T :file icepack_therm_mushy.F90
...
(icepack_warnings_aborted) ... (temperature_changes_salinity)
(temperature_changes_salinity): Picard solver non-convergence (no snow)
(icepack_warnings_aborted) ... (thermo_vertical)
(icepack_warnings_aborted) ... (icepack_step_therm1)
```

... while these errors generally indicate some sort of initial conditions mismatch, they suspiciously reported (very) bad temperatures, etc.

These work solved by enforcing double precision for all stream variables (`SST\_cpl`, `ice\_cov`, `qdp`, `hblt`).

Another (possible) solution, would be to try:

```
ice_ic = 'default'
ktherm = 1
nslyr = 1
nilyr = 4
```

in user\_nl\_cice to roll back some of the CICE physics (see [this thread](https://bb.cgd.ucar.edu/cesm/threads/force-ktherm-1-to-disable-cice-mushy-layers-scheme-in-cesm-2-3-a17b.11829/)).

### Weird SST behavior near coastal edges

The model will run if the `SSTICE_MESH` in `env_run.xml` and the DOCN `meshfile` are different, however, it may produce spurious or weird behavior in coastal regions. Recommend enforcing that all SST and ice data are on the same grid using the same mesh file.
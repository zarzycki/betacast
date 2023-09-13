# Some notes on physgrid topo
With the physgrid option in CAM, this is how files are read.

- bnd_topo lives on PG
- nudging files live on PG
- ncdata file lives on GLL

This can be tricky when performing the hydrostatic correction in atm_to_cam. This is because atm_to_cam assumes the PHIS field provided in `model_topo_file` is on the same grid as the target for the state variables (e.g., PS, T) which is GLL. However, there may only exist a PG bnd_topo file that for that resolution.

In this folder, you can store PHIS\_gll output from a PG run. If there is a stable compset already, the key is to just output one timestamp of PHIS\_gll in a fincl field. For example, in `user_nl_cam` add:

```
nhtfrq=0,-3
mfilt=1,1
fincl2='PHIS_gll:I'
```

Run a very short simulation such that the file is written out as `FHIST-ne30pg3-ndg-ERA5-Q24-N23-x101.cam.h1.2018-01-01-00000.nc`.

Then you need to average the time dimension away and compress, both with NCO

```
ncwa -a time FHIST-ne30pg3-ndg-ERA5-Q24-N23-x101.cam.h1.2018-01-01-00000.nc tmp.nc
ncks -O -4 -L1 tmp.nc tmp.nc
mv tmp.nc $BETACAST/atm_to_cam/topo_files/ne30np4_from_ne30pg3_gmted2010_modis_bedmachine_nc3000_Laplace0100_20230105.nc
```

NOTE: compression is useful because the PHIS field has so many zeros that NCO/zlib can compress this file down to the MB scale for even ultra-high-resolution meshes.

A smart strategy is just to name the file `neXXnpY_from_$PGTOPOFILENAME` so that we can guarentee a given topography file is matched to the GLL grid (e.g., if one creates a smoother PG topo field, it won't match as well to the above topo file on the GLL grid for ne30).

Once these files are stashed in `$BETACAST/atm_to_cam/topo_files/` they can be used as `model_topo_file` when generating initial conditions on the GLL grid as long as the topography file remains the same.

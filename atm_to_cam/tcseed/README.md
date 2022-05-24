## To add a vortex to an SE analysis state

In `seed-tc-in-ncdata.ncl`...

First set

```
invert_vortex = False
origfile="???"
seedfile="???_seed.nc"
```

Then define...

```
minp = 990.
target_rmw = 150000.
cen_lat = 22.75
cen_lon = 268.0
```

where minp is the target pressure of the vortex seed (in hPa), target_rmw is the target radius of maximum wind (in m), and cen_lat, cen_lon is the vortex location in the SE data.

```
$> ncl seed-tc-in-ncdata.ncl
```

should then produce an analysis with seeded TC.

Use this as ncdata in Betacast.

## To remove an existing vortex from SE analysis state

Find the rough location of the vortex's PSL center

`find-tc-fill-params.ncl` will attempt to find the correct vortex shape/settings to optimally remove the vortex from the analysis.

Edit the head of the file...

```
inic_file="/glade/u/home/zarzycki/work/sewx/INIC/sample_ne120.nc"
deltaMax=0.25            ; user-defined nominal resolution of raw data
psminlat=21.75           ; initial target for psminlat
psminlon=-53.5           ; initial target for psminlon
truePS_scan_radius=300.  ; what is radius to scan for true psmin relative to psminlat/psminlon? In km.
rad_for_corr=800.        ; what radius (in km) are we using for radial averages?
modify_q = True          ; do we modify q in addition to T, PS, U, V?
modify_q_mult = 2.5      ; multiplier on Q relative to C-C scaling (1.0 means scale mixing ratio with RH% and T).
gamma_   = 0.0065        ; lapse rate
```

Note that deltaMax is in degrees and needs to be user-specified.

Once edited, run:

```
$> ncl find-tc-fill-params.ncl
```

This should take a little while scanning all the potential settings that could "fill in" the vortex. At the end of the script, it will print...

```
(0)	BEST SETTINGS --------------------------
(0)	rp: 179761
(0)	dp: 1641.97
(0)	zp: 14305.1
(0)	gamma_: 0.0065
(0)	exppr: 1.12118
(0)	cen_lat: 21.87653196344101
(0)	cen_lat: 306.2441153488365
(0)	modify_q: True
(0)	modify_q_mult: 2.5
(0)	************* --------------------------
```

Transfer these settings over to `seed-tc-in-ncdata.ncl` in the relevant block. Also make sure you toggle `invert_vortex`.

```
invert_vortex = True
origfile="???"
seedfile="???_seed.nc"
```

Now you can run the script and it will fill in the existing vortex in the analysis.

```
$> ncl seed-tc-in-ncdata.ncl
```

Use this as ncdata in Betacast.

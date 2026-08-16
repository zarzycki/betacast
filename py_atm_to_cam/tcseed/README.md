## To add a vortex to an SE analysis state

A user can either "seed" a vortex or "deseed" a vortex. To do so, two settings need to be added to the main namelist and a "seed namelist" must be generated.

In the main namelist, add...

```
  add_vortex = true
  vortex_namelist = /glade/u/home/zarzycki/betacast/namelists/seed.default.nl
```

Then generate a "seed namelist" with the below settings.

| Namelist Variable | Type | Description | RQ_ADD | RQ_RM | Default |
| --- | --- | --- | --- | --- | --- |
| invert_vortex | Bool | **False** means we are adding (seeding), **True** means remove existing vortex ("deseed") | X | X | |
| psminlat | Numeric | Latitude of vortex center | X | X^ | |
| psminlon | Numeric | Longitude of vortex center | X | X^ | |
| deltaMax | Numeric | Nominal resolution of input data (degrees) |  | X | |
| minp | Numeric | If **positive**: target central pressure of vortex to add (hPa). If **between -19999 and 0**: magnitude is interpreted as the pressure deficit *dp* relative to ambient (Pa) (e.g., -5000. = 50 hPa deficit; note -1. means a 1 Pa deficit, *not* the default!). If **< -19999**: use default of 995 hPa. | X | | 995.0 |
| target_rmw | Numeric | Target radius of maximum wind of vortex to add (m). Values below 15000. are capped to 15 km; negative uses the default. | X* |  | 200000.0 |
| zp | Numeric | Height for calculation of P (~top of vortex) (m) | X* |  | 12000. |
| exppr | Numeric | Exponent for r dependence of p | X* |  | 1.5 |
| gamma | Numeric | TC environmental lapse rate (K/m) | X* | X* | 0.0065 |
| modify_q | Bool | Should the q profile be modified via RH%? | X | X | |
| modify_q_mult | Numeric | Multiplicative factor beyond Clausius-Clapeyron if modify_q = True | X* | X* | 1.0 |
| restart_file | Bool | Is this a CAM-SE restart file (true) or ncdata (false)? | X | X | |

An 'X' in the column RQ_ADD means that namelist variable is **required** if one is *adding* a vortex, RQ_RM means required if one is *removing* a vortex.

*Any required variable with an asterisk can be set to -1.0 for the default setting in the right column. **This does not apply to minp**: because negative minp values are interpreted as a pressure deficit in Pa (see table), setting minp=-1. produces a 1 Pa vortex. To get the 995 hPa default, set minp to any value < -19999 (e.g., -99999.).

^For vortex removal, psminlat and psminlon are "best guesses" that must be within 3 great circle degrees of the observed TC and will be overwritten by the code which finds the local minimum in PS in the model data.

### Adding a vortex example

If one is seeding (adding) a vortex, the file **must** contain the following...

```
invert_vortex=False
psminlat = 25.315
psminlon = 282.831
# 970 hPa target central pressure; use a value < -19999 for the 995 hPa default
# (do NOT use -1., which means a 1 Pa deficit -- see minp conventions above)
minp = 970.
target_rmw = -1.
gamma=-1.
zp=-1.
exppr=-1.
modify_q=True
modify_q_mult=-1.
```

### Removing a vortex example

If one is removing a vortex, the file **must** contain the following...

```
invert_vortex=True
psminlat=28.2
psminlon=312.8
deltaMax=0.25
gamma=-1.
modify_q=True
modify_q_mult=2.5
```

`find-tc-fill-params.ncl` will attempt to find the correct vortex shape/settings to optimally remove the vortex from the analysis. It will append these settings to the namelist without user intervention and then run the seeding code with the invert_vortex flag specified as True.

### Seeding/unseeding multiple TCs

A single vortex namelist may operate on multiple TCs. Separate each TC block with a line containing only `---`:

```
invert_vortex=False
psminlat = 25.315
psminlon = 282.831
minp = 970.
target_rmw = -1.
gamma=-1.
zp=-1.
exppr=-1.
modify_q=True
modify_q_mult=-1.
---
invert_vortex=True
psminlat=28.2
psminlon=312.8
deltaMax=0.25
gamma=-1.
modify_q=True
modify_q_mult=2.5
```

Each block is self-contained and follows the same rules as a single-TC namelist above (blocks may mix seeding and unseeding); keys not set in a block fall back to their defaults, except `psminlat` and `psminlon`, which are required in every block (the parser exits with an error if either is missing). A file with no `---` separators behaves exactly as a single-TC namelist. See `namelists/seed-unseed.multi.example.nl` for a working example.

Note: multi-TC namelists are only supported by the main Python pipeline (`atm_to_cam.py`, i.e., `standalone_vortex = false`). The standalone tools (`py_atm_to_cam/tcseed/py-seed-tc-in-ncdata.py` and the legacy NCL `seed-tc-in-ncdata.ncl`) handle single-TC namelists only; the Python standalone tool will exit with an error if given a file containing `---` separators.
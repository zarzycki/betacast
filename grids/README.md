# Grid information

## Creating physgrid SCRIP grids for E3SM/SCREAM

TempestRemap can generate Exodus and SCRIP files, even for "internal" grids in HOMME by specifying `--alt`.

```
GenerateCSMesh --alt --res 30 --file ./ne30.g
GenerateVolumetricMesh --in ./ne30.g --out ./ne30pg2.g --np 2 --uniform
ConvertMeshToSCRIP --in ./ne30pg2.g --out ./ne30pg2_scrip.nc
```

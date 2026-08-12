# Grid information

## Creating physgrid SCRIP grids for E3SM/SCREAM

TempestRemap can generate Exodus and SCRIP files, even for "internal" grids in HOMME by specifying `--alt`.

```
NE=30
PG=3
GenerateCSMesh --alt --res ${NE} --file ./ne${NE}.g
GenerateVolumetricMesh --in ./ne${NE}.g --out ./ne${NE}pg${PG}.g --np ${PG} --uniform
ConvertMeshToSCRIP --in ./ne${NE}pg${PG}.g --out ./ne${NE}pg${PG}_scrip.nc
```

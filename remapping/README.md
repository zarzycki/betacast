# Gridding notes

### Generate ESMF from SCRIP file

``
The ESMF_Scrip2Unstruct application is a parallel program that converts a SCRIP format grid file into an unstructured grid file in the ESMF unstructured file format or in the UGRID file format. This application program can be used together with ESMF_RegridWeightGen application for the unstructured SCRIP format grid files. An unstructured SCRIP grid file will be converted into the ESMF unstructured file format internally in ESMF_RegridWeightGen. The conversion subroutine used in ESMF_RegridWeightGen is sequential and could be slow if the grid file is very big. It will be more efficient to run the ESMF_Scrip2Unstruct first and then regrid the output ESMF or UGRID file using ESMF_RegridWeightGen. Note that a logically rectangular grid file in the SCRIP format (i.e. the dimension grid_rank is equal to 2) can also be converted into an unstructured grid file with this application.

The application usage is as follows:

ESMF_Scrip2Unstruct  inputfile outputfile dualflag [fileformat]

where
  inputfile       - a SCRIP format grid file

  outputfile      - the output file name

  dualflag        - 0 for straight conversion and 1 for dual
		    mesh.  A dual mesh is a mesh constructed
                    by putting the corner coordinates in the
                    center of the elements and using the
		    center coordinates to form the mesh
		    corner vertices.

  fileformat      - an optional argument for the output file
		    format.  It could be either ESMF or UGRID.
                    If not specified, the output file is in
		    the ESMF format.
```

```
ESMF_Scrip2Unstruct /glade/campaign/cesm/cesmdata/inputdata/lnd/clm2/mappingdata/grids/SCRIPgrid_0.125x0.125_nomask_c170126.nc tmp.nc 0 ESMF
```


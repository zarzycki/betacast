#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 4                      # yellowstone setting
#BSUB -o out.%J                  # output filename
#BSUB -e out.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J sub_ncl 
#BSUB -W 23:58                    # wall clock limit
#BSUB -P P54048000               # account number

################################################################

date

ncl digifilter.ncl machineid=1 'filtfile_name = "haiyan_48_x8.cam.h0.2013-11-03-00000.nc"' 'gridname = "haiyan_48_x8"'

date 

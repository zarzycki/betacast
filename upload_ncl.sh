#!/bin/bash

##=======================================================================
#BSUB -a poe                     # use LSF openmp elim
#BSUB -N
#BSUB -n 1                      # yellowstone setting
#BSUB -o ncl_up.%J                  # output filename
#BSUB -e ncl_up.%J                  # error filename
#BSUB -q geyser                 # queue
#BSUB -J upload_ncl
#BSUB -W 10:00                    # wall clock limit
#BSUB -P P54048000               # account number
##=======================================================================
module load ncl
###############################################################################
yearstr=2017
monthstr=11
daystr=17
cyclestr=12
cyclestrsec=43200
###############################################################################
#nclweightfile="/glade/p/work/zarzycki/maps/forecast_plot_maps/map_natlantic_30_x4_to_0.25x0.25glob_bilinear.nc"
nclweightfile=$1
runDir=$2
path_to_ncl=/glade/u/home/zarzycki/sewx-cam-forecast/plotting_ncl/
htmlFolder=/glade/u/home/zarzycki/sewx-cam-forecast/html_for_upload/
twodaysago=20120821

echo "SSHings..."
ssh zarzycki@burnt.cgd.ucar.edu "mkdir -p /web/web-data/staff/zarzycki/current/${yearstr}${monthstr}${daystr}${cyclestr} ; \ 
		 cd /web/web-data/staff/zarzycki/current/ ; \
		 cp *cfg *html ${yearstr}${monthstr}${daystr}${cyclestr} ; \
		 rm ${yearstr}${monthstr}${daystr}${cyclestr}/index.html "

## UPDATE html page
cd ${htmlFolder}
if [ ! -f index.html ]; then
        echo "copying index and sedding in this date"
	cp -v index.HOLD index.html
	sed '/<!--FORECASTHEAD-->/ r htmltemplate.html' index.html > _index1.html
	sed -e "/$twodaysago${cyclestr}/d" _index1.html > _index2.html
	sed -e "s/YYYYMMDDHH/${yearstr}${monthstr}${daystr}${cyclestr}/" _index2.html > _index3.html
	mv _index3.html index.html
	rm _index*.html
	scp index.html zarzycki@burnt.cgd.ucar.edu:/web/web-data/staff/zarzycki/current
fi
  
	filenames=`ls ${runDir}/*h0*.nc`
	numfiles=`ls ${runDir}/*h0*.nc | wc -l`
	echo $numfiles

  VARS=PRECLav,PRECCav

  for i in `seq 1 ${numfiles}`;
  do
    thisFile=`echo $filenames | cut -d" " -f${i}`
    if [ "$i" -eq 1 ] ; then
      ncks -v ${VARS} ${thisFile} ${runDir}/sum${i}.nc
      ncks -A -v ${VARS} ${runDir}/sum${i}.nc ${thisFile}
      rm ${runDir}/sum${i}.nc
    else
      iminus1=`expr $i - 1`
      lastFile=`echo $filenames | cut -d" " -f${iminus1}`
      ncrcat -v ${VARS} ${thisFile} ${lastFile} ${runDir}/tmpfile2.nc
      ncra -h -O -y ttl ${runDir}/tmpfile2.nc ${runDir}/sum${i}.nc
      ncap2 -A -s 'time=time+0.0625' ${runDir}/sum${i}.nc ${runDir}/sum${i}.nc
      ncks -A -v ${VARS} ${runDir}/sum${i}.nc ${thisFile}
      rm ${runDir}/sum${i}.nc ${runDir}/tmpfile2.nc
    fi
  done 
		
		sleep 5
		echo "Found at least one file"
		echo $filenames
		cd ${runDir}
		for i in *h0*.nc; do mv $i _$i; done
		cd ${htmlFolder}
		newfiles=`ls ${runDir}/_*h0*-00000.nc ${runDir}/_*h0*-21600.nc ${runDir}/_*h0*-43200.nc ${runDir}/_*h0*-64800.nc`
		for f in $newfiles
		do
			echo "Processing $f"
 			ncl ${path_to_ncl}/weatherplot.ncl inisec=$cyclestrsec iniday=$daystr inimon=$monthstr iniyear=$yearstr 'filename="'${f}'"' 'wgt_file="'${nclweightfile}'"' > ncl.output 2>&1
 			if [ grep FileReadVar ncl.output ]; then
 				sleep 5
 				echo "Found an error"
 				ncl ${path_to_ncl}/weatherplot.ncl inisec=$cyclestrsec iniday=$daystr inimon=$monthstr iniyear=$yearstr 'filename="'${f}'"' 'wgt_file="'${nclweightfile}'"'
 			fi
 			#rm ncl.output
		done
		
		## Trim whitespace around pngs
		genpngs=`ls *png`
		for g in $genpngs
		do
			convert -trim +repage $g $g
		done
		
		## Get file list in txt file for FLANIS viewer
		basins=(natl epac)
		outflds=(wind tmq flut prect sumprect 500vort shear850250)

		# Loop over items in outflds
		for basin in ${basins[*]}
		do
			for item in ${outflds[*]}
			do
					#printf "   %s\n" $item
					
					if [ ! -f ${item}_${basin}_files.txt ]; then
						echo "File not found!"
						ls -1 ${item}*${basin}*png > ${item}_${basin}_files_nopath.txt
					else
						mv ${item}_${basin}_files_nopath.txt orig_${item}_${basin}_files.txt
						ls -1 ${item}*${basin}*png > tocat_${item}_${basin}_files.txt
						cat orig_${item}_${basin}_files.txt tocat_${item}_${basin}_files.txt > ${item}_${basin}_files_nopath.txt
						rm tocat_${item}_${basin}_files.txt orig_${item}_${basin}_files.txt
					fi
					
					cp -v ${item}_${basin}_files_nopath.txt ${item}_${basin}_files.txt
					#sed -i.bak 's/^/##/' ${item}files.txt        
				
			done
		done

    # Copy tracker data over if around...
    mv -v /glade/u/home/zarzycki/tempest-scripts/forecast/trajs.trajectories.txt.NATLANTIC30X4-CYCLONES.png .
    cp -v /glade/u/home/zarzycki/tempest-scripts/forecast/fin-atcf/atcf.tempest.${yearstr}${monthstr}${daystr}${cyclestr} .
		
		## Move files to server
		## Google create remote directory if not existant
		## use rysnc?
		echo "Moving files to remote server"
		scp *.png *.txt atcf.tempest* zarzycki@burnt.cgd.ucar.edu:/web/web-data/staff/zarzycki/current/${yearstr}${monthstr}${daystr}${cyclestr}
		#mv $newfiles $procdir

mkdir ${yearstr}${monthstr}${daystr}${cyclestr}
mv *.png *.txt atcf.tempest* ${yearstr}${monthstr}${daystr}${cyclestr}

### Finish uploading index.html
printtime=`date -u`
sed -e 's/\"red/\"green/' index.html > _index4.html
sed -e "s/CURRENTLY UPDATING/COMPLETED AT $printtime/" _index4.html > _index5.html
mv _index5.html index.html
rm _index*html
scp index.html zarzycki@burnt.cgd.ucar.edu:/web/web-data/staff/zarzycki/current
mv -v index.html index.HOLD

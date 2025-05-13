# Remapping notes

### Interface for generating maps from command line

```
for i in $(seq 5 14); do
  IX=$(printf "%03d" $i)
  echo $IX
  python gen_analysis_to_model_wgt_file.py \
      --ANLGRID "era5_0.25x0.25" \
      --DSTGRIDNAME "mpasa3-60-tclf${IX}_scrip" \
      --DSTGRIDFILE "/glade/u/home/zarzycki/work/grids/scrip/mpasa3-60-tclf${IX}_scrip.nc" \
      --ANLGRIDPATH "../grids/anl_scrip/" \
      --WGTFILEDIR "/glade/work/zarzycki/sewx/mapping/"
done
```
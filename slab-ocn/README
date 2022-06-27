# Slab ocean w/ empirical mixing

Citation: **C. M. Zarzycki**. Tropical cyclone intensity errors associated with lack of two-way ocean coupling in high-resolution global simulations. Journal of Climate, 290 (23):0 8589--8610, 2016. [10.1175/JCLI-D-16-0273.1](http://dx.doi.org/10.1175/JCLI-D-16-0273.1).

### Creating bias correction fields

On the slab data stream, one could add a monthly bias correction field under the variable "qdp." This was used in Zarzycki (2016) to improve the correlation of the slab ocean climatology with the fixedSST climatology.

To generate:

1. Run an F simulation (>= 1 year, used 10 years in Z16) with the prescribed SST file (i.e., FIXEDSST), outputting monthly averages.
2. Run an F simulation (>= 1 year, used 10 years in Z16) with the slab ocean turned on (i.e., SLAB) (with qdp set to 0.0).
3. (optional) if being lazy and using linear interp'ed files, clean up the SST field for reasonableness with `ncap2 -O -s where(SST < 271.34 && SST > 0.0001) SST=271.34`
4. Create a monthly climatology of the difference `BIAS = SLAB - FIXEDSST`.
5. Add the monthly climatology to the SOM stream as "qdp".
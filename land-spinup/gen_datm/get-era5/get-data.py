import cdsapi
import sys

# python get-data.py 1986 02 /glade/scratch/zarzycki/ERA5-DATM/
print ("This is the name of the script: ", sys.argv[0])
print ("Number of arguments: ", len(sys.argv))
print ("The arguments are: " , str(sys.argv))

c = cdsapi.Client()

OUTYEAR=sys.argv[1]
OUTMONTH=sys.argv[2]
OUTDIR=sys.argv[3]

#OUTVAR = OUTDIR + '/out.' + OUTYEAR + '.'+ OUTMONTH + '.nc'
OUTVAR = OUTDIR + '/out.' + OUTYEAR + '.'+ OUTMONTH + '.zip'

dataset = "reanalysis-era5-single-levels"
request = {
    "product_type": ["reanalysis"],
    "variable": [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
        "2m_dewpoint_temperature",
        "2m_temperature",
        "mean_total_precipitation_rate",
        "surface_pressure",
        "surface_solar_radiation_downwards",
        "surface_thermal_radiation_downwards"
    ],
    "year": [OUTYEAR],
    "month": [OUTMONTH],
    "day": [
        "01", "02", "03",
        "04", "05", "06",
        "07", "08", "09",
        "10", "11", "12",
        "13", "14", "15",
        "16", "17", "18",
        "19", "20", "21",
        "22", "23", "24",
        "25", "26", "27",
        "28", "29", "30",
        "31"
    ],
    "time": [
            "00:00",
            "03:00",
            "06:00",
            "09:00",
            "12:00",
            "15:00",
            "18:00",
            "21:00"
    ],
    "data_format": "netcdf",
    "download_format": "unarchived"
}
target = OUTVAR

client = cdsapi.Client()
client.retrieve(dataset, request, target).download()
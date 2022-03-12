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

OUTVAR = OUTDIR + '/out.' + OUTYEAR + '.'+ OUTMONTH + '.nc'

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'format': 'netcdf',
        'variable': [
            '10m_u_component_of_wind', '10m_v_component_of_wind', '2m_dewpoint_temperature',
            '2m_temperature', 'mean_total_precipitation_rate', 'surface_pressure',
            'surface_solar_radiation_downwards', 'surface_thermal_radiation_downwards',
#            'surface_solar_radiation_downwards',
        ],
        'year': OUTYEAR,
        'month': OUTMONTH,
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
#         'time': [
#             '00:00', '01:00', '02:00',
#             '03:00', '04:00', '05:00',
#             '06:00', '07:00', '08:00',
#             '09:00', '10:00', '11:00',
#             '12:00', '13:00', '14:00',
#             '15:00', '16:00', '17:00',
#             '18:00', '19:00', '20:00',
#             '21:00', '22:00', '23:00',
#         ],
        'time': [
            '00:00',
            '03:00',
            '06:00',
            '09:00',
            '12:00',
            '15:00',
            '18:00',
            '21:00',
        ],
    },
    OUTVAR)

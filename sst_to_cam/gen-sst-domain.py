import numpy as np
import argparse
import logging
import os
from datetime import datetime
from netCDF4 import Dataset

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def parse_args():
    parser = argparse.ArgumentParser(description='Generate SST domain and SCRIP grid files')
    parser.add_argument('--inputres', type=str, required=True,
                        help='Resolution string, e.g. "180x360" (nlatxnlon)')
    return parser.parse_args()


def gen_domain_file(filename, lat, lon, inputres):
    """Write domain file with yc (lat) and xc (lon) variables."""
    ds = Dataset(filename, 'w', format='NETCDF4')

    ds.createDimension('nj', len(lat))
    ds.createDimension('ni', len(lon))

    yc = ds.createVariable('yc', 'f8', ('nj',))
    yc.long_name = 'latitude'
    yc.units = 'degrees_north'
    yc[:] = lat

    xc = ds.createVariable('xc', 'f8', ('ni',))
    xc.long_name = 'longitude'
    xc.units = 'degrees_east'
    xc[:] = lon

    ds.title = 'Betacast-generated DOCN domain file'
    ds.resolution = inputres
    ds.creation_date = datetime.now().strftime('%c')

    ds.close()
    logging.info(f"Wrote domain file: {filename}")


def gen_scrip_file(filename, lat, lon, nlat, nlon):
    """Write SCRIP grid description file for a rectilinear grid."""
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon
    grid_size = nlat * nlon

    # Build 2D center arrays (lat varies slowest, lon varies fastest)
    lon2d, lat2d = np.meshgrid(lon, lat)
    center_lat = lat2d.ravel()
    center_lon = lon2d.ravel()

    # Corner offsets (half-cell)
    dlatd2 = dlat / 2.0
    dlond2 = dlon / 2.0

    # Corners: SW, SE, NE, NW
    corner_lat = np.zeros((grid_size, 4))
    corner_lon = np.zeros((grid_size, 4))

    corner_lat[:, 0] = center_lat - dlatd2  # SW lat
    corner_lat[:, 1] = center_lat - dlatd2  # SE lat
    corner_lat[:, 2] = center_lat + dlatd2  # NE lat
    corner_lat[:, 3] = center_lat + dlatd2  # NW lat

    corner_lon[:, 0] = center_lon - dlond2  # SW lon
    corner_lon[:, 1] = center_lon + dlond2  # SE lon
    corner_lon[:, 2] = center_lon + dlond2  # NE lon
    corner_lon[:, 3] = center_lon - dlond2  # NW lon

    # Clamp corner lats to [-90, 90]
    corner_lat = np.clip(corner_lat, -90.0, 90.0)

    # Write SCRIP file
    ds = Dataset(filename, 'w', format='NETCDF4')

    ds.createDimension('grid_size', grid_size)
    ds.createDimension('grid_corners', 4)
    ds.createDimension('grid_rank', 2)

    v = ds.createVariable('grid_center_lat', 'f8', ('grid_size',))
    v.units = 'degrees'
    v[:] = center_lat

    v = ds.createVariable('grid_center_lon', 'f8', ('grid_size',))
    v.units = 'degrees'
    v[:] = center_lon

    v = ds.createVariable('grid_corner_lat', 'f8', ('grid_size', 'grid_corners'))
    v.units = 'degrees'
    v[:] = corner_lat

    v = ds.createVariable('grid_corner_lon', 'f8', ('grid_size', 'grid_corners'))
    v.units = 'degrees'
    v[:] = corner_lon

    v = ds.createVariable('grid_dims', 'i4', ('grid_rank',))
    v[:] = np.array([nlon, nlat], dtype=np.int32)

    v = ds.createVariable('grid_imask', 'i4', ('grid_size',))
    v.units = 'unitless'
    v[:] = np.ones(grid_size, dtype=np.int32)

    ds.title = f'rectilinear_to_SCRIP ({nlat},{nlon})'
    ds.Conventions = 'SCRIP'
    ds.date_created = datetime.now().strftime('%c')

    ds.close()
    logging.info(f"Wrote SCRIP file: {filename}")


def main():
    args = parse_args()
    inputres = args.inputres

    logging.info(f"GEN_SST_DOMAIN: inputres: {inputres}")

    nlat, nlon = [int(x) for x in inputres.split('x')]
    dlat = 180.0 / nlat
    dlon = 360.0 / nlon

    logging.info(f"GEN_SST_DOMAIN: nlat: {nlat}  nlon: {nlon}")
    logging.info(f"GEN_SST_DOMAIN: dlat: {dlat}  dlon: {dlon}")

    dlatd2 = dlat / 2.0
    dlond2 = dlon / 2.0

    lat = np.linspace(-90.0 + dlatd2, 90.0 - dlatd2, nlat)
    lon = np.linspace(dlond2, 360.0 - dlond2, nlon)

    # Output to grids/domains/ relative to BETACAST root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    betacast_dir = os.path.dirname(script_dir)
    storedir = os.path.join(betacast_dir, 'grids', 'domains')
    os.makedirs(storedir, exist_ok=True)

    domain_file = os.path.join(storedir, f'domain.ocn.{nlat}x{nlon}.nc')
    scrip_file = os.path.join(storedir, f'scrip.ocn.{nlat}x{nlon}.nc')

    gen_domain_file(domain_file, lat, lon, inputres)
    gen_scrip_file(scrip_file, lat, lon, nlat, nlon)

    logging.info("GEN_SST_DOMAIN: done")


if __name__ == '__main__':
    main()

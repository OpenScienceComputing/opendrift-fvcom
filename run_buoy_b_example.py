"""
Example: 3-day OpenDrift surface-drift simulation using the UMASSD FVCOM
GOM3/GOM4 hindcast (read from a public Icechunk virtual store), releasing
a 10-particle cluster at the NERACOOS Buoy B01 (Western Maine Shelf)
location.

The reference store currently spans 1978-01-01 through 2017-01-20, so the
most recent complete year is 2016; this example uses 2016-06-01 00:00
through 2016-06-04 00:00 (3 days, hourly steps). Adjust START/END (and the
icechunk bucket/prefix) for other stores or periods.
"""
import os
import time as _time

import numpy as np
import pandas as pd
import xarray as xr
import icechunk

from fvcom_icechunk_reader import Reader as FVCOMIcechunkReader

from opendrift.models.oceandrift import OceanDrift

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(HERE, 'output')
os.makedirs(OUTDIR, exist_ok=True)

# Public UMASSD FVCOM GOM3 hindcast icechunk store
BUCKET = 'umassd-fvcom'
PREFIX = 'gom3/hindcast/icechunk/gom3-whole-nochunk.icechunk'
REGION = 'us-east-1'

# NERACOOS Buoy B01 (Western Maine Shelf): site_latitude/site_longitude
BUOY_B_LAT = 43.1807
BUOY_B_LON = -70.4278

START = pd.Timestamp('2016-06-01 00:00:00')
END = pd.Timestamp('2016-06-04 00:00:00')

N_PARTICLES = 10
CLUSTER_RADIUS_M = 500

CACHE_FN = os.path.join(OUTDIR, 'surface_uv_cache.nc')


def open_store():
    config = icechunk.RepositoryConfig.default()
    config.set_virtual_chunk_container(icechunk.VirtualChunkContainer(
        url_prefix=f's3://{BUCKET}/',
        store=icechunk.s3_store(region=REGION, anonymous=True),
    ))
    credentials = icechunk.containers_credentials({
        f's3://{BUCKET}/': icechunk.s3_credentials(anonymous=True)
    })
    storage = icechunk.s3_storage(bucket=BUCKET, prefix=PREFIX, region=REGION, anonymous=True)
    repo = icechunk.Repository.open(storage, config, authorize_virtual_chunk_access=credentials)
    ds = xr.open_zarr(repo.readonly_session('main').store, consolidated=False, chunks=None)
    return ds


def load_surface_currents():
    if os.path.exists(CACHE_FN):
        print(f'Loading cached surface u/v from {CACHE_FN}')
        return xr.load_dataset(CACHE_FN)

    print('Opening icechunk store...')
    ds = open_store()

    tvals = pd.to_datetime(ds['time'].values)
    mask = (tvals >= START) & (tvals <= END)
    idx = np.where(mask)[0]
    print(f'Time window: {tvals[idx[0]]} .. {tvals[idx[-1]]}  ({len(idx)} steps)')

    sub = ds.isel(time=slice(idx[0], idx[-1] + 1))[
        ['u', 'v', 'lon', 'lat', 'lonc', 'latc']
    ]

    t0 = _time.time()
    print('Loading surface u/v into memory (streams from S3)...')
    sub = sub.isel(siglay=0).load()
    print(f'  loaded in {_time.time() - t0:.1f} s')

    sub.to_netcdf(CACHE_FN)
    print(f'  cached to {CACHE_FN}')
    return sub


def main():
    sub = load_surface_currents()

    reader = FVCOMIcechunkReader(sub, name='FVCOM-GOM-icechunk')
    print(reader)

    o = OceanDrift(loglevel=20)
    o.add_reader(reader)
    o.set_config('general:use_auto_landmask', True)
    o.set_config('drift:vertical_mixing', False)

    o.seed_elements(
        lon=BUOY_B_LON, lat=BUOY_B_LAT,
        radius=CLUSTER_RADIUS_M,
        number=N_PARTICLES,
        z=0,
        time=START.to_pydatetime(),
    )

    o.run(
        end_time=END.to_pydatetime(),
        time_step=3600,          # 1 hour, matches data frequency
        time_step_output=3600,
        outfile=os.path.join(OUTDIR, 'buoy_b_trajectories.nc'),
    )

    print(o)

    fig_fn = os.path.join(OUTDIR, 'buoy_b_trajectories.png')
    o.plot(filename=fig_fn, buffer=0.3, fast=True)
    print('Wrote', fig_fn)

    try:
        gif_fn = os.path.join(OUTDIR, 'buoy_b_animation.gif')
        o.animation(filename=gif_fn, buffer=0.3, fast=True)
        print('Wrote', gif_fn)
    except Exception as e:
        print('Animation failed (non-fatal):', e)


if __name__ == '__main__':
    main()

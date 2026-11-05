"""
A lightweight OpenDrift reader for the UMASSD FVCOM GOM3/GOM4 hindcast,
served from a public Icechunk virtual store, using surface (top sigma
layer) currents only.

Unlike opendrift.readers.reader_FVCOM_xarray.Reader (which opens its own
zarr/netcdf files and assumes 1-D siglay/siglev coordinates), this reader
is handed an already-opened, already-subset xarray.Dataset (e.g. sliced to
a short time window and loaded into memory) and uses the dataset's own
geographic lon/lat/lonc/latc coordinates directly, sidestepping the
Cartesian (NAD83 state-plane) projection and the 2-D (siglay, node) sigma
coordinate entirely.
"""
import numpy as np
import pandas as pd

from opendrift.readers.basereader.unstructured import UnstructuredReader


class Reader(UnstructuredReader):

    node_variables = []
    face_variables = ['x_sea_water_velocity', 'y_sea_water_velocity']

    variable_mapping = {
        'x_sea_water_velocity': 'u',
        'y_sea_water_velocity': 'v',
    }

    def __init__(self, ds, name='FVCOM-icechunk-surface'):
        """
        ds: xarray.Dataset already opened from the icechunk store (and
            ideally already .load()-ed into memory for the desired time
            window), with data_vars u, v already reduced to (time, nele)
            (i.e. a single, surface, sigma layer already selected) and
            coordinates lon/lat (node), lonc/latc (nele).
        """
        self.name = name
        self.dataset = ds

        # Use geographic coordinates directly - no projection needed.
        self.proj4 = '+proj=latlong +datum=WGS84'

        self.x = np.asarray(ds['lon'].values)
        self.y = np.asarray(ds['lat'].values)
        self.xc = np.asarray(ds['lonc'].values)
        self.yc = np.asarray(ds['latc'].values)

        times = pd.to_datetime(ds['time'].values).to_pydatetime()
        self.times = np.array(list(times))
        self.start_time = self.times[0]
        self.end_time = self.times[-1]
        if len(self.times) > 1:
            self.time_step = self.times[1] - self.times[0]
        else:
            self.time_step = None

        self.xmin, self.xmax = float(self.x.min()), float(self.x.max())
        self.ymin, self.ymax = float(self.y.min()), float(self.y.max())
        self.zmin, self.zmax = -np.inf, 0

        self.variables = list(self.variable_mapping.keys())

        # Plain numpy arrays, shape (time, nele), so get_variables() does
        # no further xarray indexing (ds should already be .load()-ed).
        self._u = np.asarray(ds['u'].values)
        self._v = np.asarray(ds['v'].values)

        super().__init__()

        self.boundary = self._build_boundary_polygon_(self.x, self.y)
        self.nodes_idx = self._build_ckdtree_(self.x, self.y)
        self.faces_idx = self._build_ckdtree_(self.xc, self.yc)

    def get_variables(self, requested_variables, time=None, x=None, y=None, z=None):
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        z = np.atleast_1d(z)
        if len(z) == 1:
            z = z[0] * np.ones(x.shape)

        requested_variables, time, x, y, z, _outside = self.check_arguments(
            requested_variables, time, x, y, z)

        (nearest_time, _time_before, _time_after, indx_nearest,
         _indx_before, _indx_after) = self.nearest_time(time)

        neles = self._nearest_face_(x, y)

        variables = {}
        for var in requested_variables:
            if var == 'x_sea_water_velocity':
                field = self._u
            elif var == 'y_sea_water_velocity':
                field = self._v
            else:
                raise ValueError(f'Variable not available: {var}')
            variables[var] = field[indx_nearest, neles]

        return variables

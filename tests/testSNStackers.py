from builtins import zip
import numpy as np
import unittest
import lsst.utils.tests
import stackers.coadd_stacker as sn_stacker

m5_ref = dict(
    zip('ugrizy', [23.60, 24.83, 24.38, 23.92, 23.35, 22.44]))


class TestSNStackers(unittest.TestCase):

    def testCoaddStacker(self):
        """Test the SN cadence metric """

        band = 'r'
        # Define fake data
        names = ['observationStartMJD', 'fieldRA', 'fieldDec',
                 'fiveSigmaDepth', 'visitExposureTime', 'numExposures', 'visitTime', 'season']
        types = ['f8']*len(names)
        names += ['night']
        types += ['i2']
        names += ['filter']
        types += ['O']

        day0 = 59000
        daylast = day0+0.0069*5
        cadence = 0.0069
        dayobs = np.arange(day0, daylast, cadence)
        npts = len(dayobs)
        data = np.zeros(npts, dtype=list(zip(names, types)))
        data['observationStartMJD'] = dayobs
        # data['night'] = np.floor(data['observationStartMJD']-day0)
        data['night'] = 10
        data['fiveSigmaDepth'] = m5_ref[band]
        data['visitExposureTime'] = 15.
        data['numExposures'] = 2
        data['visitTime'] = 2.*15.
        data['filter'] = band
        data['season'] = 1.

        # Run the stacker with these fake data

        stacker = sn_stacker.CoaddStacker()

        result = stacker._run(data)

        # The result should be equal to...
        result_ref = np.array([(59000.01725, 0., 0., 25.35268906,
                                90., 12., 180., 1., 10, 'r')], dtype=[('observationStartMJD', '<f8'), ('fieldRA', '<f8'), ('fieldDec', '<f8'), ('fiveSigmaDepth', '<f8'), ('visitExposureTime', '<f8'), ('numExposures', '<f8'), ('visitTime', '<f8'), ('season', '<f8'), ('night', 'i2'), ('filter', 'O')])

        for var in ['observationStartMJD', 'fieldRA', 'fieldDec', 'fiveSigmaDepth',
                    'visitExposureTime', 'numExposures', 'visitTime', 'season', 'night']:
            assert(np.isclose(np.asarray(
                result[var]), np.asarray(result_ref[var])))

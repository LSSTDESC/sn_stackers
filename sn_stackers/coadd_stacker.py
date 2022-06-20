from lsst.sims.maf.stackers import BaseStacker
import numpy as np
import numpy.lib.recfunctions as rf
import multiprocessing
from astropy.table import Table, Column
import pandas as pd
import time
import operator

__all__ = ['CoaddStacker']


class CoaddStacker(BaseStacker):
    """
    Class to coadd observations per night
    """

    colsAdded = ['sn_coadd']
    """
    def __init__(self, mjdCol='observationStartMJD', RACol='fieldRA', DecCol='fieldDec', m5Col='fiveSigmaDepth', nightCol='night', filterCol='filter', numExposuresCol='numExposures', visitTimeCol='visitTime', visitExposureTimeCol='visitExposureTime', seeingaCol='seeingFwhmEff', seeingbCol='seeingFwhmGeom', airmassCol='airmass', skyCol='sky', moonCol='moonPhase', nproc=8):
        self.colsReq = [mjdCol, RACol, DecCol, m5Col, filterCol, nightCol,
                        numExposuresCol, visitTimeCol, visitExposureTimeCol,
                        seeingaCol, seeingbCol, airmassCol, moonCol]
        self.mjdCol = mjdCol
        self.RACol = RACol
        self.DecCol = DecCol
        self.nightCol = nightCol
        self.filterCol = filterCol
        self.m5Col = m5Col
        self.numExposuresCol = numExposuresCol
        self.visitTimeCol = visitTimeCol
        self.visitExposureTimeCol = visitExposureTimeCol
        self.seeingaCol = seeingaCol
        self.seeingbCol = seeingbCol
        self.airmassCol = airmassCol
        self.skyCol = skyCol
        self.moonCol = moonCol

        self.units = ['int']
        self.nproc = 1
    """

    def __init__(self, col_sum=['numExposures', 'visitTime', 'visitExposureTime'],
                 col_mean=['observationStartMJD', 'fieldRA', 'fieldDec',
                           'fiveSigmaDepth', 'pixRA', 'pixDec', 'healpixID', 'season'],
                 col_median=['airmass', 'sky', 'moonPhase',
                             'seeingFwhmEff', 'seeingFwhmGeom'],
                 col_group=['filter', 'night'],
                 col_coadd=['fiveSigmaDepth', 'visitExposureTime']):

        self.col_sum = col_sum
        self.col_mean = col_mean
        self.col_median = col_median
        self.col_coadd = col_coadd
        self.col_group = col_group
        self.visitTimeCol = col_sum[1]
        self.exptimeCol = col_coadd[1]

    def _run(self, simData, cols_present=False):
        """Main run method

        Parameters
        ---------------

        simulation data

        Returns
        -----------

        numpy array with the following fields:
        fieldRA : RA of the field (median per night)
        fieldDec: Dec of the field (median per night)
        fiveSigmaDepth: coadded m5
        night: night number
        filter: filter
        numExposures: number of exposures (sum per night)
        visitTime: visit time (sum per night)
        visitExposureTime: visit exposure time (sum per night)

        """
        if 'note' in simData.dtype.names:
            simData = rf.drop_fields(simData, 'note')
        if cols_present:
            # Column already present in data; assume it is correct and does not need recalculating.
            return simData
        self.dtype = simData.dtype

        if self.visitTimeCol not in simData.dtype.names:
            simData = rf.append_fields(
                simData, self.visitTimeCol, [999.]*len(simData))

        r = []

        # print(type(simData))
        df = pd.DataFrame(np.copy(simData))

        # print(df)
        # time_ref = time.time()

        #keygroup = [self.filterCol, self.nightCol]
        keygroup = self.col_group
        """
        keysums =  [self.numExposuresCol,self.visitExposureTimeCol]
        if self.visitTimeCol in simData.dtype.names:
            keysums += [self.visitTimeCol]
        keymeans = [self.mjdCol, self.RACol, self.DecCol, self.m5Col]
        """

        exptime_single = df[self.exptimeCol].median()
        groups = df.groupby(keygroup)
        listref = df.columns
        tt = self.get_vals(listref, self.col_sum, groups, np.sum)
        vv = self.get_vals(listref, self.col_mean, groups, np.mean)
        if not vv.empty:
            tt = tt.merge(vv, left_on=['night', 'filter'],
                          right_on=['night', 'filter'])
        vv = self.get_vals(listref, self.col_median, groups, np.median)
        if not vv.empty:
            tt = tt.merge(vv, left_on=['night', 'filter'],
                          right_on=['night', 'filter'])

        tt = tt.sort_values(by=['night'])
        tt.loc[:, self.col_coadd[0]] += 1.25 * \
            np.log10(tt[self.col_coadd[1]]/exptime_single)
        return tt.to_records(index=False)

    """
        df.sort_values(by=keygroup, ascending=[True, True], inplace=True)
        # print('before',df[keygroup+keysums+keymeans])
        coadd_df = df.groupby(keygroup).agg({self.numExposuresCol: ['sum'],
                                             self.visitTimeCol: ['sum'],
                                             self.visitExposureTimeCol: ['sum'],
                                             self.mjdCol: ['mean'],
                                             self.RACol: ['mean'],
                                             self.DecCol: ['mean'],
                                             self.m5Col: ['mean'],
                                             self.seeingaCol: ['mean'],
                                             self.seeingbCol: ['mean'],
                                             'pixRA': ['mean'],
                                             'pixDec': ['mean'],
                                             'healpixID': ['mean'],
                                             'season': ['mean'],
                                             self.airmassCol: ['median'],
                                             # self.skyCol: ['median'],
                                             self.moonCol: ['median']}).reset_index()
        coadd_df.columns = [self.filterCol, self.nightCol, self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol, self.mjdCol, self.RACol,
                            self.DecCol, self.m5Col, self.seeingaCol, self.seeingbCol, 'pixRA', 'pixDec', 'healpixID', 'season', self.airmassCol, self.moonCol]

        # groupa = df.groupby(keygroup)[keysums].sum()[keymeans].mean()

        coadd_df.loc[:, self.m5Col] += 1.25 * \
            np.log10(coadd_df[self.visitTimeCol]/30.)

        # print(coadd_df.sort_values(by=[self.nightCol]))

        coadd_df.sort_values(by=[self.filterCol, self.nightCol], ascending=[
                             True, True], inplace=True)
        # print('coadd',coadd_df)
        # print('in stacker',coadd_df[[self.filterCol,self.m5Col,self.visitExposureTimeCol,self.visitTimeCol]])
        return coadd_df.to_records(index=False)
    """

    def get_vals(self, listref, cols, group, op):
        """
        Method to estimate a set of value using op

        Parameters
        --------------
        listref: list(str)
          list of columns of the group
        cols: list(str)
          list of columns to estimate values (op)
        group: pandas group
          data to process
        op: operator
         operator to apply

        Returns
        ----------
          pandas df with oped values

        """
        col_df = list(set(listref).intersection(cols))
        res = pd.DataFrame()
        if col_df:
            res = group[col_df].apply(lambda x: op(x)).reset_index()

        return res

    def fill(self, tab):
        """
        Field values estimation per night

        Parameters
        ---------------

        tab input table of field values (list given above)


        Returns
        -----------

        Field values per night : list
        all fields but m5, numexposures, visittime, visitexposuretime, filter: median value
        m5 : coadded (cf m5_coadd)
        numexposures, visittime, visitexposuretime: added per night
        band: unique band value

        """

        r = []

        for colname in self.dtype.names:
            # print('there lan',colname)
            if colname in ['note']:
                r.append(np.unique(tab[colname])[0])
                continue
            if colname not in [self.m5Col, self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol, self.filterCol]:
                if colname == 'sn_coadd':
                    r.append(1)
                else:
                    r.append(np.median(tab[colname]))
            if colname == self.m5Col:
                r.append(self.m5_coadd(np.copy(tab[self.m5Col])))
            if colname in [self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol]:
                r.append(np.sum(tab[colname]))
            if colname == self.filterCol:
                r.append(np.unique(tab[self.filterCol])[0])

        # print('done here',r)
        return r

    def m5_coadd(self, m5):
        """ Method to coadd m5 values

        .. math::
           \phi_{5\sigma} = 10^{-0.4*m_{5}}

           \sigma = \phi_{5\sigma}/5.

           \sigma_{tot} = 1./\sqrt(\sum(1./\sigma^2))

           \phi_{tot} = 5.*\sigma_{tot}

           m_{5}^{coadd} = -2.5*np.log10(\phi_{tot})

        Parameters
        --------------

        m5 : 5 sigma-depth values

        Returns
        ----------

        coadded m5 value


        """
        fluxes = 10**(-0.4*m5)
        sigmas = fluxes/5.
        sigma_tot = 1./np.sqrt(np.sum(1./sigmas**2))
        flux_tot = 5.*sigma_tot

        return -2.5*np.log10(flux_tot)

    def m5_coadd_grp(self, grp):
        """ Method to coadd m5 values

        .. math::
           \phi_{5\sigma} = 10^{-0.4*m_{5}}

           \sigma = \phi_{5\sigma}/5.

           \sigma_{tot} = 1./\sqrt(\sum(1./\sigma^2))

           \phi_{tot} = 5.*\sigma_{tot}

           m_{5}^{coadd} = -2.5*np.log10(\phi_{tot})

        Parameters
        --------------

        m5 : 5 sigma-depth values

        Returns
        ----------

        coadded m5 value


        """
        fluxes = 10**(-0.4*grp[self.m5Col])
        sigmas = fluxes/5.
        sigma_tot = 1./np.sqrt(np.sum(1./sigmas**2))
        flux_tot = 5.*sigma_tot

        grp[self.m5Col] = -2.5*np.log10(flux_tot)

        return grp

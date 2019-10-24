from lsst.sims.maf.stackers import BaseStacker
import numpy as np
import numpy.lib.recfunctions as rf
import multiprocessing
from astropy.table import Table, Column
import pandas as pd
import time

__all__ = ['CoaddStacker']


class CoaddStacker(BaseStacker):
    """
    Class to coadd observations per night
    """

    colsAdded = ['sn_coadd']

    def __init__(self, mjdCol='observationStartMJD',RaCol='fieldRA', DecCol='fieldDec', m5Col='fiveSigmaDepth', nightCol='night', filterCol='filter', numExposuresCol='numExposures', visitTimeCol='visitTime', visitExposureTimeCol='visitExposureTime', nproc=8):
        self.colsReq = [mjdCol,RaCol, DecCol, m5Col, filterCol, nightCol,
                        numExposuresCol, visitTimeCol, visitExposureTimeCol]
        self.mjdCol = mjdCol
        self.RaCol = RaCol
        self.DecCol = DecCol
        self.nightCol = nightCol
        self.filterCol = filterCol
        self.m5Col = m5Col
        self.numExposuresCol = numExposuresCol
        self.visitTimeCol = visitTimeCol
        self.visitExposureTimeCol = visitExposureTimeCol

        self.units = ['int']
        self.nproc = 1

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
            simData = rf.drop_fields(simData,'note')
        if cols_present:
            # Column already present in data; assume it is correct and does not need recalculating.
            return simData
        self.dtype = simData.dtype
        

        if self.visitTimeCol not in simData.dtype.names:
            simData = rf.append_fields(simData,self.visitTimeCol,[999.]*len(simData))

        r = []

        #print(type(simData))
        df = pd.DataFrame(np.copy(simData))

        
        #print(df)
        #time_ref = time.time()

        
        keygroup = [self.filterCol,self.nightCol]
        """
        keysums =  [self.numExposuresCol,self.visitExposureTimeCol]
        if self.visitTimeCol in simData.dtype.names:
            keysums += [self.visitTimeCol]
        keymeans = [self.mjdCol, self.RaCol, self.DecCol, self.m5Col]
        """
        
        

        df.sort_values(by=keygroup, ascending=[True, True], inplace=True)
        #print('before',df[keygroup+keysums+keymeans])
        coadd_df = df.groupby(keygroup).agg({self.numExposuresCol: ['sum'],
                                             self.visitTimeCol: ['sum'],
                                             self.visitExposureTimeCol: ['sum'],
                                             self.mjdCol: ['mean'],
                                             self.RaCol: ['mean'],
                                             self.DecCol: ['mean'],
                                             self.m5Col: ['mean'],
                                             'pixRa': ['mean'],
                                             'pixDec': ['mean'],
                                             'healpixID': ['mean'],
                                             'season': ['mean']}).reset_index()
        coadd_df.columns = [self.filterCol,self.nightCol,self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol,self.mjdCol, self.RaCol, self.DecCol, self.m5Col,'pixRa','pixDec','healpixID','season'] 

        #groupa = df.groupby(keygroup)[keysums].sum()[keymeans].mean()


        coadd_df.loc[:,self.m5Col] += 1.25*np.log10(coadd_df[self.visitExposureTimeCol]/30.)  
        #print(coadd_df.sort_values(by=[self.nightCol]))

        coadd_df.sort_values(by=[self.filterCol,self.nightCol], ascending=[True, True], inplace=True)
        #print('coadd',coadd_df)
    
        return coadd_df.to_records(index=False)

        
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
            #print('there lan',colname)
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

        #print('done here',r)
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

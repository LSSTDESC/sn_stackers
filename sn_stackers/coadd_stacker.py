from lsst.sims.maf.stackers import BaseStacker
import numpy as np
import numpy.lib.recfunctions as rf
import multiprocessing
from astropy.table import Table, Column
import pandas as pd

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
        
        if cols_present:
            # Column already present in data; assume it is correct and does not need recalculating.
            return simData
        self.dtype = simData.dtype
        
        r = []

        #print(type(simData))
        df = pd.DataFrame(np.copy(simData))

        #print(df)
        keygroup = [self.filterCol,self.nightCol]
        groups = df.groupby(keygroup)
        groups_m5 = groups.apply(lambda x: self.m5_coadd_grp(x))
        #print('coadded m5',groups_m5)

        groups = groups_m5.groupby(keygroup)

        means = groups.mean()
        med = groups.median()
        add = groups.sum()

        keys = np.array(list(groups.groups.keys()))
        #print('hello',keys)

        tab = Table()
        r = []
        for colname in self.dtype.names:
            if colname not in [self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol] and colname not in keygroup:
                if colname == 'sn_coadd':
                    tab.add_column(Column([1]*len(means),name='sn_coadd'))
                else:
                    tab.add_column(Column(med[colname],name=colname))
            #if colname == self.m5Col:
            #    r.append(self.m5_coadd(np.copy(tab[self.m5Col])))
            if colname in [self.numExposuresCol, self.visitTimeCol, self.visitExposureTimeCol]:
                tab.add_column(Column(add[colname],name=colname))
            if colname in keygroup:
                #print(colname,keys[:,keygroup.index(colname)])
                tab.add_column(Column(keys[:,keygroup.index(colname)],name=colname))

        
        return np.array(tab)
        #print('here',tab)
        #print('here',tab.dtype)

        #print(test)
        
        tab = Table(simData)

        groups = tab.group_by([self.filterCol,self.nightCol])

        indices = groups.groups.indices
        ngroups = len(indices)-1
        delta = ngroups
        if self.nproc > 1:
            delta = int(delta/(self.nproc))

        else:
            r = self.fillLoop(groups,0,ngroups)
            myarray = np.array(r, dtype=self.dtype)
            return myarray


        batch = range(0, ngroups, delta)

        if ngroups not in batch:
            batch = np.append(batch, ngroups)

        batch = batch.tolist()
        if batch[-1]-batch[-2]<= 2:
            batch.remove(batch[-2])

        #print('stacking',batch)
        result_queue = multiprocessing.Queue()

        r = []
        for j in range(len(batch)-1):
            ida = batch[j]
            idb = batch[j+1]
            p = multiprocessing.Process(name='Subprocess-'+str(j), target=self.fillLoop, args=(
                groups,batch[j],batch[j+1],j, result_queue))
            p.start()


        resultdict = {}
        for i in range(len(batch)-1):
            resultdict.update(result_queue.get())

        for p in multiprocessing.active_children():
            p.join()


        for key, vals in resultdict.items():
            r += vals

        """
        for ra, dec, band in np.unique(simData[[self.RaCol, self.DecCol, self.filterCol]]):
            idx = np.abs(simData[self.RaCol]-ra) < 1.e-5
            idx &= np.abs(simData[self.DecCol]-dec) < 1.e-5
            idx &= simData[self.filterCol] == band

            sel = simData[idx]
            for night in np.unique(sel[self.nightCol]):
                idxb = sel[self.nightCol] == night
                r.append(tuple(self.Fill(sel[idxb])))

        # print(r)
        """
        #print(r)
        myarray = np.array(r, dtype=self.dtype)
        return myarray

    def fillLoop(self,group, ida, idb,j=0, output_q=None):

        r = []

       
        for ind in range(ida,idb,1):
            
            grp = group.groups[ind]
            
            r.append(tuple(self.fill(np.copy(grp))))
            
        if output_q is not None:
            return output_q.put({j: r})
        else:
            return r


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

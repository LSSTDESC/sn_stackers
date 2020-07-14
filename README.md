# sn_stackers

Stackers used for SN pipeline.

Available stackers:
- CoaddStacker: "coadds" informations per filter and per night
  *  input: set of observations with the following fields: 
      * observationStartMJD: obs MJD
      * fieldRA: RA of obs.
      * fieldDec: Dec of obs.
      * fiveSigmaDepth: m5
      * visitExposureTime: visit exposure time
      * numExposures: number of exposures
      * visitTime: visit time
      * season: season number
      * night: night number
      * filter: band

  *  output: numpy array of observations with the following fields: 
      * observationStartMJD: median(obs MJD)
      * fieldRA: median(RA of obs.)
      * fieldDec: median(Dec of obs.)
      * fiveSigmaDepth: Coadd(m5)
      * visitExposureTime: sum(visit exposure time)
      * numExposures: sum(number of exposures)
      * visitTime: sum(visit time)
      * season: unique(season number)
      * night: unique(night number)
      * filter: unique(band)
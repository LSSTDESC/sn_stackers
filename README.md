# sn_stackers

Stackers used for SN pipeline.


``
This software was developed within the LSST DESC using LSST DESC resources, and so meets the criteria 
given in, and is bound by, the LSST DESC Publication Policy for being a "DESC product". 
We welcome requests to access code for non-DESC use; if you wish to use the code outside DESC 
please contact the developers.

```
## Release Status
|Release|Date|
|---|---|
|v1.0.0|2020/07/15|




## Feedback, License etc

If you have comments, suggestions or questions, please [write us an issue](https://github.com/LSSTDESC/sn_/issues).

This is open source software, available for re-use under the modified BSD license.

```
Copyright (c) 2020, the sn_stackers contributors on GitHub, https://github.com/LSSTDESC/sn_stackers/graphs/contributors.
All rights reserved.
```


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
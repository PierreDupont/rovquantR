.onAttach <- function(libname, pkgname) {
    
    suppressMessages({

        # dbinomLocal_normalCovs
        registerDistributions(
            list(
              dbinomLocal_normalCovs = list(
                    BUGSdist ='dbinomLocal_normalCovs(detNums       , detIndices    , size, p0       , p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                    Rdist = c('dbinomLocal_normalCovs(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0 = -999, p0Traps    , sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices = s, size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums = -999, detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)',
                              'dbinomLocal_normalCovs(detNums       , detIndices    , size, p0       , p0Traps = s, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, allowNoLocal, trapCovs, trapCovsIntercept, trapBetas)'
                    ),
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0Traps = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)', 'trapCovs = double(2)', 'trapCovsIntercept =  double(1)', 'trapBetas = double(1)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = F)
        
        # dcatHR
        registerDistributions(list(
          dcatHR = list(
            BUGSdist = "dcatHR(z, gamma, mhH, mhW)",
            Rdist = c( "dcatHR(z, gamma, mhH, mhW)"),
            types = c( "value = double(0)", "z = double(0)","gamma = double(0)", "mhH = double(0)", "mhW = double(0)"),
            discrete = TRUE,
            mixedSizes = TRUE,
            pqAvail = FALSE
          )))
        
    })
}

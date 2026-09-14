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
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0Traps = double(1)', 's = double(1)', 'trapCoords = double(2)',
                              'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)', 'trapCovs = double(2)', 'trapCovsIntercept =  double(1)', 'trapBetas = double(1)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )), verbose = F)
      
      
      # dbinomLocal_normalCovsResponse
      registerDistributions(
        list(
          dbinomLocal_normalCovsResponse = list(
            BUGSdist ='dbinomLocal_normalCovsResponse(detNums       , detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
            Rdist = c('dbinomLocal_normalCovsResponse(detNums       , detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums       , detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices = s, size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices = s, size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor     , habitatGrid, indicator, lengthYCombined = 0, trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices = s, size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums       , detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined    , trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums = -999, detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)',
                      'dbinomLocal_normalCovsResponse(detNums       , detIndices    , size, p0State, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor = 1 , habitatGrid, indicator, lengthYCombined = 0, trapCountries, trapCovs, trapBetas, trapResponse, betaResponse)'
            ),
            types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 'p0State = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)', 'trapCovs = double(2)', 'trapCovsIntercept =  double(1)', 'trapBetas = double(1)'),
            discrete = TRUE,
            mixedSizes = TRUE,
            pqAvail = FALSE
          )), verbose = F)
      
      
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
        
        
        # dbinomLocal_normalWolf
        registerDistributions(
          list(
            dbinomLocal_normalWolf = list(
              BUGSdist ='dbinomLocal_normalWolf(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor , habitatGrid, indicator, z, indCov, indBeta, trapCovs, trapCovsIntercept, trapBetas, lengthYCombined)',
              types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 's = double(1)', "p0=double(2)", 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)',"trapBetas = double(1)", "trapCovsIntercept =  double(1)","trapCovs =  double(2)"),
              discrete = TRUE,
              mixedSizes = TRUE,
              pqAvail = FALSE)),
          verbose = T)
        
        
        # dbinomLocal_normalWolverine
        registerDistributions(
          list(
            dbinomLocal_normalWolverine = list(
              BUGSdist = 'dbinomLocal_normalWolverine(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor, habitatGrid, indicator, indCov, indBeta, trapCovs, trapBetas, trapCovsIntercept, lengthYCombined)',
              types = c('value = double(1)', 'detNums = double(0)', 'detIndices = double(1)', 'size = double(1)', 'p0 = double(1)', 'sigma = double(0)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'resizeFactor = double(0)', 'habitatGrid = double(2)', 'indicator = double(0)', 'indCov = double(0)', 'indBeta = double(0)', 'trapCovs = double(2)', 'trapBetas = double(1)', 'trapCovsIntercept = double(1)', 'lengthYCombined = double(0)'),
              discrete = TRUE,
              mixedSizes = TRUE,
              pqAvail = FALSE)),
          verbose = T)

    })
}

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unordered_map>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector quantileCpp(NumericVector x, NumericVector q) 
  {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}


// [[Rcpp::export]]
int fastIntMode(NumericVector x, bool narm = false) 
{
  if (narm) x = x[!is_na(x)];
  int myMax = 1;
  int myMode = 0;
  std::unordered_map<int, int> modeMap;
  modeMap.reserve(x.size());
  
  for (std::size_t i = 0, len = x.size(); i < len; ++i) {
    auto it = modeMap.find(x[i]);
    
    if (it != modeMap.end()) {
      ++(it->second);
      if (it->second > myMax) {
        myMax = it->second;
        myMode = x[i];
      }
    } else {
      modeMap.insert({x[i], 1});
    }
  }
  
  return myMode;
}


// [[Rcpp::export]]
NumericVector extractUniquePositiveValues(NumericMatrix matrix) 
  {
  int nrow = matrix.nrow();
  int ncol = matrix.ncol();
  std::set<double> uniquePositiveValues;
  for (int i = 0; i < nrow; i++) {
    for (int j = 0; j < ncol; j++) {
      double val = matrix(i, j);
      if (val > 0) {
        uniquePositiveValues.insert(val);
      }
    }
  }
  return NumericVector(uniquePositiveValues.begin(), uniquePositiveValues.end());
}


// [[Rcpp::export]]
NumericMatrix createTransitionMatrix(NumericVector values) 
  {
  auto maxValues = std::max_element(values.begin(), values.end());
  int lengthValues = values.size();
  NumericMatrix result(*maxValues, *maxValues);
  int counter = 1;
  int thisI = 0;
  int thisJ = 0;
  for (int i = 0; i < lengthValues; i++) {
    for (int j = 0; j < lengthValues; j++) {
      thisI = values(i);
      thisJ = values(j);
      result(thisI - 1, thisJ - 1) = counter;
      counter++;
    }
  }
  return result;
}


// [[Rcpp::export]]
List GetDensity(NumericMatrix sx,             // X COORDINATES
                NumericMatrix sy,                // Y COORDINATES 
                NumericMatrix z,                 // Z STATE 
                NumericMatrix IDmx,              // MATRIX OF CELL IDS
                NumericVector aliveStates,       // VECTOR OF ALIVE STATES
                NumericMatrix regionID,          // MATRIX WITH REGION ID, ONE ROW PER REGION WITH 1 AND 0 WHETHER IT BELONGS TO THE REGION OR NOT NEED TO PROVIDE ROWNAMES TO IDENTIFY REGIONS.
                NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                bool display_progress = true,    // DISPLAY PROGRESS BAR
                bool returnPosteriorCells = true)
  {
  
  //==== INITIALIZE OBJECTS =====
  int niter = sx.nrow(), nind = sx.ncol();
  int sxmax = IDmx.ncol(), symax = IDmx.nrow();
  int ncells = regionID.ncol(), nregions = regionID.nrow();
  
  int sx1 = 1, sy1 = 1;
  NumericMatrix CellID(niter, nind);
  NumericMatrix Densi(ncells, niter);
  
  //==== INITIALIZE PROGRESS BAR ==== 
  Progress prog(niter, display_progress);
  
  //==== GET NUMBER OF AC WITHIN EACH STATE CONSIDERED AS ALIVE ====
  for (int ite = 0; ite < niter; ite++){
    // UPDATE PROGRESS BAR
    if (Progress::check_abort() )
      return -1.0;
    prog.increment(); 
    for (int i = 0; i < nind; i++){
      if (std::find(aliveStates.begin(), aliveStates.end(), z(ite,i))!=aliveStates.end()){
        if (sx(ite,i) >= 0 && sx(ite,i) < sxmax && sy(ite,i)>= 0 && sy(ite,i) < symax){
          sx1 = floor(sx(ite,i));
          sy1 = floor(sy(ite,i));
          CellID(ite,i) = IDmx(sy1,sx1);
          Densi(CellID(ite,i)-1,ite) += 1; 
        }
      }
    }
  }
  
  //==== CELL STATISTICS ====
  NumericVector outMean(ncells);
  NumericVector outMedian(ncells);
  NumericVector outsd(ncells);
  NumericVector outCV(ncells);
  NumericVector outCIL(ncells);
  NumericVector outCIH(ncells);
  
  for (int c = 0; c < ncells; c++){
    NumericVector tmp = Densi(c,_);
    outMean(c) = mean(tmp);
    outMedian(c) = median(tmp);
    outsd(c) = sd(tmp);
    outCV(c) = (100*outsd(c))/outMean(c);
    NumericVector outCILr = quantileCpp(tmp ,probs);
    outCIL(c) = outCILr(0);
    outCIH(c) = outCILr(1);
  }
  
  //================
  //==== RETURN ====
  //================
  
  //==== IF REGION STATISTICS ARE ASKED ==== 
  //==== INITIATE OBJECTS ====
  NumericMatrix subsetSum(nregions, niter);
  NumericMatrix summary ((nregions+1),5);
  // ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  CharacterVector Names =  rownames(regionID);
  // ADD A TOTAL ROW
  CharacterVector Names1 =  CharacterVector::create("Total");
  CharacterVector Names2(Names.size() + Names1.size());
  std::copy(Names.begin(), Names.end(), Names2.begin());
  std::copy(Names1.begin(), Names1.end(), Names2.begin() + Names.size());
  rownames(summary) = Names2;
  colnames(summary) = CharacterVector::create("mean", "median", "mode", "95%CILow","95%CIHigh");
  
  //==== SUMMARY AND POSTERIOR FOR EACH REGIONS ====
  for (int r = 0; r < nregions; r++){
    // SUBSET TO CELL WITHIN REGIONS 
    NumericVector T = regionID(r,_);
    mat Xmat(Densi.begin(), Densi.nrow(), Densi.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat subMat = Xmat.rows(find(tIdx == 1));
    int nsub = sum(T);
    // SUM THE ABUNDANCE OF ALL CELLS OF THE REGION FOR ALL ITERATIONS 
    for (int ite = 0; ite < niter; ite++){
      double total = 0;
      for (int j = 0; j < nsub; j++) {
        total += subMat(j, ite);
      }
      subsetSum(r,ite) = total;
    }
    // ==== FILL IN THE SUMMARY TABLE ====
    NumericVector tmp = subsetSum(r, _);
    summary(r,0) = mean(tmp);
    //summary(r,0) = round(summary(r,0), 2);
    summary(r,1) = median(tmp);
    //std::vector<int> tmp1 = as<std::vector<double> >(tmp); 
    summary(r,2) = fastIntMode(tmp);
    NumericVector outCILr= quantileCpp(tmp, probs);
    summary(r,3) = outCILr(0);
    summary(r,4) = outCILr(1);
  }
  
  //==== SUMMARY AND POSTERIOR FOR ALL REGIONS ====
  
  // SUBSET TO CELLS WITHIN ALL REGIONS 
  NumericVector AllRegionsID = colSums(regionID);
  // TURN VECTOR TO 0/1
  NumericVector AllRegionsID1 = ifelse(AllRegionsID>0, 1.0, 0.0);
  mat Xmat1(Densi.begin(), Densi.nrow(), Densi.ncol(), false);
  colvec tIdx1(AllRegionsID1.begin(), AllRegionsID1.size(), false); 
  mat subMat1 = Xmat1.rows(find(tIdx1 > 0));
  int nsub1 = sum(AllRegionsID1);
  // SUM THE ABUNDANCE FOR ALL REGIONS 
  NumericVector subsetSum1(niter) ;
  for (int ite = 0; ite < niter; ite++){
    double total = 0;
    for (int j = 0; j < nsub1; j++) {
      total += subMat1(j, ite);
    }
    subsetSum1(ite) = total;
  }
  rownames(subsetSum) = Names;
  
  
  // ==== FILL IN THE "Total" SUMMARY TABLE ====
  summary(nregions,0) = mean(subsetSum1);
  summary(nregions,1) = median(subsetSum1);
  summary(nregions,2) = fastIntMode(subsetSum1);
  NumericVector outCILr1= quantileCpp(subsetSum1, probs);
  summary(nregions,3) = outCILr1(0);
  summary(nregions,4) = outCILr1(1);
  
  // REGION SUMMARY
  if(returnPosteriorCells){
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorAllRegions") = subsetSum1,
                        Named("PosteriorCells") = Densi);
  }else{
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("PosteriorAllRegions") = subMat1,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorAllRegions") = subsetSum1);
  }
  
}


// [[Rcpp::export]]
List GetSpaceUse(NumericMatrix sx,                   // X COORDINATES 
                  NumericMatrix sy,           // Y COORDINATES 
                  NumericMatrix z,            // Z STATE 
                  NumericMatrix sigma,        // SIGMA VALUES 
                  NumericMatrix habitatxy,    // HABITAT COORDINATES
                  NumericVector aliveStates,  // ALIVE STATES 
                  NumericMatrix regionID,     // MATRIX WITH REGION ID, ONE ROW PER REGION WITH 1 AND 0 WHETHER IT BELONGS TO THE REGION OR NOT NEED TO PROVIDE ROWNAMES TO IDENTIFY REGIONS.
                  NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                  bool display_progress = true, // DISPLAY PROGRESS BAR
                  bool returnPosteriorCells = true) 
  {
  
  //INITIATE OBJECTS 
  int niter = sx.nrow(), nind = sx.ncol();
  int ncells = regionID.ncol(), nregions = regionID.nrow();
  double xMin = min(habitatxy(_,0)), xMax = max(habitatxy(_,0));
  double yMin = min(habitatxy(_,1)), yMax = max(habitatxy(_,1));
  
  NumericMatrix pTot(ncells, niter);
  
  //INITIATE PROGRESS BAR
  Progress prog(niter, display_progress);
  
  //LOOP FOR EACH MCMC ITERATION  
  for(int ite = 0; ite < niter; ite++){
    //UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog.increment(); 
    
    //LOOP FOR EACH INDIVIDUAL
    NumericVector p(ncells);
    NumericVector dist2(ncells);
    NumericVector pTemp(ncells);
    
    for (int i = 0; i < nind; i++) {
      //CHECK IF ALIVE OR NOT
      if(std::find(aliveStates.begin(), aliveStates.end(), z(ite,i))!=aliveStates.end()) {
        if(sx(ite,i) >= xMin && sx(ite,i) <= xMax && sy(ite,i) >= yMin && sy(ite,i) <= yMax) {
          double sigma1 = (2*sigma(ite,i)*sigma(ite,i));
          //LOOP FOR EACH HABITAT CELL
          for (int j = 0; j < ncells; j++) {
            //CALCULATE DISTANCE^2
            dist2(j) = (pow(sx(ite,i) - habitatxy(j, 0),2) + pow(sy(ite,i) - habitatxy(j, 1),2)) ;
          }
          //CALCULATE HALFNORMAL DENSITY
          pTemp = exp((-dist2)/sigma1);
          //NORMALIZE TO GET SPACE-USE
          p += pTemp/sum(pTemp);
        }
      }
    }
    pTot(_,ite) = p;
  }
  
  //CELL STATISTICS 
  NumericVector outMean(ncells);
  NumericVector outMedian(ncells);
  NumericVector outsd(ncells);
  NumericVector outCV(ncells);
  NumericVector outCIL(ncells);
  NumericVector outCIH(ncells);
  
  for (int c = 0; c < ncells; c++){
    NumericVector tmpMean = pTot(c,_);
    outMean(c) = mean(tmpMean);
    outMedian(c) = median(tmpMean);
    outsd(c) = sd(tmpMean);
    outCV(c) = (100*outsd(c))/outMean(c);
    NumericVector outCILr= quantileCpp(tmpMean ,probs);
    outCIL(c) = outCILr(0);
    outCIH(c) = outCILr(1);
  }
  
  //================
  //==== RETURN ====
  //================
  //==== IF REGION STATISTICS ARE ASKED ==== 
  //INITIATE OBJECTS
  NumericMatrix subsetSum(nregions, niter);
  NumericMatrix summary ((nregions+1),5);
  //ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  CharacterVector Names =  rownames(regionID);
  //ADD A TOTAL ROW
  CharacterVector Names1 =  CharacterVector::create("Total");
  CharacterVector Names2(Names.size() + Names1.size());
  std::copy(Names.begin(), Names.end(), Names2.begin());
  std::copy(Names1.begin(), Names1.end(), Names2.begin() + Names.size());
  rownames(summary) = Names2;
  colnames(summary) = CharacterVector::create("mean", "median", "mode", "95%CILow","95%CIHigh");
  
  //SUMMARY AND POSTERIOR FOR EACH REGIONS 
  for (int r = 0; r < nregions; r++){
    //SUBSET TO CELL WITHIN REGIONS 
    NumericVector T = regionID(r,_);
    mat Xmat(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat subMat = Xmat.rows(find(tIdx == 1));
    int nsub = sum(T);
    // SUM THE ABUNDANCE OF ALL CELLS OF THE REGION FOR ALL ITERATIONS 
    for (int ite = 0; ite < niter; ite++){
      double total = 0;
      for (int j = 0; j < nsub; j++) {
        total += subMat(j,ite);
      }
      subsetSum(r,ite) = total;
    }
    //FILL IN THE SUMMARY TABLE
    NumericVector tmp = subsetSum(r, _);
    summary(r,0) = mean(tmp);
    summary(r,1) = median(tmp);
    summary(r,2) = fastIntMode(tmp);
    NumericVector outCILr= quantileCpp(tmp, probs);
    summary(r,3) = outCILr(0);
    summary(r,4) = outCILr(1);
  }
  
  //SUMMARY AND POSTERIOR FOR ALL REGIONS
  //SUBSET TO CELLS WITHIN ALL REGIONS 
  NumericVector AllRegionsID= colSums(regionID);
  mat Xmat1(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
  colvec tIdx1(AllRegionsID.begin(), AllRegionsID.size(), false); 
  mat subMat1 = Xmat1.rows(find(tIdx1 > 0));
  int nsub1 = sum(AllRegionsID);
  
  //SUM THE ABUNDANCE FOR ALL REGIONS 
  NumericVector subsetSum1(niter) ;
  for (int ite = 0; ite < niter; ite++){
    double total = 0;
    for (int j = 0; j < nsub1; j++) {
      total += subMat1(j, ite);
    }
    subsetSum1(ite) = total;
  }
  
  //FILL IN THE "Total" SUMMARY TABLE
  summary(nregions,0) = mean(subsetSum1);
  summary(nregions,1) = median(subsetSum1);
  summary(nregions,2) = fastIntMode(subsetSum1);
  NumericVector outCILr1 = quantileCpp(subsetSum1, probs);
  summary(nregions,3) = outCILr1(0);
  summary(nregions,4) = outCILr1(1);
  
  //OUTPUT
  if(returnPosteriorCells){
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorCells") = pTot);
  }else{
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum);
  }
}


// [[Rcpp::export]]
List GetDetectability_normal( NumericMatrix p0,                 // detector-specific values of p0   
                               NumericVector sigma,              // SIGMA VALUES 
                               NumericMatrix habitatxy,          // HABITAT COORDINATES
                               NumericMatrix detectorxy,         // DETECTOR COORDINATES
                               NumericMatrix regionID,           // MATRIX WITH REGION ID, ONE ROW PER REGION WITH 1 AND 0 WHETHER IT BELONGS TO THE REGION OR NOT NEED TO PROVIDE ROWNAMES TO IDENTIFY REGIONS.
                               NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                               double localDist = 100,
                               bool display_progress = true,     // DISPLAY PROGRESS BAR
                               bool returnPosteriorCells = false)
  {
  // INITIALIZE OBJECTS 
  int niter = p0.nrow(); 
  int ncells = habitatxy.nrow();
  int ndets = detectorxy.nrow();
  int nregions = regionID.nrow();
  NumericMatrix pTot(ncells, niter);
  double maxD2 = localDist*localDist;
  
  // CALCULATE DISTANCE^2
  Progress prog1(ncells, display_progress);
  NumericMatrix dist2(ncells,ndets);
  for (int i = 0; i < ncells; i++) {
    if(Progress::check_abort())return -1.0;
    prog1.increment(); 
    for (int j = 0; j < ndets; j++) {
      dist2(i,j) = pow(habitatxy(i,0) - detectorxy(j,0),2) + pow(habitatxy(i,1) - detectorxy(j,1),2);
    }
  }
  
  
  // INITIATE PROGRESS BAR
  Progress prog2(niter, display_progress);
  
  // LOOP FOR EACH MCMC ITERATION  
  for(int ite = 0; ite < niter; ite++){
    // UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog2.increment(); 
    
    // PREPARE SIGMA
    double sigma1 = (2*sigma(ite)*sigma(ite));
    
    // LOOP FOR EACH HABITAT CELL
    for (int i = 0; i < ncells; i++){
      double prod = 1;
      // LOOP FOR EACH DETECTOR
      for (int j = 0; j < ndets; j++) {
        if(dist2(i,j) < maxD2){
          
          // CALCULATE 1-DETECTION PROBABILITY
          double OneMinusP = 1 - (p0(ite,j) * exp(-dist2(i,j)/sigma1));
          // TAKE THE PRODUCT
          prod *= OneMinusP;
        }
      }
      pTot(i,ite) = 1-prod;
    }
  }
  
  // CELL STATISTICS 
  NumericVector outMean(ncells);
  NumericVector outMedian(ncells);
  NumericVector outsd(ncells);
  NumericVector outCV(ncells);
  NumericVector outCIL(ncells);
  NumericVector outCIH(ncells);
  
  for (int i = 0; i < ncells; i++){
    NumericVector tmpMean = pTot(i,_);
    outMean(i) = mean(tmpMean);
    outMedian(i) = median(tmpMean);
    outsd(i) = sd(tmpMean);
    outCV(i) = (100*outsd(i))/outMean(i);
    NumericVector outCILr= quantileCpp(tmpMean ,probs);
    outCIL(i) = outCILr(0);
    outCIH(i) = outCILr(1);
  }
  
  //================
  //==== RETURN ====
  //================
  // INITIALIZE OBJECTS
  NumericMatrix subsetSum(nregions, niter);
  NumericMatrix summary(nregions+1,4);
  // ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  CharacterVector Names = rownames(regionID);
  // ADD A TOTAL ROW
  CharacterVector Names1 = CharacterVector::create("Total");
  CharacterVector Names2(Names.size() + Names1.size());
  std::copy(Names.begin(), Names.end(), Names2.begin());
  std::copy(Names1.begin(), Names1.end(), Names2.begin() + Names.size());
  rownames(summary) = Names2;
  colnames(summary) = CharacterVector::create("mean", "median", "95%CILow","95%CIHigh");
  
  // SUMMARY AND POSTERIOR FOR EACH REGIONS 
  for (int r = 0; r < nregions; r++){
    // SUBSET TO CELL WITHIN REGIONS 
    NumericVector T = regionID(r,_);
    colvec tIdx(T.begin(), T.size(), false); 
    mat Xmat(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
    mat subMat = Xmat.rows(find(tIdx == 1));
    int nsub = sum(T);
    // AVERAGE THE DETECTION PROBABILITY OVER ALL CELLS OF THE REGION
    for (int ite = 0; ite < niter; ite++){
      double total = 0;
      for (int j = 0; j < nsub; j++) {
        double meanP = subMat(j,ite)/nsub;
        total += meanP;
      }
      subsetSum(r,ite) = total;
    }
    // FILL IN THE SUMMARY TABLE
    NumericVector tmp = subsetSum(r,_);
    summary(r,0) = mean(tmp);
    summary(r,1) = median(tmp);
    NumericVector outCILr = quantileCpp(tmp, probs);
    summary(r,2) = outCILr(0);
    summary(r,3) = outCILr(1);
  }
  
  // SUMMARY AND POSTERIOR FOR ALL REGIONS
  // SUBSET TO CELLS WITHIN ALL REGIONS 
  NumericVector AllRegionsID = colSums(regionID);
  mat Xmat1(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
  colvec tIdx1(AllRegionsID.begin(), AllRegionsID.size(), false); 
  mat subMat1 = Xmat1.rows(find(tIdx1 > 0));
  int nsub1 = sum(AllRegionsID);
  
  // AVERAGE THE DETECTION PROBABILITY OVER ALL CELLS
  NumericVector subsetSum1(niter) ;
  for (int ite = 0; ite < niter; ite++){
    double total = 0;
    for (int j = 0; j < nsub1; j++) {
      total += subMat1(j, ite)/nsub1;
    }
    subsetSum1(ite) = total;
  }
  
  // FILL IN THE "Total" SUMMARY TABLE
  summary(nregions,0) = mean(subsetSum1);
  summary(nregions,1) = median(subsetSum1);
  NumericVector outCILr1 = quantileCpp(subsetSum1, probs);
  summary(nregions,2) = outCILr1(0);
  summary(nregions,3) = outCILr1(1);
  
  // OUTPUT
  if(returnPosteriorCells){
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorCells") = pTot);
  } else {
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum);
  }
}


// [[Rcpp::export]]
List GetDetectability_mean( NumericVector p0,             // detector-specific values of p0   
                             double sigma,                 // SIGMA VALUES 
                             NumericMatrix habitatxy,      // HABITAT COORDINATES
                             NumericMatrix detectorxy,     // DETECTOR COORDINATES
                             bool display_progress = true)
  {
  // INITIALIZE OBJECTS 
  int ncells = habitatxy.nrow();
  int ndets = detectorxy.nrow();
  NumericVector pTot(ncells);
  double sigma1 = 2*sigma*sigma;
  
  // INITIATE PROGRESS BAR
  Progress prog2(ncells, display_progress);
  
  // LOOP FOR EACH HABITAT CELL
  for (int i = 0; i < ncells; i++){
    // UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog2.increment(); 
    
    double prod = 1;
    
    // LOOP FOR EACH DETECTOR
    for (int j = 0; j < ndets; j++) {
      // CALCULATE SQUARED DISTANCE
      double d2 = pow(habitatxy(i,0) - detectorxy(j,0),2) + pow(habitatxy(i,1) - detectorxy(j,1),2);
      
      // CALCULATE 1-DETECTION PROBABILITY
      double OneMinusP = 1 - (p0(j) * exp(-d2/sigma1));
      
      // TAKE THE PRODUCT
      prod *= OneMinusP;
    }
    pTot(i) = 1-prod;
  }
  
  // OUTPUT
  return List::create(Named("MeanCell") = pTot);
}


// [[Rcpp::export]]
List GetTransitions( NumericMatrix sx1,               // X COORDINATES
                     NumericMatrix sy1,               // Y COORDINATES 
                     NumericMatrix z1,                // Z STATE 
                     NumericMatrix sx2,        // X COORDINATES t+1
                     NumericMatrix sy2,        // Y COORDINATES t+1
                     NumericMatrix z2,         // Z STATE t+1
                     NumericVector stateFrom,  // STATE FROM 
                     NumericVector stateTo,    // STATE TO 
                     NumericMatrix rgmx,       // MATRIX WITH REGION ID
                     NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                     bool display_progress = true)
  {
  //== EXTRACT USEFUL VALUES
  int niter = sx1.nrow();
  int nind = sx1.ncol();
  NumericVector regions = extractUniquePositiveValues(rgmx);
  int nrgn = regions.size();
  int ntrn = nrgn*nrgn;
  NumericMatrix trnsmx = createTransitionMatrix(regions);
  
  int sxmax = rgmx.ncol();
  int symax = rgmx.nrow();
  
  //== INITIALIZE OBJECTS
  int sx1ID = 1, sy1ID = 1;
  int sx2ID = 1, sy2ID = 1;
  int r1 = 1, r2 = 1;
  //int tID = 1;
  
  NumericMatrix numInd(ntrn, niter);
  NumericMatrix region1ID(niter,nind), region2ID(niter,nind);
  NumericMatrix trnsID(niter,nind);
  
  //== INITIALIZE PROGRESS BAR
  Progress prog(niter, display_progress);
  
  //==== GET NUMBER OF INDIVIDUALS WHO TRANSITION FROM "stateFrom" TO "stateTo"
  //==== FOR ALL COMBINATIONS OF REGIONS OF ORIGIN AND ARRIVAL
  for (int ite = 0; ite < niter; ite++){
    //== UPDATE PROGRESS BAR
    if (Progress::check_abort())return -1.0;
    prog.increment(); 
    //== LOOP THROUGH INDIVIDUALS
    for (int i = 0; i < nind; i++){
      //== IDENTIFY INDIVIDUALS DOING THE FOCAL TRANSITION
      if (std::find(stateFrom.begin(), stateFrom.end(), z1(ite,i))!=stateFrom.end() &&
          std::find(stateTo.begin(), stateTo.end(), z2(ite,i))!=stateTo.end()){
        //== IDENTIFY INDIVIDUALS STAYING IN THE FOCAL REGION(S)
        if (sx1(ite,i)>=0 && sx1(ite,i)<sxmax && sy1(ite,i)>=0 && sy1(ite,i)<symax &&
            sx2(ite,i)>=0 && sx2(ite,i)<sxmax && sy2(ite,i)>=0 && sy2(ite,i)<symax){
          //== EXTRACT REGION OF DEPARTURE
          sx1ID = floor(sx1(ite,i));
          sy1ID = floor(sy1(ite,i));
          region1ID(ite,i) = rgmx(sy1ID,sx1ID);
          //== EXTRACT REGION OF ARRIVAL
          sx2ID = floor(sx2(ite,i));
          sy2ID = floor(sy2(ite,i));
          region2ID(ite,i) = rgmx(sy2ID,sx2ID);
          if (region1ID(ite,i) > 0 && region2ID(ite,i) > 0){
            //== IDENTIFY CORRESPONDING TRANSITION
            r1 = region1ID(ite,i);
            r2 = region2ID(ite,i);
            trnsID(ite,i) = trnsmx(r1-1,r2-1);
            //== ADD ONE INDIVIDUAL TO THE CORRESPONDING TRANSITION
            numInd(trnsID(ite,i)-1,ite) += 1;
          }
        }
      }
    }
  }
  
  
  
  //================
  //==== OUTPUT ====
  //================
  
  //==== INITIATE OBJECTS ====
  NumericMatrix summary(ntrn,8);
  
  //== ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  colnames(summary) = CharacterVector::create("from","to","mean","median","sd","CV","95%CILow","95%CIHigh");
  
  //==== IDENTIFY TRANSITIONS BY THEIR DEPARTURE AND ARRIVAL REGIONS
  int counter = 0;
  for (int i = 0; i < nrgn; i++) {
    for (int j = 0; j < nrgn; j++) {
      summary(counter, 0) = regions(i);
      summary(counter, 1) = regions(j);
      counter++;
    }
  }
  
  //==== STATISTICS PER TRANSITION ====
  for (int c = 0; c < ntrn; c++){
    NumericVector tmp = numInd(c,_);
    summary(c,2) = mean(tmp);
    summary(c,3) = median(tmp);
    summary(c,4) = sd(tmp);
    summary(c,5) = (100*sd(tmp))/mean(tmp);
    NumericVector outCILr = quantileCpp(tmp, probs);
    summary(c,6) = outCILr(0);
    summary(c,7) = outCILr(1);
  }
  
  
  
  // REGION SUMMARY
  return List::create( Named("summary") = summary,
                       Named("PosteriorSamples") = numInd,
                       Named("trnsmx") = trnsmx,
                       Named("trnsID") = trnsID,
                       Named("r1") = region1ID,
                       Named("r2") = region2ID);
  
}

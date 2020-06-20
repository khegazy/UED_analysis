#include "mergeClass.h"


mergeClass::mergeClass(std::string runName_inp) : parameterClass(runName_inp) {

  runName = runName_inp;
  initializeVariables();
}


mergeClass::~mergeClass() {
  delete plt;
}


void mergeClass::initializeVariables() {

  plt = new PLOTclass("mergeCanv");
  curScan  = -1;
  NlegBins = NradLegBins*Nlegendres;
  runInd   = "";
  curDate  = "";
  curRun   = "";
  _compareReference = false;
  smearedTime       = false;
  didSMSnormalize   = false;
  didPairCorrSTD    = false;
  didSubtractT0     = false;
  SEMisBootstrap    = false;

  /////  Simulation results  /////
  /*
  atmLegDiff.resize(NradLegBins, 0.0);
  molLegDiff.resize(NradLegBins, 0.0);
  simFileNameSuffix = 
      "_Qmax-" + to_string(maxQleg)
      + "_Ieb-" + to_string(Iebeam) 
      + "_scrnD-" + to_string(screenDist) 
      + "_elE-" + to_string(elEnergy) 
      + "_Bins[" + to_string(NradLegBins) + "].dat"; 
  save::importDat<double>(atmLegDiff, simOutputDir + "/" 
              + molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);
  save::importDat<double>(molLegDiff, simOutputDir + "/" 
              + molName + "_molDiffractionPatternLineOut"
              + simFileNameSuffix);
              */

  atmAzmDiff.resize(NradAzmBins, 0.0);
  molAzmDiff.resize(NradAzmBins, 0.0);
  simFileNameSuffix = 
      "_Qmax-" + to_string(maxQazm)
      + "_Ieb-" + to_string(Iebeam) 
      + "_scrnD-" + to_string(screenDist) 
      + "_elE-" + to_string(elEnergy) 
      + "_Bins[" + to_string(NradAzmBins) + "].dat";
  save::importDat<double>(atmAzmDiff, simOutputDir + "/" 
              + molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);
  save::importDat<double>(molAzmDiff, simOutputDir + "/" 
              + molName + "_molDiffractionPatternLineOut"
              + simFileNameSuffix);


  // Calculate normalization coefficients for modified molecular scattering
  Qazm.resize(NradAzmBins);
  for (int iq=0; iq<NradAzmBins; iq++) {
    Qazm[iq] = maxQazm*(iq + 0.5)/NradAzmBins;
  }
  Qleg.resize(NradLegBins);
  for (int iq=0; iq<NradLegBins; iq++) {
    Qleg[iq] = maxQleg*(iq + 0.5)/NradLegBins;
  }
  calculateSMS();

  // Reference
  azmReference.resize(NradAzmBins, 0);
  legReference.resize(Nlegendres);
  for (auto& v : legReference) {
    v.resize(NradLegBins, 0);
  }


  // Output files and location
  fileName = "mergedScans.txt";

  legendres.resize(Nlegendres);
  legendresMs.resize(Nlegendres);
}


void mergeClass::compareReference(std::string refAddr) {

  _compareReference = true;

  compReference.resize(NradAzmBins);
  save::importDat<double>(compReference, refAddr);
}

void mergeClass::calculateSMS() {

  if (sMsLegNorm.size()) {
    std::cerr << "WARNING: Replacing sMsLegNorm values!!!\n";
  }

  sMsLegNorm.resize(NradLegBins, 0.0);
  /*
  for (int ir=0; ir<NradLegBins; ir++) {
    sMsLegNorm[ir] = Qleg[ir]/(atmLegDiff[ir]);
  }
  */

  //cerr<<"WARNING: CHANGED SMS NORM"<<endl;
  sMsAzmNorm.resize(NradAzmBins, 0.0);
  for (int ir=0; ir<NradAzmBins; ir++) {
    sMsAzmNorm[ir] = Qazm[ir]/(atmAzmDiff[ir]);
    //sMsAzmNorm[ir] /= Qazm[ir]*Qazm[ir];
  }
}


void mergeClass::compareSimulations(std::vector<std::string> radicals) {

  if (verbose) std::cout << "INFO: Entering comareFinalStates!!!\n";

  std::map< string, std::vector<double> > refAutCor, refLineOut;
  for (auto name : radicals) {
    if (molName.compare(name) != 0) {

      if (verbose) std::cout << "\tImporting " + name + " simulations\n";


      //////////////////////////////////
      /////  Import and normalize  /////
      //////////////////////////////////

      // Import data
      refLineOut[name].resize(NradLegBins, 0);
      save::importDat(refLineOut[name], simOutputDir + "/" 
                        + name + "_molDiffractionPatternLineOut"
                        + simFileNameSuffix); 
      
      if (pltVerbose) 
        delete plt->print1d(refLineOut[name], "./plots/sim_" + name + "_diffraction");

      // Normalize by atomic scattering
      /*
      for (int ir=0; ir<NradLegBins; ir++) {
        refLineOut[name][ir] *= sMsLegNorm[ir];
      } 
      */
      // Save sMs
      save::saveDat<double>(refLineOut[name], "./results/data/sim_" + name + "Diffraction["
            + to_string(NradLegBins) + "].dat");
      if (pltVerbose) 
        delete plt->print1d(refLineOut[name], "./sim_" + name + "_diffractionsMs");


      //////////////////////////////
      /////  Pair correlation  /////
      //////////////////////////////

      /*
      if (verbose) std::cout << "\tCalculating pair correlation\n";

      // Filling input and padding 
      std::vector<double> inpDiff(NradLegBins*2+1, 0);
      std::vector< std::vector<double> > autCor1d;
      int centI = (int)(inpDiff.size()/2);
      for (int ir=0; ir<NradLegBins; ir++) {
        if (ir < indHR) {
          inpDiff[centI+1+ir] = refLineOut[name][indHR]
                *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
          //inpDiff[centI+1+ir] = smearedImg[0][it][ir]
          //      *pow(sin((PI/2)*((ir+1)/((double)(indHR+1)))), 2);
          inpDiff[centI-1-ir] = inpDiff[centI+1+ir];
        }
        else {
          inpDiff[centI+1+ir] = refLineOut[name][ir];
          inpDiff[centI-1-ir] = refLineOut[name][ir];
        }
      }
      if (pltVerbose) 
        delete plt.print1d(inpDiff, name + "Sim_inpDiff");

      autCor1d = tools::fft1dRtoC(inpDiff, rMaxRat,
              NautCpadding, padDecayRat, fftB, qSpace, rSpace, false);
      refAutCor[name] = autCor1d[0];
      delete plt.print1d(refAutCor[name], name + "_pairCorrSim");
      save::saveDat<double>(refAutCor[name], "./plots/data/sim_" + name + "PairCorr["
              + to_string(refAutCor[name].size()) + "].dat");
      */
    }
  }
  if (verbose) std::cout << "\tExiting comareFinalStates\n";
}


void mergeClass::removeBadRegions(
    std::vector<double>* azmAvg,
    int64_t stagePos) {

  if (badRegions.find(stagePos) != badRegions.end()) {
    int ind1, ind2;
    for (auto const & rItr : badRegions[stagePos]) {
      ind1 = NradAzmBins*rItr.first/maxQazm;
      ind2 = NradAzmBins*rItr.second/maxQazm;
      for (int i=ind1; i<ind2; i++) {
        (*azmAvg)[i] = NANVAL;
      }
    }
  }

  return;
}


void mergeClass::addLabTimeParameter(
      long int timeStamp,
      std::string name, 
      double value) {

  labTimeParams[timeStamp][name] = value;
}


void mergeClass::addReference(int scan, int64_t stagePos, int timeStamp,
                          std::vector<double>* azmAvg, 
                          std::vector<double>* legCoeffs,
                          double imgNorm) {

  scanReferences[scan][stagePos].scale    = 1;
  scanReferences[scan][stagePos].imgNorm  = imgNorm;
  scanReferences[scan][stagePos].azmRef.resize(NradAzmBins, 0);
  for (int i=0; i<NradAzmBins; i++) {
    if ((*azmAvg)[i] != NANVAL) {
      if (mergeNormalizeImgs) {
        scanReferences[scan][stagePos].azmRef[i] = (*azmAvg)[i]/imgNorm;
      }
      else {
        scanReferences[scan][stagePos].azmRef[i] = (*azmAvg)[i];
      }
    }
    else {
      scanReferences[scan][stagePos].azmRef[i] = NANVAL;
    }
  }

  scanReferences[scan][stagePos].legRef.resize(NlegBins, 0);
  for (int i=0; i<NlegBins; i++) {
    if ((*legCoeffs)[i] != NANVAL) {
      if (mergeNormalizeImgs) {
        scanReferences[scan][stagePos].legRef[i] = (*legCoeffs)[i]/imgNorm;
      }
      else {
        scanReferences[scan][stagePos].legRef[i] = (*legCoeffs)[i];
      }
    }
    else {
      scanReferences[scan][stagePos].legRef[i] = NANVAL;
    }
  }

  labTimeMap[timeStamp].first   = scan;
  labTimeMap[timeStamp].second  = stagePos;

}


void mergeClass::addEntry(int scan, int64_t stagePos, int timeStamp,
                          std::vector<double>* azmAvg, 
                          std::vector<double>* legCoeffs,
                          double imgNorm) {
  ///////////////////////////////////////////////////////
  /////  Insert each legendre projection into a     /////
  /////  map< scan, 2d vector>  of shape            /////
  /////  [scan, stagePos, legCoeff] that is later   /////
  /////  used to remove pixel outliers by standard  /////
  /////  deviation cuts over the scan axis.         /////
  ///////////////////////////////////////////////////////

  ///  Get index of stage position  ///
  if (scanLgndrs.find(scan) == scanLgndrs.end()) {
    std::vector<double> emptyC(stagePosInds.size(), 0);
    scanScale[scan]     = emptyC;
    scanImgNorms[scan]  = emptyC; 

    std::vector< std::vector<double> > empty(stagePosInds.size());
    for (uint ipos=0; ipos<stagePosInds.size(); ipos++) {
      empty[ipos].resize(NlegBins, 0);
    }
    scanLgndrs[scan] = empty;
    
    for (uint ipos=0; ipos<stagePosInds.size(); ipos++) {
      empty[ipos].resize(NradAzmBins, 0);
    }
    scanAzmAvg[scan] = empty;
  }

  ///  Get index of stage position  ///
  auto pos = stagePosInds.find(stagePos);
  int pInd = pos->second;

  /// Add new stage position to index lookup table and all runs  ///
  if (pos == stagePosInds.end()) {
    int ind = 0;
    stagePosInds[stagePos] = -1;
    for (auto& itr : stagePosInds) {
      if (itr.second == -1) {
        pInd = ind;
        auto sCitr = scanScale.begin();
        auto sImgN = scanImgNorms.begin();
        auto sLitr = scanLgndrs.begin();
        auto sAzml = scanAzmAvg.begin();
        const std::vector<double> emptyLegVec((*legCoeffs).size(), 0.0);
        const std::vector<double> emptyAzmVec((*azmAvg).size(), 0.0);
        while (sCitr != scanScale.end()) {
          sCitr->second.insert(sCitr->second.begin()+ind, 0);
          sImgN->second.insert(sImgN->second.begin()+ind, -1);
          sLitr->second.insert(sLitr->second.begin()+ind, emptyLegVec);
          sAzml->second.insert(sAzml->second.begin()+ind, emptyAzmVec);
          sCitr++; sLitr++; sAzml++;
        }
      }
      itr.second = ind;
      ind++;
    }
  }

  ///  Add coefficients and counts to maps  ///
  scanScale[scan][pInd] = 1;
  scanImgNorms[scan][pInd] = imgNorm;
  for (uint i=0; i<(*legCoeffs).size(); i++) {
    if (mergeNormalizeImgs) {
      scanLgndrs[scan][pInd][i] = (*legCoeffs)[i]/imgNorm;
    }
    else {
      scanLgndrs[scan][pInd][i] = (*legCoeffs)[i];
    }
  }
  for (uint i=0; i<(*azmAvg).size(); i++) {
    if ((*azmAvg)[i] != NANVAL) {
      if (mergeNormalizeImgs) {
        scanAzmAvg[scan][pInd][i] = (*azmAvg)[i]/imgNorm;
      }
      else {
        scanAzmAvg[scan][pInd][i] = (*azmAvg)[i];
      }
    }
    else {
      scanAzmAvg[scan][pInd][i] = NANVAL;
    }
  }

  labTimeMap[timeStamp].first   = scan;
  labTimeMap[timeStamp].second  = stagePos;

  
/*
  for (int i=0; i<scanLgndrs[scan][pInd].size(); i++) {
    if ((std::isnan(scanLgndrs[scan][pInd][i]) || fabs(scanLgndrs[scan][pInd][i]) > 1e100 || fabs(scanLgndrs[scan][pInd][i]) < 1e-100) && (scanLgndrs[scan][pInd][i] != 0)) {
      cout<<"FOUNDNAN lgndr: "<<scan<<"  "<<pInd<<"  "<<i<<"  "<<scanLgndrs[scan][pInd][i]<<endl;
    }
  }
  for (int i=0; i<scanAzmAvg[scan][pInd].size(); i++) {
    if ((std::isnan(scanAzmAvg[scan][pInd][i]) || fabs(scanAzmAvg[scan][pInd][i]) > 1e100 || fabs(scanAzmAvg[scan][pInd][i]) < 1e-100) && (scanAzmAvg[scan][pInd][i]!=0)) {
      cout<<"FOUNDNAN azml: "<<scan<<"  "<<pInd<<"  "<<i<<"  "<<scanAzmAvg[scan][pInd][i]<<endl;
    }
  }
  */
  
}


void mergeClass::saveInitialData() {

  int scan, stagePos, pInd;

  // Saving references
  for (auto const & sItr : scanReferences) {
    scan = sItr.first;
    for (auto const & pItr : sItr.second) {
      stagePos = pItr.first;
      scanReferences_init[scan][stagePos].scale   = 
          scanReferences[scan][stagePos].scale;
      scanReferences_init[scan][stagePos].imgNorm = 
          scanReferences[scan][stagePos].imgNorm;
      scanReferences_init[scan][stagePos].azmRef.assign(
          scanReferences[scan][stagePos].azmRef.begin(),
          scanReferences[scan][stagePos].azmRef.end());
      scanReferences_init[scan][stagePos].legRef.assign(
          scanReferences[scan][stagePos].legRef.begin(),
          scanReferences[scan][stagePos].legRef.end());
    }
  }


  // Saving time dependent data
  for (auto const & sItr : scanScale) {
    scan = sItr.first;
    scanScale_init[scan].resize(stagePosInds.size());
    scanLgndrs_init[scan].resize(stagePosInds.size());
    scanAzmAvg_init[scan].resize(stagePosInds.size());
    for (auto const & pItr : stagePosInds) {
      pInd = pItr.second;

      scanScale_init[scan][pInd] = scanScale[scan][pInd];
      scanLgndrs_init[scan][pInd].assign(
          scanLgndrs[scan][pInd].begin(), 
          scanLgndrs[scan][pInd].end());
      scanAzmAvg_init[scan][pInd].assign(
          scanAzmAvg[scan][pInd].begin(),
          scanAzmAvg[scan][pInd].end());
    }
  }

  return;
}


void mergeClass::reloadInitialData() {

  int scan, stagePos, pInd;

  // Reloading references
  scanReferences.clear();
  for (auto const & sItr : scanReferences_init) {
    scan = sItr.first;
    for (auto const & pItr : sItr.second) {
      stagePos = pItr.first;
      scanReferences[scan][stagePos].scale   = 
          scanReferences_init[scan][stagePos].scale;
      scanReferences[scan][stagePos].imgNorm = 
          scanReferences_init[scan][stagePos].imgNorm;
      scanReferences[scan][stagePos].azmRef.assign(
          scanReferences_init[scan][stagePos].azmRef.begin(),
          scanReferences_init[scan][stagePos].azmRef.end());
      scanReferences[scan][stagePos].legRef.assign(
          scanReferences_init[scan][stagePos].legRef.begin(),
          scanReferences_init[scan][stagePos].legRef.end());
    }
  }


  // Reloading time dependent data
  scanScale.clear();
  scanLgndrs.clear();
  scanAzmAvg.clear();
  for (auto const & sItr : scanScale_init) {
    scan = sItr.first;
    scanScale[scan].resize(stagePosInds.size());
    scanLgndrs[scan].resize(stagePosInds.size());
    scanAzmAvg[scan].resize(stagePosInds.size());
    for (auto const & pItr : stagePosInds) {
      pInd = pItr.second;

      scanScale[scan][pInd] = scanScale_init[scan][pInd];
      scanLgndrs[scan][pInd].assign(
          scanLgndrs_init[scan][pInd].begin(), 
          scanLgndrs_init[scan][pInd].end());
      scanAzmAvg[scan][pInd].assign(
          scanAzmAvg_init[scan][pInd].begin(),
          scanAzmAvg_init[scan][pInd].end());
    }
  }

  return;
}


/*
void mergeClass::removeLabTimeOutliers() {

  float norm;
  auto tItr  = labTimeParams.begin();
  auto tmItr = tItr;
  int tInd = 0;
  int ltScan = 1;
  int NsmearSteps = 3*labParamSmear/pvSampleTimes + 1;
  smoothImgNorm.clear(); smoothImgNormSTD.clear();
  smoothImgNorm.resize(labTimeParams.size(), 0);
  smoothImgNormSTD.resize(labTimeParams.size(), 0);
  while (tItr != labTimeParams.end()) {
    ltScan = tItr->second.scan;
    norm = 0;

    tmItr = tItr;
    for (int i=0; i<=NsmearSteps; i++) {
      smoothImgNorm[tInd] += tmItr->second.imgNorm
          *exp(-1*std::pow(pvSampleTimes*i, 2)
              /(2*pow(labParamSmear, 2)));
      norm += exp(-1*std::pow(pvSampleTimes*i, 2)
                 /(2*pow(labParamSmear, 2)));

      if (tmItr == labTimeParams.begin()) break;
      tmItr--;
    
      // Only look at scan immediatly prior
      if (ltScan - tmItr->second.scan > 0) {
        if (ltScan - tmItr->second.scan == 1) {
          ltScan = tmItr->second.scan;
        }
        else {
          break;
        }
      }
    }
    smoothImgNorm[tInd] /= norm;

    tmItr = tItr;
    for (int i=0; i<=NsmearSteps; i++) {
      smoothImgNormSTD[tInd] += 
          pow(tItr->second.imgNorm -tmItr->second.imgNorm, 2)
          *exp(-1*std::pow(pvSampleTimes*i, 2)
              /(2*pow(labParamSmear, 2)));

      if (tmItr == labTimeParams.begin()) break;
      tmItr--;
    
      // Only look at scan immediatly prior
      if (ltScan - tmItr->second.scan > 0) {
        if (ltScan - tmItr->second.scan == 1) {
          ltScan = tmItr->second.scan;
        }
        else {
          break;
        }
      }
    }
    smoothImgNormSTD[tInd] = sqrt(smoothImgNormSTD[tInd]/norm);


    tInd++;
    tItr++;
  }

  tInd = 0;
  for (tItr=labTimeParams.begin(); tItr!= labTimeParams.end(); tInd++) {
    if (fabs(tItr->second.imgNorm - smoothImgNorm[tInd]) 
        > labSTDcut*smoothImgNormSTD[tInd]) {
      tItr->second.imgNorm = 0;
    }
    tInd++;
  }
}
*/
 

void mergeClass::getImgNormMeanSTD() {

  /////  Reset all vectors  /////
  // Time dependent measurements
  if (scanImgNormAzmMeans.size() != stagePosInds.size()) {
    scanImgNormAzmMeans.resize(stagePosInds.size());
    scanImgNormAzmSTDs.resize(stagePosInds.size());
  }
  std::fill(scanImgNormAzmMeans.begin(), scanImgNormAzmMeans.end(), 0);
  std::fill(scanImgNormAzmSTDs.begin(), scanImgNormAzmSTDs.end(), 0);
  scanImgNormAzmRefMean = 0;
  scanImgNormAzmRefSTD = 0;


  //////////////////////////////
  /////  Reference Images  /////
  //////////////////////////////

  /////  Calculate mean  /////
  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal mean.\n";
  double siarCount = 0;
  for (auto const & sRitr : scanReferences) {
    for (auto const & pItr : sRitr.second) {
      if (pItr.second.scale) {
        scanImgNormAzmRefMean += pItr.second.imgNorm;
        siarCount++;
      }
    }
  }
  scanImgNormAzmRefMean /= siarCount;

  /////  Calculate Standard Deviation  /////
  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal std.\n";
  siarCount = 0;
  for (auto const & sRitr : scanReferences) {
    for (auto const & pItr : sRitr.second) {
      if (pItr.second.scale) {
        scanImgNormAzmRefSTD += std::pow(pItr.second.imgNorm 
            - scanImgNormAzmRefMean, 2);
        siarCount++;
      }
    }
  }
  scanImgNormAzmRefSTD = std::sqrt(scanImgNormAzmRefSTD/siarCount);


  ///////////////////////////////////
  /////  Time Dependent Images  /////
  ///////////////////////////////////
 
  for (auto const & pItr : stagePosInds) {
    /////  Mean calculation  /////
    
    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal mean.\n";
    double siaCount = 0;
    for (auto const & sImgN : scanImgNorms) {
      if (scanScale[sImgN.first][pItr.second] > 0) {
          scanImgNormAzmMeans[pItr.second] += sImgN.second[pItr.second];
          siaCount++;
      }
    }

    scanImgNormAzmMeans[pItr.second] /= siaCount;


    /////  Standard deviation calculation  /////
    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal std.\n";
    siaCount = 0;
    for (auto const & sImgN : scanImgNorms) {
      if (scanScale[sImgN.first][pItr.second] > 0) {
        scanImgNormAzmSTDs[pItr.second]
            += std::pow(sImgN.second[pItr.second] 
                - scanImgNormAzmMeans[pItr.second], 2);
        siaCount++;
      }
    }

    scanImgNormAzmSTDs[pItr.second] = std::sqrt(scanImgNormAzmSTDs[pItr.second]/siaCount);
  }

  return;
}


   

void mergeClass::getImageMeanSTD() {

  /////  Reset all vectors  /////
  // Time dependent measurements
  if (imgAzmMeans.size() != scanAzmAvg.size()) {
    imgAzmMeans.resize(scanAzmAvg.size());
    imgAzmRefMeans.resize(scanReferences.size());
    auto scanItr = scanAzmAvg.begin();
    for (uint i=0; i<scanAzmAvg.size(); i++) {
      imgAzmMeans[i].resize(stagePosInds.size(), 0);
      scanInds[scanItr->first] = i;
      scanItr++;
    }
  }
  if (scanImgAzmMeans.size() != stagePosInds.size()) {
    scanImgAzmMeans.resize(stagePosInds.size());
    scanImgAzmSTDs.resize(stagePosInds.size());
  }
  for (uint i=0; i<scanAzmAvg.size(); i++) {
    std::fill(imgAzmMeans[i].begin(), imgAzmMeans[i].end(), 0);
  }
  for (uint i=0; i<scanReferences.size(); i++) {
    imgAzmRefMeans[i].clear();
  }
  std::fill(scanImgAzmMeans.begin(), scanImgAzmMeans.end(), 0);
  std::fill(scanImgAzmSTDs.begin(), scanImgAzmSTDs.end(), 0);
  scanImgAzmRefMean = 0;
  scanImgAzmRefSTD = 0;

  //std::vector<double>legNorms(NlegBins, 0);
  //std::vector<double>azmNorms(NradAzmBins, 0);

  double norm = 0;

  //////////////////////////////
  /////  Reference Images  /////
  //////////////////////////////

  /////  Calculate mean  /////
  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal mean.\n";
  int refInd = 0;
  double siarCount = 0;
  for (auto const & sRitr : scanReferences) {
    int rPosInd = 0;
    for (auto const &pItr : sRitr.second) {
      norm = 0;
      imgAzmRefMeans[refInd].push_back(0);
      for (int ir=0; ir<NradAzmBins; ir++) {
        if (pItr.second.azmRef[ir] != NANVAL) {
          imgAzmRefMeans[refInd][rPosInd] += pItr.second.azmRef[ir];
          norm += 1;
        }
      }
      if (norm) {
        imgAzmRefMeans[refInd][rPosInd] /= norm;
        scanImgAzmRefMean += imgAzmRefMeans[refInd][rPosInd];
        siarCount++;
      }
      else {
        imgAzmRefMeans[refInd][rPosInd] = 0;
      }
      rPosInd++;
    }
    refInd++;
  }
  scanImgAzmRefMean /= siarCount;

  /////  Calculate Standard Deviation  /////
  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal std.\n";
  siarCount = 0;
  for (auto const & sRitr : imgAzmRefMeans) {
    for (uint i=0; i<sRitr.size(); i++) {
      if (sRitr[i] != 0) {
        scanImgAzmRefSTD += std::pow(sRitr[i] - scanImgAzmRefMean, 2);
        siarCount++;
      }
    }
  }
  scanImgAzmRefSTD = std::sqrt(scanImgAzmRefSTD/siarCount);


  ///////////////////////////////////
  /////  Time Dependent Images  /////
  ///////////////////////////////////
 
  for (auto const & pItr : stagePosInds) {
    /////  Mean calculation  /////
    
    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal mean.\n";
    double siaCount = 0;
    for (auto const & sAitr : scanAzmAvg) {
      int scanInd = scanInds[sAitr.first];
      if (scanScale[sAitr.first][pItr.second] > 0) {
        double norm = 0;
        for (int i=0; i<NradAzmBins; i++) {
          if (sAitr.second[pItr.second][i] != NANVAL) {
            imgAzmMeans[scanInd][pItr.second] += sAitr.second[pItr.second][i];
            norm++;
          }
        }
        if (norm) {
          imgAzmMeans[scanInd][pItr.second] /= norm;
          scanImgAzmMeans[pItr.second] += imgAzmMeans[scanInd][pItr.second];
          siaCount++;
        }
        else {
          imgAzmMeans[scanInd][pItr.second] = 0;
        }
      }
    }

    scanImgAzmMeans[pItr.second] /= siaCount;


    /////  Standard deviation calculation  /////
    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal std.\n";
    siaCount = 0;
    for (auto const & sAitr : scanAzmAvg) {
      int scanInd = scanInds[sAitr.first];
      if (imgAzmMeans[scanInd][pItr.second]) {
        scanImgAzmSTDs[pItr.second]
            += std::pow((imgAzmMeans[scanInd][pItr.second]
                    - scanImgAzmMeans[pItr.second]), 2);
        siaCount++;
      }
    }

    scanImgAzmSTDs[pItr.second] = std::sqrt(scanImgAzmSTDs[pItr.second]/siaCount);
  }

  return;
}


void mergeClass::getRunMeanSTDSEM() {

  /////  Reset all vectors  /////
  // Time dependent measurements
  if (runAzmMeans.size() != stagePosInds.size()) {
    runLegMeans.resize(stagePosInds.size());
    runAzmMeans.resize(stagePosInds.size());
    runPCorrMeans.resize(stagePosInds.size());
    runLegSTD.resize(stagePosInds.size());
    runAzmSTD.resize(stagePosInds.size());
    runPCorrSTD.resize(stagePosInds.size());
    runLegSEM.resize(stagePosInds.size());
    runAzmSEM.resize(stagePosInds.size());
    runPCorrSEM.resize(stagePosInds.size());
    for (uint i=0; i<stagePosInds.size(); i++) {
      runLegMeans[i].resize(NlegBins, 0);
      runAzmMeans[i].resize(NradAzmBins, 0);
      runPCorrMeans[i].resize(maxRbins, 0);
      runLegSTD[i].resize(NlegBins, 0);
      runAzmSTD[i].resize(NradAzmBins, 0);
      runPCorrSTD[i].resize(maxRbins, 0);
      runLegSEM[i].resize(NlegBins, 0);
      runAzmSEM[i].resize(NradAzmBins, 0);
      runPCorrSEM[i].resize(maxRbins, 0);
    }
  }
  for (uint i=0; i<stagePosInds.size(); i++) {
    std::fill(runLegMeans[i].begin(), runLegMeans[i].end(), 0);
    std::fill(runAzmMeans[i].begin(), runAzmMeans[i].end(), 0);
    std::fill(runPCorrMeans[i].begin(), runPCorrMeans[i].end(), 0);
    std::fill(runLegSTD[i].begin(), runLegSTD[i].end(), 0);
    std::fill(runAzmSTD[i].begin(), runAzmSTD[i].end(), 0);
    std::fill(runPCorrSTD[i].begin(), runPCorrSTD[i].end(), 0);
  }
  std::vector<double>legNorms(NlegBins, 0);
  std::vector<double>azmNorms(NradAzmBins, 0);
  std::vector<double>pCorrNorms(maxRbins, 0);

  // Reference measurements
  if (!runAzmRefMean.size()) {
    runLegRefMean.resize(NlegBins, 0);
    runAzmRefMean.resize(NradAzmBins, 0);
    runLegRefSTD.resize(NlegBins, 0);
    runAzmRefSTD.resize(NradAzmBins, 0);
    runLegRefSEM.resize(NlegBins, 0);
    runAzmRefSEM.resize(NradAzmBins, 0);
  }
  std::fill(runLegRefMean.begin(), runLegRefMean.end(), 0);
  std::fill(runAzmRefMean.begin(), runAzmRefMean.end(), 0);
  std::fill(runLegRefSTD.begin(), runLegRefSTD.end(), 0);
  std::fill(runAzmRefSTD.begin(), runAzmRefSTD.end(), 0);

  double norm = 0;

  //////////////////////////////
  /////  Reference Images  /////
  //////////////////////////////

  /////  Calculate mean  /////
  /*
  // Legendre Reference
  if (verbose && (k==0)) 
    std::cout << "\tCalculating legendre mean.\n";
  for (int ir=0; ir<NlegBins; ir++) {
    norm = 0;
    for (auto sRitr : scanReferences) {
      for (auto pItr : sRitr.second) {
        if (pItr.second.legRef[ir] != NANVAL) {
          runLegRefMeans[ir] += pItr.second.legRef[ir];
          norm += 1;
        }
      }
    }
    if (norm) {
      runLegRefMeans[ir] /= norm;
    }
    else {
      runLegRefMeans[ir] = 0;
    }
  }
  */

  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal mean.\n";
  for (int ir=0; ir<NradAzmBins; ir++) {
    norm = 0;
    for (auto const & sRitr : scanReferences) {
      for (auto const & pItr : sRitr.second) {
        if (pItr.second.scale) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            runAzmRefMean[ir] += pItr.second.azmRef[ir];
            norm += 1;
          }
        }
      }
    }
    if (norm) {
      runAzmRefMean[ir] /= norm;
    }
    else {
      runAzmRefMean[ir] = 0;
    }
  }

  /////  Calculate Standard Deviation  /////
  /*
  // Legendre Reference
  if (verbose && (k==0)) 
    std::cout << "\tCalculating legendre std.\n";
  for (int ir=0; ir<NlegBins; ir++) {
    norm = 0;
    for (auto const & sRitr : scanReferences) {
      for (auto const & pItr : sRitr.second) {
        if (pItr.second.legRef[ir] != NANVAL) {
          runLegRefSTD[ir] 
              += std::pow((pItr.second.legRef[ir] 
                      - runLegRefMeans[ir]), 2);
          norm += 1;
        }
      }
    }
    if (norm) {
      runLegRefSTD[ir] = std::sqrt(runLegRefSTD[ir]/norm);
    }
    else {
      runLegRefSTD[ir] = 0;
    }
  }
  */

  // Azimuthal Reference
  if (verbose) 
    std::cout << "\tCalculating azimuthal std.\n";
  for (int ir=0; ir<NradAzmBins; ir++) {
    norm = 0;
    for (auto const & sRitr : scanReferences) {
      for (auto const & pItr : sRitr.second) {
        if (pItr.second.scale) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            runAzmRefSTD[ir] 
                += std::pow((pItr.second.azmRef[ir]
                        - runAzmRefMean[ir]), 2);
            norm += 1;
          }
        }
      }
    }
    if (norm - 1 > 0) {
      runAzmRefSTD[ir] = std::sqrt(runAzmRefSTD[ir]/(norm - 1));
      runAzmRefSEM[ir] = runAzmRefSTD[ir]/std::sqrt(norm);
    }
    else {
      runAzmRefSTD[ir] = 0;
      runAzmRefSEM[ir] = 0;
    }
  }


  ///////////////////////////////////
  /////  Time Dependent Images  /////
  ///////////////////////////////////
 
  for (auto const & pItr : stagePosInds) {
    /////  Mean calculation  /////
    std::fill(azmNorms.begin(), azmNorms.end(), 0);
    std::fill(legNorms.begin(), legNorms.end(), 0);
    
    /*
    ///  Legendres  ///
    if (verbose && (k == 0) && (pItr.second == 0)) 
      std::cout << "\tCalculating legendre mean.\n";
    for (auto sLitr : scanLgndrs) {
      if (scanScale[sLitr.first][pItr.second] > 0) {
        for (int i=0; i<NlegBins; i++) {
          if (sLitr.second[pItr.second][i] != NANVAL) {
            runLegMeans[pItr.second][i] += sLitr.second[pItr.second][i];
            legNorms[i] += 1; //scanScale[sLitr.first][pItr.second];
          }
        }
      }
    }

    for (int i=0; i<NlegBins; i++) {
      runLegMeans[pItr.second][i] /= legNorms[i];
    }
    */

    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal mean.\n";
    for (auto const & sAitr : scanAzmAvg) {
      if (scanScale[sAitr.first][pItr.second] != 0) {
        for (int i=0; i<NradAzmBins; i++) {
          if (sAitr.second[pItr.second][i] != NANVAL) {
            runAzmMeans[pItr.second][i] += sAitr.second[pItr.second][i];
            azmNorms[i] += 1; //scanScale[sAitr.first][pItr.second];
          }
        }
      }
    }

    for (int i=0; i<NradAzmBins; i++) {
      runAzmMeans[pItr.second][i] /= azmNorms[i];
      if (didSubtractT0) {
        runAzmMeans[pItr.second][i] -= runAzmRefMean[i];
      }
    }


    for (auto const & sPitr : scanAzmPCorr) {
      if (scanScale[sPitr.first][pItr.second] > 0) {
        for (int ir=0; ir<maxRbins; ir++) {
          runPCorrMeans[pItr.second][ir] += sPitr.second[pItr.second][ir];
          pCorrNorms[ir] += 1; 
        }
      }
    }

    for (int ir=0; ir<maxRbins; ir++) {
      runPCorrMeans[pItr.second][ir] /= pCorrNorms[ir];
    }

    /////  Standard deviation calculation  /////

    /*
    ///  Legendres  ///
    if (verbose && (k == 0) && (pItr.second == 0)) 
      std::cout << "\tCalculating legendre std.\n";
    for (auto const & sLitr : scanLgndrs) {
      if (scanScale[sLitr.first][pItr.second]) {
        for (int i=0; i<NlegBins; i++) {
          if (sLitr.secondpItr.second][i] != NANVAL) {
            runLegSTD[pItr.second][i]
                += std::pow((sLitr.second[pItr.second][i]
                        - runLegMeans[pItr.second][i]), 2);
            
          }
        }
      }
    }
    

    for (int i=0; i<NlegBins; i++) {
      runLegSTD[pItr.second][i] = std::sqrt(runLegSTD[pItr.second][i]/legNorms[i]);
    }
    */

    ///  Azimuthal  ///
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tCalculating azimuthal std.\n";
    for (auto const & sAitr : scanAzmAvg) {
      if (scanScale[sAitr.first][pItr.second] != 0) {
        for (int i=0; i<NradAzmBins; i++) {
          if (sAitr.second[pItr.second][i] != NANVAL) {
            if (didSubtractT0) {
              runAzmSTD[pItr.second][i]
                  += std::pow(((sAitr.second[pItr.second][i] - runAzmRefMean[i])
                        - runAzmMeans[pItr.second][i]), 2);
              if (pItr.second == 2 && i == 19) {
                //cout<<"scan "<<sAitr.first<<" : "<<std::pow(((sAitr.second[pItr.second][i] - runAzmRefMeans[i]) - runAzmMeans[pItr.second][i]), 2)<<"   "<<runAzmMeans[pItr.second][i]<<endl;
              }
            }
            else {
              runAzmSTD[pItr.second][i]
                  += std::pow((sAitr.second[pItr.second][i]
                        - runAzmMeans[pItr.second][i]), 2);
              if (pItr.second == 2 && i == 19) {
                //cout<<"scan "<<sAitr.first<<" : "<<std::pow(((sAitr.second[pItr.second][i]) - runAzmMeans[pItr.second][i]), 2)<<"   "<<runAzmMeans[pItr.second][i]<<endl;
              }

            }
          }
        }
      }
    }

    for (int i=0; i<NradAzmBins; i++) { 
      if (azmNorms[i] - 1 > 0) {
        runAzmSTD[pItr.second][i] 
            = std::sqrt(runAzmSTD[pItr.second][i]/(azmNorms[i] - 1));
        runAzmSEM[pItr.second][i] 
            = runAzmSTD[pItr.second][i]/std::sqrt(azmNorms[i]);
      }
      else {
        runAzmSTD[pItr.second][i] = 0;
        runAzmSEM[pItr.second][i] = 0;
      }
    }

    if (didPairCorrSTD) {
      for (auto const & sPitr : scanAzmPCorr) {
        if (scanScale[sPitr.first][pItr.second] != 0) {
          for (int ir=0; ir<maxRbins; ir++) {
            runPCorrSTD[pItr.second][ir]
                  += std::pow((sPitr.second[pItr.second][ir]
                          - runPCorrMeans[pItr.second][ir]), 2);
          }
        }
      }

      for (int ir=0; ir<maxRbins; ir++) {
        if (pCorrNorms[ir] - 1 > 0) {
          runPCorrSTD[pItr.second][ir] 
                = std::sqrt(runPCorrSTD[pItr.second][ir]/(pCorrNorms[ir] - 1));
          runPCorrSEM[pItr.second][ir] = runPCorrSTD[pItr.second][ir]/std::sqrt(pCorrNorms[ir]);
        }
        else {
          runPCorrSTD[pItr.second][ir] = 0;
          runPCorrSEM[pItr.second][ir] = 0;
        }

      }
    }
  }

  SEMisBootstrap = false;

  return;
}


void mergeClass::bootstrapSEM() {

  std::srand(time(NULL));
  int NtimeSteps = stagePosInds.size();
  std::vector<double> runAzmRefBstMeanCounts(NradAzmBins, 0);
  std::vector<double> runAzmRefBstSEMCounts(NradAzmBins, 0);
  std::vector< std::vector<double> > runAzmBstMeanCounts(NtimeSteps);
  std::vector< std::vector<double> > runAzmBstSEMCounts(NtimeSteps);
  std::vector< std::vector<double> > runAzmRefBstDist(mergeNbootstrap);
  std::vector< std::vector<double> > runAzmRefBstDistCounts(mergeNbootstrap);
  std::vector< std::vector< std::vector<double> > > runAzmBstDist(mergeNbootstrap);
  std::vector< std::vector< std::vector<double> > > runAzmBstDistCounts(mergeNbootstrap);
  runLegSEM.clear();     runLegSEM.resize(NtimeSteps);
  runAzmSEM.clear();     runAzmSEM.resize(NtimeSteps);
  runLegMeans.clear();   runLegMeans.resize(NtimeSteps);
  runAzmMeans.clear();   runAzmMeans.resize(NtimeSteps);
  runLegRefMean.clear(); runLegRefMean.resize(NradLegBins);
  runAzmRefMean.clear(); runAzmRefMean.resize(NradAzmBins);
  runLegRefSEM.clear();  runLegRefSEM.resize(NradLegBins);
  runAzmRefSEM.clear();  runAzmRefSEM.resize(NradAzmBins);

  std::vector<int> scans;
  for (int it=0; it<NtimeSteps; it++) {
    for (auto const & sItr : scanScale) {
      if (sItr.second[it] != 0) {
        if (std::find(scans.begin(), scans.end(), sItr.first) 
            == scans.end()) {
          scans.push_back(sItr.first);
        }
      }
    }

    runLegSEM[it].resize(NradLegBins, 0);
    runAzmSEM[it].resize(NradAzmBins, 0);
    runLegMeans[it].resize(NradLegBins, 0);
    runAzmMeans[it].resize(NradAzmBins, 0);
    runAzmBstMeanCounts[it].resize(NradAzmBins, 0);
    runAzmBstSEMCounts[it].resize(NradAzmBins, 0);
  }
  int Nscans = (int)scans.size();

  /////  Selecting Ensembles  /////

  for (int bst=0; bst<mergeNbootstrap; bst++) {
    if (bst % 500 == 0) {
      std::cout << "Bootstrap: " << bst << std::endl;
    }

    runAzmBstDist[bst].resize(NtimeSteps);
    runAzmBstDistCounts[bst].resize(NtimeSteps);
    runAzmRefBstDist[bst].resize(NradAzmBins);
    runAzmRefBstDistCounts[bst].resize(NradAzmBins);
    for (int it=0; it<NtimeSteps; it++) {
      runAzmBstDist[bst][it].resize(NradAzmBins, 0);
      runAzmBstDistCounts[bst][it].resize(NradAzmBins, 0);
    }

    for (int isc=0; isc<Nscans; isc++) {
      int scan = scans[(std::rand() % Nscans)];

      // References
      for (auto const & sRitr : scanReferences[scan]) {
        if (sRitr.second.scale) {
          for (int iazm=0; iazm<NradAzmBins; iazm++) {
            if (sRitr.second.azmRef[iazm] != NANVAL) {
              runAzmRefBstDist[bst][iazm] += sRitr.second.azmRef[iazm];
              runAzmRefBstDistCounts[bst][iazm] += 1;
            }
          }
        }
      }

      // Time Dependent Data
      for (int it=0; it<NtimeSteps; it++) {
        if (scanScale[scan][it] != 0) {
          for (int iazm=0; iazm<NradAzmBins; iazm++) {
            if (scanAzmAvg[scan][it][iazm] != NANVAL) {
              runAzmBstDist[bst][it][iazm] += scanAzmAvg[scan][it][iazm];
              runAzmBstDistCounts[bst][it][iazm] += 1;
            }
          }
        }
      }
    }

    for (int iazm=0; iazm<NradAzmBins; iazm++) {
      if (runAzmRefBstDistCounts[bst][iazm] > 0) {
        runAzmRefBstDist[bst][iazm]   /= runAzmRefBstDistCounts[bst][iazm];
      }
      else {
        runAzmRefBstDist[bst][iazm] = NANVAL;
      }

      for (int it=0; it<NtimeSteps; it++) {
        if (runAzmBstDistCounts[bst][it][iazm] > 0) {
          runAzmBstDist[bst][it][iazm] /= runAzmBstDistCounts[bst][it][iazm];
        }
        else {
          runAzmBstDist[bst][it][iazm] = NANVAL;
        }
        if (didSubtractT0) {
          if ((runAzmBstDist[bst][it][iazm] != NANVAL)
              && (runAzmRefBstDist[bst][iazm] != NANVAL)) {
            runAzmBstDist[bst][it][iazm] -= runAzmRefBstDist[bst][iazm];
          }
          else {
            runAzmBstDist[bst][it][iazm] = NANVAL;
          }
        }
      }
    }
    //plt->print1d(runAzmRefBstDist[bst], "testingBstAzmRef_"+to_string(bst));
    //plt->print1d(runAzmBstDist[bst][5], "testingBstAzm_"+to_string(bst));


    // Calculating mean
    for (int iazm=0; iazm<NradAzmBins; iazm++) {
      if (runAzmRefBstDist[bst][iazm] != NANVAL) {
        runAzmRefMean[iazm] += runAzmRefBstDist[bst][iazm];
        runAzmRefBstMeanCounts[iazm] += 1;
      }

      for (int it=0; it<NtimeSteps; it++) {
        if (runAzmBstDist[bst][it][iazm] != NANVAL) {
          runAzmMeans[it][iazm] += runAzmBstDist[bst][it][iazm];
          runAzmBstMeanCounts[it][iazm] += 1;
        }
      }
    }
  }

  // Normalizing mean
  for (int iazm=0; iazm<NradAzmBins; iazm++) {
    if (runAzmRefBstMeanCounts[iazm] > 0) {
      runAzmRefMean[iazm] /= runAzmRefBstMeanCounts[iazm];
    }
    else {
      runAzmRefMean[iazm] = NANVAL;
    }
    for (int it=0; it<NtimeSteps; it++) {
      if (runAzmBstMeanCounts[it][iazm] > 0) {
        runAzmMeans[it][iazm] /= runAzmBstMeanCounts[it][iazm];
      }
      else {
        runAzmMeans[it][iazm] = NANVAL;
      }
    }
  }
  plt->print1d(runAzmRefMean, "testingBstAzmRefMean");
  plt->print1d(runAzmMeans[5], "testingBstAzmMean");

  ////  Calculating STD of distribution of means: SEM  /////
  // Sum square differences
  for (int bst=0; bst<mergeNbootstrap; bst++) {
    for (int iazm=0; iazm<NradAzmBins; iazm++) {
      if ((runAzmRefBstDist[bst][iazm] != NANVAL)
          && (runAzmRefMean[iazm] != NANVAL)) {
        runAzmRefSEM[iazm] += std::pow(
            runAzmRefBstDist[bst][iazm] - runAzmRefMean[iazm], 2);
        runAzmRefBstSEMCounts[iazm] += 1;
      }

      for (int it=0; it<NtimeSteps; it++) {
        if ((runAzmBstDist[bst][it][iazm] != NANVAL)
            && (runAzmMeans[it][iazm] != NANVAL)) {
          runAzmSEM[it][iazm] += std::pow(
              runAzmBstDist[bst][it][iazm] - runAzmMeans[it][iazm], 2);
          runAzmBstSEMCounts[it][iazm] += 1;
        }
      }
    }
  }

  // Normalizing and taking sqrt
  for (int iazm=0; iazm<NradAzmBins; iazm++) {
    if (runAzmRefBstSEMCounts[iazm] > 0) {
      runAzmRefSEM[iazm] = std::sqrt(
          runAzmRefSEM[iazm]/runAzmRefBstSEMCounts[iazm]);
    }
    else {
      runAzmRefSEM[iazm] = NANVAL;
    }

    for (int it=0; it<NtimeSteps; it++) {
      if (runAzmBstSEMCounts[it][iazm] > 0) {
        runAzmSEM[it][iazm] = std::sqrt(
            runAzmSEM[it][iazm]/runAzmBstSEMCounts[it][iazm]);
      }
      else {
        runAzmSEM[it][iazm] = NANVAL;
      }
    }
  }

  SEMisBootstrap = true;

  return;
}


void mergeClass::testSEMbootstrap() {

  int Ntests = 4;
  std::vector<int> testN(Ntests);
  std::vector< std::vector<double> > testSEMvals(Ntests+1);

  testN[0] = mergeNbootstrap;
  testSEMvals[1].assign(
      runAzmSEM[0].begin(),
      runAzmSEM[0].end());
  getRunMeanSTDSEM();
  testSEMvals[0].assign(
      runAzmSEM[0].begin(),
      runAzmSEM[0].end());

  int NbootstrapOrig = mergeNbootstrap;
  for (int i=0; i<Ntests-1; i++) {
    mergeNbootstrap = NbootstrapOrig*(i+1)*5;
    bootstrapSEM();
    testSEMvals[i+2].assign(
        runAzmSEM[0].begin(),
        runAzmSEM[0].end());
    testN[i+1] = mergeNbootstrap;
  }

  ofstream outFile;
  outFile.open("./results/mergeNbootstrapTest.txt");

  outFile << "Ordering: \tSEM/" << testN[0];
  for (int i=1; i<Ntests; i++) {
    outFile << "\t\t" << testN[i-1] << "/" << testN[i];
  }
  outFile << "\n\n";
  for (int i=0; i<NradAzmBins; i++) {
    outFile << "Pixel " << i << ":";
    for (int j=1; j<(int)testSEMvals.size(); j++) {
      if ((testSEMvals[j-1][i] != 0) && (testSEMvals[j][i] != 0)) {
        outFile << "   \t" << testSEMvals[j-1][i]/testSEMvals[j][i];
      }
      else {
        outFile << "\t\t";
      }
    }
    outFile << std::endl;
  }
  outFile.close();

  std::cerr << "EXITING AFTER TESTING mergeNbootstrap!!!\n";
  exit(0);

  return;
}


void mergeClass::removeLowPolynomials() {

  if (!_compareReference) {
    mergeScans(true, false);
  }
  getRunMeanSTDSEM();

  int mInd, Nnans;
  int Npoly = 1;
  int size = NradAzmBins - NbinsSkip;
  double delta = maxQazm/(NradAzmBins-1);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;
  Eigen::VectorXd w;

  /////  Subtracting from references  /////
  if (verbose)
    std::cout << "\tRemoving low order polynomials from references.\n";
  for (auto& sRitr : scanReferences) {
    for (auto& pItr : sRitr.second) {
      if (pItr.second.scale) {
        Nnans = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] == NANVAL) Nnans += 1;
        }
        X.resize(size-Nnans, Npoly);
        Y.resize(size-Nnans, 1);
        w.resize(size-Nnans);
        mInd = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            X(mInd,0) = std::pow(maxQazm-ir*delta, NlowOrderPoly);
            if (_compareReference) {
              Y(mInd,0) = pItr.second.azmRef[ir] - compReference[ir];
            }
            else {
              Y(mInd,0) = pItr.second.azmRef[ir] - azmReference[ir];
            }
            w(mInd)   = std::pow(1/runAzmRefSTD[ir], 2);
            mInd += 1;
          }
        }

        Eigen::MatrixXd weights = tools::normalEquation(X, Y, w);

        for (int ir=0; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            pItr.second.azmRef[ir] -= weights(0)
                *std::pow(maxQazm-ir*delta, NlowOrderPoly);
          }
        }

      }
    }
  }

  /////  Subtracting from time dependend images  /////
  if (verbose)
    std::cout << "\tRemoving low order polynomials from td images.\n";

  if (!_compareReference) {
    mergeScans(true, false);
  }

  std::vector<PLOToptions> opts(2);
  std::vector<string> vals(2);
  opts[0] = maximum;    vals[0] = "0.3";
  opts[1] = minimum;    vals[1] = "-0.3";
  std::vector<double> plotme1(NradAzmBins);
  std::vector<double> plotme2(NradAzmBins);
  std::vector<double> fit(NradAzmBins);
  std::vector<TH1*> h(2);
  for (uint it=0; it<stagePosInds.size(); it++) {
    for (auto& sAzml : scanAzmAvg) {
      if (scanScale[sAzml.first][it] != 0) {
        Nnans = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if (sAzml.second[it][ir] == NANVAL) Nnans += 1;
        }
        X.resize(size-Nnans, Npoly);
        Y.resize(size-Nnans, 1);
        w.resize(size-Nnans);
        mInd = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if ((sAzml.second[it][ir] != NANVAL) && (azmReference[ir] != NANVAL)) {
            X(mInd,0) = std::pow(maxQazm-ir*delta, NlowOrderPoly);
            if (_compareReference) {
              Y(mInd,0) = sAzml.second[it][ir] - compReference[ir];
            }
            else {
              Y(mInd,0) = sAzml.second[it][ir] - azmReference[ir];
            }
            w(mInd)   = std::pow(1/runAzmSTD[it][ir], 2);
            mInd += 1;
          }
        }

        Eigen::MatrixXd weights = tools::normalEquation(X, Y, w);

        std::fill(fit.begin(), fit.end(), 0);
        for (int ir=0; ir<NradAzmBins; ir++) {
          if (sAzml.second[it][ir]!= NANVAL) {

            // Plot results for debugging/validation
            if (lowPolySubtractStudy) {
              if (_compareReference) {
                plotme1[ir] = sAzml.second[it][ir] - compReference[ir];
              }
              else {
                plotme1[ir] = sAzml.second[it][ir] - azmReference[ir];
              }
            }

            // Subtract fit
            if (sAzml.second[it][ir] != NANVAL) {
              sAzml.second[it][ir] -= weights(0)
                  *std::pow(maxQazm-ir*delta, NlowOrderPoly);
            }
            fit[ir] += weights(0)*std::pow(maxQazm-ir*delta, NlowOrderPoly);

            // Plot results for debugging/validation
            if (lowPolySubtractStudy) {
              if (_compareReference) {
                plotme2[ir] = sAzml.second[it][ir] - compReference[ir];
              }
              else {
                plotme2[ir] = sAzml.second[it][ir] - azmReference[ir];
              }
            }
          }
          else if (lowPolySubtractStudy) {
            plotme1[ir] = 0;
            plotme2[ir] = 0;
          }
        }

        if (lowPolySubtractStudy) {
          h[0] = plt->plot1d(plotme1, "testLowPolyRemoval-"
              + to_string(sAzml.first) + "-" + to_string(it),
              opts, vals);
          h[1] = plt->plot1d(fit, "testLowPolyRemoval-"
              + to_string(sAzml.first) + "-" + to_string(it) + "-fit",
              opts, vals);
          plt->print1d(h, "./plots/testLowPolyRemoval-"
              + to_string(sAzml.first) + "-" + to_string(it) + "-orig",
              opts, vals);
          plt->print1d(plotme2, "./plots/testLowPolyRemoval-"
              + to_string(sAzml.first) + "-" + to_string(it) + "-final",
              opts, vals);
          delete h[0];
          delete h[1];
        }
      }
    }
  }

  return;
}


void mergeClass::removeImgNormOutliers() {
  
  if (verbose)
    std::cout << "\tBeginning removing outliers in reference images.\n";

  /////  Calculate Mean and STD  /////
  getImgNormMeanSTD();

  /////  Make Reference Cut  /////
  if (verbose) 
    std::cout << "\tMaking cuts on reference images.\n";
  for (auto& sRitr : scanReferences) {
    for (auto& pItr : sRitr.second) {
      if (fabs(pItr.second.imgNorm - scanImgNormAzmRefMean) 
          > scanImgAzmRefSTDcut*scanImgNormAzmRefSTD) {
        pItr.second.scale = 0;
      }
    }
  }
    
  /////  Time Dependent  /////
  for (auto pItr : stagePosInds) {
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tMaking cuts on azimuthal images.\n";
    for (auto& sImgN : scanImgNorms) {
      if (fabs(sImgN.second[pItr.second] - scanImgNormAzmMeans[pItr.second])
          > scanImgAzmSTDcut*scanImgNormAzmSTDs[pItr.second]) {
        scanScale[sImgN.first][pItr.second] = 0;
      }
    }
  }

  if (verbose) 
    std::cout << "Finished removing outliers from time dependent images.\n";

  return;
}

 
void mergeClass::removeImageOutliers() {
  
  if (verbose)
    std::cout << "\tBeginning removing outliers in reference images.\n";

  /////  Calculate Mean and STD  /////
  getImageMeanSTD();

  /////  Make Reference Cut  /////
  if (verbose) 
    std::cout << "\tMaking cuts on reference images.\n";
  for (auto& sRitr : scanReferences) {
    int posInd = 0;
    for (auto& pItr : sRitr.second) {
      if (fabs(imgAzmRefMeans[scanInds[sRitr.first]][posInd] - scanImgAzmRefMean) 
          > scanImgAzmRefSTDcut*scanImgAzmRefSTD) {
        pItr.second.scale = 0;
      }
      posInd++;
    }
  }
     
  /////  Time Dependent  /////
  for (auto pItr : stagePosInds) {
    if (verbose && (pItr.second == 0)) 
      std::cout << "\tMaking cuts on azimuthal images.\n";
    for (auto& sAitr : scanAzmAvg) {
      if (fabs(imgAzmMeans[scanInds[sAitr.first]][pItr.second] 
          - scanImgAzmMeans[pItr.second])
          > scanImgAzmSTDcut*scanImgAzmSTDs[pItr.second]) {
        scanScale[sAitr.first][pItr.second] = 0;
        cout<<"removing image Outlier: "<<sAitr.first<<"  "<<pItr.first<<endl;
      }
    }
  }

  if (verbose) 
    std::cout << "Finished removing outliers from time dependent images.\n";

  return;
}


 

void mergeClass::removeOutliers() {

  /////  Remove time bins with few images and get mean  /////
  std::vector<int32_t> removePos;

  // End of initial NAN values from detector hole
  int legNANend = 0;
  int azmNANend = 0;
  for (auto i : scanLgndrs.begin()->second[0]) {
    if (i == NANVAL) {
      legNANend++;
    }
    else {
      break;
    }
  }
  for (auto i : scanAzmAvg.begin()->second[0]) {
    if (i == NANVAL) {
      azmNANend++;
    }
    else {
      break;
    }
  }

  
  ///////////////////////////////////////
  /////  Cleaning Reference Images  /////
  ///////////////////////////////////////
  if (verbose)
    std::cout << "\tBeginning removing outliers in reference images.\n";
  //FIX ME CHANGED k<3 to k<1
  for (int k=0; k<1; k++) {

    /////  Calculate Mean and STD  /////
    getRunMeanSTDSEM();

    /////  Make Reference Cut  /////
    /*
    // Legendre Reference
    if (verbose && (k==0)) 
      std::cout << "\tMaking cuts on legendre.\n";
    for (auto& sRitr : scanReferences) {
      for (auto& pItr : sRitr.second) {
        int imageNoise = 0;
        for (int ir=0; ir<NlegBins; ir++) {
          if (pItr.second.legRef[ir] != NANVAL) {
            if ((fabs(pItr.second.legRef[ir] - runLegRefMeans[ir])
                > mergeImageSTDScale*runLegRefSTD[ir]) 
                || ((pItr.second.legRef[ir] == NANVAL) && (ir >= legNANend))) {
              imageNoise++;
            }      
            if (fabs(pItr.second.legRef[ir] - runLegRefMeans[ir])
                > mergeSTDscale*runLegRefSTD[ir]) {
              pItr.second.legRef[ir] = NANVAL;
            }      
          }
        }
        if (imageNoise > legImageNoiseCut) {
          pItr.second.scale = 0;
        }
      }
    }
    */

    // Azimuthal Reference
    if (verbose && (k==0)) 
      std::cout << "\tMaking cuts on azimuthal.\n";
    std::vector< std::vector<double> > tempV(3);
    for (auto& sRitr : scanReferences) {
      for (auto& pItr : sRitr.second) {
        int imageNoise = 0;
        for (int ir=0; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            if ((fabs(pItr.second.azmRef[ir] - runAzmRefMean[ir])
                > mergeImageSTDScale*runAzmRefSTD[ir])
                || ((pItr.second.azmRef[ir] == NANVAL) && (ir >= azmNANend))) {
              imageNoise++;
            }      
            if (fabs(pItr.second.azmRef[ir] - runAzmRefMean[ir])
                > mergeSTDscale*runAzmRefSTD[ir]) {
              pItr.second.azmRef[ir] = NANVAL;
            }    
            else if (ir == 300) {
              tempV[0].push_back(pItr.second.azmRef[ir]);
            }
            else if (ir == 243) {
              tempV[1].push_back(pItr.second.azmRef[ir]);
            }
            else if (ir == 188) {
              tempV[2].push_back(pItr.second.azmRef[ir]);
            }
          }
        }
        if (imageNoise > azmImageNoiseCut) {
          pItr.second.scale = 0;
        }
      }
    }

    for (uint i=0; i<tempV.size(); i++) {
      //save::saveDat<double>(tempV[i], "testDistro" + to_string(i) + "["+to_string(tempV[i].size()) + "].dat");
    }
    

    if (verbose) {
      std::cout << "Finished removing reference outliers.\n";
      std::cout << "Begin removing outliers from time dependent images.\n";
    }

    for (auto pItr : stagePosInds) {
      /////  Make cuts  /////
      /*
      ///  Legendres  ///
      if (verbose && (k == 0) || (pItr.second == 0)) 
        std::cout << "\tMaking cuts on legendre images.\n";
      for (auto& sLitr : scanLgndrs) {
        int imageNoise = 0;
        for (int ir=0; ir<NlegBins; ir++) {
          if ((fabs(sLitr.second[pItr.second][ir] - runLegMeans[pItr.second][ir])
              > mergeImageSTDScale*runLegSTD[pItr.second][ir]) 
              || ((sLitr.second[pItr.second][ir] == NANVAL) && (ir >= legNANend))) {
            imageNoise++;
          }
          if (fabs(sLitr.second[pItr.second][ir] - runLegMeans[pItr.second][ir])
              > mergeSTDscale*runLegSTD[pItr.second][ir]) {
            sLitr.second[pItr.second][ir] = NANVAL;
          }
        }
        if (imageNoise > legImageNoiseCut) {
          scanScale[sLitr.first][pItr.second] = 0;
        }
      }
      */


      ///  Azimuthal  ///
      if (verbose && (k == 0) && (pItr.second == 0)) 
        std::cout << "\tMaking cuts on azimuthal images.\n";
      for (auto& sAitr : scanAzmAvg) {
        int imageNoise = 0;
        for (int ir=0; ir<NradAzmBins; ir++) {
          if ((fabs(sAitr.second[pItr.second][ir] - runAzmMeans[pItr.second][ir])
              > mergeImageSTDScale*runAzmSTD[pItr.second][ir]) 
              || ((sAitr.second[pItr.second][ir] == NANVAL) && (ir >= azmNANend))) {
            imageNoise++;
          }
          if (fabs(sAitr.second[pItr.second][ir] - runAzmMeans[pItr.second][ir])
              > mergeSTDscale*runAzmSTD[pItr.second][ir]) {
            sAitr.second[pItr.second][ir] = NANVAL;
          }
        }
        if (imageNoise > azmImageNoiseCut) {
          scanScale[sAitr.first][pItr.second] = 0;
        }
      }
    }
  }

  if (verbose) 
    std::cout << "Finished removing outliers from time dependent images.\n";

  //plt.printXY(runMeans, "testRunMeans", maximum, "5e2");

  if (verbose) {
    std::cout << "Removing " << removePos.size() << " bins\n";
    std::cout << "Before removal: " << stagePosInds.size() << "  "
      << scanLgndrs.begin()->second.size() << "  "
      << scanScale.begin()->second.size() << "  "
      << runLegMeans.size() << "  " << runLegSTD.size() << std::endl;
  }
  /*
  ///  Remove time bins  ///
  for (int i=((int)removePos.size())-1; i>=0; i--) {
    int pInd = stagePosInds[removePos[i]];
    for (auto& sLitr : scanLgndrs) {
      sLitr.second.erase(sLitr.second.begin() + pInd);
    }
    for (auto& sCitr : scanScale) {
      sCitr.second.erase(sCitr.second.begin() + pInd);
    }
    stagePosInds.erase(stagePosInds.find(removePos[i]));

    runMeans.erase(runMeans.begin() + pInd);
    runStdev.erase(runStdev.begin() + pInd);
  }
  int pCount = 0;
  for (auto& pItr : stagePosInds) {
    pItr.second = pCount;
    pCount++;
  }
  if (verbose) {
    std::cout << "After removal: " << stagePosInds.size() << "  " 
      << scanLgndrs.begin()->second.size() << "  " 
      << scanScale.begin()->second.size() << "  "
      << runMeans.size() << "  " << runStdev.size() << std::endl;
  }
  */



}


void mergeClass::scaleByFit() {

  if (verbose) {
    std::cout << "INFO: scaling by fitting patterns to mean\n";
  }

  Eigen::Matrix<double, Eigen::Dynamic, 1> X_full;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y_full;
  X_full.resize(NradAzmBins, 1);
  Y_full.resize(NradAzmBins, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;
  Eigen::VectorXd w;

  mergeScans();

  /////  Scaling References  /////
  for (auto& sRitr : scanReferences) {
    for (auto& pItr : sRitr.second) {
      if (pItr.second.scale) {
        int count = 0;
        for (int ir=100; ir<450; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            X_full(count, 0) = pItr.second.azmRef[ir];
            Y_full(count, 0) = azmReference[ir];
            count++;
          }
        }

        X.resize(count, 1);
        Y.resize(count, 1);
        for (int i=0; i<count; i++) {
          X(i, 0) = X_full(i, 0);
          Y(i, 0) = Y_full(i, 0);
        }

        w = tools::normalEquation(X, Y);

        // Scale lineout
        for (int ir=0; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            pItr.second.azmRef[ir] = w(0)*pItr.second.azmRef[ir];
          }
        }
      }
    }
  }

  /*
  cout<<"Full\n";
          for (int k=100; k<105; k++) {
            cout<<azimuthalAvg[2][k]<<"  ";
          }
          cout<<endl;
          */

  /////  Scaling Azimuthal Average  /////
  for (uint it=0; it<stagePosInds.size(); it++) {
    for (auto & sAzml : scanAzmAvg) {
      if (scanScale[sAzml.first][it] != 0) {
        int count = 0;
        for (int iazm=100; iazm<450; iazm++) {
          if (sAzml.second[it][iazm] != NANVAL)  {
            X_full(count, 0) = sAzml.second[it][iazm];
            Y_full(count, 0) = azimuthalAvg[it][iazm];
            count++;
          }
        }

        X.resize(count, 1);
        Y.resize(count, 1);
        for (int i=0; i<count; i++) {
          X(i, 0) = X_full(i, 0);
          Y(i, 0) = Y_full(i, 0);
        }

        w = tools::normalEquation(X, Y);

        /*
        if (it == 2) {
          cout<<"fitting "<<sAzml.first<<" : "<<w(1)<<"  "<<w(0)<<endl;
          for (int k=100; k<105; k++) {
            cout<<sAzml.second[it][k]<<"  ";
          }
          cout<<endl;
        }
        */
        // Scale lineout
        for (int iazm=0; iazm<NradAzmBins; iazm++) {
          if (sAzml.second[it][iazm] != NANVAL) {
            sAzml.second[it][iazm] = w(0)*sAzml.second[it][iazm];
          }
        }
      }
    }
  }

  return;
}
   


void mergeClass::mergeScans(bool refOnly, bool tdOnly) {

  /////////////////////////////////////////////////////////
  /////  Calculating time delays from stage position  /////
  /////////////////////////////////////////////////////////

  if (verbose) {
    std::cout << "INFO: Inside mergeScans\n";
    std::cout << "\tcalculating time delays: " 
      << stagePosInds.size() << std::endl;
  }

  int NtimeSteps = stagePosInds.size();
  int tInd = 1;
  timeDelays.resize(stagePosInds.size()+1);
  for (auto pItr : stagePosInds) {
    timeDelays[tInd] = 2*pItr.first/(3e3);
    tInd++;
  }
  timeDelays[0] = 2*timeDelays[1] - timeDelays[2];
  for (int i=NtimeSteps; i>=0; i--) {
    timeDelays[i] -= timeZero + timeDelays[0];
  }

  save::saveDat<double>(
      timeDelays,
      "./results/timeDelays-" + run + "_bins["
        + to_string(timeDelays.size()) + "].dat");


  ///////////////////////////
  /////  Merging scans  /////
  ///////////////////////////

  int rInd;
  double norm;

  /////  Merging reference images  /////
  if (refOnly || !tdOnly) {
    if (verbose)
      std::cout << "\tMerging legendres reference.\n";
    ///  Merging legendres  ///
    legReference.clear();
    legReference.resize(Nlegendres);
    for (int ilg=0; ilg<Nlegendres; ilg++) {
      legReference[ilg].resize(NradLegBins, 0);
      for (int ir=0; ir<NradLegBins; ir++) {
        norm = 0;
        rInd = ilg*NradLegBins + ir;
        for (auto const & sRitr : scanReferences) {
          for (auto const & pItr : sRitr.second) {
            if (pItr.second.scale) {
              if (pItr.second.legRef[rInd] != NANVAL) {
                legReference[ilg][ir] += pItr.second.legRef[rInd];
                norm += 1;
              }
            }
          }
        }
        if (norm) {
          legReference[ilg][ir] /= norm;
        }
        else {
          legReference[ilg][ir] = 0;
        }
      }
    }

    ///  Merging azimuthal averages  ///
    if (verbose)
      std::cout << "\tMerging azimuthal average reference.\n";
    azmReference.resize(NradAzmBins, 0);
    std::fill(azmReference.begin(), azmReference.end(), 0);
    std::vector< std::vector<double> > azmRefSplit(2);
    std::vector< std::vector<double> > azmRefSplitCount(2);
    for (uint k=0; k<azmRefSplit.size(); k++) {
      azmRefSplit[k].resize(NradAzmBins, 0);
      azmRefSplitCount[k].resize(NradAzmBins, 0);
    }
    std::vector<double> azmRefSplitDiff(NradAzmBins, 0);
    int splitInd = 0;
    std::map<int, int> normInd;
    azmIndReference.clear();
    for (int ir=0; ir<NradAzmBins; ir++) {
      norm = 0;
      for (auto & nItr : normInd) {
        nItr.second = 0;
      }
      for (auto const & sRitr : scanReferences) {
        for (auto const & pItr : sRitr.second) {
          if (normInd.find(pItr.first) == normInd.end()) {
            normInd[pItr.first] = 0;
            azmIndReference[pItr.first].resize(NradAzmBins, 0);
          }
          if (pItr.second.scale) {
            if (pItr.second.azmRef[ir] != NANVAL) {
              azmReference[ir] += pItr.second.azmRef[ir];
              azmIndReference[pItr.first][ir] += pItr.second.azmRef[ir];
              norm += 1;
              normInd[pItr.first] += 1;
              azmRefSplit[splitInd][ir] += pItr.second.azmRef[ir];
              azmRefSplitCount[splitInd][ir] += 1;
            }
          }
          splitInd  = (splitInd + 1)%(int)azmRefSplit.size();
        }
      }
      //cout<<"before: "<<azmReference[ir]<<"  "<<norm<<endl;
      //if (azmRefSplitCount[0][ir]) {
      if (norm) {
        azmReference[ir]    /= norm;
        //azmReference[ir] = azmRefSplit[0][ir]/azmRefSplitCount[0][ir];
      }
      else {
        azmReference[ir]    = NANVAL;
      }
      for (auto const & nItr : normInd) {
        if (nItr.second > 0) {
          azmIndReference[nItr.first][ir] /= nItr.second;
        }
        else {
          azmIndReference[nItr.first][ir] = NANVAL;
        }
      }
    }

    // Gaussian smoothing 
    if (mergeGaussSmoothRef) {
      for (int iq=0; iq<NradAzmBins; iq++) {
        if (azmReference[iq] != NANVAL) {
          azmReference[iq] *= sMsAzmNorm[iq];
        }
      }

      azmReference = imgProc::gaussianSmooth1d(
          azmReference, 
          mergeGSmoothSTD, 
          7*mergeGSmoothSTD);
      
      for (int iq=0; iq<NradAzmBins; iq++) {
        if (azmReference[iq] != NANVAL) {
          azmReference[iq] /= sMsAzmNorm[iq];
        }
      }
    }

    for (int ir=50; ir<NradAzmBins; ir++) {
        if (azmRefSplitCount[0][ir] && azmRefSplitCount[1][ir]) {
          azmRefSplitDiff[ir] = azmRefSplit[0][ir]/azmRefSplitCount[0][ir] - azmRefSplit[1][ir]/azmRefSplitCount[1][ir];
        }
      }
      plt->print1d(azmRefSplitDiff, "testingRefDiff");

    std::vector<double> azmRefCorr(0);
    if (refCorrection.compare("NULL") != 0) {
      azmRefCorr.resize(NradAzmBins);
      save::importDat<double>(azmRefCorr, refCorrection);
      for (int ir=0; ir<NradAzmBins; ir++) {
        if (azmReference[ir] != NANVAL) {
          azmReference[ir] += azmRefCorr[ir];
        }
      }
    }

  }


  /////  Merging time dependent diffraction  /////

    vector<double> sums(scanAzmAvg.size(), 0);
    vector<double> sumsC(scanAzmAvg.size(), 0);
  if (tdOnly || !refOnly) {
    // Initialize variables
    legendres.clear();
    legendresMs.clear();
    legendres.resize(Nlegendres);
    legendresMs.resize(Nlegendres);
    for (int ilg=0; ilg<Nlegendres; ilg++) {
      legendres[ilg].resize(stagePosInds.size());
      legendresMs[ilg].resize(stagePosInds.size());
      for (uint it=0; it<stagePosInds.size(); it++) {
        legendres[ilg][it].resize(NradLegBins, 0);
        std::fill(legendres[ilg][it].begin(), legendres[ilg][it].end(), 0);
        legendresMs[ilg][it].resize(NradLegBins, 0);
        std::fill(legendresMs[ilg][it].begin(), legendresMs[ilg][it].end(), 0);
      }
    }

    azimuthalAvg.resize(stagePosInds.size());
    azimuthalsMs.resize(stagePosInds.size());
    for (uint it=0; it<stagePosInds.size(); it++) {
      azimuthalAvg[it].resize(NradAzmBins, 0);
      std::fill(azimuthalAvg[it].begin(), azimuthalAvg[it].end(), 0);
      azimuthalsMs[it].resize(NradAzmBins, 0);
      std::fill(azimuthalsMs[it].begin(), azimuthalsMs[it].end(), 0);
    }

    //  Merging loop
    if (verbose)
      std::cout << "\tMerging time dependent images.\n";
    for (int it=0; it<NtimeSteps; it++) {

      ///  Merging legendres  ///
      for (int ilg=0; ilg<Nlegendres; ilg++) {
        for (int ir=0; ir<NradLegBins; ir++) {
          rInd = ilg*NradLegBins + ir;
          norm = 0;
          //for (auto sLiter : scanLgndrs) {
            //cout<<"filling: "<<ilg<<"  "<<it<<"  "<<ir<<"  "<<sLiter.first<<"  "<<sLiter.second[it][rInd]<<"  "<<sLiter.second[it][rInd] -     runMeans[it][rInd]<<"  "<< stdScale*runStdev[it][rInd]<<endl;
            //if (sLiter.second[it][rInd] != NANVAL) {
            //  legendres[ilg][it][ir] += sLiter.second[it][rInd];
            //  legNorm += scanScale[sLiter.first][it];
            //}
          //}
          //if (norm) {
          //  legendres[ilg][it][ir] /= norm;
          //}
          //else {
          //  legendres[ilg][it][ir] = 0;
          //}
          //cout<<"filled for lg: "<<ilg<<endl;
        }
      }

      ///  Merging azimuthal averages  ///
      int tInd;
      for (int iazm=0; iazm<NradAzmBins; iazm++) {
        norm = 0;
        tInd = 0;
        for (auto const & sAzml : scanAzmAvg) {
          if ((sAzml.second[it][iazm] != NANVAL) 
              && (scanScale[sAzml.first][it] != 0)) {
            azimuthalAvg[it][iazm]  += sAzml.second[it][iazm];
            norm += 1; //scanScale[sAzml.first][it];
            /*
            if ((it==27) && (iazm > 300) && (iazm < 450)) {
              sums[tInd] += sAzml.second[it][iazm];
              sumsC[tInd] += 1;
            }
            */
          }
          tInd++;
        }
        if (norm) {
          azimuthalAvg[it][iazm]  /= norm;
        }
        else {
          azimuthalAvg[it][iazm]  = NANVAL;
        }
      }
    }
    /*
    for (int i=0; i<sums.size(); i++) {
      sums[iazm] /= sumsC[iazm];
    }
    std::vector<int> inds(sums.size());
    std::iota(inds.begin(), inds.end(), 0);
    */

  }
}


void mergeClass::subtractT0() {

  didSubtractT0 = true;

  for (int ilg=0; ilg<Nlegendres; ilg++) {
    // Subtract reference or time early from legendres
    if (!subtractReference) {
      // Raw data
      std::fill(legReference[ilg].begin(), legReference[ilg].end(), 0.0);
      std::vector<int> count(NradLegBins, 0);
      //std::fill(count.begin(), count.end(), 0.0);
      for (int ir=0; ir<NradLegBins; ir++) {
        for (uint tm=0; tm<legendres[ilg].size(); tm++) {
          if (legendres[ilg][tm][ir] != NANVAL) {
            legReference[ilg][ir] += legendres[ilg][tm][ir];
            count[ir]++;
          }
        }
        if (count[ir]) {
          legReference[ilg][ir] /= count[ir];
        }
        else {
          legReference[ilg][ir] = 0;
        }
      }
    }

    for (int ir=0; ir<NradLegBins; ir++) {
      for (uint tm=0; tm<legendres[ilg].size(); tm++) {
        if (legendres[ilg][tm][ir] != NANVAL) {
          legendres[ilg][tm][ir]  -= legReference[ilg][ir];
        }
      }
    }
  }

  // Subtract reference or time early from azimuthal average
  if (!subtractReference) {
    std::vector<int> count(NradAzmBins, 0);
    count.resize(NradAzmBins,0);
    std::fill(azmReference.begin(), azmReference.end(), 0);
    for (int ir=0; ir<NradAzmBins; ir++) {
      for (uint tm=0; tm<azimuthalAvg.size(); tm++) {
        if (azimuthalAvg[tm][ir] != NANVAL) {
          azmReference[ir] += azimuthalAvg[tm][ir];
          count[ir]++;
        }
      }
      if (count[ir]) {
        azmReference[ir] /= count[ir];
      }
      else {
        azmReference[ir] = NANVAL;
      }
    }
  }


  /////  Editting azimuthal average  /////
  for (int ir=0; ir<NradAzmBins; ir++) {
    for (uint tm=0; tm<stagePosInds.size(); tm++) {
      if ((azimuthalAvg[tm][ir] != NANVAL) && (azmReference[ir] != NANVAL)) {
        azimuthalAvg[tm][ir] -= azmReference[ir];
      }
      else {
        azimuthalAvg[tm][ir] = NANVAL;
      }
    }
    for (auto & aItr : azmIndReference) {
      if ((aItr.second[ir] != NANVAL) && (azmReference[ir] != NANVAL)) {
        aItr.second[ir] -= azmReference[ir];
      }
      else {
        aItr.second[ir] = NANVAL;
      }
    }

  }

}


void mergeClass::normalizeScansResults() {

  double refReadoutNoise = 0;
  int refReadoutCount = 0;
  for (int iq=readoutAzmBinStart; iq<readoutAzmBinEnd; iq++) {
    if (azmReference[iq] != NANVAL) {
      refReadoutNoise += azmReference[iq];
      refReadoutCount++;
    }
  }
  refReadoutNoise /= refReadoutCount;
  for (int iq=0; iq<NradAzmBins; iq++) {
    if (azmReference[iq] != NANVAL) {
      azmReference[iq] -= refReadoutNoise;
    }

    for (auto & sRitr : scanReferences) {
      for (auto & pItr : sRitr.second) {
        if (pItr.second.azmRef[iq] != NANVAL) {
          pItr.second.azmRef[iq] -= refReadoutNoise;
        }
      }
    }
  }


  double readoutNoise = 0;
  int readoutCount = 0;
  for (uint it=0; it<stagePosInds.size(); it++) {
    readoutNoise = 0;
    readoutCount = 0;
    for (int iq=readoutAzmBinStart; iq<readoutAzmBinEnd; iq++) {
      if (azimuthalAvg[it][iq] != NANVAL) {
        readoutNoise += azimuthalAvg[it][iq];
        readoutCount++;
      }
    }
    readoutNoise /= readoutCount;
    for (int iq=0; iq<NradAzmBins; iq++) {
      if (azimuthalAvg[it][iq] != NANVAL) {
        azimuthalAvg[it][iq] -= readoutNoise;
      }

      for (auto & sAzml : scanAzmAvg) {
        if (sAzml.second[it][iq] != NANVAL) {
          sAzml.second[it][iq] -= readoutNoise;
        }
      }
    }
  }

  // Normalize Intensity
  int imgNormBinMin = (int)(imgNormRadMin*NradAzmBins);
  int imgNormBinMax = (int)(imgNormRadMax*NradAzmBins);

  double refNorm = 0;
  int refNormCount = 0;
  for (int iq=imgNormBinMin; iq<imgNormBinMax; iq++) {
    if (azmReference[iq] != NANVAL) {
      refNorm += azmReference[iq];
      refNormCount++;
    }
  }
  refNorm /= refNormCount;
  for (int iq=0; iq<NradAzmBins; iq++) {
    if (azmReference[iq] != NANVAL) {
      azmReference[iq] /= refNorm;
    }

    for (auto & sRitr : scanReferences) {
      for (auto & pItr : sRitr.second) {
        if (pItr.second.azmRef[iq] != NANVAL) {
          pItr.second.azmRef[iq] /= refNorm;
        }
      }
    }
  }


  double azmNorm = 0;
  int azmNormCount = 0;
  for (uint it=0; it<stagePosInds.size(); it++) {
    azmNorm = 0;
    azmNormCount = 0;
    for (int iq=imgNormBinMin; iq<imgNormBinMax; iq++) {
      if (azimuthalAvg[it][iq] != NANVAL) {
        azmNorm += azimuthalAvg[it][iq];
        azmNormCount++;
      }
    }
    azmNorm /= azmNormCount;
    if (it==2) {
      cout<<"norm range: "<<imgNormBinMin<<"  "<<imgNormBinMax<<endl;
      cout<<"norm: "<<azmNorm<<endl;
    }
    for (int iq=0; iq<NradAzmBins; iq++) {
      if (azimuthalAvg[it][iq] != NANVAL) {
        azimuthalAvg[it][iq] /= azmNorm;
      }

      for (auto & sAzml : scanAzmAvg) {
        if (sAzml.second[it][iq] != NANVAL) {
          sAzml.second[it][iq] /= azmNorm;
        }
      }
    }
  }

  return;
}



void mergeClass::sMsNormalize() {

  didSMSnormalize = true;

  /////  Normalizing legendres  /////
  for (int ilg=0; ilg<Nlegendres; ilg++) {
    for (int ir=0; ir<NradLegBins; ir++) {
      for (uint tm=0; tm<legendres[ilg].size(); tm++) {
        if (legendres[ilg][tm][ir] != NANVAL) {
          legendresMs[ilg][tm][ir] = legendres[ilg][tm][ir]*sMsLegNorm[ir];
        }
      }
    }
  }

  /////  Normalizing azimuthal average  /////
  azimuthalsMs.resize(stagePosInds.size());
  runsMsMeans.resize(stagePosInds.size());
  runsMsSEM.resize(stagePosInds.size());
  for (uint tm=0; tm<stagePosInds.size(); tm++) {
    azimuthalsMs[tm].resize(NradAzmBins, 0);
    runsMsMeans[tm].resize(NradAzmBins, 0);
    runsMsSEM[tm].resize(NradAzmBins, 0);
  }
  runsMsRefMean.resize(NradAzmBins, 0);
  runsMsRefSEM.resize(NradAzmBins, 0);


  for (int ir=0; ir<NradAzmBins; ir++) {
    for (uint tm=0; tm<stagePosInds.size(); tm++) {
      if (azimuthalAvg[tm][ir] != NANVAL) {
        azimuthalsMs[tm][ir]    = azimuthalAvg[tm][ir]*sMsAzmNorm[ir];
        runsMsMeans[tm][ir]     = runAzmMeans[tm][ir]*sMsAzmNorm[ir];
        runsMsSEM[tm][ir]       = runAzmSEM[tm][ir]*sMsAzmNorm[ir];
      }
      else {
        azimuthalsMs[tm][ir]    = NANVAL;
        runsMsMeans[tm][ir]     = NANVAL;
        runsMsSEM[tm][ir]       = NANVAL;
      }
    }
    if (runAzmRefMean[ir] != NANVAL) {
      runsMsRefMean[ir]     = runAzmRefMean[ir]*sMsAzmNorm[ir];
      runsMsRefSEM[ir]      = runAzmRefSEM[ir]*sMsAzmNorm[ir];
    }
    else {
      runsMsRefMean[ir]     = NANVAL;
      runsMsRefSEM[ir]      = NANVAL;
    }
  }

}

void mergeClass::gaussianFilterQ() {
 
  if (!didSMSnormalize) {
    std::cerr << "ERROR: cannot gaussian smooth sMs without normalizing!!!\n";
    exit(1);
  }

  int sInd;
  double scale;
  cv::Mat origLO(NradAzmBins, 1, CV_64F);
  cv::Mat smoothLO(NradAzmBins, 1, CV_64F);
  for (uint itm=0; itm<azimuthalsMs.size(); itm++) {
    sInd = 0;
    while (azimuthalsMs[itm][sInd] == NANVAL) {
      sInd++;
    }
    
    /*
    origLO.resize(NradAzmBins-sInd, 1);
    smoothLO.resize(NradAzmBins-sInd, 1);
    for (int iq=sInd; iq<NradAzmBins; iq++) {
      origLO.at<double>(iq-sInd,0) = azimuthalAvg[itm][iq];
    }
    //cv::GaussianBlur(origLO, smoothLO, cvSize(51,1), 7, 0);
    cv::medianBlur(origLO, smoothLO, 11); 
    for (int iq=sInd; iq<NradAzmBins; iq++) {
      azimuthalAvg[itm][iq] = smoothLO.at<double>(iq-sInd,0);
    }
    */

    int fsize = 50;
    std::vector<double> tmp(NradAzmBins, 0);
    for (int iq=sInd; iq<NradAzmBins; iq++) {
      double norm = 0;
      for (int iqq=-1*fsize/2; iqq<=fsize/2; iqq++) {
        if ((iq + iqq >= sInd) && (iq + iqq < NradAzmBins)) {
          scale   = std::exp(-0.5*std::pow(iqq, 2)/49);
          norm    += scale;
          tmp[iq] += scale*azimuthalsMs[itm][iq+iqq];
        }
      }
      tmp[iq] /= norm;
    }
    for (int iq=sInd; iq<NradAzmBins; iq++) {
      azimuthalsMs[itm][iq] = tmp[iq];
    }


  }
}


void mergeClass::smearTimeGaussian() {
  if (verbose) {
    std::cout << "Begin smearing legendres in time\n";
  }

  ///  Plotting time dependend legendre signal  ///
  /*
  std::vector<PLOToptions> opts(3);
  std::vector<string> vals(3);
  opts[0] = yLabel;   vals[0] = "Time [ps]";
  opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
  opts[2] = draw;     vals[2] = "COLZ";
  std::string cbarLegName = curDate + "_" + to_string(curScan) + "_Leg" + to_string(ilg);
  if (cbar.find(cbarLegName) != cbar.end()) {
    opts.push_back(minimum);
    vals.push_back(to_string(-1*cbar[cbarLegName]));
    opts.push_back(maximum);
    vals.push_back(to_string(cbar[cbarLegName]));
  }
  */

  // Find smallest difference
  double minDist = 1000;
  for (int it=0; it<(int)timeDelays.size()-1; it++) {
    if (timeDelays[it+1] - timeDelays[it] < minDist) {
      minDist = timeDelays[it+1] - timeDelays[it];
    }
  }
  minDist /= 5.;

  double totalTime = timeDelays[timeDelays.size()-1] - timeDelays[0];
  int NfineTimeDelays = totalTime/minDist;
  std::vector<double> fineTimeDelays(NfineTimeDelays);
  std::vector<int> timeMapping(NfineTimeDelays);
  int timeInd = 0;
  for (int tm=0; tm<NfineTimeDelays; tm++) {
    fineTimeDelays[tm] = timeDelays[0] + tm*minDist;
    if (timeInd + 1 != (int)azimuthalAvg.size()) {
      while (fabs(timeDelays[timeInd] - fineTimeDelays[tm]) >
          fabs(timeDelays[timeInd+1] - fineTimeDelays[tm])) {
        timeInd++;
      }
    }
    timeMapping[tm] = timeInd;
  }

  /////  Smearing  /////
  double sum, weight;
  smearedAzmAvg.resize(NfineTimeDelays);
  smearedAzmsMs.resize(NfineTimeDelays);
  for (int it=0; it<NfineTimeDelays; it++) {
    smearedAzmAvg[it].resize(NradAzmBins, 0);
    smearedAzmsMs[it].resize(NradAzmBins, 0);
  }
  for (int iq=0; iq<NradAzmBins; iq++) {
    for (int tm=0; tm<NfineTimeDelays; tm++) {
      sum = 0;
      smearedAzmAvg[tm][iq] = 0;
      for (int tmm=(int)(-1*smearTimeBinWindow/2); 
          tmm<=smearTimeBinWindow/2; tmm++) {
        if (tm + tmm < 0) continue;
        if (tm + tmm >= NfineTimeDelays) break;
        weight = std::exp(-1*std::pow(fineTimeDelays[tm] 
                     - fineTimeDelays[tm+tmm], 2)
                     /(2*std::pow(timeSmearSTD, 2)));
        smearedAzmAvg[tm][iq] += azimuthalAvg[timeMapping[tm+tmm]][iq]*weight;
        sum += weight;
      }
      smearedAzmAvg[tm][iq] /= sum;
      smearedAzmsMs[tm][iq] = smearedAzmAvg[tm][iq]*sMsAzmNorm[iq];
    }
  }



  /*
  smearedImg.resize(Nlegendres);
  std::vector<double> unPumped(legendres[0][0].size());

  for (int ilg=0; ilg<Nlegendres; ilg++) {

    double sum = 0;
    smearedImg[ilg].resize(stagePosInds.size());
    for (uint ir=0; ir<stagePosInds.size(); ir++) {
      smearedImg[ilg][ir].resize(NradLegBins, 0);
    }
    for (int ir=0; ir<NradLegBins; ir++) {
      for (int tm=0; tm<(int)stagePosInds.size(); tm++) {
        sum = 0;
        for (int tmm=0; tmm<(int)stagePosInds.size(); tmm++) {
          smearedImg[ilg][tm][ir] += legendres[ilg][tmm][ir]
                                    *std::exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(smearSTD, 2)));
          sum += exp(-1*std::pow(timeDelays[tm] - timeDelays[tmm], 2)/(2*std::pow(smearSTD, 2)));
        }
        smearedImg[ilg][tm][ir] /= sum;
      }
    }

    if (ilg == 0) {
      save::saveDat<double>(unPumped, "./results/data_unPumpedDiffractionL0Smeared["
          + to_string(NradLegBins) + "].dat");
    }

    cout<<"last entry"<<endl;
    cout<<"plotting"<<endl;

      vals[2] = "-8";
      vals[3] = "8";
      plt.printXY(img[i], timeDelays, "plotRun" + to_string(curRun) + "_" + curDate + "_" + curScan 
            + "_Leg" + to_string(i), oppts, vaals);

    ///  Plotting time dependend legendre signal  ///
    std::vector<PLOToptions> opts(3);
    std::vector<string> vals(3);
    opts[0] = yLabel;   vals[0] = "Time [ps]";
    opts[1] = xLabel;   vals[1] = "Scattering Q [arb units]";
    opts[2] = draw;     vals[2] = "COLZ";


    plt.printRC(smearedImg[ilg], timeDelays, 0, maxQ, "smeared_" + curDate + "_"
          + to_string(curScan) + "_Leg" + to_string(ilg), opts, vals);
    */

  smearedTime = true;

}

void mergeClass::smearTimeFFT() {

  /////  Rebinning into smallest bins  /////

  // Find smallest difference
  double minDist = 1000;
  for (int it=0; it<(int)timeDelays.size()-1; it++) {
    if (timeDelays[it+1] - timeDelays[it] < minDist) {
      minDist = timeDelays[it+1] - timeDelays[it];
    }
  }

  double totalTime = timeDelays[timeDelays.size()-1] - timeDelays[0];
  int NfineTimeDelays = totalTime/minDist + 1;
  std::vector<double> fineTimeDelays(NfineTimeDelays);
  for (int it=0; it<NfineTimeDelays; it++) {
    fineTimeDelays[it] = it*minDist;
  }

  ///// Filter Parameters  /////
  int filtFFToutSize = (int)(NfineTimeDelays/2 + 1);
  double* tSpace = (double*) fftw_malloc(sizeof(double)*(NfineTimeDelays));
  fftw_complex* fSpace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*filtFFToutSize);
  fftw_plan filtFFTf = fftw_plan_dft_r2c_1d(NfineTimeDelays, tSpace, fSpace, FFTW_MEASURE);
  fftw_plan filtFFTb = fftw_plan_dft_c2r_1d(NfineTimeDelays, fSpace, tSpace, FFTW_MEASURE);

  std::string filterName =
      "/reg/neh/home/khegazy/analysis/filters/" + timeFilterType
      + "Filter_Bins-" + to_string(filtFFToutSize)
      + "_WnHigh-"+ to_string(timeWnHigh) + ".dat";
  if (!tools::fileExists(filterName)) {
    cout << "INFO: Making new filter\n";
    system(("python /reg/neh/home/khegazy/analysis/filters/makeFilters.py --Nbins "
          + to_string(filtFFToutSize)
          + " --Ftype " + timeFilterType
          + " --Order " + to_string(timeFiltOrder)
          + " --WnHigh " + to_string(timeWnHigh)).c_str());
  }

  std::vector<double> bandPassFilter(filtFFToutSize);
  save::importDat<double>(bandPassFilter, filterName);


  /////  Fourier Transform  /////
  smearedAzmAvg.resize(NfineTimeDelays);
  smearedAzmsMs.resize(NfineTimeDelays);
  for (int it=0; it<NfineTimeDelays; it++) {
    smearedAzmAvg[it].resize(NradAzmBins, 0);
    smearedAzmsMs[it].resize(NradAzmBins, 0);
  }
  for (int iq=0; iq<NradAzmBins; iq++) {
    int timeInd = 0;
    for (int it=0; it<NfineTimeDelays; it++) {
      if (timeInd + 1 != (int)azimuthalAvg.size()) {
        while (fabs(timeDelays[timeInd] - fineTimeDelays[it]) >
            fabs(timeDelays[timeInd+1] - fineTimeDelays[it])) {
          timeInd++;
        }
      }
      tSpace[it] = azimuthalAvg[timeInd][iq];
    }

    fftw_execute(filtFFTf);

    for (int i=0; i<filtFFToutSize; i++) {
      fSpace[i][0] *= bandPassFilter[i]/sqrt(NfineTimeDelays);
      fSpace[i][1] *= bandPassFilter[i]/sqrt(NfineTimeDelays);
    }

    fftw_execute(filtFFTb);

    for (int it=0; it<NfineTimeDelays; it++) {
      smearedAzmAvg[it][iq] = tSpace[it]/sqrt(NfineTimeDelays);
      smearedAzmsMs[it][iq] = smearedAzmAvg[it][iq]*sMsAzmNorm[iq];
    }
  }

  smearedTime = true;

  fftw_destroy_plan(filtFFTf);
  fftw_destroy_plan(filtFFTb);
  fftw_free(tSpace);
  fftw_free(fSpace);
}



void mergeClass::makePairCorrs() {

  didPairCorrSTD = true;

  if (!didSMSnormalize) {
    sMsNormalize();
  }

  scanAzmPCorr.clear();
  int NtimeSteps = stagePosInds.size();
  std::vector< std::vector<double> > outAzm(NtimeSteps);
  std::vector< std::vector<double> > inpPC(NtimeSteps);
  for (int i=0; i<NtimeSteps; i++) {
    outAzm[i].resize(NradAzmBins, 0);
    inpPC[i].resize(maxRbins, 0);
  }
  for (auto const & sAzml : scanAzmAvg) {
    scanAzmPCorr[sAzml.first].resize(NtimeSteps);
    for (int it=0; it<NtimeSteps; it++) {
      for (int iazm=0; iazm<NradAzmBins; iazm++) {
        if (sAzml.second[it][iazm] != NANVAL) {
          outAzm[it][iazm] = sAzml.second[it][iazm]*sMsAzmNorm[iazm];
        }
        else {
          outAzm[it][iazm] = 0;
        }
      }
    }


    std::string fName = "singleScanSMS_Bins["
        + to_string(NtimeSteps) + ","
        + to_string(NradAzmBins) + "].dat";
    save::saveDat<double>(outAzm, "./results/" + fName);

    system(("./../timeDepStudies/pairCorr.exe " + runName 
        + " -Idir ./results/ -Fname "
        + fName).c_str()); 

    std::string pfName = "data-20180627_1551_pairCorrOdd["
        + to_string(NtimeSteps) + ","
        + to_string(maxRbins) + "].dat";
    save::importDat<double>(inpPC, "./results/" + pfName);

    scanAzmPCorr[sAzml.first].resize(NtimeSteps);
    for (int it=0; it<NtimeSteps; it++) {
      scanAzmPCorr[sAzml.first][it].resize(maxRbins);
      for (int ir=0; ir<maxRbins; ir++) {
        scanAzmPCorr[sAzml.first][it][ir] = inpPC[it][ir];
      }
    }

    system(("rm ./results/" + fName).c_str());
    system(("rm ./results/" + pfName).c_str());
  }
}


void mergeClass::basicLessThanCut(std::string paramName, double cut) {

  for (auto const & timeItr : labTimeParams) {
    if (labTimeParams[timeItr.first][paramName] < cut) {
      int scan      = labTimeMap[timeItr.first].first;
      int stagePos  = labTimeMap[timeItr.first].second;

      if (scanReferences[scan].find(stagePos) != scanReferences[scan].end()) {
        scanReferences[scan][stagePos].scale = 0;
      }
      else if (stagePosInds.find(stagePos) != stagePosInds.end()) {
        scanScale[scan][stagePosInds[stagePos]] = 0;
      }
      else {
        std::cerr << "ERROR: Cannot find image by scan / stagePos: "
            << scan << " / " << stagePos<<endl;
        exit(1);
      }
    }
  }

  return;
}


void mergeClass::basicGreaterThanCut(std::string paramName, double cut) {

  for (auto const & timeItr : labTimeParams) {
    if (labTimeParams[timeItr.first][paramName] > cut) {
      int scan      = labTimeMap[timeItr.first].first;
      int stagePos  = labTimeMap[timeItr.first].second;

      /*
      if (stagePos == 1542450) {
        cout<<"Cut 1542450: "<<scan<<"       "<<labTimeParams[timeItr.first][paramName]<<endl;
      }
      if (stagePos == 1542650) {
        cout<<"Cut 1542650: "<<scan<<"       "<<labTimeParams[timeItr.first][paramName]<<endl;
      }
      */

      if (scanReferences[scan].find(stagePos) != scanReferences[scan].end()) {
        scanReferences[scan][stagePos].scale = 0;
      }
      if (stagePosInds.find(stagePos) != stagePosInds.end()) {
        scanScale[scan][stagePosInds[stagePos]] = 0;
      }
    }
  }

  return;
}


void mergeClass::stdParamCut(std::string paramName, double cut) {

  int scan;
  int stagePos;

  /////  Calculating the mean  /////
  double mean   = 0;
  double count  = 0;
  for (auto const & timeItr : labTimeParams) {
    scan      = labTimeMap[timeItr.first].first;
    stagePos  = labTimeMap[timeItr.first].second;
    if (scanReferences[scan].find(stagePos) != scanReferences[scan].end()) {
      if (scanReferences[scan][stagePos].scale > 0) {
        mean += labTimeParams[timeItr.first][paramName];
        count += 1;
      }
    }
    else if (stagePosInds.find(stagePos) != stagePosInds.end()) {
      if (scanScale[scan][stagePosInds[stagePos]] > 0) {
        mean += labTimeParams[timeItr.first][paramName];
        count += 1;
      }
    }
    else {
      std::cerr << "ERROR: Cannot find image by scan / stagePos: "
          << scan << " / " << stagePos<<endl;
      exit(1);
    }
  }
  mean /= count;

  /////  Calculating the std  /////
  double stdVal   = 0;
  count  = 0;
  for (auto const & timeItr : labTimeParams) {
    scan      = labTimeMap[timeItr.first].first;
    stagePos  = labTimeMap[timeItr.first].second;
    if (scanReferences[scan].find(stagePos) != scanReferences[scan].end()) {
      if (scanReferences[scan][stagePos].scale > 0) {
        stdVal += std::pow(mean-labTimeParams[timeItr.first][paramName], 2);
        count += 1;
      }
    }
    else if (stagePosInds.find(stagePos) != stagePosInds.end()) {
      if (scanScale[scan][stagePosInds[stagePos]] > 0) {
        stdVal += std::pow(mean-labTimeParams[timeItr.first][paramName], 2);
        count += 1;
      }
    }
    else {
      std::cerr << "ERROR: Cannot find image by scan / stagePos: "
          << scan << " / " << stagePos<<endl;
      exit(1);
    }
  }
  stdVal = std::sqrt(stdVal/count);

  /////  Make cut  /////
  for (auto const & timeItr : labTimeParams) {
    scan      = labTimeMap[timeItr.first].first;
    stagePos  = labTimeMap[timeItr.first].second;
    if (fabs(labTimeParams[timeItr.first][paramName] - mean)/stdVal > cut) {
      int scan      = labTimeMap[timeItr.first].first;
      int stagePos  = labTimeMap[timeItr.first].second;

      /*
      if (stagePos == 1542450) {
        cout<<"Cut imgNorm 1542450: "<<scan<<"       "<<labTimeParams[timeItr.first][paramName]<<endl;
      }
      if (stagePos == 1542650) {
        cout<<"Cut imgNorm 1542650: "<<scan<<"       "<<labTimeParams[timeItr.first][paramName]<<endl;
      }
      */


      if (scanReferences[scan].find(stagePos) != scanReferences[scan].end()) {
        scanReferences[scan][stagePos].scale = 0;
      }
      if (stagePosInds.find(stagePos) != stagePosInds.end()) {
        scanScale[scan][stagePosInds[stagePos]] = 0;
      }
    }
  }

  return;
}



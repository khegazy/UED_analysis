#include "mergeClass.h"


mergeClass::mergeClass(std::string runName_inp) : parameterClass(runName_inp) {

  runName = runName_inp;
  initializeVariables();
}


mergeClass::~mergeClass() {
  delete plt;
  delete[] timeDelays;
}


void mergeClass::initializeVariables() {

  plt = new PLOTclass("mergeCanv");
  curScan  = -1;
  NlegBins = NradLegBins*Nlegendres;
  runInd   = "";
  curDate  = "";
  curRun   = "";
  _compareReference = false;

  timeDelays = NULL;

  /////  Simulation results  /////
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
  for (int ir=0; ir<NradLegBins; ir++) {
    sMsLegNorm[ir] = Qleg[ir]/(atmLegDiff[ir]);
  }

  sMsAzmNorm.resize(NradAzmBins, 0.0);
  for (int ir=0; ir<NradAzmBins; ir++) {
    sMsAzmNorm[ir] = Qazm[ir]/(atmAzmDiff[ir]);
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


void mergeClass::addLabTimeParameter(int timeStamp,
                          int scan, int64_t stagePos, 
                          double imgNorm) {

  labTimeParams[timeStamp].scan      = scan;
  labTimeParams[timeStamp].stagePos  = stagePos;
  labTimeParams[timeStamp].imgNorm   = imgNorm;
}


void mergeClass::addReference(int scan, int64_t stagePos,
                          std::vector<double>* azmAvg, 
                          std::vector<double>* legCoeffs,
                          double imgNorm) {

  if (!normalizeImgs) imgNorm = 1;
  scanReferences[scan][stagePos].imgNorm = imgNorm;
  scanReferences[scan][stagePos].azmRef.resize(NradAzmBins, 0);
  for (int i=0; i<NradAzmBins; i++) {
    if ((*azmAvg)[i] != NANVAL) {
      scanReferences[scan][stagePos].azmRef[i] = (*azmAvg)[i]/imgNorm;
    }
    else {
      scanReferences[scan][stagePos].azmRef[i] = NANVAL;
    }
  }

  scanReferences[scan][stagePos].legRef.resize(NlegBins, 0);
  for (int i=0; i<NlegBins; i++) {
    if ((*legCoeffs)[i] != NANVAL) {
      scanReferences[scan][stagePos].legRef[i] = (*legCoeffs)[i]/imgNorm;
    }
    else {
      scanReferences[scan][stagePos].legRef[i] = NANVAL;
    }
  }
}


void mergeClass::addEntry(int scan, int64_t stagePos, 
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
    scanCounts[scan] = emptyC;

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
        auto sCitr = scanCounts.begin();
        auto sLitr = scanLgndrs.begin();
        auto sAzml = scanAzmAvg.begin();
        const std::vector<double> emptyLegVec((*legCoeffs).size(), 0.0);
        const std::vector<double> emptyAzmVec((*azmAvg).size(), 0.0);
        while (sCitr != scanCounts.end()) {
          sCitr->second.insert(sCitr->second.begin()+ind, 0);
          sLitr->second.insert(sLitr->second.begin()+ind, emptyLegVec);
          sAzml->second.insert(sAzml->second.begin()+ind, emptyAzmVec);
          sCitr++; sLitr++; sAzml++;
        }
      }
      itr.second = ind;
      ind++;
    }
  }

  if (false && verbose) {
    std::cout << "ind/sizes: " << pInd << "  "
      << scanCounts[scan].size() << "  "
      << scanLgndrs[scan].size() << "  "
      << scanLgndrs[scan][pInd].size() << std::endl;
  }

  ///  Add coefficients and counts to maps  ///
  if (!normalizeImgs) imgNorm = 1;
  scanCounts[scan][pInd] = imgNorm;
  for (uint i=0; i<(*legCoeffs).size(); i++) {
    scanLgndrs[scan][pInd][i] = (*legCoeffs)[i]/imgNorm;
  }
  for (uint i=0; i<(*azmAvg).size(); i++) {
    if ((*azmAvg)[i] != NANVAL) {
      scanAzmAvg[scan][pInd][i] = (*azmAvg)[i]/imgNorm;
    }
    else {
      scanAzmAvg[scan][pInd][i] = NANVAL;
    }
  }

  
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
    

void mergeClass::getMeanSTD() {

  /////  Reset all vectors  /////
  // Time dependent measurements
  if (runAzmMeans.size() != stagePosInds.size()) {
    runLegMeans.resize(stagePosInds.size());
    runAzmMeans.resize(stagePosInds.size());
    runLegSTD.resize(stagePosInds.size());
    runAzmSTD.resize(stagePosInds.size());
    for (uint i=0; i<stagePosInds.size(); i++) {
      runLegMeans[i].resize(NlegBins, 0);
      runAzmMeans[i].resize(NradAzmBins, 0);
      runLegSTD[i].resize(NlegBins, 0);
      runAzmSTD[i].resize(NradAzmBins, 0);
    }
  }
  for (uint i=0; i<stagePosInds.size(); i++) {
    std::fill(runLegMeans[i].begin(), runLegMeans[i].end(), 0);
    std::fill(runAzmMeans[i].begin(), runAzmMeans[i].end(), 0);
    std::fill(runLegSTD[i].begin(), runLegSTD[i].end(), 0);
    std::fill(runAzmSTD[i].begin(), runAzmSTD[i].end(), 0);
  }
  std::vector<double>legNorms(NlegBins, 0);
  std::vector<double>azmNorms(NradAzmBins, 0);

  // Reference measurements
  if (!runAzmRefMeans.size()) {
    runLegRefMeans.resize(NlegBins, 0);
    runAzmRefMeans.resize(NradAzmBins, 0);
    runLegRefSTD.resize(NlegBins, 0);
    runAzmRefSTD.resize(NradAzmBins, 0);
  }
  std::fill(runLegRefMeans.begin(), runLegRefMeans.end(), 0);
  std::fill(runAzmRefMeans.begin(), runAzmRefMeans.end(), 0);
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
      for (auto const &pItr : sRitr.second) {
        if (pItr.second.azmRef[ir] != NANVAL) {
          runAzmRefMeans[ir] += pItr.second.azmRef[ir];
          norm += 1;
        }
      }
    }
    if (norm) {
      runAzmRefMeans[ir] /= norm;
    }
    else {
      runAzmRefMeans[ir] = 0;
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
        if (pItr.second.azmRef[ir] != NANVAL) {
          runAzmRefSTD[ir] 
              += std::pow((pItr.second.azmRef[ir]
                      - runAzmRefMeans[ir]), 2);
          norm += 1;
        }
      }
    }
    if (norm) {
      runAzmRefSTD[ir] = std::sqrt(runAzmRefSTD[ir]/norm);
    }
    else {
      runAzmRefSTD[ir] = 0;
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
      if (scanCounts[sLitr.first][pItr.second] > 0) {
        for (int i=0; i<NlegBins; i++) {
          if (sLitr.second[pItr.second][i] != NANVAL) {
            runLegMeans[pItr.second][i] += sLitr.second[pItr.second][i];
            legNorms[i] += 1; //scanCounts[sLitr.first][pItr.second];
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
      if (scanCounts[sAitr.first][pItr.second] > 0) {
        for (int i=0; i<NradAzmBins; i++) {
          if (sAitr.second[pItr.second][i] != NANVAL) {
            runAzmMeans[pItr.second][i] += sAitr.second[pItr.second][i];
            azmNorms[i] += 1; //scanCounts[sAitr.first][pItr.second];
          }
        }
      }
    }

    for (int i=0; i<NradAzmBins; i++) {
      runAzmMeans[pItr.second][i] /= azmNorms[i];
    }


    /////  Standard deviation calculation  /////

    /*
    ///  Legendres  ///
    if (verbose && (k == 0) && (pItr.second == 0)) 
      std::cout << "\tCalculating legendre std.\n";
    for (auto const & sLitr : scanLgndrs) {
      if (scanCounts[sLitr.first][pItr.second]) {
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
      if (scanCounts[sAitr.first][pItr.second]) {
        for (int i=0; i<NradAzmBins; i++) {
          if (sAitr.second[pItr.second][i] != NANVAL) {
            runAzmSTD[pItr.second][i]
                += std::pow((sAitr.second[pItr.second][i]
                        - runAzmMeans[pItr.second][i]), 2);
          }
        }
      }
    }

    for (int i=0; i<NradAzmBins; i++) {
      runAzmSTD[pItr.second][i] = std::sqrt(runAzmSTD[pItr.second][i]/azmNorms[i]);
    }
  }

  return;
}


void mergeClass::removeLowPolynomials() {

  if (!_compareReference) {
    mergeScans(true, false);
  }
  getMeanSTD();

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
      if (pItr.second.imgNorm) {
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
            X(mInd,0) = std::pow(maxQazm-ir*delta, 6);
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

        Eigen::MatrixXd weights = tools::normalEquation(X, Y);

        for (int ir=0; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            pItr.second.azmRef[ir] -= weights(0)*std::pow(maxQazm-ir*delta, 6);
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
      if (scanCounts[sAzml.first][it] != 0) {
        Nnans = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if (sAzml.second[it][ir] == NANVAL) Nnans += 1;
        }
        X.resize(size-Nnans, Npoly);
        Y.resize(size-Nnans, 1);
        w.resize(size-Nnans);
        mInd = 0;
        for (int ir=NbinsSkip; ir<NradAzmBins; ir++) {
          if (sAzml.second[it][ir] != NANVAL) {
            X(mInd,0) = std::pow(maxQazm-ir*delta, 6);
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

        Eigen::MatrixXd weights = tools::normalEquation(X, Y);

        std::fill(fit.begin(), fit.end(), 0);
        for (int ir=0; ir<NradAzmBins; ir++) {
          if (sAzml.second[it][ir]!= NANVAL) {
            if (_compareReference) {
              plotme1[ir] = sAzml.second[it][ir] - compReference[ir];
            }
            else {
              plotme1[ir] = sAzml.second[it][ir] - azmReference[ir];
            }
            sAzml.second[it][ir] -= weights(0)*std::pow(maxQazm-ir*delta, 6);
            fit[ir] += weights(0)*std::pow(maxQazm-ir*delta, 6);
            if (_compareReference) {
              plotme2[ir] = sAzml.second[it][ir] - compReference[ir];
            }
            else {
              plotme2[ir] = sAzml.second[it][ir] - azmReference[ir];
            }
          }
          else{
            plotme1[ir] = 0;
            plotme2[ir] = 0;
          }
        }
        /*
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
        */
      }
    }
  }

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
  for (int k=0; k<3; k++) {

    /////  Calculate Mean and STD  /////
    getMeanSTD();

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
          pItr.second.imgNorm = 0;
        }
      }
    }
    */

    // Azimuthal Reference
    if (verbose && (k==0)) 
      std::cout << "\tMaking cuts on azimuthal.\n";
    for (auto& sRitr : scanReferences) {
      for (auto& pItr : sRitr.second) {
        int imageNoise = 0;
        for (int ir=0; ir<NradAzmBins; ir++) {
          if (pItr.second.azmRef[ir] != NANVAL) {
            if ((fabs(pItr.second.azmRef[ir] - runAzmRefMeans[ir])
                > mergeImageSTDScale*runAzmRefSTD[ir])
                || ((pItr.second.azmRef[ir] == NANVAL) && (ir >= azmNANend))) {
              imageNoise++;
            }      
            if (fabs(pItr.second.azmRef[ir] - runAzmRefMeans[ir])
                > mergeSTDscale*runAzmRefSTD[ir]) {
              pItr.second.azmRef[ir] = NANVAL;
            }      
          }
        }
        if (imageNoise > azmImageNoiseCut) {
          pItr.second.imgNorm = 0;
        }
      }
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
          scanCounts[sLitr.first][pItr.second] = 0;
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
          scanCounts[sAitr.first][pItr.second] = 0;
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
      << scanCounts.begin()->second.size() << "  "
      << runLegMeans.size() << "  " << runLegSTD.size() << std::endl;
  }
  /*
  ///  Remove time bins  ///
  for (int i=((int)removePos.size())-1; i>=0; i--) {
    int pInd = stagePosInds[removePos[i]];
    for (auto& sLitr : scanLgndrs) {
      sLitr.second.erase(sLitr.second.begin() + pInd);
    }
    for (auto& sCitr : scanCounts) {
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
      << scanCounts.begin()->second.size() << "  "
      << runMeans.size() << "  " << runStdev.size() << std::endl;
  }
  */



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
  if (timeDelays) delete[] timeDelays;
  timeDelays = new double[NtimeSteps + 1];
  for (auto pItr : stagePosInds) {
    timeDelays[tInd] = 2*pItr.first/(3e3);
    tInd++;
  }
  timeDelays[0] = 2*timeDelays[1] - timeDelays[2];
  for (int i=NtimeSteps; i>=0; i--) {
    timeDelays[i] -= timeZero + timeDelays[0];
  }

  save::saveDat<double>(timeDelays, NtimeSteps + 1,
      "./results/timeDelays["
      + to_string(stagePosInds.size() + 1) + "].dat");


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
            if (pItr.second.imgNorm) {
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
    for (int ir=0; ir<NradAzmBins; ir++) {
      norm = 0;
      for (auto const & sRitr : scanReferences) {
        for (auto const & pItr : sRitr.second) {
          if (pItr.second.imgNorm) {
            if (pItr.second.azmRef[ir] != NANVAL) {
              azmReference[ir]    += pItr.second.azmRef[ir];
              norm += 1;
            }
          }
        }
      }
      if (norm) {
        azmReference[ir]    /= norm;
      }
      else {
        azmReference[ir]    = 0;
      }
    }
  }


  /////  Merging time dependent diffraction  /////

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
      cout<<"looking at time "<<it<<endl;

      ///  Merging legendres  ///
      for (int ilg=0; ilg<Nlegendres; ilg++) {
        for (int ir=0; ir<NradLegBins; ir++) {
          rInd = ilg*NradLegBins + ir;
          norm = 0;
          //for (auto sLiter : scanLgndrs) {
            //cout<<"filling: "<<ilg<<"  "<<it<<"  "<<ir<<"  "<<sLiter.first<<"  "<<sLiter.second[it][rInd]<<"  "<<sLiter.second[it][rInd] -     runMeans[it][rInd]<<"  "<< stdScale*runStdev[it][rInd]<<endl;
            //if (sLiter.second[it][rInd] != NANVAL) {
            //  legendres[ilg][it][ir] += sLiter.second[it][rInd];
            //  legNorm += scanCounts[sLiter.first][it];
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
      for (int iazm=0; iazm<NradAzmBins; iazm++) {
        norm = 0;
        for (auto const & sAzml : scanAzmAvg) {
          if ((sAzml.second[it][iazm] != NANVAL) 
              && (scanCounts[sAzml.first][it] != 0)) {
            azimuthalAvg[it][iazm]  += sAzml.second[it][iazm];
            norm += 1; //scanCounts[sAzml.first][it];
          }
        }
        if (norm) {
          azimuthalAvg[it][iazm]  /= norm;
        }
        else {
          azimuthalAvg[it][iazm]  = NANVAL;
        }
      }
    }
  }
}


void mergeClass::subtractT0() {

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
        azmReference[ir] = 0;
      }
    }
  }


  /////  Editting azimuthal average  /////
  for (int ir=0; ir<NradAzmBins; ir++) {
    for (uint tm=0; tm<stagePosInds.size(); tm++) {
      if (azimuthalAvg[tm][ir] != NANVAL) {
        azimuthalAvg[tm][ir] -= azmReference[ir];
      }
    }
  }

}



void mergeClass::normalize() {

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
  for (int ir=0; ir<NradAzmBins; ir++) {
    for (uint tm=0; tm<stagePosInds.size(); tm++) {
      if (azimuthalAvg[tm][ir] != NANVAL) {
        azimuthalsMs[tm][ir] = azimuthalAvg[tm][ir]*sMsAzmNorm[ir];
        runAzmMeans[tm][ir]  *= sMsAzmNorm[ir];
        runAzmSTD[tm][ir]    *= sMsAzmNorm[ir];
      }
    }
    runAzmRefMeans[ir]  *= sMsAzmNorm[ir];
    runAzmRefSTD[ir]    *= sMsAzmNorm[ir];
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

    /*
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

  }
}

/*
void mergeClass::smearTimeFFT() {

  /////  Rebinning into smallest bins  /////

  // Find smallest difference
  double minDist = 1000;
  for (int it=0; it<NtimeBins-1; it++) {
    if (timeDelays[it+1] - timeDelays[it] < minDist) {
      minDist = timeDelays[it+1] - timeDelays[it];
    }
  }

  double totalTime = timeDelays[NtimeBins] - timeDelays[0];
  int NfineTimeBins = totalTime/minDist;
  std::vector<double> fineTimeBins(NfineTimeBins);
  for (int it=0; it<NfineTimeBins; it++) {
    fineTimeBins[it] = it*minDist;
  }

  ///// Filter Parameters  /////
  int filtFFToutSize = (int)(NfineTimeBins/2 + 1);
  double* tSpace = (double*) fftw_malloc(sizeof(double)*(NfineTimeBins));
  fftw_complex* fSpace = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*filtFFToutSize);
  fftw_plan filtFFTf = fftw_plan_dft_r2c_1d(params.NradAzmBins, tSpace, fSpace, FFTW_MEASURE);
  fftw_plan filtFFTb = fftw_plan_dft_c2r_1d(params.NradAzmBins, fSpace, tSpace, FFTW_MEASURE);

  std::string filterName =
      "/reg/neh/home/khegazy/analysis/filters/" + params.timeFilterType
      + "Filter_Bins-" + to_string(filtFFToutSize)
      + "_WnHigh-"+ to_string(params.timeWnHigh) + ".dat";
  if (!tools::fileExists(filterName)) {
    cout << "INFO: Making new filter\n";
    system(("python /reg/neh/home/khegazy/analysis/filters/makeFilters.py --Nbins "
          + to_string(filtFFToutSize)
          + " --Ftype " + params.timeFilterType
          + " --Order " + to_string(params.timeOrder)
          + " --WnHigh " + to_string(params.timeWnHigh)).c_str());
  }

  std::vector<double> bandPassFilter(filtFFToutSize);
  save::importDat<double>(bandPassFilter, filterName);


  /////  Fourier Transform  /////
  std::vector< std::vector<double> > smearedAzmAvg(NfineTimeBins);
  for (int it=0; it<NfineTimeBins; it++) 
    smearedAzmAvg[it].resize(params.NradAzmAvg,0);

  for (int iq=0; iq<NradAzmBins; iq++) {
    int timeInd = 0;
    for (int it=0; it<NfineTimeBins; it++) {
      if (timeInd != NtimeBins-1) {
        while (fabs(timeDelays[timeInd] - fineTimeBins[it]) >
            fabs(timeDelays[timeInd+1] - fineTimeBins[it])) {
          timeInd++;
        }
      }
      tSpace[it] = azimuthalAvg[timeInd][iq];
    }

    fftw_execute(filtFFTf);

    for (int ir=0; ir<filtFFToutSize; ir++) {
      fSpace[ir][0] *= bandPassFilter[ir]/sqrt(NfineTimeBins);
      fSpace[ir][1] *= bandPassFilter[ir]/sqrt(NfineTimeBins);
    }

    fftw_execute(filtFFTb);

    for (int it=0; it<NfineTimeBins; it++) {
      smearedAzmAvg[it][iq] = tSpace[it]/sqrt(NfineTimeBins);
    }
  }


  save::saveDat<double>(smearedAzmAvg, 
     "./results/data_smearedAzmAvg[" 
     + to_string(NradAzmBins) + "].dat");

  plt.printRC(smearedAzmAvg, "smearedAzmAvg");

  fftw_destroy_plan(filtFFTf);
  fftw_destroy_plan(filtFFTb);
  fftw_free(tSpace);
  fftw_free(tSpace);
}
  */
  /*
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
  */












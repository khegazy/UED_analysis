#include "alignClass.h"


alignClass::alignClass(std::string runName) : parameterClass(runName) {

  initializeVariables();
}


alignClass::~alignClass() {
  delete plt;
  delete[] timeDelays;
}


void alignClass::initializeVariables() {

  plt = new PLOTclass("alignCanv");
  legSize  = -1;
  curScan  = -1;
  NrefBins = -1;
  runInd   = "";
  curDate  = "";
  curRun   = "";

  timeDelays = NULL;

  /////  Simulation results  /////
  atmLegDiff.resize(NradLegBins, 0.0);
  molLegDiff.resize(NradLegBins, 0.0);
  simFileNameSuffix = "_Bins-" + to_string(NradLegBins) 
                      + "_Qmax-" + to_string(maxQ)
                      + "_Ieb-" + to_string(Iebeam) 
                      + "_scrnD-" + to_string(screenDist) 
                      + "_elE-" + to_string(elEnergy) + ".dat";
  save::importDat<double>(atmLegDiff, simReferenceDir + "/" 
              + molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);
  save::importDat<double>(molLegDiff, simReferenceDir + "/" 
              + molName + "_molDiffractionPatternLineOut"
              + simFileNameSuffix);

  atmAzmDiff.resize(NradAzmBins, 0.0);
  molAzmDiff.resize(NradAzmBins, 0.0);
  simFileNameSuffix = "_Bins-" + to_string(NradAzmBins) 
                      + "_Qmax-" + to_string(NradAzmBins*QperPix)
                      + "_Ieb-" + to_string(Iebeam) 
                      + "_scrnD-" + to_string(screenDist) 
                      + "_elE-" + to_string(elEnergy) + ".dat";
  save::importDat<double>(atmAzmDiff, simReferenceDir + "/" 
              + molName + "_atmDiffractionPatternLineOut"
              + simFileNameSuffix);
  save::importDat<double>(molAzmDiff, simReferenceDir + "/" 
              + molName + "_molDiffractionPatternLineOut"
              + simFileNameSuffix);


  plt->print1d(atmAzmDiff, "TestingAzmAtm");
  // Calculate normalization coefficients for modified molecular scattering
  calculateSMS();


  // Output files and location
  fileName = "alignment.txt";
  outputDir = "output/data/";

  legendres.resize(Nlegendres);
}


void alignClass::calculateSMS() {

  if (sMsLegNorm.size()) {
    std::cerr << "WARNING: Replacing sMsLegNorm values!!!\n";
  }

  sMsLegNorm.resize(NradLegBins, 0.0);
  for (int ir=0; ir<NradLegBins; ir++) {
    sMsLegNorm[ir] = (maxQ*(ir+0.5)/NradLegBins)/(1e20*atmLegDiff[ir]);
  }

  sMsAzmNorm.resize(NradAzmBins, 0.0);
  for (int ir=0; ir<NradAzmBins; ir++) {
    sMsAzmNorm[ir] = ((NradAzmBins*QperPix)*(ir+0.5)/NradAzmBins)
                        /(1e20*atmAzmDiff[ir]);
  }
}


void alignClass::compareSimulations(std::vector<std::string> radicals) {

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
      save::importDat(refLineOut[name], simReferenceDir + "/" 
                        + name + "_molDiffractionPatternLineOut"
                        + simFileNameSuffix); 
      
      if (pltVerbose) 
        delete plt->print1d(refLineOut[name], "./plots/sim_" + name + "_diffraction");

      // Normalize by atomic scattering
      for (int ir=0; ir<NradLegBins; ir++) {
        refLineOut[name][ir] *= sMsLegNorm[ir];
      } 
      // Save sMs
      save::saveDat<double>(refLineOut[name], "./plots/data/sim_" + name + "Diffraction["
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


void alignClass::addEntry(int scan, int stagePos, 
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

  legSize = (*legCoeffs).size();

  cout<<"INFO: sizes: "<<(*azmAvg).size()<<"   "<<(*legCoeffs).size()<<endl;
  ///  Get index of stage position  ///
  if (scanLgndrs.find(scan) == scanLgndrs.end()) {
    std::vector<double> emptyC(stagePosInds.size(), 0);
    scanCounts[scan] = emptyC;

    std::vector< std::vector<double> > empty(stagePosInds.size());
    for (uint ipos=0; ipos<stagePosInds.size(); ipos++) {
      empty[ipos].resize(legSize, 0);
    }
    scanLgndrs[scan] = empty;
    
    for (uint ipos=0; ipos<stagePosInds.size(); ipos++) {
      empty[ipos].resize((*azmAvg).size(), 0);
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
  scanCounts[scan][pInd] = imgNorm;
  for (uint i=0; i<(*legCoeffs).size(); i++) {
    scanLgndrs[scan][pInd][i] = (*legCoeffs)[i];
  }
  //cout<<"INFO: "<<scan<<"  "<<pInd<<"  "<<(*azmAvg).size()<<"  "<<scanAzmAvg[scan][pInd].size()<<endl;
  for (uint i=0; i<(*azmAvg).size(); i++) {
    scanAzmAvg[scan][pInd][i] = (*azmAvg)[i];
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


void alignClass::removeOutliers() {

  if (legSize == -1) {
    std::cerr << "ERROR: Must add entries before removing outliers!!!\n";
  }

  /////  Remove time bins with few images and get mean  /////
  std::vector<int32_t> removePos;
  runMeans.clear(); runStdev.clear();
  runMeans.resize(stagePosInds.size());
  runStdev.resize(stagePosInds.size());
  for (uint i=0; i<runMeans.size(); i++) {
    runMeans[i].resize(legSize, 0);
    runStdev[i].resize(legSize, 0);
  }

  // looping over time bins
  for (auto pItr : stagePosInds) {
    int count = 0;
    double norm = 0;

    ///  mean calculation  ///
    for (auto sLitr : scanLgndrs) {
      // mean
      for (uint i=0; i<legSize; i++) {
        runMeans[pItr.second][i] += sLitr.second[pItr.second][i];
      }

      // counting non zero bins
      if (scanCounts[sLitr.first][pItr.second] > 0) {
        count++;
      }
      norm += 1; //scanCounts[sLitr.first][pItr.second];
    }

    // calculate mean
    for (uint i=0; i<legSize; i++) {
      runMeans[pItr.second][i] /= norm;
    }

    ///  Standard deviation calculation  ///
    norm = 0;
    for (auto sLitr : scanLgndrs) {
      for (uint i=0; i<legSize; i++) {
        runStdev[pItr.second][i]
            += std::pow((sLitr.second[pItr.second][i]
                    - runMeans[pItr.second][i]), 2);
        /*
            += scanCounts[sLitr.first][pItr.second]
                  *std::pow((sLitr.second[pItr.second][i]
                    /scanCounts[sLitr.first][pItr.second]
                    - runMeans[pItr.second][i]), 2);
                    */
      }
      norm += 1; //scanCounts[sLitr.first][pItr.second];
    }

    cout<<"\n\nTIME: "<<pItr.first<<endl<<endl;
    // calculate standard deviation
    for (uint i=0; i<legSize; i++) {
      runStdev[pItr.second][i] = std::sqrt(runStdev[pItr.second][i]/norm);
      cout<<runMeans[pItr.second][i]<<"/"<<runStdev[pItr.second][i]<<"\t"<<endl;
    }

    for (auto sLitr : scanLgndrs) {
      for (uint rInd=0; rInd<legSize; rInd++) {
        if (fabs(sLitr.second[pItr.second][rInd] - runMeans[pItr.second][rInd])
            < stdScale*runStdev[pItr.second][rInd]) {
          sLitr.second[pItr.second][rInd] = NANVAL;
        }
      }
    }

  }

  //plt.printXY(runMeans, "testRunMeans", maximum, "5e2");

  if (verbose) {
    std::cout << "Removing " << removePos.size() << " bins\n";
    std::cout << "Before removal: " << stagePosInds.size() << "  "
      << scanLgndrs.begin()->second.size() << "  "
      << scanCounts.begin()->second.size() << "  "
      << runMeans.size() << "  " << runStdev.size() << std::endl;
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



void alignClass::mergeScans() {
  ///  Calculating time delays from stage position  ///
  if (verbose) {
    std::cout << "calculating time delays: " << runMeans.size() << std::endl;
  }

  int NtimeSteps = stagePosInds.size();
  int tInd = 1;
  if (timeDelays) delete[] timeDelays;
  timeDelays = new double[NtimeSteps + 1];
  for (auto pItr : stagePosInds) {
    cout<<"stage: "<<pItr.first<<"   "<<pItr.second<<endl;
    timeDelays[tInd] = 2*pItr.first/(3e3);
    tInd++;
  }
  cout<<"timeZ: "<<timeZero<<endl;
  timeDelays[0] = 2*timeDelays[1] - timeDelays[2];
  for (int i=NtimeSteps; i>=0; i--) {
    cout<<"time: "<<timeDelays[i]<<"     ";
    timeDelays[i] -= timeZero + timeDelays[0];
    cout<<timeDelays[i]<<endl;
  }

  save::saveDat<double>(timeDelays, NtimeSteps + 1,
      "./plots/data/timeDelays["
      + to_string(stagePosInds.size() + 1) + "].dat");


  /////  Initialize variables  /////
  legendres.clear();
  legendres.resize(Nlegendres);
  for (uint ilg=0; ilg<Nlegendres; ilg++) {
    legendres[ilg].resize(stagePosInds.size());
    for (uint it=0; it<stagePosInds.size(); it++) {
      legendres[ilg][it].resize(NradLegBins, 0);
    }
  }

  azimuthalAvg.clear();
  azimuthalAvg.resize(stagePosInds.size());
  for (uint it=0; it<stagePosInds.size(); it++) {
    azimuthalAvg[it].resize(NradAzmBins, 0);
  }


  ///////////////////////////
  /////  Merging scans  /////
  ///////////////////////////
  int rInd;
  double legNorm, azmNorm;
  auto pItr = stagePosInds.begin();
  for (uint it=0; it<NtimeSteps; it++) {

    /////  Merging legendres  /////
    for (uint ilg=0; ilg<Nlegendres; ilg++) {
    cout<<"starting to fill"<<endl;
      cout<<"it: "<<it<<endl;
      cout<<"resized"<<endl;
      for (uint ir=0; ir<NradLegBins; ir++) {
        rInd = ilg*NradLegBins + ir;
        legNorm = 0;
        //for (auto sLiter : scanLgndrs) {
          //cout<<"filling: "<<ilg<<"  "<<it<<"  "<<ir<<"  "<<sLiter.first<<"  "<<sLiter.second[it][rInd]<<"  "<<sLiter.second[it][rInd] -     runMeans[it][rInd]<<"  "<< stdScale*runStdev[it][rInd]<<endl;
          //if (sLiter.second[it][rInd] != NANVAL) {
          //  legendres[ilg][it][ir] += sLiter.second[it][rInd];
          //  legNorm += scanCounts[sLiter.first][it];
          //}
        //}
        //legendres[ilg][it][ir] /= legNorm;
        //cout<<"filled for lg: "<<ilg<<endl;
      }
    }

    for (uint iazm=0; iazm<NradAzmBins; iazm++) {
      azmNorm = 0;
      for (auto sAzml : scanAzmAvg) {
        if ((std::isnan(azimuthalAvg[it][iazm]) || fabs(azimuthalAvg[it][iazm]) > 1e100 || fabs(azimuthalAvg[it][iazm]) < 1e-100) && (azimuthalAvg[it][iazm]!=0)) {
          cout<<"FOUNDNAN azml merging: "<<sAzml.first<<"  "<<it<<"  "<<iazm<<"  "<<azimuthalAvg[it][iazm]<<endl;
        }
        if (sAzml.second[it][iazm] != NANVAL) {
          azimuthalAvg[it][iazm] += sAzml.second[it][iazm];
          azmNorm += 1; //scanCounts[sAzml.first][it];
        }
      }
      if (azmNorm) {
        azimuthalAvg[it][iazm] /= azmNorm;
      }
      else {
        azimuthalAvg[it][iazm] = NANVAL;
      }
    }
    pItr++;
  }
  for (int i=0; i<azimuthalAvg[12].size(); i++) {
    cout<<"dbl check "<<i<<"  "<<azimuthalAvg[12][i]<<endl;
    if ((std::isnan(azimuthalAvg[12][i]) || fabs(azimuthalAvg[12][i]) > 1e100 || fabs(azimuthalAvg[12][i]) < 1e-100) && (azimuthalAvg[12][i]!=0)) {
      cout<<"FOUNDNAN azml merge: "<<"12 "<<"  "<<i<<"  "<<azimuthalAvg[12][i]<<endl;
    }
  }

  for (uint i=0; i<15; i++) {
    cout<<"CHECK: "<<i<<"  "<<azimuthalAvg[12][i]<<endl;
  }
}


void alignClass::getNrefBins() {
  ///  Number of images to average before T0  ///
  std::vector<double> unPumped(NradLegBins, 0);
  NrefBins = 0;
  if (subtractT0) {
    while (timeDelays[NrefBins+1] < 0) {
      NrefBins++;
    }
  }
  else {
    NrefBins = stagePosInds.size();
  }
  cout<<"NrefBins: "<<NrefBins<<endl;
}


void alignClass::subtractT0andNormalize() {

  if (NrefBins == -1) getNrefBins();

  std::vector<double> unPumped(legendres[0][0].size());
  std::vector<double> count(legendres[0][0].size());

  for (uint ilg=0; ilg<Nlegendres; ilg++) {
    // Mean/t0 subtraction and normalizing by atomic scattering
    // Raw data
    std::fill(unPumped.begin(), unPumped.end(), 0.0);
    std::fill(count.begin(), count.end(), 0.0);
    for (int ir=0; ir<NradLegBins; ir++) {
      for (int tm=0; tm<NrefBins; tm++) {
        if (legendres[ilg][tm][ir] != NANVAL) {
          unPumped[ir] += legendres[ilg][tm][ir];
          count[ir]++;
        }
      }
      if (count[ir]) {
        unPumped[ir] /= count[ir];
      }
      else {
        unPumped[ir] = 0;
      }
    }
    cout<<"unp filled"<<endl;
    for (int ir=0; ir<NradLegBins; ir++) {
      for (uint tm=0; tm<legendres[ilg].size(); tm++) {
        if (legendres[ilg][tm][ir] != NANVAL) {
          legendres[ilg][tm][ir] -= unPumped[ir];
          legendres[ilg][tm][ir] *= sMsLegNorm[ir];
        }
      }
    }
    cout<<"scaled"<<endl;
    if (ilg==0) {
      plt->print1d(unPumped, "testUnp");
    }

    ///  Saving data  ///
    /*
    save::saveDat<double>(legendres[ilg], "./plots/data/data_sMsL"
        + to_string(ilg) + "Diff["
        + to_string(legendres[ilg].size()) + ","
        + to_string(legendres[ilg][0].size()) + "].dat");
    save::saveDat<double>(legendres[ilg][legendres[ilg].size()-1], "./plots/data/data_sMsFinalL"
        + to_string(ilg) + "Diff["
        + to_string(legendres[ilg][0].size()) + "].dat");


    // Save raw unpumped data to fit background
    if (ilg == 0) {
      save::saveDat<double>(unPumped, "./qScale/results/unPumpedDiffractionL0["
          + to_string(unPumped.size()) + "].dat");
      plt.print1d(unPumped, "unPumpedRaw");
    }
    */
  }


  /////  Editting azimuthal average  /////
  count.clear();
  unPumped.clear();
  count.resize(NradAzmBins, 0);
  unPumped.resize(NradAzmBins, 0);
  for (int ir=0; ir<NradAzmBins; ir++) {
    for (int tm=0; tm<NrefBins; tm++) {
      if (azimuthalAvg[tm][ir] != NANVAL) {
        unPumped[ir] += azimuthalAvg[tm][ir];
        count[ir]++;
      }
    }
    if (count[ir]) {
      unPumped[ir] /= count[ir];
    }
    else {
      unPumped[ir] = 0;
    }
  }

  for (int ir=0; ir<NradAzmBins; ir++) {
    for (int tm=0; tm<stagePosInds.size(); tm++) {
      if (azimuthalAvg[tm][ir] != NANVAL) {
        azimuthalAvg[tm][ir] -= unPumped[ir];
        azimuthalAvg[tm][ir] *= sMsAzmNorm[ir];
      }
    }
  }

}


void alignClass::smearTime() {
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

  if (NrefBins == -1) getNrefBins();

  smearedImg.resize(Nlegendres);
  std::vector<double> unPumped(legendres[0][0].size());

  for (uint ilg=0; ilg<Nlegendres; ilg++) {

    double sum = 0;
    smearedImg[ilg].resize(stagePosInds.size());
    for (uint ir=0; ir<(int)stagePosInds.size(); ir++) {
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

    ///  Subtract before t0  ///
    std::fill(unPumped.begin(), unPumped.end(), 0.0);
    for (int ir=0; ir<NradLegBins; ir++) {
      for (int tm=0; tm<NrefBins; tm++) {
        unPumped[ir] += smearedImg[ilg][tm][ir];
      }
      unPumped[ir] /= (double)NrefBins;
    }
    for (int ir=0; ir<NradLegBins; ir++) {
      for (uint tm=0; tm<stagePosInds.size(); tm++) {
        smearedImg[ilg][tm][ir] -= unPumped[ir];
      }
    }

    if (ilg == 0) {
      save::saveDat<double>(unPumped, "./plots/data/data_unPumpedDiffractionL0Smeared["
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



















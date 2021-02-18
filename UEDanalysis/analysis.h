#include <cmath>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <ctime>

//From ROOT
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>

//From OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>

//Home Grown
#include "/cds/home/k/khegazy/baseTools/tools/tools.h"
#include "/cds/home/k/khegazy/baseTools/tools/imageProcessing.h"
#include "/cds/home/k/khegazy/baseTools/tools/plotClass.h"
//#include "/cds/home/k/khegazy/baseTools/tools/saveClass.h"
#include "/cds/home/k/khegazy/baseTools/tools/saving.h"
#include "/cds/home/k/khegazy/baseTools/tools/parameters.h"



using namespace std;

//////////////////////////////////////
//   Declaration of useful classes  //
//////////////////////////////////////

PLOTclass plt;


//////////////////////////////////////////
//   Declaration of analysis variables  //
//////////////////////////////////////////

std::vector<int>           scans;
std::vector<int>           *imgNums;
std::vector<int>           *imgIsRefs;
std::vector<int>           *timeStamps;
std::vector<int>           *stagePos;
std::vector<float>         *t0StagePos;
std::vector<float>         *t0Times;
std::vector<float>         *throttles;
std::vector<int>           *centerCs;
std::vector<int>           *centerRs;
std::vector<int>           *I0centerCs;
std::vector<int>           *I0centerRs;
std::vector<float>         *centerCstdRatios;
std::vector<float>         *centerRstdRatios;
std::vector<float>         *imgNorms;
std::vector<float>         *I0norms;
std::vector<float>         *readoutNoises;
std::vector<double>         *UVcounts;
std::vector<double>         *UVcountDers;
std::vector<double>         *bunkerTemps;
std::vector<double>         *bunkerTempDers;
std::vector<double>         *highBayTemps;
std::vector<double>         *highBayTempDers;
std::vector<double>         *pressures;
std::vector<double>         *pressureDers;
std::vector<vector<double> >  *legCoeffs;
std::vector<vector<int> >     *legCoeffs_nanMaps;
std::vector<vector<double> >  *rawAzmAvgs;
std::vector<vector<double> >  *azmAvgs;
std::vector<vector<int> >     *azmAvg_nanMaps;
std::vector<vector<double> >  *imgRadSTDs;
std::vector<vector<vector<double> > > *imgOrigs;
std::vector<vector<vector<double> > > *imgSubBkgs;


//////////////////////////////////////
//  analysisClass for loading data  // 
//////////////////////////////////////

class analysisClass : public parameterClass {

   public:
        analysisClass(std::string input_run) : parameterClass(input_run) {
          start             = clock();
          imgNums           = NULL;
          imgIsRefs         = NULL;
          timeStamps        = NULL;
          stagePos          = NULL;
          t0StagePos        = NULL;
          t0Times           = NULL;
          throttles         = NULL;
          centerCs          = NULL;
          centerRs          = NULL;
          I0centerCs        = NULL;
          I0centerRs        = NULL;
          centerCstdRatios  = NULL;
          centerRstdRatios  = NULL;
          imgNorms          = NULL;
          I0norms           = NULL;
          readoutNoises     = NULL;
          UVcounts          = NULL;
          UVcountDers       = NULL;
          bunkerTemps       = NULL;
          bunkerTempDers    = NULL;
          highBayTemps      = NULL;
          highBayTempDers   = NULL;
          pressures         = NULL;
          pressureDers      = NULL;
          legCoeffs         = NULL;
          legCoeffs_nanMaps = NULL;
          rawAzmAvgs        = NULL;
          azmAvgs           = NULL;
          azmAvg_nanMaps    = NULL;
          imgRadSTDs        = NULL;
          imgOrigs          = NULL;
          imgSubBkgs        = NULL;
        }


        void init_preProcesses_data(std::vector<std::string> variables);
        ~analysisClass();
        uint64_t setupEnvironment();
        int loadEvent(uint64_t entry);

	clock_t start, stop;

  private:
        std::vector<int>           _imgNums;
        std::vector<int>           _imgIsRefs;
        std::vector<int>           _timeStamps;
        std::vector<int>           _stagePos;
        std::vector<float>         _t0StagePos;
        std::vector<float>         _t0Times;
        std::vector<float>         _throttles;
        std::vector<int>           _centerCs;
        std::vector<int>           _centerRs;
        std::vector<int>           _I0centerCs;
        std::vector<int>           _I0centerRs;
        std::vector<float>         _centerCstdRatios;
        std::vector<float>         _centerRstdRatios;
        std::vector<float>         _imgNorms;
        std::vector<float>         _I0norms;
        std::vector<float>         _readoutNoises;
        std::vector<double>        _UVcounts;
        std::vector<double>        _UVcountDers;
        std::vector<double>        _bunkerTemps;
        std::vector<double>        _bunkerTempDers;
        std::vector<double>        _highBayTemps;
        std::vector<double>        _highBayTempDers;
        std::vector<double>        _pressures;
        std::vector<double>        _pressureDers;
        std::vector<vector<double> >  _legCoeffs;
        std::vector<vector<int> >     _legCoeffs_nanMaps;
        std::vector<vector<double> >  _rawAzmAvgs;
        std::vector<vector<double> >  _azmAvgs;
        std::vector<vector<int> >     _azmAvg_nanMaps;
        std::vector<vector<double> >  _imgRadSTDs;
        std::vector<vector<vector<double> > > _imgOrigs;
        std::vector<vector<vector<double> > > _imgSubBkgs;

};      


void analysisClass::init_preProcesses_data(std::vector<std::string> variables) {


  std::vector<int> sizes, _scans;
  std::vector<std::string> fileNames, var_types;
  std::vector<int> _scans_unique;
  // Get Scan Numbers
  std::string data_dir = preProcOutputDir + "Run-" + run + "/";
  DIR* dir = opendir(data_dir.c_str());
  struct dirent* ent;
  while ((ent = readdir(dir)) != NULL) {
    string curFileName(ent->d_name);
   
    if (curFileName.size() < 4) continue;

    if (curFileName.substr(0,4).compare("Scan") == 0) {
      auto spos = curFileName.find("_");
      int scan_num = stoi(curFileName.substr(5, curFileName.find("_") - 5));
      std::string variable_type = curFileName.substr(spos+1,
          curFileName.rfind("_") - (spos + 1));
      auto itr = std::find(variables.begin(), variables.end(),
          variable_type);

      if (itr != variables.end()) {
        fileNames.push_back(data_dir+curFileName);
        _scans.push_back(scan_num);
        var_types.push_back((*itr));
        if (std::find(_scans_unique.begin(), _scans_unique.end(), scan_num)
            == _scans_unique.end()) {
          _scans_unique.push_back(scan_num);
        }
        
        // Get Data Shape
        std::vector<int> shape = save::getShape(
            data_dir, curFileName.substr(0, curFileName.find("hape[")));
        sizes.push_back(shape[0]);
      }
    }
    else {
      std::cerr << "WARNING!!! Cannot handle file " 
          + curFileName << " ... Skipping!!!\n";
    }
  }
  closedir(dir);

  std::sort(_scans_unique.begin(), _scans_unique.end());

  // Order Indices
  for (uint isc=0; isc<_scans_unique.size(); isc++) {

    // Append entries to scan vector
    int scan_size = 0;
    {
      int ind=0;
      while (_scans_unique[isc] != _scans[ind]) ind++;
      scan_size = sizes[ind];
    }
        
    scans.resize(scans.size()+scan_size, _scans_unique[isc]);

    for (uint ifl=0; ifl<fileNames.size(); ifl++) {

      if (_scans[ifl] == _scans_unique[isc]) {
        std::string fName = fileNames[ifl];

        // Get Data Shape
        std::vector<int> shape = save::getShape(data_dir,
            fName.substr(data_dir.size(), fName.find("hape[")-data_dir.size()));

        cout<<fName<<endl;
        if (var_types[ifl].compare("imgNum") == 0) {
          int pSize = _imgNums.size();
          _imgNums.resize(pSize + shape[0]);
          save::importDat<int>(_imgNums, fName, pSize);
         
          if (!imgNums) imgNums = &_imgNums;
        }
        else if (var_types[ifl].compare("imgIsRef") == 0) {
          int pSize = _imgIsRefs.size();
          _imgIsRefs.resize(pSize + shape[0]);
          save::importDat<int>(_imgIsRefs, fName, pSize);
          
          if (!imgIsRefs) imgIsRefs = &_imgIsRefs;
        }
        else if (var_types[ifl].compare("timeStamp") == 0) {
          int pSize = _timeStamps.size();
          _timeStamps.resize(pSize + shape[0]);
          save::importDat<int>(_timeStamps, fName, pSize);
          
          if (!timeStamps) timeStamps = &_timeStamps;
        }
        else if (var_types[ifl].compare("stagePos") == 0) {
          int pSize = _stagePos.size();
          _stagePos.resize(pSize + shape[0]);
          save::importDat<int>(_stagePos, fName, pSize);
          
          if (!stagePos) stagePos = &_stagePos;
        }
        else if (var_types[ifl].compare("t0StagePos") == 0) {
          int pSize = _t0StagePos.size();
          _t0StagePos.resize(pSize + shape[0]);
          save::importDat<float>(_t0StagePos, fName, pSize);
          
          if (!t0StagePos) t0StagePos = &_t0StagePos;
        }
        else if (var_types[ifl].compare("t0Time") == 0) {
          int pSize = _t0Times.size();
          _t0Times.resize(pSize + shape[0]);
          save::importDat<float>(_t0Times, fName, pSize);
          
          if (!t0Times) t0Times = &_t0Times;
        }
        else if (var_types[ifl].compare("throttle") == 0) {
          int pSize = _throttles.size();
          _throttles.resize(pSize + shape[0]);
          save::importDat<float>(_throttles, fName, pSize);
          
          if (!throttles) throttles = &_throttles;
        }
        else if (var_types[ifl].compare("centerC") == 0) {
          int pSize = _centerCs.size();
          _centerCs.resize(pSize + shape[0]);
          save::importDat<int>(_centerCs, fName, pSize);
          
          if (!centerCs) centerCs = &_centerCs;
        }
        else if (var_types[ifl].compare("centerR") == 0) {
          int pSize = _centerRs.size();
          _centerRs.resize(pSize + shape[0]);
          save::importDat<int>(_centerRs, fName, pSize);
          
          if (!centerRs) centerRs = &_centerRs;
        }
        else if (var_types[ifl].compare("I0centerC") == 0) {
          int pSize = _I0centerCs.size();
          _I0centerCs.resize(pSize + shape[0]);
          save::importDat<int>(_I0centerCs, fName, pSize);
          
          if (!I0centerCs) I0centerCs = &_I0centerCs;
        }
        else if (var_types[ifl].compare("I0centerR") == 0) {
          int pSize = _I0centerRs.size();
          _I0centerRs.resize(pSize + shape[0]);
          save::importDat<int>(_I0centerRs, fName, pSize);
          
          if (!I0centerRs) I0centerRs = &_I0centerRs;
        }
        else if (var_types[ifl].compare("centerCstdRatio") == 0) {
          int pSize = _centerCstdRatios.size();
          _centerCstdRatios.resize(pSize + shape[0]);
          save::importDat<float>(_centerCstdRatios, fName, pSize);
          
          if (!centerCstdRatios) centerCstdRatios = &_centerCstdRatios;
        }
        else if (var_types[ifl].compare("centerRstdRatio") == 0) {
          int pSize = _centerRstdRatios.size();
          _centerRstdRatios.resize(pSize + shape[0]);
          save::importDat<float>(_centerRstdRatios, fName, pSize);
          
          if (!centerRstdRatios) centerRstdRatios = &_centerRstdRatios;
        }
        else if (var_types[ifl].compare("imgNorm") == 0) {
          int pSize = _imgNorms.size();
          _imgNorms.resize(pSize + shape[0]);
          save::importDat<float>(_imgNorms, fName, pSize);
          
          if (!imgNorms) imgNorms = &_imgNorms;
        }
        else if (var_types[ifl].compare("I0norm") == 0) {
          int pSize = _I0norms.size();
          _I0norms.resize(pSize + shape[0]);
          save::importDat<float>(_I0norms, fName, pSize);
          
          if (!I0norms) I0norms = &_I0norms;
        }
        else if (var_types[ifl].compare("readoutNoise") == 0) {
          int pSize = _readoutNoises.size();
          _readoutNoises.resize(pSize + shape[0]);
          save::importDat<float>(_readoutNoises, fName, pSize);
         
          if (!readoutNoises) readoutNoises = &_readoutNoises;
        }
        else if (var_types[ifl].compare("UVcounts") == 0) {
          int pSize = _UVcounts.size();
          _UVcounts.resize(pSize + shape[0]);
          save::importDat<double>(_UVcounts, fName, pSize);
         
          if (!UVcounts) UVcounts = &_UVcounts;
        }
        else if (var_types[ifl].compare("UVcountsDer") == 0) {
          int pSize = _UVcountDers.size();
          _UVcountDers.resize(pSize + shape[0]);
          save::importDat<double>(_UVcountDers, fName, pSize);
         
          if (!UVcountDers) UVcountDers = &_UVcountDers;
        }
        else if (var_types[ifl].compare("bunkerTemp") == 0) {
          int pSize = _bunkerTemps.size();
          _bunkerTemps.resize(pSize + shape[0]);
          save::importDat<double>(_bunkerTemps, fName, pSize);
         
          if (!bunkerTemps) bunkerTemps = &_bunkerTemps;
        }
        else if (var_types[ifl].compare("bunkerTempDer") == 0) {
          int pSize = _bunkerTempDers.size();
          _bunkerTempDers.resize(pSize + shape[0]);
          save::importDat<double>(_bunkerTempDers, fName, pSize);
         
          if (!bunkerTempDers) bunkerTempDers = &_bunkerTempDers;
        }
        else if (var_types[ifl].compare("highBayTemp") == 0) {
          int pSize = _highBayTemps.size();
          _highBayTemps.resize(pSize + shape[0]);
          save::importDat<double>(_highBayTemps, fName, pSize);
         
          if (!highBayTemps) highBayTemps = &_highBayTemps;
        }
        else if (var_types[ifl].compare("highBayTempDer") == 0) {
          int pSize = _highBayTempDers.size();
          _highBayTempDers.resize(pSize + shape[0]);
          save::importDat<double>(_highBayTempDers, fName, pSize);
         
          if (!highBayTempDers) highBayTempDers = &_highBayTempDers;
        }
        else if (var_types[ifl].compare("pressure") == 0) {
          int pSize = _pressures.size();
          _pressures.resize(pSize + shape[0]);
          save::importDat<double>(_pressures, fName, pSize);
         
          if (!pressures) pressures = &_pressures;
        }
        else if (var_types[ifl].compare("pressureDer") == 0) {
          int pSize = _pressureDers.size();
          _pressureDers.resize(pSize + shape[0]);
          save::importDat<double>(_pressureDers, fName, pSize);
         
          if (!pressureDers) pressureDers = &_pressureDers;
        }
        else if (var_types[ifl].compare("legCoeffs") == 0) {
          int pSize = _legCoeffs.size();
          _legCoeffs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_legCoeffs.size(); i++) {
            _legCoeffs[i].resize(shape[1]);
          }
          save::importDat<double>(_legCoeffs, fName, pSize);
         
          if (!legCoeffs) legCoeffs = &_legCoeffs;
        }
        else if (var_types[ifl].compare("legCoeffs_nanMap") == 0) {
          int pSize = _legCoeffs_nanMaps.size();
          _legCoeffs_nanMaps.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_legCoeffs_nanMaps.size(); i++) {
            _legCoeffs_nanMaps[i].resize(shape[1]);
          }
          save::importDat<int>(_legCoeffs_nanMaps, fName, pSize);
         
          if (!legCoeffs_nanMaps) legCoeffs_nanMaps = &_legCoeffs_nanMaps;
        }
        else if (var_types[ifl].compare("rawAzmAvg") == 0) {
          int pSize = _rawAzmAvgs.size();
          _rawAzmAvgs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_rawAzmAvgs.size(); i++) {
            _rawAzmAvgs[i].resize(shape[1]);
          }
          save::importDat<double>(_rawAzmAvgs, fName, pSize);
         
          if (!rawAzmAvgs) rawAzmAvgs = &_rawAzmAvgs;
        }
        else if (var_types[ifl].compare("azmAvg") == 0) {
          int pSize = _azmAvgs.size();
          _azmAvgs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_azmAvgs.size(); i++) {
            _azmAvgs[i].resize(shape[1]);
          }
          save::importDat<double>(_azmAvgs, fName, pSize);
         
          if (!azmAvgs) azmAvgs = &_azmAvgs;
        }
        else if (var_types[ifl].compare("azmAvg_nanMap") == 0) {
          cout<<"wtf1"<<endl;
          int pSize = _azmAvg_nanMaps.size();
          _azmAvg_nanMaps.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_azmAvg_nanMaps.size(); i++) {
            _azmAvg_nanMaps[i].resize(shape[1]);
          }
          save::importDat<int>(_azmAvg_nanMaps, fName, pSize);
         
          if (!azmAvg_nanMaps) azmAvg_nanMaps = &_azmAvg_nanMaps;
          cout<<"wtf1"<<endl;
        }
        else if (var_types[ifl].compare("imgRadSTD") == 0) {
          int pSize = _imgRadSTDs.size();
          _imgRadSTDs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_imgRadSTDs.size(); i++) {
            _imgRadSTDs[i].resize(shape[1]);
          }
          save::importDat<double>(_imgRadSTDs, fName, pSize);
         
          if (!imgRadSTDs) imgRadSTDs = &_imgRadSTDs;
        }
        else if (var_types[ifl].compare("imgOrig") == 0) {
          int pSize = _imgOrigs.size();
          _imgOrigs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_imgOrigs.size(); i++) {
            _imgOrigs[i].resize(shape[1]);
            for (int j=0; j<shape[1]; j++) {
              _imgOrigs[i][j].resize(shape[2]);
            }
          }
          save::importDat<double>(_imgOrigs, fName, pSize);
         
          if (!imgOrigs) imgOrigs = &_imgOrigs;
        }
        else if (var_types[ifl].compare("imgSubBkg") == 0) {
          int pSize = _imgSubBkgs.size();
          _imgSubBkgs.resize(pSize + shape[0]);
          for (int i=pSize; i<(int)_imgSubBkgs.size(); i++) {
            _imgSubBkgs[i].resize(shape[1]);
            for (int j=0; j<shape[1]; j++) {
              _imgSubBkgs[i][j].resize(shape[2]);
            }
          }
          save::importDat<double>(_imgSubBkgs, fName, pSize);
         
          if (!imgSubBkgs) imgSubBkgs = &_imgSubBkgs;
        }
        else {
          std::cerr << "WARNING: Cannot handle variable type "
              << var_types[ifl] << " ... Skipping!!!\n";
        }
      }
    }
  }
 
  if (std::find(variables.begin(), variables.end(), "legCoeffs") != variables.end()) {
    std::cerr << "WARNING: Setting legCoeffs to azmAvg since preProce can't calculate\n";

    for (auto a : _azmAvgs) {
      _legCoeffs.push_back(a);
    }
    legCoeffs = &_legCoeffs;
  }
  if (std::find(variables.begin(), variables.end(), "legCoeffs_nanMap") != variables.end()) {
    std::cerr << "WARNING: Setting legCoeffs to azmAvg since preProce can't calculate\n";

    for (auto a : _azmAvg_nanMaps) {
      _legCoeffs_nanMaps.push_back(a);
    }
    legCoeffs_nanMaps = &_legCoeffs_nanMaps;
  }


}


analysisClass::~analysisClass() {

  //delete fChain;
  stop = clock();
  cout<<endl<<"Total running time: "<<double(stop-start)/CLOCKS_PER_SEC<<endl<<endl;
}


#ifndef PREPROCESSING_H
#define PREPROCESSING_H

#include "/reg/neh/home/khegazy/baseScripts/imageProcessing.h"
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"
#include "/reg/neh/home/khegazy/baseScripts/tools.h"
#include "../../parameters.h"
#include <fftw3.h>
#include <dirent.h>

namespace ppFunct {

  bool getScanRunInfo(std::vector<imgProc::imgInfoStruct> &imgINFO, 
            std::string runListName, bool verbose);
  void makeRunLists(std::vector<imgProc::imgInfoStruct> &imgINFO,
            std::string runName, std::string preProcFolder);

  extern std::map<int, int> monthLengths;

}

/*
ppFunct::monthLengths[1] = 31;
ppFunct::monthLengths[2] = 28;
ppFunct::monthLengths[3] = 31;
ppFunct::monthLengths[4] = 30;
ppFunct::monthLengths[5] = 31;
ppFunct::monthLengths[6] = 30;
ppFunct::monthLengths[7] = 31;
ppFunct::monthLengths[8] = 31;
ppFunct::monthLengths[9] = 30;
ppFunct::monthLengths[10] = 31;
ppFunct::monthLengths[11] = 30;
ppFunct::monthLengths[12] = 31;
*/
#endif

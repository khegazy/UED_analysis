#include "preProcessing.h"



/////////////////////////////
/////  Days in a Month  /////
/////////////////////////////

std::map<int, int> ppFunct::monthLengths = {
  {1 , 31},
  {2 , 28},
  {3 , 31},
  {4 , 30},
  {5 , 31},
  {6 , 30},
  {7 , 31},
  {8 , 31},
  {9 , 30},
  {10 , 31},
  {11 , 30},
  {12 , 31}};


////////////////////////////////
////  Retrieving File Info  ////
////////////////////////////////

bool ppFunct::getScanRunInfo(std::vector<imgProc::imgInfoStruct> &imgINFO, std::string runListName, bool verbose) {

  /*
monthLengths[1] = 31;
monthLengths[2] = 28;
monthLengths[3] = 31;
monthLengths[4] = 30;
monthLengths[5] = 31;
monthLengths[6] = 30;
monthLengths[7] = 31;
monthLengths[8] = 31;
monthLengths[9] = 30;
monthLengths[10] = 31;
monthLengths[11] = 30;
monthLengths[12] = 31;
*/


  ifstream fileNames(runListName.c_str());
  if (!fileNames.is_open()) {
    cerr << "ERROR: Cannot open file " << runListName << endl;
    exit(0);
  }

  size_t ipos, fpos;
  string fileName;
  std::map<float, imgProc::imgInfoStruct> imgInfoMap;

  std::string curRun="";
  std::string curDate="";
  std::string stgP1, stgP2;
  int curScan=-1;

  std::vector<int> months;
  std::vector<int> days;
  std::vector<int> times;
  int time;

  while (getline(fileNames, fileName)) {

    if (verbose) std::cout << "Now looking at file " << fileName << endl;

    imgProc::imgInfoStruct imgInfo;

    // Finding date (2016)
    ipos = fileName.find("201",0);
    imgInfo.date = fileName.substr(ipos, 8);
    if (verbose) std::cout << "\tDATE: " << imgInfo.date;

    // Finding Scan name
    ipos = fileName.find("scan",0);
    if (ipos == string::npos) {
      ipos = fileName.find("Scan",0);
      if (ipos == string::npos) {
        ipos = fileName.find("Range",0);
        if (ipos == string::npos) {
          cerr << "ERROR: Cannot find 'scan' or 'Scan' in image name!!!\n\n";
        }
      }
    }
    if (ipos != string::npos) {
      if (imgInfo.date.find("2016") != string::npos) {
        while (fileName.at(ipos-1) != '/') ipos--;
        fpos = fileName.find("/", ipos);
        imgInfo.run = fileName.substr(ipos, fpos-ipos);
      }
      else {
        while (fileName.at(ipos+3) != '/') ipos++;
        fpos = fileName.find("/", ipos);
        imgInfo.scan = std::stoi(fileName.substr(ipos, fpos-ipos));
      }
    }
    else {
      imgInfo.scan = 1;
    }

    // Finding run name
    ipos = fileName.find("run",0);
    if (ipos == string::npos) {
      ipos = fileName.find("Run",0);
      if (ipos == string::npos) {
        if (fileName.find("20180629_1630",0) != string::npos) {
          ipos = fileName.find("Background",0);
          ipos += 7;
        }
      }
    }
    imgInfo.runType = "Run";
    if (ipos == string::npos) {
      ipos = fileName.find("Background",0);
      ipos += 7;
      imgInfo.runType = "Background";
    }
    if (ipos == string::npos) {
      std::cerr << "ERROR: Cannot find Run or Background indication!!!\n";
      exit(0);
    }

    if (imgInfo.date.find("2016") != string::npos) {
      fpos = fileName.find("/", ipos);
      imgInfo.scan = stoi(fileName.substr(ipos+3, fpos-ipos-3));
    }
    else {
      fpos = fileName.find("/", ipos+4);
      imgInfo.run = fileName.substr(ipos+4, fpos-ipos-4);
    }
    if (verbose) cout << "\tRUN: " << imgInfo.run;
    if (verbose) cout << "\tSCAN: " << imgInfo.scan;

    // Finding file name
    ipos = fileName.rfind("/");
    imgInfo.fileName = fileName.substr(ipos, fileName.length()-ipos);
    imgInfo.path = fileName.substr(0, ipos);
    if (verbose) cout << "\tFILEPATH: " << imgInfo.path;
    if (verbose) cout << "\tFILENAME: " << imgInfo.fileName;

    // Finding image number
    ipos = imgInfo.fileName.find("-");
    imgInfo.imgNum = stoi(imgInfo.fileName.substr(ipos+1, 3));
    if (verbose) cout << "\tIMGNUM: " << imgInfo.imgNum;

    // Finding stage position
    stgP1 = imgInfo.fileName.substr(imgInfo.fileName.length()-17, 3);
    stgP2 = imgInfo.fileName.substr(imgInfo.fileName.length()-13, 4);
    imgInfo.stagePos = stoi(stgP1+stgP2);
    if (verbose) cout << "\tSTAGE POSITION: " << imgInfo.stagePos;


    // Fill imgINFO with ordered events
    if ((curScan != imgInfo.scan) || (curRun != imgInfo.run) || (curDate != imgInfo.date)) {
     
      // Get time image was taken
      times.clear();
      DIR* dir = opendir(imgInfo.path.c_str());
      struct dirent* ent;
      while ((ent = readdir(dir)) != NULL) {
        string txtName(ent->d_name);
        if (txtName.length() < 10) continue;
        if (txtName.substr(txtName.length()-4, 4).compare(".txt") == 0) {
          time = 0;
          for (int i=1; i<stoi(txtName.substr(26,2)); i++) {
            time += monthLengths[i]*24*3600;
          }
          time += stoi(txtName.substr(28,2))*24*3600;
          time += stoi(txtName.substr(31,2))*3600;
          time += stoi(txtName.substr(33,2))*60;
          time += stoi(txtName.substr(35,2));
          times.push_back(time);
        }
      }

      sort(times.begin(), times.end(),
            [](long int i1, long int i2)
            {return i1 < i2;});

      if (verbose) cout << "\n\nFilling imgINFO\n";

      for (auto itr: imgInfoMap) {
        itr.second.time = times[itr.second.imgNum-1];
        imgINFO.push_back(itr.second);
      }
      imgInfoMap.clear();

      curScan = imgInfo.scan;
      curRun = imgInfo.run;
      curDate = imgInfo.date;
    }

    if (verbose) cout << endl;
    imgInfoMap[imgInfo.stagePos] = imgInfo;
  }

  for (auto itr: imgInfoMap) {
    itr.second.time = times[itr.second.imgNum-1];
    cout<<itr.second.time<<endl;
    imgINFO.push_back(itr.second);
  }

  fileNames.close();
  if (verbose) cout << "\nINFO: Finished retrieving info from runList!\n\n";

  // Check stagePos is ordered
  if (verbose) {
    cout << "Checking images are in order of stage position\n";
    for (uint i=0; i<imgINFO.size(); i++) {
      cout << imgINFO[i].scan << "  " << imgINFO[i].run 
        << "  " << imgINFO[i].stagePos << endl;
    }
  }


  return true;
}



//////////////////////////////////////////////////////
/////  Creating smaller runList files if needed  /////
//////////////////////////////////////////////////////

void ppFunct::makeRunLists(std::vector<imgProc::imgInfoStruct> &imgINFO,
      std::string runName, std::string preProcFolder) {
  ofstream outList;
  std::string curRun   = ""; 
  std::string curDate  = ""; 
  std::string fileName = "";

  int maxScan = 0;
  int curScan = -1;
  for (uint ifl=0; ifl<imgINFO.size(); ifl++) {

    if ((imgINFO[ifl].run.compare("20180629_1630") == 0) 
        && (imgINFO[ifl].stagePos == 1545800)) continue;

    if ((curScan != imgINFO[ifl].scan) || (curRun != imgINFO[ifl].run)
          || (curDate != imgINFO[ifl].date)) {

      curRun = imgINFO[ifl].run;
      curDate = imgINFO[ifl].date;
      curScan = imgINFO[ifl].scan;
      if (curScan > maxScan) {
        maxScan = curScan;
      }

      if (outList.is_open()) {
        outList.close();
      }

      fileName = "runLists/runList_" + imgINFO[ifl].runType + "-";
      if (curDate.find("2016") != string::npos) {
        fileName += imgINFO[ifl].date + "_";
      }
      fileName += imgINFO[ifl].run + "_Scan-" + to_string(imgINFO[ifl].scan) + ".txt";

      outList.open(fileName.c_str());
    }
    cout<<fileName<<endl;
    cout<<imgINFO[ifl].path + imgINFO[ifl].fileName<<endl;
    outList << imgINFO[ifl].path + imgINFO[ifl].fileName << endl;
  }
  outList.close();


  /////  Runlist of preprocessed files  /////
  if (imgINFO[0].runType.compare("Background") != 0) {
    outList.open(("../mergeScans/runLists/run-" + imgINFO[0].run + ".txt").c_str());
    for (int i=1; i<=maxScan; i++) {
      outList << preProcFolder + "Run-" 
                  + imgINFO[0].run + "_Scan-"
                  + to_string(i) + ".root" << endl;
    }
    outList.close();
  }

  // Exiting program so you can use correct runLists
  exit(1);
}




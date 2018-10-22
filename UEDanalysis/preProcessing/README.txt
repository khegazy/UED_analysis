


Creating Runlist:
  Make a single list of all the files
      ls /reg/ued/ana/data/nitrobenzene/20180627/Run/20180627_1551/scan*/image*/*tif > list.txt

  Set doFunLists in preProcessing.cpp to true
      bool doRunLists = true;

  Compile code and run
      ./preProcessing.exe list.txt

  Reset doRunLists to false and delete list.txt


Finding course center:
  Set imgProc::centerSearchCOM::verbose to true, it will output 
    plots/data/centering_selectedRange[1024,1024].dat
    plots/data/centering_smoothed[1024,1024].dat
    
  run plotting script in plots/scripts
    python courseCOMcenterFinding.py

  View files in plots folder and change parameters accordingly

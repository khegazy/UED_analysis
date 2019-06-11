#include "../analysis.h"
#include "/reg/neh/home/khegazy/baseTools/tools/parameters.h"


int main(int argc, char* argv[]) {

  string runName(argv[1]);
  cout<<"RunName: "<<runName<<endl;

  ///  Get Parameters  ///
  parameterClass params(runName);

  int fitIndBegin = (int)(params.NradAzmBins*params.fitQbegin/params.maxQazm);
  int fitIndEnd   = (int)(params.NradAzmBins*params.fitQend/params.maxQazm);


  ///  Simulation/data data structs  ///
  std::vector<double> atmDiff(params.NradAzmBins);
  std::vector<double> initStateMolDiff(params.NradAzmBins);
  std::vector<double> finalStateMolDiff(params.NradAzmBins);
  std::vector<double> dataFinalState(params.NradAzmBins);
  std::vector<double> dataReference(params.NradAzmBins);

  // Getting initial state simulation
  save::importDat<double> (atmDiff, 
      params.simOutputDir + params.molName 
      + "_atmDiffractionPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(params.NradAzmBins) + "].dat");
 
  save::importDat<double> (initStateMolDiff, 
      params.simOutputDir + params.molName 
      + "_molDiffractionPatternLineOut_Qmax-"
      + to_string(params.maxQazm) + "_Ieb-"
      + to_string(params.Iebeam) + "_scrnD-"
      + to_string(params.screenDist) + "_elE-"
      + to_string(params.elEnergy) + "_Bins["
      + to_string(params.NradAzmBins) + "].dat");

  // Getting data 
  save::importDat<double> (dataFinalState,
      params.mergeScansOutputDir 
      + "data-" + runName + 
      + "_diffFinalState[" 
      + to_string(params.NradAzmBins) + "].dat");
 
  save::importDat<double> (dataReference,
      params.mergeScansOutputDir 
      + "data-" + runName + 
      + "-referenceAzm[" 
      + to_string(params.NradAzmBins) + "].dat");
 

  ///////////////////////////////////////
  /////  Loop Through Final States  /////
  ///////////////////////////////////////

  for (int ifs=0; ifs<params.finalStates.size(); ifs++) {

    std::cout << "Starting Final State: " + params.finalStates[ifs] + "\n";

    ///  Get molecular diffraction of final state  ///
    save::importDat<double> (finalStateMolDiff, 
        params.simOutputDir + params.finalStates[ifs] 
        + "_molDiffractionPatternLineOut_Qmax-"
        + to_string(params.maxQazm) + "_Ieb-"
        + to_string(params.Iebeam) + "_scrnD-"
        + to_string(params.screenDist) + "_elE-"
        + to_string(params.elEnergy) + "_Bins["
        + to_string(params.NradAzmBins) + "].dat");

    ///  Calculate percentage of transient signal  ///
    double atmDiffSum     = 0;
    double initMolDiffSum = 0;
    double diffMolDiffSum = 0;
    double finalStateDiffSum = 0;
    double referenceDiffSum = 0;
    for (int iq=fitIndBegin; iq<fitIndEnd; iq++) {
      atmDiffSum      += std::fabs(atmDiff[iq]);
      initMolDiffSum  += std::fabs(initStateMolDiff[iq]);
      diffMolDiffSum  += std::fabs(finalStateMolDiff[iq] - initStateMolDiff[iq]);

      finalStateDiffSum += std::fabs(dataFinalState[iq]);
      referenceDiffSum  += std::fabs(dataReference[iq]);
    }

    std::cout << "\t Percent of initial state diffraction is molecular: "
      << 100*initMolDiffSum/(initMolDiffSum + atmDiffSum) << "% \n";

    std::cout << "\t Percent of molecular signal that changes: "
      << 100*diffMolDiffSum/initMolDiffSum << "% \n";

    std::cout << "\t Percent of total molecular signal that changes: "
      << 100*diffMolDiffSum/(initMolDiffSum + atmDiffSum) << "% \n";

    std::cout << "\t Percent of total data that changes: "
      << 100*finalStateDiffSum/referenceDiffSum << "% \n";

    std::cout << "\t Percent of molecules excited: "
      << 100*(finalStateDiffSum/referenceDiffSum)
        /(diffMolDiffSum/(initMolDiffSum + atmDiffSum)) << "% \n";
  }


  return 0;
}

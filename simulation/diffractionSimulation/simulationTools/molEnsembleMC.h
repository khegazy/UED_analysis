#ifndef MOLENSEMBLEMC_H
#define MOLENSEMBLEMC_H

#include <string>
#include <vector>
#include <TRandom3.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <fstream>
#include <sstream>
#include "Rtypes.h"

//From OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/core/types_c.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>

//Home Grown
#include "../../../tools/tools.h"
#include "../../../tools/constants.h"
#include "../../../tools/plotClass.h"
#include "atomClass.h"
#include "moleculeClass.h"
#include "pdfClass.h"


using namespace std;


class MOLENSEMBLEMCclass : public PDFclass{

   public:
	MOLENSEMBLEMCclass(long int seed, std::string pdfPath, std::string pdfNames);
	MOLENSEMBLEMCclass(long int seed, std::string xyzPath);
	MOLENSEMBLEMCclass(long int seed);
	~MOLENSEMBLEMCclass();

	int Nmols;
	int NmolAtoms;
	MOLECULEclass* molecules;
	vector<ATOMS> atomTypes;
        bool verbose;
        
	void positionMolecule(MOLECULEclass &molecule);
	bool usePositionMC;
        double stdDevPositionLength;
	string XpositionPDF;
	string YpositionPDF;
	string ZpositionPDF;

	void orientMolecule(MOLECULEclass &molecule);
	bool useOrientationMC;
	string orientationPDF;

	TH2F* orientDist;
	TGraph2D* orientGraph;

	virtual void buildMolecule(MOLECULEclass &molecule, 
              std::map<std::string, double> inpVals);
	TGraph2D* testBuildMolecule();
	TGraph2D* testBuildMolecule(std::map<std::string, double> inpVals);

	// creates molecular ensemble with only position of the base atoms specified
	void makeMolEnsemble();
	void makeMolEnsemble(std::map<std::string, double> inpVals);

        std::map< std::string, std::vector<double> >* getBonds(
            std::map< std::string, std::vector<double> >* bonds=NULL);
        std::vector<double> simulatePairCorr(
            int Nbins, 
            float maxR,
            bool smear=false,
            PLOTclass* pltVerbose=NULL);

	void reset();

   private:
	TTree* saveTree;

	TH1F* testAzm;
	TH1F* testPol;

 	int Ngr;

        std::string xyzFile;
        std::map< std::string, int > atomCount;
        std::map< std::string, Eigen::Vector3d > atomPos;

        std::map< ATOMS, int > charges;

	void initialize();
        void importXYZfile(std::string fileName);
        ATOMS getATOMtype(std::string atm);
};

#endif

#ifndef DIFFRACTIONCLASS_H
#define DIFFRACTIONCLASS_H

#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <stdlib.h>
#include <string>
#include <TH2F.h>
#include <fftw3.h>

#include "../../../tools/constants.h"
#include "moleculeClass.h"
#include "molEnsembleMC.h"

//#include <boost/math/special_functions/bessel.hpp>

using namespace std;

class DIFFRACTIONclass {

   public:
	DIFFRACTIONclass(double Iebeam_in, double elEnergy, string scatAmpPath);
	DIFFRACTIONclass(MOLENSEMBLEMCclass* molEnsemble, double sMax_in, double Iebeam_in, double screenDist_in, double elEnergy, int Nbins_in, string scatAmpPath);
	DIFFRACTIONclass(string name, MOLENSEMBLEMCclass* molEnsemble, double sMax_in, double Iebeam_in, double screenDist_in, double elEnergy, int Nbins_in, string scatAmpPath);
	~DIFFRACTIONclass();

	vector< vector<double> > diffPattern;
	vector< vector<double> > diffMolPattern;
	vector< vector<double> > diffAtmPattern;
	vector< vector<double> > sPattern;

	TH2F* FT2D_samplePattern;
	TH2F* FT2D_FTpattern;

	double Iebeam;
	double screenDist;                      // m
	double screenSize;                      // m
	double lambda;				// AU
	double k0;				// AU
 	double sMax;
	int Nmols;

	MOLECULEclass* molecules;    

	map<ATOMS, vector<float> > scatAmps;
	//map<ATOMS, boost::math::barycentric_rational<double>* > scatAmps;
        double interpScatAmp(ATOMS atmT, double sInp);
 
        void sPatternCalc();
	void diffPatternCalc();
        std::vector< std::vector<double> > lineOut_uniform();
        std::vector< std::vector<double> > azimuthalAvg_uniform();
        std::vector< std::vector<double> > diffPatternCalc_uniform();
	//void diffPatternCalc(MOLECULEclass* molecules_in);
	void FT2d();
        Eigen::Vector3d sCalc(double xpos, double ypos, double zpos);
	void reset();

	void Initialize(string name, MOLENSEMBLEMCclass* molEnsemble, 
            double sMax_in, double Iebeam_in, double screenDist_in, 
            double elEnergy, int Nbins_in, string scatAmpPath);

   private:
  	int Nbins;
        bool fillLineOut;
	double xScreenInit, zScreenInit, xScreenFin, zScreenFin;
	double xSampleInit, zSampleInit, xSampleFin, zSampleFin;
	vector<float> scatSInterp;
	vector<ATOMS> atomTypes;

        void importScatteringAmps(std::string scatAmpPath);
};

#endif

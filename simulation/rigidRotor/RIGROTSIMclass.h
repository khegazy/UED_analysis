#ifndef RIGIDROTORCLASS_H
#define RIGIDROTORCLASS_H

#include <Sparse>
#include <Core>
#include <Eigenvalues>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <dirent.h>
#include <Math/SpecFuncMathMore.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include "/reg/neh/home/khegazy/baseTools/tools/tools.h"

#include "distributions.h"

using namespace std;


class RIGROTSIMclass {

  
   public:
        double startTime;            	/* ps */        
        double endTime;            	/* ps */        
                                                        
        double temperature;         	/* K */         
        double laserIntensity;        	/* W/cm^2 */    
                                                        
        double rotConstB;        	/* cm^(-1) */   
        double rotConstD;          	/* cm^(-1) */   
        double rotConstH;          	/* cm^(-1) */   
        double deltaAlpha;   		/* m^3 */       
                                                        
        double startSampleTime;      	/* ps */       
        double endSampleTime;       	/* ps */       
        double sampleStep;        	/* ps */       
        double dTimeEvolStep;     	/* ps */       
        double pulseTlength;       	/* ps */       
        uint Npulses;                                   
        double pulseTspacing;   	/* ps */       
        
        bool doCosSqEV;
        bool doSinCosEV;
        uint cosOrder, sinCosOrder;
        bool makePDFs;                                 
        std::string savePDFformat;
        double startPDFTime;		/* ps */
        double endPDFTime;          	/* ps */       
        int NYlmEVs;
        bool evenOnlyAxisDist;
        std::string axisDistWFPrefix;
        std::vector<double> sampleTimes;
        std::vector< std::vector<double> > cosSqEVals;
        std::vector< std::vector<double> > sinCosEVals;
        std::vector< std::vector<double> > YlmEVals;
        std::vector< std::vector< std::vector< std::vector< complex<double> > > > > axisDistWFcoeffs;
                                                       
        int MAXj;                                      
	bool hasQuarterRev; 
        double indexRefr;                              
        string vibKey;

        std::string fileName;
        std::string outputDir;
        std::string PDFfileName;
        std::string PDFoutputDir;

        std::vector<double> population0;
        std::vector<double> jPopulation;
        std::vector<double> vibThermalDist;	
        std::vector< std::vector<double> > rotThermalDist;	
        std::map< int, std::map<int, std::map<int, 
        Eigen::SparseVector< complex<double> > > > > eigenAmps;      // Eigenvectors of |VJm>
        map<int, vector<double> > cosSqEigVal;
        std::map<int, Eigen::SparseMatrix< complex<double> > > sinCos;
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSq;
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSqExpDiag;
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSqTran;
        std::vector< std::map<int, Eigen::SparseMatrix< complex<double> > > > Ylm;

        RIGROTSIMclass();
        void runSimulation();
        std::vector<double> getPopulationDistribution();
        std::map<int, std::vector< complex<double> > > projectDensitySH();



    private:
  	class prunefnctr {

    	    public:
        	double cutoff = 1e-10;
        	complex<double> complexNumber;
        	bool operator() (const Eigen::SparseMatrix< complex<double> >::Index& row, 
			const Eigen::SparseMatrix< complex<double> >::Index& col, 
			const Eigen::SparseMatrix< complex<double> >::Scalar& value) const {
          	  //complexNumber = value;
          	  return (fabs(value.real()) > cutoff) || (fabs(value.imag()) > cutoff);
        	}
	};

        prunefnctr prunefx;

	//PLOTclass pltRR;
	double E0;
	double startTimeAU;
  	double endTimeAU;
	double pulseTlengthAU;
	double pulseTspacingAU;
	double startSampleTimeAU;
	double endSampleTimeAU;
	double dTimeEvolStepAU;
	double rotConstB_AU;
	double rotConstD_AU;
	double rotConstH_AU;
	double deltaAlphaAU;

        bool doVib;
        int NvibStates;
        double vibPopThresh;

	void checkAmps();
        std::vector< complex<double> > EVcos();
        std::vector< complex<double> > EVsinCos();
        std::vector< complex<double> > EVYlm();
        std::vector< std::vector<double> > anglePDFCalc();
	double laserPotential(double time);
	double timeInt(double fTime, double iTime);
	void sampleSim(double time_fs);

        void calculateThermalDistribution();
        void calculateCosSqMatrices();
        void calculateSinCosMatrices();
        void calculateYlmMatrices();

        std::vector<uint> pdfSampleInds;
        std::vector<float> pdfSampleTimes;
        std::vector<double> pulseTimes;

  	uint ist = 0;         // index of sampleTimes
  	uint ipt = 0;         // index of pulseTimes
  	uint ipi = 0;         // index of pdfSampleInds
	uint NsampleSteps;

        std::vector< std::vector<double> > rotEnergyJ;
        std::map<std::string, std::vector<double>* > vibEns;
        std::map<std::string, std::vector<double>* > rotVibBconst;
        std::map<std::string, std::vector<double>* > rotVibDconst;
        std::map<std::string, std::vector<double>* > rotVibHconst;
        std::map<int, std::map<int, std::map<int, std::map<int, double> > > > lgndrInts;

        std::string baseCodeDir;

};


#endif

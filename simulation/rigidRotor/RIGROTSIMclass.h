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
#include <Math/SpecFuncMathMore.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "/reg/neh/home/khegazy/baseScripts/plotClass.h"

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
                                                      
        uint cosOrder;
        bool makePDFs;                                 
        std::string savePDFformat;
        double startPDFTime;		/* ps */
        double endPDFTime;          	/* ps */       
        std::vector<double> sampleTimes;
                                                       
        int MAXj;                                      
	bool hasQuarterRev; 
        double indexRefr;                              
        string vibKey;

        std::string fileName;
        std::string outputDir;

        std::vector<double> jPopulation;
        std::vector<double> vibThermalDist;	
        std::vector< std::vector<double> > rotThermalDist;	
        std::map< int, std::map<int, std::map<int, 
          Eigen::SparseVector< complex<double> > > > > eigenAmps;      // Eigenvectors of |VJm>
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSq;
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSqExpDiag;
        std::map<int, Eigen::SparseMatrix< complex<double> > > cosSqTran;

        RIGROTSIMclass();
        std::vector< std::vector<double> > runSimulation();



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
        std::vector< std::vector<double> > anglePDFCalc();
	double laserPotential(double time);
	double timeInt(double fTime, double iTime);
	void sampleSim(std::string fileName);

        std::vector<uint> pdfSampleInds;
        std::vector<double> pulseTimes;
        std::vector< std::vector<double> > cosEVals;

  	uint ist = 0;         // index of sampleTimes
  	uint ipt = 0;         // index of pulseTimes
  	uint ipi = 0;         // index of pdfSampleInds
	uint NsampleSteps;

        std::map<std::string, std::vector<double>* > vibEns;
        std::map<std::string, std::vector<double>* > rotVibBconst;
        std::map<std::string, std::vector<double>* > rotVibDconst;
        std::map<std::string, std::vector<double>* > rotVibHconst;

};


#endif

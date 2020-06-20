#include "RIGROTSIMclass.h"

RIGROTSIMclass::RIGROTSIMclass() {

  // Vibration variables
  doVib = false;
  vibKey = "NULL";
  NvibStates = 1;     // Count ground state
  vibPopThresh = 5e-3;
  vibEns.emplace("n2o", &n2o_VibEns);
  rotVibBconst.emplace("n2o", &n2o_rotVibBconst);
  rotVibDconst.emplace("n2o", &n2o_rotVibDconst);
  rotVibHconst.emplace("n2o", &n2o_rotVibHconst);

  savePDFformat = "binary";

  cosOrder = 1;
  cosSqEVals.resize(1);
  NYlmEVs = 0;

  baseCodeDir = "/reg/neh/home/khegazy/baseTools/simulation/rigidRotor/";
}


void RIGROTSIMclass::checkAmps() {

  int mSize = rotThermalDist[0].size();

  for (int iv=0; iv<NvibStates; iv++) {
    cout<<"\nnorm vib="<<iv<<":  ";
    for (int im=0; im<=mSize; im++) {
      cout<<"\nnorm m="<<im<<":  ";
      for (int j=0; j<=mSize-im; j++) {
        cout<<"j="<<j<<": "<<eigenAmps[iv][j][im].squaredNorm()<<"   ";
      }
    }
  }
}


vector< complex<double> > RIGROTSIMclass::EVcos() {

  std::vector< complex<double> > ExpVals(cosOrder);
  for (uint i=0; i<ExpVals.size(); i++) {
    ExpVals[i] = std::complex<double>(0,0);
  }
  int mSize = cosSq.size();
  int m = 0;

  // Looping over m
  for (auto& cosSqItr : cosSq) {
    m = cosSqItr.first;
    
    // Looping over cos power
    for (uint ip=0; ip<cosOrder; ip++) {
      
      Eigen::SparseMatrix< complex<double> > cosMat = cosSqItr.second;
      for (uint j=0; j<ip; j++) {
        cosMat = cosMat*cosSqItr.second;
      }

      //if (ip==5) {
        //cout<<"No loop"<<endl<<cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second<<endl;
        //cout<<"\n\n\n\n"<<"Loop"<<cosMat<<endl;
        //cout<<"\n\n\n\n"<<"Diff"<<cosMat - cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second<<endl;
      //}
      // Looping over J for each m in cosSq
      for (int j=0; j<mSize-m; j++) {
        // Looping over vibrational states
        for (int iv=0; iv<NvibStates; iv++) {
          // if/else accounts for 2j+1 multiplicity (only one m=0 state)
          if (cosSqItr.first) {
            ExpVals[ip] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosMat*eigenAmps[iv][j][m]);
            /*
            ExpVals[0] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second*eigenAmps[iv][j][m]);
            
            ExpVals[1] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[2] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[3] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[4] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[5] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
                        */
          }
          else {
            ExpVals[ip] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosMat*eigenAmps[iv][j][m]);
            /*
            ExpVals[0] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[1] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[2] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[3] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
            ExpVals[4] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*cosSqItr.second*eigenAmps[iv][j][m]);          
            ExpVals[5] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*cosSqItr.second
                        *cosSqItr.second*eigenAmps[iv][j][m]);
                        */
          }
        }
      }
    }
  }

  return ExpVals;
}


vector< complex<double> > RIGROTSIMclass::EVsinCos() {

  cout<<"start EVsincCos\n";
  std::vector< complex<double> > ExpVals(cosOrder);
  for (uint i=0; i<ExpVals.size(); i++) {
    ExpVals[i] = std::complex<double>(0,0);
  }

  // Looping over m
  //for (auto& sinCosItr : sinCos) {
  for (int im=0; im<=MAXj; im++) {
    if (sinCos.find(im) == sinCos.end()) {
      std::cerr << "ERROR: sinCos does not have matrix for m=" 
        << im << std::endl;
      exit(1);
    }
    if (sinCos.find(-1*im) == sinCos.end()) {
      std::cerr << "ERROR: sinCos does not have matrix for m=" 
        << -1*im << std::endl;
      exit(1);
    }

    auto sinCosMatP = sinCos[im];
    auto sinCosMatM = sinCos[-1*im];

    
    // Looping over sincos power
    for (uint ip=0; ip<sinCosOrder; ip++) {
      
      for (uint j=0; j<ip; j++) {
        sinCosMatP = sinCosMatP*sinCosMatP;
        sinCosMatM = sinCosMatM*sinCosMatM;
      }

      //if (ip==5) {
        //cout<<"No loop"<<endl<<cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second<<endl;
        //cout<<"\n\n\n\n"<<"Loop"<<cosMat<<endl;
        //cout<<"\n\n\n\n"<<"Diff"<<cosMat - cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second*cosSqItr.second<<endl;
      //}
      // Looping over J for each m in cosSq
      for (int j=0; j<=MAXj-im; j++) {
        // Looping over vibrational states
        for (int iv=0; iv<NvibStates; iv++) {
          // if/else accounts for 2j+1 multiplicity (only one m=0 state)
            cout<<"EV: "<<im<<"  "<<j<<"  "<<iv<<"  "<<endl;
          if (im) {
            //cout<<"mat: "<<sinCosMatM.rows()<<"  "<<sinCosMatM.cols()<<endl;
            //cout<<"vec: "<<eigenAmps[iv][j][im].rows()<<"  "<<eigenAmps[iv][j][im].cols()<<endl;
            if (im==1) {
            cout<<"vec: "<<eigenAmps[iv][j][im]<<endl;
            cout<<"P: "<<rotThermalDist[iv][j+im]
                        *eigenAmps[iv][j][im].dot(sinCosMatP*eigenAmps[iv][j][im])<<endl;
            cout<<"M: "<<rotThermalDist[iv][j+im]
                        *eigenAmps[iv][j][im].dot(sinCosMatM*eigenAmps[iv][j][im])<<endl;
            }
            ExpVals[ip] += 2*rotThermalDist[iv][j+im]
                        *eigenAmps[iv][j][im].dot(sinCosMatP*eigenAmps[iv][j][im]);
          }
          else {
            ExpVals[ip] += rotThermalDist[iv][j+im]
                        *eigenAmps[iv][j][im].dot(sinCosMatP*eigenAmps[iv][j][im]);
          }
        }
      }
    }
  }

  cout<<"end EVsinCos\n";
  return ExpVals;
}


vector< complex<double> > RIGROTSIMclass::EVYlm() {

  int mSize = cosSq.size();
  std::vector< complex<double> > ExpVals(NYlmEVs+1);

  for (int J_=0; J_<=NYlmEVs; J_++) {
    ExpVals[J_] = std::complex<double>(0,0);
    // Looping over vibrational states
    for (int iv=0; iv<NvibStates; iv++) {
      // Looping over m
      for (int m=0; m<=MAXj; m++) {
        // Looping over J for each m in cosSq
        for (int j=0; j<mSize-m; j++) {
          // if/else accounts for 2j+1 multiplicity (only one m=0 state)
          if (m) {
            ExpVals[J_] += 2*rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(Ylm[J_][m]*eigenAmps[iv][j][m]);
          }
          else {
            ExpVals[J_] += rotThermalDist[iv][j+m]
                        *eigenAmps[iv][j][m].dot(Ylm[J_][m]*eigenAmps[iv][j][m]);
          }
        }
      }
    }
  }

  return ExpVals;
}



std::vector< std::vector<double> > RIGROTSIMclass::anglePDFCalc() {

// Coordinates are in molecular frame, where 
//      Z is along laser polarization
//      THETA is the polar angle from the laser polarization [0,PI]
//      PHI is the azimuthal angle about the laser polarization [0,2*PI)
//      sin(theta)*(PI/Nitr)*(2*PI/Nitr) in probability summation is the Jacobian

  cout << "\nINFO: Making angular PDF!!!\n";

  int iphi, itheta;
  int Nitr = 200;
  double sum = 0;
  double theta, phi, prob;
  double angItr = PI/((double)Nitr);
  double binJacobian = (2*PI*PI)/(Nitr*Nitr);

  //complex<double>* Ylm = new complex<double>[MAXj*MAXj*Nitr*Nitr];
  complex<double> wavefxn0, wavefxnM, wavefxnNM;
  std::vector< std::vector<double> > pdf(Nitr);
  for (uint i=0; i<pdf.size(); i++) {
    pdf[i].resize(Nitr, 0);
  }

  // Open or make Ylm references
  FILE* file;
  int im, ij;
  std::vector< std::vector< std::vector< std::vector< complex<double> > > > > YlmVec(MAXj+1);
  for (ij=0; ij<=MAXj; ij++) {
    int index = 0;
    //string fName = "/reg/d/psdm/amo/amoi0314/scratch/YlmRef/YlmRef_J-" 
    //string fName = "/reg/ued/ana/scratch/n2o/YlmRef/YlmRef_J-" 
    string fName = "/reg/neh/home/khegazy/baseTools/simulation/rigidRotor/YlmRef/YlmRef_J-" 
      + to_string(ij) + "_Nitr-" + to_string(Nitr) + ".dat";

    // Create file of Ylm if it doesn't exist
    if (!tools::fileExists(fName)) {
      cout << "INFO: Creating new Ylm reference file " + fName + "!!!\n";
      size_t nwritten = 0;
      vector< complex<double> > YlmJM((ij+1)*Nitr*Nitr);

      file = fopen(fName.c_str(), "wb");

      index = 0;
      for (im=0; im<=ij; im++) {
        for (itheta=0; itheta<Nitr; itheta++) {
          theta = itheta*angItr + angItr/2.0;
          for (iphi=0; iphi<Nitr; iphi++) {
            phi = iphi*2*angItr + angItr;
            YlmJM[index] = (boost::math::spherical_harmonic<double, double>(ij, im, theta, phi));
            index++;
          }
        }
      }

      nwritten += fwrite(&YlmJM[0], sizeof(complex<double>), YlmJM.size(), file);
      fclose(file);
    }

    // Load Ylm reference
    complex<double>* YlmArr = new complex<double>[(ij+1)*Nitr*Nitr];
    file = fopen(fName.c_str(), "rb");
    fread(YlmArr, sizeof(complex<double>), (ij+1)*Nitr*Nitr, file);
    fclose(file); 

    index = 0;
    YlmVec[ij].resize(ij+1);
    for (im=0; im<=ij; im++) {
      YlmVec[ij][im].resize(Nitr);
      for (itheta=0; itheta<Nitr; itheta++) {
        YlmVec[ij][im][itheta].resize(Nitr);
        for (iphi=0; iphi<Nitr; iphi++) {
          YlmVec[ij][im][itheta][iphi] = YlmArr[index];
          index++;
        }
      }
    }
    delete[] YlmArr;

  }



/*
int ij, im;
    for (itheta=0; itheta<Nitr; itheta++) {
      theta = itheta*angItr + angItr/2.0;
      for (iphi=0; iphi<Nitr; iphi++) {
        phi = iphi*2*angItr + angItr;

        for (ij=0; ij<MAXj; ij++) {
          for (im=0; im<MAXj; im++) {
	    if (boost::math::spherical_harmonic<double, double>(ij, im, theta, phi) != Ylm[itheta*(Nitr*MAXj*MAXj) + iphi*(MAXj*MAXj) + ij*MAXj + im]) cout << "Incorrec vals: "<<itheta<<"  "<<iphi<<"  "<<ij<<"  "<<im<<endl;
	  }
	}
}}
*/

  for (itheta=0; itheta<Nitr; itheta++) {
    //theta = itheta*angItr + angItr/2.0;
    for (iphi=0; iphi<Nitr; iphi++) {
      //phi = iphi*2*angItr + angItr;
      prob = 0;

      //YlmTP = Ylm + itheta*Titr + iphi*Pitr;
      for (auto& ampVItr : eigenAmps) {
        for (auto& ampJItr : ampVItr.second) {
          for (auto& ampMItr : ampJItr.second) {
            wavefxn0 = wavefxnM = wavefxnNM = complex<double>(0,0);
  
            // For m!=0
            if (ampMItr.first) {
              for (Eigen::SparseVector< complex<double> >::InnerIterator JinnerItr(ampMItr.second); JinnerItr; ++JinnerItr) {
                if (ampMItr.first > JinnerItr.index()) {
                  continue;
                }
            //cout<<"size: "<<YlmVec[JinnerItr.index()].size()<<endl;
            //cout<<"size: "<<YlmVec[JinnerItr.index()][0].size()<<endl;
            //cout<<"size: "<<YlmVec[JinnerItr.index()][0][0].size()<<endl;
            //cout<<"fill wavefxn: "<<JinnerItr.index()<<"  "<<ampMItr.first<<endl;
            //cout<<"val: "<<JinnerItr.value()<<endl;
            //cout<<YlmVec[JinnerItr.index()][ampMItr.first][itheta][iphi]<<endl;
                wavefxnM += JinnerItr.value()*YlmVec[JinnerItr.index()][ampMItr.first][itheta][iphi];
                //cout<<"wm: "<<wavefxnM<<endl;
                wavefxnNM += conj(JinnerItr.value()*YlmVec[JinnerItr.index()][ampMItr.first][itheta][iphi]);
                //cout<<"wmn: "<<wavefxnM<<endl;
              }

              //prob += thermalDist[ampJItr.first]*(std::norm(wavefxnM) + std::norm(wavefxnNM))*sin(theta);  //For pdf to be norm
              prob += rotThermalDist[ampVItr.first][ampJItr.first]*(std::norm(wavefxnM) + std::norm(wavefxnNM));
            }
            // For m==0
            else {
              for (Eigen::SparseVector< complex<double> >::InnerIterator JinnerItr(ampMItr.second); JinnerItr; ++JinnerItr) {
                if (ampMItr.first > JinnerItr.index()) {
                  continue;
                }
                //if (JinnerItr.index() != 2) {
                //  continue; //DELETE ME
                //}
            //cout<<"fill wavefxn: "<<JinnerItr.index()<<"  "<<ampMItr.first<<endl;
                wavefxn0 += JinnerItr.value()*YlmVec[JinnerItr.index()][0][itheta][iphi];
              }

              //prob += thermalDist[ampJItr.first]*std::norm(wavefxn0)*sin(theta);     //For pdf to be norm
              prob += rotThermalDist[ampVItr.first][ampJItr.first]*std::norm(wavefxn0);
            }


          }
        }
        
      }

      prob *= binJacobian;
      pdf[itheta][iphi] = prob;
      sum += prob;

    }
  }

  
  // Normalize the pdf
  for (itheta=0; itheta<(int)pdf.size(); itheta++) {
    for (iphi=0; iphi<(int)pdf[itheta].size(); iphi++) {
      pdf[itheta][iphi] /= sum;
    }
  }

  return pdf;
}


double RIGROTSIMclass::laserPotential(double time) {

  return -pow(E0*sin(PI*time/(pulseTlengthAU)), 2)*deltaAlphaAU/4;
}


double RIGROTSIMclass::timeInt(double fTime, double iTime) {

  // Integral of sin^2
  double A = PI/pulseTlengthAU;
  return (fTime-iTime)/2.0 - (sin(2*fTime*A) - sin(2*iTime*A))/(4*A);
}


void RIGROTSIMclass::sampleSim(double time_fs) {

      // Save PDF of the angular distribution
      if (ist == pdfSampleInds[ipi]) {
        std::string filePrefix = PDFoutputDir + PDFfileName;
        cout<<"FNAME: " << filePrefix 
                    + "_Time-" + to_string(pdfSampleTimes[ipi]).substr(0,7)<<endl;
        std::vector< std::vector<double> > pdf = anglePDFCalc();
        if (savePDFformat == "ROOT") {
          std::string fileName = filePrefix 
                    + "_Time-" + to_string(pdfSampleTimes[ipi]).substr(0,7)
                    + "_ThetaBins-" + to_string(pdf.size()) 
                    + "_PhiBins-" + to_string(pdf[0].size()) 
                    + ".root";
          cout<<"FILENAME: "<<fileName<<endl;
          TFile* pdfFile;
          pdfFile = TFile::Open(fileName.c_str(), "RECREATE");
          TH2D* pdfHist = new TH2D("thetaPhiPDF", "thetaPhiPDF", 
                    pdf.size(), 0, PI, pdf[0].size(), 0, 2*PI);
          for (uint itheta=0; itheta<pdf.size(); itheta++) {
            for (uint iphi=0; iphi<pdf[itheta].size(); iphi++) {
              pdfHist->SetBinContent(itheta+1, iphi+1, pdf[itheta][iphi]);
            }
          }
          pdfHist->Write();
          pdfFile->Close();
	  cout << "INFO: Angular distribution saved in " << 
		filePrefix + ".root!!!" << endl;
        }
        //else if (savePDFformat == "binary") {
        FILE* binFile = fopen((filePrefix 
                  + "_Time-" + to_string(pdfSampleTimes[ipi]).substr(0,7)
                  + "_ThetaBins-" + to_string(pdf.size()) 
                  + "_PhiBins-" + to_string(pdf[0].size()) 
                  + ".dat").c_str(), "wb");
        for (uint itheta=0; itheta<pdf.size(); itheta++) {
          fwrite(&pdf[itheta][0], sizeof(double), pdf[itheta].size(), binFile);
        }
        fclose(binFile);
	cout << "INFO: Angular distribution saved in " << 
	    filePrefix + ".dat!!!" << endl;
        //}

	ipi++;
      }

      // Calculate expectation values
      if (doCosSqEV) {
        std::vector< std::complex<double> > expVals = EVcos();
        for (uint i=0; i<cosOrder; i++) {
          cosSqEVals[i][ist] = expVals[i].real();
        }
      }
      if (doSinCosEV) {
        std::vector< std::complex<double> > expVals = EVsinCos();
        for (uint i=0; i<sinCosOrder; i++) {
          sinCosEVals[i][ist] = expVals[i].real();
        }
      }

      // Calculate wavefunction coefficients
      if (NYlmEVs) {
        std::vector< std::complex<double> > expVals = EVYlm();
        for (int J_=0; J_<=NYlmEVs; J_++) {
          YlmEVals[J_][ist] = expVals[J_].real();
        }
      }

      // Increment to the next sample time
      ist++;
	
      return;
}


std::map<int, std::vector< complex<double> > > RIGROTSIMclass::projectDensitySH() {

  std::map<int, std::vector< complex<double> > > coeffs;
  std::map<int, std::vector< complex<double> > > waveFxn;
  int jSize = rotThermalDist[0].size();
  complex<double> init(0,0);


  for (int m=0; m<=jSize; m++) {
    for (int j=m; j<=jSize; j++) {
      if (waveFxn[j].size() == 0) {
        waveFxn[j].resize(j+1, init);
      }
      //waveFxn[j][m] = init;
      for (int jj=m; jj<=jSize; jj++) {
        for (int iv=0; iv<NvibStates; iv++) {
          waveFxn[j][m] += rotThermalDist[iv][jj]*eigenAmps[iv][jj][m].coeffRef(j);
        }
      }
    }
  }

  /*
  for (int J=0; J<=jSize; J++) {
    coeffs[J].resize(2*J+1, init);
    for (int M=-1*J; M<=J; M++) {

      for (int j1=0; j1<=jSize; j1++) {
        for (int m1_=-1*j1; m1_<=j1; m1_++) {
          m1 = -1*m1_;
          m2 = m1 + M;
          j2Max = jSize < j1 + J ? jSize : j1 + J;
          for (int j2=abs(m2); j2<=j2Max; j2++) {

            coeffs[J][M] += conj(waveFxn[j1][m1])*waveFxn[j2][m2]
                *sqrt((2*j1 + 1)*(2*J + 1)*(2*j2 + 1)/(4*PI))
                *ROOT::Math::wigner_3j(1,2,4,5,6,7);
          }
        }
      }
    }
  }
  */

  return coeffs;
}



std::vector<double> RIGROTSIMclass::getPopulationDistribution() {

  int J;
  std::vector<double> population(MAXj, 0);
  for (int im=0; im<=MAXj; im++) {
    for (int j=0; j<=MAXj-im; j++) {
      for (int iv=0; iv<NvibStates; iv++) {
        for (int jj=0; jj<MAXj-im; jj++) {
          J = im + jj;
          if (im) {
            population[J] += std::real(rotThermalDist[iv][j+im]
                *eigenAmps[iv][j][im].coeffRef(jj)*std::conj(eigenAmps[iv][j][im].coeffRef(jj)));
            population[J] += std::real(rotThermalDist[iv][j+im]
                *eigenAmps[iv][j][im].coeffRef(jj)*std::conj(eigenAmps[iv][j][im].coeffRef(jj)));
          }
          else {
            population[J] += std::real(rotThermalDist[iv][j+im]
                *eigenAmps[iv][j][im].coeffRef(jj)*std::conj(eigenAmps[iv][j][im].coeffRef(jj)));
          }
        }
      }
    }
  }

  return population;
}


void RIGROTSIMclass::calculateThermalDistribution() {

  double Zvib = 0;
  if (doVib) {

    // Only considering most important vibrational states
    double bwVal;
    vibThermalDist.clear();
    for (int iv=0; iv<(int)(*vibEns[vibKey]).size(); iv++) {
      Zvib += exp(-(*vibEns[vibKey])[iv]*icm_to_au/(KBOLTZ_AU*temperature));
    }
    for (int iv=0; iv<(int)(*vibEns[vibKey]).size(); iv++) {
      bwVal = exp(-(*vibEns[vibKey])[iv]*icm_to_au/(KBOLTZ_AU*temperature))/Zvib;
      if (bwVal > vibPopThresh) {
        vibThermalDist.push_back(bwVal);
      }
      else {
        break;
      }
    }
    NvibStates = (int)vibThermalDist.size();

    // Rotational energies 
    rotEnergyJ.resize(NvibStates);
    for (int iv=0; iv<NvibStates; iv++) {
      rotEnergyJ[iv].resize(MAXj+1);
      for (int j=0; j<=MAXj; j++) {
        rotEnergyJ[iv][j] = (*rotVibBconst[vibKey])[iv]*j*(j + 1)*icm_to_au
                      - (*rotVibDconst[vibKey])[iv]*pow(j*(j + 1), 2)*icm_to_au
                      + (*rotVibHconst[vibKey])[iv]*pow(j*(j + 1), 3)*icm_to_au;
      }
    }
  }
  else {
    cout<<"in not dovib"<<endl;
    rotEnergyJ.resize(1);
    rotEnergyJ[0].resize(MAXj+1);
    for (int j=0; j<=MAXj; j++) {
      rotEnergyJ[0][j] = rotConstB_AU*j*(j + 1)
                      - rotConstD_AU*pow(j*(j + 1), 2)
                      + rotConstH_AU*pow(j*(j + 1), 3);
    }
  }
  cout<<"NvibStates: "<<NvibStates<<endl;

  complex<double> complexNumber;
  complexNumber.imag(0);
  complexNumber.real(1);
  for (int iv=0; iv<NvibStates; iv++) {
    for (int im=0; im<=MAXj; im++) {
      for (int j=0; j<=MAXj-im; j++) {
        eigenAmps[iv][j][im].resize(MAXj+1-im);
        eigenAmps[iv][j][im].insert(j) = complexNumber;
        //if (im <= j) eigenAmps[j][im].insert(j) = complexNumber;
      }
    }
  }

  cout<<"partition"<<endl;
  // Partition function
  double Zrot = 0;
  double qrtrRv = 1;
  rotThermalDist.resize(NvibStates);
  if (doVib) {
    for (int iv=0; iv<NvibStates; iv++) {
      rotThermalDist[iv].resize(MAXj+1, 0);
      for (int j=0; j<=MAXj; j++) {
        if (hasQuarterRev) {
          qrtrRv = ((j + 1)%2) + 1;
        }
        Zrot += qrtrRv*(2*j + 1)
          *exp(-(rotEnergyJ[iv][j] + (*vibEns[vibKey])[iv]*icm_to_au)
              /(KBOLTZ_AU*temperature));
      }
    }
    for (int iv=0; iv<NvibStates; iv++) {
      for (int j=0; j<=MAXj; j++) {
        if (hasQuarterRev) {
          qrtrRv = ((j + 1)%2) + 1;
        }
        rotThermalDist[iv][j] = qrtrRv
          *exp(-(rotEnergyJ[iv][j] + (*vibEns[vibKey])[iv]*icm_to_au)
              /(KBOLTZ_AU*temperature))/Zrot;
      }
    }
  }
  else {
    rotThermalDist[0].resize(MAXj+1, 0);
    for (int j=0; j<=MAXj; j++) {
      if (hasQuarterRev) {
        qrtrRv = ((j + 1)%2) + 1;
      }
      Zrot += qrtrRv*(2*j + 1)
        *exp(-rotEnergyJ[0][j]/(KBOLTZ_AU*temperature));
    }
    for (int j=0; j<=MAXj; j++) {
      if (hasQuarterRev) {
        qrtrRv = ((j + 1)%2) + 1;
      }
      rotThermalDist[0][j] = qrtrRv
        *exp(-rotEnergyJ[0][j]/(KBOLTZ_AU*temperature))/Zrot;
    }
  }
}



void RIGROTSIMclass::calculateCosSqMatrices() {
  // Setup CosSq matrix and tranform matrix for time evolution in pulse
  for (int im=0; im<=MAXj; im++) {
    cosSq[im].resize(MAXj+1-im, MAXj+1-im);
    cosSqExpDiag[im].resize(MAXj+1-im, MAXj+1-im);
    cosSqTran[im].resize(MAXj+1-im, MAXj+1-im);
    cosSqEigVal[im].resize(MAXj+1-im, MAXj+1-im);

    double jj;
    complex<double> complexNumber;
    complexNumber.imag(0);
    for (int j=0; j<=MAXj-im; j++) {
      complexNumber.real(0);
      cosSqExpDiag[im].insert(j, j) = complexNumber;

      jj = j + im;
      //if (im > j) continue;

      // <J,M| cos^2 |J+2,M>    
      if ((j <= MAXj-2-im)) {
        complexNumber.real((sqrt((double)((2*jj + 1)*(2*jj + 5)*(jj+1 - im)))
              /((double)((2*jj + 1)*(2*jj + 3)*(2*jj + 5))))
              *sqrt((double)((jj + 2 - im)*(jj + 1 + im)*(jj + 2 + im))));
        cosSq[im].insert(j, j+2) = complexNumber;
      }

      // <J,M| cos^2 |J,M>      
      complexNumber.real(1.0/3.0 + 2.0/3.0*(((double)(jj*(jj + 1) - 3*im*im))
              /((double)((2*jj + 3)*(2*jj - 1)))));
      cosSq[im].insert(j, j) = complexNumber;

      // <J,M| cos^2 |J-2,M>    
      if (j >1) {
        complexNumber.real((sqrt((double)((2*jj - 3)*(2*jj + 1)*(jj - 1 - im)))
              /((double)((2*jj - 3)*(2*jj - 1)*(2*jj + 1))))
              *sqrt(((double)((jj - im)*(jj - 1 + im)*(jj + im)))));
        cosSq[im].insert(j, j-2) = complexNumber;
      }
    }

    // Diagonalizing CosSq matrix and defining transfer matrix
    Eigen::MatrixXcd cosSqDense = Eigen::MatrixXcd(cosSq[im]);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigSolver(cosSqDense);

    for (int jr=0; jr<=MAXj-im; jr++) {
      cosSqEigVal[im][jr] = eigSolver.eigenvalues()[jr];
      //coeffCSD[jr] = T(jr, jr, exp(eigSolver.eigenvalues()[jr]));
      for (int jc=0; jc<=MAXj-im; jc++) {
        cosSqTran[im].insert(jr, jc) = eigSolver.eigenvectors()(jr,jc);
      }
    }

    cosSqTran[im].prune(prunefx);
    cosSq[im].prune(prunefx);
  }
}


void RIGROTSIMclass::calculateSinCosMatrices() {
  // Creating sin*cos matrix
  if (doSinCosEV) {
    int im = 0;
    for (int im=0; im<=MAXj; im++) {
      sinCos[im].resize(MAXj+1-im, MAXj+1-im);
      sinCos[im].setZero();
      sinCos[-1*im].resize(MAXj+1-im, MAXj+1-im);
      sinCos[-1*im].setZero();
    }

    int ir, ic;
    std::string line, fileName;
    for (int ij=1; ij<=MAXj; ij++) {
      for (int ijj=ij; ijj<=MAXj; ijj++) {
        if ((ijj - ij) % 2 == 1) {
          fileName = 
            "/reg/d/psdm/amo/amoi0314/scratch/lgndrInts/mANDmp1/assLegendreIntegral_l["
            + to_string(ij) + "," + to_string(ijj) + "].txt";
         
          ifstream file(fileName.c_str());
          if (file.is_open()) {
            im = -1*ij;
            while (getline(file, line)) {
              //cout<<"m/val: "<<im<<"  "<<ij<<"/"<<ijj<<"  "<<sinCosSq[im].rows()<<"/"<<sinCosSq[im].cols()<<"   "+line<<endl;
              lgndrInts[ij][im][ijj][im+1] = std::stod(line);
              im++;
            }
          }
          else {
            std::cerr << "ERROR: Cannot find file " << fileName << std::endl;
          }
          file.close();
          
          fileName = 
            "/reg/d/psdm/amo/amoi0314/scratch/lgndrInts/mp1ANDm/assLegendreIntegral_l["
            + to_string(ij) + "," + to_string(ijj) + "].txt";
         
          file.open(fileName.c_str(), std::ifstream::in);
          if (file.is_open()) {
            im = -1*ij;
            while (getline(file, line)) {
              ir = ij - fabs(im);
              ic = ijj - fabs(im);
              lgndrInts[ij][im][ijj][im-1] = std::stod(line);
              im++;
            }
          }
          else {
            std::cerr << "ERROR: Cannot find file " << fileName << std::endl;
          }
          file.close();

        }
      }
    }

    cout<<"starting to fill sincos"<<endl;
    int mRange;
    double lgInt1, lgInt2, lgInt3, lgInt4;
    double cJpM, cJmM, cJpMp, cJmMp;
    complex<double> complexNumber;
    complexNumber.real(0);
    complexNumber.imag(0);
    for (int ij=0; ij<=MAXj; ij++) {
      for (int ijj=0; ijj<=MAXj; ijj++) {
        if ((ijj - ij) % 2 == 1) {
          if (ij > ijj) {
            mRange = ijj;
          }
          else {
            mRange = ij;
          }
cout<<"J: "<<ij<<" / "<<ijj<<endl;

          // Filling sinCos matrix where m index is for m and m+1 ME
          // <ij,im| sinCos |ijj,im>
          for (int im=-1*mRange; im<=mRange; im++) {
            cout<<"m: "<<im<<endl;
            cout<<"111"<<endl;
            cJpM  = 0;
            cJpMp = 0;
            cJmM  = 0;
            cJmMp = 0;
            if (ij + 1 >= fabs(im)) {
              cJpM = std::sqrt(std::pow(ij-im+1, 2)*boost::math::factorial<double>(ij-im)
                  /(4*PI*(2*ij+1)*boost::math::factorial<double>(ij+im)));
            }
            cout<<"111"<<endl;
            if (ijj + 1 >= fabs(im) + 1) {
              cJpMp = std::sqrt(boost::math::factorial<double>(ijj-im)
                  /(4*PI*(2*ijj+1)*boost::math::factorial<double>(ijj+im)));
            }
            cout<<"111"<<endl;
            if ((ij > 0) && (ij - 1 >= fabs(im))) {
              cJmM = std::sqrt((ij+im)*boost::math::factorial<double>(ij-im)
                  /(4*PI*(2*ij+1)*boost::math::factorial<double>(ij+im-1)));
            }
            cout<<"111"<<endl;
            if ((ijj > 0) && (ijj - 1 > fabs(im) + 1)) {
              cJmMp = std::sqrt((ijj-im)*(ijj-im-1)*boost::math::factorial<double>(ijj-im-2)
                  /(4*PI*(2*ijj+1)*boost::math::factorial<double>(ijj+im)));
            }


            // WHAT DOES IT MEAN IF IJ-1 OR IJJ-1 < 0? ARE THE INTEGRALS 0? LOOK AT RECURSIVE RELATIIONS

            cout<<"starting integrals"<<endl;
            lgInt1 = 0;
            lgInt2 = 0;
            lgInt3 = 0;
            lgInt4 = 0;

            if ((ij + 1 <= MAXj) && (ijj + 1 <= MAXj)) {
              if (lgndrInts.find(ij+1) == lgndrInts.end()) {
                cout<<"Missing: "<<ij+1<<endl;
              }
              else {
                if (lgndrInts[ij+1].find(im) == lgndrInts[ij+1].end()) {
                  cout<<"Missing: "<<ij+1<<"  "<<im<<endl;
                }
                else {
                  if (lgndrInts[ij+1][im].find(ijj+1) == lgndrInts[ij+1][im].end()) {
                    cout<<"Missing: "<<ij+1<<"  "<<im<<"  "<<ijj+1<<endl;
                  }
                  else {
                    if (lgndrInts[ij+1][im][ijj+1].find(im+1) == lgndrInts[ij+1][im][ijj+1].end()) {
                      cout<<"Missing: "<<ij+1<<"  "<<im<<"  "<<ijj+1<<"  "<<im+1<<endl;
                    }
                  }
                }
              }

              lgInt1 = lgndrInts[ij+1][im][ijj+1][im+1];
                cout<<lgInt1<<endl;
            }
             
            if ((ij + 1 <= MAXj) && (ijj - 1 >= 0)) {
              if (lgndrInts.find(ij+1) == lgndrInts.end()) {
                cout<<"Missing: "<<ij+1<<endl;
              }
              else {
                if (lgndrInts[ij+1].find(im) == lgndrInts[ij+1].end()) {
                  cout<<"Missing: "<<ij+1<<"  "<<im<<endl;
                }
                else {
                  if (lgndrInts[ij+1][im].find(ijj-1) == lgndrInts[ij+1][im].end()) {
                    cout<<"Missing: "<<ij+1<<"  "<<im<<"  "<<ijj-1<<endl;
                  }
                  else {
                    if (lgndrInts[ij+1][im][ijj-1].find(im+1) == lgndrInts[ij+1][im][ijj-1].end()) {
                      cout<<"Missing: "<<ij+1<<"  "<<im<<"  "<<ijj-1<<"  "<<im+1<<endl;
                    }
                  }
                }
              }

              lgInt2 = lgndrInts[ij+1][im][ijj-1][im+1];
                cout<<lgInt2<<endl;
            }

            if ((ij - 1 >= 0) && (ijj + 1 <= MAXj)) {
              if (lgndrInts.find(ij-1) == lgndrInts.end()) {
                cout<<"Missing: "<<ij-1<<endl;
              }
              else {
                if (lgndrInts[ij-1].find(im) == lgndrInts[ij-1].end()) {
                  cout<<"Missing: "<<ij-1<<"  "<<im<<endl;
                }
                else {
                  if (lgndrInts[ij-1][im].find(ijj+1) == lgndrInts[ij-1][im].end()) {
                    cout<<"Missing: "<<ij-1<<"  "<<im<<"  "<<ijj+1<<endl;
                  }
                  else {
                    if (lgndrInts[ij-1][im][ijj+1].find(im+1) == lgndrInts[ij-1][im][ijj+1].end()) {
                      cout<<"Missing: "<<ij-1<<"  "<<im<<"  "<<ijj+1<<"  "<<im+1<<endl;
                    }
                  }
                }
              }

              lgInt3 = lgndrInts[ij-1][im][ijj+1][im+1];
                cout<<lgInt3<<endl;
            }

            if ((ij - 1 >= 0) && (ijj - 1 >= 0)) {
            if (lgndrInts.find(ij-1) == lgndrInts.end()) {
                cout<<"Missing: "<<ij-1<<endl;
              }
              else {
                if (lgndrInts[ij-1].find(im) == lgndrInts[ij-1].end()) {
                  cout<<"Missing: "<<ij-1<<"  "<<im<<endl;
                }
                else {
                  if (lgndrInts[ij-1][im].find(ijj-1) == lgndrInts[ij-1][im].end()) {
                    cout<<"Missing: "<<ij-1<<"  "<<im<<"  "<<ijj-1<<endl;
                  }
                  else {
                    if (lgndrInts[ij-1][im][ijj-1].find(im+1) == lgndrInts[ij-1][im][ijj-1].end()) {
                      cout<<"Missing: "<<ij-1<<"  "<<im<<"  "<<ijj-1<<"  "<<im+1<<endl;
                    }
                  }
                }
              }
                lgInt4 = lgndrInts[ij-1][im][ijj-1][im+1];
                cout<<lgInt4<<endl;
            }

            complexNumber.real(
                cJpM*cJpMp*lgInt1
                - cJpM*cJmMp*lgInt2
                + cJmM*cJpMp*lgInt3
                - cJmM*cJmMp*lgInt4);

            ir = ij - fabs(im);
            ic = ijj - fabs(im);
            sinCos[im].insert(ir, ic) = complexNumber; 
          }
        }
      }
    }

    // Prune small values to speed up simulation
    double prevPruneVal = prunefx.cutoff;
    prunefx.cutoff = 1e-7;
    for (auto & sinCosItr : sinCos) {
      sinCosItr.second.prune(prunefx);
    }
    prunefx.cutoff = prevPruneVal;
    cout<<"sinCos Mat\n"<<sinCos[-1]<<endl;
  }
}


void RIGROTSIMclass::calculateYlmMatrices() {

  int J;
  Ylm.resize(NYlmEVs+1);
  for (int J_=0; J_<=NYlmEVs; J_++) {
    if (evenOnlyAxisDist) {
      J = 2*J_;
    }
    else {
      J = J_;
    }
    if (J > MAXj) {
      break;
    }

    std::string wigner3J_fileName = baseCodeDir 
        + "YlmME/YlmME_L-" + to_string(J)
        + "_M-" + to_string(0)
        + "_maxJ-" + to_string(MAXj) + ".dat"; 
    if (!tools::fileExists(wigner3J_fileName)) {
      std::cerr << "Cannot find file " + wigner3J_fileName
          + " now exiting!!!\n";
    }

    
    std::string wigner3J_inds_fileName =  
        + "YlmME_indices_L-" + to_string(J)
        + "_M-" + to_string(0)
        + "_maxJ-" + to_string(MAXj); 
    std::vector<int> shape = save::getShape(
        baseCodeDir + "YlmME",
        wigner3J_inds_fileName);
    wigner3J_inds_fileName +=
        "_bins[" + to_string(shape[0]) + "," + to_string(shape[1]) + "].dat";
    wigner3J_inds_fileName = baseCodeDir + "YlmME/" + wigner3J_inds_fileName;

    std::vector<double> wgn3j(shape[0]);
    std::vector< std::vector<int> > indices(shape[0]);
    for (int i=0; i<shape[0]; i++) {
      indices[i].resize(shape[1], 0);
    }

    save::importDat<int32_t>(indices, wigner3J_inds_fileName);
    save::importDat<double>(wgn3j, wigner3J_fileName);


    for (int im=0; im<=MAXj; im++) {
      Ylm[J_][im].resize(MAXj+1-im, MAXj+1-im);
    }

    int l1, l2, m1, m2;
    for (int i=0; i<shape[0]; i++) {
      l1 = indices[i][0];  m1 = abs(indices[i][1]);
      l2 = indices[i][2];  m2 = abs(indices[i][3]);
      Ylm[J_][m1].coeffRef(l1-m1, l2-m1) = wgn3j[i];
    }

  }
}
      
      

void RIGROTSIMclass::runSimulation() {

// Each |jm> eigenstate is propogated independently, when 
//      expectation values or norms are calculated each 
//      state is weighted by the amplitude of that state at 
//      the t=0 given by the thermal distribution

  if (MAXj > MAXj) {
    cerr << "ERROR: Quantum number m cannot be larger than j (MAXj > MAXj)!!!\n\n\n";
    exit(0);
  }
  if (MAXj < MAXj) {
    cerr << "WARNING: Quantum number m does not reach the rull range for all j states (MAXj < MAXj)!!!\n\n\n";
  }
  if (temperature == 0) temperature = 1e-20;
  if (vibKey != "NULL") {
    doVib = true;
    if (vibEns.find(vibKey) == vibEns.end()) {
      cerr << "ERROR: Cannot find distrubution " << vibKey << "!!!\n\n";
      exit(0);
    }
  }


  // Printing values for simulation
  cout << "Start Time:\t\t"       <<  startTime       <<  endl;
  cout << "End Time:  \t\t"       <<  endTime         <<  endl;
  cout << "Temperature:\t\t"      <<  temperature     <<  endl;
  cout << "Laser Power:\t\t"      <<  laserIntensity  <<  endl;
  cout << "Rot Const B:\t\t"      <<  rotConstB       <<  endl;
  cout << "Rot Const D:\t\t"      <<  rotConstD       <<  endl;
  cout << "Rot Const H:\t\t"      <<  rotConstH       <<  endl;
  cout << "deltaAlpha:\t\t"       <<  deltaAlpha      <<  endl;
  cout << "Sample Start Time:\t"  <<  startSampleTime <<  endl;
  cout << "Sample End Time:\t"    <<  endSampleTime   <<  endl;
  cout << "Sample Step:\t\t"      <<  sampleStep      <<  endl;
  cout << "Do PDFs:\t\t"          <<  makePDFs        <<  endl;
  cout << "PDF Start Time:\t\t"   <<  startPDFTime    <<  endl;
  cout << "PDF End Time:\t\t"     <<  endPDFTime      <<  endl;
  cout << "Number Axis Dist:\t\t" <<  NYlmEVs  <<  endl;
  cout << "dTime Evol:\t\t"       <<  dTimeEvolStep   <<  endl;
  cout << "Pulse Length:\t\t"     <<  pulseTlength    <<  endl;
  cout << "Npulses:   \t\t"       <<  Npulses         <<  endl;
  cout << "Max J:    \t\t"        <<  MAXj            <<  endl;
  cout << "Vibration Dist:\t"     <<  vibKey          <<  endl;

  if (!opendir(PDFoutputDir.c_str())) {
    mkdir(PDFoutputDir.c_str(), 0777);
  }


  // Filling time vectors here to avoid double precision errors
  ist = ipt = ipi = 0;
  sampleTimes.clear();
  pulseTimes.clear();
  pdfSampleInds.clear();

  uint k = 0;
  pulseTimes.resize(Npulses);
  for (k=0; k<Npulses; k++) {
    pulseTimes[k] = k*pulseTspacing*1e3*fs_to_au;
  }
  k = 0;
  while (k*sampleStep + startSampleTime <= endSampleTime) {
    sampleTimes.push_back((k*sampleStep + startSampleTime)*1e3*fs_to_au);
    if (makePDFs 
	&& (k*sampleStep + startSampleTime >= startPDFTime*0.99999) 
	&& (k*sampleStep + startSampleTime <= endPDFTime*1.00001)) {
      pdfSampleInds.push_back(k);
      int tm = (int)((k*sampleStep + startSampleTime)*1000) % 10;
      if (tm == 9) {
        pdfSampleTimes.push_back(k*sampleStep + startSampleTime + 1);
      }
      else if (tm == 1) {
        pdfSampleTimes.push_back(k*sampleStep + startSampleTime - 1);
      }
      else {
        pdfSampleTimes.push_back(k*sampleStep + startSampleTime);
      }

    }
    k++;
  }
  if (!pdfSampleInds.size()) {
    pdfSampleInds.push_back(-1);
  }

cout<<"Num sample Points: "<<sampleTimes.size()<<endl;


  // Converting to AU
  startTimeAU = startTime*1e3*fs_to_au;
  endTimeAU = endTime*1e3*fs_to_au;
  pulseTlengthAU = pulseTlength*1e3*fs_to_au;
  pulseTspacingAU = pulseTspacing*1e3*fs_to_au;
  startSampleTimeAU = startSampleTime*1e3*fs_to_au;
  endSampleTimeAU = endSampleTime*1e3*fs_to_au;
  dTimeEvolStepAU = dTimeEvolStep*1e3*fs_to_au;
  rotConstB_AU = rotConstB*icm_to_au;
  rotConstD_AU = rotConstD*icm_to_au;
  rotConstH_AU = rotConstH*icm_to_au;
  deltaAlphaAU = deltaAlpha*pow(m_to_au,3);

  
  /////  Declaring Variables  /////

  clock_t clockBegin, clockEnd;
  clockBegin=clock();

  int j, jr, jc, im, iv;
  complex<double> complexNumber;
  vector<double> population(MAXj);

  eigenAmps.clear();
  cosSq.clear();
  cosSqExpDiag.clear();
  cosSqTran.clear();
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > pulsedTgenL;
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > pulsedTgenR;
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > pulseEvolMat;
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > H0dTevol;
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > H0dTgen;
  map<int, map<int, Eigen::SparseMatrix< complex<double> > > > H0Tevol;

  prunefx.cutoff = 1e-4;



  ////////////////////////////////////////////
  /////  Defining simulation parameters  /////
  ////////////////////////////////////////////

  // Electric field
  bool firstPulse = (Npulses > 1);
  double E0 = sqrt(2*laserIntensity*1.0e4/(C_SI*indexRefr*EPS0_SI))*Vpm_to_au;

  /////  Making thermal distributions  /////
  calculateThermalDistribution();

  /////  Calculate CosSq, Ylm, ... Matrices  /////
  calculateCosSqMatrices();
  calculateSinCosMatrices();
  calculateYlmMatrices();

  // Setup H0 dTimeEvolStep evolution first, for timeEvolStep
  for (iv=0; iv<NvibStates; iv++) {
    for (im=0; im<=MAXj; im++) {
      H0dTevol[iv][im].resize(MAXj+1-im, MAXj+1-im);
      H0dTgen[iv][im].resize(MAXj+1-im, MAXj+1-im);
      H0Tevol[iv][im].resize(MAXj+1-im, MAXj+1-im);
      for (int j=0; j<=MAXj-im; j++) {
        complexNumber.real(0);
        complexNumber.imag(-rotEnergyJ[iv][j+im]*dTimeEvolStepAU/2.0);
        H0dTevol[iv][im].insert(j, j) = exp(complexNumber);
        H0dTgen[iv][im].insert(j, j) = exp(complexNumber);
      }
    }
  }

  // Defining commonly used variables
  complexNumber.real(1);
  complexNumber.imag(0);
  for (int iv=0; iv<NvibStates; iv++) {
    for (int im=0; im<=MAXj; im++) {
      pulsedTgenR[iv][im].resize(MAXj+1-im, MAXj+1-im);
      pulsedTgenL[iv][im].resize(MAXj+1-im, MAXj+1-im);
      pulseEvolMat[iv][im].resize(MAXj+1-im, MAXj+1-im);
      for (int j=0; j<=MAXj-im; j++) {
        pulseEvolMat[iv][im].insert(j, j) = complexNumber;
        //if (im <= j) eigenAmps[j][im].insert(j) = complexNumber;
      }
    }
  }


  for (int im=0; im<=MAXj; im++) {
    for (iv=0; iv<NvibStates; iv++) {
      pulsedTgenL[iv][im] = H0dTgen[iv][im]*cosSqTran[im];
      pulsedTgenR[iv][im] = cosSqTran[im].adjoint()*H0dTgen[iv][im];
      pulsedTgenL[iv][im].prune(prunefx);
      pulsedTgenR[iv][im].prune(prunefx);
    }
  }

  cout<<"start sim"<<endl;
  ////////////////////////////////////
  /////  Simulating starts here  /////
  ////////////////////////////////////

  double curTime = startTimeAU;
  bool doPulse, doSampling;
  bool regTimeStep;
  Eigen::SparseMatrix< complex<double> > timeEvolMat;

  double TnextSample, TnextPulse, endPulseSim;
  double timeEvolStep, dTimePulseEvol;
  NsampleSteps = sampleTimes.size(); 

  if (doCosSqEV) {
    cosSqEVals.resize(cosOrder);
    for (uint i=0; i<cosOrder; i++) {
      cosSqEVals[i].resize(NsampleSteps, 0);
    }
  }
  if (doSinCosEV) {
    sinCosEVals.resize(cosOrder);
    for (uint i=0; i<sinCosOrder; i++) {
      sinCosEVals[i].resize(NsampleSteps, 0);
    }
  }
  if (NYlmEVs) {
    NYlmEVs = min(NYlmEVs, MAXj);
    YlmEVals.resize(NYlmEVs + 1);
    for (uint j=0; j<YlmEVals.size(); j++) {
      YlmEVals[j].resize(NsampleSteps);
    }
  }

  // Get initial distribution
  population0 = getPopulationDistribution();

  while ((curTime <= endTimeAU) && (ist != NsampleSteps)) {

    /////  Decide how to evolve through time  /////
    TnextSample = sampleTimes[ist] - curTime;
    TnextPulse = pulseTimes[ipt] - curTime;
    if (((TnextSample < TnextPulse) || ipt == pulseTimes.size())
                && (ist != sampleTimes.size())) {
      timeEvolStep = TnextSample;
      doSampling = true;
      doPulse = false;
    }
    else if (ipt != pulseTimes.size()) {
      timeEvolStep = TnextPulse;
      doPulse = true;
      doSampling = false;
    }
    else {
      break;
    }

    /////  Time Evolution Through Free Space  /////
    for (iv=0; iv<NvibStates; iv++) {
      for (im=0; im<=MAXj; im++) {
        for (j=0; j<=MAXj-im; j++) {
         complexNumber.real(0);
         complexNumber.imag(-rotEnergyJ[iv][j+im]*timeEvolStep);
         H0Tevol[iv][im].coeffRef(j, j) = exp(complexNumber);
        }
        for (j=0; j<=MAXj-im; j++) {
         eigenAmps[iv][j][im] = H0Tevol[iv][im]*eigenAmps[iv][j][im];
        }
      }
    }

    // Update curTime
    if (doSampling) {
      curTime = sampleTimes[ist];
    }
    else {
      curTime = pulseTimes[ipt];
    }

    // Sampling  
    if (doSampling) {
      sampleSim(curTime/fs_to_au);
    }

    ///// Time Evolution Through Pulse  /////
    while (doPulse) {

      if (((pulseTimes[ipt] + pulseTlengthAU) > sampleTimes[ist])) {
        endPulseSim = sampleTimes[ist];
        doSampling = true;
      }
      else {
        endPulseSim = pulseTimes[ipt] + pulseTlengthAU;
        doSampling = false;
        if (!firstPulse && (Npulses > 1) && (curTime == pulseTimes[ipt])) {
          for (iv=0; iv<NvibStates; iv++) {
            for (im=0; im<=MAXj; im++) {
              for (j=0; j<=MAXj-im; j++) {
                eigenAmps[iv][j][im] = pulseEvolMat[iv][im]*eigenAmps[iv][j][im];
              }
            }
          }
        curTime = endPulseSim;
        }
      }


      // Increment until next sampling or pulse end
      while (curTime < endPulseSim) {

        if (curTime + dTimeEvolStepAU < endPulseSim) {
          dTimePulseEvol = dTimeEvolStepAU;
          regTimeStep = true;
        }
        else {
          dTimePulseEvol = endPulseSim - curTime;
          regTimeStep = false;
        }

        // Time evolution of amplitudes
        double t0 = fmod(curTime, pulseTspacingAU);
        t0 > pulseTlengthAU ? t0 = 0 : t0 = t0;
        double E2_t = E0*E0*deltaAlphaAU
                                //*pow(sin(PI*fmod(curTime, pulseTspacingAU)/pulseTlengthAU),2)
                        *timeInt(t0 + dTimePulseEvol, t0)/4.0;
        for (im=0; im<=MAXj; im++) {
          for (j=0; j<=MAXj-im; j++) {
            complexNumber.real(0);
            complexNumber.imag(E2_t*cosSqEigVal[im][j]);
            cosSqExpDiag[im].coeffRef(j, j) = exp(complexNumber);
          }

          for (iv=0; iv<NvibStates; iv++) {
            if (regTimeStep) {
              timeEvolMat = pulsedTgenL[iv][im]*cosSqExpDiag[im]*pulsedTgenR[iv][im];
            }
            else {
              for (j=0; j<=MAXj-im; j++) {
                complexNumber.real(0);
                complexNumber.imag(-rotEnergyJ[iv][j+im]*dTimePulseEvol/2.0);
                H0dTgen[iv][im].coeffRef(j, j) = exp(complexNumber);
              }
              timeEvolMat = H0dTgen[iv][im]*cosSqTran[im]*cosSqExpDiag[im]*cosSqTran[im].adjoint()*H0dTgen[iv][im];
            }
            timeEvolMat.prune(prunefx);

            // evolution
            for (j=0; j<=MAXj-im; j++) {
              eigenAmps[iv][j][im] = timeEvolMat*eigenAmps[iv][j][im];
            }
            if (firstPulse) {
              pulseEvolMat[iv][im] = timeEvolMat*pulseEvolMat[iv][im];
            }
          }
        }

        // Update curTime
        if (regTimeStep) {
          curTime += dTimeEvolStepAU;
        }
        else {
          curTime = endPulseSim;
        }
      }

      if (doSampling) {
        sampleSim(curTime/fs_to_au);
      }
      else {
        curTime = pulseTimes[ipt] + pulseTlengthAU;
        doPulse = false;
        firstPulse = false;
        ipt++;
      }
    }
  }

/*
  vector<PLOToptions> opts(3);
  vector<string> optVals(3);
  opts[0] = xSpan;              optVals[0] = to_string(startSampleTimeAU/fs_to_au/1e3)+","+to_string(endSampleTime);
  opts[1] = xLabel;             optVals[1] = "Time [ps]";
  opts[2] = yLabel;             optVals[2] = "<Cos^{2}#theta>";

  //pltRR.print1d(cosSqEVals, "n2oAlignment", opts, optVals);
  //pltRR.print1d(cosSqEVals, "n2oAlignment");
*/

  /*
  jPopulation.resize(MAXj+1, 0);
  qrtrRv = 1;
  double prob, mult;
  for (int j=0; j<=MAXj+1; j++) {
    prob = 0;
    if (hasQuarterRev) {
      qrtrRv = ((j + 1)%2) + 1;
    }
    for (im=0; im<=MAXj+1; im++) {
      if (im > j) {
        break;
      }
      mult = im%1 + 1;   // Account for +/- m values
      for (iv=0; iv<NvibStates; iv++) {
        for (int j0=im; j0<MAXj+1; j0++) {
          //prob += mult*qrtrRv*rotThermalDist[iv][j0]
          //  *std::norm(eigenAmps[iv][j0][im].coeffRef(j-im));
        }
      }
    }
    jPopulation[j] = prob;
  }
*/

  clockEnd = clock();

  cout << "Time to simulate: " << double(clockEnd - clockBegin)/CLOCKS_PER_SEC << endl;

  return;
}





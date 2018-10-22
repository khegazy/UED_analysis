#include "diffractionClass.h"

using namespace std;


DIFFRACTIONclass::DIFFRACTIONclass(double Iebeam_in, double elEnergy, string scatAmpPath) {

  Initialize("", NULL, 0, Iebeam_in, 0, elEnergy, 0, scatAmpPath);
}


DIFFRACTIONclass::DIFFRACTIONclass(MOLENSEMBLEMCclass* molEnsemble, double sMax_in, double Iebeam_in, double screenDist_in, double elEnergy, int Nbins_in, string scatAmpPath) {

  Initialize("", molEnsemble, sMax_in, Iebeam_in, screenDist_in, elEnergy, Nbins_in, scatAmpPath);
}


DIFFRACTIONclass::DIFFRACTIONclass(string name, MOLENSEMBLEMCclass* molEnsemble, double sMax_in, double Iebeam_in, double screenDist_in, double elEnergy, int Nbins_in, string scatAmpPath) {

  Initialize(name, molEnsemble, sMax_in, Iebeam_in, screenDist_in, elEnergy, Nbins_in, scatAmpPath);
}


void  DIFFRACTIONclass::Initialize(string name, MOLENSEMBLEMCclass* molEnsemble, double sMax_in, double Iebeam_in, double screenDist_in, double elEnergy, int Nbins_in, string scatAmpPath) {

  if (molEnsemble) {
    // Check molEnsemble->atomTypes
    for (int iatm=0; iatm<molEnsemble->NmolAtoms; iatm++) {
      if (molEnsemble->atomTypes.end() == find(	molEnsemble->atomTypes.begin(),
          					molEnsemble->atomTypes.end(),
          					molEnsemble->molecules->atoms[iatm]->atomType)) {
        cerr << "ERROR: Not all types of atoms are included in MOLENSEMBLEMCclass::atomTypes\n\n";
        exit(0);
      }
    }

    fillLineOut = false;

    molecules = molEnsemble->molecules;
    Nmols = molEnsemble->Nmols;
    atomTypes = molEnsemble->atomTypes;
  }
  else {
    Nmols = 0;
    molecules = NULL;
    atomTypes.resize(NUMBER_ATOMS);
    for (int ia=0; ia<(int)NUMBER_ATOMS; ia++) {
      atomTypes[ia] = (ATOMS)ia;
    }
  }

  // Max s inv(Angs)            (|s| = 4*PI*sin(theta/2)/lambda)
  sMax = sMax_in;

  // De Broglie wavelength angs
  lambda = 2*PI*C_AU/sqrt(pow(elEnergy*eV_to_au + C_AU*C_AU,2) - pow(C_AU,4)); //au
  lambda /= angs_to_au;  // angs
  k0 = 2*PI/lambda;

  // Scattering amplitudes
  importScatteringAmps(scatAmpPath);

  // Detector geometry
  Iebeam = Iebeam_in;
  screenDist = screenDist_in;
  xScreenInit = -screenDist*tan(2*asin(sMax*lambda/(4*PI)));
  zScreenInit = xScreenInit;
  xScreenFin = -xScreenInit;
  zScreenFin = -zScreenInit;
  screenSize = xScreenFin - xScreenInit;
  Nbins = Nbins_in;


  //lambda = H_EV/sqrt(2*Mel*elEnergy*(1+elEnergy/(2*Mel*C_CMpS*C_CMpS)));



  diffPattern.resize(Nbins);
  diffAtmPattern.resize(Nbins);
  diffMolPattern.resize(Nbins);
  sPattern.resize(Nbins);
  for (int ir=0; ir<Nbins; ir++) {
    diffPattern[ir].resize(Nbins, 0);
    diffAtmPattern[ir].resize(Nbins, 0);
    diffMolPattern[ir].resize(Nbins, 0);
    sPattern[ir].resize(Nbins, 0);
  }


  //FT2D_samplePattern = new TH2F((name+"FT2D_samplePattern").c_str(),(name+"FT2D_samplePattern").c_str(),Nbins,xSampleInit,xSampleFin,Nbins,ySampleInit,ySampleFin);
  //FT2D_FTpattern = new TH2F((name+"FT2D_FTpattern").c_str(),(name+"FT2D_FTpattern").c_str(),Nbins,xSampleInit,xSampleFin,Nbins,ySampleInit,ySampleFin);
}


DIFFRACTIONclass::~DIFFRACTIONclass() {

  //delete FT2D_samplePattern;
  //delete FT2D_FTpattern;
}


void DIFFRACTIONclass::importScatteringAmps(std::string scatAmpPath) {

  string filename, line;
  std::vector<double> inpScatAmps;
  scatSInterp.resize(20000);

  filename = scatAmpPath+"/"+"interpolationAngs.dat";
  ifstream angFile(filename.c_str());
  int ind = 0;
  while ( getline(angFile,line) ) {
    scatSInterp[ind] = 4*PI*sin(atof(line.c_str())*PI/(180*2))/lambda;
    ind++;
  }

  for (uint iatm=0; iatm<atomTypes.size(); iatm++) {
    inpScatAmps.clear();

    //filename = scatAmpPath+"/"+atomNames[atomTypes[iatm]]+"_dcs.dat";
    filename = scatAmpPath+"/"+atomNames[atomTypes[iatm]]+"_interpolation.dat";
    if (access(filename.c_str(), F_OK) == -1) {
      cerr << "ERROR: Cannot find file " << filename << endl << endl;
      exit(0);
    }


    ind = 0;
    ifstream dcsFile(filename.c_str());
    scatAmps[atomTypes[iatm]].resize(20000);
    while ( getline(dcsFile,line) ) {
      scatAmps[atomTypes[iatm]][ind] = atof(line.c_str());
      ind++;
    }


    /*
    // Open file and remove header text
    ifstream dcsFile(filename.c_str());
    for (int il=0; il<31; il++) getline(dcsFile,line);

    // Import scattering amplitudes
    cout<<"begin loading"<<endl;
    scatSInterp.resize(20000);
    FILE* angFile = fopen("/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/interpolationAngs.dat", "rb");
    fread(&scatSInterp[0], sizeof(double), 20000, angFile);
    fclose(angFile);

    while ( getline(dcsFile,line) ) {
      scatAmps[atomTypes[iatm]].resize(20000);
      cout<<"loading: "<<atomNames[atomTypes[iatm]]<<endl;
      FILE* scatFile = fopen(("/reg/neh/home/khegazy/simulations/scatteringAmplitudes/3.7MeV/"
          + atomNames[atomTypes[iatm]] + "_interpolationAngs.dat").c_str(), "rb");
      fread(&scatAmps[atomTypes[iatm]][0], sizeof(double), 20000, scatFile);
      fclose(scatFile);

      if (iatm == 0) {
        scatSInterp.push_back(
            4*PI*sin(atof((line.substr(2,9)).c_str())*PI/(180*2))/lambda);
        cout<<"scatAng: "<<scatSInterp[ind]<<endl;
        ind++;
      }

      //inpScatAmps.push_back(
      //|s| = 4*PI*sin(theta/2)/lambda
      scatAmps[atomTypes[iatm]].push_back(
      		//sqrt(atof((line.substr(26,11)).c_str())));  //cm
         	sqrt(atof((line.substr(39,11)).c_str())));  //a0
    }
    */

    /*
    scatAmps[atomTypes[iatm]] = 
        new boost::math::barycentric_rational<double>(
              scatSInterp.data(), inpScatAmps.data(), inpScatAmps.size());
              */
  }
}


void DIFFRACTIONclass::reset() {

  for (uint ir=0; ir<diffPattern.size(); ir++) {
    for (uint ic=0; ic<diffPattern[ir].size(); ic++) {
      diffPattern[ir][ic] = 0;
      diffAtmPattern[ir][ic] = 0;
      diffMolPattern[ir][ic] = 0;
      sPattern[ir][ic] = 0;
    }
  }

  //diffPattern->Reset();
  //diffMolPattern->Reset();
  //diffAtmPattern->Reset();
} 


Eigen::Vector3d DIFFRACTIONclass::sCalc(double xpos, double ypos, double zpos) {
// s is in angs

  Eigen::Vector3d vInit(0,k0,0);
  Eigen::Vector3d vDir(xpos, ypos, zpos);
  Eigen::Vector3d vFin = k0*vDir.normalized();

  return vInit - vFin;
}


double DIFFRACTIONclass::interpScatAmp(ATOMS atmT, double sInp) {
// theta must be in degrees


  /////////////////////////////////////////////////////////////////
  /////  Find index of the closest angle to theta in scatAng  /////
  /////////////////////////////////////////////////////////////////

  int ind;
  while (fabs(sInp - scatSInterp[ind]) >= fabs(sInp - scatSInterp[ind+1])) 
    ind++;

  return scatAmps[atmT][ind];

  ///////////////////////////
  /////  Interpolating  /////
  ///////////////////////////
  
  //return scatAmps[atmT][strt] + (scatAmps[atmT][end] - scatAmps[atmT][strt])/(scatSInterp[end] - scatSInterp[strt])*(sInp - scatSInterp[strt]);

  /*
  double angIntNorm = 1 - (scatSInterp[ind+1] - theta)/(scatSInterp[ind+1] - scatSInterp[ind]);


  return ((-0.5*scatAmps[atmT][fabs(ind-1)]+1.5*scatAmps[atmT][ind]
		-1.5*scatAmps[atmT][ind+1]+0.5*scatAmps[atmT][ind+2])*pow(angIntNorm,3)
	+(scatAmps[atmT][fabs(ind-1)]-2.5*scatAmps[atmT][ind]
                +2*scatAmps[atmT][ind+1]-0.5*scatAmps[atmT][ind+2])*pow(angIntNorm,2)
	+(-0.5*scatAmps[atmT][fabs(ind-1)] + 0.5*scatAmps[atmT][ind+1])*angIntNorm
	+ scatAmps[atmT][ind]);
        */
}


void DIFFRACTIONclass::diffPatternCalc() {

  double scatLength=0;
  int rflX, rflZ;
  double dX, dZ;
  double sumAtmDiff, sumMolDiff;
  Eigen::Vector3d s;  

  map<ATOMS, double> scatAmpInterp;

  for (int ix=0; ix<Nbins; ix++) {
    dX = xScreenInit + (xScreenFin-xScreenInit)*(ix + 0.5)/((double)Nbins);
    for (int iz=0; iz<Nbins/2; iz++) {
//double rad = sqrt(pow(ix-(Nbins/2-0.5), 2) + pow(iz-(Nbins/2-0.5), 2));
//if (rad<((8./40.)*Nbins/2) || rad>(14./40.)*Nbins/2) continue;
      dZ = zScreenInit + (zScreenFin - zScreenInit)*(iz + 0.5)/((double)Nbins);
      scatLength = 100*sqrt(dX*dX + dZ*dZ + screenDist*screenDist);	// Convert to cm
      s = sCalc(dX, screenDist, dZ);

      for (uint ia=0; ia<atomTypes.size(); ia++) {
        scatAmpInterp[atomTypes[ia]] = interpScatAmp(atomTypes[ia], s.norm());
      }

      sumAtmDiff=sumMolDiff=0;
      for (int imol=0; imol<Nmols; imol++) {
        for (auto& atomI : molecules[imol].atoms) {
	  sumAtmDiff +=	(Iebeam/pow(scatLength,2))
			  *pow(scatAmpInterp[atomI->atomType],2);

	  for (auto& atomJ : molecules[imol].atoms) {
	    if (atomI->index == atomJ->index) continue;
            
	    sumMolDiff += (Iebeam/pow(scatLength,2))
			    *scatAmpInterp[atomI->atomType]*scatAmpInterp[atomJ->atomType]
			    *cos(-s.dot(atomI->location - atomJ->location));
	  }
        }
      }

      diffAtmPattern[iz][ix] = sumAtmDiff;
      diffMolPattern[iz][ix] = sumMolDiff;
      diffPattern[iz][ix] = sumAtmDiff+sumMolDiff;
      sPattern[iz][ix] = s.norm();

      rflX = Nbins - 1 - ix;
      rflZ = Nbins - 1 - iz;
      diffAtmPattern[rflZ][rflX] = sumAtmDiff;
      diffMolPattern[rflZ][rflX] = sumMolDiff;
      diffPattern[rflZ][rflX] = sumAtmDiff+sumMolDiff;
      sPattern[rflZ][rflX] = s.norm();
    }
  }

  return;
}


std::vector< std::vector<double> > DIFFRACTIONclass::lineOut_uniform() {

  fillLineOut = true;
  std::vector< std::vector<double> > result = diffPatternCalc_uniform();
  fillLineOut = false;

  return result;
}


std::vector< std::vector<double> > DIFFRACTIONclass::azimuthalAvg_uniform() {

  fillLineOut = true;
  std::vector< std::vector<double> > result = diffPatternCalc_uniform();
  fillLineOut = false;

  return result;
}


std::vector< std::vector<double> > DIFFRACTIONclass::diffPatternCalc_uniform() {

  if (Nmols == 0) {
    cerr << "ERROR: Molecule ensemble must have at least one entry to "
		<< "make reference diffraction pattern!!!\n\n";
    exit(0);
  } 

  Eigen::Vector3d s;  
  double dX, dZ;
  double scatLength=0;
  double sumAtmDiff, sumMolDiff;
  double cntrBins = (1 - Nbins%2)*0.5;    // Evaluate at center of bins

  map<ATOMS, double> scatAmpInterp;

  int NbinsX, NbinsZ;
  // Update Nbins for lineouts
  std::vector< std::vector<double> > lineOuts(4);
  if (fillLineOut) {
    NbinsX = 1;
    NbinsZ = (int)(Nbins/2. + Nbins%2);
    for (uint il=0; il<lineOuts.size(); il++) {
      lineOuts[il].resize(NbinsZ, 0);
    }
  }
  else {
    NbinsX = Nbins;
    NbinsZ = Nbins;
  }


  /////////////////////////////////////////
  /////  Making diffraction patterns  /////
  /////////////////////////////////////////

  int ix, iz, xItr, zItr;
  double arg;
  for (ix=0; ix<NbinsX; ix++) {
    xItr = ix + cntrBins - (int)(Nbins/2);
    dX = (xScreenFin - xScreenInit)*xItr/Nbins;

    // Adjust parameters for line outs
    if (fillLineOut) {
      dX = 0;
    }

    for (iz=0; iz<NbinsZ; iz++) {
      zItr = iz + cntrBins - (int)(Nbins/2);
      dZ = (zScreenFin - zScreenInit)*zItr/Nbins;

      // Scattering length and S
      scatLength = 100*sqrt(dX*dX+dZ*dZ+screenDist*screenDist);	
      s = sCalc(dX, screenDist, dZ);

      // Interpolating scattering amplitude 
      for (uint ia=0; ia<atomTypes.size(); ia++) {
        scatAmpInterp[atomTypes[ia]] = interpScatAmp(atomTypes[ia], s.norm());
      }

      /////  Building diffraction patterns  /////
      sumAtmDiff = 0;
      sumMolDiff = 0;
      for (auto& atomI : molecules[0].atoms) {
        // Atomic Diffraction
        sumAtmDiff += (Iebeam/pow(scatLength,2))
      			*pow(scatAmpInterp[atomI->atomType],2);

        // Molecular Diffraction
        for (auto& atomJ : molecules[0].atoms) {
          if (atomI->index != atomJ->index) {
 
            // Limit sinc(x) as x->0
	    if (s.norm()) {     
              arg         = s.norm()*(atomI->location - atomJ->location).norm();
              sumMolDiff += (Iebeam/pow(scatLength,2))
      	  		      *scatAmpInterp[atomI->atomType]*scatAmpInterp[atomJ->atomType]
      			      *sin(arg)/arg;
	    }
	    else {
              sumMolDiff += (Iebeam/pow(scatLength,2))
      	            	      *scatAmpInterp[atomI->atomType]*scatAmpInterp[atomJ->atomType];
	    }
          }
        }
      }
     

      if (fillLineOut) {
        lineOuts[1][NbinsZ - 1 - iz] = sumMolDiff;
        lineOuts[2][NbinsZ - 1 - iz] = sumAtmDiff;
        lineOuts[3][NbinsZ - 1 - iz] = s.norm();
        lineOuts[0][NbinsZ - 1 - iz] = sumMolDiff + sumAtmDiff;
      }
      else {
        diffAtmPattern[iz][ix] = sumAtmDiff;
        diffMolPattern[iz][ix] = sumMolDiff;
        diffPattern[iz][ix] = sumAtmDiff+sumMolDiff;
        sPattern[iz][ix] = s.norm();

      }
    }
  }

  if (fillLineOut) {
    return lineOuts;
  }

  return diffPattern;
}


void DIFFRACTIONclass::sPatternCalc() {

  double dX, dZ;
  double scatLength;
  Eigen::Vector3d s;  

  for (int ix=0; ix<Nbins; ix++) {
    dX = xScreenInit + (xScreenFin-xScreenInit)*ix/Nbins;
    for (int iz=0; iz<Nbins; iz++) {
      dZ = zScreenInit + (zScreenFin-zScreenInit)*iz/Nbins;
      scatLength = 100*sqrt(dX*dX+dZ*dZ+screenDist*screenDist);	// Convert to cm
      s = sCalc(dX, screenDist, dZ);

      sPattern[iz][ix] = s.norm();
    }
  }

  return;
}

/*
double DIFFRACTIONclass::gaussXsect(double x, double y, double xs, double stdev) {

  return (x/(stdev*sqrt(3.14159)))*exp(-(x^2+y^2)/(stdev^2));
}
*/

/*
void DIFFRACTIONclass::FT2d() {

  TVector3 s;
  double dX,dZ;
  int iang;
  double binSizeX = (xSampleFin-xSampleInit)/Nbins;
  double binSizeY = (ySampleFin-ySampleInit)/Nbins;
  double* in;
  fftw_complex* out;
  in = (double*) fftw_malloc(sizeof(double)*Nbins*Nbins);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nbins*Nbins);
  fftw_plan p = fftw_plan_dft_r2c_2d(Nbins, Nbins, in, out, FFTW_ESTIMATE);
for (int i=0; i<Nbins*Nbins; i++) in[i]=0; //if(in[i]!=0) cout<<"Entry "<<i<<":  "<< in[i]<<endl;

cout<<"111"<<endl;
  for (int imol=0; imol<Nmols; imol++) {
    for (int iatm=0; iatm<molecules[imol].Natoms; iatm++) {
      //cout<<"index y,x:  "<<int((molecules[imol].atoms[iatm]->location.Y()-zInitFt)/binSizeY)<<"   "<<int((molecules[imol].atoms[iatm]->location.X()-xInitFt)/binSizeX)<<endl<<molecules[imol].atoms[iatm]->location.Y()<<"   "<<molecules[imol].atoms[iatm]->location.X()<<endl<<binSizeX<<"   "<<binSizeY<<endl;
      in[Nbins*int((molecules[imol].atoms[iatm]->location.Y()-ySampleInit)/binSizeY) + int((molecules[imol].atoms[iatm]->location.X()-xSampleInit)/binSizeX)]++;
    }
  }

for (int i=0; i<Nbins*Nbins; i++) { if(in[i]!=0) cout<<"Nonzero: "<<in[i]<<"   "<<i<<endl;}

cout<<"222"<<endl;

  
cout<<"333"<<endl;
  fftw_execute(p);

cout<<"444"<<endl;
  iang=0;

  for (int ix=0; ix<Nbins; ix++) {
    dX = xScreenInit+(xScreenFin-xScreenInit)*ix/Nbins;
    for (int iz=0; iz<Nbins; iz++) {
      dZ = zScreenInit+(zScreenFin-zScreenInit)*iz/Nbins;
      sumAtmDiff=sumMolDiff=0;
      s = sCalc(dX, dZ, screenDist);
      iang = angleIndex(atan(sqrt(dX*dX+dZ*dZ)/screenDist)*180/PI,iang);
      for (int imol=0; imol<Nmols; imol++) {
        for (int iatm=0; iatm<molecules[imol].Natoms; iatm++) {
          sumAtmDiff += scatDCS[molecules[imol].atoms[iatm]->atomType][iang];
        }
      }
      sPattern->SetBinContent(ix+1,iz+1, s.norm());
      diffAtmPattern->SetBinContent(ix+1, iz+1, sumAtmDiff);
      FT2D_samplePattern->SetBinContent(ix+1, iz+1, in[Nbins*iz+ix]);
//cout<<"Bin:  "<<ix<<"  "<<iz<<"      diffPattern: "<<out[iz*Nbins+ix][0]<<"  "<<out[iz*Nbins+ix][1]<<"      input: "<<in[Nbins*iz+ix]<<endl; 
    //diffPattern->SetBinContent(ix+1, iz+1, sumMolDiff+sumAtmDiff);
    }
  }

  for (int ix=0; ix<Nbins/2+1; ix++) for (int iz=0; iz<Nbins/2; iz++) {
    FT2D_FTpattern->SetBinContent(ix+1, Nbins/2+iz+1, (pow(out[iz*(Nbins/2+1)+ix][0],2) + pow(out[iz*(Nbins/2+1)+ix][1],2)));
    FT2D_FTpattern->SetBinContent(Nbins-ix, Nbins/2+iz+1, FT2D_FTpattern->GetBinContent(ix+1,Nbins/2+iz+1));
    FT2D_FTpattern->SetBinContent(ix+1, iz+1,(pow(out[(Nbins/2)*(Nbins/2+1)-iz*(Nbins/2+1)+ix][0],2) + pow(out[(Nbins/2)*(Nbins/2+1)-iz*(Nbins/2+1)+ix][1],2))); 
    FT2D_FTpattern->SetBinContent(Nbins-ix, iz+1, FT2D_FTpattern->GetBinContent(ix+1,iz+1));
    diffPattern->SetBinContent(ix+1, Nbins/2+iz+1, diffAtmPattern->GetBinContent(ix+1,Nbins/2+iz+1)*(pow(out[iz*(Nbins/2+1)+ix][0],2) + pow(out[iz*(Nbins/2+1)+ix][1],2)));
    diffPattern->SetBinContent(Nbins-ix, Nbins/2+iz+1, diffPattern->GetBinContent(ix+1,Nbins/2+iz+1));
    diffPattern->SetBinContent(ix+1, iz+1, diffAtmPattern->GetBinContent(ix+1,iz+1)*(pow(out[(Nbins/2)*(Nbins/2+1)-iz*(Nbins/2+1)+ix][0],2) + pow(out[(Nbins/2)*(Nbins/2+1)-iz*(Nbins/2+1)+ix][1],2)));
    diffPattern->SetBinContent(Nbins-ix, iz+1, diffPattern->GetBinContent(ix+1,iz+1));
    diffMolPattern->SetBinContent(ix+1, Nbins/2+iz+1, diffPattern->GetBinContent(ix+1,Nbins/2+iz+1));
    diffMolPattern->SetBinContent(Nbins-ix, Nbins/2+iz+1, diffPattern->GetBinContent(Nbins-ix,Nbins/2+iz+1));
    diffMolPattern->SetBinContent(ix+1, iz+1, diffPattern->GetBinContent(ix+1,iz+1));
    diffMolPattern->SetBinContent(Nbins-ix, iz+1, diffPattern->GetBinContent(Nbins-ix,iz+1));
 }

 diffMolPattern->Add(diffAtmPattern,-1);
//for (int iz=0; iz<Nbins/2; iz++) cout<<"Bin "<<iz+1<<"  "<<Nbins-iz<<" :   "<<diffPattern->GetBinContent(25,iz+1)<<"      "<<diffPattern->GetBinContent(25,Nbins-iz)<<endl;

cout<<"555"<<endl;
  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);

  return;
}
*/ 
  




/*
void DIFFRACTIONclass::FT2d(MOLECULEclass* mol) {

  double Iatm;
  iang=0;

  for (int row=0; row<=Nrows; row++) {
    dX = xInit+(xFin-xInit)*row/Nrows;
    for (int col=0; col<=Ncols; col++) {
      dZ = zInit+(zFin-zInit)*col/Ncols;
      sumAtmDiff=sumMolDiff=0;
      s = sCalc(dX, dZ, screenDist);
      iang = angleIndex(atan(sqrt(dX*dX+dZ*dZ)/screenDist)*180/PI,iang);
      for (int imol=0; imol<Nmol; imol++) {
        for (int iatm=0; iatm<mol->Natoms; iatm++) {
          Iatm = scatDCS[mol->atoms[iatm]->atomType][iang];
          cout<<"Angle: "<<scatSInterp[iang]<<"        scatDCS: "<<scatDCS[mol->atoms[iatm]->atomType][iang]<<endl;
          sumAtmDiff += Iatm;
          sumMolAmp += sqrt(Iatm)*exp(-complex<double>(0,s.X()*(mol->atoms[iatm]->location.X())+s.Y()*(mol->atoms[iatm]->location.Y())));
        }
      sumMolDiff += norm(sumMolAmp) - sumAtmDiff;
      }
    sPattern->SetBinContent(row+1,col+1, s.norm());
    diffAtmPattern->SetBinContent(row+1, col+1, sumAtmDiff);
    diffMolPattern->SetBinContent(row+1, col+1, sumMolDiff);
    diffPattern->SetBinContent(row+1, col+1, sumMolDiff+sumAtmDiff);
    }
  }

  return;
} 
*/  


/*
void diffMolCalc(MOLECULEclass* mol, double* s, double &diffAtmBin, double &diffMolBin, double diffBin) {

  int iatm=0;
  diffAtmBin += IatmCalc(mol->baseAtom);
  diffBin += diffAtmBin + diffMolBin;

  for (int ib=0; ib<mol->baseAtom->Nbonds, ib++) diffMolCalcRec(mol->baseAtom, mol->baseAtom->atoms[ib], s, diffAtmBin, diffMolBin, diffBin);
  return;
}


void diffMolCalcRec(ATOMclass* atomI, ATOMclass atomF, double*, double &diffAtmBin, double &diffMolBin, double diffBin) {

  diffAtmBin += IatmCalc(atomF);
  diffMolBin += ImolCalc(atomI,);
  diffBin += diffAtmBin + diffMolBin;

  for ...  diffMolCalcRec(iatm++, mol, molStr, s, diffAtmBin, diffMolBin, diffBin); 
}
*/

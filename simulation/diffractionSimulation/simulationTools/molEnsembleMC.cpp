#include "molEnsembleMC.h"

using namespace std;


MOLENSEMBLEMCclass::MOLENSEMBLEMCclass(long int seed, std::string xyzPath, std::string pdfPath, std::string pdfNames) 
      : PDFclass(seed, pdfPath, pdfNames) {

  initialize();
  importXYZfile(xyzPath);
}


MOLENSEMBLEMCclass::MOLENSEMBLEMCclass(long int seed, std::string pdfPath, std::string pdfNames) 
      : PDFclass(seed, pdfPath, pdfNames) {

  initialize();
}


MOLENSEMBLEMCclass::MOLENSEMBLEMCclass(long int seed, std::string xyzPath) 
      : PDFclass(seed, "NULL", "") {

  initialize();
  importXYZfile(xyzPath);
}


MOLENSEMBLEMCclass::MOLENSEMBLEMCclass(long int seed) 
      : PDFclass(seed, "NULL", "") {

  initialize();
}


void MOLENSEMBLEMCclass::initialize() {

  // Debugging
  verbose = 1;

  molecules   = NULL;
  NmolSamples = 1;
  NmolAtoms   = -1;
  Nmols       = -1;
  Ngr = 0;
  
  usePositionMC = false;
  useOrientationMC = false;

  charges[H] = 1;
  charges[C] = 6;
  charges[N] = 7;
  charges[O] = 8;
  charges[I] = 53;

  testAzm = new TH1F("testAzm", "testAzm", 50, -5, 5);
  testPol = new TH1F("testPol", "testPol", 50, 0 , 12);
  orientDist = new TH2F("orientDist", "orientDist", 50, 0 , PI, 50, 0, 2*PI);
  orientGraph = new TGraph2D();
  orientGraph->SetTitle("orientGraph");
}


void MOLENSEMBLEMCclass::importXYZfile(std::string fileName) {

  NmolAtoms = -1;
  ifstream XYZ(fileName.c_str());
  std::string str;

  // Get atom positions
  int pInd, fInd;
  std::string xPos, yPos, zPos, name;
  std::map< std::string, Eigen::Vector3d > curAtomPos;
  while (std::getline(XYZ, str)) {
  
    atomCount.clear();
    curAtomPos.clear();

    // Number of atoms
    int curNmolAtoms = stoi(str);
    if (NmolAtoms == -1) {
      NmolAtoms = curNmolAtoms;
    }
    else if (NmolAtoms != curNmolAtoms) {
      std::cerr << "ERROR: Number of atoms in the 'same' molecule are not equal ["
        + to_string(NmolAtoms) + ", " + to_string(curNmolAtoms) + "!!!\n";
      exit(0);
    }

    // Skipping extra line
    std::getline(XYZ, str);

    for (int i=0; i<NmolAtoms; i++) {
      if (std::getline(XYZ, str)) {
        pInd = tools::stringRemove(str, ' ', 0);
        fInd = str.find(' ', pInd);
        std::string atm = str.substr(pInd, fInd - pInd);
        if (std::find(atomTypes.begin(), atomTypes.end(), getATOMtype(atm.c_str())) 
            == atomTypes.end()) {
          atomTypes.push_back(getATOMtype(atm));
        }

        pInd = tools::stringRemove(str, ' ', 3);
        fInd = str.find(' ', pInd);
        xPos = str.substr(pInd, fInd - pInd);
        pInd = fInd;

        pInd = tools::stringRemove(str, ' ', str.find(' ', pInd));
        fInd = str.find(' ', pInd);
        yPos = str.substr(pInd, str.find(' ', pInd) - pInd);
        pInd = fInd;

        pInd = tools::stringRemove(str, ' ', str.find(' ', pInd));
        fInd = str.find(' ', pInd);
        zPos = str.substr(pInd, str.find(' ', pInd) - pInd);
        pInd = fInd;

        name = atm + "_" + to_string(atomCount[atm]);
        atomCount[atm]++;

        Eigen::Vector3d pos(stof(xPos), stof(yPos), stof(zPos));
        curAtomPos[name] = pos;
      }
      else {
        cerr << "ERROR: Ran out of lines in the XYZ file before filling all atom positions!!!\n\n";
        exit(0);
      }
    }

    atomPos.push_back(curAtomPos);
  }

}


MOLENSEMBLEMCclass::~MOLENSEMBLEMCclass() {

  if (molecules) {
    delete[] molecules;
  }

  TFile* f = TFile::Open("testAng.root", "RECREATE");
  testAzm->Write();
  testPol->Write();
  f->Close();

  TCanvas* MyC = new TCanvas("MyCC", "MyCC", 800, 600);
  orientDist->Draw();
  MyC->Print("orientDist_test.png");
  delete MyC;

  delete testAzm;
  delete testPol;
  delete orientDist;
  delete orientGraph;
}

void MOLENSEMBLEMCclass::reset() {

  delete[] molecules;
  molecules = NULL;
}


ATOMS MOLENSEMBLEMCclass::getATOMtype(std::string atm) {

  if (atm.compare("H") == 0)      return H;
  if (atm.compare("C") == 0)      return C;
  if (atm.compare("N") == 0)      return N;
  if (atm.compare("O") == 0)      return O;
  else {
    cerr << "ERROR: Cannot find atom type (" << atm << ")!!!\n";
    exit(0);
  }
}


void MOLENSEMBLEMCclass::makeMolEnsemble() {

  std::map<std::string, double> empty;
  makeMolEnsemble(empty);
}

void MOLENSEMBLEMCclass::makeMolEnsemble(std::map<string, double> inpVals) {

  /*
  if (Nmols == -1) {
    cerr << "ERROR: Did not set the number of molecules desired (MOLECULEclass::Nmols)!!!\n\n";
    exit(0);
  }
  */
  if (NmolAtoms == -1) {
    cerr << "ERROR: Did not set the number of atoms per molecule (MOLECULEclass::NmolAtoms)!!!\n\n";
    exit(0);
  }
  if (!atomTypes.size()) {
    cerr << "ERROR: Did not set the types of atoms (MOLECULEclass::atomTypes)!!!\n\n";
    exit(0);
  }
  if (molecules) {
    cerr << "ERROR: Array of molecules is already set, you must reset first!!!\n\n"; 
    exit(0);
  }

  if (useOrientationMC) {
    if ((orientationPDF == "Uniform") || (orientationPDF == "uniform") 
		|| (orientationPDF == "UNIFORM")) {
      orientationPDF = "uniform";
    }
    if (verbose) {
      cout << "INFO:    Orientation MC will be applied with the following distributions:\n";
      cout << "\t\tOrientation distribution: " << orientationPDF << endl;
    }
  }
  else if (verbose) {
    cout << "INFO:    Will not apply orientation MC!\n";
  }
  if (verbose) {
    if (usePositionMC) {
      cout << "INFO:    Position MC will be applied with the following distributions:\n";
      cout << "\t\tStandard deviation of distribution (if needed): " << stdDevPositionLength << endl;
      cout << "\t\tX position PDF: " << XpositionPDF << endl;
      cout << "\t\tY position PDF: " << YpositionPDF << endl;
      cout << "\t\tZ position PDF: " << ZpositionPDF << endl;
    }
    else {
      cout<< "INFO:   Will not apply position MC!\n";
    }
  }

  NdiffMols = atomPos.size();
  cout<<"Number of diff mols: "<<NdiffMols<<endl;
  if (Nmols == -1) {
    Nmols = NmolSamples*NdiffMols;
  }
  else if (Nmols != NmolSamples*NdiffMols) {
    int NcurMols = (int)atomPos.size();
    atomPos.resize(Nmols);
    int ind=0;
    for (int i=0; i<Nmols; i++) {
      ind = i%NcurMols;
      atomPos[i] = atomPos[ind];
    }
  }
  NdiffMols = atomPos.size();


  molecules = new MOLECULEclass[Nmols];
 

  int cmol = 0;
  for (int imol=0; imol<NdiffMols; imol++) {
    for (int ism =0; ism<NmolSamples; ism++) {

      // Build an individual molecule
      buildMolecule(molecules[cmol], imol, inpVals);

      // Orientation MC
      if (useOrientationMC) orientMolecule(molecules[cmol]);

      // Set the position of the molecule
      if (usePositionMC) positionMolecule(molecules[cmol]);

      cmol++;
    }
  }
}


void MOLENSEMBLEMCclass::orientMolecule(MOLECULEclass &molecule) {

  pair<double, double> angs = sampleOrientPDF(orientationPDF);

  // Switch pol/azm angles from simulation since beam axis is 
  //     not z unlike laser polarization in simulation
  double randPolAng = angs.first;
  double randAzmAng = angs.second; 

//randAzmAng = 0;
  // Correcting for simulation switching theta and phi
  //	Z axis in simulation is Y in here
/*
  double randPolAng = samplePDF("uniform");
  double randAzmAng = samplePDF("uniform");
  if (orientationPDF == "uniform") {
cout<<"OLDrands: "<<randAzmAng<<"   "<<randPolAng<<endl;
    randAzmAng = acos(2*randAzmAng - 1);	// norm integral of sin 0->PI
    randPolAng *= 2*PI;
cout<<"OLDangs: "<<randAzmAng<<"   "<<randPolAng<<endl<<endl;
  }
*/

  testAzm->Fill(randAzmAng);
  testPol->Fill(randPolAng);
  orientDist->Fill(randPolAng, randAzmAng);
  // orig orientDist->Fill(randAzmAng, randPolAng);
  orientGraph->SetPoint(Ngr,
               sin(randPolAng)*cos(randAzmAng),
               sin(randPolAng)*sin(randAzmAng),
               cos(randPolAng));
  Ngr++;

  Eigen::Matrix3d azimRotate;
  azimRotate  = Eigen::AngleAxisd(randAzmAng, Eigen::Vector3d::UnitY());
  Eigen::Matrix3d polRotate;
  polRotate   = Eigen::AngleAxis<double>(randPolAng, Eigen::Vector3d::UnitX());
  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    if ((molecule.atoms[iatm]->location(0) == 0)
        && (molecule.atoms[iatm]->location(1) == 0) 
        && (molecule.atoms[iatm]->location(2) == 0)) {
      continue;
    }

    // Axis of rotation: Rotate XZ plane projection by 90 degrees and normalize
    //Eigen::Vector3d polAxis(
    //    -1*molecule.atoms[iatm]->location(2),
    //    0,
    //    molecule.atoms[iatm]->location(0));
    //polAxis /= polAxis.norm();
    //polRotate = Eigen::AngleAxis<double>(randPolAng, polAxis);
    molecule.atoms[iatm]->location = polRotate*molecule.atoms[iatm]->location;
    
    molecule.atoms[iatm]->location = azimRotate*molecule.atoms[iatm]->location;


    /*
    rad = molecule.atoms[iatm]->location.Mag();
    polAng = molecule.atoms[iatm]->location.Theta();
    azmAng = molecule.atoms[iatm]->location.Phi();
    molecule.atoms[iatm]->setLocation( 	
		//molecule.atoms[iatm]->location.X() + rad*sin(polAng)*cos(azmAng),
		//molecule.atoms[iatm]->location.Y() + rad*sin(polAng)*sin(azmAng),
		//molecule.atoms[iatm]->location.Z() + rad*cos(polAng) 	
		rad*sin(polAng+randPolAng)*cos(azmAng+randAzmAng),
		rad*sin(polAng+randPolAng)*sin(azmAng+randAzmAng),
		rad*cos(polAng+randPolAng) 	
    );
    */
  }
}


void MOLENSEMBLEMCclass::positionMolecule(MOLECULEclass &molecule) {

  double Xpos = stdDevPositionLength*samplePDF(XpositionPDF);
  double Ypos = stdDevPositionLength*samplePDF(YpositionPDF);
  double Zpos = stdDevPositionLength*samplePDF(ZpositionPDF);

  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    molecule.atoms[iatm]->setLocation(
        molecule.atoms[iatm]->location(0) + Xpos,
        molecule.atoms[iatm]->location(1) + Ypos,
        molecule.atoms[iatm]->location(2) + Zpos);
  }
}


std::vector<double> MOLENSEMBLEMCclass::simulatePairCorr(int Nbins, float maxR, bool smear, PLOTclass* pltVerbose) {

  MOLECULEclass &molecule = molecules[0]; 
  int ind, chargeI;
  float dist;
  std::vector<double> counts(Nbins, 0);
  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    chargeI = charges[molecule.atoms[iatm]->atomType];
    for (int jatm=iatm+1; jatm<NmolAtoms; jatm++) {
      dist = (molecule.atoms[iatm]->location 
          - molecule.atoms[jatm]->location).norm();
      if (dist > maxR) {
        continue;
      }
      ind = (int)(Nbins*(dist/maxR));
      counts[ind] += chargeI*charges[molecule.atoms[jatm]->atomType];
    }
  }
  
  if (pltVerbose)
    pltVerbose->print1d(counts, "./plots/simPairCorrCounts");

  if (smear) {
    float var = std::pow(Nbins*(0.1/maxR), 2);
    std::vector<double> pairCorr(Nbins, 0);
    for (int i=0; i<Nbins; i++) {
      for (int ii=0; ii<Nbins; ii++) {
        pairCorr[i] += exp(-1*std::pow(i-ii, 2)/(2*var))*counts[ii];
      }
    }

    if (pltVerbose)
      pltVerbose->print1d(pairCorr, "./plots/simPairCorr");
   
    return pairCorr;
  }

  return counts;
}
 

std::map< std::string, std::vector<double> >* MOLENSEMBLEMCclass::getBonds(
      std::map< std::string, std::vector<double> >* bonds) {

  MOLECULEclass &molecule = molecules[0]; 
  float dist;
  std::string bondType;
  if (!bonds) {
    bonds = new std::map< std::string, std::vector<double> >;
  }
  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    for (int jatm=iatm+1; jatm<NmolAtoms; jatm++) {
      if (molecule.atoms[iatm]->atomType < molecule.atoms[jatm]->atomType) {
        bondType  = atomNames[molecule.atoms[iatm]->atomType] + "-"
                      + atomNames[molecule.atoms[jatm]->atomType];
      }
      else {
        bondType  = atomNames[molecule.atoms[jatm]->atomType] + "-"
                      + atomNames[molecule.atoms[iatm]->atomType];
      }
      dist = (molecule.atoms[iatm]->location 
          - molecule.atoms[jatm]->location).norm();
      (*bonds)[bondType].push_back(dist);
    }
  }
 
  return bonds;
}
   
TGraph2D* MOLENSEMBLEMCclass::testBuildMolecule() {
  std::map<std::string, double> empty;
  return testBuildMolecule(empty);
}


TGraph2D* MOLENSEMBLEMCclass::testBuildMolecule(std::map<std::string, double> inpVals) {

  TFile* f = TFile::Open("testplot.root", "RECREATE");

  TGraph2D* testMolecule = new TGraph2D();
  //TGraph2D* testMolecule = new TGraph2D("testbuild", "testbuild", 51, -10, 10, 51, -10, 10, 51, -10, 10);
  MOLECULEclass molecule;
  buildMolecule(molecule, 0, inpVals);
cout<<"sizes: "<<NmolAtoms<<"   "<<molecule.atoms.size()<<endl;
  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    testMolecule->SetPoint(
        iatm,	
        molecule.atoms[iatm]->location(0), 
	molecule.atoms[iatm]->location(1),
	molecule.atoms[iatm]->location(2));
  }

  testMolecule->Write();
  f->Close();
  return testMolecule;
}


void MOLENSEMBLEMCclass::buildMolecule(
    MOLECULEclass &molecule,
    int imol,
    std::map<std::string, double> inpVals) {

    ATOMS aType;

    for (auto mItr : atomPos[imol]) {
      aType = getATOMtype(mItr.first.substr(0, mItr.first.find("_")));
      molecule.addAtom(new ATOMclass(mItr.first, aType, mItr.second));
    }
}




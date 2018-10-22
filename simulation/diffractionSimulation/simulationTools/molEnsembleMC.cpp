#include "molEnsembleMC.h"

using namespace std;



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

  molecules = NULL;
  Nmols = -1;
  NmolAtoms = -1;
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

  ifstream XYZ(fileName.c_str());
  std::string str;

  // Number of atoms
  std::getline(XYZ, str);
  NmolAtoms = stoi(str);
  cout<<NmolAtoms<<endl;

  // Get atom positions
  std::getline(XYZ, str);
  int pInd;
  std::string xPos, yPos, zPos, name;
  while (std::getline(XYZ, str)) {
    pInd = 1;

    std::string atm = str.substr(1, str.find(' ', 1) - 1);
    if (std::find(atomTypes.begin(), atomTypes.end(), getATOMtype(atm.c_str())) 
        == atomTypes.end()) {
      atomTypes.push_back(getATOMtype(atm));
    }

    pInd = tools::stringRemove(str, ' ', 3);
    xPos = str.substr(pInd, str.find(' ', pInd) - pInd);

    pInd = tools::stringRemove(str, ' ', str.find(' ', pInd));
    yPos = str.substr(pInd, str.find(' ', pInd) - pInd);

    pInd = tools::stringRemove(str, ' ', str.find(' ', pInd));
    zPos = str.substr(pInd, str.find(' ', pInd) - pInd);

    name = atm + "_" + to_string(atomCount[atm]);
    atomCount[atm]++;

    Eigen::Vector3d pos(stof(xPos), stof(yPos), stof(zPos));
    atomPos[name] = pos;
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

  if (Nmols == -1) {
    cerr << "ERROR: Did not set the number of molecules desired (MOLECULEclass::Nmols)!!!\n\n";
    exit(0);
  }
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


  molecules = new MOLECULEclass[Nmols];
  
  for (int imol=0; imol<Nmols; imol++) {

    // Build an individual molecule
    buildMolecule(molecules[imol], inpVals);

    // Orientation MC
    if (useOrientationMC) orientMolecule(molecules[imol]);

    // Set the position of the molecule
    if (usePositionMC) positionMolecule(molecules[imol]);
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
  orientDist->Fill(randAzmAng, randPolAng);
  orientGraph->SetPoint(Ngr,
               sin(randPolAng)*cos(randAzmAng),
               sin(randPolAng)*sin(randAzmAng),
               cos(randPolAng));
  Ngr++;

  double rad, azmAng, polAng;

  for (int iatm=0; iatm<NmolAtoms; iatm++) {
    cout<<"FIX THIS: switched ot eigen, look into anglaxis or calculate theta and phi"<<endl;
    exit(0);
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
   
TGraph2D* MOLENSEMBLEMCclass::testBuildMolecule() {
  std::map<std::string, double> empty;
  return testBuildMolecule(empty);
}


TGraph2D* MOLENSEMBLEMCclass::testBuildMolecule(std::map<std::string, double> inpVals) {

  TFile* f = TFile::Open("testplot.root", "RECREATE");

  TGraph2D* testMolecule = new TGraph2D();
  //TGraph2D* testMolecule = new TGraph2D("testbuild", "testbuild", 51, -10, 10, 51, -10, 10, 51, -10, 10);
  MOLECULEclass molecule;
  buildMolecule(molecule, inpVals);
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


void MOLENSEMBLEMCclass::buildMolecule(MOLECULEclass &molecule,
    std::map<std::string, double> inpVals) {

    ATOMS aType;

    for (auto mItr : atomPos) {
      aType = getATOMtype(mItr.first.substr(0, mItr.first.find("_")));
      molecule.addAtom(new ATOMclass(mItr.first, aType, mItr.second));
    }
}




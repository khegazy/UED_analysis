#include "pdfClass.h"
#include <sstream>

using namespace std;



PDFclass::PDFclass(long int seed, string pdfPath, string pdfList) {

  rand = TRandom3(seed);

  appliedJacobian = false;

  pdfFile = NULL;
  if (pdfPath.length() && (pdfPath != "NULL")) {
    cout<<"Trying to get PDF"<<endl;
    pdfFile = TFile::Open(pdfPath.c_str());
    cout<<"Got PDF"<<endl;
    if (!pdfFile->IsOpen()) {
      cerr<<"ERROR: Cannot open PDF file "<<pdfPath<<"!!!"<<endl;
      exit(0);
    }
    else {
      cout << "INFO: Using pdf file " << pdfPath << "!\n";
    }

    istringstream ss(pdfList);
    string name;

    while (getline(ss, name, ',')) {
      cout << "\tINFO: Getting pdf: " << name << endl;
      pdfs[name] = (TH1*)pdfFile->Get(name.c_str());
      pdfNames.push_back(name);

      if (!pdfs[name]) {
	cerr << "ERROR: Cannot find pdf " + name 
	    + " in file " + pdfPath + "!!!\n\n";
      }
    }   
  }
    
}


PDFclass::~PDFclass() {

  for (auto& ih : pdfs) {
    delete ih.second;
  }

  if (pdfFile) {
    pdfFile->Close();
  }
}


double PDFclass::samplePDF(const string pdfName) {

  if (std::find(pdfNames.cbegin(), pdfNames.cend(), pdfName) != pdfNames.end()) {
    /////  Sampling PDF distribution  /////
    int ibx;
    TH1D* pdf = (TH1D*)pdfs[pdfName];
 
    double sample = rand.Uniform();
 
    for (ibx=1; ibx<=pdf->GetNbinsX() && sample>0; ibx++) {
      sample -= pdf->GetBinContent(ibx);
    }
  
    // FIX (MAKE BETTER): Currently choosing bin and picking place randomly within bin. Consider interpolation.
     return pdf->GetXaxis()->GetXmin() 
        + (pdf->GetXaxis()->GetXmax() - pdf->GetXaxis()->GetXmin())
        /pdf->GetNbinsX()*(ibx+rand.Uniform());
  } 
  else if ((pdfName=="Uniform") || (pdfName=="uniform") || (pdfName=="UNIFORM")) {
    return rand.Uniform();
  }
  else if ((pdfName=="Gaus") || (pdfName=="gaus") || (pdfName=="GAUS")) {
    return rand.Gaus();
  }
  else {
    cerr << "ERROR: PDF " + pdfName 
	+ " is not in list of pdfnames and thus not accounted for!!!\n\n";
    exit(0);
  }
}


pair<double, double> PDFclass::sampleOrientPDF(const string pdfName) {

  pair<double, double> res;

  if (std::find(pdfNames.cbegin(), pdfNames.cend(), pdfName) != pdfNames.end()) {
    /////  Sampling PDF distribution  /////
    int ibx, iby;
    TH2D* pdf = (TH2D*)pdfs[pdfName];

double norm; 
    // Apply sin(theta) jacobian to distribution
    if (!appliedJacobian) {
      double jacobian;
      for (ibx=0; ibx<pdf->GetNbinsX(); ibx++) {
        jacobian = -(cos((ibx+1)*PI/((double)pdf->GetNbinsX())) 
		   - cos(ibx*PI/((double)pdf->GetNbinsX())));
        for (iby=0; iby<pdf->GetNbinsY(); iby++) {
	  pdf->SetBinContent(ibx + 1, iby + 1,
			pdf->GetBinContent(ibx + 1, iby + 1)*jacobian);
	}
      }
      norm = pdf->Integral();
      pdf->Scale(1.0/norm);

      appliedJacobian = true;
    }

    // Sample distribution
    double sample = rand.Uniform();
    double theta;
 
    for (ibx=0; ibx<pdf->GetNbinsX() && sample>0; ibx++) {
      theta = (ibx - 1)*PI/((double)pdf->GetNbinsX()) 
		+ PI/((double)2*(pdf->GetNbinsX()));
      for (iby=0; iby<pdf->GetNbinsY() && sample>0; iby++) {
        sample -= pdf->GetBinContent(ibx + 1, iby + 1);
      }
    }
  
    // FIX (MAKE BETTER): Currently choosing bin and picking place randomly within bin. Consider interpolation.
     res.first = pdf->GetXaxis()->GetXmin() 
        	+ (pdf->GetXaxis()->GetXmax()-pdf->GetXaxis()->GetXmin())
        	/pdf->GetNbinsX()*(ibx+rand.Uniform());
     res.second = pdf->GetYaxis()->GetXmin() 
        	+ (pdf->GetYaxis()->GetXmax()-pdf->GetYaxis()->GetXmin())
        	/pdf->GetNbinsY()*(iby+rand.Uniform());
     return res;
  } 
  else if ((pdfName=="Uniform") || (pdfName=="uniform") || (pdfName=="UNIFORM")) {
    res.first = acos(1 - 2*rand.Uniform());
    res.second = 2*PI*rand.Uniform();
    return res;
  }
//  else if ((pdfName=="Gaus") || (pdfName=="gaus") || (pdfName=="GAUS")) {
//    return rand.Gaus();
//  }
  else {
    cerr << "ERROR: PDF " + pdfName 
	+ " is not in list of pdfnames and thus not accounted for!!!\n\n";
    exit(0);
  }
}


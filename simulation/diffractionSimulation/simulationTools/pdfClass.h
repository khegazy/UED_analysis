#ifndef PDFCLASS_H
#define PDFCLASS_H

#include <string>
#include <cmath>
#include <map>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <algorithm>

#include "../../../tools/constants.h"

using namespace std;


class PDFclass {

   public:
	PDFclass(long int seed, string pdfPath, string pdfList);
	~PDFclass();

	map<string, TH1*> pdfs;
	vector<string> pdfNames;

	double samplePDF(string pdfName);
        pair<double, double> sampleOrientPDF(const string pdfName);
	//double* sample2dPDF(PDFenum pdf, double sample);
	//double* sample3dPDF(PDFenum pdf, double sample);

   private:
	TFile* pdfFile;
	TRandom3 rand;

	bool appliedJacobian;
};



#endif

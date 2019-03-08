#ifndef PLOTCLASS_H
#define PLOTCLASS_H


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <TPaveText.h>
#include <TStyle.h>
#include <TColor.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPad.h>

//HomeGrown
#include "constants.h"
//#include "tools.h"

using namespace std;

enum PLOToptions {xLabel,yLabel,zLabel,xBinSpan,yBinSpan,zBinSpan,xSpan,ySpan,zSpan,markerStyle,markerColor,markerSize,maximum,minimum,logx,logy,logz,pads,draw,fileType};

class PLOTclass {

   public:
	PLOTclass(std::string canvasName = "MyC");
	~PLOTclass();

	TCanvas* MyC;
	
	void setStyle();
 	void setFolder(string fld);

	TH1* print1d(TH1* h, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH1* print1d(TH1* h, string name, PLOToptions opt, string val, string canv="orig");
	TH1* print1d(TH1* h, string name, string canv="orig");
	void print1d(vector<TH1*> hists, string name, string canv="orig");
	void print1d(vector<TH1*> hists, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH1* print1d(vector<double> data, string name, string canv="orig");
	TH1* print1d(vector<double> data, string name, PLOToptions opts, string vals, string canv="orig");
	TH1* print1d(vector<double> data, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* print2d(TH2* h, string name, PLOToptions opt, string val, string canv="orig");
	TH2* print2d(TH2* h, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* print2d(TH2* h, string name, string option, string canv="orig");
	void print2d(vector<TH2*> vect, PLOToptions opt, string val, string canv="orig");
	TH2* print2d(TH2* h, string name, string canv="orig");
	void print2d(vector<TH2*> vect, string canv="orig");
	TH2* printRC(vector< vector<double> > data, string name, string canv="orig");
	TH2* printRC(vector< vector<double> > data, string name, PLOToptions opt, string val, string canv="orig");
	TH2* printRC(vector< vector<double> > data, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* printRC(vector< vector<double> > data, double* rbins, double cMin, double cMax, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* printXY(vector< vector<double> > data, string name, string canv="orig");
	TH2* printXY(vector< vector<double> > data, string name, PLOToptions opt, string val,  string canv="orig");
	TH2* printXY(vector< vector<double> > data, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* printXY(vector< vector<double> > data, double* xbins, double yMin, double yMax, string name, vector<PLOToptions> opts, vector<string> vals, string canv="orig");
	TH2* printXY(map<int, vector<double> > &data, string name, string canv="orig");
	TH2* printXY(map<int, vector<double> > &data, string name, string xtitle, string ytitle, string canv="orig");
	void print(vector<TH1*> hists, string name, PLOToptions opt, string canv="orig");
	TH1* print(TH1* hist, string name, string canv);
	TH2* print(TH2* hist, string name, string canv);
	void printPads(vector<TH1*> hists, string name, string canv);


	TH1F* plot1d(vector<double> data, string name);
	TH1F* plot1d(vector<double> data, string name, PLOToptions opt, string val);
	TH1F* plot1d(vector<double> data, string name, vector<PLOToptions> opts, vector<string> vals);
	TH2F* plot2d(vector< vector<double> > data, string name, bool rowCol=true);
	TH2F* plotRC(vector< vector<double> > data, string name);
	TH2F* plotRC(vector< vector<double> > data, string name, PLOToptions opt, string val); 
	TH2F* plotRC(vector< vector<double> > data, string name, vector<PLOToptions> opts, vector<string> vals); 
	TH2F* plotRC(vector< vector<double> > data, string name, double xbins[], double ybins[], vector<PLOToptions> opts, vector<string> vals);
	TH2F* plotRC(vector< vector<double> > data, string name, double xbins[], double ybins[]);
	TH2F* plotRC(vector< vector<double> > data, double* rbins, double cMin, double cMax, string name, vector<PLOToptions> opts, vector<string> vals);
	TH2F* plotXY(vector< vector<double> > data, string name);
	TH2F* plotXY(vector< vector<double> > data, string name, vector<PLOToptions> opts, vector<string> vals); 
	TH2F* plotXY(vector< vector<double> > data, double* xbins, double yMin, double yMax, string name, vector<PLOToptions> opts, vector<string> vals); 
	TH2F* plotXY(map<int, vector<double> > &data, string name);
	TH2F* plotXY(map<int, vector<double> > &data, string name, vector<PLOToptions> opts, vector<string> vals);

	void doubleCanvas(string name, int width, int height, char type, float perc);
	void executeCommand(TH1* h, PLOToptions opt, string inp);
	void executeCommand(TH1* h, vector<PLOToptions> opts, vector<string> inps);
	void setAxisSpan(TAxis* h, string range);
	void setAxisBinSpan(TAxis* h, string range);
	void drawHist(TH1* h, string opt);

   private:
	bool DEFAULT;
	string folder;

	bool setlogx, setlogy, setlogz;

	struct canvas {
	    public:
		~canvas();

		TCanvas* can;
		vector<TPad*> pads;
	};
	
	string fType;
	TStyle* ATLASstyle();
	map<string, canvas*> CANVS;
};



#endif

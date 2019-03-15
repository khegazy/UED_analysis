#include "plotClass.h"

using namespace std;


PLOTclass::PLOTclass(std::string canvasName) {

  setStyle();
 //gStyle->SetPalette(kRainBow);
  DEFAULT = true;
  folder = "";
  fType = ".png";
  setlogx = setlogy = setlogz = false;
  CANVS["orig"] = new canvas();
  CANVS["orig"]->can = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
  MyC = CANVS["orig"]->can;
}

PLOTclass::~PLOTclass() { 

std::cout<<"Start Plot dec"<<std::endl;
  for (std::map<std::string, canvas*>::iterator it=CANVS.begin(); it!=CANVS.end(); it++) delete it->second; 
std::cout<<"End Plot dec"<<std::endl;
}


PLOTclass::canvas::~canvas() {

std::cout<<"Start Canv dec"<<std::endl;
  delete can;
std::cout<<"Start Canv dec"<<std::endl;
  for (auto& ip : pads) {std::cout<<ip<<std::endl; delete ip;}
std::cout<<"Start Canv dec"<<std::endl;
}


void PLOTclass::setStyle() {

  ATLASstyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();

  return;
}


void PLOTclass::setFolder(std::string fld) {

  folder = fld+"/";
  return;
}


TH1F* PLOTclass::plot1d(
    std::vector<double> data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint k;
  TH1F* h = new TH1F(
      name.c_str(), name.c_str(), 
      data.size(), 1, data.size());
  for (k=0; k<data.size(); k++) {
    h->SetBinContent(k+1, data[k]);
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(h, opts[k], vals[k]);
  }

  return h;
}


TH1F* PLOTclass::plot1d(
    std::vector<double> data, 
    std::string name, 
    PLOToptions opt, 
    std::string val) {

  TH1F* h = new TH1F(name.c_str(), name.c_str(), data.size(), 1, data.size());
  for (uint k=0; k<data.size(); k++) {
    h->SetBinContent(k+1, data[k]);
  }
  executeCommand(h, opt, val);

  return h;
}


TH1F* PLOTclass::plot1d(
    std::vector<double> data, 
    std::string name) {

  uint k;
  TH1F* h = new TH1F(
      name.c_str(), name.c_str(), 
      data.size(), 1, data.size());
  for (k=0; k<data.size(); k++) {
    h->SetBinContent(k+1, data[k]);
  }

  return h;
}


TH1* PLOTclass::print1d(
    TH1* h, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  CANVS[canv]->can->cd();
  int id=-1;
  for (int k=0; k<(int)opts.size(); k++) {
    if (opts[k] == draw) {
      id=k; 
      continue;
    }
    executeCommand(h, opts[k], vals[k]);
  }
  if (id == -1) h->Draw();
  else executeCommand(h, draw, vals[id]);

  return print(h, name, canv);
}


TH1* PLOTclass::print1d(
    TH1* h, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  CANVS[canv]->can->cd();
  executeCommand(h, opt, val);
  if (opt != draw) h->Draw();

  return print(h, name, canv);
}


TH1* PLOTclass::print1d(
    TH1* h, 
    std::string name, 
    std::string canv) {

  return print1d(h, name, draw, "", canv);
}


void PLOTclass::print1d(
    std::vector<TH1*> hists, 
    std::string name, 
    std::string canv) {

  CANVS[canv]->can->cd();
  for (uint ih=0; ih<hists.size(); ih++) {
    if (hists[ih]->GetMaximum() > hists[0]->GetMaximum()) {
      hists[0]->SetMaximum();
    }
  }

  hists[0]->Draw();
  for (uint ih=0; ih<hists.size(); ih++) hists[ih]->Draw("SAME");
  print(hists[0], name, canv); 
}


void PLOTclass::print1d(
    std::vector<TH1*> hists, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  CANVS[canv]->can->cd();
  for (uint ih=0; ih<hists.size(); ih++) {
    if (hists[ih]->GetMaximum() > hists[0]->GetMaximum()) {
      hists[0]->SetMaximum();
    }
  }

  for (int k=0; k<(int)opts.size(); k++) {
    if (opts[k] == draw) {
      cerr << "Implement Draw here/n"; 
      exit(0);
    }
    for (uint ih=0; ih<hists.size(); ih++) {
      executeCommand(hists[ih], opts[k], vals[k]);
    }
  }

  hists[0]->Draw();
  for (uint ih=0; ih<hists.size(); ih++) {
    hists[ih]->Draw("SAME");
  }
  print(hists[0], name, canv); 
}


void PLOTclass::print1d(
    std::vector<TH1*> hists, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  CANVS[canv]->can->cd();
  for (uint ih=0; ih<hists.size(); ih++) {
    if (hists[ih]->GetMaximum() > hists[0]->GetMaximum()) {
      hists[0]->SetMaximum();
    }
  }

  for (uint ih=0; ih<hists.size(); ih++) {
    executeCommand(hists[ih], opt, val);
  }

  hists[0]->Draw();
  for (uint ih=0; ih<hists.size(); ih++) {
    hists[ih]->Draw("SAME");
  }
  print(hists[0], name, canv); 
}


TH1* PLOTclass::print1d(
    std::vector<double> data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  return print1d(plot1d(data, name, opts, vals), name, draw, "", canv);
}


TH1* PLOTclass::print1d(
    std::vector<double> data, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  return print1d(plot1d(data, name, opt, val), name, draw, "", canv);
}


TH1* PLOTclass::print1d(
    std::vector<double> data, 
    std::string name, 
    std::string canv) {

  return print1d(plot1d(data, name), name, draw, "", canv);
}


TH2F* PLOTclass::plot2d(
    std::vector< std::vector<double> > data, 
    std::string name, 
    bool rowCol) {

  if (rowCol) return plotRC(data, name);
  return plotXY(data, name);
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    double xbins[], 
    double ybins[], 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint ir, ic, k;
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data[0].size(), xbins, 
      data.size(), ybins);
  for (ir=0; ir<data.size(); ir++) {
    for (ic=0; ic<data[0].size(); ic++) {
      hist->SetBinContent(ic+1, ir+1, data[ir][ic]);
    }
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  return hist;
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    double xbins[], 
    double ybins[]) {

  uint ir, ic;
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data[0].size(), xbins, 
      data.size(), ybins);
  for (ir=0; ir<data.size(); ir++) {
    for (ic=0; ic<data[0].size(); ic++) {
      hist->SetBinContent(ic+1, ir+1, data[ir][ic]);
    }
  }

  hist->SetContour(80);
  return hist;
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    PLOToptions opt, 
    std::string val) {

  uint ir, ic;
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data[0].size(), 1, 
      data[0].size(), data.size(), 
      1, data.size());
  for (ir=0; ir<data.size(); ir++) {
    for (ic=0; ic<data[0].size(); ic++) {
      hist->SetBinContent(ic+1, ir+1, data[ir][ic]);
    }
  }
  executeCommand(hist, opt, val);

  hist->SetContour(80);
  return hist;
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint ir, ic, k;
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data[0].size(), 1, 
      data[0].size(), data.size(), 
      1, data.size());
  for (ir=0; ir<data.size(); ir++) {
    for (ic=0; ic<data[0].size(); ic++) {
      hist->SetBinContent(ic+1, ir+1, data[ir][ic]);
    }
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  return hist;
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    std::string name) {

  
  std::vector<PLOToptions> o;
  std::vector<std::string> s;
  return plotRC(data, name, o, s);
}


TH2F* PLOTclass::plotRC(
    std::vector< std::vector<double> > data, 
    double* rbins, 
    double cMin, 
    double cMax, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint ic, ir, k;
  double* cbins = new double[data[0].size()+1];
  double cItr = (cMax-cMin)/((double)(data[0].size()+1));
  for (ic=0; ic<data[0].size()+1; ic++) {
    cbins[ic] = cMin + cItr*ic;
  }
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data[0].size(), cbins, 
      data.size(), rbins);
  for (ic=0; ic<data[0].size(); ic++) {
    for (ir=0; ir<data.size(); ir++) {
      hist->SetBinContent(ic+1, ir+1, data[ir][ic]);
    }
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  delete[] cbins;
  return hist;
}


TH2* PLOTclass::printRC(
    std::vector< std::vector<double> > data, 
    double* rbins, 
    double cMin, 
    double cMax, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  TH2F* h = plotRC(data, rbins, cMin, cMax, name, opts, vals);
  h->SetContour(80);
  print2d(h, name, opts, vals, canv);

  return h;
}


TH2F* PLOTclass::plotXY(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint ix, iy, k;
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data.size(), 1, 
      data.size(), data[0].size(), 
      1, data[0].size());
  for (ix=0; ix<data.size(); ix++) {
    for (iy=0; iy<data[0].size(); iy++) {
      hist->SetBinContent(ix+1, iy+1, data[ix][iy]);
    }
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  return hist;
}


TH2F* PLOTclass::plotXY(
    std::vector< std::vector<double> > data, 
    double* xbins, 
    double yMin, 
    double yMax, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  uint ix, iy, k;
  double* ybins = new double[data[0].size()+1];
  double yItr = (yMin-yMax)/((double)(data[0].size()+1));
  for (iy=0; iy<data[0].size()+1; iy++) {
    ybins[iy] = yMin + yItr*iy;
  }
  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data.size(), xbins, 
      data[0].size(), ybins);
  for (ix=0; ix<data.size(); ix++) {
    for (iy=0; iy<data[0].size(); iy++) {
      hist->SetBinContent(ix+1, iy+1, data[ix][iy]);
    }
  }
  for (k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  delete[] ybins;
  return hist;
}

TH2F* PLOTclass::plotXY(
    std::vector< std::vector<double> > data, 
    std::string name) {

  std::vector<PLOToptions> o;
  std::vector<std::string> s;
  return plotXY(data, name, o, s);
}


TH2F* PLOTclass::plotXY(
    std::map< int, std::vector<double> > &data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals) {

  TH2F* hist = new TH2F(
      name.c_str(), name.c_str(), 
      data.size(), 0, 
      data.size(), data[0].size(), 
      0, data[0].size());

  int ix=0;
  for (std::map<int, std::vector<double> >::iterator m_it=data.begin(); 
      m_it!=data.end(); 
      ++m_it) {

    for (uint iy=0; iy<m_it->second.size(); iy++) {
      hist->SetBinContent(ix+1, iy+1, m_it->second[iy]);
    }
    ++ix;
  }
  for (uint k=0; k<opts.size(); k++) {
    executeCommand(hist, opts[k], vals[k]);
  }

  hist->SetContour(80);
  return hist;
}


TH2* PLOTclass::printRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  TH2F* h = plotRC(data, name);
  print2d(h, name, opt, val, canv);

  return h;
}


TH2* PLOTclass::printRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  TH2F* h = plotRC(data, name);
  print2d(h, name, opts, vals, canv);

  return h;
}


TH2* PLOTclass::printRC(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::string canv) {

  TH2F* h = plotRC(data, name);
  print2d(h, name, draw, "COLZ", canv);

  return h;
}

 
TH2* PLOTclass::printXY(
    std::vector< std::vector<double> > data, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  TH2F* h = plotXY(data, name);
  h->SetContour(80);
  print2d(h, name, opt, val, canv);

  return h;
}


TH2* PLOTclass::printXY(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  TH2F* h = plotXY(data, name);
  h->SetContour(80);
  print2d(h, name, opts, vals, canv);

  return h;
}


TH2* PLOTclass::printXY(
    std::vector< std::vector<double> > data, 
    double* xbins, 
    double yMin, 
    double yMax, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  TH2F* h = plotXY(
      data, xbins, 
      yMin, yMax, 
      name, opts, vals);
  h->SetContour(80);
  print2d(h, name, opts, vals, canv);

  return h;
}


TH2* PLOTclass::printXY(
    std::vector< std::vector<double> > data, 
    std::string name, 
    std::string canv) {

  TH2F* h = plotXY(data, name);
  h->SetContour(80);
  print2d(h, name, draw, "COLZ", canv);

  return h;
}


TH2* PLOTclass::printXY(
    std::map<int, std::vector<double> > &data, 
    std::string name, 
    std::string canv) {

  std::vector<PLOToptions> o;  
  std::vector<std::string> s;
  TH2F* h = plotXY(data, name, o, s);
  h->SetContour(80);
  print2d(h, name, draw, "COLZ", canv);

  return h;
}


TH2* PLOTclass::print2d(
    TH2* h, 
    std::string name, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> vals, 
    std::string canv) {

  int id=-1;

  h->SetContour(80);
  CANVS[canv]->can->cd();
  for (int k=0; k<(int)opts.size(); k++) {
    if (opts[k] == draw) {id=k; continue;}
    executeCommand(h, opts[k], vals[k]);
  }
  if (id == -1) h->Draw("COLZ");
  else executeCommand(h, draw, vals[id]);

  return print(h, name, canv);
}


TH2* PLOTclass::print2d(
    TH2* h, 
    std::string name, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  h->SetContour(80);
  CANVS[canv]->can->cd();
  executeCommand(h, opt, val);
  if (opt != draw) h->Draw("COLZ");

  return print(h, name, canv);
}


void PLOTclass::print2d(
    std::vector<TH2*> vect, 
    PLOToptions opt, 
    std::string val, 
    std::string canv) {

  for (uint k=0; k<vect.size(); k++) {
    print2d(vect[k], std::string(vect[k]->GetName()), opt, val, canv);
  }
}


TH2* PLOTclass::print2d(
    TH2* h, 
    std::string name, 
    std::string canv) {

  return print2d(h, name, draw, "COLZ", canv);
}


void PLOTclass::print2d(
    std::vector<TH2*> vect, 
    std::string canv) {

  print2d(vect, draw, "COLZ", canv);
}


TH1* PLOTclass::print(
    TH1* hist, 
    std::string name, 
    std::string canv) {

  if (DEFAULT) {
    CANVS[canv]->can->Print((folder+name+fType).c_str());
  }
  else {
    if (setlogx) {
      CANVS[canv]->can->SetLogx(1);
      setlogx = false;
    }
    if (setlogy) {
      CANVS[canv]->can->SetLogy(1);
      setlogy = false;
    }
    if (setlogz) {
      CANVS[canv]->can->SetLogz(1);
      setlogz = false;
    }

    if (fType.find("ascii")!=std::string::npos 
        || fType.find("ASCII")!=std::string::npos) {

      ofstream file;
      file.open(name+fType);

      for (int ix=0; ix<hist->GetNbinsX(); ix++) {
        file<<hist->GetBinContent(ix+1)<<"\t";
      }

      file.close();
    }
    else CANVS[canv]->can->Print((folder+name+fType).c_str());

    fType=".png";
    CANVS[canv]->can->SetLogx(0);
    CANVS[canv]->can->SetLogy(0);
    CANVS[canv]->can->SetLogz(0);

    DEFAULT=true;
  }

  return hist;
}


TH2* PLOTclass::print(
    TH2* hist, 
    std::string name, 
    std::string canv) {

  if (DEFAULT) {
    CANVS[canv]->can->Print((folder+name+fType).c_str());
  }
  else {
    if (setlogx) {
      CANVS[canv]->can->SetLogx(1);
      setlogx = false;
    }
    if (setlogy) {
      CANVS[canv]->can->SetLogy(1);
      setlogy = false;
    }
    if (setlogz) {
      CANVS[canv]->can->SetLogz(1);
      setlogz = false;
    }

    if (fType.find("ascii")!=std::string::npos 
        || fType.find("ASCII")!=std::string::npos) {

      ofstream file;
      file.open(name+fType);

      int ix, iy;
      for (iy=0; iy<hist->GetNbinsY(); iy++) {
        for (ix=0; ix<hist->GetNbinsX(); ix++) {
          file<<hist->GetBinContent(ix+1,iy+1)<<"\t";
        }
        file<<"\n";
      }

      file.close();
    }
    else CANVS[canv]->can->Print((folder+name+fType).c_str());

    fType=".png";
    CANVS[canv]->can->SetLogx(0);
    CANVS[canv]->can->SetLogy(0);
    CANVS[canv]->can->SetLogz(0);

    DEFAULT=true;
  }

  return hist;
}


void PLOTclass::print(
    std::vector<TH1*> hists, 
    std::string name, 
    PLOToptions opt, 
    std::string canv) {

  switch (opt) {
    case pads: printPads(hists, name, canv);  return;
    default:
      cerr << "WARNING: Cannot find plotting type for enum "
        << opt << " !!!\n";
      return;
  }
}


void PLOTclass::printPads(
    std::vector<TH1*> hists, 
    std::string name, 
    std::string canv) {

std::cout<<"starting pads"<<std::endl;
  CANVS[canv]->can->cd();
  CANVS[canv]->pads[0]->Draw();
std::cout<<"changed canv"<<std::endl;
  for (uint k=0; k<hists.size(); k++) {
  CANVS[canv]->pads[0]->cd();
std::cout<<"111"<<std::endl;
    CANVS[canv]->pads[k]->Draw();
std::cout<<"222"<<std::endl;
    CANVS[canv]->pads[k]->cd();
std::cout<<"333"<<std::endl;
    hists[k]->Draw();
    CANVS[canv]->pads[k]->Update();
std::cout<<"444"<<std::endl;
  }
std::cout<<"555"<<std::endl;
  CANVS[canv]->pads[0]->cd();
  CANVS[canv]->pads[0]->Update();
std::cout<<"printing"<<std::endl;
    CANVS[canv]->pads[0]->SaveAs((name+fType).c_str());
  //CANVS[canv]->can->SaveAs((name+fType).c_str());
std::cout<<"changing back canv"<<std::endl;
  CANVS["orig"]->can->cd();
}
  
  

////////////////////////////////////
/////  Histogram manipulation  /////
////////////////////////////////////

void PLOTclass::setAxisSpan(TAxis* ax, std::string range) {

  std::size_t spot;
  spot = range.find(", ");
  if (spot != std::string::npos) {
    ax->SetLimits(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+2,(range.length()-(spot+2)))));
    return;
  }
  spot = range.find(",");
  if (spot != std::string::npos) {
    ax->SetLimits(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+1,(range.length()-(spot+1)))));
    return;
  }
  spot = range.find(" ");
  if (spot != std::string::npos) {
    ax->SetLimits(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+1,(range.length()-(spot+1)))));
    return;
  }

  cerr << "WARNING: Did not have the correct deliniation between min and max ("
    << range << "), must use ',' ' ' ', '!!!\n";
}


void PLOTclass::setAxisBinSpan(TAxis* ax, std::string range) {

  std::size_t spot;
  spot = range.find(", ");
  if (spot != std::string::npos) {
    ax->SetRange(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+2,(range.length()-(spot+2)))));
    return;
  }
  spot = range.find(",");
  if (spot != std::string::npos) {
    ax->SetRange(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+1,(range.length()-(spot+1)))));
    return;
  }
  spot = range.find(" ");
  if (spot != std::string::npos) {
    ax->SetRange(
        stof(range.substr(0,spot)), 
        stof(range.substr(spot+1,(range.length()-(spot+1)))));
    return;
  }

  cerr << "WARNING: Did not have the correct deliniation between min and max ("
    << range << "), must use ',' ' ' ', '!!!\n";
}


void PLOTclass::doubleCanvas(
    std::string name, 
    int width, int height, 
    char type, float perc) {

  CANVS[name] = new canvas();
  CANVS[name]->can = new TCanvas(name.c_str(), name.c_str(), width, height);
  CANVS[name]->pads.resize(3);
  CANVS[name]->pads[0] = new TPad("main", "main", 0.0, 0.0, 1.0, 1.0);
  if (type == 'x' || type == 'X') {
    CANVS[name]->pads[0] = new TPad(
        (name+"pad1").c_str(), 
        (name+"pad1").c_str(), 
        0.0, perc, 1.0, 1.0);
    CANVS[name]->pads[1] = new TPad(
        (name+"pad2").c_str(), 
        (name+"pad2").c_str(), 
        0.0, 0.0, 1.0, perc);
  }
  if (type == 'y' || type == 'Y') {
    CANVS[name]->pads[0] = new TPad(
        (name+"pad1").c_str(), 
        (name+"pad1").c_str(), 
        perc, 0.0, 1.0, 1.0);
    CANVS[name]->pads[1] = new TPad(
        (name+"pad2").c_str(), 
        (name+"pad2").c_str(), 
        0.0, 0.0, perc, 1.0);
  }
}


void PLOTclass::executeCommand(TH1* h, PLOToptions opt, std::string inp) {

  if (inp.compare("null") == 0) return;

  switch(opt) {
    case xLabel: 	h->GetXaxis()->SetTitle(inp.c_str());	return;
    case yLabel: 	h->GetYaxis()->SetTitle(inp.c_str());	return;
    case zLabel: 	h->GetZaxis()->SetTitle(inp.c_str());	return;
    case xSpan: 	setAxisSpan(h->GetXaxis(), inp);   	return;
    case ySpan: 	setAxisSpan(h->GetYaxis(), inp);   	return;
    case zSpan: 	setAxisSpan(h->GetZaxis(), inp);   	return;
    case xBinSpan: 	setAxisBinSpan(h->GetXaxis(), inp);   	return;
    case yBinSpan: 	setAxisBinSpan(h->GetYaxis(), inp);   	return;
    case zBinSpan: 	setAxisBinSpan(h->GetZaxis(), inp);   	return;
    case markerStyle: 	h->SetMarkerStyle(stoi(inp));  		return; 
    case markerColor: 	h->SetMarkerColor(stoi(inp));   	return;
    case markerSize: 	h->SetMarkerSize(stof(inp));   		return;
    case maximum: 	h->SetMaximum(stof(inp));   		return;
    case minimum: 	h->SetMinimum(stof(inp));   		return;
    case logx: 		setlogx=true; DEFAULT=false;   		return; 
    case logy: 		setlogy=true; DEFAULT=false;   		return; 
    case logz: 		setlogz=true; DEFAULT=false;   		return; 
    case draw: 		h->Draw(inp.c_str());			return;
    case fileType: 	fType="."+inp; DEFAULT=false; 		return;
    default: 
      cerr << "WARNING: Cannot find plotting function for enum " 
        << opt << " !!!\n";
      return;
  }
}


void PLOTclass::executeCommand(
    TH1* h, 
    std::vector<PLOToptions> opts, 
    std::vector<std::string> inps) {

  for (uint iopt=0; iopt<opts.size(); iopt++) {
    executeCommand(h, opts[iopt], inps[iopt]);
  }
}


//////////////////////////////
/////  ATLAS Plot Style  /////
//////////////////////////////

TStyle* PLOTclass::ATLASstyle() {

  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
    Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  //atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadTopMargin(0.04);
  //atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadRightMargin(0.15);
  atlasStyle->SetPadBottomMargin(0.14);
  atlasStyle->SetPadLeftMargin(0.14);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.4);
  atlasStyle->SetTitleYOffset(1.4);
  
  // use large fonts
  // Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.04;
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");

  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth(2.);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  atlasStyle->SetPalette(55);
  //atlasStyle->SetPalette(kDarkBodyRadiator);
  //atlasStyle->SetPalette(kVisibleSpectrum);

  return atlasStyle;
}




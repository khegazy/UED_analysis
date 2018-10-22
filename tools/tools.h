#ifndef TOOLS_H
#define TOOLS_H

#include <cmath>
#include <ctime>
#include <stdlib.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <sys/stat.h>

//From ROOT
#include <TH1.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>

//From Eigen
#include <Sparse>
#include <Core>
#include <Eigenvalues>

//From fftw
#include <fftw3.h>

//From Boost
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel.hpp>

//Home Grown
#include "constants.h"
#include "plotClass.h"

using namespace std;


namespace tools {


  inline bool fileExists(const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
  }

  inline int round(double a) {
    return fmod(a,(int)a) >= 0.5 ? (int)(a+1) : (int)a;
  }
  inline void nrShift2(double& a, double& b, const double c) {
    a=b; b=c;
  }
  inline void nrShift3(double& a, double& b, double& c, const double d) {
    a=b; b=c; c=d;
  }
  inline void nrShift3(vector<double>& a, vector<double>& b, vector<double>& c, const vector<double> d) {
    a=b; b=c; c=d;
  }
  inline void nrSign(double& a, double& b) {
   a = b >= 0.0 ? fabs(a) : -fabs(a); 
  }
  inline void swap(double& a, double& b) {
    double c=a; a=b; b=c;
  }   
  inline void swap(vector<double>& a, vector<double>& b) {
    vector<double> c=a; a=b; b=c;
  }   

  int stringRemove(std::string inp, char remove, int ind);
  
  // Matrix/Vector manipulation
  double vNorm(vector<double>& vect);
  double mag(vector<double>& vect);
  double vDotP(vector<double>& v1, vector<double> v2);
  void gramSchmidt(vector< vector<double>* > &gsOrth, bool normalize);
  Eigen::MatrixXd SVDinvert(Eigen::MatrixXd &mat);

  template <typename T>
  void vScale(std::vector<T> &vec, T scale);

  // Legendre
  double legendreIntegral(int order, double init, double fin);

  // FFT
  std::vector< std::vector<double> > fft1dRtoC(const std::vector<double> &input, 
      double outRat, int Npadding, double padDecayRat, 
      const fftw_plan &fftPlan, double* in, fftw_complex* out, 
      PLOTclass* pltVerbose=NULL, std::string pltName="fftInput"); 
  std::vector< std::vector<double> > fft1dRtoC(const std::vector<double> &input, 
      const fftw_plan &fftPlan, double* in, fftw_complex* out, 
      bool forwardFFT, bool centeredOrigin, 
      PLOTclass* pltVerbose=NULL, std::string pltName="fftInput"); 

  // Machine Learning
  Eigen::MatrixXd normalEquation(Eigen::MatrixXd X, Eigen::MatrixXd Y);
  std::pair<Eigen::MatrixXcd, Eigen::VectorXcd>
    PCA(Eigen::MatrixXd inpArrays, bool zeroMean = true, bool unitVariance = true);

  // Minimizing functions
  template <typename fType>
  double powellMin(fType fxn, vector<double>& p, double scale, double powellTol, double fracTol1d, bool verbose = false);
  template <typename fType>
  double powellMin(fType fxn, vector<double>& p, double scale, double minScale, double powellTol, double fracTol1d, bool verbose = false);
  template <typename fType>
  double powellMin(fType fxn, vector<double>& p, vector<double> scale, vector<double> minScale, vector<double> powellTol, double fracTol1d, bool verbose = false);
  template <typename fType>
  vector<double> goldenMin(fType fxn, vector<double>& ap, vector<double>& bp, double fracTol, double difTol, double& fmin, bool verbose);
  template <typename fType>
  vector<double> goldenMinSubR(fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double fracTol, double difTol, double& fmin, bool verbose);
  template <typename fType>
  void mnbrak(fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double& fa, double& fb, double& fc, bool verbose);
 
 
  template <class fType>
  class F1DIMclass {

    public:
	F1DIMclass(fType fxn_i, vector<double> slope_i, vector<double> p0_i);

 	double eval(double x);
 

   private:
    	fType fxn;
	vector<double> p0;
	vector<double> slope;
	int Nvar;
  };


  // Other
  double getAzmAngle(double X, double Y);
}






template <typename T>
void tools::vScale(std::vector<T> &vec, T scale) {
  std::for_each(vec.begin(), vec.end(), [scale](T &d) 
      { d *= scale;});
}











template <typename fType>
double tools::powellMin(fType fxn, vector<double>& p, double scale, double powellTol, double fracTol1d, bool verbose) {

  double minScale = scale/1000.0 > 1e-6 ? scale/1000.0 : 1e-6;
  return tools::powellMin<fType>(fxn, p, scale, minScale, powellTol, fracTol1d, verbose); 
}


template <typename fType>
double tools::powellMin(fType fxn, vector<double>& p, vector<double> scale, vector<double> minScale, vector<double> powellTol, double fracTol1d, bool verbose) {

  if (verbose) cout<<"ENTERING POWELLMIN"<<endl<<endl;
  const int ITMAX = 200;
  int Nvar = p.size();
  double fVal, fValP, fValPo, del, vmag;
  int iter, im, k, iBig;
  bool pass;
  vector<double> pP(Nvar), pE(Nvar), pf(Nvar), pi(Nvar), pdif(Nvar), avgDir(Nvar);
  vector< vector<double> > varDir;

  varDir.resize(Nvar);
  for (int i=0; i<Nvar; i++) {
    varDir[i].resize(Nvar);
    for (int k=0; k<Nvar; k++) varDir[i][k] = 0;
    varDir[i][i] = 1;
  }


  fVal = fxn(p);
  for (iter=0; ; ++iter) {
    fValPo=fVal;
    pP=p;
    del=0;
    iBig=0;

    // minimize along each variable
    vmag = tools::mag(p);
    for (im=0; im<Nvar; im++) {
      if (verbose) cout<<"minimizing var "<<im<<endl;
      fValP = fVal;
      for (k=0; k<Nvar; k++) {
        pf[k] = p[k] + scale[k]*varDir[im][k];
        pi[k] = p[k] - scale[k]*varDir[im][k];
      }
      p = goldenMin(fxn, pi, pf, fracTol1d, 0, fVal, verbose);
      if (fValP - fVal > del) {
        del = fValP - fVal;
      }
    }

    // Minimizing along the avgDir
    if (verbose) cout<<endl<<"minimizing along avgDir"<<endl;
    fValP = fVal;
    for (k=0; k<Nvar; k++) {
      avgDir[k] = p[k] - pP[k];
      pdif[k] = avgDir[k];
      pf[k] = p[k] + 0.5*avgDir[k];
      pi[k] = p[k];
    }
    if (tools::mag(pdif)) p = goldenMin(fxn, pi, pf, fracTol1d, 0, fVal, verbose);

    // Check to see if avgDir should be added (needs to be predominently one direction)
    if (2.0*(fValPo-2.0*fValP+fVal)*pow(fValPo-fValP-del, 2) < del*pow(fValPo-fVal, 2)) {
      if (verbose) cout<<"adding new vector to variables"<<endl;
      vmag = tools::mag(avgDir);
      for (k=0; k<Nvar; k++) {
        varDir[iBig][k] = varDir[Nvar-1][k];
        varDir[Nvar-1][k] = avgDir[k]/vmag;
      }
    }

    if (verbose) {
      cout<<endl<<"powell iter: "<<iter<<endl;
      cout<<"P0 P"<<endl;
      for (k=0; k<Nvar; k++) cout<<pP[k]<<"  "<<p[k]<<endl;
      cout<<endl<<"Checking powell condition (each var < powellTol): "<<endl;
      for (k=0; k<Nvar; k++) cout<<fabs(pP[k]-p[k])<<" < "<<powellTol[k]<<endl;
      cout<<endl;
    }

    // Check if close enough (Powell condition)
    //if ( 2.0*(fValPo - fVal) < tol*(fabs(fValPo) + fabs(fVal)) + TINY) return 1; //return p;
    //if (fabs(tools::mag(pP)-tools::mag(p)) < sqrt(Nvar)*powellTol ) break;
    pass = true;
    for (k=0; k<Nvar; k++) pass = pass && (fabs(pP[k]-p[k])<powellTol[k]);
    if (pass) break;
    if (iter > ITMAX) {
      cerr << "ERROR: Powell method reached the maximum iterations without satisfying tolerance!!!" << endl;
      exit(0);
    }

    // Recalculating the scale of the window around min guess
    for (k=0; k<Nvar; k++) {
      pdif[k] = p[k] - pP[k];
      scale[k] = fabs(pdif[k]/pi[k])*scale[k] > minScale[k] ? 
			fabs(pdif[k]/pi[k])*scale[k] : minScale[k];
      if (verbose) cout << "scale" << k << ": " << scale[k] << endl;
    }
  }

return fVal;
}





template <typename fType>
double tools::powellMin(fType fxn, vector<double>& p, double scale, double minScale, double powellTol, double fracTol1d, bool verbose) {

  if (verbose) cout<<"ENTERING POWELLMIN"<<endl<<endl;
  const int ITMAX = 200;
  int Nvar = p.size();
  double fVal, fValP, fValPo, del, vmag;
  int iter, im, k, iBig;
  bool pass;
  vector<double> pP(Nvar), pE(Nvar), pf(Nvar), pi(Nvar), pdif(Nvar), avgDir(Nvar);
  vector< vector<double> > varDir;

  varDir.resize(Nvar);
  for (int i=0; i<Nvar; i++) {
    varDir[i].resize(Nvar);
    for (int k=0; k<Nvar; k++) varDir[i][k] = 0;
    varDir[i][i] = 1;
  }


  fVal = fxn(p);
  for (iter=0; ; ++iter) {
    fValPo=fVal;
    pP=p;
    del=0;
    iBig=0;

    // minimize along each variable
    vmag = tools::mag(p);
    for (im=0; im<Nvar; im++) {
      if (verbose) cout<<"minimizing var "<<im<<endl;
      fValP = fVal;
      for (k=0; k<Nvar; k++) {
        pf[k] = p[k] + scale*varDir[im][k];
        pi[k] = p[k] - scale*varDir[im][k];
      }
      p = goldenMin(fxn, pi, pf, fracTol1d, powellTol, fVal, verbose);
      if (fValP - fVal > del) {
        del = fValP - fVal;
      }
    }

    // Minimizing along the avgDir
    if (verbose) cout<<endl<<"minimizing along avgDir"<<endl;
    fValP = fVal;
    for (k=0; k<Nvar; k++) {
      avgDir[k] = p[k] - pP[k];
      pdif[k] = avgDir[k];
      pf[k] = p[k] + 0.5*avgDir[k];
      pi[k] = p[k];
    }
    if (tools::mag(pdif)) p = goldenMin(fxn, pi, pf, fracTol1d, powellTol, fVal, verbose);

    // Check to see if avgDir should be added (needs to be predominently one direction)
    if (2.0*(fValPo-2.0*fValP+fVal)*pow(fValPo-fValP-del, 2) < del*pow(fValPo-fVal, 2)) {
      if (verbose) cout<<"adding new vector to variables"<<endl;
      vmag = tools::mag(avgDir);
      for (k=0; k<Nvar; k++) {
        varDir[iBig][k] = varDir[Nvar-1][k];
        varDir[Nvar-1][k] = avgDir[k]/vmag;
      }
    }

    if (verbose) {
      cout<<endl<<"powell iter: "<<iter<<endl;
      cout<<"P0 P"<<endl;
      for (k=0; k<Nvar; k++) cout<<pP[k]<<"  "<<p[k]<<endl;
      //cout<<endl<<"Checking powell condition: "<<fabs(tools::mag(pP)-tools::mag(p))<<" < "<<sqrt(Nvar)*powellTol<<endl<<endl;
      cout<<endl<<"Checking powell condition (each var < powellTol): "<<endl;
      for (k=0; k<Nvar; k++) cout<<fabs(pP[k]-p[k])<<" < "<<powellTol<<endl;
      cout<<endl;
    }

    // Check if close enough (Powell condition)
    //if ( 2.0*(fValPo - fVal) < tol*(fabs(fValPo) + fabs(fVal)) + TINY) return 1; //return p;
    //if (fabs(tools::mag(pP)-tools::mag(p)) < sqrt(Nvar)*powellTol ) break;
    pass = true;
    for (k=0; k<Nvar; k++) pass = pass && (fabs(pP[k]-p[k])<powellTol);
    if (pass) break;
    if (iter > ITMAX) {
      cerr<<"ERROR: Powell method reached the maximum iterations without satisfying tolerance!!!"<<endl;
      exit(0);
    }

    // Recalculating the scale of the window around min guess
    for (k=0; k<Nvar; k++) pdif[k] = p[k] - pi[k];
    vmag = tools::mag(pdif);
    scale = vmag/5.0 > minScale ? vmag/5.0 : minScale;
    for (k=0; k<Nvar; k++) pdif[k] = p[k] - pP[k];
    vmag = tools::mag(pdif);
    if (verbose) cout<<"scale: "<<scale<<endl;
  }

return fVal;
}





template <typename fType>
vector<double> tools::goldenMin(fType fxn, vector<double>& ap, vector<double>& bp, double fracTol, double difTol, double& fmin, bool verbose) {

  double fa = fxn(ap);
  double fb = fxn(bp);
  double fc;
  vector<double> cp(ap.size()), result(ap.size());

  if (verbose) cout<<"ENTERING GOLDENMIN"<<endl<<endl;
  tools::mnbrak(fxn, ap, bp, cp, fa, fb, fc, verbose);
  result = tools::goldenMinSubR(fxn, ap, bp, cp, fracTol, difTol, fmin, verbose);
  bp = cp;
  if (verbose) cout<<endl<<"EXITING GOLDENMIN"<<endl;
  return result;
}





template <typename fType>
void tools::mnbrak( fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double& fa, double& fb, double& fc, bool verbose) {
  ////// This function will find a minimum between the two endpoints. If ap and bp do not contain a minimum 
  //////        it will travel downhill until ap and cp contain a minimum. 
  // ap and bp are the endpoints of the initial guess (brackets) cp should be empty.
  // ap, bp, and cp will be altered so that the final result will be ap and cp contain
  //    a minimum where fxn(ap)>fxn(bp)<fxn(cp) and the point bp is between ap and cp.

 if (verbose) cout<<endl<<"ENTERING MNBRAK"<<endl<<endl;
  double GLIMIT = 200;
  double ulim,u,r,q,fu;
  int k;
  int Nvar = ap.size();
  const vector<double> p0 = ap;
  vector<double> slope(Nvar), pi(Nvar), pf(Nvar);
  double ax = 0.0;
  double bx = 1.0;

  for (k=0; k<Nvar; k++) slope[k] = bp[k]-ap[k];

  fa = fxn(ap);
  fb = fxn(bp);

  // Checking if initial guess contains minimum 
  //    cp is in the middle just for this check
  for (k=0; k<Nvar; k++) cp[k] = (ap[k]+bp[k])/2.0;
  fc = fxn(cp);
  if (fc < fa && fc < fb) {
    tools::swap(bp, cp);
    tools::swap(fb, fc);
    if (verbose) cout<<"EXITING MNBRAK: initial guess works"<<endl<<endl;
    return;
  }

  // Find new endpoints (ap, cp) to bracket a minimum
  if (fb > fa) {
    tools::swap(ap,bp);
    tools::swap(ax,bx);
    tools::swap(fb,fa);
  }

  tools::F1DIMclass<fType> f1d(fxn, slope, p0); // turning fxn into a 1d function given an initial point and slope
  double cx = bx+GOLD*(bx-ax);
  fc = f1d.eval(cx);

  if (verbose) {
    cout<<"ap  "<<"bp  "<<"cp  "<<"slope  "<<endl;
    for (k=0; k<Nvar; k++) cout<<ap[k]<<"  "<<bp[k]<<"  "<<cp[k]<<"  "<<slope[k]<<endl;
  }

  while (fb > fc) {
    if (verbose) cout<<"Parameters: "<<ax<<"  "<<bx<<"  "<<cx<<endl;
    r = (bx-ax)*(fb-fc);
    q = (bx-cx)*(fb-fa);
    u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*( (q-r)>=0.0 ? fabs(max(fabs(q-r),TINY)) : -fabs(max(fabs(q-r),TINY))));
    ulim = bx+GLIMIT*(cx-bx);

    if ((bx-u)*(u-cx) > 0.0) {
      fu = f1d.eval(u);
      if (fu < fc) {
        ax = bx;
        bx = u;
        fa = fb;
        fb = fu;
        break;
      }
      else if (fu > fb) {
        cx = u;
        fc = fu;
        break;
      }
      u = cx+GOLD*(cx-bx);
      fu = f1d.eval(u);
    }
    else if ((cx-u)*(u-ulim) > 0.0) {
      fu = f1d.eval(u);
      if (fu < fc) {
        tools::nrShift3(bx,cx,u,(u+GOLD*(u-cx)));
        cout<<f1d.eval(u)<<endl;
        tools::nrShift3(fb,fc,fu,f1d.eval(u));
      }
    }
    else if ((u-ulim)*(ulim-cx) >= 0.0) {
      u = ulim;
      fu = f1d.eval(u);
    }
    else {
      u = cx+GOLD*(cx-bx);
      fu = f1d.eval(u);
    }
    tools::nrShift3(ax,bx,cx,u);
    tools::nrShift3(fa,fb,fc,fu);
  }

  for (k=0; k<Nvar; k++) {
    ap[k] = p0[k] + ax*slope[k];
    bp[k] = p0[k] + bx*slope[k];
    cp[k] = p0[k] + cx*slope[k];
  }

  if (verbose) {
    cout<<endl<<"final vectors"<<endl;
    cout<<"parameters: "<<ax<<"  "<<bx<<"  "<<cx<<endl;
    cout<<"ap  "<<"bp  "<<"cp  "<<endl;
    for (k=0; k<Nvar; k++) cout<<ap[k]<<"  "<<bp[k]<<"  "<<cp[k]<<endl;
    cout<<endl<<"EXITING MNBRAK"<<endl<<endl;
  }

  return;
}





template <typename fType>
vector<double> tools::goldenMinSubR(fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double fracTol, double difTol, double& fmin, bool verbose) {
  ////// This function will find the minimum along the n-dimensional path between ap and cp using the golden ratio search
  //////        the exit condition can be triggered if the percent change in search signal fabs(x3-x0)/(fabs(x1)+fabs(x2)) 
  //////        is small or if the difference between the endpoints is small enough. Usually calling function will desire
  //////        difference condition as the final condition and keep lowering fracTol as the calling function converges on 
  //////        the minimum. This is because earlier search points do not need to find precise points until they are around min.
  //  REQUIREMENTS: ap and cp must be the endpoints of the search segment and contain bp
  //  fracTol: fractional tolerance, if fractional change in a search step is smaller then exit routine
  //  difTol: difference tolerance, if the difference between search endpoints is smaller, then exit routine.
  //  	difTol<2*minScale since difference will be 2*minScale near end.

  if (verbose) cout<<"STARTING GOLDENMINSUBR"<<endl<<endl;

  int Nvar = ap.size();
  int k;
  int iter = 0;
  double difMag = difTol + 1;
  double x0, x1, x2, x3, xb, f1, f2;
  const vector<double> p0 = ap;
  vector<double> pret(Nvar), slope(Nvar);
  int ITERM = 100;
  double GOLDN = GOLD-1.0;
  double GOLDC = 1.0-GOLDN;

  // calculating the slope and parameratizing fxn to be 1 dimensional
  for (k=0; k<Nvar; k++) slope[k] = cp[k] - ap[k];
  tools::F1DIMclass<fType> f1d(fxn, slope, p0);

  // Initial conditions
  x0 = 0.0;
  x3 = 1.0;
  xb = -99999;
  for (k=0; k<Nvar; k++) {
    if (fabs(cp[k]-ap[k])>0) {
      xb = (bp[k]-ap[k])/(cp[k]-ap[k]); 
      break;
    }
  }
  if (xb == -99999) cerr<<"WARNING: Cannot find a middle point between bracket points in goldenMinSubR!!!"<<endl;
  if (fabs(x3-xb) > fabs(x0-xb)) {
    x1 = xb;
    x2 = x1 + GOLDC*(x3-x1);
  }
  else {
    x2 = xb;
    x1 = x2 - GOLDC*(x2-x0);
  }
  f1 = f1d.eval(x1);
  f2 = f1d.eval(x2);

  while (fabs(x3-x0)/(fabs(x1)+fabs(x2)) > fracTol && difMag > difTol) {
    if (f2 < f1) {
      tools::nrShift3(x0, x1, x2, GOLDN*x2+GOLDC*x3);
      tools::nrShift2(f1,f2,f1d.eval(x2));
    }
    else {
      tools::nrShift3(x3, x2, x1, GOLDN*x1+GOLDC*x0);
      tools::nrShift2(f2,f1,f1d.eval(x1));
    }
    iter++;
    if (iter > ITERM || fabs(x0-x3) < 1e-7) {
      cerr<<"ERROR: iter > ITERMAX or cannot resolve further due to double precision in goldenSubR!!! "<<iter<<"  "<<fabs(x0-x3)<<endl;
      exit(0);
    }
    difMag=0;
    for (k=0; k<Nvar; k++) difMag += pow((x0-x3)*slope[k],2);
    difMag = sqrt(difMag);
    if (verbose) {
      cout<<"iter f1 f2: "<<iter<<"  "<<f1<<"  "<<f2<<endl;
      cout<<"Parameters: "<<x0<<"  "<<x1<<"  "<<x2<<"  "<<x3<<endl;
      cout<<"p0  "<<"p1  "<<"p2  "<<"p3  "<<endl;
      for (k=0; k<Nvar; k++) cout<<p0[k]+x0*slope[k]<<"  "<<p0[k]+x1*slope[k]<<"  "<<p0[k]+x2*slope[k]<<"  "<<p0[k]+x3*slope[k]<<endl;
      cout<<"Checking golden condition: "<<fabs(x0-x3)/(fabs(x1)+fabs(x2))<<" < "<<fracTol<<endl;
      cout<<"Checking golden condition: "<<difMag<<" < "<<difTol<<endl<<endl;
    }
  }

  // Updating vectors with final values
  for (k=0; k<Nvar; k++) {
    ap[k] = p0[k]+x0*slope[k];
    cp[k] = p0[k]+x3*slope[k];
  }

  // Choosing point with the smallest value to return
  if (f1 < f2) {
    fmin = f1;
    for (k=0; k<Nvar; k++) {
      bp[k] = p0[k]+x2*slope[k];
      pret[k] = p0[k]+x1*slope[k];
    }
  }
  else {
    fmin = f2;
    for (k=0; k<Nvar; k++) {
      bp[k] = p0[k]+x1*slope[k];
      pret[k] = p0[k]+x2*slope[k];
    }
  }

  if (verbose) cout<<"EXITING GOLDENMINSUBR"<<endl;
  return pret;
}






template <class fType>
tools::F1DIMclass<fType>::F1DIMclass(fType fxn_i, const vector<double> slope_i, const vector<double> p0_i) {

  fxn = fxn_i;
  slope = slope_i;
  p0 = p0_i;
  Nvar = p0.size();
}





template <class fType>
double tools::F1DIMclass<fType>::eval(double x) {

  vector<double> p(Nvar);
  for (int k=0; k<Nvar; k++) {
    p[k] = p0[k] + x*slope[k];
  }
  return fxn(p);
}



#endif

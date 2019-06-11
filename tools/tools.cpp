#include "tools.h"

using namespace std;


double tools::vNorm(vector<double>& vect) {

  double sqrsum = 0;
  for (uint k=0; k<vect.size(); k++) sqrsum += pow(vect[k],2);
  return sqrsum; 
}


double tools::mag(vector<double>& vect) {

  return sqrt(vNorm(vect));
}


double tools::vDotP(vector<double>& v1, vector<double> v2) {

  if (v1.size() != v2.size()) {
    std::cerr << "ERROR: Vectors in dot product do not have the same size!!!\n"; 
    exit(0);
  }
  double sum=0;
  for (uint k=0; k<v1.size(); k++) {
    sum += v1[k]*v2[k];
  }

  return sum;
}


int tools::stringRemove(std::string inp, char remove, int ind) {

  while (inp[ind] == remove) {
    ind++;
  }
  return ind;
}


void tools::gramSchmidt(vector< vector<double>* > &gsOrth, bool normalize) {

  uint l,j,i;
  double nrm, dp;
  for (l=0; l<gsOrth.size(); l++) {
    for (j=0; j<l; j++) {
      nrm = vNorm((*gsOrth[j]));
      dp = vDotP((*gsOrth[j]), (*gsOrth[l]));
      for (i=0; i<(*gsOrth[l]).size(); i++) {
        (*gsOrth[l])[i] -= (*gsOrth[j])[i]*dp/nrm;
      }
    }
  }

  if (normalize) {
    for (l=0; l<gsOrth.size(); l++) {
      nrm = tools::mag((*gsOrth[l]));
      for (i=0; i<(*gsOrth[l]).size(); i++) {
        (*gsOrth[l])[i] /= nrm;                        
      }
    }
  }

  return;
}


Eigen::MatrixXd tools::SVDinvert(Eigen::MatrixXd &mat) {

  uint svSize = std::min(mat.rows(), mat.cols());
  Eigen::SparseMatrix<double> svMat(svSize, svSize);
  Eigen::BDCSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

  for (uint i=0; i<svSize; i++) {
    svMat.coeffRef(i,i) = 1./svd.singularValues()(i);
  } 
  
  return svd.matrixV()*svMat*svd.matrixU().transpose();
}


Eigen::MatrixXd tools::normalEquation(
    Eigen::MatrixXd X, 
    Eigen::MatrixXd Y) {

  auto w = Eigen::VectorXd::Ones(X.rows());

  return normalEquation(X, Y, w);
}


Eigen::MatrixXd tools::normalEquation(
    Eigen::MatrixXd X, 
    Eigen::MatrixXd Y, 
    Eigen::VectorXd w) {

  auto W = w.asDiagonal();
  Eigen::MatrixXd nrmInp = X.transpose()*W*X;
  Eigen::MatrixXd norm = tools::SVDinvert(nrmInp);

  return norm*(X.transpose()*W*Y);
}


std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> 
  tools::PCA(Eigen::MatrixXd inpArrays, bool zeroMean, bool unitVariance) {


  if (zeroMean) {
    inpArrays.colwise() -= inpArrays.rowwise().mean();
  }
  // FIX
  if (unitVariance) {
    cerr << "ERROR: Must write variance normalization for PCA!" << endl;
    exit(0);
  }

  /////  Make covariance matrix  /////
  Eigen::MatrixXd covMat = (inpArrays.row(0)).transpose()*inpArrays.row(0);

  for (int ir=1; ir<inpArrays.rows(); ir++) {
    covMat += (inpArrays.row(ir)).transpose()*inpArrays.row(ir);
  }

  // normalize
  covMat /= inpArrays.rows();

  /////  Solve for Eigenvectors  /////
  Eigen::EigenSolver<Eigen::MatrixXd> solver(covMat);

  cout<<"eigenvectors pca"<<endl<<solver.eigenvectors()<<endl<<endl;
  cout<<"eigenvalues pca"<<endl<<solver.eigenvalues()<<endl<<endl;

  std::pair<Eigen::MatrixXcd, Eigen::VectorXcd> 
    results(solver.eigenvectors(), solver.eigenvalues()); 
          
  return results;
}

std::vector< std::vector<double> > tools::pythonFit(
    std::vector< std::vector<double> >* X,
    std::vector<double>* Y,
    std::vector<double>* W,
    std::vector<double>* p0,
    bool globalMin,
    std::vector< std::vector<double> >* searchRange) {

  std::string XfileName = "./fitValuesY_Bins["
      + to_string(X->size()) 
      + "," + to_string((*X)[0].size()) + "].dat";
  save::saveDat<double>((*X), XfileName);
  std::string YfileName = "./fitValuesY_Bins["
      + to_string(Y->size()) + "].dat";
  save::saveDat<double>((*Y), YfileName);

  std::string WfileName = "./fitValuesW_Bins["
      + to_string(W->size()) + "].dat";
  if (W != NULL) {
    save::saveDat<double>((*W), WfileName);
  }

  std::string p0RangeFileName;
  if (globalMin) {
    p0RangeFileName += "./globalMinP0ranges["
      + to_string(searchRange->size())
      + "," + to_string((*searchRange)[0].size());
    save::saveDat<double>((*searchRange), p0RangeFileName);
  }


  std::string codeDir = "/reg/neh/home/khegazy/baseTools/tools/";
  std::string command = "python " + codeDir
        + "fitLinCurveFit.py --X " + XfileName
        + " --Y " + YfileName;
  if (W != NULL) {
    command += " --W " + WfileName;
  }
  if (globalMin) {
    command += " --globalMin " + p0RangeFileName;
  }
  
  std::vector< std::vector<double> > fitCoeffs(2);
  fitCoeffs[0].resize(X->size());
  fitCoeffs[1].resize(X->size());
  std::string coeffsFileName = "./fitCoefficients_Bins["
      + to_string(X->size()) + "].dat";
  std::string coeffsErrFileName = "./fitCoefficientErrors_Bins["
      + to_string(X->size()) + "].dat";

  save::importDat<double>(fitCoeffs[0], coeffsFileName);
  save::importDat<double>(fitCoeffs[1], coeffsErrFileName);

  if (W != NULL) {
    std::remove(WfileName.c_str());
  }
  if (globalMin) {
    std::remove(p0RangeFileName.c_str());
  }
  std::remove(YfileName.c_str());
  std::remove(XfileName.c_str());
  std::remove(coeffsFileName.c_str());
  std::remove(coeffsErrFileName.c_str());

  return fitCoeffs;
}
   


double tools::getAzmAngle(double X, double Y) {
  double ang = atan(Y/X);
  if ((Y >= 0) && (X <= 0)) ang += PI;
  if ((Y < 0)  && (X < 0))  ang += PI;
  if ((Y <= 0) && (X >= 0)) ang += 2*PI;

  return ang;
}



/*
template <typename fType>
double tools::powellMin(fType fxn, vector<double>& p, double scale, double powellTol, double difTol1d, double fracTol1d, bool verbose) {

  if (verbose) cout<<"ENTERING POWELLMIN"<<endl<<endl;
  const int ITMAX = 200;
  int Nvar = p.size();
  double minScale = scale/1000.0 > 1e-6 ? scale/1000.0 : 1e-6;
  double fVal, fValP, fValPo, del, vmag;
  int iter, im, k, iBig;
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
      p = goldenMin(fxn, pi, pf, difTol1d, fracTol1d, fVal, verbose);
      if (fValP - fVal > del) {
        del = fValP - fVal;
      }
    }

    // Minimizing along the avgDir
    if (verbose) cout<<endl<<"minimizing along avgDir"<<endl;
    fValP = fVal;
    for (k=0; k<Nvar; k++) {
      avgDir[k] = p[k] - pP[k];
      pf[k] = p[k] + 0.5*avgDir[k];
      pi[k] = p[k];
    }
    p = goldenMin(fxn, pi, pf, difTol1d, fracTol1d, fVal, verbose);

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
      cout<<endl<<"Checking powell condition: "<<fabs(tools::mag(pP)-tools::mag(p))<<" < "<<sqrt(Nvar)*powellTol<<endl<<endl;
    }

    // Check if close enough (Powell condition)
    //if ( 2.0*(fValPo - fVal) < tol*(fabs(fValPo) + fabs(fVal)) + TINY) return 1; //return p;
    if (fabs(tools::mag(pP)-tools::mag(p)) < sqrt(Nvar)*powellTol ) break;
    if (iter > ITMAX) {
      cerr<<"ERROR: Powell method reached the maximum iterations without satisfying tolerance!!!"<<endl;
      exit(0);
    }

    // Recalculating the scale of the window around min guess
    for (k=0; k<Nvar; k++) pdif[k] = p[k] - pi[k];
    vmag = tools::mag(pdif);
    scale = 3.0*scale*vmag > minScale ? 3.0*scale*vmag : minScale;
    for (k=0; k<Nvar; k++) pdif[k] = p[k] - pP[k];
    vmag = tools::mag(pdif);
    fracTol1d = fracTol1d*vmag/tools::mag(pi);
  }

return fVal;
}
*/





/*
template <typename fType>
vector<double> tools::goldenMin(fType fxn, vector<double>& ap, vector<double>& bp, double difTol, double fracTol, double& fmin, bool verbose) {

  double fa = fxn(ap);
  double fb = fxn(bp);
  double fc;
  vector<double> cp(ap.size()), result(ap.size());

  if (verbose) cout<<"ENTERING GOLDENMIN"<<endl<<endl;
  tools::mnbrak(fxn, ap, bp, cp, fa, fb, fc, verbose);
  result = tools::goldenMinSubR(fxn, ap, bp, cp, difTol, fracTol, fmin, verbose);
  bp = cp;
  if (verbose) cout<<endl<<"EXITING GOLDENMIN"<<endl;
  return result;
}





template <typename fType>
void tools::mnbrak( fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double& fa, double& fb, double& fc, bool verbose) {
  ////// This function will find a minimum between the two endpoints. If ap and bp do not contain a minimum 
  ////// 	it will travel downhill until ap and cp contain a minimum. 
  // ap and bp are the endpoints of the initial guess (brackets) cp should be empty.
  // ap, bp, and cp will be altered so that the final result will be ap and cp contain 
  // 	a minimum where fxn(ap)>fxn(bp)<fxn(cp) and the point bp is between ap and cp.

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

  tools::F1DIMclass<fType> f1d(fxn, slope, p0);	// turning fxn into a 1d function given an initial point and slope
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
vector<double> tools::goldenMinSubR(fType fxn, vector<double>& ap, vector<double>& bp, vector<double>& cp, double difTol, double fracTol,  double& fmin, bool verbose) {
  ////// This function will find the minimum along the n-dimensional path between ap and cp using the golden ratio search
  //////	the exit condition can be triggered if the percent change in search signal fabs(x3-x0)/(fabs(x1)+fabs(x2)) 
  //////	is small or if the difference between the endpoints is small enough. Usually calling function will desire
  //////	difference condition as the final condition and keep lowering fracTol as the calling function converges on 
  //////	the minimum. This is because earlier search points do not need to find precise points until they are around min.
  //  REQUIREMENTS: ap and cp must be the endpoints of the search segment and contain bp
  //  fracTol: fractional tolerance, if fractional change in a search step is smaller then exit routine
  //  difTol: difference tolerance, if the difference between search endpoints is smaller then exit routine


  if (difTol <= 1e-7) {
    cerr<<"WARNING: due to double precision cannot find window less than 1e-8 precision, 1e-7 or smaller will take a long time (See NR 10.1)!!!"<<endl;
  }

  if (verbose) cout<<"STARTING GOLDENMINSUBR"<<endl<<endl;

  int Nvar = ap.size();
  int k;
  int iter = 0;
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
  for (k=0; k<Nvar; k++) { 
    if (fabs(cp[k]-ap[k])>0) {xb = (bp[k]-ap[k])/(cp[k]-ap[k]); break;} 
    else cerr<<"WARNING: Cannot find a middle point between bracket points in goldenMinSubR!!!"<<endl;
  }
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

  while (fabs(x3-x0)/(fabs(x1)+fabs(x2)) > fracTol && fabs(x3-x0) > difTol) {
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
    if (verbose) {
      cout<<"iter f1 f2: "<<iter<<"  "<<f1<<"  "<<f2<<endl;
      cout<<"Parameters: "<<x0<<"  "<<x1<<"  "<<x2<<"  "<<x3<<endl;
      cout<<"p0  "<<"p1  "<<"p2  "<<"p3  "<<endl;
      for (k=0; k<Nvar; k++) cout<<p0[k]+x0*slope[k]<<"  "<<p0[k]+x1*slope[k]<<"  "<<p0[k]+x2*slope[k]<<"  "<<p0[k]+x3*slope[k]<<endl;
      cout<<"Checking golden condition: "<<fabs(x0-x3)<<" < "<<difTol<<endl;
      cout<<"Checking golden condition: "<<fabs(x0-x3)/(fabs(x1)+fabs(x2))<<" < "<<fracTol<<endl<<endl;
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
*/





double tools::legendreIntegral(int order, double init, double fin) {
// -1 <= init,fin <= 1

  if (init < -1 || fin > 1) {
    cerr<<"ERROR: legendreIntegral requires -1 <= init,fin <= 1 !!!"<<endl;
    exit(0);
  }

  switch(order) {
    case 0:
	return (fin - init);
	break;
    case 1:
	return (pow(fin,2)-pow(init,2))/2;
	break;
    case 2:
	return ((pow(fin,3)-pow(init,3)) - (fin-init))/2;
	break;
    case 3:
	return ((pow(fin,4)-pow(init,4))*5/4 - (pow(fin,2)-pow(init,2))*3/2)/2;
	break;
    case 4:
	return ((pow(fin,5)-pow(init,5))*7 - (pow(fin,3)-pow(init,3))*10 + (fin-init)*3)/8;
	break;
    case 5:
	return ((pow(fin,6)-pow(init,6))*63/6 - (pow(fin,4)-pow(init,4))*70/4 + (pow(fin,2) - pow(init,2))*15/2)/8;
	break;
    case 6:
	return ((pow(fin,7)-pow(init,7))*231/7 - (pow(fin,5)-pow(init,5))*63 + (pow(fin,3)-pow(init,3))*35 - (fin-init)*5)/16;
	break;
    case 7:
	return ((pow(fin,8)-pow(init,8))*429/8 - (pow(fin,6)-pow(init,6))*693/6 + (pow(fin,4)-pow(init,4))*315/4 - (pow(fin,2)-pow(init,2))*35/2)/16;
	break;
    case 8:
	return ((pow(fin,9)-pow(init,9))*6435/9 - (pow(fin,7)-pow(init,7))*12012/7 + (pow(fin,5)-pow(init,5))*1386 - (pow(fin,3)-pow(init,3))*420 + (fin-init)*35)/128;
	break;
    case 9:
	return ((pow(fin,10)-pow(init,10))*1215.5 - (pow(fin,8)-pow(init,8))*25740/8 + (pow(fin,6)-pow(init,6))*18018/6 - (pow(fin,4)-pow(init,4))*4620/4 + (pow(fin,2)-pow(init,2))*315/2)/128;
	break;
    case 10:
	return ((pow(fin,11)-pow(init,11))*46189/11 - (pow(fin,9)-pow(init,9))*109395/9 + (pow(fin,7)-pow(init,7))*90090/7 - (pow(fin,5)-pow(init,5))*30030/5 + (pow(fin,3)-pow(init,3))*3465/3 - (fin-init)*63)/256;
	break;
    case 11:
	return ((pow(fin,12)-pow(init,12))*88179/12 - (pow(fin,10)-pow(init,10))*230945/10 + (pow(fin,8)-pow(init,8))*218790/8 - (pow(fin,6)-pow(init,6))*90090/6 + (pow(fin,4)-pow(init,4))*15015/4 - (pow(fin,2)-pow(init,2))*693/2)/256;
    case 12:
        return ((pow(fin,13)-pow(init,13))*676039/13 - (pow(fin,11)-pow(init,11))*1939938/11 + (pow(fin,9)-pow(init,9))*2078505/9 - (pow(fin,7)-pow(init,7))*1021020/7 + (pow(fin,5)-pow(init,5))*225225/5 - (pow(fin,3)-pow(init,3))*18018/3 + (fin-init)*231)/1024;
    case 13:
	return ((pow(fin,14)-pow(init,14))*1300075/14 - (pow(fin,12)-pow(init,12))*4056234/12 + (pow(fin,10)-pow(init,10))*4849845/10 - (pow(fin,8)-pow(init,8))*2771340/8 + (pow(fin,6)-pow(init,6))*765765/6 - (pow(fin,4)-pow(init,4))*90090/4 + (pow(fin,2)-pow(init,2))*3003/2)/1024;
    case 14:
	return (-429*(fin-init) + 45045*(pow(fin,3)-pow(init,3))/3 - 765765*(pow(fin,5)-pow(init,5))/5 + 4849845*(pow(fin,7)-pow(init,7))/7 - 14549535*(pow(fin,9)-pow(init,9))/9 + 22309287*(pow(fin,11)-pow(init,11))/11 - 16900975*(pow(fin,13)-pow(init,13))/13 + 5014575*(pow(fin,15)-pow(init,15))/15)/2048;
    case 15:
	return (9694845*(pow(fin,16)-pow(init,16))/16 - 35102025*(pow(fin,14)-pow(init,14))/14 + 50702925*(pow(fin,12)-pow(init,12))/12 - 37182145*(pow(fin,10)-pow(init,10))/10 + 14549535*(pow(fin,8)-pow(init,8))/8 - 2909907*(pow(fin,6)-pow(init,6))/6 + 255255*(pow(fin,4)-pow(init,4))/4 - 6435*(pow(fin,2)-pow(init,2))/2)/2048;
    case 16:
	return (300540195*(pow(fin,17)-pow(init,17))/17 - 1163381400*(pow(fin,15)-pow(init,15))/15 + 1825305300*(pow(fin,13)-pow(init,13))/13 - 1487285800*(pow(fin,11)-pow(init,11))/11 + 669278610*(pow(fin,9)-pow(init,9))/9 - 162954792*(pow(fin,7)-pow(init,7))/7 + 19399380*(pow(fin,5)-pow(init,5))/5 - 875160*(pow(fin,3)-pow(init,3))/3 + 6435*(fin-init))/32768;
    case 17:
	return (583401555*(pow(fin,18)-pow(init,18))/18 - 2404321560*(pow(fin,16)-pow(init,16))/16 + 4071834900*(pow(fin,14)-pow(init,14))/14 - 3650610600*(pow(fin,12)-pow(init,12))/12 + 1859107250*(pow(fin,10)-pow(init,10))/10 - 535422888*(pow(fin,8)-pow(init,8))/8 + 81477396*(pow(fin,6)-pow(init,6))/6 - 5542680*(pow(fin,4)-pow(init,4))/4 + 109395*(pow(fin,2)-pow(init,2))/2)/32768;
    case 18:
	return (2268783825*(pow(fin,19)-pow(init,19))/19 - 9917826435*(pow(fin,17)-pow(init,17))/17 + 18032411700*(pow(fin,15)-pow(init,15))/15 - 17644617900*(pow(fin,13)-pow(init,13))/13 + 10039179150*(pow(fin,11)-pow(init,11))/11 - 3346393050*(pow(fin,9)-pow(init,9))/9 + 624660036*(pow(fin,7)-pow(init,7))/7 - 58198140*(pow(fin,5)-pow(init,5))/5 + 2078505*(pow(fin,3)-pow(init,3))/3 - 12155*(fin-init))/65536;
    case 19:
	return (4418157975*(pow(fin,20)-pow(init,20))/20 - 20419054425*(pow(fin,18)-pow(init,18))/18 + 39671305740*(pow(fin,16)-pow(init,16))/16 - 42075627300*(pow(fin,14)-pow(init,14))/14 + 26466926850*(pow(fin,12)-pow(init,12))/12 - 10039179150*(pow(fin,10)-pow(init,10))/10 + 2230928700*(pow(fin,8)-pow(init,8))/8 - 267711444*(pow(fin,6)-pow(init,6))/6 + 14549535*(pow(fin,4)-pow(init,4))/4 - 230945*(pow(fin,2)-pow(init,2))/2)/65536;
    case 20:
 	return (34461632205*(pow(fin,21)-pow(init,21))/21 - 167890003050*(pow(fin,19)-pow(init,19))/19 + 347123925225*(pow(fin,17)-pow(init,17))/17 - 396713057400*(pow(fin,15)-pow(init,15))/15 + 273491577450*(pow(fin,13)-pow(init,13))/13 - 116454478140*(pow(fin,11)-pow(init,11))/11 + 30117537450*(pow(fin,9)-pow(init,9))/9 - 4461857400*(pow(fin,7)-pow(init,7))/7 + 334639305*(pow(fin,5)-pow(init,5))/5 - 9699690*(pow(fin,3)-pow(init,3))/3 + 46189*(fin-init))/262144;
    case 21:
	return (67282234305*(pow(fin,22)-pow(init,22))/22 - 344616322050*(pow(fin,20)-pow(init,20))/20 + 755505013725*(pow(fin,18)-pow(init,18))/18 - 925663800600*(pow(fin,16)-pow(init,16))/16 + 694247850450*(pow(fin,14)-pow(init,14))/14 - 328189892940*(pow(fin,12)-pow(init,12))/12 + 97045398450*(pow(fin,10)-pow(init,10))/10 - 17210021400*(pow(fin,8)-pow(init,8))/8 + 1673196525*(pow(fin,6)-pow(init,6))/6 - 74364290*(pow(fin,4)-pow(init,4))/4 + 969969*(pow(fin,2)-pow(init,2))/2)/262144;
    case 22:
 	return (263012370465*(pow(fin,23)-pow(init,23))/23 - 1412926920405*(pow(fin,21)-pow(init,21))/21 + 3273855059475*(pow(fin,19)-pow(init,19))/19 - 4281195077775*(pow(fin,17)-pow(init,17))/17 + 3471239252250*(pow(fin,15)-pow(init,15))/15 - 1805044411170*(pow(fin,13)-pow(init,13))/13 + 601681470390*(pow(fin,11)-pow(init,11))/11 - 124772655150*(pow(fin,9)-pow(init,9))/9 + 15058768725*(pow(fin,7)-pow(init,7))/7 - 929553625*(pow(fin,5)-pow(init,5))/5 + 22309287*(pow(fin,3)-pow(init,3))/3 - 88179*(fin-init))/524288;
    case 23:
	return (514589420475*(pow(fin,24)-pow(init,24))/24 - 2893136075115*(pow(fin,22)-pow(init,22))/22 + 7064634602025*(pow(fin,20)-pow(init,20))/20 - 9821565178425*(pow(fin,18)-pow(init,18))/18 + 8562390155550*(pow(fin,16)-pow(init,16))/16 - 4859734953150*(pow(fin,14)-pow(init,14))/14 + 1805044411170*(pow(fin,12)-pow(init,12))/12 - 429772478850*(pow(fin,10)-pow(init,10))/10 + 62386327575*(pow(fin,8)-pow(init,8))/8 - 5019589575*(pow(fin,6)-pow(init,6))/6 + 185910725*(pow(fin,4)-pow(init,4))/4 - 2028117*(pow(fin,2)-pow(init,2))/2)/524288;
    case 24:
	return (8061900920775*(pow(fin,25)-pow(init,25))/25 - 47342226683700*(pow(fin,23)-pow(init,23))/23 + 121511715154830*(pow(fin,21)-pow(init,21))/21 - 178970743251300*(pow(fin,19)-pow(init,19))/19 + 166966608033225*(pow(fin,17)-pow(init,17))/17 - 102748681866600*(pow(fin,15)-pow(init,15))/15 + 42117702927300*(pow(fin,13)-pow(init,13))/13 - 11345993441640*(pow(fin,11)-pow(init,11))/11 + 1933976154825*(pow(fin,9)-pow(init,9))/9 - 194090796900*(pow(fin,7)-pow(init,7))/7 + 10039179150*(pow(fin,5)-pow(init,5))/5 - 202811700*(pow(fin,3)-pow(init,3))/3 + 676039*(fin-init))/4194304;
    case 25:
	return (15801325804719*(pow(fin,26)-pow(init,26))/26 - 96742811049300*(pow(fin,24)-pow(init,24))/24 + 260382246760350*(pow(fin,22)-pow(init,22))/22 - 405039050516100*(pow(fin,20)-pow(init,20))/20 + 402684172315425*(pow(fin,18)-pow(init,18))/18 - 267146572853160*(pow(fin,16)-pow(init,16))/16 + 119873462177700*(pow(fin,14)-pow(init,14))/14 - 36100888223400*(pow(fin,12)-pow(init,12))/12 + 7091245901025*(pow(fin,10)-pow(init,10))/10 - 859544957700*(pow(fin,8)-pow(init,8))/8 + 58227239070*(pow(fin,6)-pow(init,6))/6 - 1825305300*(pow(fin,4)-pow(init,4))/4 + 16900975*(pow(fin,2)-pow(init,2))/2)/4194304;
    case 26:
	return (61989816618513*(pow(fin,27)-pow(init,27))/27 - 395033145117975*(pow(fin,25)-pow(init,25))/25 + 1112542327066950*(pow(fin,23)-pow(init,23))/23 - 1822675727322450*(pow(fin,21)-pow(init,21))/21 + 1923935489951475*(pow(fin,19)-pow(init,19))/19 - 1369126185872445*(pow(fin,17)-pow(init,17))/17 + 667866432132900*(pow(fin,15)-pow(init,15))/15 - 222622144044300*(pow(fin,13)-pow(init,13))/13 + 49638721307175*(pow(fin,11)-pow(init,11))/11 - 7091245901025*(pow(fin,9)-pow(init,9))/9 + 601681470390*(pow(fin,7)-pow(init,7))/7 - 26466926850*(pow(fin,5)-pow(init,5))/5 + 456326325*(pow(fin,3)-pow(init,3))/3 - 1300075*(fin-init))/8388608;
    case 27:
	return (121683714103007*(pow(fin,28)-pow(init,28))/28 - 805867616040669*(pow(fin,26)-pow(init,26))/26 + 2370198870707850*(pow(fin,24)-pow(init,24))/24 - 4079321865912150*(pow(fin,22)-pow(init,22))/22 + 4556689318306125*(pow(fin,20)-pow(init,20))/20 - 3463083881912655*(pow(fin,18)-pow(init,18))/18 + 1825501581163260*(pow(fin,16)-pow(init,16))/16 - 667866432132900*(pow(fin,14)-pow(init,14))/14 + 166966608033225*(pow(fin,12)-pow(init,12))/12 - 27577067392875*(pow(fin,10)-pow(init,10))/10 + 2836498360410*(pow(fin,8)-pow(init,8))/8 - 164094946470*(pow(fin,6)-pow(init,6))/6 + 4411154475*(pow(fin,4)-pow(init,4))/4 - 35102025*(pow(fin,2)-pow(init,2))/2)/8388608;
    case 28:
	return (956086325095055*(pow(fin,29)-pow(init,29))/29 - 6570920561562378*(pow(fin,27)-pow(init,27))/27 + 20146690401016725*(pow(fin,25)-pow(init,25))/25 - 36343049350853700*(pow(fin,23)-pow(init,23))/23 + 42832879592077575*(pow(fin,21)-pow(init,21))/21 - 34630838819126550*(pow(fin,19)-pow(init,19))/19 + 19624141997505045*(pow(fin,17)-pow(init,17))/17 - 7823578204985400*(pow(fin,15)-pow(init,15))/15 + 2170565904431925*(pow(fin,13)-pow(init,13))/13 - 408140597414550*(pow(fin,11)-pow(init,11))/11 + 49638721307175*(pow(fin,9)-pow(init,9))/9 - 3610088822340*(pow(fin,7)-pow(init,7))/7 + 136745788725*(pow(fin,5)-pow(init,5))/5 - 2035917450*(pow(fin,3)-pow(init,3))/3 + 5014575*(fin-init))/33554432;
    case 29:
	return (1879204156221315*(pow(fin,30)-pow(init,30))/30 - 13385208551330770*(pow(fin,28)-pow(init,28))/28 + 42710983650155457*(pow(fin,26)-pow(init,26))/26 - 80586761604066900*(pow(fin,24)-pow(init,24))/24 + 99943385714847675*(pow(fin,22)-pow(init,22))/22 - 85665759184155150*(pow(fin,20)-pow(init,20))/20 + 51946258228689825*(pow(fin,18)-pow(init,18))/18 - 22427590854291480*(pow(fin,16)-pow(init,16))/16 + 6845630929362225*(pow(fin,14)-pow(init,14))/14 - 1447043936287950*(pow(fin,12)-pow(init,12))/12 + 204070298707275*(pow(fin,10)-pow(init,10))/10 - 18050444111700*(pow(fin,8)-pow(init,8))/8 + 902522205585*(pow(fin,6)-pow(init,6))/6 - 21037813650*(pow(fin,4)-pow(init,4))/4 + 145422675*(pow(fin,2)-pow(init,2))/2)/33554432;
    case 30:
	return (7391536347803839*(pow(fin,31)-pow(init,31))/31 -54496920530418135*(pow(fin,29)-pow(init,29))/29 +180700315442965395*(pow(fin,27)-pow(init,27))/27 - 355924863751295475*(pow(fin,25)-pow(init,25))/25 + 463373879223384675*(pow(fin,23)-pow(init,23))/23 - 419762220002360235*(pow(fin,21)-pow(init,21))/21 +271274904083157975*(pow(fin,19)-pow(init,19))/19 - 126155198555389575*(pow(fin,17)-pow(init,17))/17 + 42051732851796525*(pow(fin,15)-pow(init,15))/15 - 9888133564634325*(pow(fin,13)-pow(init,13))/13 + 1591748329916745*(pow(fin,11)-pow(init,11))/11 - 166966608033225*(pow(fin,9)-pow(init,9))/9 + 10529425731825*(pow(fin,7)-pow(init,7))/7 - 347123925225*(pow(fin,5)-pow(init,5))/5 + 4508102925*(pow(fin,3)-pow(init,3))/3 - 9694845*(fin-init))/67108864;
    default:
	cerr<<"ERROR: Cannot compute the integral for Legendre of order "<<order<<", you must add this!!!"<<endl;
	exit(0);
  }
}



std::vector< std::vector<double> > tools::fft1dRtoC(const std::vector<double> &input, 
      double outRat, int Npadding, double padDecayRat, const fftw_plan &fftPlan, 
      double* in, fftw_complex* out, PLOTclass* pltVerbose, std::string pltName) {

  int ir, centShift;
  int centI = (int)(input.size()/2);
  int fftSize = input.size() + Npadding;
  int outSize = (int)(fftSize*outRat/2);
  const int indData = (int)(input.size()/2-1);
  std::vector< std::vector<double> > results(4);
  for (auto& v : results) {
    v.resize(outSize, 0);
  }
  results[3].resize(fftSize, 0);
  if (input.size()%2 == 1) {
    centShift = 1;
    in[0] = input[centI];
  }
  else {
    centShift = 0;
  }
  for (ir=0; ir<indData; ir++) {
    in[centShift+ir] = input[centI+1+ir];
    in[fftSize-1-ir] = input[centI-1-ir];
  }

  const int NpadDecay = Npadding*padDecayRat/2;
  for (ir=0; ir<Npadding/2; ir++) {
    if (ir<NpadDecay) {
      in[indData + ir] = in[indData]*pow(sin((PI/2)*(ir/NpadDecay)), 2);
      in[fftSize-indData-ir] = in[fftSize-indData]*pow(sin((PI/2)*(ir/NpadDecay)), 2);
    } 
    else {
      in[indData+ir] = 0;
      in[fftSize-indData-ir] = 0;
    }
  }
  for (ir=0; ir<fftSize; ir++) {
    results[3][ir] = in[ir];
  }

  if (pltVerbose) {
    std::vector<double> plotMe(fftSize,0);
    for (int ir=0; ir<fftSize; ir++) {
      plotMe[ir] = in[ir];
    }
    pltVerbose->print1d(plotMe, "./plots/" + pltName);
  }

  fftw_execute(fftPlan);

  for (ir=0; ir<outSize; ir++) {
      results[0][ir] = out[ir][0]/(fftSize);
      results[1][ir] = out[ir][1]/(fftSize);
      results[2][ir] = sqrt(pow(results[0][ir], 2) + pow(results[1][ir], 2));
  }

  return results;
}




std::vector< std::vector<double> > tools::fft1dRtoC(const std::vector<double> &input, 
      const fftw_plan &fftPlan, double* in, fftw_complex* out, 
      bool forwardFFT, bool centeredOrigin, 
      PLOTclass* pltVerbose, std::string pltName) {

  /////  Variables  /////
  int inpSize = input.size();
  int outSize = inpSize/2 + 1;
  int centI = inpSize/2;

  /////  Format input for fftw  /////
  if (centeredOrigin) {
    for (int i=0; i<inpSize; i++) {
      in[i] = input[(centI+i)%inpSize];
    }
  }
  else { 
    for (int i=0; i<inpSize; i++) {
      in[i] = input[i];
    }
  }

  //  Plot input
  if (pltVerbose) {
    std::vector<double> plotMe(inpSize,0);
    for (int i=0; i<inpSize; i++) {
      plotMe[i] = in[i];
    }
    pltVerbose->print1d(plotMe, "./plots/" + pltName);
  }

  /////  FFT  /////
  fftw_execute(fftPlan);

  /////  Return results  /////
  double scale = std::sqrt(inpSize);
  float FFTdirection = 1;
  if (!forwardFFT) {
    FFTdirection = -1;
  }
  std::vector< std::vector<double> > results(4);
  for (auto& v : results) {
    v.resize(outSize, 0);
  }
  for (int ir=0; ir<outSize; ir++) {
      results[0][ir] = out[ir][0]/scale;
      results[1][ir] = out[ir][1]*FFTdirection/scale;
      results[2][ir] = sqrt(pow(results[0][ir], 2) + pow(results[1][ir], 2));
  }

  return results;
}

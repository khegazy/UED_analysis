#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#include <dirent.h>

//From OpenCV
#include <opencv2/core/core.hpp>
#include <opencv2/core/types_c.h>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/objdetect/objdetect.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/contrib/contrib.hpp>

//HomeGrown
#include "constants.h"
#include "tools.h"
#include "saving.h"
#include "plotClass.h"

using namespace std;


namespace imgProc {


  struct imgInfoStruct {
 
    string date;
    string path;
    string fileName;
    string run;
    string runType;
    int scan;
    int imgNum;
    int throttle;
    long int time;
    int32_t stagePos;
 
    imgInfoStruct () {
      scan = -999;
    }
  };


  template<typename type>
  std::vector< std::vector<double> > getImgVector(cv::Mat imgMat,
      int roiSize, int centerR, int centerC, bool subMedian = false,
      std::vector< std::vector<double> > *background = NULL);
      //int rmRad=-1, int rmCenterR=-1, int rmCenterC=-1,
      //std::vector< std::vector<int> >* nanMap=NULL);
      //bool eadoutSubtraction=false, 
      //int roR=970, int roC=30, int window=30);

  void saveAscii(cv::Mat mat, string fileName);
  void saveAscii(TH2F *hist, string fileName);

  void threshold(vector< vector<double> > &img, int hotpixel);
  void threshold(cv::Mat &img, int hotpixel);
  template<typename T>
  void thresholdT(cv::Mat &img, int hotpixel);
  std::vector< std::vector<int> > removeXrayHits(
    std::vector<std::vector<double> >* img,
    std::vector<std::vector<int> > nanMap,
    double highCut, double lowCut, 
    double stdCut, int window,
    TH1F** plots=NULL);
  void averagePixel(cv::Mat &img, int row, int col);
  template<typename T>
  void averagePixelT(cv::Mat &img, int row, int col);

  void medianPixel(cv::Mat &mat, int row, int col, int hotpixel);
  template<typename T>
  void medianPixelT(cv::Mat &mat, int row, int col, int hotpixel);

  std::vector<double> gaussianSmooth1d(
      std::vector<double> inp, 
      double sigma, 
      int totalRange);

  double cubicInterpolate(double p[4], double x);
  double bicubicInterpolate(double p[4][4], double row, double col);
  template<typename T>
  double interpolateT(cv::Mat mat, double row, double col);
  double interpolate(cv::Mat mat, double row, double col);

  void histMeanStd(vector< vector<double> >& hist, double& mean, double& std);

  double imageNorm(
      std::vector< std::vector<double> > &img, 
      std::vector< std::vector<int> > nanMap,
      float radMin=0.7, float radMax=0.72);


  ////////////////////////////
  /////  Image Cleaning  /////
  ////////////////////////////

  /////  Radial processing  /////
  class radProcTool {
    public:

      //radProcTool(string inpIndPath);
      radProcTool(int shellWidth, int maxRad);

      double getMean(
          vector<double> vals, int leftCut=-1, int rightCut=-1);
      double getMMean(
          vector<double> vals, vector<int> ord,
          double range, bool verbose=false);
      void   getSmearedDist(std::map<int, double> &smearedDist, 
                  vector<double> &vals, double stdev, bool verbose=false);
      double getSTDev(
          vector<double> vals, double mean,
          int leftCut=-1, int rightCut=-1);
      double getLeftSTDev(
          vector<double> vals, vector<int> ord, double mean, int leftCut=-1);

      int    getNrightEntries(vector<double> vals, vector<int> ord, double mean);
      int    getNleftEntries(vector<double> vals, vector<int> ord, double mean);
      int    getLeftOutliers(vector<double> &vals, vector<int> &orderedInds, double mean, double stdev, double Nstdev);
      int    getRightOutliers(vector<double> &vals, vector<int> &orderedInds, double mean, double stdev, double Nstdev);
      std::vector< std::vector<double> > removeOutliers(
              vector< vector<double> > &image,
              vector< vector<int> > &nanMap,
              int centerR, int centerC, int buffer,
              int maxRad, int NshellOutlierLoops,
              int shellWidth, int Npoly,
              double stdOutlierCutLeft, double stdOutlierCutRight,
              double stdCutLeft, double stdCutRight,
              int stg, double outlierMapSTDcut,   
              bool getOutlierImage, bool verbose,
              PLOTclass* pltVerbose, TH1F** radPixHistos=NULL);

      std::vector< std::vector<double> > removeOutliers_stdRatio(
              vector< vector<double> > &image,
              vector< vector<int> > &nanMap,
              int centerR, int centerC, int buffer,
              int maxRad, int NshellOutlierLoops,
              int shellWidth, int Npoly,
              double stdChangeRatioLeft, double stdChangeRatioRight,
              double stdAccRatioLeft, double stdAccRatioRight,
              double stdCutLeft, double stdCutRight,
              double fracShellSTDcutLeft, double fracShellSTDcutRight,
              int stg, double outlierMapSTDcut,   
              bool getOutlierImage, bool verbose,
              PLOTclass* pltVerbose, TH1F** radPixHistos=NULL);

      std::vector< std::vector<double> > removeOutliers_stdInclude(
              vector< vector<double> > &image,
              vector< vector<int> > &nanMap,
              int centerR, int centerC, int buffer,
              int maxRad, int NshellOutlierLoops,
              int shellWidth, int Npoly,
              double stdIncludeLeft, double distSTDratioLeft,
              double stdCutLeft, double fracShellSTDcutLeft, 
              double stdIncludeRight, double distSTDratioRight,
              double stdCutRight, double fracShellSTDcutRight,
              double stdChangeRatio, int stg, double outlierMapSTDcut,   
              bool getOutlierImage, bool verbose,
              PLOTclass* pltVerbose, TH1F** radPixHistos=NULL);
      std::vector< std::vector<double> > removeOutliersSimple(
              vector< vector<double> > &image,
              vector< vector<int> > &nanMap,
              float centerR_f, float centerC_f, int buffer,
              int maxRad, int shellWidth, 
              bool removeLowPoly, int Npoly,
              double stdCut, int stg, double outlierMapSTDcut,
              bool getOutlierImage, bool verbose, PLOTclass* pltVerbose);
      vector<double> getPolarLineOut(vector< vector<double> >* image, 
              int centerR, int centerC, int rad, 
              int shellWidth, int NangleBins, bool verbose=false);

      double radialSliceVar(vector< vector<double> >* image,
              int centerR, int centerC, 
              int rad, int shellWidth, bool verbose=false);

    private:
      string indPath;
      int makeMaxRad;

      std::map<int, std::map<int, std::vector< std::vector<int> > > > allIndices;
      void importIndices(int shellWidth, int rad);
      void makeIndices(int shellWidth, int maxRad);
      bool checkForIndices(int shellWidth, int rad);
  };

  std::vector< std::vector<double> > asymmetrize(
          std::vector< std::vector<double> > inpImg,
          int centerR, int centerC,
          int height, int width,
          std::vector< std::vector<double> > oddImgImgn,
          fftw_plan &fftFref, fftw_complex* fftIn,
          fftw_plan &fftBref, fftw_complex* fftOut,
          PLOTclass* plt=NULL);

  std::vector< std::pair<uint, uint> > findClusters(
          std::vector< std::vector<double> > &inpImg,
          int centerR, int centerC, double minCenterRad,
          double coreValThresh, int coreRad,
          int minClusterSize, int minPixelSize, 
          double minDensity, int clusterRad,
          double borderValThresh, int borderDistLimit,
          int borderRad, int padRad, 
          PLOTclass* pltVerbose = NULL);

  std::vector< std::pair<uint, uint> > findClusters(
          std::vector< std::vector<double> > &inpImg, 
          int centerR, int centerC, double minCenterRad,
          double coreValThresh, int coreRad,
          int minClusterSize, int minPixelSize, 
          double minDensity, int clusterRad,
          double borderValThresh, int borderDistLimit,
          int borderRad, int padRad,
          std::vector< std::vector<double> > &clusters_coreBorder, 
          PLOTclass* pltVerbose = NULL); 

  std::vector< std::pair<uint, uint> > findClusters(
          std::vector< std::vector<double> > &inpImg, 
          int centerR, int centerC, double minCenterRad,
          double coreValThresh, int coreRad,
          int minClusterSize, int minPixelSize, 
          double minDensity, double shapeVarRatio, double shapeEdgeRatio, int clusterRad,
          double borderValThresh, int borderDistLimit,
          int borderRad, int padRad,
          float rMinScale, float rMaxScale, float cMinScale, float cMaxScale,
          std::vector< std::vector<double> > &clusters_coreBorder, 
          PLOTclass* pltVerbose = NULL); 


  double removeReadOutNoise(
          std::vector< std::vector<double> > &img,
          std::vector< std::vector<int> > &nanMap,
          int roiR=950, int roiC=30, int window=50, int buffer=20);
  double removeAvgReadOutNoise(
          std::vector< std::vector<double> > &img,
          std::vector< std::vector<int> > &nanMap,
          int centerR, int centerC, double minRad, double maxRad, int buffer=20);
  double removeMedianReadOutNoise(std::vector< std::vector<double> > &img,
          int centerR, int centerC, double minRad, double maxRad, int buffer=20, 
          std::vector< std::vector<int> >* nanMap=NULL);


  /////// Center Finding ///////
  class centerfnctr {
      public:
          int fxnType;
          int minRadBin;
          int centShellW;
          double meanVal, valSTD;
          radProcTool* radProc;
          std::vector<int>* indsR;
          std::vector<int>* indsC;
          std::vector<double>* vals;
          std::vector< std::vector<double> >* img;
          double operator() (std::vector<double> vect);
  };

  std::vector<int> centerSearchCOM(
      std::vector<imgInfoStruct> &imgINFO,
      std::vector< std::vector<int> > nanMap,
      double hotPixel, double sigmaX,
      int blockCentR, int blockCentC, 
      float minRad, float maxRad,
      int meanInd, double stdScale,  
      bool verbose, PLOTclass* pltVerbose=NULL);
  std::vector<int> centerSearchCOM(std::vector< std::vector<double> > imgCent,
      int blockCentR, int blockCentC, 
      float minRad, float maxRad,
      int meanInd, double stdScale,  
      bool verbose, PLOTclass* pltVerbose=NULL);

  // Teiring image and finding contours
  cv::Mat buildTeirGrayImg(cv::Mat &inp, long int modVal, long int imScale=1);
  template<typename T>
  cv::Mat buildTeirGrayImgT(cv::Mat &inp, long int modVal, long int imScale=1);

  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, int cannyLow=200, int cannyHigh=400);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int inRad, int outRad);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, int Ncontours);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow=200, int cannyHigh=400);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int inRad, int outRad);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad);
  vector<vector<cv::Point> > findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int innerRad, int outterRad, cv::Mat &cannyEdges, vector<cv::Vec4i> &hierarchy, int Ncontours);

  cv::Point findContourCenters(vector< vector<cv::Point> > contours);
  cv::Point findContourCenters(vector< vector<cv::Point> > contours, vector<cv::Point> &centers);

  cv::Point findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int inRad, int outRad);
  cv::Point findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad);
  cv::Point findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int inRad, int outRad);
  cv::Point findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad);
  cv::Point findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Point> &centers, int Ncontours);
  cv::Point findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Point> &centers);
  cv::Point findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Vec4i> &hierarchy, vector<cv::Point> &centers);

  // Summing opposite sections of 2d shell and comparing
  //double centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW, double angle=0);
  double centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW, bool nanBins);
  double centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW);
  double centerOddLeg(vector< vector<double> >* img, double centerR, double centerC, int Nshells, double minR, double shellW, bool nanBins);

/*
  class centerfnctr {
  
      public:
          int fxnType;
          int minRadBin;
          int centShellW;
          imgProc::radProcTool* radProc = NULL;
          std::vector< std::vector<double> >* img;
  
          centerfnctr();
          
          double operator() (std::vector<double> vect);
   };
*/



  ////// Rebin to polar binning //////
  vector< vector<double> > polarBinning(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins);
  vector< vector<double> > polarBinning(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan);
  template<typename T>
  vector< vector<double> > polarBinningT(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan);
  vector< vector<double> > polarBinning(vector< vector<double> > mat, double centerR, double centerC, int Nslices, double NmatBins, double NhistBins);
  vector< vector<double> > polarBinning(vector< vector<double> > mat, double centerR, double centerC, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan);


  void makeMaskRad(double ang, double radI, double radF, int rows, int cols, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask);
  void makeMaskRC(double ang, double radI, double radF, int rads, int angs, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask);
  void makeMask(bool rad, double ang, double radI, double radF, vector< vector<double> > &vect, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask);



  ////// Legendre Coefficients //////
  std::vector<double> legendreFit(std::vector< std::vector<double> > img, int rebin, const int Nlg, const int NradBins, const int Nrows, const int Ncols);
  std::vector<double> legendreFit(std::vector< std::vector<double> > img, int rebin, const int Nlg, const int NradBins, const int Nrows, const int Ncols, const Eigen::MatrixXd g);
  double legendre1dCoeff(TH1F* hist, int order);
  double legendre1dCoeff(vector<double> &vect, int order);
  double legendre1dCoeff(vector<double> &vect, int order, double nan);
  double legendre1dCosCoeff(vector<double> &vect, int order);
  double legendre1dCosCoeff(vector<double> &vect, int order, double nan);
  TH2F* polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order);
  TH2F* polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order, double &max);
  TH2F* polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order, bool maxes, double &max);
  std::vector<double> fitImgLegendre(
          std::vector< std::vector<double> > &img,
          const int Nrows, const int rebinRow,
          const int Ncols, const int rebinCol,
          const int Nlg, const int NradBins,
          bool verbose);


  ////// Others //////
  vector< vector<double> > imgRotatePId2(vector< vector<double> > origImg);

  double Aparam(vector<double> vect);
  double Aparam(vector< vector<double> > vect, vector< pair<int,int> > pmask, vector< pair<int,int> > nmask);

  vector< vector<cv::Point> > centContours(cv::Mat &img, vector<cv::Vec4i> &hierarchy, int low=200, int high=400);
  cv::Point findCenter(vector< vector<cv::Point> > &contours);
  cv::Point findCenter(cv::Mat &img, vector< vector<cv::Point> > &contours, vector<float> &radii, vector<cv::Vec4i> &hierarchy, int low=200, int high=400);

} 







template<typename type>
std::vector< std::vector<double> > imgProc::getImgVector(cv::Mat imgMat,
    int roiSize, int centerR, int centerC, bool subMedian,
    std::vector< std::vector<double> > *background) {



  int icol, irow, index;
  int Nrows = (int)imgMat.rows;
  int Ncols = (int)imgMat.cols;
  type* pv = imgMat.ptr<type>(0);
  double readoutNoise = 0;
  double rad = 0;
  std::vector< std::vector<double> > img(roiSize);


  if ((centerR < roiSize/2) || (centerR + roiSize/2 > Nrows)) {
    cerr << "ERROR: Cannot center image with center (" 
      << centerR << ", " << centerC <<") and " << Nrows << " rows!!!\n";
    exit(0);
  }
  if ((centerC < roiSize/2) || (centerC + roiSize/2 > Ncols)) {
    cerr << "ERROR: Cannot center image with center (" 
      << centerR << ", " << centerC <<") and " << Ncols << " columns!!!\n";
    exit(0);
  }


  /*
  if (removeReadout) {
    std::vector<double> readoutWindow;
    for (int ir=roR; ir<roR+window; ir++) {
      if (ir >= Nrows) {
        cerr << "ERROR: readout noise window exceeds rows ("
          << ir << ">=" << Nrows <<"!!!\n" << endl;
        exit(0);
      }
      for (int ic=roC; ic<roC+window; ic++) {
        if (ic >= Ncols) {
          cerr << "ERROR: readout noise window exceeds cols ("
            << ic << ">=" << Ncols << "!!!\n" << endl;
          exit(0);
        }

        index = ir*Ncols + ic;
        readoutWindow.push_back(pv[index]);
      }
    }

    std::vector<int> orderedInds(readoutWindow.size(), 0);
    iota(orderedInds.begin(), orderedInds.end(), 0);
    sort(orderedInds.begin(), orderedInds.end(),
          [&readoutWindow](uint i1, uint i2)
          {return readoutWindow[i1] < readoutWindow[i2];});

    readoutNoise = readoutWindow[orderedInds[(uint)(readoutWindow.size()/2.)]];
  }
  */

  double median = 0;
  if (subMedian) {
    int sbInd = 0;
    std::vector<int> sbMedVec(100);
    for (int ir=0; ir<10; ir++) {
      for (int ic=0; ic<10; ic++) {
        index = ir*Nrows + ic;
        sbMedVec[sbInd] = pv[index];
        sbInd++;
      }
    }

    std::sort(sbMedVec.begin(), sbMedVec.end());

    if (sbMedVec.size() % 2 == 0) {
      median = (sbMedVec[(int)(sbMedVec.size()/2)]
                + sbMedVec[(int)(sbMedVec.size()/2)-1])/2.;
    }
    else {
      median = sbMedVec[(int)(sbMedVec.size()/2)];
    }
  }

  for (int ir=0; ir<roiSize; ir++) {
    img[ir].resize(roiSize);
    irow = ir - roiSize/2 + centerR;
    for (int ic=0; ic<roiSize; ic++) {
      icol = ic - roiSize/2 + centerC;
      index = irow*Ncols + icol;
      img[ir][ic] = ((double)(pv[index])) - median;
      if (background) {
        img[ir][ic] -= (*background)[irow][icol];
      }
    }
  }

  return img;
}

template<typename T>
double imgProc::interpolateT(cv::Mat mat, double row, double col) {

  int clow = int(col);
  int rlow = int(row);
  double p[4][4] = {0};

  if (((clow+1)%mat.cols >= mat.cols) || (clow%mat.cols < 0) || ((rlow+1)%mat.rows >= mat.rows) || (rlow%mat.rows < 0)) {
    cerr<<"ERROR: Cannot interpolate point between col("<<clow<<", "<<clow+1<<") and row("<<rlow<<", "<<rlow+1<<") which is out of range of "<<mat.rows<<"X"<<mat.cols<<" matrix!!!"<<endl;
    exit(0);
  }

  bool padLeft, padRight, padTop, padBottom;
  T* val;

  if (mat.isContinuous()) {

    val = mat.ptr<T>(0);
    padLeft = (0 < clow%mat.cols);
    padRight = (clow%mat.cols < (mat.cols-2));
    padBottom = (rlow%mat.rows < (mat.rows-2));
    padTop = ( 0 < rlow%mat.rows);

    p[1][1] = val[clow + mat.rows*rlow];      	p[2][1] = val[clow + mat.rows*(rlow+1)];
    p[1][2] = val[clow+1 + mat.rows*rlow]; 	p[2][2] = val[clow+1 + mat.rows*(rlow+1)];

//cout<<"Bools: "<<padTop<<"  "<<padBottom<<"  "<<padLeft<<"  "<<padRight<<endl;
    if (padLeft && padRight && padTop && padBottom) {
      p[0][0] = val[clow-1 + mat.rows*(rlow-1)];        p[1][0] = val[clow-1 + mat.rows*rlow];
      p[2][0] = val[clow-1 + mat.rows*(rlow+1)];        p[3][0] = val[clow-1 + mat.rows*(rlow+2)];
      p[0][1] = val[clow + mat.rows*(rlow-1)];          p[3][1] = val[clow + mat.rows*(rlow+2)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[3][2] = val[clow+1 + mat.rows*(rlow+2)];
      p[0][3] = val[clow+2 + mat.rows*(rlow-1)];        p[1][3] = val[clow+2 + mat.rows*rlow];
      p[2][3] = val[clow+2 + mat.rows*(rlow+1)];        p[3][3] = val[clow+2 + mat.rows*(rlow+2)];
    }
        else if (!padLeft && padRight && padTop && padBottom) {
      p[0][1] = val[clow + mat.rows*(rlow-1)];          p[3][1] = val[clow + mat.rows*(rlow+2)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[3][2] = val[clow+1 + mat.rows*(rlow+2)];
      p[0][3] = val[clow+2 + mat.rows*(rlow-1)];        p[1][3] = val[clow+2 + mat.rows*rlow];
      p[2][3] = val[clow+2 + mat.rows*(rlow+1)];        p[3][3] = val[clow+2 + mat.rows*(rlow+2)];
      p[0][0] = p[0][1] - (p[0][2]-p[0][1]);            p[1][0] = p[1][1] - (p[1][2]-p[1][1]);
      p[2][0] = p[2][1] - (p[2][2]-p[2][1]);            p[3][0] = p[3][1] - (p[3][2]-p[3][1]);
    }
    else if (padLeft && !padRight && padTop && padBottom) {
      p[0][1] = val[clow + mat.rows*(rlow-1)];          p[3][1] = val[clow + mat.rows*(rlow+2)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[3][2] = val[clow+1 + mat.rows*(rlow+2)];
      p[0][0] = val[clow-1 + mat.rows*(rlow-1)];        p[1][0] = val[clow-1 + mat.rows*rlow];
      p[2][0] = val[clow-1 + mat.rows*(rlow+1)];        p[3][0] = val[clow-1 + mat.rows*(rlow+2)];
      p[0][3] = p[0][2] + (p[0][2]-p[0][1]);            p[1][3] = p[1][2] + (p[1][2]-p[1][1]);
      p[2][3] = p[2][2] + (p[2][2]-p[2][1]);            p[3][3] = p[3][2] + (p[3][2]-p[3][1]);
    }
    else if (padLeft && padRight && !padTop && padBottom) {
      p[1][0] = val[clow-1 + mat.rows*rlow];            p[1][3] = val[clow+2 + mat.rows*rlow];
      p[2][0] = val[clow-1 + mat.rows*(rlow+1)];        p[2][3] = val[clow+2 + mat.rows*(rlow+1)];
      p[3][0] = val[clow-1 + mat.rows*(rlow+2)];        p[3][1] = val[clow + mat.rows*(rlow+2)];
      p[3][2] = val[clow+1 + mat.rows*(rlow+2)];        p[3][3] = val[clow+2 + mat.rows*(rlow+2)];
      p[0][0] = p[1][0] - (p[2][0]-p[1][0]);            p[0][1] = p[1][1] - (p[2][1]-p[1][1]);
      p[0][2] = p[1][2] - (p[2][2]-p[1][2]);            p[0][3] = p[1][3] - (p[2][3]-p[1][3]);
    }
    else if (padLeft && padRight && padTop && !padBottom) {
      p[1][0] = val[clow-1 + mat.rows*rlow];            p[1][3] = val[clow+2 + mat.rows*rlow];
      p[2][0] = val[clow-1 + mat.rows*(rlow+1)];        p[2][3] = val[clow+2 + mat.rows*(rlow+1)];
      p[0][0] = val[clow-1 + mat.rows*(rlow-1)];        p[0][1] = val[clow + mat.rows*(rlow-1)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[0][3] = val[clow+2 + mat.rows*(rlow-1)];
      p[3][0] = p[2][0] + (p[2][0]-p[1][0]);            p[3][1] = p[2][1] + (p[2][1]-p[1][1]);
      p[3][2] = p[2][2] + (p[2][2]-p[1][2]);            p[3][3] = p[2][3] + (p[2][3]-p[1][3]);
    }
    else if (!padLeft && padRight && padTop && !padBottom) {
      p[1][3] = val[clow+2 + mat.rows*rlow];            p[2][3] = val[clow+2 + mat.rows*(rlow+1)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[0][3] = val[clow+2 + mat.rows*(rlow-1)];
      p[0][1] = val[clow + mat.rows*(rlow-1)];          p[0][0] = p[0][1] - (p[0][2]-p[0][1]);
      p[1][0] = p[1][1] - (p[1][2]-p[1][1]);            p[2][0] = p[2][1] - (p[2][2]-p[2][1]);
      p[3][0] = p[2][0] + (p[2][0]-p[1][0]);            p[3][1] = p[2][1] + (p[2][1]-p[1][1]);
      p[3][2] = p[2][2] + (p[2][2]-p[1][2]);            p[3][3] = p[2][3] + (p[2][3]-p[1][3]);
    }
    else if (!padLeft && padRight && !padTop && padBottom) {
      p[1][3] = val[clow+2 + mat.rows*rlow];            p[2][3] = val[clow+2 + mat.rows*(rlow+1)];
      p[3][1] = val[clow + mat.rows*(rlow+2)];          p[3][2] = val[clow+1 + mat.rows*(rlow+2)];
      p[3][3] = val[clow+2 + mat.rows*(rlow+2)];        p[1][0] = p[1][1] - (p[1][2]-p[1][1]);
      p[2][0] = p[2][1] - (p[2][2]-p[2][1]);            p[3][0] = p[3][1] - (p[3][2]-p[3][1]);
      p[0][0] = p[1][0] - (p[2][0]-p[1][0]);            p[0][1] = p[1][1] - (p[2][1]-p[1][1]);
      p[0][2] = p[1][2] - (p[2][2]-p[1][2]);            p[0][3] = p[1][3] - (p[2][3]-p[1][3]);
    }
    else if (padLeft && !padRight && !padTop && padBottom) {
      p[1][0] = val[clow-1 + mat.rows*rlow];            p[2][0] = val[clow-1 + mat.rows*(rlow+1)];
      p[3][0] = val[clow-1 + mat.rows*(rlow+2)];        p[3][1] = val[clow + mat.rows*(rlow+2)];
      p[3][2] = val[clow+1 + mat.rows*(rlow+2)];        p[1][3] = p[1][2] + (p[1][2]-p[1][1]);
      p[2][3] = p[2][2] + (p[2][2]-p[2][1]);            p[3][3] = p[2][3] + (p[2][3]-p[1][3]);
      p[0][0] = p[1][0] - (p[2][0]-p[1][0]);            p[0][1] = p[1][1] - (p[2][1]-p[1][1]);
      p[0][2] = p[1][2] - (p[2][2]-p[1][2]);            p[0][3] = p[1][3] - (p[2][3]-p[1][3]);
    }
    else if (padLeft && !padRight && padTop && !padBottom) {
      p[1][0] = val[clow-1 + mat.rows*rlow];            p[2][0] = val[clow-1 + mat.rows*(rlow+1)];
      p[0][0] = val[clow-1 + mat.rows*(rlow-1)];        p[0][1] = val[clow + mat.rows*(rlow-1)];
      p[0][2] = val[clow+1 + mat.rows*(rlow-1)];        p[0][3] = p[0][2] + (p[0][2]-p[0][1]);
      p[1][3] = p[1][2] + (p[1][2]-p[1][1]);            p[2][3] = p[2][2] + (p[2][2]-p[2][1]);
      p[3][0] = p[2][0] + (p[2][0]-p[1][0]);            p[3][1] = p[2][1] + (p[2][1]-p[1][1]);
      p[3][2] = p[2][2] + (p[2][2]-p[1][2]);            p[3][3] = p[2][3] + (p[2][3]-p[1][3]);
    }
    else {
      cerr<<"ERROR: Cannot find correct conditions to fill array of neighboring cells!!!"<<endl;
      exit(0);
    }
  }
  else {
    cerr<<"ERROR: This function has not been finished to accomodate non continuous matrices, now it's your job!"<<endl;
    exit(1);
  }

  //for (int i=0;i<4;i++) {
  //    for (int j=0;j<4;j++)     cout<<p[i][j]<<" ";
  //    cout<<endl;
  //}

  return bicubicInterpolate(p, row-rlow, col-clow);
}

#endif

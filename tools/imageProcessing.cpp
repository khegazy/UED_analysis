#include "imageProcessing.h"
#include <TCanvas.h>

using namespace std;


//enum { CV_8U=0, CV_8S=1, CV_16U=2, CV_16S=3, CV_32S=4, CV_32F=5, CV_64F=6 };


/*
template<typename type>
std::vector< std::vector<double> > imgProc::getImgVector(cv::Mat imgMat, 
    int roih, int roiw, int centerR, int centerC, 
    bool removeReadout, int roR, int roC, int window) {



  int ir, ic, icol, irow, index;
  int Nrows = (int)imgMat.rows;
  int Ncols = (int)imgMat.cols;
  type* pv = imgMat.ptr<type>(0);
  double readoutNoise = 0;
  std::vector< std::vector<double> > img(roih);

  if (removeReadout) {
    std::vector<double> readoutWindow;
    for (ir=roR; ir<roR+window; ir++) {
      if (ir >= Nrows) {
        cerr << "ERROR: readout noise window exceeds rows (" 
          << ir << ">=" << Nrows <<"!!!\n" << endl;
        exit(0);
      }
      for (ic=roC; ic<roC+window; ic++) {
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

  for (ir=0; ir<roih; ir++) {
    img[ir].resize(roiw);
    irow = ir - roih/2 + centerR;
    for (ic=0; ic<roiw; ic++) {
      icol = ic - roiw/2 + centerC;
      index = irow*Ncols + icol;
      img[ir][ic] = ((double)(pv[index])) - readoutNoise;
    }
  }

  return img;
}
*/


template<typename T>
void imgProc::medianPixelT(cv::Mat &mat, int row, int col, int hotpixel) {

  vector<int> vect(25);
  int count=0;
  int size=0;
  int index;
  T* val = mat.ptr<T>(0);

  while (count != (int)vect.size()) {
    size += 2;
    count = 0;
    for (int ir=-size/2; ir<=size/2; ir++) {
      for (int ic=-size/2; ic<=size/2; ic++) {
	index = min(max(0,row+ir),mat.rows-1)*mat.cols + min(max(col+ic,0),mat.cols-1);
	if (val[index] < hotpixel) {
	  vect[count] = (val[index]);
	  count++;
	}
	if (count == (int)vect.size()) break;
      }
      if (count == (int)vect.size()) break;
    }
  }

  sort(vect.begin(), vect.end());
  val[row*mat.cols + col] = vect[vect.size()/2];
}

void imgProc::medianPixel(cv::Mat &mat, int row, int col, int hotpixel) {

  switch (mat.depth()) {
    case 0:
	medianPixelT<uint8_t>(mat, row, col, hotpixel);
	break;
    case 1:
   	medianPixelT<int8_t>(mat, row, col, hotpixel);
	break;
    case 2:
	medianPixelT<uint16_t>(mat, row, col, hotpixel);
	break;
    case 3:
	medianPixelT<int16_t>(mat, row, col, hotpixel);
	break;
    case 4:
	medianPixelT<float>(mat, row, col, hotpixel);
	break;
    case 5:
	medianPixelT<double>(mat, row, col, hotpixel);
	break;
    default:
	cerr<<"ERROR: Do not recognize bit depth enum "<<mat.depth()<<"!!!"<<endl;
	exit(0);
  }
}


template<typename T>
void imgProc::averagePixelT(cv::Mat &mat, int row, int col) {
  
  T* val; 
  T xSum=0;
  T ySum=0;
  int Nrows = mat.rows;
  int Ncols = mat.cols;

  if (mat.isContinuous()) {
    val = mat.ptr<T>(0);
    if (col > Ncols && col < (Nrows-1)*Ncols) {ySum = val[col-Ncols]+val[col+Ncols];}
    if ((col+1)%Ncols!=0 && col%Ncols!=0) xSum = val[col+1]+val[col-1];
  
    if (xSum && ySum) val[row*col] = (xSum+ySum)/4;
    else if (xSum) val[row*col] = xSum/2;
    else val[row*col] = ySum/2;
  } 
  
  else {
    if ((col+1)<Ncols && col > 0) {
      val = mat.ptr<T>(row);
      xSum = val[col+1]+val[col-1];
    } 
    if ((row+1)<Nrows && row > 0) {
      val = mat.ptr<T>(row-1);
      ySum = val[col];
      val = mat.ptr<T>(row+1);
      ySum += val[col]; 
    } 
        
    if (xSum && ySum) val[col] = (xSum+ySum)/4;
    else if (xSum) val[col] = xSum/2;
    else val[col] = ySum/2;
  }       
}


void imgProc::averagePixel(cv::Mat &mat, int row, int col) {

  switch (mat.depth()) {
    case 0:
	averagePixelT<uint8_t>(mat, row, col);
	break;
    case 1:
   	averagePixelT<int8_t>(mat, row, col);
	break;
    case 2:
	averagePixelT<uint16_t>(mat, row, col);
	break;
    case 3:
	averagePixelT<int16_t>(mat, row, col);
	break;
    case 4:
	averagePixelT<float>(mat, row, col);
	break;
    case 5:
	averagePixelT<double>(mat, row, col);
	break;
    default:
	cerr<<"ERROR: Do not recognize bit depth enum "<<mat.depth()<<"!!!"<<endl;
	exit(0);
  }
}


std::vector<double> imgProc::gaussianSmooth1d(
      std::vector<double> inp, 
      double sigma, 
      int totalRange) {

  double norm, weight;
  std::vector<double> out(inp.size(), 0);

  int range = totalRange/2 + 1;

  for (int i=0; i<(int)inp.size(); i++) {
    if (inp[i] != NANVAL) {
      norm = 0;
      for (int j=-1*range; j<=range; j++) {
        if ((i + j < 0) || (i + j >= (int)inp.size())) {
          continue;
        }

        if (inp[i+j] != NANVAL) {
          weight = std::exp(-1*j*j/(2*std::pow(sigma, 2)));
          norm += weight;
          out[i] += inp[i+j]*weight;
        }
      }

      out[i] /= norm;
    }
    else {
      out[i] = NANVAL;
    }
  }

  return out;
}

    
///////////////////////////////////////////////////////////////////////////////////////////////////////
void imgProc::threshold(vector< vector<double> >  &img, int hotpixel) {

  for (auto &ir : img) {
    for (auto &ic: ir) if (fabs(ic) >= hotpixel) ic=0; //medianPixelT<T> (mat, ir, ic, hotpixel);
  }
}


template<typename T>
void imgProc::thresholdT(cv::Mat &mat, int hotpixel) {

  int Nrows = mat.rows; 
  int Ncols = mat.cols;
  if (mat.isContinuous()) {
    int pixels = Ncols*Nrows;
    T* valhp = mat.ptr<T>(0);
    for (int ip=0; ip<pixels; ip++) {
      if (fabs(valhp[ip]) >= hotpixel) medianPixelT<T> (mat, ip/Ncols, ip%Ncols, hotpixel);
    } 
  } 
  else {
    T* valhp;
    for (int ir=0; ir<Nrows; ir++)  {
      valhp = mat.ptr<T>(ir);
      for (int ic=0; ic<Ncols; ic++) {
        if (fabs(valhp[ic]) >= hotpixel) medianPixelT<T> (mat, ir, ic, hotpixel);
      } 
    }
  } 
}


void imgProc::threshold(cv::Mat &mat, int hotpixel) {

  switch (mat.depth()) {
    case 0:
	thresholdT<uint8_t>(mat, hotpixel);
	break;
    case 1:
   	thresholdT<int8_t>(mat, hotpixel);
	break;
    case 2:
	thresholdT<uint16_t>(mat, hotpixel);
	break;
    case 3:
	thresholdT<int16_t>(mat, hotpixel);
	break;
    case 4:
	thresholdT<float>(mat, hotpixel);
	break;
    case 5:
	thresholdT<uint32_t>(mat, hotpixel);
	break;
    default:
	cerr<<"ERROR: Do not recognize bit depth enum "<<mat.depth()<<"!!!"<<endl;
	exit(0);
  }
}



void imgProc::removeXrayHits(
    std::vector<std::vector<double> >* img,
    double highCut, double lowCut, 
    double stdCut, int window,
    TH1F** plots) {

  std::vector<double> binStats;
  for (int ir=0; ir<(int)(*img).size(); ir++) {
    for (int ic=0; ic<(int)(*img)[ir].size(); ic++) {
      if ((*img)[ir][ic] == NANVAL) {
        continue;
      }

      if ((*img)[ir][ic] > highCut) {
        for (int icc = ic; icc<(int)(*img)[ir].size(); icc++) {
          (*img)[ir][icc] = NANVAL;
        }
      }
      else if (plots) {
        plots[0]->Fill((*img)[ir][ic], 1);
      }
    }
  }

  for (int ir=0; ir<(int)(*img).size(); ir++) {
    for (int ic=0; ic<(int)(*img)[ir].size(); ic++) {
      if ((*img)[ir][ic] == NANVAL) {
        continue;
      }

      if ((*img)[ir][ic] > lowCut) {
        int irow, icol;
        binStats.clear();

        // Find non X-ray bins
        for (int irr=-1*window; irr<=window; irr++) {
          irow = ir+irr;
          if (irow >= 1024 || irow < 0) {
            continue;
          }
          for (int icc=-1*window; icc<=window; icc++) {
            icol = ic+icc;
            if (icol >= 1024 || icol < 0) {
              continue;
            }

            if (((*img)[irow][icol] != NANVAL) 
                && ((*img)[irow][icol] < lowCut)) {
              binStats.push_back((*img)[irow][icol]);
            }
          }
        }

        double mean = std::accumulate(binStats.begin(), binStats.end(), 0.0);
        mean /= binStats.size();

        double stdev = 0;
        for (uint i=0; i<binStats.size(); i++) {
          stdev += std::pow(binStats[i] - mean, 2);
        }
        stdev = std::sqrt(stdev/binStats.size());

        // Find non X-ray bins
        for (int irr=-1*window; irr<=window; irr++) {
          irow = ir+irr;
          if (irow >= 1024 || irow < 0) {
            continue;
          }
          for (int icc=-1*window; icc<=window; icc++) {
            icol = ic+icc;
            if (icol >= 1024 || icol < 0) {
              continue;
            }

            if (((*img)[irow][icol] != NANVAL)
                && (fabs((*img)[irow][icol]-mean) > stdCut*stdev)) {
              (*img)[irow][icol] = NANVAL;
            }
          }
        }
      }
      else if (plots) {
        plots[1]->Fill((*img)[ir][ic], 1);
      }
    }
  }
}



void imgProc::histMeanStd(vector< vector<double> >& hist, double& mean, double& std) {

  mean = std = 0;
  for (uint ir=0; ir<hist.size(); ir++) {
    for (uint ic=0; ic<hist[ir].size(); ic++) {
      mean += hist[ir][ic];
      std += pow(hist[ir][ic], 2);
    }
  }

  mean /= (double)(hist.size()*hist[0].size());
  std = sqrt(std/((double)(hist.size()*hist[0].size())) - mean*mean);

  return;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////


double imgProc::removeReadOutNoise(std::vector< std::vector<double> > &img, 
    int roiR, int roiC, int window, int buffer) {

  /// Calculate image norm ///
  int rows = (int)img.size();
  int cols = (int)img[0].size();
  int ir, ic;
  std::vector<double> vals;

  if (roiR+window > rows) {
    cerr << "ERROR: readout noise window exceeds rows!!!\n" << endl;
    exit(0);
  }
  if (roiC+window > cols) {
    cerr << "ERROR: readout noise window exceeds cols!!!\n" << endl;
    exit(0);
  }

  for (ir=roiR; ir<roiR+window; ir++) {
    if ((ir < buffer) || (ir + buffer >= rows)) continue;
    for (ic=roiC; ic<roiC+window; ic++) {
      if ((ic < buffer) || (ic + buffer >= cols)) continue;
      if (img[ir][ic] != NANVAL) {
        vals.push_back(img[ir][ic]);
      }
    }
  }

  std::vector<int> orderedInds(vals.size(), 0);
  iota(orderedInds.begin(), orderedInds.end(), 0);
  sort(orderedInds.begin(), orderedInds.end(),
        [&vals](uint i1, uint i2)
        {return vals[i1] < vals[i2];});

  double readoutNoise = vals[orderedInds[(uint)(vals.size()/2.)]];

  for (ir=0; ir<rows; ir++) {
    for (ic=0; ic<cols; ic++) {
      if (img[ir][ic] != NANVAL) {
        img[ir][ic] -= readoutNoise;
      }
    }
  }

  return readoutNoise;
}


double imgProc::centerfnctr::operator() (std::vector<double> vect) {
  if (fxnType == 0)  {
    if (!img) {
      std::cerr 
        << "ERROR: Must specifiy img in centerfnctr to use this option!!!\n";
      exit(1);
    }
    return imgProc::centerSymXsqr(
        img, vect[0], vect[1], 
        8, 1, minRadBin, centShellW, true);
  }
  if (fxnType == 1)  {
    if (!img) {
      std::cerr 
        << "ERROR: Must specifiy img in centerfnctr to use this option!!!\n";
      exit(1);
    }
    return imgProc::centerOddLeg(
        img, vect[0], vect[1], 
        1, minRadBin, centShellW, true);
  }
  if (fxnType == 2)  {
    if (!radProc) {
      std::cerr 
        << "ERROR: Must specifiy radProc in centerfnctr to use this option!!!\n";
      exit(1);
    }    
    if (!img) {
      std::cerr 
        << "ERROR: Must specifiy img in centerfnctr to use this option!!!\n";
      exit(1);
    }

    return radProc->radialSliceVar(
        img, vect[0], vect[1], 
        minRadBin, centShellW);
  }
  if (fxnType == 3)  {
    double loss = 0;
    std::vector<float> dists(indsR->size());
    for (uint i=0; i<indsR->size(); i++) {
      dists[i] = std::sqrt(
          std::pow((*indsR)[i] - vect[0], 2)
          + std::pow((*indsC)[i] - vect[1], 2));
    }

    double meanR = std::accumulate(dists.begin(), dists.end(), 0);
    meanR /= (float)dists.size();

    for (uint i=0; i<indsR->size(); i++) {
      loss += std::pow(meanR - dists[i], 2)
                *std::exp(-0.5*std::pow(((*vals)[i] - meanVal)/valSTD, 2));
    }
    loss /= (float)dists.size();

    return loss;
  }
  else {
    cerr << "ERROR: Did not select a correct center finding alogrithm, now exiting!!!\n\n";
    exit(0);
  }
  //N2 return imgProc::centerOddLeg(img, vect[0], vect[1], 1, 70, 40, true, NANVAL);
}


/*
///////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> imgProc::findCenterSingleImg(imgInfoStruct info, 
    double hotPixel, double sigma,
    int blockCentR, int blockCentC, float minRad, float maxRad, 
    int meanInd, double stdScale, bool verbose, PLOTclass* pltVerbose) {

  // Create initial matrix
  cv::Mat imgSum = cv::imread((imgINFO[0].path + imgINFO[0].fileName).c_str(),
      CV_LOAD_IMAGE_ANYDEPTH);
  imgSum.convertTo(imgSum, 5);
  int Nrows = imgSum.rows;
  int Ncols = imgSum.cols;

  /////  Average images  /////
  for (uint ifl=1; ifl<imgINFO.size(); ifl++) {

    string imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH);
    imgProc::threshold(imgMat, hotPixel);
    cv::add(imgSum, imgMat, imgSum, cv::Mat(), 5);
    if (verbose) cout << "\tpassed!\n";
  }
  imgSum /= imgINFO.size();

  /////  Smooth averaged image  /////
  if (verbose) cout << "Smoothing  ...  ";
  cv::Mat smoothMat;
  int size = (int)(sigma*6);
  size += 1 - size%2;
  cv::GaussianBlur(imgSum, smoothMat, cvSize(size,size), sigma, sigma);

  // Convert to vector
  if (verbose) cout << "Vectorizing  ...  ";
  std::vector< std::vector<double> > imgCent 
        = imgProc::getImgVector<float>(smoothMat, Nrows, Nrows/2, Ncols/2);
  if (pltVerbose) save::saveDat<double>(imgCent, "./results/centering_smoothed[" 
      + to_string(Nrows) + "," + to_string(Ncols) + "].dat");

  // Remove readout noise
  if (verbose) cout << "Remove readout  ...  ";
  imgProc::removeReadOutNoise(imgCent);
  if (verbose) cout << "done\n";

  /////  Find center by averaging points  /////

  std::vector< std::vector<double> > testMap(Nrows);
  if (verbose || pltVerbose) {
    testMap.resize(Nrows);
    for (uint ir=0; ir<Nrows; ir++) {
      testMap[ir].resize(Ncols,0);
    }
  }

  double weight, rad;
  double count = 0;
  double centerR = 0;
  double centerC = 0;
  double scale = imgCent[(int)(Nrows/2)][(int)(Ncols/2)];
  double val = imgCent[meanInd][meanInd]*scale;
  for (int ir=0; ir<(int)imgCent.size(); ir++) {
    for (int ic=0; ic<(int)imgCent[ir].size(); ic++) {
      imgCent[ir][ic] *= scale;
      rad = std::sqrt(std::pow(ir-blockCentR, 2) + std::pow(ic-blockCentC, 2));
      if ((rad < minRad) || (rad > maxRad)) continue;

      weight = std::exp(-1*std::pow((val - imgCent[ir][ic])/val, 2)/stdScale);
      //if (verbose) cout << "vals: " << val << "/" << imgCent[ir][ic] <<"    "<<std::pow((val - imgCent[ir][ic])/val, 2)/stdScale
       // << "\tweight: "<<weight<<endl;

      if (weight > 0.1) {
        centerR += ir*weight;
        centerC += ic*weight;
        count += weight;

        if (pltVerbose) {
          testMap[ir][ic] = weight;
        }
      }
    }
  }

  if (pltVerbose) {
    save::saveDat<double>(testMap, "./results/centering_selectedRange["
      + to_string(Nrows) + "," + to_string(Ncols) + "].dat");
    pltVerbose->printRC(testMap, "plots/centering_selectedRange");
  }
  std::vector<int> result = {(int)(centerR/count), (int)(centerC/count)};

  if (pltVerbose) {
    cout<<"here"<<endl;
    TH2F* cImg = pltVerbose->plotRC(imgCent, "plots/roughCenteredImage");
    for (int ir=(result[0])-5; ir<=result[0]+5; ir++) {
      for (int ic=(result[1])-5; ic<=result[1]+5; ic++) {
        if (sqrt(pow(ir-result[0],2)+pow(ic-result[1],2)) < 5)
          cImg->SetBinContent(ic, ir, -100000);
      }
    }
    cImg->SetMinimum(-1);
    pltVerbose->print2d(cImg, "plots/roughCenteredImage");
    delete cImg;
  }

  return result;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////


double imgProc::removeAvgReadOutNoise(std::vector< std::vector<double> > &img, 
    int centerR, int centerC, double minRad, double maxRad, int buffer, 
    std::vector< std::vector<double> >* nanMap) {

  /// Calculate image norm ///
  int Nrows = (int)img.size();
  int Ncols = (int)img[0].size();
  int ir, ic;
  double count = 0;
  double sum = 0;
  double rad;

  for (ir=0; ir<Nrows; ir++) {
    if ((ir < buffer) || (ir + buffer >= Nrows)) continue;
    for (ic=0; ic<Ncols; ic++) {
      if ((ic < buffer) || (ic + buffer >= Ncols)) continue;
      rad = std::sqrt(std::pow(ir-centerR, 2) + std::pow(ic-centerC, 2));
      if (rad > minRad && rad<maxRad ) {
        if ((img[ir][ic] != NANVAL)) {
          if (nanMap) {
            if ((*nanMap)[ir][ic] == NANVAL) {
              continue;
            }
          }
          sum += img[ir][ic];
          count++;
        }
      }
    }
  }

  double readoutNoise = sum/count;

  sum = 0;
  for (ir=0; ir<Nrows; ir++) {
    for (ic=0; ic<Ncols; ic++) {
      if (img[ir][ic] != NANVAL) {
        img[ir][ic] -= readoutNoise;
      }
    }
  }

  return readoutNoise;
}


double imgProc::removeMedianReadOutNoise(std::vector< std::vector<double> > &img, 
    int centerR, int centerC, double minRad, double maxRad, int buffer, 
    std::vector< std::vector<double> >* nanMap) {

  /// Find image norm ///
  std::vector<double> vals;
  int Nrows = (int)img.size();
  int Ncols = (int)img[0].size();
  int ir, ic;
  double rad;

  for (ir=0; ir<Nrows; ir++) {
    if ((ir < buffer) || (ir + buffer >= Nrows)) continue;
    for (ic=0; ic<Ncols; ic++) {
      if ((ic < buffer) || (ic + buffer >= Ncols)) continue;
      rad = std::sqrt(std::pow(ir-centerR, 2) + std::pow(ic-centerC, 2));
      if (rad > minRad && rad<maxRad ) {
        if ((img[ir][ic] != NANVAL)) {
          if (nanMap) {
            if ((*nanMap)[ir][ic] == NANVAL) {
              continue;
            }
          }
          vals.push_back(img[ir][ic]);
        }
      }
    }
  }

  vector<int> orderedInds(vals.size());
  iota(orderedInds.begin(), orderedInds.end(), 0);
  sort(orderedInds.begin(), orderedInds.end(),
      [&vals](int i1, int i2)
      {return vals[i1] < vals[i2];});

  double readoutNoise = vals[orderedInds[vals.size()/2]];

  for (ir=0; ir<Nrows; ir++) {
    for (ic=0; ic<Ncols; ic++) {
      if (img[ir][ic] != NANVAL) {
        img[ir][ic] -= readoutNoise;
      }
    }
  }

  return readoutNoise;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<int> imgProc::centerSearchCOM(std::vector<imgInfoStruct> &imgINFO, 
    double hotPixel, double sigma,
    int blockCentR, int blockCentC, float minRad, float maxRad, 
    int meanInd, double stdScale, bool verbose, PLOTclass* pltVerbose) {

  // Create initial matrix
  cv::Mat imgSum = cv::imread((imgINFO[0].path + imgINFO[0].fileName).c_str(),
      CV_LOAD_IMAGE_ANYDEPTH);
  imgSum.convertTo(imgSum, 5);
  int Nrows = imgSum.rows;
  int Ncols = imgSum.cols;

  /////  Average images  /////
  for (uint ifl=1; ifl<imgINFO.size(); ifl++) {

    string imgAddr = imgINFO[ifl].path + imgINFO[ifl].fileName;
    if (verbose) cout << "INFO: Trying to open " << imgAddr << "\t .....";
    cv::Mat imgMat = cv::imread(imgAddr.c_str() ,CV_LOAD_IMAGE_ANYDEPTH);
    imgProc::threshold(imgMat, hotPixel);
    cv::add(imgSum, imgMat, imgSum, cv::Mat(), 5);
    if (verbose) cout << "\tpassed!\n";
  }
  imgSum /= imgINFO.size();

  /////  Smooth averaged image  /////
  if (verbose) cout << "Smoothing  ...  ";
  cv::Mat smoothMat;
  int size = (int)(sigma*6);
  size += 1 - size%2;
  cv::GaussianBlur(imgSum, smoothMat, cvSize(size,size), sigma, sigma);

  // Convert to vector
  if (verbose) cout << "Vectorizing  ...  ";
  std::vector< std::vector<double> > imgCent 
        = imgProc::getImgVector<float>(smoothMat, Nrows, Nrows/2, Ncols/2);
  if (pltVerbose) save::saveDat<double>(imgCent, "./results/centering_smoothed[" 
      + to_string(Nrows) + "," + to_string(Ncols) + "].dat");

  // Remove readout noise
  if (verbose) cout << "Remove readout  ...  ";
  imgProc::removeReadOutNoise(imgCent);
  if (verbose) cout << "done\n";

  /////  Find center by averaging points  /////

  std::vector< std::vector<double> > testMap(Nrows);
  if (verbose || pltVerbose) {
    testMap.resize(Nrows);
    for (int ir=0; ir<Nrows; ir++) {
      testMap[ir].resize(Ncols,0);
    }
  }

  double weight, rad;
  double count = 0;
  double centerR = 0;
  double centerC = 0;
  double scale = imgCent[(int)(Nrows/2)][(int)(Ncols/2)];
  double val = imgCent[meanInd][meanInd]*scale;
  for (int ir=0; ir<(int)imgCent.size(); ir++) {
    for (int ic=0; ic<(int)imgCent[ir].size(); ic++) {
      imgCent[ir][ic] *= scale;
      rad = std::sqrt(std::pow(ir-blockCentR, 2) + std::pow(ic-blockCentC, 2));
      if ((rad < minRad) || (rad > maxRad)) continue;

      weight = std::exp(-1*std::pow((val - imgCent[ir][ic])/val, 2)/stdScale);
      //if (verbose) cout << "vals: " << val << "/" << imgCent[ir][ic] <<"    "<<std::pow((val - imgCent[ir][ic])/val, 2)/stdScale
       // << "\tweight: "<<weight<<endl;

      if (weight > 0.1) {
        centerR += ir*weight;
        centerC += ic*weight;
        count += weight;

        if (pltVerbose) {
          testMap[ir][ic] = weight;
        }
      }
    }
  }

  if (pltVerbose) {
    save::saveDat<double>(testMap, "./results/centering_selectedRange["
      + to_string(Nrows) + "," + to_string(Ncols) + "].dat");
    pltVerbose->printRC(testMap, "plots/centering_selectedRange");
  }

  // Can use for messier data
  /*
  int shift=15;
  double rad;
  for (uint k=0; k<4; k++) {
    for (ir=shift; ir<imgCent.size()-shift; ir++) {
      for (ic=shift; ic<imgCent[ir].size()-shift; ic++) {
        rad = std::sqrt(std::pow(ic-Ncols/2, 2) + std::pow(ir-Nrows/2, 2));
        if (rad < 100) continue;

        count = 0;

        for (int irr=-1*shift; irr<=shift; irr++) {
          for (int icc=-1*shift; icc<=shift; icc++) {
            rad = std::sqrt(irr*irr + icc*icc);
            if (rad < shift) {
              count += testMap[ir+irr][ic+icc]*
                std::exp(-1*rad*rad/(shift/4));
            }
          }
        }
        if (count > 0.5 + k*0.15) {
          finMap[ir][ic] = count;
        }
      }
    }

            plt.printRC(testMap, "testmap"+to_string(k));
    for (ir=0; ir<Nrows; ir++) {
      for (ic=0; ic<Ncols; ic++) {
        testMap[ir][ic] = finMap[ir][ic];
      }
    }
  }
  */

  std::vector<int> result = {(int)(centerR/count), (int)(centerC/count)};

  if (pltVerbose) {
    cout<<"here"<<endl;
    TH2F* cImg = pltVerbose->plotRC(imgCent, "plots/roughCenteredImage");
    for (int ir=(result[0])-5; ir<=result[0]+5; ir++) {
      for (int ic=(result[1])-5; ic<=result[1]+5; ic++) {
        if (sqrt(pow(ir-result[0],2)+pow(ic-result[1],2)) < 5)
          cImg->SetBinContent(ic, ir, -100000);
      }
    }
    cImg->SetMinimum(-1);
    pltVerbose->print2d(cImg, "plots/roughCenteredImage");
    delete cImg;
  }

  return result;
}


std::vector<int> imgProc::centerSearchCOM(
    std::vector< std::vector<double> > imgCent, 
    int blockCentR, int blockCentC, float minRad, float maxRad, 
    int meanInd, double stdScale, bool verbose, PLOTclass* pltVerbose) {


  /////  Find center by averaging points  /////
  uint Nrows = imgCent.size();
  uint Ncols = imgCent[0].size();
  double weight, rad;
  double count = 0;
  double centerR = 0;
  double centerC = 0;
  double scale = imgCent[(int)(Nrows/2)][(int)(Ncols/2)];
  double val = imgCent[meanInd][meanInd]*scale;
  
  std::vector< std::vector<double> > testMap(Nrows);
  if (verbose || pltVerbose) {
    testMap.resize(Nrows);
    for (uint ir=0; ir<Nrows; ir++) {
      testMap[ir].resize(Ncols,0);
    }
  }

  for (int ir=0; ir<(int)imgCent.size(); ir++) {
    for (int ic=0; ic<(int)imgCent[ir].size(); ic++) {
      imgCent[ir][ic] *= scale;
      rad = std::sqrt(std::pow(ir-blockCentR, 2) + std::pow(ic-blockCentC, 2));
      if ((rad < minRad) || (rad > maxRad)) continue;

      weight = std::exp(-1*std::pow((val - imgCent[ir][ic])/val, 2)/stdScale);

      if (weight > 0.1) {
        centerR += ir*weight;
        centerC += ic*weight;
        count += weight;

        if (pltVerbose) {
          testMap[ir][ic] = weight;
        }
      }
    }
  }

  if (pltVerbose) {
    save::saveDat<double>(testMap, "./results/centering_selectedRange["
      + to_string(Nrows) + "," + to_string(Ncols) + "].dat");
    pltVerbose->printRC(testMap, "plots/centering_selectedRange");
  }

  std::vector<int> result = {(int)(centerR/count), (int)(centerC/count)};

  if (pltVerbose) {
    TH2F* cImg = pltVerbose->plotRC(imgCent, "plots/roughCenteredImage");
    for (int ir=(result[0])-5; ir<=result[0]+5; ir++) {
      for (int ic=(result[1])-5; ic<=result[1]+5; ic++) {
        if (sqrt(pow(ir-result[0],2)+pow(ic-result[1],2)) < 5)
          cImg->SetBinContent(ic, ir, -100000);
      }
    }
    cImg->SetMinimum(-1);
    pltVerbose->print2d(cImg, "plots/roughCenteredImage");
    delete cImg;
  }

  return result;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////


double imgProc::imageNorm(std::vector< std::vector<double> > &img, 
    float radMin, float radMax) {

  double imgNorm = 0;
  double count = 0;
  double rad;
  float rows = (int)img.size();
  float cols = (int)img[0].size();
  for (int ir=0; ir<rows; ir++) { 
    for (int ic=0; ic<cols; ic++) {
      rad = sqrt(pow(ir-rows/2,2) + pow(ic-cols/2,2))/(rows/2);
      if ((rad > radMin) && (rad < radMax) && (img[ir][ic] != NANVAL)) {
        imgNorm += img[ir][ic];
        count++;
      }
    }
  }
  
  return imgNorm /= count;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
//  Use 2D sine transformation to fine asymmetric part of the image  //
///////////////////////////////////////////////////////////////////////
std::vector< std::vector<double> > imgProc::asymmetrize(
        std::vector< std::vector<double> > inpImg,
        int centerR, int centerC,
        int height, int width,
        std::vector< std::vector<double> > oddImgImgn,
        fftw_plan &fftFref, fftw_complex* fftIn,
        fftw_plan &fftBref, fftw_complex* fftOut,
        PLOTclass* plt) {

  int ir, ic;
  int irow, icol;

  // Real and imaginary parts of final image
  oddImgImgn.resize(height);
  std::vector< std::vector<double> > oddImgReal(height);
  for (ir=0; ir<height; ir++) {
    oddImgReal[ir].resize(width, 0);
    oddImgImgn[ir].resize(width, 0);
  }

  // Filling fftIn array
  std::vector< std::vector<double> > testImg(height);
  for (ir=0; ir<height; ir++) {
    testImg[ir].resize(width);
    irow = centerR + ir;
    if (irow >= centerR + height/2) irow = irow%(centerR + height/2) - height/2 + centerR;
    for (ic=0; ic<width; ic++) {
      icol = centerC + ic; 
      if (icol >= centerC + width/2) icol = icol%(centerC + width/2) - width/2 + centerC;
      fftIn[ir*width+ic][0] = inpImg[irow][icol];
      fftIn[ir*width+ic][1] = 0;

      if (plt) {
        testImg[ir][ic] = fftIn[ir*width+ic][0];
      }
    }
  }

  if (plt) {
    plt->printRC(testImg,"plots/testImg", maximum, "1000");
  }

  // Execute forward fft
  fftw_execute(fftFref);

  // Only keep imaginary terms for sine transform
  for (ir=0; ir<height*width; ir++) {
    fftIn[ir][1] = fftOut[ir][1];
    fftIn[ir][0] = 0;
  }

  // Execute backwards fft
  fftw_execute(fftBref);

  // Fill results
  for (ir=0; ir<height; ir++) {
    for (ic=0; ic<width; ic++) {
      irow = ((int)(height/2)+ir)%height;
      icol = ((int)(width/2)+ic)%width;

      oddImgImgn[irow][icol] = fftOut[ir*width+ic][1]/(height*width);
      if (fftOut[ir*width+ic][0] > 0) {
        oddImgReal[irow][icol] = 2*fftOut[ir*width+ic][0]/(height*width);
      }
    }
  }

  return oddImgReal;
}


/*
///////////////////////////////////
//  Filter diffraction pattern.  //
///////////////////////////////////
std::vector< std::vector<double> > imgProc::filter2Ddiffraction(
        std::vector< std::vector<double> > inpImg,
        int centerR, int centerC,
        int height, int width,
        std::vector<double>&  filter,
        fftw_plan &fftFref, fftw_complex* fftIn,
        fftw_plan &fftBref, fftw_complex* fftOut) {

  int ir, ic;
  int irow, icol;

  // Real and imaginary parts of final image
  oddImgImgn.resize(height);
  std::vector< std::vector<double> > oddImgReal(height);
  for (ir=0; ir<height; ir++) {
    oddImgReal[ir].resize(width, 0);
    oddImgImgn[ir].resize(width, 0);
  }

  // Filling fftIn array
  for (ir=0; ir<height; ir++) {
    irow = centerR + ir;
    if (irow >= centerR + height/2) irow = irow%(centerR + height/2) - height/2 + centerR;
    for (ic=0; ic<width; ic++) {
      icol = centerC + ic;
      if (icol >= centerC + width/2) icol = icol%(centerC + width/2) - width/2 + centerC;
      fftIn[ir*width+ic][0] = inpImg[irow][icol];
      fftIn[ir*width+ic][1] = 0;
    }
  }

  // Execute forward fft
  fftw_execute(fftFref);

  // Only keep imaginary terms for sine transform
  for (ir=0; ir<height; ir++) {
    irow = centerR + ir;
    if (irow >= centerR + height/2) irow = irow%(centerR + height/2) - height/2 + centerR;
    irSq = std::pow(ir%(int(height/2)), 2);
    for (ic=0; ic<width; ic++) {
      icol = centerC + ic;
      if (icol >= centerC + width/2) icol = icol%(centerC + width/2) - width/2 + centerC;
      filterInd = filter.size()*(sqrt(irSq + std::pow(icol - centerC, 2))/(height/2 + 0.5));
      fftIn[ir*width+ic][0] = filter[filterInd]*fftOut[ir*width+ic][0];
      fftIn[ir*width+ic][1] = 0;
    }
  }

  // Execute backwards fft
  fftw_execute(fftBref);

  // Fill results
  for (ir=0; ir<height; ir++) {
    for (ic=0; ic<width; ic++) {
      irow = ((int)(height/2)+ir)%height;
      icol = ((int)(width/2)+ic)%width;

      oddImgImgn[irow][icol] = fftOut[ir*width+ic][1]/(height*width);
      if (fftOut[ir*width+ic][0] > 0) {
        oddImgReal[irow][icol] = 2*fftOut[ir*width+ic][0]/(height*width);
      }
    }
  }

  return oddImgReal;
}
*/


///////////////////////////////////////////////////////////////////////////////////
//  Use varient of density clustering casted as a 1/r potential to find clusters //
//                                                                               //
//    1) First find initial core points s.t. "grav" potential is > coreValThresh //
//    2) Look at core pixels and trace out clusters s.t. each pixel < clusterRad //
//          from another pixel in the cluster                                    //
//          if size > minCluster pixel then pixels belong in a cluster           //
//          if minPixelSize < size < minCluster size => not a cluster but still  //
//                consider pixels (label as border pixels)                       //
//          if size < minPixelSize => Do not consider these pixels               //
//    3) Starting from only cluster core pixels, calculate neighboring pixels'   //
//          potential and if > borderValThresh then add to cluster and look      //
//          at neighbors of new border pixel.                                    //  
//    4) Add padding out to padRad from clusters.                                //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////
std::vector< std::pair<uint, uint> > imgProc::findClusters(
        std::vector< std::vector<double> > &inpImg, 
        int centerR, int centerC, double minCenterRad,
        double coreValThresh, int coreRad,
        int minClusterSize, int minPixelSize, 
        double minDensity,
        int clusterRad, double borderValThresh, 
        int borderDistLimit, int borderRad, int padRad, 
        PLOTclass* pltVerbose) {

 std::vector< std::vector<double> > clusters_coreBorder;
 return findClusters(inpImg, centerR, centerC, minCenterRad, coreValThresh, coreRad, 
                minClusterSize, minPixelSize, minDensity, clusterRad,
                borderValThresh, borderDistLimit, borderRad, padRad, clusters_coreBorder, pltVerbose);
}


std::vector< std::pair<uint, uint> > imgProc::findClusters(
        std::vector< std::vector<double> > &inpImg, 
        int centerR, int centerC, double minCenterRad,
        double coreValThresh, int coreRad,
        int minClusterSize, int minPixelSize, 
        double minDensity,
        int clusterRad, double borderValThresh, 
        int borderDistLimit, int borderRad, int padRad,
        std::vector< std::vector<double> > &clusters_coreBorder,
        PLOTclass* pltVerbose) {

 return findClusters(inpImg, centerR, centerC, minCenterRad, coreValThresh, coreRad, 
                minClusterSize, minPixelSize, minDensity, 0, 0, clusterRad,
                borderValThresh, borderDistLimit, borderRad, padRad, 0, 0, 0, 0, clusters_coreBorder, pltVerbose);
}

  
std::vector< std::pair<uint, uint> > imgProc::findClusters(
        std::vector< std::vector<double> > &inpImg, 
        int centerR, int centerC, double minCenterRad,
        double coreValThresh, int coreRad,
        int minClusterSize, int minPixelSize, 
        double minDensity, double shapeVarRatio, double shapeEdgeRatio,
        int clusterRad, double borderValThresh, 
        int borderDistLimit, int borderRad, int padRad,
        float rMinScale, float rMaxScale, float cMinScale, float cMaxScale,
        std::vector< std::vector<double> > &clusters_coreBorder,
        PLOTclass* pltVerbose) {


  map<int, map<int, double> > testVals;
  for (uint ir=970; ir<990; ir++) {
    for (uint ic=260; ic<300; ic++) {
      if (inpImg[ir][ic])
      testVals[ir][ic] = inpImg[ir][ic];
    }
  }

  double irSqrd, energy, rad, sum, Nchecks;
  double rMin, rMax, cMin, cMax;
  int row, col;
  int ir, ic, irr, icc, irrr, iccc;
  int Nrows = inpImg.size();
  int Ncols = inpImg[0].size();

  // Value of cores before cluster finding
  int initCoreVal = 10;
  // Value of cores after cluster finding
  int coreVal = 3;
  // Value of border pionts and non-cluster cores
  int borderVal = 2;
  // Value of padding for clusters
  int padVal = 1;

  // Store initial core points
  std::vector< std::pair<uint, uint> > corePairs;
  // Store all points that are in clusters (core points + border)
  std::vector< std::pair<uint, uint> > clusterPairs;
  // Store all points that need to be removed
  std::vector< std::pair<uint, uint> > removePairs;
  // Map showing core, border, and padding points
  std::vector< std::vector<double> > energyMap;
  clusters_coreBorder.resize(Nrows);
  energyMap.resize(Nrows);
  for (ir=0; ir<Nrows; ir++) {
    clusters_coreBorder[ir].resize(Ncols, 0);
    energyMap[ir].resize(Ncols, 0);
  }


  /////  Searching for initial cores  /////
  for (ir=0; ir<Nrows; ir++) {
    for (ic=0; ic<Ncols; ic++) {
      sum = 0;
      if (std::sqrt(pow(ir-centerR, 2) + pow(ic-centerC, 2)) < minCenterRad) continue; 

      for (irr=0; irr<coreRad; irr++) {
        irSqrd = irr*irr;

      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<"Looking at "<<ir<<" "<<ic<<"   "<<endl;
      //}
        // Looping columns
        for (icc=0; icc<coreRad; icc++) {
          if ((icc == 0) && (irr == 0)) continue;
          rad = std::sqrt(irSqrd + icc*icc);

          if (rad < coreRad) {
            // Looking at rows above
            if (ir + irr < Nrows) {
              // Check if within radius
              if (ic + icc < Ncols) {
                sum += inpImg[ir+irr][ic+icc]/rad;
      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<"\tar1"<<ir<<" "<<ic<<"   "<<sum<<"  "<<inpImg[ir+irr][ic+icc]<<"  "<<rad<<endl;
      //}
              }
              if (ic - icc >= 0) {
                sum += inpImg[ir+irr][ic-icc]/rad;
      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<"\tar2"<<ir<<" "<<ic<<"   "<<sum<<"  "<<inpImg[ir+irr][ic-icc]<<"  "<<rad<<endl;
      //}
              }
            }
            // Looking at rows below
            if (ir - irr >= 0) {
              // Check if within radius
              if (ic + icc < Ncols) {
                sum += inpImg[ir-irr][ic+icc]/rad;
      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<"\tab1"<<ir<<" "<<ic<<"   "<<sum<<"  "<<inpImg[ir-irr][ic+icc]<<"  "<<rad<<endl;
      //}
              }
              if (ic - icc >= 0) {
                sum += inpImg[ir-irr][ic-icc]/rad;
      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<"\tab1"<<ir<<" "<<ic<<"   "<<sum<<"  "<<inpImg[ir-irr][ic-icc]<<"  "<<rad<<endl;
      //}
              }
            }
          }
        }
      }

      energy = sum*inpImg[ir][ic];
      energyMap[ir][ic] = energy;
      //if ((ir>420) && (ir<450) && (ic>320) && (ic<350)) {
      //  cout<<ir<<" "<<ic<<"   "<<energy<<"  "<<sum<<"  "<<inpImg[ir][ic]<<endl<<endl<<endl;
     // }

      //if ((ic>630) && (ic<660) && (ir>380) && (ir<410)) {
      //  cout<<ir<<" "<<ic<<"   "<<energy<<"  "<<sum<<"  "<<inpImg[ir][ic]<<endl;
      //}

      if (energy > coreValThresh) {
        clusters_coreBorder[ir][ic] = initCoreVal;
        std::pair<uint, uint> coreInds(ir, ic);
        corePairs.push_back(coreInds);
      }
      else {
        clusters_coreBorder[ir][ic] = 0;
      }
    }
  }
  removePairs.insert(removePairs.end(), corePairs.begin(), corePairs.end());


  if (pltVerbose) {
    pltVerbose->printRC(energyMap, "imgEnergyMap", maximum, to_string(1));
    pltVerbose->printRC(clusters_coreBorder, "imgCoreMap");
  }


  /////  Classify as initial core points based upon size of cluster  /////
  /////       clusters that are too small -> label as border points  /////
  vector< std::pair<uint, uint> > tempPairs;
  int tempVal = 20;
  double rMean, cMean, rVar, cVar, count;
  for (uint ip=0; ip<corePairs.size(); ip++) {
    sum = 0;
    Nchecks = 0;
    tempPairs.clear();
    row = (int)corePairs[ip].first;
    col = (int)corePairs[ip].second;
    rMax = rMin = (int)corePairs[ip].first;
    cMax = cMin = (int)corePairs[ip].second;

    if (clusters_coreBorder[row][col] != initCoreVal) continue;
    // Count initial pixel
    sum++;
    Nchecks++;
    std::pair<uint, uint> coreInds(row, col);
    tempPairs.push_back(coreInds);

    for (uint ipp=0; ipp<tempPairs.size(); ipp++) {
      row = (int)tempPairs[ipp].first;
      col = (int)tempPairs[ipp].second;
      for (ir=0; ir<=clusterRad; ir++) {
        irSqrd = ir*ir;
        if ((row + ir >= Nrows) || (row + ir < 0)) continue;

        for (ic=0; ic<=clusterRad; ic++) {
          if ((col + ic >= Ncols) || (col + ic < 0)) continue;

          // Initial point is already in cluster
          if ((ir == 0) && (ic == 0)) continue;

          // Must be in cluster radius
          rad = sqrt(irSqrd + ic*ic);
          if (rad <= clusterRad) {
            // Looking at rows above
            if (row + ir < Nrows) {
              if (col + ic < Ncols) {
                // Avoid center region with large values
                if (sqrt(pow(row+ir-centerR, 2) + pow(col+ic-centerC, 2))
                    > minCenterRad) {
                  // Fill vector tempPairs: contains pairs in minicluster
                  if (clusters_coreBorder[row+ir][col+ic] == initCoreVal) {
                    sum++;
                    Nchecks++;
                    clusters_coreBorder[row+ir][col+ic] = tempVal;
                    std::pair<uint, uint> coreInds(row+ir, col+ic);
                    tempPairs.push_back(coreInds);

                    if (row + ir > rMax) rMax = row + ir;
                    if (col + ic > cMax) cMax = col + ic;
                  }
                  else if (clusters_coreBorder[row+ir][col+ic] == 0) {
                    Nchecks++;
                  }
                }
                //else {
                //  clusters_coreBorder[row+ir][col+ic] = coreVal;
                //}
              }
              if (col - ic >= 0) {
                // Avoid center region with large values
                if (sqrt(pow(row+ir-centerR, 2) + pow(col-ic-centerC, 2))
                    > minCenterRad) {
                  // Fill vector tempPairs: contains pairs in minicluster
                  if (clusters_coreBorder[row+ir][col-ic] == initCoreVal) {
                    sum++;
                    Nchecks++;
                    clusters_coreBorder[row+ir][col-ic] = tempVal;
                    std::pair<uint, uint> coreInds(row+ir, col-ic);
                    tempPairs.push_back(coreInds);

                    if (row + ir > rMax) rMax = row + ir;
                    if (col - ic < cMin) cMin = col - ic;
                  }
                  else if (clusters_coreBorder[row+ir][col-ic] == 0) {
                    Nchecks++;
                  }
                }
                //else {
                //  clusters_coreBorder[row+ir][col-ic] = coreVal;
                //}
              }
            }
            // Looking at rows below
            if (row - ir >= 0) {
              if (col + ic < Ncols) {
                // Avoid center region with large values
                if (sqrt(pow(row-ir-centerR, 2) + pow(col+ic-centerC, 2))
                    > minCenterRad) {
                  // Fill vector tempPairs: contains pairs in minicluster
                  if (clusters_coreBorder[row-ir][col+ic] == initCoreVal) {
                    sum++;
                    Nchecks++;
                    clusters_coreBorder[row-ir][col+ic] = tempVal;
                    std::pair<uint, uint> coreInds(row-ir, col+ic);
                    tempPairs.push_back(coreInds);

                    if (col + ic > cMax) cMax = col + ic;
                    if (row - ir < rMin) rMin = row - ir;
                  }
                  else if (clusters_coreBorder[row-ir][col+ic] == 0) {
                    Nchecks++;
                  }
                }
                //else {
                //  clusters_coreBorder[row-ir][col+ic] = coreVal;
                //}
              }
              if (col - ic >= 0) {
                // Avoid center region with large values
                if (sqrt(pow(row-ir-centerR, 2) + pow(col-ic-centerC, 2))
                    > minCenterRad) {
                  // Fill vector tempPairs: contains pairs in minicluster
                  if (clusters_coreBorder[row-ir][col-ic] == initCoreVal) {
                    sum++;
                    Nchecks++;
                    clusters_coreBorder[row-ir][col-ic] = tempVal;
                    std::pair<uint, uint> coreInds(row-ir, col-ic);
                    tempPairs.push_back(coreInds);

                    if (col - ic < cMin) cMin = col - ic;
                    if (row - ir < rMin) rMin = row - ir;
                  }
                  else if (clusters_coreBorder[row-ir][col-ic] == 0) {
                    Nchecks++;
                  }
                }
                //else {
                //  clusters_coreBorder[row-ir][col-ic] = coreVal;
                //}
              }
            }
          }
        }
      }
    }

    rMean = 0;
    cMean = 0;
    count = 0;
    for (auto inds : tempPairs) {
      rMean += inds.first*energyMap[inds.first][inds.second];
      cMean += inds.second*energyMap[inds.first][inds.second];
      count += energyMap[inds.first][inds.second];
    }
    rMean /= count;
    cMean /= count;
    for (auto inds : tempPairs) {
      rVar += pow(rMean - inds.first, 2)*energyMap[inds.first][inds.second];
      cVar += pow(cMean - inds.second, 2)*energyMap[inds.first][inds.second];
    }
    rVar /= count;
    cVar /= count;

    ///  If cluster is large enough label points as core, otherwise border  ///
    int pixelMarker;

    sum > minClusterSize      ? pixelMarker = coreVal : pixelMarker = borderVal;
    sum < minPixelSize        ? pixelMarker = 0 : pixelMarker = pixelMarker;
    sum/Nchecks < minDensity  ? pixelMarker = 0 : pixelMarker = pixelMarker;
    if ((cVar/rVar < shapeVarRatio) && ((cMax - cMin)/(rMax - rMin) < shapeEdgeRatio)) {
      pixelMarker = 0;
    }
    for (uint ipp=0; ipp<tempPairs.size(); ipp++) {
      //if (sum > minClusterSize) {
      if (pixelMarker == coreVal) {
        clusterPairs.push_back(tempPairs[ipp]);
      }
      clusters_coreBorder[tempPairs[ipp].first][tempPairs[ipp].second] = pixelMarker;
    }
  }


  if (pltVerbose) {
    pltVerbose->printRC(clusters_coreBorder, "imgCoreClusterMap");
  }
  /////  Add border points to core points  /////
  for (uint ip=0; ip<clusterPairs.size(); ip++) {
    for (ir=-1*borderDistLimit; ir<=borderDistLimit; ir++) {
      row = (int)clusterPairs[ip].first + ir;
      if ((row >= Nrows) || (row < 0)) continue;

      for (ic=-1; ic<=1; ic++) {
        col = (int)clusterPairs[ip].second + ic;
        if ((col >= Ncols) || (col < 0)) continue;

        // Initial point is already a core point
        if ((ir == 0) && (ic == 0)) continue;

        // If it's a core point continue
        if (clusters_coreBorder[row][col]) continue;

        //cout<<"inpImg: "<<row<<"  "<<col<<"  "<<inpImg[row][col]<<endl;
        if (testVals.find(row) != testVals.end()) {
          if (testVals[row].find(col) != testVals[row].end()) {
            if (inpImg[row][col] == 0 && testVals[row][col]) cout<<"NOT ZERO"<<endl;
          }
        }

        sum = 0;
        for (irr=0; irr<borderRad; irr++) {
          irSqrd = irr*irr;
          // Looping columns
          for (icc=0; icc<borderRad; icc++) {
            if ((icc == 0) && (irr == 0)) continue;
            rad = std::sqrt(irSqrd + icc*icc);

            if (rad < borderRad) {
              //cout<<"row/col: "<<row<<"/"<<col<<"   "<<irr<<"/"<<icc<<"   "<<ratioInd<<endl;
              // Looking at rows above
              if (row + irr < Nrows) {
                //cout<<"11"<<endl;
                // Check if within radius
                if (col + icc < Ncols) {
                  //cout<<"summing: "<<row+irr<<"  "<<col+icc<<"  "<<inpImg[row+irr][col+icc]<<"  "<<inpImg[row+irr][col+icc]/rad<<endl;
                  sum += inpImg[row+irr][col+icc]/rad;//[ratioInd+irr*Ncols+icc]/rad;
                }
                if (col - icc >= 0) {
                  sum += inpImg[row+irr][col-icc]/rad;//[ratioInd+irr*Ncols-icc]/rad;
                }
                //cout<<"out11"<<endl;
              }
              // Looking at rows below
              if (row - irr >= 0) {
                // Check if within radius
                if (col + icc < Ncols) {
                  sum += inpImg[row-irr][col+icc]/rad;//[ratioInd-irr*Ncols+icc]/rad;
                }
                if (col - icc >= 0) {
                  sum += inpImg[row-irr][col-icc]/rad;//[ratioInd-irr*Ncols-icc]/rad;
                }
              }
            }
          }
        }

        energy = sum*inpImg[row][col]; //[ratioInd];

      //cout<<"energy: "<<energy<<"  "<<sum<<"  "<<inpImg[row][col]<<endl;
      //if (inpImg[row][col]) cout<<"NOT ZERO"<<endl;
        if (energy > borderValThresh) {
          clusters_coreBorder[row][col] = borderVal;

          std::pair<uint, uint> coreInds(row, col);
          clusterPairs.push_back(coreInds);
          removePairs.push_back(coreInds);
        }
      }
    }
  }

  if (pltVerbose) {
    pltVerbose->printRC(clusters_coreBorder, "imgCoreClusterBorderMap");
  }


  /////////////////////////////////////
  /////  Add padding to clusters  /////
  /////////////////////////////////////

  vector< std::pair<uint, uint> > padPairs;
  for (uint ip=0; ip<clusterPairs.size(); ip++) {
    for (ir=-1; ir<=1; ir++) {
      if (((int)clusterPairs[ip].first + ir >= Nrows)
        || ((int)clusterPairs[ip].first + ir < 0)) continue;

      for (ic=-1; ic<=1; ic++) {
        if (((int)clusterPairs[ip].second + ic >= Ncols)
          || ((int)clusterPairs[ip].second + ic < 0)) continue;

        // Initial point is already a core point
        if ((ir == 0) && (ic == 0)) continue;

        std::pair<uint, uint> padPair(row, col);
        padPairs.push_back(padPair);

        row = (int)clusterPairs[ip].first + ir;
        col = (int)clusterPairs[ip].second + ic;
        // If it's a cluster point continue
        if (clusters_coreBorder[row][col]) continue;

        ///  Initial pad pixel, now search for more around it  ///
        padPairs.clear();
        std::pair<uint, uint> padInds(row, col);
        padPairs.push_back(padInds);
        clusters_coreBorder[row][col] = padVal;
        rMin = rMax = row;
        cMin = cMax = col;
        for (uint ipp=0; ipp<padPairs.size(); ipp++) {
          for (irr=-1; irr<=1; irr++) {
            row = (int)padPairs[ipp].first + irr;
            if ((row >= Nrows) || (row < 0)) continue;

            for (icc=-1; icc<=1; icc++) {
              col = (int)padPairs[ipp].second + icc;
              if ((col >= Ncols) || (col < 0)) continue;

              if (clusters_coreBorder[row][col]) continue;

              // If there is a border point in padRad exit loop
              try {
                for (irrr=0; irrr<padRad; irrr++) {
                  irSqrd = irrr*irrr;
                  // Looping columns
                  for (iccc=0; iccc<padRad; iccc++) {
                    if ((iccc == 0) && (irrr == 0)) continue;
                    rad = std::sqrt(irSqrd + iccc*iccc);

                    if (rad < padRad) {
                      // Looking at rows above
                      if (row + irrr < Nrows) {
                        // Check if within radius
                        if (col + iccc < Ncols) {
                          if (clusters_coreBorder[row+irrr][col+iccc] > padVal) {
                            padInds.first = row;
                            padInds.second = col;
                            throw padInds;
                          }
                        }
                        if (col - iccc >= 0) {
                          if (clusters_coreBorder[row+irrr][col-iccc] > padVal) {
                            padInds.first = row;
                            padInds.second = col;
                            throw padInds;
                          }
                        }
                      }
                      // Looking at rows below
                      if (row - irrr >= 0) {
                        // Check if within radius
                        if (col + iccc < Ncols) {
                          if (clusters_coreBorder[row-irrr][col+iccc] > padVal) {
                            padInds.first = row;
                            padInds.second = col;
                            throw padInds;
                          }
                        }
                        if (col - iccc >= 0) {
                          if (clusters_coreBorder[row-irrr][col-iccc] > padVal) {
                            padInds.first = row;
                            padInds.second = col;
                            throw padInds;
                          }
                        }
                      }
                    }
                  }
                }
              }
              catch (std::pair<uint, uint> inds) {
                if (clusters_coreBorder[inds.first][inds.second]) cout<<"PADERR"<<endl;

                clusters_coreBorder[inds.first][inds.second] = padVal;
                padPairs.push_back(inds);
                if (inds.first < rMin) rMin = inds.first;
                if (inds.first > rMax) rMax = inds.first;
                if (inds.second < cMin) cMin = inds.second;
                if (inds.second > cMax) cMax = inds.second;
              }
            }
          }
        }

        // Extend padding based on shape
        int rDiff = rMax - rMin;
        int cDiff = cMax - cMin;
        for (int ir=(int)(rMin+rMinScale*rDiff); ir<(int)(rMin+rMaxScale*rDiff); ir++) {
          if ((ir >= Nrows) || (ir < 0)) continue;
          for (int ic=(int)(cMin+cMinScale*cDiff); ic<(int)(cMin+cMaxScale*cDiff); ic++) {
            if ((ic >= Ncols) || (ic < 0)) continue;
            if (clusters_coreBorder[ir][ic] == 0) {
              padInds.first = ir; 
              padInds.second = ic;
              padPairs.push_back(padInds);
              clusters_coreBorder[ir][ic] = padVal;
            }
          }
        }

        removePairs.insert(removePairs.end(),
        padPairs.begin(), padPairs.end());
      }
    }
  }

  if (pltVerbose) {
    pltVerbose->printRC(clusters_coreBorder, "imgCoreClusterBorderPadMap");
  }
  return removePairs;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
cv::Mat imgProc::buildTeirGrayImgT(cv::Mat &inp, long int modVal, long int imScale) {

  cv::Mat out(inp.rows, inp.cols, CV_8U);

  if (modVal > 256) {cerr<<"ERROR: modVal must be <=256 since gray scale is 8 bit!!!"<<endl; exit(0);}
  int Nrows = inp.rows;
  int Ncols = inp.cols;

  if (inp.isContinuous()) {
    Ncols *= Nrows;
    Nrows = 1;
  }

  T* val;

  uint8_t* valG;
  for (int ir=0; ir<Nrows; ir++)  {
    val = inp.ptr<T>(ir);
    valG = out.ptr<uint8_t>(ir);
    for (int ic=0; ic<Ncols; ic++) valG[ic] = (uint8_t) fmod(((long int)val[ic])/imScale, 180.0);
  } 

  return out;
}


cv::Mat imgProc::buildTeirGrayImg(cv::Mat &inp, long int modVal, long int imScale) {

  switch (inp.depth()) {
    case 0:
        return buildTeirGrayImgT<uint8_t>(inp, modVal, imScale);
        break;
    case 1:
        return buildTeirGrayImgT<int8_t>(inp, modVal, imScale);
        break;
    case 2:
        return buildTeirGrayImgT<uint16_t>(inp, modVal, imScale);
        break;
    case 3:
        return buildTeirGrayImgT<int16_t>(inp, modVal, imScale);

        break;
    case 4:
        return buildTeirGrayImgT<float>(inp, modVal, imScale);
        break;
    case 5:
        return buildTeirGrayImgT<double>(inp, modVal, imScale);
        break;
    default:
        cerr<<"ERROR: Do not recognize bit depth enum "<<inp.depth()<<" in function buildTeirGrayImg!!!"<<endl;
        exit(0);
  }
}



///////////////////////////////////////////////////////////////////////////////////////////////////////
vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, cv::Mat &cannyEdges, vector<cv::Vec4i> &hierarchy, int Ncontours) {

  vector< vector<cv::Point> > allContours;
  vector< vector<cv::Point> > centContours;

  if (!(mat.depth() == 0 || mat.depth() == 1)) {
    cerr<<"ERROR: findCenteredContours can only be used for 8 bit matrices!!!"<<endl;
    return centContours;
  }

  Canny(mat, cannyEdges, cannyLow, cannyHigh, cannyKernel);

  findContours( cannyEdges, allContours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0) );

  bool passx, passy;
  int xdist, ydist;
  vector<cv::Point> centers;
  for(unsigned int ic=0; ic<allContours.size(); ic++ ) {
    passx=passy=false;
    xdist = (int)(center.x-allContours[ic][0].x);
    ydist = (int)(center.y-allContours[ic][0].y);
    if (outRad && (fabs(xdist)>outRad && fabs(ydist)>outRad)) continue;
    if (inRad && (fabs(xdist)<inRad && fabs(ydist)<inRad)) continue;
    if ((abs(allContours[ic][0].x-allContours[ic][allContours[ic].size()-1].x)>3) || (abs(allContours[ic][0].y-allContours[ic][allContours[ic].size()-1].y)>3)) continue;
    for (unsigned int ip=0; ip<allContours[ic].size(); ip++) {
      if (passx && passy) break;
      if (!passx && (abs((center.x-allContours[ic][ip].x)+xdist)<20)) passx=true;
      if (!passy && (abs((center.y-allContours[ic][ip].y)+ydist)<20)) passy=true;
    }
    
    if (passx && passy) centContours.push_back(allContours[ic]);
  }

  if (int(centContours.size())<Ncontours) cerr<<"WARNING: Only find "<<centContours.size()<<" contours instead of required "<<Ncontours<<", increase scaling for imgProc::buildTeirGrayImg!!"<<endl;
  return centContours;
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad) {

  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0);
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad) {

  cv::Point center;
  center.x = mat.cols/2;
  center.y = mat.rows/2;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0);
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, int Ncontours) {

  cv::Point center;
  center.x = mat.cols/2;
  center.y = mat.rows/2;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, Ncontours);
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh, int inRad, int outRad) {

  cv::Point center;
  center.x = mat.cols/2;
  center.y = mat.rows/2;
  int cannyKernel = 3;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0);
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int inRad, int outRad) {

  int cannyKernel = 3;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh) {

  int cannyKernel = 3;
  int inRad = 0;
  int outRad = 0;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
}


vector<vector<cv::Point> > imgProc::findCenteredContours(cv::Mat &mat, int cannyLow, int cannyHigh) {

  cv::Point center;
  center.x = mat.cols/2;
  center.y = mat.rows/2;
  int cannyKernel = 3;
  int inRad = 0;
  int outRad = 0;
  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  return findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
}


double imgProc::centerOddLeg(vector< vector<double> >* img, double centerR, double centerC, int Nshells, double minR, double shellW, bool nanBins) {
  ////// This function takes divides a 2d shell into 2*Nsections and sums up all the contributions 
  //////	in each section. It then averages the square fractional difference between opposite 
  //////	sections and returns this value as the chi^2.
  // Nsections: number of sections per quadrant
  // Nshells: how many shells to consider
  // shellW: bin width of each shell
  // minR: bin radii to start at 
 
  double angle;
  double goodness=0;
  double prad, sum1L0, sum2L0, sum1L1, sum2L1;
  double rads[2];
  int col, row, row0;
  bool pass1, pass2, pass3, pass4;

  // Summing sections on opposite sides
  for (int ish=0; ish<(int)Nshells; ish++) {
    rads[0] = minR + shellW*ish;
    rads[1] = rads[0] + shellW;
    sum1L0 = sum2L0 = 0;
    sum1L1 = sum2L1 = 0;
    col = centerC+ (int)rads[1]; 
    row = row0 = centerR+1;
    if (!nanBins) {
      while (col >= centerC) {
        row = row0;
        prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        while (prad <= rads[0]) {
          row++;
          row0 = row;
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        }
        while (prad <= rads[1]) {
          angle = atan(((double)(row-centerR))/((double)(col-centerC)));
          // Legendre 1
          sum1L0 += (*img)[row][col];
          sum1L0 += (*img)[row][2*centerC-col];
          sum2L0 += (*img)[2*centerR-row][2*centerC-col];
          sum2L0 += (*img)[2*centerR-row][col];
          // Legendre 1
          sum1L1 += (*img)[row][col]*boost::math::legendre_p(1, cos(angle));
          sum1L1 += (*img)[row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI/2.0));
          sum2L1 += (*img)[2*centerR-row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI));
          sum2L1 += (*img)[2*centerR-row][col]*boost::math::legendre_p(1, cos(angle+PI*0.75));
          row++;
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        }
        col--;
      }


/*
      row = centerR;
      col = centerC + rad[0];
      while (col <= centerC+rads[1]) 
        // Legendre 1
        sum1L1 += (*img)[row][col]*boost::math::legendre_p(1, 1));
        sum1L1 += (*img)[row][2*centerC-col]*boost::math::legendre_p(1, cos(PI/2.0));
        sum2L1 += (*img)[2*centerR-row][2*centerC-col]*boost::math::legendre_p(1, cos(PI));
        sum2L1 += (*img)[2*centerR-row][col]*boost::math::legendre_p(1, cos(PI*0.75);
*/


    }
    else {
      while (col >= centerC) {
        row = row0;
        prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        while (prad <= rads[0]) {
          row++;
          row0 = row;
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        }
        while (prad <= rads[1]) {
          angle = atan(((double)(row-centerR))/((double)(col-centerC)));
 
          pass1 = ((*img)[row][col] != NANVAL);
          pass2 = ((*img)[row][2*centerC-col] != NANVAL);
          pass3 = ((*img)[2*centerR-row][2*centerC-col] != NANVAL);
          pass4 = ((*img)[2*centerR-row][col] != NANVAL);
          if (pass1 && pass3) {
            // Legendre 0
            sum1L0 += (*img)[row][col];
            sum2L0 += (*img)[2*centerR-row][2*centerC-col];
            // Legendre 1
            sum1L1 += (*img)[row][col]*boost::math::legendre_p(1, cos(angle));
            sum2L1 += (*img)[2*centerR-row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI));
          }
          else if (pass1 && !pass3) {
            // Legendre 0
            sum1L0 += (*img)[row][col];
            sum2L0 += (*img)[row][col];
            // Legendre 1
            sum1L1 += (*img)[row][col]*boost::math::legendre_p(1, cos(angle));
            sum2L1 += (*img)[row][col]*boost::math::legendre_p(1, cos(angle+PI));
          }
          else if (!pass1 && pass3) {
            // Legendre 0
            sum1L0 += (*img)[2*centerR-row][2*centerC-col];
            sum2L0 += (*img)[2*centerR-row][2*centerC-col];
            // Legendre 1
            sum1L1 += (*img)[2*centerR-row][2*centerC-col]*boost::math::legendre_p(1, cos(angle));
            sum2L1 += (*img)[2*centerR-row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI));
          }
 
          if (pass2 && pass4) {
            // Legendre 0
            sum1L0 += (*img)[row][2*centerC-col];
            sum2L0 += (*img)[2*centerR-row][col];
            // Legendre 1
            sum1L1 += (*img)[row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI/2.0));
            sum2L1 += (*img)[2*centerR-row][col]*boost::math::legendre_p(1, cos(angle+PI*0.75));
          }
          else if (pass2 && !pass4) {
            // Legendre 0
            sum1L0 += (*img)[row][2*centerC-col];
            sum2L0 += (*img)[row][2*centerC-col];
            // Legendre 1
            sum1L1 += (*img)[row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI*2.0));
            sum2L1 += (*img)[row][2*centerC-col]*boost::math::legendre_p(1, cos(angle+PI*0.75));
          }
          else if (!pass2 && pass4) {
            // Legendre 0
            sum1L0 += (*img)[2*centerR-row][col];
            sum2L0 += (*img)[2*centerR-row][col];
            // Legendre 1
            sum1L1 += (*img)[2*centerR-row][col]*boost::math::legendre_p(1, cos(angle+PI/2.0));
            sum2L1 += (*img)[2*centerR-row][col]*boost::math::legendre_p(1, cos(angle+PI*0.75));
          }

          row++;
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
        }
        col--;
      }
    }
  
    // function to measure goodness of center
    //goodness += sqrt(pow(sum1L1,2)+pow(sum2L1,2))*fabs(sum1L0-sum2L0)/(pow(sum1L0,2)+pow(sum2L0,2));
    goodness += pow(fabs(sum1L1) + fabs(sum2L1),3)*fabs(sum1L0-sum2L0)/(pow(sum1L0,2)+pow(sum2L0,2));
  }

//cout<<centerR<<"/"<<centerC<<"   "<<(pow(sum1L1,2)+pow(sum2L1,2))<<"/"<<(pow(sum1L0,2)+pow(sum2L0,2))<<"    "<<goodness<<endl;
  return goodness;
}



double imgProc::centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW, bool nanBins) {
  ////// This function takes divides a 2d shell into 2*Nsections and sums up all the contributions 
  //////	in each section. It then averages the square fractional difference between opposite 
  //////	sections and returns this value as the chi^2.
  // Nsections: number of sections per quadrant
  // Nshells: how many shells to consider
  // shellW: bin width of each shell
  // minR: bin radii to start at 
 
  double angle[2];
  double angIter = (PI/2.0)/Nsections;
  double goodness=0;
  double entries=0;
  double weight;
  double prad, pang, sum1, sum2, sum3, sum4;
  double pol[4];
  double rads[2];
  int maxr, minr, maxc, minc, col, row;

  for (int ish=0; ish<(int)Nshells; ish++) {
    rads[0] = minR + shellW*ish;
    rads[1] = rads[0] + shellW;
    angle[0] = 0;
    for (int isct=0; isct<(int)Nsections; isct++) {
      angle[0] = isct*angIter;  
      angle[1] = angle[0] + angIter;
      pol[0] = sin(angle[0]); pol[1] = sin(angle[1]);
      pol[2] = cos(angle[0]); pol[3] = cos(angle[1]);

      // Making opposing sections to search for bins within opposite sections
      maxr = rads[1]*pol[1] + centerR + 1;
      minr = rads[0]*pol[0] + centerR - 1;
      maxc = rads[1]*pol[2] + centerC + 1;
      minc = rads[0]*pol[3] + centerC - 1;

      // Summing sections on opposite sides
      sum1 = sum2 = sum3 = sum4 = 0;
      row = minr; col = minc;
      if (!nanBins) {
        while (col <= maxc) {
          while (row <= maxr) {
            prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
            pang = atan(((double)(row-centerR))/((double)(col-centerC)));
            if ((prad>=rads[0] && prad<=rads[1]) && (angle[0]<=pang && pang<=angle[1])) {
              // Quadrant 1
              sum1 += (*img)[row][col];
              // Quadrant 2
              sum2 += (*img)[row][2*centerC-col];
              // Quadrant 3
              sum3 += (*img)[2*centerR-row][2*centerC-col];
              // Quadrant 4
              sum4 += (*img)[2*centerR-row][col];
            }
            row++;
          }
          col++;
          row = minr;
        }
      }
      else {
        while (col <= maxc) {
          while (row <= maxr) {
            prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
            pang = atan(((double)(row-centerR))/((double)(col-centerC)));
            if ((prad>=rads[0] && prad<=rads[1]) && (angle[0]<=pang && pang<=angle[1])) {

	      if (((*img)[row][col] != NANVAL) && ((*img)[2*centerR-row][2*centerC-col] != NANVAL)) {
                // Quadrant 1
                sum1 += (*img)[row][col];
                // Quadrant 3
                sum3 += (*img)[2*centerR-row][2*centerC-col];
	      }
	      if (((*img)[row][2*centerC-col] != NANVAL) && ((*img)[2*centerR-row][col] != NANVAL)) {
                // Quadrant 2
                sum2 += (*img)[row][2*centerC-col];
                // Quadrant 4
                sum4 += (*img)[2*centerR-row][col];
	      }
            }
            row++;
          }
          col++;
          row = minr;
        }
      }
 
      // function to measure goodness of center
      //goodness += pow((fabs(sum1-sum2)/(sum1+sum2)),2);
      (pol[0]+pol[1])/2 <= (pol[2]+pol[3])/2 ? weight = pow((pol[2]+pol[3])/2,2) : weight = pow((pol[0]+pol[1])/2,2);
      goodness += fabs(sum1-sum3)*weight;
      goodness += fabs(sum2-sum4)*weight;
      entries += 1.0;
    }
  }

  //cout<<"RESULT: "<<goodness<<endl;
  //return goodness/entries;
  return goodness;
}


double imgProc::centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW) {

  
  return centerSymXsqr(img, centerR, centerC, Nsections, Nshells, minR, shellW, false); 
}
/*
double imgProc::centerSymXsqr(vector< vector<double> >* img, double centerR, double centerC, double Nsections, double Nshells, double minR, double shellW, double angle) {
  ////// This function takes divides a 2d shell into 2*Nsections and sums up all the contributions 
  //////	in each section. It then averages the square fractional difference between opposite 
  //////	sections and returns this value as the chi^2.
  // Nsections: number of sections per each hemishell
  // Nshells: how many shells to consider
  // shellW: bin width of each shell
  // minR: bin radii to start at 
  // angle: where to define the angle 0   
 
  double angIter = PI/Nsections;
  double goodness=0;
  double entries=0;
  double ang, rad, prad, pang, sum1, sum2, A0, A1;
  double pol[4];
  double rads[2];
  int maxr1, maxr2, minr1, minr2, maxc1, maxc2, minc1, minc2, col, row;

  for (int ish=0; ish<(int)Nshells; ish++) {
    rad = minR + shellW*ish;
    rads[0] = rad;    rads[1] = rad + shellW;
    for (int isct=0; isct<(int)Nsections; isct++) {
      ang = angle + isct*angIter;
      pol[0] = sin(ang); pol[1] = sin(ang+angIter);
      pol[2] = cos(ang); pol[3] = cos(ang+angIter);

      // Making opposing sections to search for bins within opposite sections
      maxr1=minr1=maxr2=minr2=rads[0]*pol[0];
      maxc1=minc1=maxc2=minc2=rads[0]*pol[2];
      for (int ir=0; ir<2; ir++) {
        for (int ip=0; ip<2; ip++) {
          maxr1 = maxr1>(rads[ir]*pol[ip]) ? maxr1 : (rads[ir]*pol[ip]);
          minr1 = minr1<(rads[ir]*pol[ip]) ? minr1 : (rads[ir]*pol[ip]);
          maxr2 = maxr2>(-rads[ir]*pol[ip]) ? maxr2 : (-rads[ir]*pol[ip]);
          minr2 = minr2<(-rads[ir]*pol[ip]) ? minr2 : (-rads[ir]*pol[ip]);
          maxc1 = maxc1>(rads[ir]*pol[ip+2]) ? maxc1 : (rads[ir]*pol[ip+2]);
          minc1 = minc1<(rads[ir]*pol[ip+2]) ? minc1 : (rads[ir]*pol[ip+2]);
          maxc2 = maxc2>(-rads[ir]*pol[ip+2]) ? maxc2 : (-rads[ir]*pol[ip+2]);
          minc2 = minc2<(-rads[ir]*pol[ip+2]) ? minc2 : (-rads[ir]*pol[ip+2]);
        }
      }
      maxr1+=centerR+1; minr1+=centerR-1;
      maxr2+=centerR+1; minr2+=centerR-1;
      maxc1+=centerC+1; minc1+=centerC-1;
      maxc2+=centerC+1; minc2+=centerC-1;

      // Defining boundary angles since atan = [-PI/2,PI/2]
      A0 = ang>(PI/2.0) ? ang-PI : ang;
      A1 = (ang+angIter)>(PI/2.0) ? (ang+angIter)-PI : (ang+angIter);

      // Summing sections on opposite sides
      sum1 = 0;
      row = minr1; col = minc1;
      while (col <= maxc1) {
        while (row <= maxr1) {
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
          pang = atan(((double)(row-centerR))/((double)(col-centerC)));
          if (fabs(A0-A1)<1.2*angIter) {
            if ((prad>=rad && prad<=(rad+shellW)) && (pang>=A0 && pang<=A1)) sum1 += (*img)[row][col];
          }
          else {
            if ((prad>=rad && prad<=(rad+shellW)) && (pang>=A0 || pang<=A1)) sum1 += (*img)[row][col];
          }
          row++;
        }
        col++;
        row = minr1;
      }

      sum2 = 0;
      row = minr2; col = minc2;
      while (col <= maxc2) {
        while (row <= maxr2) {
          prad = sqrt(pow((double)(row-centerR),2) + pow((double)(col-centerC),2));
          pang = atan(((double)(row-centerR))/((double)(col-centerC)));
          if (fabs(A0-A1)<1.2*angIter) {
            if ((prad>=rad && prad<=(rad+shellW)) && (pang>=A0 && pang<=A1)) sum2 += (*img)[row][col];
          }
          else {
            if ((prad>=rad && prad<=(rad+shellW)) && (pang>=A0 || pang<=A1)) sum2 += (*img)[row][col];
          }
          row++;
        }
        col++;
        row = minr2;
      }

      // function to measure goodness of center
      //goodness += pow((fabs(sum1-sum2)/(sum1+sum2)),2);
      goodness += fabs(sum1-sum2);
      entries += 1.0;
    }
  }

  //cout<<"RESULT: "<<goodness<<endl;
  //return goodness/entries;
  return goodness;
}
*/




//int imgProc::findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector<vector<cv::Point> > &contours, cv::Mat &cannyEdges, vector<cv::Vec4i> &hierarchy) 
  //return findCenterTC(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, contours, cannyEdges, hierarchy) 
///////////////////////////////////////////////////////////////////////////////////////////////////////
cv::Point imgProc::findContourCenters(vector< vector<cv::Point> > contours, vector<cv::Point> &centers) {

  int sumx, sumy;
  cv::Point avgCent;

  if (contours.size() == 0) {
    cerr<<"WARNING: Cannot find the center of 0 contours, now exiting with false value!!!"<<endl;
    return avgCent;
  }

  for(unsigned int ic=0; ic<contours.size(); ic++ ) {
    sumx=sumy=0;
    for (unsigned int ip=0; ip<contours[ic].size(); ip++) {
      sumx += contours[ic][ip].x;
      sumy += contours[ic][ip].y;
    }
    cv::Point p(sumx/contours[ic].size(), sumy/contours[ic].size());
    avgCent.x += p.x;
    avgCent.y += p.y;
    centers.push_back(p);
    //cv::Scalar color = cv::Scalar( 255,255, 255 );
    //drawContours( mat, contours, ic, color, 2, 8, hierarchy, 0, cv::Point() );
  }

  avgCent.x /= contours.size();
  avgCent.y /= contours.size();

  return avgCent;
}


cv::Point imgProc::findContourCenters(vector< vector<cv::Point> > contours) {

  vector<cv::Point> centers;

  return findContourCenters(contours, centers);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
cv::Point imgProc::findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Vec4i> &hierarchy, vector<cv::Point> &centers) {

  cv::Mat cannyEdges;

  contours = findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
  
  return findContourCenters(contours, centers);
}

 
cv::Point imgProc::findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Point> &centers) {

  cv::Mat cannyEdges;
  vector<cv::Vec4i> hierarchy;

  contours = findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
  
  return findContourCenters(contours, centers);
}

 
cv::Point imgProc::findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad, vector< vector<cv::Point> > &contours, vector<cv::Point> &centers, int Ncontours) {

  contours = findCenteredContours(mat, cannyLow, cannyHigh, cannyKernel, inRad, outRad, Ncontours); 

  return findContourCenters(contours, centers);
}

 
cv::Point imgProc::findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad) {

  cv::Mat cannyEdges;
  vector< vector<cv::Point> > contours;
  vector<cv::Vec4i> hierarchy;

  contours = findCenteredContours(mat, center, cannyLow, cannyHigh, cannyKernel, inRad, outRad, cannyEdges, hierarchy, 0); 
  
  return findContourCenters(contours);
}

  
cv::Point imgProc::findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int cannyKernel, int inRad, int outRad) {

  vector< vector<cv::Point> > contours;

  contours = findCenteredContours(mat, cannyLow, cannyHigh, cannyKernel, inRad, outRad);
  
  return findContourCenters(contours);
}
 

cv::Point imgProc::findCenterTC(cv::Mat &mat, cv::Point &center, int cannyLow, int cannyHigh, int inRad, int outRad) {

  vector< vector<cv::Point> > contours;

  contours = findCenteredContours(mat, center, cannyLow, cannyHigh, inRad, outRad);
  
  return findContourCenters(contours);
}
 

cv::Point imgProc::findCenterTC(cv::Mat &mat, int cannyLow, int cannyHigh, int inRad, int outRad) {

  cv::Mat cannyEdges;
  vector< vector<cv::Point> > contours;
  vector<cv::Vec4i> hierarchy;

  contours = findCenteredContours(mat, cannyLow, cannyHigh, inRad, outRad);
  
  return findContourCenters(contours);
}
 



void imgProc::saveAscii(cv::Mat mat, string fileName) {

  ofstream output;
  output.open((fileName+".dat").c_str());

  int Nrows = mat.rows;
  int Ncols = mat.cols;

  if (mat.isContinuous()) {
    Ncols *= Nrows;
    Nrows = 1;
  }

  uint16_t *val;
  for (int ir=0; ir<Nrows; ir++)  {
    val = mat.ptr<uint16_t>(ir);
    for (int ic=0; ic<Ncols; ic++) {
      output << val[ic] << "\t";
      if (ic%mat.cols == 0) output << "\n";
    }
    if ((ir+1) != mat.rows) output << "\n";
  }

  output.close();
}
 

void imgProc::saveAscii(TH2F *hist, string fileName) {

  ofstream output;
  output.open((fileName+".dat").c_str());

  int Nrows = hist->GetNbinsY();
  int Ncols = hist->GetNbinsX();

  for (int ir=0; ir<Nrows; ir++)  {
    for (int ic=0; ic<Ncols; ic++) output << hist->GetBinContent(ic,ir) << "\t";
    if ((ir+1) != Nrows) output << "\n";
  }

  output.close();
}
    


 
///////////////////////////////////////////////////////////////////////////////////////////////////////
double imgProc::cubicInterpolate (double p[4], double x) {
  return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}


double imgProc::bicubicInterpolate (double p[4][4], double row, double col) {

  if (col<0 || row<0 || row>1 || col>1) {
    cerr<<"ERROR: BicubicInterpolation requires (x,y) bin values to be between 0-1, instead given ("<<col<<","<<row<<")!!!"<<endl;
    exit(0);
  }

  double arr[4];
  arr[0] = cubicInterpolate(p[0], col);
  arr[1] = cubicInterpolate(p[1], col);
  arr[2] = cubicInterpolate(p[2], col);
  arr[3] = cubicInterpolate(p[3], col);
  return cubicInterpolate(arr, row);
}


double imgProc::interpolate(cv::Mat mat, double row, double col) {

  switch (mat.depth()) {
    case 0:
        return interpolateT<uint8_t>(mat, row, col);
        break;
    case 1:
        return interpolateT<int8_t>(mat, row, col);
        break;
    case 2:
        return interpolateT<uint16_t>(mat, row, col);
        break;
    case 3:
        return interpolateT<int16_t>(mat, row, col);
        break;
    case 4:
        return interpolateT<float>(mat, row, col);
        break;
    case 5:
        return interpolateT<double>(mat, row, col);
        break;
    default:
        cerr<<"ERROR: Do not recognize bit depth enum "<<mat.depth()<<" in function interpolate!!!"<<endl;
        exit(0);
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////

vector< vector<double> > imgProc::polarBinning(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins) {
  return polarBinning(mat, center, Nslices, NmatBins, NhistBins, false, 0);
}

vector< vector<double> > imgProc::polarBinning(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan) {

  switch (mat.depth()) {
    case 0:
        return polarBinningT<uint8_t>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    case 1:
        return polarBinningT<int8_t>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    case 2:
        return polarBinningT<uint16_t>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    case 3:
        return polarBinningT<int16_t>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    case 4:
        return polarBinningT<float>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    case 5:
        return polarBinningT<double>(mat, center, Nslices, NmatBins, NhistBins, doNAN, nan);
        break;
    default:
        cerr<<"ERROR: Do not recognize bit depth enum "<<mat.depth()<<" in function polarBinning!!!"<<endl;
        exit(0);
  }
}


template<typename T>
vector< vector<double> > imgProc::polarBinningT(cv::Mat mat, cv::Point center, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan) {
  ////// Plots the radial projection as a function of angle. Each angle has 4 hists in cross shape. 
  //////	hists 0,1,2,3 start at angle 0,PI/2,PI,3*PI/2.
  //  Nslices: how many angular steps from 0 to PI/2, total of 4*Nslices hists
  //  NmatBins: number of bins to go out to from center of mat
  //  NhistBins: number of bins in the projection histograms
  //  nan: any bin with a contribution of a cartesian bin with this value will be nan
  //  combined: plot all slices on a 2D histo in order as a functin of angle
  //  normalization of binning is done where one mat bin has area 1
  //  OUTPUT: vector of the 4 slices in cross shape for each angle

  cerr<<"WARNING: This version of polarBinning has not been thoroughly tested but changed accordering the tested version!"<<endl;

  if (!mat.isContinuous()) {
    cerr<<"ERROR: Radial Plots only works with continuous matrices, can easily be modified!!!"<<endl;
    exit(0);
  }
  if (NhistBins > NmatBins) {
    cerr<<"ERROR: Cannot make radial plots with NhistBins > NmatBins, cannot pull information out of nowhere"<<endl;
    exit(0);
  }

  vector< vector<double> > radHists;
  bool combineBins = true;

  int NradSamp = 20;
  int NangSamp = 10;
  int NtotAngSamp = NangSamp/Nslices;
  ostringstream strStm;

  double angle = 0;
  double iterA;
  double iterR = ((double)NmatBins)/((double)(NhistBins*NradSamp));
  double iterS = (PI/2.0)/((double)Nslices);
  double binRat = NmatBins/NhistBins;

  radHists.resize((int)NhistBins);
  for (uint irad=0; irad<radHists.size(); irad++) radHists[irad].resize(4*Nslices);
  T* matVals = mat.ptr<T>(0);
  int Nrows = mat.rows;

  double sampCosA, sampSinA, sinA, cosA, sampRad, sampAng, sampArea;
  double binVal1, binVal2, binVal3, binVal4;
  cv::Point samp1,samp2,samp3,samp4;
  CvPoint2D64f point1,point2,point3,point4;

  for (int isl=0; isl<Nslices; isl++) {
    angle = iterS*isl;
    strStm << angle;

    if (combineBins) {  	// Fills each bin by sampling bins within the polar bin's area
      for (int ib=0; ib<(int)NhistBins; ib++) {
        binVal1=binVal2=binVal3=binVal4=0;
        for (int irad=0; irad<NradSamp; irad++) {
          sampRad = ib*binRat + irad*iterR + iterR/2.0;
	  NtotAngSamp = ((int)sampRad+1)*NangSamp*iterS + 1;
	  iterA = iterS/((double)NtotAngSamp);
          sampArea = (iterA/2)*(pow(sampRad+iterR-iterR/2, 2) - pow(sampRad-iterR/2, 2));
          for (int iang=0; iang<NtotAngSamp; iang++) {
            sampAng = angle + iang*iterA + iterA/2.0; 
            sampCosA = cos(sampAng);	sampSinA = sin(sampAng);
            samp1.x = (int)(center.x + sampRad*sampCosA);	samp1.y = (int)(center.y - sampRad*sampSinA);
            samp2.x = (int)(center.x - sampRad*sampSinA);	samp2.y = (int)(center.y - sampRad*sampCosA);
            samp3.x = (int)(center.x - sampRad*sampCosA);	samp3.y = (int)(center.y + sampRad*sampSinA);
            samp4.x = (int)(center.x + sampRad*sampSinA);	samp4.y = (int)(center.y + sampRad*sampCosA);
 
            binVal1 += sampArea*matVals[samp1.y*Nrows + samp1.x];
            binVal2 += sampArea*matVals[samp2.y*Nrows + samp2.x];
            binVal3 += sampArea*matVals[samp3.y*Nrows + samp3.x];
            binVal4 += sampArea*matVals[samp4.y*Nrows + samp4.x];
          }
        }

	radHists[ib][isl] = binVal1;
	radHists[ib][isl+Nslices] = binVal2;
	radHists[ib][isl+2*Nslices] = binVal3;
	radHists[ib][isl+3*Nslices] = binVal4;
      }
    }
    else {	// Takes an interpolation at the center of each bin, treats each bin as the same size
      sinA = sin(angle);
      cosA = cos(angle);
      point1.x=point2.x=point3.x=point4.x=center.x;
      point1.y=point2.y=point3.y=point4.y=center.y;
      for (int ib=0; ib<NhistBins; ib++) {
	radHists[ib][isl] = binVal1;
	radHists[ib][isl+Nslices] = binVal2;
	radHists[ib][isl+2*Nslices] = binVal3;
	radHists[ib][isl+3*Nslices] = binVal4;

        //      COLUMNS			 ROWS
        point1.x += iterR*cosA; 	point1.y += -iterR*sinA;
        point2.x += -iterR*sinA; 	point2.y += -iterR*cosA;
        point3.x += -iterR*cosA; 	point3.y += iterR*sinA;
        point4.x += iterR*sinA; 	point4.y += iterR*cosA;
      }
    }  
  
    strStm.str("");
    strStm.clear();
  }

  return radHists;
}


vector< vector<double> > imgProc::polarBinning(vector< vector<double> > mat, double centerR, double centerC, int Nslices, double NmatBins, double NhistBins) {
  return polarBinning(mat, centerR, centerC, 
      Nslices, NmatBins, NhistBins, false, 0);
}


vector< vector<double> > imgProc::polarBinning(vector< vector<double> > mat, double centerR, double centerC, int Nslices, double NmatBins, double NhistBins, bool doNAN, double nan) {
  ////// Plots the radial projection as a function of angle. Each angle has 4 hists in cross shape. 
  //////	hists 0,1,2,3 start at angle 0,PI/2,PI,3*PI/2.
  //  Nslices: how many angular steps from 0 to PI/2, total of 4*Nslices hists
  //  NmatBins: number of bins to go out to from center of mat
  //  NhistBins: number of bins in the projection histograms
  //  nan: any polar bin with a contribution from a nan cartesian bin will become nan
  //  combined: plot all slices on a 2D histo in order as a functin of angle
  //  normalization of binning is done where one mat bin has area 1
  //  OUTPUT: vector of the 4 slices in cross shape for each angle

  if (NhistBins > NmatBins) {
    cerr<<"ERROR: Cannot make radial plots with NhistBins > NmatBins, cannot pull information out of nowhere"<<endl;
    exit(0);
  }


  vector< vector<double> > radHists;

  int NradSamp = 20;
  int NangSamp = 20;
  int NtotAngSamp;
  ostringstream strStm;

  double angle = 0;
  double iterA;
  double iterR = NmatBins/(NhistBins*NradSamp);
  double iterS = (PI/2.0)/((double)Nslices);
  double binRat = NmatBins/NhistBins;

  radHists.resize((int)NhistBins);
  for (uint irad=0; irad<radHists.size(); irad++) {
    radHists[irad].resize(4*Nslices);
  }

  double sampCosA, sampSinA, sampRad, sampAng, sampArea;
  double binVal1, binVal2, binVal3, binVal4;
  double TbinVal1, TbinVal2, TbinVal3, TbinVal4;
  int samp1R, samp2R, samp3R, samp4R;
  int samp1C, samp2C, samp3C, samp4C;

  // Used for nan bins
  bool binVal1_nan;
  bool binVal2_nan;
  bool binVal3_nan;
  bool binVal4_nan;


  /////  Looping over angular bins  /////
  for (int isl=0; isl<Nslices; isl++) {
    angle = iterS*isl;
    strStm << angle;

      /////  Looping over radial bins  /////
      for (int ib=0; ib<(int)NhistBins; ib++) {
        binVal1=binVal2=binVal3=binVal4=0;
//cout<<"angle isl Bin: "<<angle<<"  "<<isl<<"  "<<ib<<endl;
        if (doNAN) {
          binVal1_nan = false;
          binVal2_nan = false;
          binVal3_nan = false;
          binVal4_nan = false;
        }
        ///  Looping over small portions of rad bin  ///
        for (int irad=0; irad<NradSamp; irad++) {
          sampRad = ib*binRat + irad*iterR + iterR/2.0;
	  NtotAngSamp = ((int)sampRad+1)*NangSamp*iterS + 1;
	  iterA = iterS/((double)NtotAngSamp);
          sampArea = (iterA/2)*(pow(sampRad+iterR-iterR/2, 2) - pow(sampRad-iterR/2, 2));
          TbinVal1=TbinVal2=TbinVal3=TbinVal4=0;
//cout<<"  radSamp: "<<sampRad<<endl;
          ///  Looping over small portions of angular bins  ///
          for (int iang=0; iang<NtotAngSamp; iang++) {
            sampAng = angle + iang*iterA + iterA/2.0; 
            sampCosA = cos(sampAng);	sampSinA = sin(sampAng);
            samp1C = (int)(centerC + sampRad*sampCosA);	samp1R = (int)(centerR + sampRad*sampSinA);
            samp2C = (int)(centerC - sampRad*sampSinA);	samp2R = (int)(centerR + sampRad*sampCosA);
            samp3C = (int)(centerC - sampRad*sampCosA);	samp3R = (int)(centerR - sampRad*sampSinA);
            samp4C = (int)(centerC + sampRad*sampSinA);	samp4R = (int)(centerR - sampRad*sampCosA);
 
            TbinVal1 += mat[samp1R][samp1C];
            TbinVal2 += mat[samp2R][samp2C];
            TbinVal3 += mat[samp3R][samp3C];
            TbinVal4 += mat[samp4R][samp4C];

            if (doNAN) {
              binVal1_nan = (binVal1_nan || (mat[samp1R][samp1C] == nan));
              binVal2_nan = (binVal2_nan || (mat[samp2R][samp2C] == nan));
              binVal3_nan = (binVal3_nan || (mat[samp3R][samp3C] == nan));
              binVal4_nan = (binVal4_nan || (mat[samp4R][samp4C] == nan));
            }
//if (mat[samp4R][samp4C]) cout<<"    nonzero isl ang row col: "<<isl<<"  "<<sampAng<<"  "<<samp4R<<"  "<<samp4C<<endl;
          }
          binVal1 += sampArea*TbinVal1;
          binVal2 += sampArea*TbinVal2;
          binVal3 += sampArea*TbinVal3;
          binVal4 += sampArea*TbinVal4;
//cout<<"  Binvalue: "<<binVal4<<endl;
        }

	radHists[ib][isl] = binVal1;
	radHists[ib][isl+Nslices] = binVal2;
	radHists[ib][isl+2*Nslices] = binVal3;
	radHists[ib][isl+3*Nslices] = binVal4;

        if (doNAN) {
          if (binVal1_nan) radHists[ib][isl] = nan;
          if (binVal2_nan) radHists[ib][isl+Nslices] = nan;
          if (binVal3_nan) radHists[ib][isl+2*Nslices] = nan;
          if (binVal4_nan) radHists[ib][isl+3*Nslices] = nan;
        }
          
//cout<<"END angle isl Bin value: "<<angle<<"  "<<isl<<"  "<<ib<<"  "<<binVal4<<endl;
      }

    strStm.str("");
    strStm.clear();
  }

  return radHists;
}


std::vector<double> imgProc::legendreFit(std::vector< std::vector<double> > img, int rebin, const int Nlg, const int Nrad, const int Nrows, const int Ncols) {

  // Check if g matrix already exixsts, else make new one
  string matrix_folder = "/reg/neh/home/khegazy/analysis/legendreFitMatrices/";
  string matrix_fileName = "gMatrix_row-" + to_string(Nrows)
      + "_col-" + to_string(Ncols) + "_Nrad-" + to_string(Nrad)
      + "_Nlg-" + to_string(Nlg) + ".dat";
  if (access((matrix_folder + matrix_fileName).c_str(), F_OK) == -1) {
    cout << "INFO: Making new g matrix\n";
    system(("python " + matrix_folder + "makeLgMatrix.py --NradBins="
          + to_string(Nrad) + " --Ncols=" + to_string(Ncols)
          + " --Nrows=" + to_string(Nrows)
          + " --Nlg=" + to_string(Nlg)).c_str());
  }

  // Import the g matrix
  const int NgMat = Nrows*Ncols*Nrad*Nlg;
  const int Npix = Nrows*Ncols;
  double* gInp = new double[NgMat];
  FILE* inpFile = fopen((matrix_folder + matrix_fileName).c_str(), "rb");
  fread(gInp, sizeof(double), NgMat, inpFile);

  Eigen::MatrixXd gOrig(Npix, Nlg*Nrad);
  for (int ir=0; ir<Npix; ir++) {
    for (int ic=0; ic<Nlg*Nrad; ic++) {
      gOrig(ir,ic) = gInp[ir*Nlg*Nrad + ic];
    }
  }
  //Eigen::Map< Eigen::Matrix<double, Npix, Nlg*Nrad, Eigen::RowMajor> > gOrig(gInp);

  return legendreFit(img, rebin, Nlg, Nrad, Nrows, Ncols, gOrig);
}


std::vector<double> imgProc::legendreFit(std::vector< std::vector<double> > img, int rebin, const int Nlg, const int Nrad, const int Nrows, const int Ncols, const Eigen::MatrixXd g) {

  assert(rebin*Nrows == (int)img.size());
  assert(rebin > 0);

  const int Npix = Nrows*Ncols;
  Eigen::MatrixXd g_inv(Nlg*Nrad, Npix);
  uint svSize = std::min(Npix, Nrad*Nlg);
  Eigen::SparseMatrix<double> svMat(svSize, svSize);
  for (uint i=0; i<svSize; i++) {
    svMat.insert(i,i) = 0;
  }

  // Remove nans
  int Nnans = 0;
  int sum = 0;
  std::pair<int, int> indPair(0,0);
  std::vector< std::pair<int, int> > indPairs;
  std::vector<double> imgFlatVec;
  bool foundNan = false;
  bool inNan = false;
    
  int irr, icc, indR, indC;
  for (int ir=0; ir<Nrows; ir++) {
    for (int ic=0; ic<Ncols; ic++) {
      indR = ir*rebin;
      indC = ic*rebin;

      sum = 0;
      foundNan = false;
      for (irr=0; irr<rebin; irr++) {
        for (icc=0; icc<rebin; icc++) {
          if (img[indR+irr][indC+icc] == NANVAL) {
            foundNan = true;
            break;
          }
          sum += img[indR+irr][indC+icc];
        }

        if (foundNan) {
          break;
        }
      }

      // In in nan and out of nan then start new pair 
      if (inNan && !foundNan) {
        inNan = false;
        indPair.first = ir*Ncols + ic;
      }
      // If found nan then start new pair if it's the first, else append
      if (foundNan) {
        if (!inNan) {
          indPair.second = ir*Ncols + ic;
          indPairs.push_back(indPair);
        }
        inNan = true;
        Nnans++;
      }
      else {
        imgFlatVec.push_back(sum);
      }

    }
  }
  if (!inNan) {
    indPair.second = Npix;
    indPairs.push_back(indPair);
  }

  Eigen::MatrixXd g_masked(Npix - Nnans, Nlg*Nrad);
  int sInd = 0;
  for (uint i=0; i<indPairs.size(); i++) {
    int indDiff = indPairs[i].second - indPairs[i].first;
    g_masked.block(sInd, 0, indDiff, g.cols())
        = g.block(indPairs[i].first, 0, indDiff, g.cols());
    sInd += indDiff;
  }

  Eigen::VectorXd imgFlat(Npix - Nnans);
  for (uint i=0; i<imgFlatVec.size(); i++) {
    imgFlat(i) = imgFlatVec[i];
  }

  // Inverting
  /*
  Eigen::BDCSVD<Eigen::MatrixXd> svd(g_masked, Eigen::ComputeThinU | Eigen::ComputeThinV);
  for (uint i=0; i<svSize; i++) {
    svMat.coeffRef(i,i) = 1./svd.singularValues()(i);
  }
  g_inv = svd.matrixV()*svMat*svd.matrixU().transpose();
  */
  g_inv = tools::SVDinvert(g_masked);

  Eigen::VectorXd lgC = g_inv*imgFlat;

  std::vector<double> legCoeffs(lgC.rows());
  for (uint i=0; i<lgC.rows(); i++) {
    legCoeffs[i] = lgC(i);
  }

  return legCoeffs;
}


TH2F* imgProc::polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order) {

  double m;
  return imgProc::polarLegenderCoeffs(pimgs, imgNums, order, false, m);
}



TH2F* imgProc::polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order, double &max) {

  return imgProc::polarLegenderCoeffs(pimgs, imgNums, order, true, max);
}



TH2F* imgProc::polarLegenderCoeffs(vector< vector< vector< double> > > &pimgs, vector<int> imgNums, int order, bool maxes, double &max) {

  //  All images must have the same row column dimensions 

  double coeff;
  max = -500;
  int min = 500;
  for (uint im=0; im<imgNums.size(); im++) min = imgNums[im]<min ? imgNums[im] : min;

  // Creating histograms for each order
  cout<<"N ybins: "<<pimgs[0].size()<<endl;
  TH2F* legHist = new TH2F(("legHist"+to_string(order)).c_str(), ("legHist"+to_string(order)).c_str(), imgNums.size(), 1, imgNums.size()+1, pimgs[0].size(), 0, pimgs[0].size()+1);

//cout<<"111"<<endl;
  // Calculating legender coeffs and filling the histograms
  for (uint im=0; im<imgNums.size(); im++) {
//cout<<"222"<<endl;
    imgNums[im] += -min + 1; 		// Normalize so the the first number is 1
    for (uint ir=0; ir<pimgs[im].size(); ir++) {
      coeff = imgProc::legendre1dCoeff(pimgs[im][ir], order);
//cout<<"333"<<endl;
      legHist->SetBinContent(imgNums[im], ir+1, coeff);
//cout<<"444"<<endl;
      if (maxes & (coeff > max)) max = coeff;
    }
  }

  return legHist;
}


/*
std::vector<double> imgProc::fitImgLegendre(
        std::vector< std::vector<double> > &img,
        const int Nrows, const int rebinRow,
        const int Ncols, const int rebinCol,
        const int Nlg, const int NradBins, 
        double NANVAL, bool verbose) {

  if (rebinRow != rebinCol) {
    cerr << "ERROR: Currently not setup for rebinRow != rebinCol!!!" << endl;
    exit(0);
  }

  auto multiplyConst = [&](int opt) {
    if (opt==0) {
      return Nrows*Ncols/(rebinRow*rebinCol);
    }
    else if (opt==1) {
      return Nlg*NradBins;
    }
    else {
      cerr << "ERROR: I do not support this option" << endl;
    }
  };

  //assert(Nrows%rebinRow == 0);
  //assert(Ncols%rebinCol == 0);
  const int lgFit_Rows = Nrows/rebinRow;
  const int lgFit_Cols = Ncols/rebinCol;
  // Check if g matrix already exists, else make new one
  string matrix_folder = "/reg/neh/home/khegazy/analysis/legendreFitMatrices/";
  string matrix_fileName = "gMatrix_row-" + to_string(lgFit_Rows)
      + "_col-" + to_string(lgFit_Cols) + "_Nrad-" + to_string(NradBins)
      + "_Nlg-" + to_string(Nlg) + ".dat";
  if (access((matrix_folder + matrix_fileName).c_str(), F_OK) == -1) {
    cout << "INFO: Making new g matrix\n";
    system(("python " + matrix_folder + "makeLgMatrix.py ----NradBins="
          + to_string(NradBins) + " --Ncols=" + to_string(lgFit_Cols)
          + " --Nrows=" + to_string(lgFit_Rows)
          + " --Nlg=" + to_string(Nlg)).c_str());
  }

  // Import the g matrix
  const int NgMat = lgFit_Rows*lgFit_Cols*NradBins*Nlg;
  const int Npix = lgFit_Rows*lgFit_Cols;
  double* gInp = new double[NgMat];
  FILE* inpFile = fopen((matrix_folder + matrix_fileName).c_str(), "rb");
  fread(gInp, sizeof(double), NgMat, inpFile);
  //Eigen::Map< Eigen::Matrix<double, const_cast<int>(Nrows*Ncols/(rebinRow*rebinCol)), Nlg*NradBins, Eigen::RowMajor> > gOrig(gInp);    
  const int m0 = multiplyConst(0);
  Eigen::Map< Eigen::Matrix<double, m0, multiplyConst(1), Eigen::RowMajor> > gOrig(gInp);    

  clock_t begin = clock();
  std::vector<double> legCoeffs = legendreFit(img, rebinRow, Nlg, NradBins, lgFit_Rows, lgFit_Cols, NANVAL, gOrig);
  clock_t end = clock();

  if (verbose) {
    cout << "Time to fit Legendres: " << double(end - begin) << endl;
  }

  return legCoeffs;
}
*/




/*
template<typename T>
vector< vector<double> > imgProc::polarBinningT(cv::Mat mat, cv::Point center, int Nslices, int NmatBins, int NhistBins, float rad) {
  ////// Plots the radial projection as a function of angle. Each angle has 4 hists in cross shape. 
  //////	hists 0,1,2,3 start at angle 0,PI/2,PI,3*PI/2.
  //  Nslices: how many angular steps from 0 to PI/2, total of 4*Nslices hists
  //  NmatBins: number of bins to go out to from center of mat
  //  NhistBins: number of bins in the projection histograms
  //  rad: number value of the last bin for labelling histos
  //  combined: plot all slices on a 2D histo in order as a functin of angle
  //  normalization of binning is done where one mat bin has area 1
  //  OUTPUT: vector of the 4 slices in cross shape for each angle

  if (!mat.isContinuous()) {
    cerr<<"ERROR: Radial Plots only works with continuous matrices, can easily be modified!!!"<<endl;
    exit(0);
  }
  if (NhistBins > NmatBins) {
    cerr<<"ERROR: Cannot make radial plots with NhistBins > NmatBins, cannot pull information out of nowhere"<<endl;
    exit(0);
  }

  vector< vector<double> > radHists;
  bool combineBins = true;

  double angle=0;
  double iterA = PI/(2.0*Nslices);
  double iterR = ((double)NmatBins)/((double)NhistBins);
  ostringstream strStm;
  radHists.resize(NhistBins);
  for (uint irad=0; irad<radHists.size(); irad++) radHists[irad].resize(4*Nslices);
  T* matVals = mat.ptr<T>(0);
  int Nrows = mat.rows;
  if (rad == -1) rad = (2*((double)NhistBins)/((double)Nrows));

  int NradSamp = (int)(3*iterR+1);
  int NangSamp = 1000;
  while (NangSamp%Nslices != 0) NangSamp++;
  NangSamp /= Nslices;
  double sampRadIter = iterR/NradSamp;
  double sampAngIter = iterA/NangSamp;
  double sampCosA, sampSinA, sinA, cosA, sampRad, sampAng, sampArea;
  double binVal1, binVal2, binVal3, binVal4;
  cv::Point samp1,samp2,samp3,samp4;
  CvPoint2D64f point1,point2,point3,point4;

  for (int isl=0; isl<Nslices; isl++) {
    angle = isl*iterA;
    strStm << angle;

    if (combineBins) {  	// Fills each bin by sampling bins within the polar bin's area
      for (int ib=0; ib<NhistBins; ib++) {
        binVal1=binVal2=binVal3=binVal4=0;
        for (int irad=0; irad<NradSamp; irad++) {
          sampRad = ib + irad*sampRadIter + sampRadIter/2;
          sampArea = sampAngIter*(pow(sampRad+sampRadIter/2,2) - pow(sampRad-sampRadIter/2,2));
          for (int iang=0; iang<NangSamp; iang++) {
            sampAng = angle + iang*sampAngIter + sampAngIter/2; 
            sampCosA = cos(sampAng);	sampSinA = sin(sampAng);
            samp1.x = (int)(center.x + sampRad*sampCosA);	samp1.y = (int)(center.y - sampRad*sampSinA);
            samp2.x = (int)(center.x - sampRad*sampSinA);	samp2.y = (int)(center.y - sampRad*sampCosA);
            samp3.x = (int)(center.x - sampRad*sampCosA);	samp3.y = (int)(center.y + sampRad*sampSinA);
            samp4.x = (int)(center.x + sampRad*sampSinA);	samp4.y = (int)(center.y + sampRad*sampCosA);
 
            binVal1 += sampArea*matVals[samp1.y*Nrows + samp1.x];
            binVal2 += sampArea*matVals[samp2.y*Nrows + samp2.x];
            binVal3 += sampArea*matVals[samp3.y*Nrows + samp3.x];
            binVal4 += sampArea*matVals[samp4.y*Nrows + samp4.x];
          }
        }

	radHists[ib][isl] = binVal1;
	radHists[ib][isl+Nslices] = binVal2;
	radHists[ib][isl+2*Nslices] = binVal3;
	radHists[ib][isl+3*Nslices] = binVal4;
      }
    }
    else {	// Takes an interpolation at the center of each bin, treats each bin as the same size
      sinA = sin(angle);
      cosA = cos(angle);
      point1.x=point2.x=point3.x=point4.x=center.x;
      point1.y=point2.y=point3.y=point4.y=center.y;
      for (int ib=0; ib<NhistBins; ib++) {
	radHists[ib][isl] = binVal1;
	radHists[ib][isl+Nslices] = binVal2;
	radHists[ib][isl+2*Nslices] = binVal3;
	radHists[ib][isl+3*Nslices] = binVal4;

        //      COLUMNS			 ROWS
        point1.x += iterR*cosA; 	point1.y += -iterR*sinA;
        point2.x += -iterR*sinA; 	point2.y += -iterR*cosA;
        point3.x += -iterR*cosA; 	point3.y += iterR*sinA;
        point4.x += iterR*sinA; 	point4.y += iterR*cosA;
      }
    }  
  
    strStm.str("");
    strStm.clear();
  }

  return radHists;
}
*/

double imgProc::legendre1dCoeff(TH1F* hist, int order) {
// Calculates the coefficients of legendre polynomial. 
// Function fits top half and lower half, then averages two results
// ASSUMES: legendre polynomial repeats, that is hist goes from 0-2PI

  double coeff=0;
  int binsTot = hist->GetNbinsX();
  int bins = binsTot/2;
  double iter = 2.0/((double)bins);
  double legendre;
  for (int ib=1; ib<=bins; ib++) {
    legendre = tools::legendreIntegral(order, 1-ib*iter, 1-(ib-1)*iter);
    coeff += legendre*hist->GetBinContent(ib);
    coeff += legendre*hist->GetBinContent(binsTot-(ib-1));
  }
  if (fabs(coeff) < 0.0000000001) coeff = 0;
  else coeff /= (2.0*((2.0*(double)order+1.0)/2.0));           // (2*io+1)/2 is normalization of legendre, 2 is average

  return coeff;
}


double imgProc::legendre1dCoeff(vector<double> &vect, int order) {
// Calculates the coefficient of legendre polynomial in this single vector. 
// Function fits top half and lower half, then averages two results
// ASSUMES: legendre polynomial repeats, that is vect goes from 0-2PI

  double coeff=0;
  int binsTot = vect.size();
  int bins = binsTot/2;
  double iter = 2.0/((double)bins);
  if (binsTot%2 == 1) iter = 2.0/(((double)bins)+0.5); 
  double legendre;

  for (int ib=0; ib<bins; ib++) {
    legendre = tools::legendreIntegral(order, 1.0-((double)ib+1.0)*iter, 1.0-(double)ib*iter);
    coeff += legendre*vect[ib];
    coeff += legendre*vect[binsTot-(ib+1)];
  }

  if (binsTot%2 == 1) coeff += vect[bins]*2.0*tools::legendreIntegral(order, -1.0, 1.0-((double)bins)*iter); 

  if (fabs(coeff) < 0.0000000001) coeff = 0;
  else coeff /= (2.0*((2.0*(double)order+1.0)/2.0));           // (2*order+1)/2 is normalization of legendre, 2 is average

  return coeff;
}


double imgProc::legendre1dCoeff(vector<double> &vect, int order, double nan) {
// Calculates the coefficient of legendre polynomial in this single vector. 
// Function fits legendre polynomial and normalizes based upon bins used
// ASSUMES: legendre polynomial repeats, that is vect goes from 0-2PI

  double coeff = 0;
  double norm = 0;
  int bins = vect.size();
  double iter = 2.0/((double)bins);
  double legendre;

  for (int ib=0; ib<bins; ib++) {
    if (vect[ib] != nan) {
      legendre = tools::legendreIntegral(order, 1.0-((double)ib+1.0)*iter, 1.0-(double)ib*iter);
      coeff += legendre*vect[ib];
      norm += legendre*legendre;
    }
  }

  if ((fabs(coeff) < 1e-10) || (norm < 1e-10)) coeff = 0;
  else {
    // (2*order+1)/2 is normalization of legendre, 2 is average
    coeff /= (2.0*((2.0*(double)order+1.0)/2.0));
    coeff /= norm;
  }

  return coeff;
}


double imgProc::legendre1dCosCoeff(vector<double> &vect, int order) {
// Calculates the coefficient of legendre polynomial in this single vector. 
// Function fits top half and lower half, then averages two results
// ASSUMES: legendre polynomial repeats, that is vect goes from 0-2PI

  double coeff=0;
  int binsTot = vect.size();
  int bins = binsTot/2;
  double w = PI/((double)binsTot/2.);
  double legendre;

  for (int ib=0; ib<bins; ib++) {
    legendre = tools::legendreIntegral(order, cos((ib + 1)*w), cos(ib*w));
    coeff += legendre*vect[ib];
    coeff += legendre*vect[binsTot-(ib+1)];
  }

  if (binsTot%2 == 1) coeff += vect[bins]*2.0*tools::legendreIntegral(order, -1.0, cos(bins*w));

  if (fabs(coeff) < 0.0000000001) coeff = 0;
  else coeff /= (2.0*((2.0*(double)order+1.0)/2.0));           // (2*order+1)/2 is normalization of legendre, 2 is average

  return coeff;
}


double imgProc::legendre1dCosCoeff(vector<double> &vect, int order, double nan) {
// Calculates the coefficient of legendre polynomial in this single vector. 
// Function fits legendre polynomial and normalizes based on bins used
// ASSUMES: legendre polynomial repeats, that is vect goes from 0-2PI

  double coeff = 0;
  double norm = 0;
  int bins = vect.size();
  double w = PI/((double)bins/2.);
  double legendre;

  for (int ib=0; ib<bins; ib++) {
    if (vect[ib] != nan) {
      legendre = tools::legendreIntegral(order, cos((ib + 1)*w), cos(ib*w));
      coeff += legendre*vect[ib];
      norm += legendre*legendre;
    }
  }

  if ((fabs(coeff) < 1e-10) || (norm < 1e-10)) coeff = 0;
  // (2*order+1)/2 is normalization of legendre, 2 is average
  else {
    coeff /= (2.0*((2.0*(double)order+1.0)/2.0)); 
    coeff /= norm;
  }

  return coeff;
}


void imgProc::makeMaskRC(double ang, double radI, double radF, int rows, int cols, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask) {

  if (radI>1 || radF>1) {
    cerr<<"ERROR: radI and radF must be less than 1, radii is normalized to 1!!!"<<endl;
    exit(0);
  }

  double r,a,norm,centC,centR;
  centC=cols/2; centR=rows/2;
  norm = (rows+cols)/2/2;
  pair<int, int> p;
  for (int ir=0; ir<rows; ir++) {
    for (int ic=0; ic<cols; ic++) {
      r = sqrt(pow(ic-centC,2)+pow(ir-centR,2))/norm;
      if (r<radF && r>radI) {
        a = atan(fabs((ir-centR)/(ic-centC)));
        if (a<=ang/2.0) {
          p.first = ir;
          p.second = ic;
          pmask.push_back(p);
        }
        if (a>=(PI-ang)/2.0) {
          p.first = ir;
          p.second = ic;
          nmask.push_back(p);
        }
      }
    }
  }
  return;
}


void imgProc::makeMaskRad(double ang, double radI, double radF, int rads, int angs, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask) {

  if (radI>1 || radF>1) {
    cerr<<"ERROR: radI and radF must be less than 1, radii is normalized to 1!!!"<<endl;
    exit(0);
  }

  double r,a;
  pair<int, int> p;
  for (int ir=1; ir<=rads; ir++) {
    r = ((double)ir)/((double)rads);
    if (r>radF || r<radI) continue;
    for (int ia=1; ia<=angs; ia++) {
      a = fmod(((double)ia)*(2*PI)/((double)angs), PI);
      if (a<=ang/2 || a>(PI-ang/2)) {
        p.first = ir-1;
        p.second = ia-1;
        pmask.push_back(p);
      }
      if (ang/2>=fabs(PI/2-a)) {
        p.first = ir-1;
        p.second = ia-1;
        nmask.push_back(p);
      }
    }
  }
  return;
}


void imgProc::makeMask(bool rad, double ang, double radI, double radF, vector< vector<double> > &vect, vector< pair<int,int> > &pmask, vector< pair<int,int> > &nmask) {

  if (rad) makeMaskRad(ang, radI, radF, vect.size(), vect[0].size(), pmask, nmask);
  else makeMaskRC(ang, radI, radF, vect.size(), vect[0].size(), pmask, nmask);
  return;
}


double imgProc::Aparam(vector<double> vect) {

  int quart = vect.size()/4;
  double param1=0;
  double param2=0;
  for (int k=0; k<(int)vect.size(); k++) {
    if ((k+1>quart/2 && k<quart*1.5) || (k+1>quart*2.5 && k<quart*3.5)) {param1 += vect[k];}
    else {param2 += vect[k];}
  }

  //cout<<"parameters: "<<param1<<"  "<<param2<<endl;
  return (param1-param2);///(param1+param2);
  //return param1/param2;
}


vector< vector<double> > imgProc::imgRotatePId2(vector< vector<double> > origImg) {

  int rpix = origImg[0].size();
  int cpix = origImg.size();
  vector< vector<double> > rotImg(rpix);
  for (int ir=0; ir<rpix; ir++) {
    rotImg[ir].resize(cpix, 0);
    for (int ic=0; ic<cpix; ic++) {
      rotImg[ir][ic] = origImg[ic][rpix - ir - 1];
    }
  }

  return rotImg;
}


double imgProc::Aparam(vector< vector<double> > vect, vector< pair<int,int> > pmask, vector< pair<int,int> > nmask) {

  double param1=0;
  double param2=0;
  for (uint k=0; k<pmask.size(); k++) param1 += vect[pmask[k].first][pmask[k].second];
  for (uint k=0; k<nmask.size(); k++) param2 += vect[nmask[k].first][nmask[k].second];

  cout<<"parameters: "<<param1<<"  "<<param2<<endl;
  //return (param1-param2);///(param1+param2);
  return param1/param2;
}




















//imgProc::radProcTool::radProcTool(string inpIndPath) {
//  indPath = inpIndPath;
//}

imgProc::radProcTool::radProcTool(int shellWidth, int maxRad) {
  indPath = "constructor";
  makeIndices(shellWidth, maxRad);
}

void imgProc::radProcTool::makeIndices(int shellWidth, int maxRad) {
  int rad;
  makeMaxRad = maxRad;
  for (int ir=-1*maxRad; ir<=maxRad; ir++) {
    for (int ic=-1*maxRad; ic<=maxRad; ic++) {
      rad = std::round(std::sqrt(ir*ir + ic*ic));

      if (rad <= maxRad) {
        std::vector<int> vInds{ir, ic};
        allIndices[shellWidth][rad].push_back(vInds);
      }
    }
  }
}

bool imgProc::radProcTool::checkForIndices(int shellWidth, int rad) {
  if (allIndices.find(shellWidth) != allIndices.end()) {
    if (allIndices[shellWidth].find(rad) != allIndices[shellWidth].end()) {
      return true;
    }
  }

  if (indPath.compare("constructor") == 0) {
    makeIndices(shellWidth, rad);
  }
  else {
    importIndices(shellWidth, rad);
  }
  return true;
}

void imgProc::radProcTool::importIndices(int shellWidth, int rad) {

  string indFileName;
  string fileName = "NULL";
  bool found = false;
  bool made = false;
  vector<int16_t> inpInds;

  /////  Find file containing indices  /////
  found = false;
  made = false;
  indFileName = "indicesInRange[" + to_string(rad) + "," + to_string(rad+shellWidth) + "]";
  while (!found) {
    DIR* dir = opendir(indPath.c_str());
    struct dirent* ent;
    while ((ent = readdir(dir)) != NULL) {
      string curFileName(ent->d_name);
      if (indFileName.compare(curFileName.substr(0,indFileName.length())) == 0) {
        fileName = curFileName;
        found = true;
        break;
      }
    }
    closedir(dir);

    // Make a new file of indices if needed
    if (!found) {
      if (made) break;

      cout << "INFO: Making new inds\n";
      cout<<"command: "<<("python " + indPath + "getIndices.py --rMin " + to_string(rad) + " --rMax " + to_string(rad+shellWidth))<<endl;
      system(("python " + indPath + "getIndices.py --rMin "
            + to_string(rad) + " --rMax " + to_string(rad+shellWidth)).c_str());
      made = true;
    }
  }
  if (!found && made) {
    std::cerr << "ERROR: Cannot find file beginning with "
      + indFileName + " even though it was made!!!\n";
    exit(0);
  }

  // Get number of indices
  auto iPos = fileName.find("Size") + 5;
  int size = std::atoi(fileName.substr(iPos, fileName.find("]", iPos) - iPos).c_str());

  // Import indices
  inpInds.resize(2*size);
  save::importDat<int16_t>(inpInds, indPath + fileName);

  
  for (uint i=0; i<inpInds.size(); i+=2) {
    std::vector<int> vInds{inpInds[i], inpInds[i+1]};
    allIndices[shellWidth][rad].push_back(vInds);
  }
}


void imgProc::radProcTool::getSmearedDist(std::map<int, double> &smearedDist,
    vector<double> &vals, double stdev, bool verbose) {

  /////  Initialize variables  /////
  std::map<int, double> dist;
  smearedDist.clear();

  // Scale stdev
  double var = std::pow(stdev/5, 2);

  for (uint i=0; i<vals.size(); i++) {
    if (vals[i] == NANVAL) continue;
    if (dist.find((int)vals[i]) == dist.end()) {
      dist[(int)vals[i]] = 1;
    }
    else {
      dist[(int)vals[i]]++;
    }
  }

  int begin = dist.begin()->first - 4*stdev;
  int end = dist.rbegin()->first + 4*stdev;
  auto dItr = dist.begin();
  for (int ind=begin; ind<=end; ind++) {
    smearedDist[ind] = 0;
    dItr = dist.begin();
    while ((ind - dItr->first)/stdev > 4) {
      dItr++;
    }

    smearedDist[ind] = 0;
    while (fabs(ind - dItr->first)/stdev <= 4) {
      smearedDist[ind] += dItr->second
              *exp(-1*pow(ind-dItr->first, 2)/(2*var));

      dItr++;
      if (dItr == dist.end()) break;
    }
  }
}




double imgProc::radProcTool::getMMean(vector<double> vals, 
    vector<int> ord, double range, bool verbose) {
  // Returns the value that has the largest bin height as the mean

  /*
  if (ord.size() < 250) {
    return getMean(vals);
  }
  */
  //cout<<"start getmmean"<<endl;
  int maxCenter = vals[ord[ord.size()*2/3]];
  double weight = 0;
  // Start at median and go to lower values
  int center = (int)vals[ord[0]];
  if (verbose) cout<<"ORIGINFO: "<<ord.size()<<"  "<<maxCenter<<"  "<<vals[ord[ord.size()-1]]<<endl;

  if ((center == maxCenter) || (center == (int)vals[ord[ord.size()-1]])) {
    return center;
  }

  int ind = 0;
  for (ind=0; center - vals[ord[ind]] > range; ind++) {
      if (verbose) cout<<"BAD: "<<ind<<"  "<<center<<"  "<<vals[ord[ind]]<<"  "<<range<<endl;
      cerr << "ERROR: Should not be here Reached end of indices without exiting!!!\n";
  }
  if (verbose)cout<<"h2"<<endl;

  for (; fabs(center - vals[ord[ind]]) <= range; ind++){
    weight++;
    if (ind + 1 == (int)ord.size()) break;
  }

  if (verbose) cout<<"beginning vals: "<<center<<"  "<<weight<<endl;
  double maxWeight = weight;
  int bestCenter = center;

 //cout<<"h3"<<endl;
  // Search smaller values until we certainly passed the maximum

  while ((center < maxCenter )
      && (center != vals[ord[ord.size()-1]])) {
    center++;
    weight = 0;
    for (ind=0; center - vals[ord[ind]] > range; ind++){
    if (verbose) cout<<"test: "<<ind<<"  "<<center<<" "<<vals[ord[ind]]<<"  "<<range<<"  "<<ord.size()<<endl;
      if (ind + 1 == (int)ord.size()) {
        cerr << "ERROR: Reached end of indices without exiting!!!\n";
        break;
      }
    }
    if (verbose) cout<<"\tcheck1: "<<ind<<" "<<vals[ord[ind]]<<endl;

    for (; fabs(center - vals[ord[ind]]) <= range; ind++) {
      weight++;
      if (ind + 1 == (int)ord.size()) break;
    }
    //cout<<"\tcheck2: "<<ind<<" "<<vals[ord[ind]]<<endl;
  //cout<<"h6"<<endl;
    if (verbose) cout<<"\tlooping: "<<center<<"  "<<weight<<"  "<<maxWeight<<endl;
    if (weight > maxWeight) {
      maxWeight = weight;
      bestCenter = center;
    }
  }
  //cout<<"done"<<endl;

  return bestCenter;
}


double imgProc::radProcTool::getMean(vector<double> vals) {
  double sum = 0;
  double count = 0;
  for_each(vals.begin(), vals.end(),
      [&sum, &count] (double d)
      { if (d != NANVAL) {
          sum += d;
          count += 1;
        }
      });

  return sum/count;
}


int imgProc::radProcTool::getNleftEntries(vector<double> vals, vector<int> ord, double mean) {
  int count = 0;
  for_each(ord.begin(), ord.end(),
      [&count, vals, mean] (int ind)
      { if (vals[ind] > mean) {
          return;
          }
        if (vals[ind] != NANVAL) {
          count += 1;
        }
      });
  return count;
}


double imgProc::radProcTool::getSTDev(vector<double> vals, double mean) {
  double var = 0;
  double count = 0;
  for_each(vals.begin(), vals.end(),
      [&var, &count, mean] (double val)
      { if (val != NANVAL) {
          var += pow(val - mean, 2);
          count += 1;
        }
      });
  return sqrt(var/count);
}


double imgProc::radProcTool::getLeftSTDev(vector<double> vals, vector<int> ord, double mean) {
  double var = 0;
  double count = 0;
  for_each(ord.begin(), ord.end(),
      [&var, &count, vals, mean] (int ind)
      { if (vals[ind] > mean) {
          return;
          }
        if (vals[ind] != NANVAL) {
          var += pow(vals[ind] - mean, 2);
          count += 1;
        }
      });
  return sqrt(var/count);
}


int imgProc::radProcTool::getNrightEntries(vector<double> vals, vector<int> ord, double mean) {
  int count = 0;
  for_each(ord.rbegin(), ord.rend(),
      [&count, vals, mean] (int ind)
      { if (vals[ind] < mean) {
          return;
          }
        if (vals[ind] != NANVAL) {
          count += 1;
        }
      });
  return count;
}


int getNleftEntries(vector<double> vals, vector<int> ord, double mean) {
  int count = 0;
  for_each(ord.begin(), ord.end(),
      [&count, vals, mean] (int ind)
      { if (vals[ind] > mean) {
          return;
          }
        if (vals[ind] != NANVAL) {
          count += 1;
        }
      });
  return count;
}


int imgProc::radProcTool::getLeftOutliers(vector<double> &vals, vector<int> &orderedInds,
    double mean, double stdev, double Nstdev) {

  int Nfront = 0;
  for (auto ind : orderedInds) {
    if (vals[ind] < mean - stdev*Nstdev) {
      Nfront += 1;
    }
    else break;
  }
  return Nfront;
}


int imgProc::radProcTool::getRightOutliers(vector<double> &vals, vector<int> &orderedInds,
    double mean, double stdev, double Nstdev) {

  int Nfront = 0;
  auto ind = orderedInds.rbegin();
  for (; vals[*ind] > mean + Nstdev*stdev; ind++) {
    Nfront += 1;
  }
  return Nfront;
}

std::vector< std::vector<double> > imgProc::radProcTool::removeOutliersSimple(
    vector< vector<double> > &image, 
    float centerR_f, float centerC_f, int buffer,
    int maxRad, int shellWidth, 
    bool removeLowPoly, int Npoly,
    double stdCut, int stg, double outlierMapSTDcut, 
    bool getOutlierImage,  
    bool verbose, PLOTclass* pltVerbose) {

  int centerR = (int)centerR_f;
  int centerC = (int)centerC_f;
  /*
  for (int ir=400; ir<405; ir++) {
    for (int ic=400; ic<405; ic++) {
      cout<<ir<<" / "<<ic<<" :  "<<image[ir][ic]<<endl;
    }
  }
  */

  if (verbose) cout << "Entered removeOutliers\n";
  std::map<int, double> smearedDist;
  vector< vector<double> > removedOutliers(image.size());
  vector< vector<double> > outlierSTDmap(image.size());
  for (uint ir=0; ir<removedOutliers.size(); ir++) {
    removedOutliers[ir].resize(image[ir].size(), 0);
    outlierSTDmap[ir].resize(image[ir].size(), 0);
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;

  // Plotting variables
  vector<double> pltPixDist;

  vector<double> vals, angles;
  vector< vector<int> > indices;
  vector< vector<int> >* INDS;
    cerr<<"WARNING: not subtracting low poly"<<endl;
  for (int rad=0; rad<maxRad; rad+=shellWidth) {

    //if (rad != 64) continue;
    // Check if indices have been imported
    if (verbose) cout << "\tChecking for indices " 
                      << shellWidth << " " << rad << "  ...  ";
    checkForIndices(shellWidth, rad);
    if (verbose) cout << "found!\n";

    // Filter out bins that are off the image or contain nans
    if (verbose) cout << "\tFilling indices of non-NANVAL bins on image\n";
    int ir, ic;
    vals.clear();
    angles.clear();
    indices.clear();
    INDS = &allIndices[shellWidth][rad];
    for (uint i=0; i<(*INDS).size(); i++) {
      ir = centerR + (*INDS)[i][0];
      ic = centerC + (*INDS)[i][1];
      if ((ir + buffer < (int)image.size()) && (ir >= buffer)
          && (ic + buffer < (int)image[0].size()) && (ic >= buffer)) {
        if (image[ir][ic] != NANVAL) {
          vals.push_back(image[ir][ic]);
          angles.push_back(tools::getAzmAngle((*INDS)[i][1], (*INDS)[i][0]));
          std::vector<int> vInds{ir, ic};
          indices.push_back(vInds);
        }
      }
    }

    // Skip if vals is full of NANVAL
    if (vals.size() < 2) continue;

 
    if (removeLowPoly) {
      // Subtract background fluctuations
      X.resize(indices.size(), Npoly+1);
      Y.resize(indices.size(), 1);
      for (uint i=0; i<angles.size(); i++) {
        for (int p=0; p<Npoly; p++) {
          X(i,p+1) = std::pow(angles[i], 2*p+1);
        }
        X(i,0) = 1;
        Y(i,0) = vals[i];
      }
      Eigen::MatrixXd weights = tools::normalEquation(X, Y);


      for (uint i=0; i<angles.size(); i++) {
        for (int p=0; p<Npoly; p++) {
          vals[i] -= weights(p+1)*std::pow(angles[i], 2*p+1);
          image[indices[i][0]][indices[i][1]] 
              -= weights(p+1)*std::pow(angles[i], 2*p+1);
        }
      }
    }

    // Order indices based on pixel values
    if (verbose) cout << "\tOrdering indices\n";
    std::vector<int> removedInds;
    std::vector<double> removedVals;
    vector<int> orderedInds(vals.size());
    iota(orderedInds.begin(), orderedInds.end(), 0);
    sort(orderedInds.begin(), orderedInds.end(),
        [&vals](int i1, int i2)
        {return vals[i1] < vals[i2];});


    ///////////////////////////////////
    /////  Remove pixel outliers  /////
    ///////////////////////////////////

    /////  Remove left outliers  /////
    //cout<<"111"<<endl;
    //cout<<"CHECKING: "<<stg<<"  "<<rad<<"  "<<vals.size()<<endl;

    double mean = getMean(vals);
    double stdev = getSTDev(vals, mean);
    /*
    if (rad %25 == 0) {
      cout<<"mean/std "<<rad<<" : "<<mean<<" / "<<stdev<<endl;
    }
    */

    int ind = 0;
    vector<double> rms;
    while (fabs(vals[orderedInds[ind]] - mean) > stdCut*stdev) {
      removedInds.push_back(orderedInds[ind]);
      removedVals.push_back(vals[orderedInds[ind]]);
      rms.push_back(vals[orderedInds[ind]]);
      vals[orderedInds[ind]] = NANVAL;
      image[indices[orderedInds[ind]][0]][indices[orderedInds[ind]][1]] = NANVAL;
      ind++;
    }
 
    ind = orderedInds.size() - 1;
    while (fabs(vals[orderedInds[ind]] - mean) > stdCut*stdev) {
      removedInds.push_back(orderedInds[ind]);
      removedVals.push_back(vals[orderedInds[ind]]);
      rms.push_back(vals[orderedInds[ind]]);
      vals[orderedInds[ind]] = NANVAL;
      image[indices[orderedInds[ind]][0]][indices[orderedInds[ind]][1]] = NANVAL;
      ind--;
    }

    //if (rad % 25 == 0) {
    //  cout<<"Count: "<<rms.size()<<endl;
    /*
      for (uint k=0; k<rms.size(); k++) {
        cout<<rms[k]<<"   ";
      }
      cout<<endl;
    */
    //}


    /////  Outlier images  /////
    
    // Image of outliers cut on right end
    if (getOutlierImage || pltVerbose) {
      for (uint i=0; i<removedInds.size(); i++) {
        removedOutliers[indices[removedInds[i]][0]]
                       [indices[removedInds[i]][1]] = removedVals[i];
      }
    }

    // Map of hot pixel stdevs using left mean/stdev with lower thresholds
    if (getOutlierImage || pltVerbose) {
      double oMean = getMean(vals);
      double oSTDev = getSTDev(vals, oMean);
      /*
      for (uint i=0; i<indices.size(); i++) {
        outlierSTDmap[indices[i][0]]
                     [indices[i][1]] = pow((vals[i] - meanLeft)/stdevLeft, 4);
      }
      */

      for (uint i=0; i<removedInds.size(); i++) {
        outlierSTDmap[indices[removedInds[i]][0]]
                     [indices[removedInds[i]][1]] 
                        = pow((removedVals[i] - oMean)/oSTDev, 4);
      }
      /*
      auto inds = orderedInds.rbegin();
      while ((vals[*inds] - oMean)/oSTDev > outlierMapSTDcut) {
        outlierSTDmap[indices[*inds][0]]
                     [indices[*inds][1]] = pow((vals[*inds] - oMean)/oSTDev, 4);
        inds++;
      }
      */
    }

    /////  Distribution of pixel values for current ring  /////
    if ((verbose && pltVerbose)) {
      //cout << "\tSaving value distributions ... ";
      std::vector<double> outP;
      for_each(vals.begin(), vals.end(), [&outP](double d)
          { if (d != NANVAL) {
              outP.push_back(d);
            }
          });
      save::saveDat<double>(outP,  
          "./plots/radialPixelDist/data/vals_" 
          + to_string(stg) + "_" 
          + to_string(rad) + "_"
          + to_string(stdev) + ".dat");
      //cout << "saved\n";
    }

    /////  Polar lineout of ring  /////
    if (verbose && pltVerbose) {
      //cout << "\tSaving polar lineout ... ";
      vector<double> outP = getPolarLineOut(
          &image, centerR, centerC, rad, shellWidth, -1);
      save::saveDat<double>(outP, 
                  "plots/data/polLO_" + to_string(stg) 
                  + "_" + to_string(rad) + 
                  + "[" + to_string(outP.size()) + "].dat");
    }
  }

  return outlierSTDmap;
}


 
std::vector< std::vector<double> > imgProc::radProcTool::removeOutliers(
    vector< vector<double> > &image, 
    int centerR, int centerC, int buffer,
    int maxRad, int shellWidth, int Npoly,
    double stdIncludeLeft, double distSTDratioLeft,
    double stdCutLeft, int meanBinSize,
    double stdIncludeRight, double distSTDratioRight,
    double stdChangeRatio, double stdCutRight,
    int stg, double outlierMapSTDcut, 
    bool getOutlierImage, bool verbose, 
    PLOTclass* pltVerbose, TH1F** radPixHistos) {


  if (verbose) cout << "Entered removeOutliers\n";
  std::map<int, double> smearedDist;
  vector< vector<double> > removedOutliers(image.size());
  vector< vector<double> > outlierSTDmap(image.size());
  for (uint ir=0; ir<removedOutliers.size(); ir++) {
    removedOutliers[ir].resize(image[ir].size(), 0);
    outlierSTDmap[ir].resize(image[ir].size(), 0);
  }
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> X;
  Eigen::Matrix<double, Eigen::Dynamic, 1> Y;

  // Plotting variables
  vector<double> pltPixDist;
  int pdInd;

  vector<double> vals, angles;
  vector< vector<int> > indices;
  vector< vector<int> >* INDS;
  for (int rad=0; rad<maxRad; rad+=shellWidth) {

    //if (rad != 64) continue;
    // Check if indices have been imported
    if (verbose) cout << "\tChecking for indices " 
                      << shellWidth << " " << rad << "  ...  ";
    checkForIndices(shellWidth, rad);
    if (verbose) cout << "found!\n";

    for (int rep=0; rep<NshellOutlierLoops; rep++) { 

      // Filter out bins that are off the image or contain nans
      if (verbose) cout << "\tFilling indices of non-NANVAL bins on image\n";
      int ir, ic;
      vals.clear();
      angles.clear();
      indices.clear();
      INDS = &allIndices[shellWidth][rad];
      for (uint i=0; i<(*INDS).size(); i++) {
        ir = centerR + (*INDS)[i][0];
        ic = centerC + (*INDS)[i][1];
        if ((ir + buffer < image.size()) && (ir >= buffer)
            && (ic + buffer < image[0].size()) && (ic >= buffer)) {
          if (image[ir][ic] != NANVAL) {
            vals.push_back(image[ir][ic]);
            angles.push_back(tools::getAzmAngle((*INDS)[i][1], (*INDS)[i][0]));
            std::vector<int> vInds{ir, ic};
            indices.push_back(vInds);
          }
        }
      }

      // Skip if vals is full of NANVAL
      if (vals.size() < 2) continue;

   
      // Subtract background fluctuations
      X.resize(indices.size(), Npoly);
      Y.resize(indices.size(), 1);
      for (uint i=0; i<angles.size(); i++) {
        for (int p=0; p<Npoly; p++) {
          X(i,p) = std::pow(angles[i], p);
        }
        Y(i,0) = vals[i];
      }
      Eigen::MatrixXd weights = tools::normalEquation(X, Y);

      for (uint i=0; i<angles.size(); i++) {
        for (int p=0; p<Npoly; p++) {
          vals[i] -= weights(p)*std::pow(angles[i], p);
        }
      }


      // Order indices based on pixel values
      if (verbose) cout << "\tOrdering indices\n";
      std::vector<int> removedInds;
      std::vector<double> removedVals;
      vector<int> orderedInds(vals.size());
      iota(orderedInds.begin(), orderedInds.end(), 0);
      sort(orderedInds.begin(), orderedInds.end(),
          [&vals](int i1, int i2)
          {return vals[i1] < vals[i2];});




      ///////////////////////////////////
      /////  Remove pixel outliers  /////
      ///////////////////////////////////

      /////  Remove left outliers  /////
      //cout<<"111"<<endl;
      //cout<<"CHECKING: "<<stg<<"  "<<rad<<"  "<<vals.size()<<endl;

      double mean = getMean(vals);
      double stdev = getSTDev(vals, mean);
      getSmearedDist(smearedDist, vals, stdev);

      //double meanLeft = getMMean(vals, orderedInds, stdev/4, ((rad==28) && verbose && false));
      double meanLeft;
      double maxVal = 0;
      if (pltVerbose) {
        pltPixDist.resize(smearedDist.size(),0);
        pdInd = 0;
      }

      for (auto dItr : smearedDist) {
        if (dItr.second > maxVal*1.1) {
          maxVal = dItr.second;
          meanLeft = dItr.first;
        }

        if (pltVerbose) {
          pltPixDist[pdInd] = dItr.second;
          pdInd++;
        }
      }

      if (pltVerbose && (rep == NshellOutlierLoops-1)) {
        pltVerbose->MyC->cd();
        pltVerbose->MyC->SetLogy(1);
        std::string pltName = 
            "./plots/radialPixelDist/vals_" 
            + to_string(stg) + "_" 
            + to_string(rad) + "_"
            + to_string(meanLeft) + "_"
            + to_string(stdev) + "_"
            + to_string(mean);

        TH1F* h = pltVerbose->plot1d(pltPixDist, pltName,
            xSpan, 
            to_string(smearedDist.begin()->first) + ","
            + to_string(smearedDist.rbegin()->first));
        h->SetMinimum(0.5);

        h->Draw();
        TLine *lineL = new TLine(meanLeft, 0, meanLeft, 1000);
        lineL->SetLineColor(kRed);
        lineL->SetLineWidth(3);
        lineL->Draw("SAME");
        TLine *lineM = new TLine(mean, 0, mean, 1000);
        lineM->SetLineColor(kBlue);
        lineM->SetLineWidth(3);
        lineM->Draw("SAME");
        pltVerbose->MyC->Print((pltName + ".png").c_str());

        delete lineL;
        delete lineM;
        delete h;
      }
   
      //cout<<"222"<<endl;
      double stdevLeft = getLeftSTDev(vals, orderedInds, meanLeft);
      //cout<<"333"<<endl;
      double Nentries = getNleftEntries(vals, orderedInds, meanLeft);
      //cout<<"444"<<endl;
      double Noutliers = getLeftOutliers(vals, orderedInds,
                        meanLeft, stdevLeft, stdIncludeLeft);
      //cout<<"555"<<endl;

      if (verbose) cout << "\tRemoving left outliers\n\t" 
                        << meanLeft << "\t" << stdevLeft << "\t" 
                        << (1 - Noutliers/Nentries) << endl;
      while (((1 - Noutliers/Nentries) < distSTDratioLeft)
              && (orderedInds.size() > 100)) {
        for (int i=0; i<max((int)(orderedInds.size()*fracShellSTDcutLeft), 1); i++) {
          removedInds.push_back(orderedInds[0]);
          removedVals.push_back(vals[orderedInds[0]]);
          vals[orderedInds[0]] = NANVAL;
          image[indices[orderedInds[0]][0]][indices[orderedInds[0]][1]] = NANVAL;
          orderedInds.erase(orderedInds.begin());
        }

        stdevLeft = getLeftSTDev(vals, orderedInds, meanLeft);
        Nentries = getNleftEntries(vals, orderedInds, meanLeft);
        Noutliers = getLeftOutliers(vals, orderedInds,
                      meanLeft, stdevLeft, stdIncludeLeft);
      }
      while (vals[orderedInds[0]] < meanLeft - stdCutLeft*stdevLeft) {
        removedVals.push_back(vals[orderedInds[0]]);
        removedInds.push_back(orderedInds[0]);
        vals[orderedInds[0]] = NANVAL;
        image[indices[orderedInds[0]][0]][indices[orderedInds[0]][1]] = NANVAL;
        orderedInds.erase(orderedInds.begin());
      }


      /////  Remove right outliers  /////
      double stdevPrev;
      double meanRight = getMean(vals);
      double stdevRight = getSTDev(vals, meanRight);
      Nentries = getNrightEntries(vals, orderedInds, meanRight);
      Noutliers = getRightOutliers(vals, orderedInds,
                    meanRight, stdevRight, stdIncludeRight);
      
      if (verbose) cout << "\tRemoving right outliers\n\t" 
                        << meanRight << "\t" << stdevRight << "\t" 
                        << (1 - Noutliers/Nentries) << endl;
      while (((1 - Noutliers/Nentries) < distSTDratioRight) 
              && (orderedInds.size() > 100)) {

        for (int i=0; i<(int)(orderedInds.size()*fracShellSTDcutRight); i++) {
          int ind = orderedInds.size() - 1;
          removedVals.push_back(vals[orderedInds[ind]]);
          removedInds.push_back(orderedInds[ind]);
          vals[orderedInds[ind]] = NANVAL;
          image[indices[orderedInds[ind]][0]][indices[orderedInds[ind]][1]] = NANVAL;
          orderedInds.erase(orderedInds.begin() + ind);
          //cout<<"REMOVING!!!!!!!!!!!"<<endl;
        }

        // Break out of the the loop if the stdev does not change appreciably
        meanRight = getMean(vals);
        stdevPrev = stdevRight;
        stdevRight = getSTDev(vals, meanRight);
        if ((stdevPrev-stdevRight)/stdevRight < stdChangeRatio) {
          if (verbose) cout << "\tbreaking by ratio: " + to_string((stdevPrev-stdevRight)/stdevRight) << endl;
          break;
        }

        Noutliers = getRightOutliers(vals, orderedInds,
                      meanRight, stdevRight, stdIncludeRight);
        Nentries = getNrightEntries(vals, orderedInds, meanRight);
      }
      int ind = orderedInds.size() - 1;
      while (vals[orderedInds[ind]] > meanRight + stdCutRight*stdevRight) {
        removedVals.push_back(vals[orderedInds[ind]]);
        removedInds.push_back(orderedInds[ind]);
        vals[orderedInds[ind]] = NANVAL;
        image[indices[orderedInds[ind]][0]][indices[orderedInds[ind]][1]] = NANVAL;
        orderedInds.erase(orderedInds.begin() + ind);
        ind = orderedInds.size() - 1;
      }


      //////////////////////////////
      /////  Plotting results  /////
      //////////////////////////////

      if (rep == NshellOutlierLoops-1) {
        /////  Outlier images  /////
        if (radPixHistos) {
          for (uint i=0; i<vals.size(); i++) {
            if (vals[i] != NANVAL) {
              radPixHistos[rad]->Fill(vals[i], 1);
            }
          }
        }

        // Image of outliers cut on right end
        if (getOutlierImage || pltVerbose) {
          for (int i=0; i<removedInds.size(); i++) {
            removedOutliers[indices[removedInds[i]][0]]
                           [indices[removedInds[i]][1]] = removedVals[i];
          }
        }

        // Map of hot pixel stdevs using left mean/stdev with lower thresholds
        if (getOutlierImage || pltVerbose) {
          double oMean = getMean(vals);
          double oSTDev = getSTDev(vals, oMean);
          /*
          for (uint i=0; i<indices.size(); i++) {
            outlierSTDmap[indices[i][0]]
                         [indices[i][1]] = pow((vals[i] - meanLeft)/stdevLeft, 4);
          }
          */

          for (int i=0; i<removedInds.size(); i++) {
            outlierSTDmap[indices[removedInds[i]][0]]
                         [indices[removedInds[i]][1]] = pow((removedVals[i] - oMean)/oSTDev, 4);
          }
          auto inds = orderedInds.rbegin();
          while ((vals[*inds] - oMean)/oSTDev > outlierMapSTDcut) {
            outlierSTDmap[indices[*inds][0]]
                         [indices[*inds][1]] = pow((vals[*inds] - oMean)/oSTDev, 4);
            inds++;
          }
        }

        /////  Distribution of pixel values for current ring  /////
        if ((verbose && pltVerbose)) {
          //cout << "\tSaving value distributions ... ";
          std::vector<double> outP;
          for_each(vals.begin(), vals.end(), [&outP](double d)
              { if (d != NANVAL) {
                  outP.push_back(d);
                }
              });
          save::saveDat<double>(outP,  
              "./plots/radialPixelDist/data/vals_" 
              + to_string(stg) + "_" 
              + to_string(rad) + "_"
              + to_string(meanLeft) + "_"
              + to_string(stdev) + "_"
              + to_string(meanRight) + ".dat");
          //cout << "saved\n";
        }

        /////  Polar lineout of ring  /////
        if (verbose && pltVerbose) {
          //cout << "\tSaving polar lineout ... ";
          vector<double> outP = getPolarLineOut(&image, centerR, centerC, rad, shellWidth, 180);
          save::saveDat<double>(outP, 
                      "./results/polLO_" + to_string(stg) 
                      + "_" + to_string(rad) + ".dat");

          std::fill(outP.begin(), outP.end(), 0);
          std::vector<double> outPcount(outP.size(), 0);
          for (auto oInd : orderedInds) {
            outP[(int)(angles[oInd]*180/(2*PI))] += vals[oInd];
            outPcount[(int)(angles[oInd]*180/(2*PI))] ++;
          }
          for (uint i=0; i<outP.size(); i++) {
            if (outPcount[i]) {
              outP[i] /= outPcount[i];
            }
            else {
              outP[i] = 0;
            }
          }
          save::saveDat<double>(outP, 
                      "./results/polNormLO_" + to_string(stg) 
                      + "_" + to_string(rad) + ".dat");
          //cout << "saved\n";
        }
      }
    }
  }


  ////////////////////
  /////  Saving  /////
  ////////////////////

  if (pltVerbose) {
    save::saveDat<double>(removedOutliers, "./results/outlierRemoved_"
                    + to_string(stg) + "_["
                    + to_string(removedOutliers.size()) + ","
                    + to_string(removedOutliers[0].size()) + "].dat");
    save::saveDat<double>(outlierSTDmap, "./results/outlierSTD_"
                    + to_string(stg) + "_["
                    + to_string(outlierSTDmap.size()) + ","
                    + to_string(outlierSTDmap[0].size()) + "].dat");
  }

  return outlierSTDmap;
}


vector<double> imgProc::radProcTool::getPolarLineOut(vector< vector<double> >* image, 
    int centerR, int centerC, int rad, 
    int shellWidth, int NangleBins, bool verbose) {

  // Check if indices have been imported
  if (verbose) cout << "\tChecking for indices " 
                    << shellWidth << " " << rad <<endl;
  checkForIndices(shellWidth, rad);

  // Filter out bins that are off the image or contain nans
  if (verbose) cout << "\tFilling indices of non-NANVAL bins on image\n";
  double ang;
  int ir, ic, angInd;
  vector<double> lineOut;
  if (NangleBins != -1) {
    lineOut.resize(NangleBins, 0);
  }
  vector<double> lineOutCount;
  if (NangleBins != -1) {
    lineOutCount.resize(NangleBins, 0);
  }
  vector< vector<int> >* INDS = &allIndices[shellWidth][rad];
  for (uint k=0; k<(*INDS).size(); k++) {
    ir = centerR+(*INDS)[k][0];
    ic = centerC+(*INDS)[k][1];
    if ((ir < (*image).size()) && (ir >= 0)
        && (ic < (*image)[ir].size()) && (ic >= 0)) {
      if ((*image)[ir][ic] != NANVAL) {

        if (NangleBins == -1) {
          lineOut.push_back((*image)[ir][ic]);
        }
        else {
          ang = tools::getAzmAngle((*INDS)[k][1], (*INDS)[k][0]);
          /*
          ang = atan(((double)(*INDS)[k][0])/((double)(*INDS)[k][1]));
          if (((*INDS)[k][0] >= 0) && ((*INDS)[k][1] <= 0)) ang += PI;
          if (((*INDS)[k][0] < 0)  && ((*INDS)[k][1] < 0))  ang += PI;
          if (((*INDS)[k][0] <= 0) && ((*INDS)[k][1] >= 0)) ang += 2*PI;
          */
          angInd = ang*NangleBins/(2*PI);
          lineOut[angInd] += (*image)[ir][ic];
          lineOutCount[angInd] += 1;
        }
      }
    }
  }

  // Averaging
  if (NangleBins != -1) {
    for (uint i=0; i<NangleBins; i++) {
      if (lineOutCount[i]) {
        lineOut[i] /= lineOutCount[i];
      }
      else {
        lineOut[i] = 0;
      }
    }
  }

  return lineOut;
}


double imgProc::radProcTool::radialSliceVar(vector< vector<double> >* image,
    int centerR, int centerC, int rad, int shellWidth, bool verbose) {

  // Check if indices have been imported
  if (verbose) cout << "\tChecking for indices " 
                    << shellWidth << " " << rad <<endl;
  checkForIndices(shellWidth, rad);

  // Filter out bins that are off the image or contain nans
  if (verbose) cout << "\tFilling indices of non-NANVAL bins on image\n";
  int ir, ic;
  vector<double> vals;
  vector< vector<int> >* INDS = &allIndices[shellWidth][rad];
  for (uint i=0; i<(*INDS).size(); i++) {
    ir = centerR+(*INDS)[i][0];
    ic = centerC+(*INDS)[i][1];
    if ((ir < (*image).size()) && (ir >= 0)
        && (ic < (*image)[0].size()) && (ic > 0)) {
      if ((*image)[ir][ic] != NANVAL) {
        vals.push_back((*image)[ir][ic]);
      }
    }
  }

  // Get mean
  double mean = getMean(vals);

  // Get Variance
  return std::pow(getSTDev(vals, mean), 2);
}



////////////////////////////////////////////////////////////////////////////
/*
imgProc::centerfnctr::centerfnctr() {
  radProc = NULL;
  verbose = false;
}

double imgProc::centerfnctr::operator() (std::vector<double> vect) {
  if (fxnType == 0)  return imgProc::centerSymXsqr(img, vect[0], vect[1], 
                                8, 1, minRadBin, centShellW, true);
  if (fxnType == 1)  return imgProc::centerOddLeg(img, vect[0], vect[1], 
                                1, minRadBin, centShellW, true); 
  if (fxnType == 2)  {
    if (!radProc) {
      std::cerr << 
          "ERROR: Must specifiy radProc in centerfnctr to use this option!!!\n";
      exit(0);
    }  
    return radProc->radialSliceVar(img, vect[0], vect[1], 
              minRadBin, centShellW, verbose);
  }
  else { 
      cerr << 
          "ERROR: Did not select a correct center finding alogrithm, now exiting!!!\n\n";
      exit(0);
    }
    //N2 return imgProc::centerOddLeg(img, vect[0], vect[1], 1, 70, 40, true);
  } 
} 
*/

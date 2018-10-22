#include "hdf5.h"

using namespace std;


template<typename T>
struct dataVector {

  T* contig1Darray;
  int rank;
  uint* dims;

  dataVector(int inpRank) {
    rank = inpRank;
    contig1Darray = NULL;
    dims = NULL;
  };
  ~dataVector() {
    clear1DCV();
  };

  virtual void initArray(hsize_t* inpDims, int shift) = 0;
  virtual void initArray(hsize_t* inpDims) = 0;
  void init(hsize_t* inpDims, int shift) {
    clear();

    uint size=1;
    dims = new uint[rank]; 	// Change if dataScalar
    for (int ir=0; ir<rank; ir++) {
      dims[ir] = (uint)inpDims[ir+shift]; 
      size *= dims[ir];
    }

    contig1Darray = new T[size];
  };

  uint size(int ir) {
    return dims[ir];
  };

  bool compareDims(hsize_t* compD, int shift) {
    if (!dims) return false;

    bool same = true;
    for (int ir=0; ir<rank; ir++) {
      if (dims[ir] != (uint)compD[ir+shift]) {
	same = false;
	break;
      }
    }

    return same;
  }

  virtual void clear() = 0;
  void clear1DCV() {
    if (!contig1Darray) return;
    delete[] contig1Darray;
    delete[] dims;		// Change if dataScalar
    contig1Darray = NULL;
    dims = NULL;
  };
};


/*
// Scalar
template<typename T>
struct dataScalar: public dataVector<T> {

  T* array;

  dataScalar() : dataVector<T>(1) {
    array = NULL;
  };

  T operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims) {
    if (array) return;
    clear();
    this->dataVector<T>::init(inpDims);
    array = this->dataVector<T>::contig1Darray;
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    array = NULL;
  };
};
*/


// 1D array
template<typename T>
struct data1Dvec: public dataVector<T> {

  T* array;

  data1Dvec() : dataVector<T>(1) {
    array = NULL;
  };

  T operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims, int shift) {
    if (this->dataVector<T>::compareDims(inpDims, shift)) return;
    this->dataVector<T>::init(inpDims, shift);
    array = this->dataVector<T>::contig1Darray;
  };
  void initArray(hsize_t* inpDims) {
    this->initArray(inpDims, 0);
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    array = NULL;
  };
};


template<typename T>
struct data2Dvec: public dataVector<T> {

  T** array;

  data2Dvec() : dataVector<T>(2) {
    array = NULL;
  };

  T* operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims, int shift) {
    if (this->dataVector<T>::compareDims(inpDims, shift)) return;
    this->dataVector<T>::init(inpDims, shift);
    uint* p2dims = this->dims;

    array = new T*[p2dims[0]];
    for (uint i0=0; i0<p2dims[0]; i0++) {
      array[i0] = &(this->dataVector<T>::contig1Darray[i0*p2dims[1]]);
    }
  };
  void initArray(hsize_t* inpDims) {
    this->initArray(inpDims, 0);
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    delete[] array;
    array = NULL;
  };
};


template<typename T>
struct data3Dvec: public dataVector<T> {

  T*** array;

  data3Dvec() : dataVector<T>(3) {
    array = NULL;
  };

  T** operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims, int shift) {
    if (this->dataVector<T>::compareDims(inpDims, shift)) return;
    this->dataVector<T>::init(inpDims, shift);
    uint* p2dims = this->dims;

    array = new T**[p2dims[0]];
    for (int i0=0; i0<p2dims[0]; i0++) {
      array[i0] = new T*[p2dims[1]];
      for (uint i1=0; i1<p2dims[1]; i1++) {
	array[i0][i1] = &(this->dataVector<T>::contig1Darray[i0*(p2dims[1]*p2dims[2]) + i1*(p2dims[2])]);
      } 
    } 
  };
  void initArray(hsize_t* inpDims) {
    this->initArray(inpDims, 0);
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    uint* p2dims = this->dims;
    for (uint i0=0; i0<p2dims[0]; i0++) delete[] array[i0];
    delete[] array;
    array = NULL;
  };
};


template<typename T>
struct data4Dvec: public dataVector<T> {

  T**** array;

  data4Dvec() : dataVector<T>(4) {
    array = NULL;
  };

  T*** operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims, int shift) {
    if (this->dataVector<T>::compareDims(inpDims, shift)) return;
    this->dataVector<T>::init(inpDims, shift);
    uint* p2dims = this->dims;

    array = new T***[p2dims[0]];
    for (uint i0=0; i0<p2dims[0]; i0++) {
      array[i0] = new T***[p2dims[1]];
      for (uint i1=0; i1<p2dims[1]; i1++) {
        array[i0][i1] = new T**[p2dims[2]];
      	for (uint i2=0; i2<p2dims[2]; i2++) {
	  array[i0][i1][i2] = &(this->dataVector<T>::contig1Darray[i0*(p2dims[1]*p2dims[2]*p2dims[3]) + i1*(p2dims[2]*p2dims[3]) + i2*(p2dims[3])]);
	}
      } 
    } 
  };
  void initArray(hsize_t* inpDims) {
    this->initArray(inpDims, 0);
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    uint* p2dims = this->dims;
    for (uint i0=0; i0<p2dims[0]; i0++) {
      for (uint i1=0; i1<p2dims[1]; i1++) {
	delete[] array[i0][i1];
      }
      delete[] array[i0];
    }
    delete[] array;
    array = NULL;
  };
};


template<typename T>
struct data5Dvec: public dataVector<T> {

  T***** array;

  data5Dvec() : dataVector<T>(5) {
    array = NULL;
  };

  T**** operator[](int n) {
    return this->array[n];
  };
 
  void initArray(hsize_t* inpDims, int shift) {
    if (this->dataVector<T>::compareDims(inpDims, shift)) return;
    this->dataVector<T>::init(inpDims, shift);
    array[0][0][0][0] = this->dataVector<T>::contig1Darray;
    uint* p2dims = this->dims;

    array = new T****[p2dims[0]];
    for (uint i0=0; i0<p2dims[0]; i0++) {
      array[i0] = new T****[p2dims[1]];
      for (uint i1=0; i1<p2dims[1]; i1++) {
        array[i0][i1] = new T***[p2dims[2]];
      	for (uint i2=0; i2<p2dims[2]; i2++) {
          array[i0][i1][i2] = new T***[p2dims[3]];
      	  for (uint i3=0; i3<p2dims[3]; i3++) {
	    array[i0][i1][i2][i3] = &(this->dataVector<T>::contig1Darray[i0*(p2dims[1]*p2dims[2]*p2dims[3]*p2dims[4]) + i1*(p2dims[2]*p2dims[3]*p2dims[4]) + i2*(p2dims[3]*p2dims[4]) + i3*(p2dims[4])]);
	  }
	}
      } 
    } 
  };
  void initArray(hsize_t* inpDims) {
    this->initArray(inpDims, 0);
  };

  void clear() {
    if (!array) return;
    this->dataVector<T>::clear1DCV();
    uint* p2dims = this->dims;
    for (uint i0=0; i0<p2dims[0]; i0++) {
      for (uint i1=0; i1<p2dims[1]; i1++) {
        for (uint i2=0; i2<p2dims[1]; i2++) {
 	  delete[] array[i0][i1][i2];
	}
	delete[] array[i0][i1];
      }
      delete[] array[i0];
    }
    delete[] array;
    array = NULL;
  };
};










/*
template<typename T>
struct data2Dvec: public dataVector {

  T** array;
  T* operator[](int n) {
    return this->array[n];
  };

  data2Dvec() {
    rank = 2;
    array = NULL;
  };
  void init(T** inpArray, hsize_t* inpDims) {
    this->dataVector::init(2, inpDims);
    array = inpArray;
  };

  ~data2Dvec(){
    clear();
  };
  void clear() {
    if (! array) return;
    delete[] array[0];
    array = NULL;
  };
};

template<typename T>
struct data3Dvec: public dataVector {

  T*** array;
  T** operator[](int n) {
    return this->array[n];
  };

  data3Dvec() {
    rank = 3;
    array = NULL;
  };
  void init(T*** inpArray, hsize_t* inpDims) {
    this->dataVector::init(3, inpDims);
    array = inpArray;
  };

  ~data3Dvec(){
    clear();
  };
  void clear() {
    if (! array) return;
    delete[] array[0][0];
    array = NULL;
  };
};


template<typename T>
struct data4Dvec: public dataVector {

  T**** array;
  T*** operator[](int n) {
    return this->array[n];
  };

  data4Dvec() {
    rank = 4;
    array = NULL;
  };
  void init(T**** inpArray, hsize_t* inpDims) {
    this->dataVector::init(4, inpDims);
    array = inpArray;
  };

  ~data4Dvec(){
    clear();
  };
  void clear() {
    if (! array) return;
    delete[] array[0][0][0];
    array = NULL;
  };
};


template<typename T>
struct data5Dvec: public dataVector {

  T***** array;
  T**** operator[](int n) {
    return this->array[n];
  };

  data5Dvec() {
    rank = 5;
    array = NULL;
  };
  void init(T***** inpArray, hsize_t* inpDims) {
    this->dataVector::init(5, inpDims);
    array = inpArray;
  };

  ~data5Dvec(){
    clear();
  };
  void clear() {
    if (! array) return;
    delete[] array[0][0][0][0];
    array = NULL;
  };
};
*/

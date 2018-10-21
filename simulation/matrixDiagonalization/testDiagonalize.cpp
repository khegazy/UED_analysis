#include <Sparse>
#include <iostream>
#include <vector>
#include <ctime>
#include <Math/SpecFuncMathMore.h>

using namespace Eigen;
using namespace std;

/*
 * Test case is for stacked pulses
 */

int main() {

  typedef Triplet<double> T;

  clock_t cB, cE;
  double B = 10.5e9;
  double h = 6.62607004e-34;
  int rank = 100;
  int m=0;

  SparseMatrix<double> interH(rank,rank);
  SparseMatrix<double> Ap(rank,rank);
  SparseMatrix<double> V(rank,rank);
  SparseMatrix<double> Q,R,A,Vt;
  Matrix<int,50,50> Pd,Pdi;
  PermutationMatrix<Dynamic, Dynamic, int> P;
  vector<T> coeffs(3*rank);
  vector<T> coeffsI(3*rank);

cout<<"start to fill"<<endl;
  for (int j=0; j<rank; j++) {
    if (j>1) coeffs.push_back(T(j-2,j,sqrt((2*(j-2)+1)*(2*j+1))*(ROOT::Math::wigner_3j((j-2),2,j,-m,0,m)*ROOT::Math::wigner_3j((j-2),2,j,0,0,0)+ROOT::Math::wigner_3j((j-2),0,j,-m,0,m)*ROOT::Math::wigner_3j((j-2),0,j,0,0,0)/3)));
    coeffs.push_back(T(j,j,(2*j+1)*(ROOT::Math::wigner_3j(j,2,j,-m,0,m)*ROOT::Math::wigner_3j(j,2,j,0,0,0)+ROOT::Math::wigner_3j(j,0,j,-m,0,m)*ROOT::Math::wigner_3j(j,0,j,0,0,0)/3)));
    coeffsI.push_back(T(j,j,1));
    if (j<rank-2) coeffs.push_back(T(j+2,j, sqrt((2*(j+2)+1)*(2*j+1))*(ROOT::Math::wigner_3j((j+2),2,j,-m,0,m)*ROOT::Math::wigner_3j((j+2),2,j,0,0,0)+ROOT::Math::wigner_3j((j+2),0,j,-m,0,m)*ROOT::Math::wigner_3j((j+2),0,j,0,0,0)/3)));
  }
  interH.setFromTriplets(coeffs.begin(),coeffs.end());
  interH.makeCompressed();
  V.setFromTriplets(coeffsI.begin(), coeffsI.end());
  V.makeCompressed();

  cout<<interH<<endl;

  //for (int ir=0; ir<rank; ir++) {
  //  for (SparseMatrix<double>::InnerIterator it(interH,ir); it; ++it) cout<<it.row()<<"   "<<it.col()<<"   "<<it.value()<<endl;
  //  cout<<endl;
  //}




  cB=clock();
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > decompQR;
  decompQR.analyzePattern(interH);
  decompQR.factorize(interH); 
  V = decompQR.matrixQ();

  for (int k=0; k<7000; k++) {
    P = decompQR.colsPermutation();
    Q = decompQR.matrixQ();
    R = decompQR.matrixR();
    //Ap = decompQR.matrixR()*(decompQR.colsPermutation()).transpose()*decompQR.matrixQ();
    //Ap = decompQR.matrixR()*(decompQR.colsPermutation()).transpose();
    Ap = R*P.transpose();
    Ap = Ap*Q;

    decompQR.analyzePattern(Ap);
    decompQR.factorize(Ap); 
    Q = decompQR.matrixQ();
    V=V*Q;

    if(k%20) Ap.prune(1e-18);
/*    if (k%10==0) {
      for (int ir=0; ir<rank; ir++) {
    	for (SparseMatrix<double>::InnerIterator it(Ap,ir); it; ++it){
      	  if (abs(it.value()) <1e-10) Ap.coeffRef(it.row(), it.col()) = 0;  
    	}
      }
      for (int ir=0; ir<rank; ir++) {
    	for (SparseMatrix<double>::InnerIterator it(R,ir); it; ++it){
      	  if (abs(it.value()) <1e-10) R.coeffRef(it.row(), it.col()) = 0;  
    	}
      }
      for (int ir=0; ir<rank; ir++) {
    	for (SparseMatrix<double>::InnerIterator it(Q,ir); it; ++it){
      	  if (abs(it.value()) <1e-10) Q.coeffRef(it.row(), it.col()) = 0;  
    	}
      }
//      cout<<Ap<<endl<<R<<endl<<Q<<endl<<endl<<"------------------------------------------------------"<<endl<<endl;  

    }
*/
  } 
  cE=clock();

  for (int ir=0; ir<rank; ir++) {
    for (SparseMatrix<double>::InnerIterator it(Ap,ir); it; ++it){
      if (abs(it.value()) <1e-15) Ap.coeffRef(it.row(), it.col()) = 0;  
    }
  }

  //for (int ir=0; ir<rank; ir++) {
  //  for (SparseMatrix<double>::InnerIterator it(V,ir); it; ++it){
  //	  if (abs(it.value()) <1e-10) V.coeffRef(it.row(), it.col()) = 0;  
  //  }
  //}
  Ap.prune(1e-15);
  cout<<Ap<<endl;
  //cout<<V<<endl;
  //cout<<V.transpose()*interH*V<<endl;

  cout<<"Time to diagonalize: "<<double(cE-cB)/CLOCKS_PER_SEC<<endl;





/*
  SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > decompQR;
  decompQR.analyzePattern(interH);
  decompQR.factorize(interH); 

  P = decompQR.colsPermutation();
  Pd = P.toDenseMatrix();
  Pdi = (P.inverse()).toDenseMatrix(); 
  Q = decompQR.matrixQ();
  R = decompQR.matrixR();

  for (int ir=0; ir<rank; ir++) {
    for (SparseMatrix<double>::InnerIterator it(Q,ir); it; ++it){
      if (abs(it.value()) <1e-10) Q.coeffRef(it.row(), it.col()) = 0;  
    }
  }
  for (int ir=0; ir<rank; ir++) {
    for (SparseMatrix<double>::InnerIterator it(R,ir); it; ++it){
      if (abs(it.value()) <1e-10) R.coeffRef(it.row(), it.col()) = 0;  
    }
  }
  A = Q*R*P.transpose();
  Ap = R*P.transpose();
  Ap = A*Q;
  cout<<Q*(R*P.transpose())<<endl;
  cout<<(R*P.transpose())*Q<<endl;

  cout<<Q<<endl;
  cout<<R<<endl;
*/
return 1;
}

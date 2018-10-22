#ifndef ATOMCLASS_H
#define ATOMCLASS_H

#include <iostream>
#include <complex>
#include <vector>
#include <string>

#include "../../../tools/constants.h"
#include "../../../tools/tools.h"

using namespace std;

///////////////////////////////////////
// Atoms to be used in the molecules //
//  must be changed per expirement   //
///////////////////////////////////////

enum ATOMS {H,C,N,O,I,NUMBER_ATOMS};
static string atomNames[NUMBER_ATOMS] = {"hydrogen","carbon","nitrogen","oxygen","iodine"};


/////////////////////////////////////////////////
// ATOMclass stores information about each atom//
// Every atom used is its own instantiation    //
// of this class.			       //
// 					       //	
// There will generally be a base atom from  //
// which we will orientate the molecule, this  //
// atom will be kept track of in the molecule. //
// Starting from the base atom we count the  //
// bonds extending from it and link the atoms  //
// together by storing the pointers to the     //
// atoms which are bonded to the atom of       //
// interest. For the next atom we will do the  //
// same, however we will omit the bond that    //
// we have already counted by travelling from  //
// the base atom to the atom of interest.    //
// 					       //
// Let's say each atom (a) is bonded to N_a    //
// atoms, then the base atom will store N_a  //
// bonds in its class, and all other atoms     //
// will have (N_a'-1) bonds.                   //
// 					       //
// Atoms will generally be made with new and   //
// stored within molecules, for this reason    //
// they will by default be deleted by their    //
// molecule in MOLECULEclass unless otherwise  //
// changed via deleteAtoms.		       //
/////////////////////////////////////////////////

class ATOMclass {

   public:
	ATOMclass();
	ATOMclass(ATOMS atm);
	ATOMclass(ATOMS atm, Eigen::Vector3d location_in);
	ATOMclass(string index, ATOMS atm);
	ATOMclass(string index, ATOMS atm, Eigen::Vector3d location_in);
	~ATOMclass();

	ATOMS atomType;			// Type of atom
	string index;
        Eigen::Vector3d location;		// Location of atom (X,Y,Z)
	
	void initialize(ATOMS atm);
	void initialize(ATOMS atm, Eigen::Vector3d location_in);
	void initialize(string index, ATOMS atm);
	void initialize(string index, ATOMS atm, Eigen::Vector3d location_in);
	void setLocation(double x, double y, double z);
};

#endif


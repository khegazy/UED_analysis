#ifndef MOLECULECLASS_H
#define MOLECULECLASS_H

#include <iostream>

#include "../../../tools/constants.h"
#include "atomClass.h"

using namespace std;


//////////////////////////////////////////////////////////////
// MOLECULEclass is primarily a means of storing the atoms  //
// that belong to a single molecule.                        //
//                                                          //
// Molecules are the primary entity the code will interact  //
// with, the atoms are thus managed by the molecule. For    //
// this reason the molecule will delete the atoms in the    //
// destructor by default. This can be changed by changing   //
// deleteAtoms to false.                                    //
//                                                          //
// baseAtom is the primary atoms from which our perception  //
// of the molecule is orientated, it is the starting point  //
// from where you can traverse the molecule to other atoms. //
// atoms is simply an array of pointers to all the atoms in //
// the molecule and Natoms is the number of atoms in the    //
// molecule.                                                //
//                                                          //
// One will generally interact with the class itself rather //
// than a pointer when first making the molecules from the  //
// MC. Afterwards pointers will be used when accessing from //
// root files.                                              //
//////////////////////////////////////////////////////////////

class MOLECULEclass {

   public:
	MOLECULEclass();
	MOLECULEclass(bool deleteAtoms_in);
	~MOLECULEclass();

	vector<ATOMclass*> atoms;

	void initialize(bool deleteAtoms_in);
	void addAtom(ATOMclass* atom);

   private:
	bool deleteAtoms;
};

#endif

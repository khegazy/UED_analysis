#include "atomClass.h"

using namespace std;


ATOMclass::ATOMclass() {}


ATOMclass::ATOMclass(ATOMS atm) {

  initialize(atm);
}


ATOMclass::ATOMclass(ATOMS atm, Eigen::Vector3d location_in) {

  initialize(atm, location_in);
}


ATOMclass::ATOMclass(string index_in, ATOMS atm) {

  initialize(index_in, atm);
}


ATOMclass::ATOMclass(string index_in, ATOMS atm, Eigen::Vector3d location_in) {

  initialize(index_in, atm, location_in);
}


ATOMclass::~ATOMclass() {
}


void ATOMclass::initialize(ATOMS atm) {

  initialize("NULL", atm);
}


void ATOMclass::initialize(ATOMS atm, Eigen::Vector3d location_in) {

  initialize("NULL", atm);
  location = location_in;
}


void ATOMclass::initialize(string index_in, ATOMS atm) {

  index = index_in;
  atomType = atm;
}


void ATOMclass::initialize(string index_in, ATOMS atm, Eigen::Vector3d location_in) {

  initialize(index_in, atm);
  location = location_in;
}


void ATOMclass::setLocation(double x, double y, double z) {
 
  location(0) = x;
  location(1) = y;
  location(2) = z;
}


#include "Photon.h"
#include <sstream>

Photon::Photon(double e, double tx, double ty, double x, double y,double z, double p) {
	_ephot=e;
	_thetax=tx;
	_thetay=ty;
	_xemis = x;
	_yemis = y;
	_zemis=z;
	_poids=p;
}

string Photon::outputFlow() const{
 ostringstream sortie;

//   sortie <<  " ehpot: " <<  _ephot <<  " thetax: " <<  _thetax << " thetay: " <<  _thetay <<  " zemis: " << _zemis <<  " poids: " <<  _poids << endl;
//  sortie <<  "  " <<  _ephot <<  "  " <<  _thetax << "  " <<  _thetay <<  "  " << _zemis <<  "  " <<  _poids << endl;
 return sortie.str();

}


void Photon::print() const{
  cout << outputFlow();
}

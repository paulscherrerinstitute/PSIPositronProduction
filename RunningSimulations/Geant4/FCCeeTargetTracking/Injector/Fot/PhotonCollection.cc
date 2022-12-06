#include "PhotonCollection.h"
#include <fstream>
#include "GlobalConstants.h"
#include <iostream>



void PhotonCollection::saveOnFile(string filename)
{
	ofstream output_file;
output_file.open(filename.c_str());
  list<Photon>::const_iterator it;
	int numero = 1;
	
	  for (it = _photon.begin(); it != _photon.end(); it++)
    {
      //		output_file << numero << it->print() ;
      output_file << numero << it->outputFlow() ;

      numero++;
    }
	output_file.close();
	_currentIterator = _photon.begin();
}


const Photon* PhotonCollection::getNextPhoton(){
  if(_currentIterator != _photon.end()){
    const Photon* p = &(*_currentIterator);
    _currentIterator++;
		
    return p;
  }
  else{
    return NULL;
  }
}




void PhotonCollection::printPhotons() const
{
  list<Photon>::const_iterator it;
  int compteur = 1;
  for (it = _photon.begin(); it != _photon.end(); it++)
    {
      std::cout <<"photon n° " << compteur;
      it->print();
      //      it->print(); 
      compteur++;
    }
}

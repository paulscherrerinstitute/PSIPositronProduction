#ifndef _PHOTONCOLLECTION_H
#define _PHOTONCOLLECTION_H

#include "Photon.h"
#include <list>

/**
 * \file PhotonCollection.h
 * \brief The PhotonCollection class provides a useful way to create collection of photons
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

using namespace std;


/** \class PhotonCollection
 *
 *  The PhotonCollection class is the object that will contain all photons emitted by a particle in the crystal
 */
class PhotonCollection {

  list<Photon> _photon;
  list<Photon>::const_iterator _currentIterator;

 public:

	
  /**
   *  \brief Constructor
   *
   *  Default constructor of a collection of photons
   *
   */

 PhotonCollection() { _currentIterator = _photon.begin(); }
  
  const list<Photon>& getPhotonList() const
  {
    return _photon;
  }

  //list<Photon>& getPhotonList()
  //  {
  //    return _photon;
  //}


void clear() 
{
  _photon.clear();
  _currentIterator = _photon.begin();
}

	
/**
 *  \brief return the lenght of the collection of photons
 *
 *  return the number of photons in the collection
 *
 * \return int: return the number of photons
 */

int lenght() const
{
return _photon.size();
}
	
	
/**
 *  \brief fill the list
 *
 *  Add a new photon emitted in the list _photon
 *
 *  \param double :
 *  \param double : X photon emission angle
 *  \param double : Y photon emission angle
 *  \param double : Z position of the photon emission
 *  \param double : weight of the photon
 */

	
void fill(double e, double tx, double ty, double x, double y, double z, double p){
  _photon.push_back(Photon(e,tx,ty,x,y,z,p));
}

	
/**
 *  \brief Register all photons of the collection in a file
 *
 *  Fill a file with data from all photons of the list _photons
 *
 *  \param string : the name of registration file
 */
void saveOnFile(string filename);
	
	
/**
 *  \brief return the next photon
 *
 *  return a pointer to the next photon in the collection
 *
 * \return Photon*: return the next photon
 */
const Photon* getNextPhoton();

	
/**
 *  \brief return the lenght of the collection of photons
 *
 *  return the number of photons in the collection
 *
 * \return int: return the number of photons
 */

int getNbPhotons() const
{
  return _photon.size();
}


const Photon& getLastPhoton() const
{
  return _photon.back();
}	
/**
 *  \brief Print all photons which are in the collection
 *
 *  Edit data concerning all photons of the collection 
 *
 */
void printPhotons() const;

};

#endif

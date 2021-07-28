#ifndef _FOT_H
#define _FOT_H

#include <iostream>
#include <math.h>

#include "Particle.h"
//#include "Photon.h"
#include "Snake.h"
#include "Crystal.h"
#include "Lindhard.h"
#include "ParticleCollection.h"
#include "RunParameters.h"
//#include "Bremsstrahlung.h"
#include "statistiques.h"


/**
 * \file Fot.h
 * \brief The Fot class provides a useful way to implement the phenomenon of channeling
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

using namespace std;


/** \class Fot
 *
 *  The Fot class is the object that will contain all data concerning the phenomenon of channeling 
 */
class Fot {

  PhotonCollection _photons;
  Snake* _snak;
  ParticleInCrystal* _partCrys;
  Lindhard _lind;
  statistiques* _stat;
  const RunParameters& _runPar;

  double _etmax, _vtmax, _zexit;


  int _nevnt;


 public:

  /**
   *  \brief Constructor
   *
   *  Construct the Fot object from data 
   *
   *  \param RunParameters& : parameters concerning the crystal and the Kuma photon
   *
   */
  Fot(RunParameters& rp);
	
	
  /**
   *  \brief Destructor
   *
   *  Destroy the Fot object
   *
   *
   */
  ~Fot()
    {
      if ( _snak != NULL ) delete _snak;
      if ( _partCrys != NULL ) delete _partCrys;
      if ( _stat != NULL ) delete _stat;
    }


	
  /**
   *  \brief Main method where the algorithm of channeling is implemented
   *
   *  Main method where the algorithm of channeling is implemented
   *
   */
  //   void make();
  void makeKumakhov(ParticleCollection& partColl);

  const ParticleInCrystal& makeSingleParticleKumakhov(const Particle* part);

  void poirot();


  const PhotonCollection& getPhotonCollection() const
  {
    return _photons;
  }

  void finir() 
  {
    cout << " end of the bunch " << endl;
    cout << endl;
#ifdef DO_STATS
    _stat->tsatis();
#endif
  }
};

#endif

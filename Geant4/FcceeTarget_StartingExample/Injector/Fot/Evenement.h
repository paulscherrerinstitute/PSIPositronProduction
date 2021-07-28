#ifndef _EVENEMENT_H
#define _EVENEMENT_H

#include "Particle.h"
#include "Photon.h"
#include "Snake.h"
#include "Bremsstrahlung.h"
#include "Crystal.h"
#include "Lindhard.h"

#include <iostream>
#include <cmath>
#include "statistiques.h"

class Evenement 
{
  // const RunParameters* _runPar;
  ParticleInCrystal* _partCrys;
  Snake* _snak;
  BremsStrahlung* _bremse;
  PhotonCollection* _photons;

  double _Zexit;

  double _vtmax;
  double _etmax;
  
  double _pt2max;

  double _vt, _ecin, _dtinv1, _dtinv2, _gamma;

  bool _restartSnake;

  double _zj;

 public:

  Evenement( Snake* snak, BremsStrahlung* brms, PhotonCollection* photons, double etmax, double vtmax, double Zexit);

  bool reInitSnake();
  bool makeStep();

};


#endif

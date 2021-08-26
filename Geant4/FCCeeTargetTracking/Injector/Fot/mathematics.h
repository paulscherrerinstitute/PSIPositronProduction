#ifndef _MATHEMATICS_H
#define _MATHEMATICS_H
#include "rndm.h"
#include "GlobalConstants.h"

#include <iostream>
#include <math.h>

/**
 * \file Particle.h
 * \brief The mathematics class provides a useful way to do mathematical calculation
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

using namespace std;


/** \class mathematics
 *
 *  The mathematics class is the object that contain two mathematical functions
 */

class mathematics  {

  static  RNDM _alea;

 public:

  static void set_seed(int seed){
    _alea.set_seed(seed);
  }

  static double rndm()
  {
    return _alea.rndm();
  } 

  static double gauss( double xm, double s )
  {
    return _alea.gaussian(xm, s);
  }
  
  static  RNDM* getRandomGenerator() 
  {
    return &_alea;
  }	

};

#endif

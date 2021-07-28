#ifndef _LINDHARD_H
#define _LINDHARD_H

#include <iostream>


// C   DCOLL,DCOLTA: inverse du libre parcours moyen de diff. incoherente.



/**
 * \file Particle.h
 * \brief The Lindhard class provides a useful way to make a linear interpolation of potentials to calculate the Lindhard force
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

using namespace std;



/** \class  Lindhard
 *
 *  The Lindhard class is the object that will calculate the Lindhard force
 */
typedef struct  {
    double _epot;
	
	double _f;

    double _fx;

    double _fy;

    double _dcoll;

} Lindhard;
#endif

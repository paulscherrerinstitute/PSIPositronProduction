#ifndef _BREMSSTRAHLUNG_H
#define _BREMSSTRAHLUNG_H


/**
 * \file BremsStrahlung.h
 * \brief The BremsStrahlung class provides a useful way to control the multiple scattering and the emission of photons of bremsstrahlung
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

#include <vector>
#include <math.h>
using namespace std;

#include "Particle.h"
#include "PhotonCollection.h"
#include "GlobalConstants.h"
#include "statistiques.h"


/** \class  BremsStrahlung
 *
 *  The BremsStrahlung class is the object that will contain all data concerning the multiple scattering
 */
class BremsStrahlung 
{
  //  RNDM& _hasard;
  PhotonCollection& _photons;

  ParticleInCrystal* _part;

  double _phomin;
  double _poimin;
  double _xmilog;

  double _ancoll;

  statistiques* _stat;

 public :

	/**
     *  \brief Constructor
     *
     *  Construct the photon of bremsstrahlung from data 
     *
     *  \param PhotonCollection& : the collection of the photons emitted
     *
     */
  BremsStrahlung(PhotonCollection& fcol, double phomin, double poimin, statistiques* stat);

	/**
     *  \brief 
     *
     *
     *
	 *  \param ParticleInCrystal& : the particle in the crystal
	 *  \param double : random number
     *
     */
  void reInit(ParticleInCrystal* p);

	
	/**
     *  \brief 
     *
     *
     *
	 *  \param double : weight of the crystal
	 *  \param double : time step
	 *  \param double : density of inelastic collision
     *
	 *  \return bool :
     */
  bool multipleDiffusion( double dt,  double Dcoll);

 private: 

	
	/**
     *  \brief 
     *
     *
     *
	 *  \param double : weight of the crystal
	 *  \param double : step of proper time
	 *  \param double : local cut-off
	 *  \param int : number of collision
	 *  \param double& : 
	 *  \param double& : 
	 *  \param double& : 
     *
     */
  void difmul(double DTO,double UC,int ncol, double& Q, double& QX, double& QY);

	
	/**
     *  \brief 
     *
     *
     *
	 *  \param double : step of proper time
	 *  \param double : local cut-off
	 *  \param double : 
	 *  \param double : 
	 *  \param double : 
	 *  \param double& : 
	 *  \param double& : 
	 *  \param double& : 
	 *  \param double& : 
     *
     */
  void abrems(double DTO, double UC, double qq, double qx, double qy, double& AA, double& ZETA, double& GATX, double& GATY);



};


#endif

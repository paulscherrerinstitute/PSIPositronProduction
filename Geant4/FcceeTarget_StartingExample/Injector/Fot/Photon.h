#ifndef _PHOTON_H
#define _PHOTON_H

#include <iostream>



/**
 * \file Photon.h
 * \brief The Photon class provides a useful way to create a photon with a predifined specifications
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

using namespace std;


/** \class  Photon
 *
 *  The Photon class is the object that will contain all data concerning the photon
 */
class Photon {

  double _ephot;

    double _thetax;

    double _thetay;


    double _xemis, _yemis;
    double _zemis;

    double _poids;

  public:
	
	/**
     *  \brief Constructor
     *
     *  Construct the photon from data 
     *
     *  \param double :
	 *  \param double : X photon emission angle
	 *  \param double : Y photon emission angle
	 *  \param double : Z position of the photon emission
	 *  \param double : weight of the photon
     *
     */
    Photon(double e, double tx, double ty, double x, double y,double z, double p);	
	


    double getEnergy() const 
    {
      return _ephot;
    }


    double getXemis() const 
    {
      return _xemis;
    }

    double getYemis() const 
    {
      return _yemis;
    }

    double getZemis() const 
    {
      return _zemis;
    }

    double getThetax() const
    {
      return _thetax;
    }

    double getThetay() const
    {
      return _thetay;
    }

    double getWeight() const
    {
      return _poids;
    }


	/**
     *  \brief Stock data of the photon
     *
     *  Stock data of the photon in a variable
     *
     * \return bool: return data concerning the photon
     */
	string  outputFlow() const;

	
	/**
     *  \brief Print the photon
     *
     *  Edit data concerning the photon 
	 *
     */
	void print() const;
};
#endif

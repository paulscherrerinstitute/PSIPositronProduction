#ifndef _CRYSTAL_H
#define _CRYSTAL_H



// c p3 ------------------------------------------------
// C    UTHERM: amplitude de vibration du reseau cristalin (en X ou en Y)
// C  F,FTABLE: force de Lindhard en m(elec)/angstrom
// C  Les variables VCOR,FCOR sont specifiques de <110> axial type diamant
// C   DCOLL,DCOLTA: inverse du libre parcours moyen de diff. incoherente.
// C     1/BCOH: recul minimum au cours d'une telle diffusion.


/**
 * \file Crystal.h
 * \brief The Crystal class provides a useful way to create a crystal with a predifined structure
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */



#include <vector>
#include "Lindhard.h"



using namespace std;

  typedef struct 
  {
    double V, F, deltaV;
  } Vdecanal ;


/** \class Crystal
 *
 *  The Crystal class is the object that will contain all data concerning the crystal in which the particle will move.
 */
class Crystal {

    double _utherm;
    double _amaille;
    double _xlib, _ylib; // "decanalisation" point

    double _distatom;

    double _bcoh;

    double _qqcoh;

    double _rbohr2;

    double _zreco;

    double _barrier;

    vector<double> _potable;

    vector<double> _ftable;

    vector<double> _dcolta;

    vector<double> _cristal;


    //    Lindhard _lind;

  public:
	
	
    /**
     *  \brief Constructor
     *
     *  Default constructor of the crystal
     *
     */
    Crystal() {}

	
	
    /**
     *  \brief Constructor
     *
     *  Constructor of a crystal from an atom and an axe
     *
     *  \param atom : string 
     *  \param axe : index of the crital axis
     *
     */
    Crystal(string atom, int axe);

	
	/**
     *  \brief Destructor
     *
     *  Destroy the crystal 
     *
     *
     */
    ~Crystal();

	
	
	/**
     *  \brief return the amplitude of vibration
     *
     *  Allows access to the amplitude of vibration of the crystal lattice (X or Y) 
     *
     * \return double: return the amplitude of vibration _utherm
     */
       
double getUtherm() const{
	return _utherm;
}
	
	
	/**
     *  \brief return the maximum impact parameter
     *
     * Allows access to the maximum impact parameter for incoherent scattering on nucleus
     *
     * \return double: return the maximum impact parameter _bcoh
     */

double getBcoh() const{
	return _bcoh;
}

	
	/**
     *  \brief
     *
     * 
     *
     * \return double:
     */

double getQqcoh() const{
	return _qqcoh;
}
	
	
	/**
     *  \brief
     *
     * 
     *
     * \return double:
     */

double getRbohr2() const{
	return _rbohr2;
}
	
	
	/**
     *  \brief return the potentiel 
     *
     * Allows access to the potentiel of the crystal 
     *
     * \return vector<double>: return 
     */

const vector<double>& getPotable() const{
	return _potable;
}
	
	
	/**
     *  \brief return the Lindhard force 
     *
     * Allows access to the Lindhard force of the crystal 
     *
     * \return vector<double>: return 
     */

const vector<double>& getFtable()  const{
	return _ftable;
}
	
	
	/**
     *  \brief return the inverse average free path
     *
     * Allows access to the inverse average free path of incoherent scattering 
     *
     * \return vector<double>: return 
     */
const vector<double>& getDcolta() const{
	return _dcolta;
}
	
	
	/**
     *  \brief return axis
     *
     * Allows access to the 16 axis of the crystal
     *
     * \return vector<double>: return 
     */

const vector<double>& getCristal() const{
	return _cristal;
}

	
	/**
     *  \brief print geometric parameters
     *
     * Print geometric parametres of the crystal.
     *
     */
	void printGeometricParameters() const;
	
	
	/**
     *  \brief calculate fields
     *
     * Method to calculate the distance to the nearest axis, the potential energy and the component of x and y of the Lindhard force  
     *
	 * \param double: X position of the particle
     * \param double: Y position of the particle
	 * \param double: charge of the particle
	 *
     * \return Lindhard: return fields
     */
	Lindhard calculationOfFields(double xx, double yy, double rayon, double charge) const;

 private :

	void setTungstene(int axe);
	void struc111(string typeOfNetwork, double zatom);
};
#endif

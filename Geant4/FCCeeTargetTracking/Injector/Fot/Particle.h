#ifndef _PARTICLE_H
#define _PARTICLE_H


#include "Crystal.h"
#include "GlobalConstants.h"
#include "Lindhard.h"
#include <iostream>

#include <math.h>

/**
 * \file Particle.h
 * \brief The Particle class provides a useful way to create a particle with a predifined specifications
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */



/** \class  Particle
 *
 *  The Particle class is the object that will contain all data concerning the particle
 */
class Particle {
 protected :

  double _charge;
  double _xPosition;
  double _yPosition;
  double _rayon;
  double _zPosition;
  double _px;
  double _py;
  double _gamma;

  bool  _isPosit;

 public:

  Particle() {}
  
  /**
   *  \brief Copy constructor
   *
   *  Construct the particle from data of another particle
   *
   *  \param Particle& : particle for the copy
   *
   */
  Particle(const Particle& part);
	
  /**
   *  \brief Constructor
   *
   *  Construct the particle from data 
   *
   *  \param double : charge of the particle
   *  \param double : X position of the particle
   *  \param double : Y position of the particle
   *  \param double : Z position of the particle
   *  \param double : X moment of the particle
   *  \param double : Y moment of the particle
   *  \param double : photon energy
   *
   */
  Particle(double ch, double x, double y, double z, double px, double py, double g);
	

  /**
   *  \brief return the charge of the particle
   *
   *  Allows access to the charge of the particle
   *
   * \return double: return the charge _charge
   */
  inline double getCharge() const
  {
    return _charge; 
  }


  double getXPosition() const{
    return _xPosition;
  }

  double getYPosition() const{
    return _yPosition;
  }

  double getRadius() const{
    return _rayon;
  }

  /*  void getTransversePosition(double& xx, double& yy, double& rayon) const */
  /*  { */
  /*    xx = _xPosition; */
  /*    yy = _yPosition; */
  /*    rayon = _rayon; */
  /*  } */

  double getZPosition() const{
    return _zPosition;
  }

  double getPx() const{
    return _px;
  }

  double getPy() const{
    return _py;
  }


  double getGamma() const{
    return _gamma;
  }

  bool isPositron() const
  {
    return _isPosit;
  }

  void printPoirot() const
  {
    cout << " X= " << _xPosition << " Y= " << _yPosition << " Z= " << _zPosition << " PX= " << _px << " PY= " << _py << " GAMMA= " <<  _gamma << endl;
  }
	
	
  /**
   *  \brief return the X position of the particle
   *
   *  Allows access to the X position of the particle
   *
   * \return double: return the X position _xPosition
   */
  //	double getXPosition() const;
	
	
  /**
   *  \brief return the Y position of the particle
   *
   *  Allows access to the Y position of the particle
   *
   * \return double: return the Y position _yPosition
   */
  //	double getYPosition() const;
	
	
  /**
   *  \brief return the Z position of the particle
   *
   *  Allows access to the Z position of the particle
   *
   * \return double: return the Z position _zPosition
   */
  //	double getZPosition() const;
	
	
	
  /**
   *  \brief return the X moment of the particle
   *
   *  Allows access to the X moment of the particle
   *
   * \return double: return the X moment _px
   */
  //	double getPx() const;
	
	
  /**
   *  \brief return the Y moment of the particle
   *
   *  Allows access to the Y moment of the particle
   *
   * \return double: return the Y moment _py
   */
  //	double getPy() const;

	
  /**
   *  \brief return the photon energy
   *
   *  Allows access to the photon energy
   *
   * \return double: return the photon energy _gamma
   */
  //	double getGamma() const;

	
  /**
   *  \brief return the module of angular moment
   *
   *  Calculate the module of angular moment of the particle
   *
   * \return double: return the module of angular moment
   */
  //	double moduleOfAngularMomentum() const;

	
	
  /**
   *  \brief return the sum of squared moment
   *
   *  Calculate the sum of squared moment
   *
   * \return double: return the sum of squared moment
   */
  //	double Ptsquare() const;

	
  /**
   *  \brief Check the charge of the particle
   *
   *  Check whether the particle is a positron or not
   *
   * \return bool: return true if the particle is a positron and false otherwise
   */
  //	bool isPositron() const;
	

	
  /**
   *  \brief move the particle
   *
   *  Move the particle in the crystal using a time step
   *
   * \param double: time step
   */
  //   void move(double dt);
	
	
	
  /**
   *  \brief acceleration of the particle
   *
   *  Calculate transverse acceleration of the particle
   *
   * \param double: time step
   * \param double: X Lindhard force
   * \param double: Y Lindhard force
   */
  //	void acceleration(double dt,double fx, double fy);

	
  /**
   *  \brief change the values of _xPosition and _yPosition
   *
   *  Replace _xPosition and _yPosition by theirs new values xnew and ynew  
   *
   * \param double: new value xnew of _xPosition 
   * \param double: new value ynew of _yPosition
   */
  //	void derapeXY( double xnew, double ynew);

	
  /**
   *  \brief change the value of _zPosition
   *
   *  Replace _zPosition by its new values z  
   *
   * \param double: new value z of _zPosition
   */
  //	void setZ(double z);

	
  /**
   *  \brief change the value of energy moment
   *
   *  Replace _px, _py and _gamma by theirs new values px, py and gamma  
   *
   * \param double: new value px of the X moment of the particle
   * \param double: new value py of the Y moment of the particle
   * \param double: new value gamma of the photon energy
   */
  //	void setEnergyMomentum( double px, double py, double gamma);

	
  /**
   *  \brief increment the transverse moment
   *
   *  Increment the X and Y moments of the particles of dpx and dpy   
   *
   * \param double: X step for X moment
   * \param double: Y step for Y moment
   */
  //	void incrementTransverseMomentum( double dpx, double dpy);

	
  /**
   *  \brief Multiply the photon energy
   *
   *  Multiply the photon energy by factor  
   *
   * \param double: the factor of multiplication
   */
  //	void multiplyGamma(double factor);



  double moduleOfAngularMomentum() const
  {

    return fabs(_xPosition*_py - _yPosition*_px);
  }

  double Ptsquare() const
  {
    return  _px*_px + _py*_py;
  }



  //déplace la particule
  void move(double dt) {
    _xPosition=_xPosition+_px*dt/_gamma;
    _yPosition=_yPosition+_py*dt/_gamma;
    _zPosition=_zPosition+dt;
    _rayon = hypot(_xPosition, _yPosition);
  }


  //acceleration transverse de la particule
  void acceleration(double dt,double fx, double fy) {
    _px=_px+fx*dt;
    _py=_py+fy*dt;
  }

  void setZ(double z)
  {
    _zPosition = z;
  }

  void derapeXY( double xnew, double ynew)
  {
    _xPosition = xnew;
    _yPosition = ynew;
    _rayon = hypot(_xPosition, _yPosition);
  }

  void setEnergyMomentum( double px, double py, double gamma) 
  {
    _px = px;
    _py = py;
    _gamma = gamma;
  } 

  void incrementTransverseMomentum( double dpx, double dpy)
  {
    _px += dpx;
    _py += dpy;

  }

  void multiplyGamma(double factor)
  {
    _gamma *= factor;
  }


	
  /**
   *  \brief Print the particle
   *
   *  Edit data concerning the particle  
   *
   */
  void printParticle() const;


};



class ParticleInCrystal : public Particle
{

 private :

  bool _iemis;

  bool _ientree;

  bool _iexit;

  // c p5  -  X0CELL,Y0CELL: position du centre de la cellule

  double _x0cell;
	
  double _y0cell;

  int _nsaut;

	
  const Crystal & _crys;


  //bigJump corespond a la variable miggg de fortran
  bool _bigJump;

 public :

	
	
  /**
   *  \brief Constructor
   *
   *  Construct the particle from a crystal 
   *
   *  \param Crystal& : the crystal where the particle move
   *
   */
  // ParticleInCrystal(Crystal& crys);

	
  /**
   *  \brief Constructor
   *
   *  Construct the particle from a crystal and data of another particle 
   *
   *  \param Crystal& : the crystal where the particle move
   *  \param Particle& : the particle for the copy
   *
   */
  ParticleInCrystal(const Crystal& crys, const Particle& part);
  ParticleInCrystal(const Crystal& crys);
 ParticleInCrystal(const Crystal& crys, double ch, double x, double y, double z, double px, double py, double g) : Particle( ch,  x,  y,  z,  px, py, g), _crys(crys) {}
  
  /**
   *  \brief Constructor
   *
   *  Construct the particle from a crystal and data 
   *
   *  \param Crystal& : the crystal where the particle move
   *  \param double : charge of the particle
   *  \param double : X position of the particle
   *  \param double : Y position of the particle
   *  \param double : Z position of the particle
   *  \param double : X moment of the particle
   *  \param double : Y moment of the particle
   *  \param double : photon energy
   *
   */
	
  /**
   *  \brief return the crystal of the particle
   *
   *  Allows access to the crystal in which the particle is
   *
   * \return Crystal: return the crystal _c
   */

  const Crystal& getCrystal() const {
    return _crys;
  }


  /**
   *  \brief Check whether the particle emitted
   *
   *  Check whether the particle emitted a photon in the crystal
   *
   * \return bool: return true is the particle emitted a photon and false otherwise
   */




  bool getIemis() const{
    return _iemis;
  }

	
  /**
   *  \brief Check whether the particle entered in the crystal
   *
   *  Check whether the particle entered in the crystal
   *
   * \return bool: return true if the particle entered in the crystal and false otherwise
   */

  bool getIentree() const
  {
    return _ientree;
  }

	
  /**
   *  \brief Check whether the particle exited
   *
   *  Check whether the particle exited of the crystal
   *
   * \return bool: return true if the particle exited of the crystal and false otherwise
   */

  bool getIexit() const
  {
    return _iexit;
  }

	
  /**
   *  \brief Check whether the particle made a big jump
   *
   *  Check whether the particle made a big jump in the crystal
   *
   * \return bool: return true if the particle made a big jump and false otherwise
   */

  bool bigJump() const{
    return _bigJump;
  }

	
  /**
   *  \brief change the value of _iemis
   *
   *  Replace _iemis by its new values b  
   *
   * \param double: new value b of _iemis
   */

  void setIemis(bool b) {
    _iemis=b;
  }

	
  /**
   *  \brief change the value of _ientree
   *
   *  Replace _ientree by its new values b  
   *
   * \param double: new value b of _ientree
   */

  void setIentree(bool b)
  {
    _ientree = b;
  }
	
	
  /**
   *  \brief change the value of _iexit
   *
   *  Replace _iexit by its new values b  
   *
   * \param double: new value b of _iexit
   */

  void setIexit(bool b)
  {
    _iexit = b;
  }
	
	
  /**
   *  \brief change the value of _bigJump
   *
   *  Replace _bigJump by its new values b  
   *
   * \param double: new value b of _bigJump
   */

  void setBigJump(bool b) {
    _bigJump=b;
  }
	
       

	
  /**
   *  \brief return the X center position of the cell 
   *
   *  Allows access to the X center position of the cell in which the particle is
   *
   * \return double: return the X center position _x0cell
   */

  double getX0cell() const{
    return _x0cell;
  }
	
	
  /**
   *  \brief return the Y center position of the cell 
   *
   *  Allows access to the Y center position of the cell in which the particle is
   *
   * \return double: return the Y center position _y0cell
   */

  double getY0cell() const{
    return _y0cell;
  }
	
	
  /**
   *  \brief change the value of _x0cell
   *
   *  Replace _x0cell by its new values x0  
   *
   * \param double: new value x0 of _x0cell
   */

  void setX0cell(double x0){
    _x0cell=x0;
  }

	
	
  /**
   *  \brief change the value of _y0cell
   *
   *  Replace _y0cell by its new values y0  
   *
   * \param double: new value y0 of _y0cell
   */

  void setY0cell(double y0){
    _y0cell=y0;
  }

	
  /**
   *  \brief Check whether the particle migrated
   *
   *  Check whether the particle migrated to another cell in the crystal
   *
   * \return bool: return true if the particle migrated to another cell and false otherwise
   */
  bool migration();
	
	
  /**
   *  \brief change the values of particle and cell coordinates
   *
   *  Replace _xPosition, _yPosition, _x0cell and _y0cell by theirs new values x, y, xCelle and yCelle  
   *
   * \param double: new value x of _xPosition
   * \param double: new value y of _yPosition
   * \param double: new value xCelle of _x0cell
   * \param double: new value yCelle of _y0cell
   */
  void reIniCoorAndCell(double x, double y, double xCelle, double yCelle);




  Lindhard calculationOfFields() const
  {
    return _crys.calculationOfFields(_xPosition,_yPosition,_rayon,_charge);
  }


  void printPoirot() const
  {

    cout << " PARTICLE : " << endl;

    Particle::printPoirot();

    cout << " IENTREE= " << _ientree << " LEXIT= " << _iexit << " NSAUT= " << _nsaut << endl;

  }


 private: 


	
  /**
   *  \brief update coordinates of the particle and its cell
   *
   *  Update coordinates of the particle and its cell when the particle migrated
   *
   * \param double:
   * \param double:
   * \param double:
   */
  void updateCell(double h1, double h2, double h3);
	
};

#endif

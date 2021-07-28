#ifndef _SNAKE_H
#define _SNAKE_H


/**
 * \file Snake.h
 * \brief The Snake class provides a useful way to model the trajectory of the particle in the crystal and manage emissions of photons.
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */


#include <vector>
#include <math.h>
using namespace std;

#include "Particle.h"
#include "RunParameters.h"
#include "Crystal.h"
#include "PhotonCollection.h"
#include "statistiques.h"




/** \class  Snake
 *
 *  The Snake class is the object that will contain all data concerning the trajectory of the particle in the crystal
 */
class Snake {

  typedef struct 
  {
    double xl, yl, zl, px, py, dtau, ucl, chic, dp, gamdv2;
  } SnakeStep ;


  PhotonCollection& _photons;
  const RunParameters & _runPar;
  ParticleInCrystal & _part;

  statistiques* _stat;

  vector<SnakeStep> _xyz;
  int _snakeBegin;

  double _ancoll;

  double _ra;

  int _lpas;

  int _la;

  int _l1;

  int _lb;

  int _l2;

  int _lmax;
    
  ////int _nkuma;
	
	
  double _dto;
	
  double _uc;

  bool _partIsPositron;
  bool _mature;
	
  bool _jzaug;
  bool _eloigne;

  // parametres pour l'integration kumakhov

  double _ximin, _xisup, _usup, _umin, _xmilog;

  // statistiques 

  /**
     double _fotk, _wk;
     double _promax;
     int _mcoll;
     double _favort, _wavort;
     vector<int> _nproba;

     vector< vector<double> > _efficacite;
	
	
  **/

 public:

	
  /**
   *  \brief Constructor
   *
   *  Construct the snake from data 
   *
   *  \param RNDM& : random number
   *  \param RunParameters& : parameters concerning the crystal and the Kuma photon
   *  \param ParticleInCrystal& : the particle in the crystal
   *  \param PhotonCollection& : the collection of photons emitted
   *
   */
  Snake(const RunParameters & rp, ParticleInCrystal & p, PhotonCollection& fcol, statistiques* stat);

  /**
   *  \brief Check whether the snake is mature
   *
   *  Check whether the distance to the nearest axis is maximal
   *
   * \return bool: return true if the snake is mature and false otherwise
   */

  inline bool isMature() const{
    return _mature;
  }

	
  /**
   *  \brief Initialize the snake 
   *
   *  Initializes the snake by removing all elements of vectors and by statistical variables to 0
   *
   */
  void initSnake();
	
	
  /**
   *  \brief Increment the step _lpas
   *
   *  Increase _lpas by 1 and check whether _lpas = _lmax. In this case we can say that the snake is mature
   *
   */
  //  void incrementLpas();
	
  //  void reinitL2();

	
	
  /**
   *  \brief return the space step 
   *
   *  Allows access to the space step
   *
   * \return int : return the space step _lpas
   */



  inline  ParticleInCrystal& getParticle() const
  {
    return _part;
  }

  inline  ParticleInCrystal* getParticlePtr() const
  {
    return &_part;
  }


	
	
	
  /**
   *  \brief  
   *
   *  Allows access to the 
   *
   * \return double:
   */

  inline double getXiMin() const
  {
    return _ximin;
  }

	
  /**
   *  \brief Update vectors of proper time step and local cut-off
   *
   *  Add the new proper time step, the new cut-off in vectors _dtau and _ucl. They are calculated at the step _lpas  
   *
   * \param double: time step
   * \param double: Lindhard force
   */
  inline void updateDtauChicUcl(double dt, double F)
  {
    _dto = dt/_part.getGamma();
    _xyz[_lpas].dtau = _dto;
    double value=5.00*fabs(F)*_part.getGamma();
    _xyz[_lpas].chic = value; // [p14]
    _uc = CONSTANT_p15/_dto; // [p15]
    _xyz[_lpas].ucl = _uc;
  }

	
  /**
   *  \brief Update the snake for the new step
   *
   *  update at every new step vectors calling addLinePositrons, updateGamdv2 and check if the snake is mature calling maturation  
   *
   * \param bool: true if the particule migrated
   * \param double: potential energy
   * \param double: kinetic energy
   * \param double: angular momentum
   */
  void updateNewStep(bool migrated,  double PotentialEnergy, double KineticEnergy, double angularMomentum);


  /**
   *  \brief Register the photon emitted
   *
   *  Register the photon emitted by the particle in a string variable  
   *
   *  \return string : data concerning the photon emitted by the particle
   */

  int makeStepKumakhov();

  //  bool multipleDiffusion(double dt,  double Dcoll);

	
  /**
   *  \brief Update moments
   *
   *  update moments calling addLineMoments and updateDp 
   *
   * \param double: 
   * \param double: 
   */
 
  inline void updateMoments(double pxa, double pya)
  {
    _xyz[_lpas].px = _part.getPx();
    _xyz[_lpas].py = _part.getPy();

    //intensite locale d'emission  [p18]  
    updateDp(pxa,pya);
  }


  void printPoirot() const
  {
    cout << " SNAKE : " << endl;

    cout << " LA= " << _la - _snakeBegin << " L1= " << _l1 - _snakeBegin << " LB= " << _lb - _snakeBegin << " L2= " << _l2 - _snakeBegin << " LPAS= " << _lpas - _snakeBegin << " LMAX= " << _lmax << endl;
    cout << " MATURE= " << _mature << endl;
    cout << "  UMIN= " << _umin << " USUP= " << _usup << endl;
    cout << " L,X/Y/Z/CHIC for last 10 steps " << endl;
    int ll = max(0, _lpas - 10);
    int k;

    for ( k = ll;  k <= _lpas; k++)
      {
	cout << " L= " << k  - _snakeBegin << " XL(L)= " << _xyz[k].xl << " YL(L)= " << _xyz[k].yl << " ZL(L)= " << _xyz[k].zl << " CHIC(L) " << _xyz[k].chic  << endl;
      }

    cout << endl;
    cout << " L,PXL/PYL/DTAU/dp" << endl;

    for ( k = ll;  k <= _lpas-1; k++)
      {
	cout << " L= " << k - _snakeBegin << " PXL(L)= " << _xyz[k].px << " PYL(L)= " << _xyz[k].py << " DTAU(L)= " << _xyz[k].dtau << " dp(L) " << _xyz[k].dp  << endl;
      }

    cout << " DTO= " << _dto << " UC= " << _uc << endl;

  }

 private: 

	
  inline int snakeIndex(int relativeIndex) const
  {
    //   cout << " index relatif " << relativeIndex << " absolu " << relativeIndex + _snakeBegin << endl;
    return relativeIndex + _snakeBegin;
  }



  inline void updateDp(double pxa,double pya) {
    double px = _part.getPx();
    double py = _part.getPy();
    _xyz[_lpas].dp = sqrt( (px - pxa)*(px - pxa) + (py - pya)*(py - pya) );
  }


  /**
   *  \brief Check conditions of maturation of the snake
   *
   *  Modify the variable _mature depending on maturation conditions  
   *
   * \param bool: true if the particule migrated
   * \param double: potential energy
   * \param double: kinetic energy
   * \param double: angular momentum
   */
  void maturation(bool migrated,  double PotentialEnergy, double KineticEnergy, double angularMomentum );

	
  /**
   *  \brief Control the firing of photons of KUMA
   *
   *   Control the firing of photons of KUMA and in case of sucessful, essais() defines a new trajectory 
   *
   */
  void essais();

	
  /**
   *  \brief Cut the tail of the snake
   *
   *  Cut the tail of the snake by removing all values of snake vectors, that correspond to an index between 0 and l1 
   *
   */
  void cutTale();

	
  /**
   *  \brief Draw at random the impulse of the photon
   *
   *  Draw at random the impulse of the photon and calculate the amplitude at square of the emission of the photon 
   *
   * \param int:
   * \param double&: 
   * \param double&: 
   * \param double&: 
   * \param double&: 
   * \param int&:
   * \param double&:
   * \param int&: 
   * \param int&:
   */
  void kuma(int lessai, double & U, double& OX, double& OY,double& AA,  int& LM, double& GM, int& IHU, int& IHD);


  // void difmul (int ncol, double& Q, double& QX, double& QY);


  // void abrems( double qq, double qx, double qy, double& AA, double& ZETA, double& GATX, double& GATY);

};



// c p14 - chi=gamma*F/m^2 ; 
// c       CHIC = (chi/249) * [coef.num. devant DABS(F) dans DTINV   [p12,p13] 
// c       1/259 = hbar/mc = lambdabar Compton, en angstrom.

// c p15 - UC est un cut-off local pour U = E(photon)/E(electron final)
// C  - En champ uniforme (le premier terme de DTINV1ou2 domine),
// C    UC est proportionel. a  CHI:
// C    UC=(0.0243)*(coef.num. de |F| dans DTINV1ou2)*259*CHI.
// C  - En champ rapidement variable (deuxieme terme dominant),
// C    UC est inversement. prop. a l'echelle des inhomogeneites 
// C    de la trajectoire.
// C  Dans l'amplitude d'un photon a U voisin de UC, 
// C  la variation de la phase pendant une inhomogeneite doit etre
// C  de l'ordre de 1 radian au plus.
// C  Une condition SUFFISANTE est: coef.num. de UC > 2/259.
// c  * Verifier le coef.num. a l'aide d'un test WCLASS :
// C  Si WCLASS simule < WCLASS calcule par la regle de Lienard, augmenter
// C  le coef.

// c p18 - le nombre de photons par pas est estime a priori a |dp|
// c      Cette quantite interviendra dans la facteur DENOR de ESSAI et KUMA,
// c      ainsi que dans la recherche du "point d'emission" de KUMA.
// c      La version ancienne utilisait chi*dtau (d'apres Landau "QED" par.90)


#endif


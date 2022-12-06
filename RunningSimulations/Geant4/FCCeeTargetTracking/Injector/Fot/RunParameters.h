#ifndef _RUNPARAMETERS_H
#define  _RUNPARAMETERS_H


/**
 * \file RunParameters.h
 * \brief The RunParameters class provides a useful way to gather some parameters of the crystal and photons of Kuma
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */

#include "Crystal.h"

// c p2 ---------------------------------------------
// C   ETMAX et VTMAX sont les cut-offs d'ejection (si l'electron 
// C   franchit l'un d'eux, il est considere comme en amorphe et ejecte
// C   du programme).
// C   Si  ETMAX(eV)>0.5*10**9*E(GeV)*VTMAX**2, l'ejection depend de VT
// C   Si  ETMAX(eV)<0.5*10**9*E(GeV)*VTMAX**2, l'ejection depend de ET
// C  FREKUMA gouverne la frequence des essais de photons KUMA  [es6]
// C  Parametre POIMIN: mode normnal, POIMIN=1
// C  Mode cristal mince, POIMIN < 1             [es10]
// C  En mode "cristal mince", augmenter FREKUMA.


/** \class RunParameters
 *
 *  The RunParameters class is the object that will contain some data concerning the crystal and photons of Kuma
 */
class RunParameters
{


  Crystal _crys;

		
  double _Zexit;

  //
  double _vtmax;
  double _etmax;

  double _phomin; // ENERGIE MINIMUM DU PHOTON EN GEV
  double _frekuma; // controle la frequence des tirs (voir note es7)
		
		
  //  Parametre POIMIN: mode normnal, POIMIN=1
  //  Mode cristal mince, POIMIN < 1             [es10]
  double _poimin; 
 public:
		
  /**
   *  \brief Constructor
   *
   *  Construct the RunParameters object from a crystal 
   *
   *  \param Crystal : the crystal where the particle move
   *
   */
		
  RunParameters(Crystal crystal) 
    {
			
      _crys=crystal;
			
    }

  void  set(double phomin, double etmax, double vtmax, double poimin, double frekuma)
  {
    _phomin = phomin;
    _etmax = etmax;
    _vtmax = vtmax;
    //  Parametre POIMIN: mode normnal, POIMIN=1
    //  Mode cristal mince, POIMIN < 1             [es10]
    _poimin = poimin;
    _frekuma = frekuma;
			
  }
		
  /**
   *  \brief return the Z position of emmission of the photon
   *
   *  Allows access to the Z position of emmission of the photon
   *
   * \return double: return the Z position of emmission of the photon _zexit
   */

  inline void setZexit(double z) 
  {
    _Zexit = z;
  }

		
  inline double getZexit() const
  {
    return _Zexit;
  }
		
		
  /**
   *  \brief return the minimum energy of the photon in GeV
   *
   *  Allows access to the minimum energy of the photon in GeV
   *
   * \return double: return the minimum energy of the photon in GeV _phomin
   */
		
		
  inline double getPhomin() const
  {
    return _phomin;
  }
		
  /**
   *  \brief return the testing frequency of KUMA
   *
   *  Allows access to the testing frequency of KUMA
   *
   * \return double: return the testing frequency of KUMA _frekuma
   */
		
  inline  double getFrekuma() const
  {
    return _frekuma;
  }
		
		
  /**
   *  \brief return the weight of the crystal
   *
   *  Allows access to the weight of the crystal
   *
   * \return double: return the weight of the crystal _poimin
   */
		
  inline   double getPoimin() const
  {
    return _poimin;
  }
		
		
  /**
   *  \brief return the crystal of the particle
   *
   *  Allows access to the crystal in which the particle is
   *
   * \return Crystal&: return the crystal _crys
   */
		
  inline const Crystal& getCrys() const {
    return _crys;
  }

  inline void getEjectionCut_off(double& etmax, double& vtmax) const
  {
    etmax = _etmax;
    vtmax = _vtmax;
  }		
};
#endif

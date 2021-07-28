#include "Evenement.h"

#ifdef TEST_WCLASS
#include "wclass.h"
#endif

Evenement::Evenement( Snake* snak, BremsStrahlung* brms, PhotonCollection* photons, double etmax, double vtmax, double Zexit) : _snak(snak), _bremse(brms), _photons(photons), _Zexit(Zexit), _vtmax(vtmax), _etmax(etmax)
{
   

//   _vtmax = 1e-2;
//   _etmax = 100000.0;

  _zj = 0.0;

  _partCrys = _snak->getParticlePtr();

  bool migre = _partCrys->migration();
  //Modification de MIGGG
  _partCrys->setBigJump(false);
  if (migre) _zj = _partCrys->moduleOfAngularMomentum();
  _restartSnake = true;

}

bool Evenement::makeStep()
{
  

  if ( _partCrys->getIexit()  )  return false;
    
	

  if ( _restartSnake && !reInitSnake() ) return false;

  double dt=2/(_dtinv2 + _dtinv1); // [p9]

  //avancement de la particule et actualisation de la cellule
  _partCrys->move(dt);
  if ( _partCrys->getZPosition() >= _Zexit   )
    {
      _partCrys->setIexit(true);
    }
	


	 
  //migration??
		
  _zj = _partCrys->moduleOfAngularMomentum();
  bool migre = _partCrys->migration();
		

  if (_partCrys->bigJump()){
    cerr << " jump is to big " << endl;
    throw string("poirot");
  }
		  
		
  //calcul des potentiels
  //   double xx,yy, rayon;
  //   _partCrys->getTransversePosition(xx, yy, rayon);
  //   double charge = _partCrys->getCharge();

  //  Lindhard lind = _partCrys->getCrystal().calculationOfFields(xx,yy,rayon,charge);

  Lindhard lind = _partCrys->calculationOfFields();

  double epot = lind._epot;
		

  _snak->updateNewStep(migre,epot,_ecin,_zj);


  int continuer = _snak->makeStepKumakhov();

  if ( continuer != 2 ) 
    {
      if ( continuer == 1) 
	{
	  _restartSnake = true;
	  return true;
	}
      else return false;
    }
  else _restartSnake = false;

  //première actualiation du pas : dt1, ou plutot DTINV1 = 1./dt1  [p12]

  double rayon = _partCrys->getRadius();

  // PAS TRANSVERSE MAXIMUM   ?? a augmenter ??
  double drmax=0.2*(rayon + _partCrys->getCrystal().getUtherm());
  _dtinv1=5*fabs(lind._f) + _vt/drmax;  // [p13a]
  dt=2.0/(_dtinv2+_dtinv1);

  _snak->updateDtauChicUcl(dt, lind._f);
		
  //acceleration transverse
  double fx=lind._fx;
  double fy=lind._fy;
  double pxa = _partCrys->getPx();
  double pya = _partCrys->getPy();
  _partCrys->acceleration(dt,fx,fy);

#ifdef TEST_WCLASS
  WCLASS += lind._f * lind._f * dt; // (seulement POUR TEST REGLE DE SOMME)
#endif

  bool bremsReussi = false;

#ifndef TEST_WCLASS
  bremsReussi = _bremse->multipleDiffusion( dt, lind._dcoll);
#endif

  if ( bremsReussi ) 
    {
      _restartSnake = true;
      return true;
    }
		

  _snak->updateMoments( pxa, pya);

  double pt2 = _partCrys->Ptsquare();
  _ecin = 0.5*pt2*ELECTRON_MASS_EV/_gamma;
  if (pt2 > _pt2max )
    {
      return false;
    }
  else
    {		
      //seconde actualisation du demi pas 'DT2'
      _vt = sqrt(pt2)/_gamma;
      _dtinv2 = 5.00*fabs(lind._f) + _vt/drmax; // [p13b] 
      return true;
    }
  
}

bool Evenement::reInitSnake()
{

  _gamma = _snak->getParticle().getGamma();
		
  _pt2max = min((_vtmax*_gamma)*(_vtmax*_gamma) , 2.0*_gamma*_etmax/ELECTRON_MASS_EV); // [p8]

  //initialisation du serpent
  _snak->initSnake();
	    
	
  // si l'electron est trop lent, evenement suivant 	
  if ( _snak->getXiMin() >= 0.9999 )  return false;


  double pt2 = _snak->getParticle().Ptsquare();

  _vt = sqrt(pt2)/_gamma;
	      
  _ecin = 0.5*pt2*ELECTRON_MASS_EV/_gamma;

#ifndef TEST_WCLASS
  _bremse->reInit(&_snak->getParticle());
#endif
  //préparation des variables pour calculer dt dans le move(dt)
  _dtinv2 = 1.e30;
  _dtinv1=_dtinv2;
  return true;
}



// NOTES DE XAVIER ARTRU
// c *************  NOTES du programme principal **********************
// c p1 ----------------------------------------------------
// C     Unification relativiste temps=longueur, masse=energie, etc. (c=1). 
// c     Unite de temps ou distance : ANGSTROM
// C     PT=(PX,PY)=GAMMA*(VIT.TRANSVERSE)
// C     =(impulsion transverse)/m(elec)
// C     Par contre, l'energie totale E est en GEV
// C     DTO ou DTAU(L) = pas en temps propre = DT/GAMMA
// C     Potentiel POT et energie transverse ET en eV

// c p8  --  PT^2 > PT2MAX: particule consideree comme en amorphe. Voir [p2]

// c p9 --------------------------------------------------------------
// C  DT1 et DT2 sont pratiquement nuls initialement.
// C  Le premier passage (LPAS=0) ne changera pratiquement pas X,Y,Z,PX,PY





// c p12 -------------------------------------------------------------
// C   POUR DIMINUER LE PAS: 
// c   - augmenter le coefficient de DABS(F) dans DTINV1 et DTINV2
// c   - diminuer le coefficient de DRMAX
// C  (ATTENTION! CECI ENTRAINE UNE AUGMENTATION DE UC=UCL).
// C  Coefficients conseilles: 
// C    0.2 pour DRMAX,  
// C    5.0 pour DABS(F).
// C   Verifier que le pas est assez petit avec un test de conservation 
// C   de l'energie transverse.

// c p13a,b ---------------------------------------------------------------
// C   Le premier terme sert a LIMITER LA VARIATION DE PT/M  (=VT*GAMMA).
// c   Le deuxieme sert a limiter le pas transverse a DRMAX, qui est d'autant
// c   plus petit que R est petit (on veut que le vecteur force varie peu
// c   le long du pas). 
// c   a) la premiere actualisation (DTINV1) a lieu apres celle de F
// c   b) la deuxieme actualisation (DTINV2) a lieu apres celle de VT




// c p19 - coef.num. 0.96E-8 = ?? 

// C-----------------------------------------------------------------------------

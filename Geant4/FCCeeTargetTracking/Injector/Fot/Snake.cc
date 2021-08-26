
#include "Snake.h"
#include "mathematics.h"

#ifdef TEST_WCLASS
#include "wclass.h"
#endif

// c p4 - L1MAX (voir REPTIL) est choisi = 0.7 fois la dimension de XL, etc.

#define _SNAKE_SIZE 2048
#define _L1MAX int(0.7*(_SNAKE_SIZE-1))

Snake::Snake(const RunParameters & rp, ParticleInCrystal & p, PhotonCollection& fcol, statistiques* stat): _photons(fcol), _runPar(rp), _part(p), _stat(stat), _mature(false),_jzaug(false) {


  _eloigne = false;

  _ra = -1.0;
  _partIsPositron = _part.isPositron();

  _xyz.resize(_SNAKE_SIZE);
  _snakeBegin = 0;

  _ancoll = 0.0;
  //_mcoll = 0;

  _dto = 0.0;
	
  _uc = 0.0;

}

void Snake::cutTale() {

#ifdef DO_STATS
  _stat->incrementMrept();
#endif
  _snakeBegin = _l1;


  _la=_lb;
  _l1=_lpas-1;
  _lmax=std::min(_SNAKE_SIZE-1, _l1 + _L1MAX);


  _mature = false;
}
    

void Snake::updateNewStep(bool migrated,  double PotentialEnergy, double KineticEnergy, double angularMomentum)
{
  _lpas++;
  if ( _lpas == _SNAKE_SIZE) 
    {
      _xyz.erase( _xyz.begin(), _xyz.begin()+ _snakeBegin); 
      _xyz.resize(_SNAKE_SIZE);
 
      _l1 -= _snakeBegin;
      _l2 -= _snakeBegin;
      _la -= _snakeBegin;
      _lb -= _snakeBegin;
      _lpas -= _snakeBegin;
      _snakeBegin = 0;
    }


#ifdef DO_STATS
  _stat->incrementMpas();
#endif

  if (_lpas - _snakeBegin == _lmax){
#ifdef DO_STATS
    _stat->incrementNlmax();
#endif
    _mature=true;
  } // [p10]

  SnakeStep& aux = _xyz.at(_lpas);

  aux.xl = _part.getXPosition();
  aux.yl = _part.getYPosition();
  aux.zl = _part.getZPosition();


  //   premier cut-off angulaire GAMDV2 du kumakhov   [p7]
  double toto=_part.getGamma()*fabs(PotentialEnergy)/ELECTRON_MASS_EV;
  aux.gamdv2 = min(toto,toto*fabs(PotentialEnergy)/(KineticEnergy+1e-9));
  maturation(migrated, PotentialEnergy, KineticEnergy,angularMomentum );
}




// la valeur de retour peut etre 0, 1 ou 2
// 0 : interrompre le traitement de l'evenement
// 1 : emission
// 2 : rien ( pas d'emission, serpent non mur) : on continue
int Snake::makeStepKumakhov()
{
  if ( !_mature ) return 2; 
  else 
    {
      _mature = false;
      _l2 =  _lpas - 1;
      if ( _lpas - _snakeBegin != 1 )
	{
	  essais();
	  if (_part.getIexit() ) return 0;
	  else
	    {
	      if ( _part.getIemis() ) return 1;
	      else
		{
		  cutTale(); // [p11a]
		}
	    }
	}
      if (_part.getIexit() ) return 0;
      else return 2;
	
    } // fin mature
}



void Snake::maturation(bool migrated,  double PotentialEnergy, double KineticEnergy, double angularMomentum ){

  //  updateGamdv2(PotentialEnergy,KineticEnergy);
  if ( _part.getIexit() ) 
    {
      _mature = true;
    }
  //actualise la variable _mature lorsque le serpent est mature EN CAS DE MIGRATION
  if (migrated){
		
    bool jzaugav = _jzaug;
    _jzaug = _part.moduleOfAngularMomentum() < angularMomentum;
    if(jzaugav && !_jzaug) 
      {
	_mature = true;	
      }
  }
	
  // double Px = _part.getPx();
  // double Py = _part.getPy();
  double r = _part.getRadius();

  bool eloigna = _eloigne;
  _eloigne = r > _ra;
  _ra = r;
  // double Et = PotentialEnergy + KineticEnergy;
  bool canalise = (!_partIsPositron)&&(!_eloigne);
  bool lapogee = eloigna && (!_eloigne);
  if (lapogee && canalise)  _mature = true; // [p11]

}


	

void Snake::initSnake() {



  _snakeBegin = 0;
	
  _l1=0;
  _la=0;
  _lmax = _L1MAX;	
  _lpas=-1;

  _ximin = _runPar.getPhomin()/(_part.getGamma()*ELECTRON_MASS_GEV);

  _umin = _ximin/(1.0 - _ximin);

#ifdef TEST_WCLASS  
  _umin = _ximin; // pour test WCLASS
#endif

  _xmilog = - log (_ximin); 
  _mature = false;
  _jzaug = false;
  _ra = _part.getRadius();	
  _eloigne = false;
  _ancoll = 0.0;
}




//C  commande les tirs des photons KUMA et, en cas de reussite,
//C  definit la nouvelle trajectoire  
void Snake::essais()
{

  double poid; // je declare ici cette variable, car dans le fortran elle n'est pas
  // toujours initialisee. Je ne sais pas trop ce que ca fait. GLM
  double epho;
  double proba;
	
  _part.setIemis(false);
  if ( _l2 <= _l1 )
    {
      cerr << "     essais : ERROR;  L2 <= L1 , l1= " << _l1 << " l2 = " << _l2 << " lpas " << _lpas << endl;
      throw string("poirot");
    }


  //  CALCUL DE LB = FIN DE LA ZONE DE TIR (LA= DEBUT) [es1]

  double sup = -1.0;
  int lessai;
  double pxxl1, pxxl2,pyyl1, pyyl2;
  double distance;
  double pxless, pyless;

  pxxl1 = _xyz.at(_l1).px;
  pxxl2 = _xyz.at(_l2).px;
  pyyl1 = _xyz.at(_l1).py;
  pyyl2 = _xyz.at(_l2).py;

 
  for ( lessai = _l1; lessai < _l2; lessai++)
    {

      pxless = _xyz.at(lessai).px;
      pyless = _xyz.at(lessai).py;

      distance = (pxless - pxxl1) * (pxless - pxxl1) + (pyless - pyyl1) * (pyless - pyyl1);
      distance *= (pxless - pxxl2) * (pxless - pxxl2) +  (pyless - pyyl2) * (pyless - pyyl2);
      if ( distance > sup )
	{
	  sup = distance;
	  _lb = lessai;
	}
    }       //  [es2]


  if (_part.getIentree() ) _la = snakeIndex(0);
  if (_part.getIexit() ) _lb = _l2;

  if ( _l1 == _snakeBegin && !_part.getIexit() ) return; // [es3]
  
  if ( _la > _lb ) 
    {
      cerr << " ERROR :  LA > LB in Snake::essais  " << endl;
      throw string("essais");
    }

  //   CALCUL DE USUP POUR SERPENT        ![es4]
  _xyz.at(snakeIndex(0)).ucl = _xyz.at(snakeIndex(1)).ucl;  // [es5]
  _usup = 0.0;
  for (lessai = _la; lessai <= _lb; lessai++) 
    {
      _usup = max(_usup, _xyz.at(lessai).ucl);
    }
  if ( _usup <= _umin )  //  (trajectoire trop douce)
    {
      _part.setIentree(false); 
      return;
    }

  double sx = _part.getXPosition();
  double sy = _part.getYPosition();
  double sxCell =  _part.getX0cell();
  double syCell = _part.getY0cell();


  //  BOUCLE DE TIRAGE

  _xisup = _usup/(1.0 + _usup);

#ifdef TEST_WCLASS
  _xisup = _usup;
#endif



  double h0 = log(_xisup/_ximin)*_runPar.getFrekuma();  // [es6]

  double gamma = _part.getGamma();
  double poimin = _runPar.getPoimin();
  bool newTraj = false;
  double Ulocal;
  double OX, OY;
  int LM;

  for ( lessai = _la; lessai <= _lb; lessai++)
    {
      //   nombre de tirs au point(LESSAI)   [es7]
      double tirs = h0 * _xyz.at(lessai).dp; // [p18]
      //     (on a enleve /259 en redefinissant FREKUMA)

      double hlocal = mathematics::rndm();

      if ( tirs >= hlocal ) // si tirs < h, on passe au point suivant
	{
	  int ntirs = floor(tirs);
	  tirs -= ntirs;
	  if ( tirs >= hlocal ) ntirs++; // [es7a] 
	  //          - et s'il y a une diffusion a angle >> 1/gamma ???
	  //          - j'essaye d'en tenir compte en prenant DP plutot que CHI
	  int IHU = 0; // initialisation (n'importe quelle valeur entre 0 et 10)
	  int IHD = 0;
	  int ijk;
	  for (ijk = 0; ijk < ntirs; ijk++)
	    {

	      double AA, GM;

	      kuma(lessai, Ulocal, OX, OY, AA, LM, GM, IHU, IHD);

	      proba = max( AA, 0.0)/_runPar.getFrekuma();
	      epho = gamma * ELECTRON_MASS_GEV * Ulocal/(1.0 + Ulocal);


#ifdef TEST_WCLASS
	      epho = gamma * ELECTRON_MASS_GEV * Ulocal;
#endif


	      //   sttatistiques SUR LES PHOTONS KUMAKHOV:
#ifdef DO_STATS
	      _stat->addEnergyKuma(proba, epho, IHU, IHD, AA);
#endif

#ifndef TEST_WCLASS
	      //   decision aleatoire d'emettre le photon (sous reserve de passer ETAFIN)

	      hlocal = mathematics::rndm();
	      if ( proba >= hlocal * poimin )    // [es10]
		{

		  //   position initiale pour la nouvelle trajectoire 
 		  
		  //  "derapage"   [ef3]
		  double dpx = _xyz.at(LM).px - _xyz.at(LM - 1).px;
		  double dpy = _xyz.at(LM).py - _xyz.at(LM-1).py;
		  double toto = 0.5* _xyz.at(LM).dtau * Ulocal * GM /(dpx*dpx + dpy*dpy);
		  _part.derapeXY( _xyz.at(LM).xl + toto*dpx, _xyz.at(LM).yl + toto*dpy ); // [ef3]

		  //   condition de non-migration [ef4]
		  _part.setIemis(!_part.migration());
		  //Modification de MIGGG
		  _part.setBigJump(false); //  (si on veut continuer malgre un MIGGG)
	      
		  //       cas ETAFIN reussi
		  if ( _part.getIemis() ) 
		    {
		      //            enregistrement du photon 
		      double thetx = OX / gamma;
		      double thety = OY / gamma;
		      poid = max( proba, poimin);    // [es10]
		      _photons.fill(epho, thetx, thety, _part.getXPosition(), _part.getYPosition(), _part.getZPosition(), poid);
#ifdef DO_STATS
		      _stat->incrementMfot();
#endif			
		      //            decision de prendre la nouvelle trajectoire [es10]
		      if ( proba > hlocal)   // decide oui 
			{
			  newTraj = true;
			  break;
			}
		    }
#ifdef DO_STATS
		  else // cas ETAFIN rate - a supprimer si ca n'arrive presque jamais ???
		    {
		      _stat->incrementAvort(poid, epho);
		    }
#endif
		}
#endif
	    } // fin boucle ijk (etiq. fortran 52)
	  if ( newTraj ) break;



	} // fin if tirs
    } // fin lessai

  if (!newTraj ) 
    {
      //  cas ou il n'y a pas de photons emis - du moins, cas ou on ne modifie
      //  pas la trajectoire (aucun photon avec POID > 1.)
      _part.setIemis(false);
      _part.setIentree(false);
      //        rappel de X,Y sauvegardes (efface le derapage)
      _part.reIniCoorAndCell(sx,sy, sxCell, syCell);
    }
  else // Le photon est enregistre ET on modifie la trajectoire
    {
      _part.setZ(_xyz.at(LM).zl);
      _part.setIentree(false);
      //    nouvelle energie-impulsion   [ef2]
      _part.setEnergyMomentum(_xyz.at(LM).px - OX*Ulocal/(1.0 + Ulocal),
			      _xyz.at(LM).py - OY*Ulocal/(1.0 + Ulocal),
			      gamma/(1.0 + Ulocal));
		
		
#ifdef DO_STATS
      _stat->addPhotKuma(proba, epho);
#endif
    }
}






// C  TIRE AU HASARD L'IMPULSION DU PHOTON,
// C  CALCULE L'AMPLITUDE AU CARRE D'EMISSION DU PHOTON.
void Snake::kuma(int lessai, double & U, double& OX, double& OY,double& AA,  int& LM, double& GM, int& IHU, int& IHD)
{

  //   try
  //     {
  vector<double> alocal(6,0.0);
  vector<double> gatx(_SNAKE_SIZE);
  vector<double> gaty(_SNAKE_SIZE);
  vector<double> glocal(_SNAKE_SIZE);
  

  double phi = 0.0;
  double aux;



  //   1. tirage de l'energie et de l'angle du photon   [k1]

  //  TIRAGE  DE  U= hbar*OMEGA/(E-hbar*OMEGA)
  double h = mathematics::rndm();
  double xi = _ximin*pow(_xisup/_ximin, h);
  U = xi/(1.0 - xi);

#ifdef TEST_WCLASS
  U = xi;
#endif

  if ( U > _xyz.at(lessai).ucl ) 
    {
      AA  = 0.0;
      return;
    }

#ifdef DO_STATS
  _stat->incrementNkuma();
#endif
  // ici MODIF DE MAI 2010
  double xmax = _xyz.at(lessai).ucl/(_xyz.at(lessai).ucl + 1.0);
  xmax = min(xmax, _xisup);
  double hprime = h * log(_xisup/_ximin) / log(xmax/_ximin);
  IHU = floor( 10.0 * hprime);

  //  TIRAGE DE LA DIRECTION (OX,OY) DU PHOTON

  double dlta02 = pow( _usup/(100.0 * U), 0.67); 
  // 0.67 : 2/3, valeur theorique (voir documentation)
  dlta02 = min( dlta02, _xyz.at(lessai).gamdv2) + 1.0;

  h = mathematics::rndm();
  IHD = floor(10.0 * h); //  POUR EFFIC

  double gat = dlta02 * ( 1.0/h - 1.0 );
  gat = sqrt( gat );

  {
    const double phi_ = 2.0 * M_PI * mathematics::rndm();
    OX = _xyz.at(lessai).px - gat * cos(phi_);
    OY = _xyz.at(lessai).py - gat * sin(phi_);
  }

  //   2. CALCUL DE LA PROBABILITE D'EMISSION     [k4]

  double omegag = U*MC_OVER_HBAR/(8.0*M_PI);

  // C   PREPARATION . LIMITES D'INTEGRATION LU,LV [k6],
  // c   POINT D'EMISSION LM,
  // c   FACTEUR NORMALISANT DENOR [k,k7].

  double denor = 0.0;
  // double fsg2ma = 0.0;

  //   tabulation de GATX,GATX,G, recherche des bornes d'integration LU,LV
  int lu = 0, lv;
  double fsg2 = 0.0;
  double fsg2mi = 1.0e30; // infinity
  int l;
  for (l=_snakeBegin; l <= _l1; l++)
    {
      gatx.at(l) = _xyz.at(l).px - OX;
      gaty.at(l) = _xyz.at(l).py - OY;
      glocal.at(l) = 1.0 + gatx.at(l)*gatx.at(l) + gaty.at(l)*gaty.at(l);

      //       recherche de LU

      fsg2 = _xyz.at(l).chic/(glocal.at(l)*glocal.at(l)); //  FSG2 = "f sur g^2"
      if ( fsg2 < fsg2mi ) 
	{
	  lu = l; // LU entre 0 ET L1
	  fsg2mi = fsg2;
	}
      if ( l >= _la && U < _xyz.at(l).ucl )
	{
	  aux = dlta02 + glocal.at(l) - 1.0;
	  denor += _xyz.at(l).dp/( aux*aux );
	}
    }

  if (_part.getIentree() ) lu = 0;

  //   - recherche de LV
  lv = _l1;
  fsg2mi = fsg2;
  for (l= _l1+1; l <= _l2; l++)
    {
      gatx.at(l) = _xyz.at(l).px - OX;
      gaty.at(l) = _xyz.at(l).py - OY;
      glocal.at(l) = 1.0 + gatx.at(l)*gatx.at(l) + gaty.at(l)*gaty.at(l);

      fsg2 = _xyz.at(l).chic/(glocal.at(l)*glocal.at(l)); 

      if ( fsg2 < fsg2mi ) 
	{
	  lv = l; // LV entre L1 ET L2
	  fsg2mi = fsg2;
	}

      if ( l <= _lb && U < _xyz.at(l).ucl )
	{
	  aux = dlta02 + glocal.at(l) - 1.0;
	  denor += _xyz.at(l).dp/( aux*aux );
	}
    }
  if (_part.getIexit() ) lv = _l2;



  if ( lu == lv ) 
    {
      AA = 0.0;
      return;
    }

  if ( denor == 0.0 ) 
    {
      cerr << " KUMA : ERREUR denor = 0  "  << endl;
      cerr << " lessai= " << lessai <<  " lu= " << lu << " lm= " << LM <<  " lv= " << lv << " l= " << l   << endl;
      cerr << " denor= " << denor << " g(l)= " <<   glocal.at(l) << endl;
      int kk;
      cerr << " DTAU " ;
      for (kk=_snakeBegin; kk <= _l2; kk++) cerr << " " << _xyz.at(l).dtau;
      cerr << endl;
      throw string("poirot");
    }

  //   FIN DE LA PREPARATION;INTEGRATION

  double vtau = 1.0/glocal.at(lu);
  double vx = gatx.at(lu) * vtau; // vitesse apparente   [k8]
  double vy = gaty.at(lu) * vtau; 


  double slocal = 0.0; // initialisations
  double clocal = 1.0; 
  double wlocal = 0.0; 
  double bsup = 0.0;

  //   boucle d'integration de l'amplitude

  for ( l = lu + 1; l <= lv; l++)
    {


      double dw = (glocal.at(l) + glocal.at(l-1)) * _xyz.at(l).dtau * omegag; // phase [k9]
      if ( dw == 0.0 ) 
	{
	  cerr << " KUMA : ERREUR dw = 0 " << endl;
	  cerr << " lessai= " << lessai <<  " lu= " << lu << " lm= " << LM <<  " lv= " << lv << " l= " << l   << endl;
	  cerr << " denor= " << denor << " g(l)= " <<   glocal.at(l) << " dw= " << dw << endl;
	  int kk;
	  cerr << " DTAU " ;
	  for (kk=_snakeBegin; kk <= _l2; kk++) cerr << " " << _xyz.at(l).dtau;
	  cerr << endl;
	  throw string("poirot");	      
	}

      wlocal += dw;
      double wwlocal = 2.0 * M_PI * (wlocal - floor(wlocal));
      double ca = clocal;
      double sa = slocal;
      
      clocal = cos(wwlocal);
      slocal = sin(wwlocal);

      //       calcul des cosinus et sinus moyens    [k10]
      double cm, sm;
      if ( dw > 1.0e-3 ) 
	{
	  cm = ( slocal - sa) / dw;
	  sm = ( ca - clocal) / dw;
	}
      else
	{
	  cm = M_PI * (clocal + ca);
	  sm = M_PI * (slocal + sa);
	}
      double vtaua = vtau; // antecedents
      double vxa = vx;
      double vya = vy;

      vtau = 1.0/glocal.at(l); 
      vx = gatx.at(l) * vtau; 
      vy = gaty.at(l) * vtau; 
      double dvx = vx - vxa; // accroissements
      double dvy = vy - vya;
      double dvtau = vtau - vtaua;

      //       recherche du point d'emission
      aux = max( vtaua, vtau);
      double blocal = max(  _xyz.at(l).dp, _xyz.at(l).chic*_xyz.at(l).dtau) * aux * aux;
      if ( blocal > bsup ) 
	{
	  bsup = blocal;
	  LM = l;
	  phi = wlocal;
	} // [k14]

      alocal[0] += cm * dvx;
      alocal[1] += sm * dvx;
      alocal[2] += cm * dvy;
      alocal[3] += sm * dvy;
      alocal[4] += cm * dvtau;
      alocal[5] += sm * dvtau;

    } // fin de la boucle d'integration

  double a0 = alocal[0];
  double a1 = alocal[1];
  double a2 = alocal[2];
  double a3 = alocal[3];
  double a4 = alocal[4];
  double a5 = alocal[5];

  AA = a0*a0 + a1*a1 + a2*a2 + a3*a3;
  AA *=  1.0 + (1.0 + U)*(1.0 + U);
  AA += U * U * (  a4*a4 + a5*a5 );

  //   AA = alocal[0]*alocal[0] + alocal[1]*alocal[1] + alocal[2]*alocal[2] + alocal[3]*alocal[3];
  //   AA *=  1.0 + (1.0 + U)*(1.0 + U);
  //   AA += U * U * (  alocal[4]*alocal[4] + alocal[5]*alocal[5] );

  AA /= (1.0 + U)*(1.0 + U) * denor * dlta02 * EIGHT_PI_CUB_OV_ALPHA;



  GM = glocal.at(LM);



#ifdef TEST_WCLASS
  AA = 2.0 * ( alocal[0] * alocal[0] + alocal[1] * alocal[1] + alocal[2] * alocal[2] + alocal[3] * alocal[4] );
  AA /= denor * dlta02 * EIGHT_PI_CUB_OV_ALPHA;
#endif

  //   REJET POUR PHASE DEBILE ?
  phi -= GM * _xyz.at(LM).dtau * omegag;

  if ( phi < 0.08 && !_part.getIentree() ) AA = -AA;
  phi = wlocal - phi;
  if ( phi < 0.08 && !_part.getIexit() ) AA = -AA;


}

// C p7 -----------------------------------------------------------    
// c   GAMDV2 = (gamma*angle typique de deflection)**2 au cours d'une 
// c   ondulation de la trajectoire.
// c   Cet "angle typique" est inferieur a - ou de l'ordre de l'angle de Lindhard.
// c   On admettra que c'est un cut-off angulaire pour KUMA


// c p10 ---------------------------------------------------------------
// C    Voir REPTIL; Si trop frequent, augmenter la dim. de XL, etc.
// C    On doit avoir NLMAX << NREPT


// C p11 ------------------------------------------------------------------
// C     Cas canalise: Le serpent est MATURE si la distance a l'axe le + proche
// C     VIENT DE PASSER par un maximum. 
// C     Cas non canalise: voir MIGRE


// c p11a  - LE SERPENT N'EST PAS MUR OU LES ESSAIS KUMA ONT ECHOUE


// c p18 - le nombre de photons par pas est estime a priori a |dp|
// c      Cette quantite interviendra dans la facteur DENOR de ESSAI et KUMA,
// c      ainsi que dans la recherche du "point d'emission" de KUMA.
// c      La version ancienne utilisait chi*dtau (d'apres Landau "QED" par.90)


// c ------------------------------------------------------------------------
// C LABEL DES POINTS DU SERPENT AVANT ET APRES REPTIL:
// C                      
// C   0             LA           L1            LB                L2   L(=LPAS ??
// C   .  .  .  . .  .  .  . . .  .  .  . . .   .  .  . .  .  .   .    .
// C                              0             LA                L1   L(=LPAS ??
// c ------------------------------------------------------------------------


// c es1 ------------------------------         
// C   Je veux que LB soit loin a la fois de L1 et de L2 dans l'espace
// C   des vitesses transverses (Idem pour LA,L0,L1, grace a REPTIL).

// c es2 ----------------------------
// C   LB entre L1+1 ET L2-1. Apres REPTIL : LA entre 1 ET L1-1,
// C   sauf SORTIE (LB=L2) , ENTREE (LA=0) , ou cas rare L2-L1=1

// c es3 - cas ou le serpent n'a qu'une revolution.

// c es4 - USUP est un cut-off global pour l'intervalle [LA,LB]
// c      On tiendra compte du cut-off local dans KUMA [k2]

// c es5 -----------------------------------
// C  MOTIF: A L'INITIALISATION D'UN NOUVEAU SERPENT, DT EST PRIS A PEU
// C  PRES NUL, D'OU UCL(0) EXAGEREMENT GRAND

// c es6 -  FREKUMA trop grand: programme trop lent
// C        FREKUMA trop petit: beaucoup de photons oublies (PROBA >1)

// c es7 - Un essai a chaque pas serait trop long. 
// c       Si on ignore la diffusion multiple, le nombre moyen "TIRS" 
// c       de tirs effectues au pas(LPAS) est proportionnel a la
// c       longueur du pas DTAU et a la valeur locale de CHI.
// c  Avec la diff. mul. on le generalise CHI*DTAU par DP(LESSAI)/259.
// c  (le facteur 259 etant omis, on a divise l'ancien FREKUMA par 259)  
// c  On tient compte de ce nombre moyen d'essais par le facteur DENOR dans KUMA

// c es7a - Le nombre (entier) de tirs realises, NTIR, est choisi au hazard:
// c  - soit l'entier juste au dessus, avec la probabilite = partie frac.(TIRS)
// c  - soit l'entier juste en dessous, avec la probabilite complementaire.
      


// c es10 - simulation normale : POIMIN = 1.  
// c       simulation forcee  : POIMIN < 1  (pour cible mince). Dans ce cas
// c       on enregistre plus de photons.
// c       Cependant on ne leur attribue qu'un POID normal,
// c       et on ne modifie la trajectoire que si le photon aurait ete accepte
// c       en simulation normale.  Voir [p2]



// c - NOTES/ETAFIN ---------------------------------------------------------

// c  e.f.1  On fait un retour en arriere au "point d'emission"
// c         LM defini dans KUMA pour initier la nouvelle trajectoire.
// c         (On devrait en principe aussi reprendre les X0CELL et Y0CELL
// c         correspondants, mais nous ne les avons pas gardes en memoire)  

// c  e.f.2  Dans le cas ou il n'y a pas de diffusion multiple,
// c         on suppose la conservation de l'impulsion dans le processus
// c         electron --> electron + photon. Ceci est fait dans ESSAI  

// c  e.f.3. La conservation de l'energie-impulsion est impossible 
// c         pour le processus electron --> electron + photon en l'absence
// c         de champ exterieur. Dans le champ cristalin,  elle est possible 
// c         si l'electron descend dans le potentiel d'une denivelee 
// c         = m_e * GM * U/(2*GAMMA)
// c         Ce "derapage" est de nature quantique, mais, pour un grand nombre
// c         de photons mous, le cumul des derapages equivaut a la force 
// c         d'Abraham-Lorentz. 
// c         Dans le cas ou il y a diffusion multiple, le derapage n'est pas
// c         necessaire, car il peut y avoir transfert d'impulsion longitudinal.
// c         La recette adoptee tient compte du premier cas et donne un derapage
// c         negligable dans le second.

// c  e.f.4  le derapage ne doit pas faire changer de cellule. 
// c         S'il y a changement (MIG), le photon est avorte.


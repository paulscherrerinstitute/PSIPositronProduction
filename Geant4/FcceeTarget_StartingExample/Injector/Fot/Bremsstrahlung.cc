#include "Bremsstrahlung.h"
#include "mathematics.h"
BremsStrahlung::BremsStrahlung(PhotonCollection& fcol, double phomin, double poimin,  statistiques* stat) : _photons(fcol), _phomin(phomin), _poimin(poimin), _stat(stat)
{
  _part = NULL;
  _ancoll = 0.0;
}

void BremsStrahlung::reInit(ParticleInCrystal* p)
{
    _part = p;
  _ancoll = 0.0;
 double xmin = _phomin/(_part->getGamma()*ELECTRON_MASS_GEV);
 _xmilog = - log (xmin); 
}


bool BremsStrahlung::multipleDiffusion( double dt,  double Dcoll)
{

  _ancoll += dt* Dcoll;  // [p17]
  double DTO = dt/_part->getGamma();
  double UC = CONSTANT_p15/DTO;


  // diffusion multiple
  if ( _ancoll > 1.0 ) 
    {
      int ncol = floor( _ancoll );
      _ancoll -=  ncol;  // reliquat
      double Q, QX, QY;
      difmul( DTO,  UC,ncol, Q, QX, QY);
      if ( _part->getIemis() )  return  true; //  bremsstrahlung reussi !
      else _part->incrementTransverseMomentum( QX, QY);
    }
  return(false);
}



void BremsStrahlung::abrems(double DTO, double UC, double qq, double qx, double qy, double& AA, double& ZETA, double& GATX, double& GATY)
{
  //   tirage de l'energie du photon    [b2]
  ZETA = 1.0 - exp( - _xmilog * mathematics::rndm() ); // = 1.-hbar*OMEGA/E
  //   tirage de l'angle                [b3]
  double g1 = 1.0/mathematics::rndm(); //  1+GAT1**2
  double gat1 = sqrt( g1 - 1.0); // GAMMA * ANGLE( ELEC.INCID./PHOTON )
  double phi =  2.0 * M_PI * mathematics::rndm(); //  (AZIMUT DE GAT1) * 1024/2PI
  
  GATX = gat1 * cos(phi);
  GATY = gat1 * sin(phi);
  
  // 1+GAT2**2	
  double g2 = 1.0 + (GATX + qx) * (GATX + qx) + (GATY + qy) * (GATY + qy);
  //   GAT2=GAT1+Q(VECTEURS):GAM FINAL * ANGLE(ELECTRON FINAL/PHOTON)
  
  double ulocal = 1.0 / ZETA - 1.0;
  if ( ulocal * DTO * (g1 + g2) <  CST_BIZARRE && ulocal < UC ){
    AA = 0;
    return; 
  }
  AA = ( (1.0 + ZETA * ZETA) * qq * g1 * g2 -2.0 * ZETA * (g2 - g1) * (g2 - g1) ) / (g1 * g1 + g2 * g2);
  
  if (AA * (qq -AA) < 0) {
    cout << " ERROR BREMS AA= " << AA << " QQ= " << qq << endl;
    throw string("poirot");
  }
  
  //   Une fois sur deux, on intervertit les vecteurs GAT1 ET -GAT2:	
  if (phi<M_PI){
    GATX = - (GATX + qx);
    GATY = - (GATY + qy);
  }
}

void BremsStrahlung::difmul(double DTO, double UC, int ncol, double& Q, double& QX, double& QY){

  //  tirage de Q^2                    [dif1]
  double qq = 0.0;

  for (int m = 1; m <= ncol; m++){
    qq +=  1.0 / (mathematics::rndm() + 1e-6);
  }
  qq = - qq * log ( 1e-9 + mathematics::rndm() );


  double rayon = _part->getRadius();


  if ( rayon <  2.0*_part->getCrystal().getUtherm() ) {
    qq *=  _part->getCrystal().getQqcoh(); // DIFFUSION SUR NOYAU
  }  
  else{
    //    il y a un petit ecart entre 1.5e-5 et 1/MC_OVER_HBAR
    // la ligne commentee est identique au fortran; il eput etre utile 
    // de la retablir pour des comparaisons
    //    qq *=  1.5e-5/ ( r * r + __part.getC().getRbohr2() ); //  DIFFUSION SUR ELECTRON
    qq /=  MC_OVER_HBAR * MC_OVER_HBAR * ( rayon * rayon + _part->getCrystal().getRbohr2() ); //  DIFFUSION SUR ELECTRON
  }  // [dif2]

#ifdef DO_STATS
  _stat->addDiffMul(ncol);
#endif

	Q = sqrt ( qq );
	double phi = 2.0 * M_PI * mathematics::rndm(); //  angle azimutal aleatoire *1024/2P
	QX = Q * cos(phi);
	QY = Q * sin(phi);

// C   TIRE SELON DQX*DQY*(1.D0-(1.D0+QQ)*DEXP(-QQ))/(PI*Q**4)
// C   ET FAIT LA CONVOLUTION DES NCOL DIFFUSIONS
	
	_part->setIemis(false);
	
	//   PREMIERE EPREUVE :                      [dif3]

	double proba = ALPHA_OV_PI * _xmilog * qq;
// C   Ce premier PROBA = majorant de la probabilite de bremsstrahlung:
// C   dN(Brems)/dLog(EPHO)=(ALPHA/PI)*AA
// C   et AA < QQ  .   (Ce premier PROBA peut depasser 1)

	double hlocal = mathematics::rndm();
	if ( proba < hlocal * _poimin ) return ; // aucune emission - [dif4]  

      
	//  DEUXIEME EPREUVE :
	double AA, ZETA, GATX, GATY;
	
#ifdef DO_STATS
	_stat->incrementNbrems();
#endif

	abrems( DTO, UC, qq, QX, QY, AA, ZETA, GATX, GATY);

	double GAMMA = _part->getGamma();
	if ( (GATX+QX)*(GATX+QX) + (GATY+QY)*(GATY+QY) + 1 > (GAMMA*ZETA)*(GAMMA*ZETA) ) return;
	
	proba = max ( proba, _poimin ) * AA/qq;
	
	double epho = ELECTRON_MASS_GEV * _part->getGamma() * (1.0 - ZETA );
	
	//statistique Sandrine
// 	st.setWb(proba * epho);
// 	st.setFotb(proba);
	
#ifdef DO_STATS
	_stat->addEnergyBrems(proba, epho);
#endif

	// DECISION ALEATOIRE D'ENREGISTRER LE PHOTON
	hlocal = mathematics::rndm();
	if ( proba < hlocal * _poimin ) return ;
	
	double thetx = ( _part->getPx() - GATX ) / _part->getGamma();
	double thety = ( _part->getPy() - GATY ) / _part->getGamma();
	double poid = max ( proba, _poimin );
	
	_photons.fill(epho , thetx , thety , _part->getXPosition(), _part->getYPosition(), _part->getZPosition() , poid);

#ifdef DO_STATS
	_stat->incrementMfot();
#endif	 
	// decision de modifier la trajectoire de l'electron
	if ( proba < hlocal ) return ; // (traj. non modifiee)
	
	_part->setIemis(true); // photon emis, traj. modifiee (ci-dessous)

	//   ETAT FINAL DE L'ELECTRON APRES EMISSION BREMSSTRAHLUNG

	double deltaPx = QX - ( 1.0 - ZETA ) * _part->getGamma() * thetx;
	double deltaPy = QY - ( 1.0 - ZETA ) * _part->getGamma() * thety;
	_part->incrementTransverseMomentum( deltaPx, deltaPy);
	_part->multiplyGamma(ZETA);

	_part->setIentree(false);
	
#ifdef DO_STATS
	_stat->addPhotBrems(proba, epho);
#endif	
}



// c p17 -----------------------------------------------------------------
// C    DCOLL=1/(L.P.M. de collision); ANCOLL=nombre moyen de collision
// c    a faire.

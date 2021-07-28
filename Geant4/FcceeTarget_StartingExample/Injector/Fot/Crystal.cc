
#include "Crystal.h"
#include "GlobalConstants.h"
#include "Lindhard.h"
#include <cmath>
#include <cstdlib>
#include <iostream>

#define LONG_POT 1024

using namespace std;

static Vdecanal formul(double r)
{

  // BAIER-KATKOV-STRAKHOVENKO  AXIAL 
      
  //  PARAMETRES DE LA FORMULE DU POTENTIEL

  double as = 0.215;
  double eta = 0.115;
  double v0 = 417.;
  Vdecanal vdeca;
  double x = r*r/(as*as);
  vdeca.V = v0*log1p(1.0/(x + eta));
  double denom = (x + eta ) * ( x + eta + 1.0 );
  vdeca.F = 2.0 * r * v0 / (as * as * denom);
  vdeca.deltaV = 4.0 * v0 * ((eta - x) * (eta + x) + eta ) / (as*as*denom*denom);
  return vdeca;
}

Crystal::Crystal(string atom, int axe)
{

  _potable = vector<double>(LONG_POT);
  _ftable = vector<double>(LONG_POT);
  _dcolta = vector<double>(LONG_POT);
	

  if ( atom == string("W") ) 
    {
      setTungstene(axe);
    }
  else 
    {
      cerr << " Crystal:: type of atom not programmed : " << atom << endl;
 
    }
  _qqcoh = 1.0/(MC_OVER_HBAR*MC_OVER_HBAR * _bcoh *   _bcoh);
}

Crystal::~Crystal() {
}

void Crystal::setTungstene(int axe)
{
  if ( axe != 111) 
    {
      cerr << " Crystal::setTunstene : unknown axe direction " << axe << endl;
      cerr << " only axe <111> programmed " << endl;
    }
  double zatom = 74.;

  _amaille = 3.165;
  //  _distatom = _amaille*SQROOT3*0.5; // (for W, diamant, Si, Ge)
  _distatom = _amaille*sqrt(3.0)*0.5; // (for W, diamant, Si, Ge)
  _utherm = 0.05; // at regular temperature

  // Thomas-Fermi radius
  double atf;
  atf = FACTOR_THOMAS_FERMI*RBOHR1/pow(zatom, 1.0/3.);
  struc111(string("BCC"), zatom);

  // parameters for coulombian diffusion
  // c 5 ------------------------------------------------------------  
  // c     BCOH= PARAMETRE IMPACT MAXIMUM POUR DIFFUSION INCOHERENTE SUR NOYAU
  // C  SIGMANU= SECTION EFFICACE TOTALE DE DIFFUSION COULOMBIENNE SUR NOYAU :
  // C    N(COLLISION)=DT*SIGMANU*EXP(-0.5(R/UTHERM)**2)
  // C    /(2PI*UTHERM**2)/(DISTANCE DES NOYAUX)

  double toto = exp(-1.0)/(_utherm*_utherm) + exp(EULER)/(atf*atf);
  _bcoh = 1.0/sqrt(toto);
  //    cout << " _bcoh " << _bcoh << " atf " << atf << endl;
  double sigmanu = 4.*M_PI*(zatom*EM_ALPHA)*(zatom*EM_ALPHA)*_bcoh*_bcoh;
  _rbohr2 = ( RBOHR1 * RBOHR1) /( zatom * zatom);

  // calculation of the decanalisation potential vlib
  // c 6 ----------------------------------------------------
  // c     le "point de decanalisation " (XLIB,YLIB) est situe sur la frontiere 
  // c     entre 2 cellules 
  // C     a cause de la soustraction de VLIB, le potentiel
  // c     tabule sera AUTOMATIQUEMENT NUL AU POINT (XLIB,YLIB) 
  // c    (DU MOINS POUR LE CAS UNIAXIAL)

  double r = hypot(_xlib, _ylib);

  Vdecanal decalib = formul(r);
  //    cout << " r= " << r << " vlib= " << decalib.V << endl;
  // tabulation loop

  double ffmax = 0.0;
  _zreco = 0.0;
  double dr = 0.0025; // discretization step for r
  Vdecanal decapot;
  int ir;
	
  for ( ir= 0; ir < LONG_POT; ir++)
    {
      r = dr * ir;
      decapot = formul(r);
      _potable[ir] = decapot.V - decalib.V;

      ffmax = max ( ffmax, fabs(decapot.F) ); // maximum field (eV/angstrom)
      _ftable[ir] = decapot.F/ELECTRON_MASS_EV; // changement d'unite eV --> masse(electron)
      // computation of electronic and nuclear densities

      double ether = 0.5 / (_utherm * _utherm);
      double rhonu = exp(-r*r*ether)*ether/(M_PI*_distatom);
      double rhoel = zatom*rhonu- decapot.deltaV/PI_4_HBAR_ALPHA;

      // c 7 --------------------------------------
      // c  RHONU = densite des noyaux 
      // c  RHOEL: densite des electrons
      // c                 = (densite des noyaux)/Z - (densite de charge totale)/e 
      // c
      // c  Densite de charge totale (electrons+noyaux) = DELTAV/e   
      // c
      // c                181 eV*Angstrom = e^2 = 4*M_PI*HBAR/137  


      if ( rhoel < 0.0 ) 
	{
	  cerr << " negative electronic density !, rhoel = " << rhoel << endl;
	}


      //	cout << " ir= " << ir << " rhonu " << rhonu << " sigmanu " << sigmanu << endl;
      // DCOLNU = inverse du libre parcour moyen de collision avec les noyaux.
      double dcolnu = rhonu*sigmanu;



      // c 9 --------------------------------
      // c    DCOLEL = inverse du libre parcour moyen de collision avec les electrons.
      // c    le cut-off en parametre d'impact de la collision e-e est suppose
      // c    etre environ max{R,Rbohr}.   
      // c    6.695D-4 = 4*M_PI/137**2 

      double dcolel = rhoel*PI_4_ALPHA_SQ*( r*r + _rbohr2);


      // c 10 -----------------------------------------------------------
      // c      DCOLTA(IR)=DCOLNU+DCOLEL  = "densite de collision" = inverse du l.p.m. 
      // c      de collision incoherente
      // c      Ici, on ne separe pas les collisions sur noyaux de celles sur electrons.
      // c       Dans "DIFMUL", on assumera que 
      // c  - si R < 2*UTHERM, les collisions sont sur des noyaux
      // c  - si R > 2*UTHERM, les collisions sont sur des electrons

      //	cout << " ir= " << ir << " dcolnu " << dcolnu << " dcolel " << dcolel << endl;
      _dcolta[ir] = dcolnu + dcolel;

      // c  11 - ZRECO = "Z reconstitue" = nombre d'electrons par atome situes
      // c      a l'interieur du cylindre de rayon R.
      // c      A la fin de la boucle DO, ZRECO doit etre egal a ZATOM

      _zreco += rhoel*2.0*M_PI*_distatom*r*dr;


      // c 12 - on montre les resultats pour r < 0.1 Angstrom et pour des r espaces
      // c      de 50 pas.

      if (r <= 0.1 || ir%50 == 0) 
	{
	  double courbu = 99.;
	  if (ir >= 2) 
	    {
	      courbu = _ftable[ir-1] - 0.5*( _ftable[ir] + _ftable[ir-2]);
	      courbu /= _ftable[ir-1];
	    }
	  //	    	    cout << r << " " << _potable[ir] << " " << _dcolta[ir] << " " << rhonu << " " << rhoel << " " << decapot.F << " " << courbu << endl;
	}

    }
  _barrier = _potable[0];
  //    cout << " BARRIER(eV)= " << _barrier << " FMAX(eV/angstrom)= " << ffmax << endl;
  //    cout << " disgtatom " << _distatom << " _zreco " << _zreco << endl;
}


// geometric parameters of the axial cell <111>
void Crystal::struc111(string typeOfNetwork, double zatom)
{
  _cristal = vector<double>(6);
  _cristal[0] = zatom;  // atomic number
  _cristal[1] = double(111.); 
  double oi;
  if ( typeOfNetwork == string("BCC") ) oi = _amaille/sqrt(6.0);
  else if ( typeOfNetwork == string("FCC") ) oi = _amaille/sqrt(24.0);
  else { cout << " Crystal::struc111 : typeOfNetwork " << typeOfNetwork << " not implemented " << endl; exit(1); }
  
  //  double ai = oi/SQROOT3;
  double ai = oi/sqrt(3.0);
  _cristal[2] = 2.0*oi;
  _cristal[3] = 1.0/_cristal[2];
  _cristal[4] = ai; // abscisse of B
  _cristal[5] = oi; // twice the ordinate of the point J


  // decanalisation point
  _xlib = 0.0;
  _ylib = oi;

  // C ------------------------------------------
  // C                                           |
  // C sheme           A   I   B                 | 
  // C of the                                    |  
  // C cell                      J               |
  // C <111>                                     |
  // C              C      O       D             |
  // C                                           |
  // C                                           |
  // C                                           |
  // C                 E       F                 |
  // C                                           |
  // C------------------------------------------
  // C  note : un point de decanalisation (par exemple I OU J) est un "col"
  // c  de potentiel sur la frontiere entre 2 cellules.

}

void Crystal::printGeometricParameters() const
{
  cout << " ----------------------------------------------- " << endl;
  cout << " geometric parameters : " << endl;
  for (size_t k=0; k < _cristal.size(); k++)
    {
      cout << " cristal[" << k << "] = " << _cristal[k] << endl;
    }
  cout << " ----------------------------------------------- " << endl;
}




Lindhard Crystal::calculationOfFields(double xx, double yy, double rayon, double charge) const {
  double ur=rayon/DR_LINDHARD;
  int ir=floor(ur);
  ur=ur-ir;
	
  if (ir>=(int)_potable.size() || rayon > 2.0*_cristal[4]){
    cerr<<" too big index ir in Crystal::calculationOfFields() "<< endl;
    cerr << " CRISTAL ";
    for (size_t k=0; k < _cristal.size(); k++) cerr << " " << _cristal[k];
    cerr << endl;
    throw string("poirot");
  }
	

  Lindhard lind;

  lind._epot=charge*(_potable[ir]+(_potable[ir+1] - _potable[ir])*ur);
  lind._f=charge*(_ftable[ir] + (_ftable[ir+1]-_ftable[ir])*ur);

  lind._fx = lind._f*xx/(rayon+0.1*1.0e-20);
  lind._fy = lind._f*yy/(rayon+0.1*1.0e-20);

  lind._dcoll = _dcolta[ir];

	
  //	Lindhard l(epot,f,fx,fy,dcoll);
  //	_lind.init(epot,f,fx,fy,dcoll);
  return lind;
	
}

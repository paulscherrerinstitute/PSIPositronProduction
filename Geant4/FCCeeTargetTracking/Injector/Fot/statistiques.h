#ifndef _STATISTIQUES_H
#define _STATISTIQUES_H



#include <vector>
#include <iostream>
#include <sstream>
#include "GlobalConstants.h"

using namespace std;


class statistiques {

 private:

  int _ndebut;	
  // int _numero;	
  int _nevnt;	
  int _mpas;	
  int _npas; 	
  int _mrept; 	
  int _nlmax;	
  int _nkuma; 	
  int _mfotk; 	
  int _nfotk;	
  int _ndifmu; 	
  int _mcoll; 	
  int _ncoll; 	
  int _neject; 	
  int _nbrems; 	
  int _mfotb; 	
  int _nfotb; 	
  int _nrept;	
  int _mfot;

	
  vector<int> _nproba; 
	
  double _fotk; 	
  double _wwk; 	
  double _wk; 	
  double _fotkou; 	
  double _wkou; 	
  double _favort; 	
  double _wavort;	
  double _fotb; 	
  double _wwb; 	
  double _wb;	
  double _fotbou; 	
  double _wbou; 	
  double _promax; 
	
  vector< vector<double> >  _effic; 
	
	
 public:
	
  statistiques();
  ~statistiques();
	

  void incrementNeject()
  {
    _neject++;
  } 

  void incrementMfot() 
  {
    _mfot++;
    if ( _mfot > MAX_PHOTONS )
      {
	cerr << " statistiques : to much photons MFOT   > " << MAX_PHOTONS << endl;
	throw string("poirot");
      }
  }


  void addEvent(); 

  void addEnergyKuma(double proba, double epho, int IHU, int IHD, double AA); 

  void addPhotKuma(double proba, double epho);


void addEnergyBrems(double proba, double epho)
{
  _wb += proba * epho;
  _fotb += proba;
}

  void addPhotBrems(double proba, double epho);
void addDiffMul(int ncol)
{
  _mcoll += ncol;
  _ndifmu++;
}

	
	
	
void incrementMrept() {
	
	_mrept++;
}
	             	
void incrementMpas() {
	
	_mpas++;
}
	

void incrementNlmax() {	
	_nlmax++;
}
	
	
void incrementNkuma() {	
	_nkuma++;
}
	
	
void incrementNbrems() {
	
	_nbrems++;
}
	
	
       
	
void incrementNevnt() {
	
	_nevnt++;
}
	

	
	
	

void incrementAvort( double poid, double epho) {
	
	_favort += poid;
	_wavort += poid*epho;
}

	
  void zeros();
	
  void tsatis();
	
  void reset();
	
	
  void printPoirot() const	
  {
    cout << " MPAS=  " << _mpas << " MREPT= " << _mrept  << " MCOLL= " << _mcoll << endl;
    cout << " MFOTB(bremse)= " << _mfotb << " MFOTK(kumakhov)= " << _mfotk << endl;
  }
	
};
#endif

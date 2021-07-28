/*
 *  statistiques.cpp
 *  
 *
 *  Created by Berte on 02/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include <math.h>

#include "statistiques.h"



using namespace std;

statistiques::statistiques() {
	
	_nproba = vector<int>(11);
	_effic = vector< vector<double> >(11);
	
	int k;
	for ( k=0; k < 11; k++) 
    {
		_effic[k] = vector<double>(11);
    }
	
	
	//initialisation des variables qui sont dans la principale
	_ndebut = 1;
	
	_nkuma = 0;
	_fotk = 0;
	_wk = 0;
	_promax = 0.0;
	_favort = 0.0;
	_wavort = 0.0;
	_mcoll = 0;
	_neject = 0;
	
}


statistiques::~statistiques() {
}


void statistiques::addEnergyKuma(double proba, double epho, int IHU, int IHD, double AA) 
{
  _fotk += proba;
  _wk += proba*epho;  // [es8]

  _effic.at(IHU).at(IHD) += proba;
  //_effic[IHU][IHD] += proba;
 _promax = max(proba, _promax);

 if (proba != 0.0 ) 
   {
     int klocal = 5 + round(log(proba) +0.5 );
     klocal = max(1, min(klocal, 10) );
     _nproba[klocal]++;
   }

 //  PHOTONS REJETES POUR PHASE DEBILE ( AVANT TEST etat final )  [es9]

 if(AA < 0.0) _nproba[0]++;
}


void statistiques::addPhotKuma(double proba, double epho)
{
  _mfotk++;
  _wwk += epho;     // [es11]
  //photons oublies
  double aux = max(proba, 1.0) - 1.0;
  _fotkou += aux;
  _wkou += aux * epho;   // [es12]
}



void statistiques::addPhotBrems(double proba, double epho)
{
  _mfotb++;
  _wwb += epho;
  if (proba > 1.0){
    double aux = proba - 1.0;	
    _fotbou += aux;
    _wbou += aux * epho;
  }
}



void statistiques::addEvent() 
{
	_npas += _mpas;
	_nrept += _mrept;
	_ncoll += _mcoll;
	_nfotk += _mfotk;
	_nfotb +=  _mfotb;
}







void statistiques::tsatis() 
{
	
	std::cout << "pas: " << _npas << " rept: " <<  _nrept <<  " DIFMU " << _ndifmu << " collis " << _ncoll << std::endl;
	
	cout << "n(lmax): " << _nlmax << " ejectes: " <<  _neject << endl;
	


	cout << "appels BREMS: " << _nbrems << " reussis: " <<  _nfotb << " " << _fotb << " Wbrems " << _wwb << " " << _wb << endl;

	
	//photons brems oublies
	cout << "photons brems. oublies: " << _fotbou << " energie= " <<  _wbou << endl;
	
	cout << "appels KUMA: " << _nkuma << " reussis: " <<  _nfotk << " " << _fotk << " Wkuma " << _wwk << " " << _wk << endl;
	
	
	//photons kuma oublies
	cout << "photons kuma. oublies: " << _fotkou << " energie= " <<  _wkou << endl;
	
	cout <<  _nproba[0] << " rejets pour phase débile " << endl;
	
	cout << " avortes dans ETAFIN: " << _favort << " energie= " <<  _wavort << endl;
	
	
	//histogramme des proba des photons kuma
	
	cout << "hist. proba(kuma) < 1: " << endl;
	  cout <<  _nproba[1] << " " << _nproba[2] << " " << _nproba[3] << " " << _nproba[4] << " " << _nproba[5] << endl;
	
	  cout << "hist. proba(kuma) > 1: "  << endl;
	  cout <<  _nproba[6] << " " << _nproba[7] << " " << _nproba[8] << " " << _nproba[9] << " " << _nproba[10] << endl;
	
	cout << "promax: " << _promax << endl;
	
	
	//tableau efficacu-ite de kuma
	
	
	cout << "efficacite du tirage en U et GAT dans kuma" << endl;
	
	double toto = 100.0000 / (_fotk + 1e-25);
	
	cout << "vers le haut: GAT croissant " << endl; 
	
	cout << "vers la droite: U croissant. Derniere col. : cumul/U " << endl; 
	
	ostringstream s;
	
	ostringstream ss; 
	
	for (int ihd =0; ihd < 10; ihd++){
		
		for (int ihu =0; ihu < 10; ihu++){

			
			_effic[ihu][ihd] = _effic[ihu][ihd] * toto;

			//somme sur GAT
			_effic[ihu][10] = _effic[ihu][10] + _effic[ihu][ihd] * 0.1;

			
			//somme sur U 
			_effic[10][ihd] = _effic[10][ihd] + _effic[ihu][ihd] * 0.1;

			
		}
		
		for (int ihu =0; ihu < 10; ihu++){
			
			s << " " <<  _effic[ihu][ihd];
		}
		
		cout << s.str() << ":" <<  _effic[10][ihd] << endl;
		
		s.str("");
		
		_effic [10][10] = _effic[10][10] + _effic[10][ihd] * 0.1;
		
	}
	
	cout << "cumul sur GAT: " << endl; 
	
	for (int ihu =0; ihu < 10; ihu++){
		
		ss << " " << _effic[ihu][10];
	}
	
	cout << ss.str() << " " << _effic[10][10] << endl;
	
	ss.str("");
	
	
	
}


void statistiques::zeros() {
	
	for (int i1 = 0; i1 < 11; i1++){
		
		for (int i2 = 0; i2 < 11; i2++){
			
			_effic[i1][i2] = 0;
		}
	}
	
	for (int i1 = 0; i1 < 11; i1++){
		
		_nproba[i1]=0;
	}
	
	_npas = 0;
	_nrept = 0;
	_ncoll = 0;
	_nevnt = 0;
	_ndifmu = 0;
	_nkuma = 0;
	_nbrems = 0;
	_fotk = 0;
	_wk = 0;
	_nfotk = 0;
	_wwk = 0;
	_fotkou = 0;
	_wkou = 0;
	_fotbou = 0;
	_wbou = 0;
	_favort = 0;
	_wavort = 0;
	_fotb = 0;
	_wb = 0;
	_nfotb = 0;
	_wwb = 0;
	_promax = 0;
	_nlmax = 0;
	_neject = 0;
}

// mise a zero des statistiques   [p6]
void statistiques::reset(){
	_mpas = 0;
	_mrept = 0;
	_mcoll = 0;
	_mfotb = 0;
	_mfotk = 0;
	_mfot = 0;
}
	
// c p6 ------------------------------------------------------
// c   Nuance entre "Mtoto" et "Ntoto" :  
// c  "Mtoto" compte les "totos" d'un meme evenement, 
// c  tandis que  "Ntoto" compte les "totos" en cumulant tous les evenements.

// es8 - WK=energie KUMA "CONTINUE", avant test ETAFIN (evenements cumules).

// c es9 - Si ces photons sont trop nombreux, e.g. NPROBA(0)/NKUMA>0.5, 
// c        augmenter PHOMIN.

// c es11 -  WWK: energie KUMA discrete (evenements cumules)

// c es12 -  maintenir WKOU << WK, FOTKOU << FOTK, au besoin augmenter FREKUMA

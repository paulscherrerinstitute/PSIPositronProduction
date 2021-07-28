
#include "Fot.h"
#include "Evenement.h"

#ifdef TEST_WCLASS
#include "wclass.h"
#endif

Fot::Fot(RunParameters& rp) : _snak(NULL),  _partCrys(NULL), _runPar(rp)
{
#ifdef DO_STATS
 _stat = new statistiques;
_stat->zeros();
#else
 _stat = NULL;
#endif

      _zexit = _runPar.getZexit();
      _runPar.getEjectionCut_off(_etmax, _vtmax);
 _nevnt  = 0;

 
}

void Fot::makeKumakhov(ParticleCollection& partColl)
{

  cout<<"run makeKumakhov "<<endl;

  try
    {
#ifdef DO_STATS
      _stat->zeros();	
#endif
      const Particle* part = partColl.initParticleIteration();      
      // BremsStrahlung* bremse = NULL;

      _nevnt  = 0;

#ifdef TEST_WCLASS
      double gamma;
      WCLASS = 0.0; // TEST REGLE DE SOMME
#endif


      while(part != NULL)
	{
	  _nevnt++;
	  makeSingleParticleKumakhov(part);
	  part = partColl.getNextParticle();
	} //   fin boucle particule


      _photons.saveOnFile("photons.dat");

#ifdef TEST_WCLASS
      WCLASS *= gamma* gamma *0.96e-8;  //  cas du test "W classique"  ![p19]
      cout << "  Wclassique " << WCLASS << endl;
#endif
	
      finir();
      cout << " nombre de photons " << _photons.getNbPhotons() << endl;
    }

  catch (string erreur)
    {
      if ( erreur == string("poirot") ) poirot();
      else
	{
	  cerr << " error " << erreur << " not managed " << endl;
	}
    }
}


const ParticleInCrystal& Fot::makeSingleParticleKumakhov(const Particle* part)
{


//   try
//     {

      BremsStrahlung* bremse = NULL;


#ifdef TEST_WCLASS
      double gamma;
      WCLASS = 0.0; // TEST REGLE DE SOMME
#endif



//       double etmax, vtmax, zexit;
//       zexit = _runPar.getZexit();
//       _runPar.getEjectionCut_off(etmax, vtmax);

      if ( _partCrys != NULL ) delete _partCrys;

      _partCrys = new ParticleInCrystal(_runPar.getCrys(),*part);

      _snak = new Snake(_runPar, *_partCrys, _photons, _stat);
#ifdef DO_STATS
      _stat->incrementNevnt();
      _stat->reset();
#endif
#ifndef TEST_WCLASS
      bremse = new BremsStrahlung( _photons, _runPar.getPhomin(), _runPar.getPoimin(), _stat);
#endif
      Evenement eve( _snak, bremse, &_photons, _etmax, _vtmax, _zexit);
      bool suite = true;

      // calcul proprement dit
      while (suite) 
	{
	  suite =   eve.makeStep();  
	} // fin boucle 30

#ifdef DO_STATS
      _stat->addEvent();
#endif
      delete _snak;
      _snak = NULL;

#ifdef TEST_WCLASS
      gamma = _partCrys->getGamma();
#endif

      if ( bremse != NULL ) delete bremse;

#ifdef TEST_WCLASS
      WCLASS *= gamma* gamma *0.96e-8;  //  cas du test "W classique"  ![p19]
      cout << "  Wclassique " << WCLASS << endl;
#endif
	
      //      cout << " nombre de photons " << _photons.getNbPhotons() << endl;
      return *_partCrys;
//     }
//   catch (string erreur)
//     {
//       if ( erreur == string("poirot") ) poirot();
//       else
// 	{
// 	  cerr << " error " << erreur << " not managed " << endl;
// 	}
//     }
}



void Fot::poirot() 
{
  cout << " ***************************************************************** " << endl;
  cout << " FOTPP Blague!... POIROT :   " << endl;
  _snak->printPoirot();
#ifdef DO_STATS
  _stat->printPoirot();
#endif
  _partCrys->printPoirot();
  cout << " EPOT= " << _lind._epot << " FX= " << _lind._fx << " FY= " << _lind._fy  << " F= " << _lind._f  << endl;
  cout << " ***************************************************************** " << endl;
  finir();
}






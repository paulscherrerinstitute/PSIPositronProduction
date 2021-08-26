#include "ParticleCollection.h"
#include <fstream>
#include "GlobalConstants.h"
#include "mathematics.h"
#include <iostream>



ParticleCollection::ParticleCollection(double einit, int ntir, double xmx, 
   double xmy, double xmvx, double xmvy, double sx, 
    double sy, double svx, double svy) {
	double charge;
	
	for (int i=1; i <= ntir; i++){
	  double x = mathematics::gauss(xmx , sx);
		double y = mathematics::gauss(xmy , sy);
		double z = 0.0;
		bool lposit = false;
		double vx = mathematics::gauss(xmvx , svx);
		double vy = mathematics::gauss(xmvy , svy);
		
		double gamma = einit/ ELECTRON_MASS_GEV;
		
		double px = gamma*vx;
		double py = gamma*vy;
		if ( lposit) charge = 1.0;
		else charge = -1.0;
		_particles.push_back(Particle(charge,x,y,z,px,py,gamma));
	}
	_currentIterator = _particles.begin();
}



void ParticleCollection::fill(string filename, int nbegin, int nend)
{
  ifstream input_file;
  int numero, lposit;
  double x,y,z,EGeV,vx,vy, gamma;
  double px, py, charge;
  input_file.open(filename.c_str());
  if ( !input_file )
    {
      cout << " ParticleCollection::fill : could not open the file " << filename << endl;
    }
  //  cout << " numero " << numero << endl;
  while ( input_file >> numero >>  x >> y >> z >> EGeV >> vx >> vy >> lposit)
    {
      if (numero >= nbegin && numero <= nend) 
	{

	  gamma = EGeV/ ELECTRON_MASS_GEV;

	  px = gamma*vx;
	  py = gamma*vy;
	  if ( lposit) charge = 1.0;
	  else charge = -1.0;
	  _particles.push_back(Particle(charge,x,y,z,px,py,gamma));
	}
    }
  input_file.close();
  _currentIterator = _particles.begin();
}


const Particle* ParticleCollection::getNextParticle(){
  if(_currentIterator != _particles.end()){
    const Particle* p = &(*_currentIterator);
    _currentIterator++;
		
    return p;
  }
  else{
    return NULL;
  }
}


void ParticleCollection::printParticles() const
{
  list<Particle>::const_iterator it;
  // double x,y,z,EGeV,vx,vy, gamma;
  // double px, py;
  int compteur = 1;
  for (it = _particles.begin(); it != _particles.end(); it++)
    {
      std::cout <<"particle n° "<< compteur <<std::endl;
      it->printParticle();
      compteur++;
    }
}


void ParticleCollection::saveOnFile(string filename)
{
  ofstream output_file;
  int numero, lposit;
  // double px, py, charge;
  output_file.open(filename.c_str());
  if ( !output_file )
    {
      cout << " ParticleCollection::saveOnFile : could not open the file " << filename << endl;
    }

  list<Particle>::const_iterator itr;
  numero = 1;
  for ( itr = _particles.begin(); itr != _particles.end(); itr++ ) 
    {
      double x,y,z,EGeV,vx,vy, gamma;
      x = itr->getXPosition();
      y = itr->getYPosition();
      z = itr->getZPosition();
      gamma = itr->getGamma();
      EGeV = gamma * ELECTRON_MASS_GEV;
      vx = itr->getPx()/gamma;
      vy = itr->getPy()/gamma;
      if ( itr->isPositron() ) lposit = 1;
      else lposit = 0;
      output_file <<  numero << "  " <<  x << "  " << y << "  " << z  << "  " <<  EGeV << "  " << vx << "  " << vy << "  " <<  lposit << endl;
      numero++;
    }

  output_file.close();
}



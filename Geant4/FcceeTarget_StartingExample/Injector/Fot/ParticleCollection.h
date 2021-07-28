#ifndef _PARTICLECOLLECTION_H
#define _PARTICLECOLLECTION_H

#include "Particle.h"

#include <list>
#include <math.h>

/**
 * \file ParticleCollection.h
 * \brief The ParticleCollection class provides a useful way to create collection of particles
 * \author Guy LE MEUR & Sandrine BERTE
 * \date 01.09.2010
 */



/** \class ParticleCollection
 *
 *  The ParticleCollection class is the object that will contain particles set randomly
 */
class ParticleCollection
{

 private : 

  list<Particle> _particles;
  list<Particle>::const_iterator _currentIterator;
	
	
public :

	
	/**
     *  \brief Constructor
     *
     *  Default constructor of a collection of particles
     *
     */

 ParticleCollection() { _currentIterator = _particles.begin(); }
	
	/**
     *  \brief Constructor
     *
     *  Construct a collection of particles from data 
     *
     *  \param double :
	 *  \param double : nombre de tirs ntir
	 *  \param double :
	 *  \param double :
	 *  \param double :
	 *  \param double :
	 *  \param double :
	 *  \param double :
	 *  \param double :
	 *  \param double :
     *
     */
  ParticleCollection(double einit, int ntir, double xmx, double xmy, double xmvx, double xmvy, double sx, double sy, double svx, double svy);
  
	
	/**
     *  \brief Destructor
     *
     *  Destroy the collection of particles 
     *
     *
     */
  ~ParticleCollection()  {;}
	
	
	/**
     *  \brief return the lenght of the collection of particles
     *
     *  return the number of particles in the collection
     *
     * \return int: return the number of particles
     */

  void addParticle(double charge, double x, double y, double z, double px, double py, double gamma)
  {
    _particles.push_back(Particle(charge,x,y,z,px,py,gamma));
  }


int lenght(){
  return _particles.size();
}
	
	
	/**
     *  \brief fill a file
     *
     *  Fill a file with data from particle number nbegin-nend
     *
     *  \param string : the name of file
	 *  \param int : the number of the first particle in the file
	 *  \param int : the number of the last particle in the file
     */
  void fill(string filename, int nbegin, int nend);
	
	
	/**
     *  \brief return the next particle
     *
     *  return a pointer to the next particle in the collection
     *
     * \return Particles*: return the next particle
     */
  const Particle* getNextParticle();
  const Particle* initParticleIteration()
  {
    if ( _particles.size() != 0 ) 
      {
	_currentIterator = _particles.begin();
	const Particle* p =  &(*_currentIterator);
	_currentIterator++;
	return p;
      }
    else return NULL;    
  }

	
	
  Particle getFirstParticle()
  {
    return  *_particles.begin();
  }

	/**
     *  \brief Print all particles which are in the collection
     *
     *  Edit data concerning all particles of the collection 
	 *
     */
  void printParticles() const;
  void saveOnFile(string filename);

};

#endif

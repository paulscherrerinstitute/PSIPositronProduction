
#include "Particle.h"

Particle::Particle(const Particle& part)
{
  _charge = part._charge;
  _xPosition=part._xPosition;
  _yPosition=part._yPosition;
  _rayon = hypot(_xPosition, _yPosition);
  _zPosition=part._zPosition;
  _px=part._px;
  _py=part._py;
  _gamma=part._gamma;
  _isPosit = part._isPosit;
}



Particle::Particle(double ch, double x, double y, double z, double px, double py, double g) {
  _charge = ch;
  _xPosition=x;
  _yPosition=y;
  _rayon = hypot(_xPosition, _yPosition);
  _zPosition=z;
  _px=px;
  _py=py;
  _gamma=g;
  if ( _charge > 0.0 ) _isPosit = true;
  else _isPosit = false;
	//cout << "::::::))))))  "  <<  _xPosition <<  "  " << _yPosition << endl;
}

// double Particle::getCharge() const
// {
//   return _charge; 
// }





void Particle::printParticle() const
{
  double EGeV = _gamma * ELECTRON_MASS_GEV;
  // double vx = _px/_gamma;
  // double vy = _py/_gamma;
  int lposit= _charge > 0.0;
  printf(" x= %20.17e y= %20.17e z= %20.17e \n", _xPosition,  _yPosition, _zPosition);
  std::cout <<  " egev: " << EGeV << " px: " << _px << "py: " << _py << " lposit= "<< lposit << std::endl;
}




ParticleInCrystal::ParticleInCrystal(const Crystal& crys, const Particle& part) : Particle(part), _iemis(false), _ientree(true), _iexit(false),  _x0cell(0.0), _y0cell(0.0), _crys(crys), _bigJump(false)
{_nsaut = 0;}

ParticleInCrystal::ParticleInCrystal(const Crystal& crys) : Particle(),  _iemis(false), _ientree(true), _iexit(false), _x0cell(0.0), _y0cell(0.0), _crys(crys), _bigJump(false)
{_nsaut = 0;}


//Vraie si la particule a migre vers une autre cellule
bool ParticleInCrystal::migration() {
  double hmax = (_crys.getCristal())[2];
  double h = _xPosition*SQROOT3;
  double h1 = _yPosition + h;
  double h2 = _yPosition - h;
  double h3 = - _yPosition  - _yPosition;
  h = max(max(fabs(h1),fabs(h2)),fabs(h3));
  bool mig = h > hmax;
	
	
  if(!mig) 
    return false;
  else
    {
      _nsaut++;
      // _bigJump corespond a la variable miggg de fortran que l'on actualise dans migre
      _bigJump = h > 2*hmax;
      updateCell(h1, h2, h3);
      return true;
    }
}



  void ParticleInCrystal::reIniCoorAndCell(double x, double y, double xCell, double yCell)
{
  _xPosition = x;
  _yPosition = y;
  _rayon = hypot(_xPosition, _yPosition);
  setX0cell(xCell);
  setY0cell(yCell);
}




//met a jour les coordonnees de la cellule dans laquelle la particle se trouve en cas de migration
void ParticleInCrystal::updateCell(double h1, double h2, double h3) {
  int ithree=3;
	
	
  int n1=round(h1*(_crys.getCristal())[3]+0.5);
  int n2=round(h2*(_crys.getCristal())[3]+0.5);
  int n3=round(h3*(_crys.getCristal())[3]+0.5);
  n3=n1+n2+n3-1;
  int n12=(n1-n2)%ithree;
	
  if (n12 == 0){
    n1=n1-n3;
    n2=n2-n3;
  }
  else{
    if(n12 == 1 || n12 == -2) {
      n1=n1-1;
    }
    else{
      n2=n2-1;
    }
  }
	
  double h=(n1-n2)*(_crys.getCristal())[4];
  _xPosition=_xPosition-h;
  _x0cell=_x0cell+h;
  h=(n1+n2)*(_crys.getCristal())[5];
  _yPosition=_yPosition-h;
  _y0cell=_y0cell+h;
  _rayon = hypot(_xPosition, _yPosition);
}

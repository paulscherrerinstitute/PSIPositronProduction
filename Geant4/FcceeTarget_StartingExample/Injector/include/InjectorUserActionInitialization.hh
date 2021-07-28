#ifdef G4MULTITHREADED

#ifndef InjectorUserActionInitialization_h
#define InjectorUserActionInitialization_h 1

#endif

#include "G4VUserActionInitialization.hh"

class G4GeneralParticleSource;
class InjectorDetectorConstruction;

class InjectorUserActionInitialization : public G4VUserActionInitialization{

public:

  InjectorUserActionInitialization(InjectorDetectorConstruction*);

  ~InjectorUserActionInitialization();

  void Build() const;
  void BuildForMaster() const;

private:

  G4GeneralParticleSource* masterGPS;

  InjectorDetectorConstruction* fDetector;

};

#endif


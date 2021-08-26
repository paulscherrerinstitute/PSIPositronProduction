#ifdef G4MULTITHREADED

#include "InjectorUserActionInitialization.hh"
#include "InjectorPrimaryGeneratorAction.hh"
#include "InjectorTrackingAction.hh"
#include "InjectorSteppingAction.hh"
#include "InjectorEventAction.hh"
#include "InjectorRunAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"

InjectorUserActionInitialization::InjectorUserActionInitialization(InjectorDetectorConstruction* det)
:fDetector(det){
}

InjectorUserActionInitialization::~InjectorUserActionInitialization() {
}

void InjectorUserActionInitialization::Build() const {
    SetUserAction(new InjectorPrimaryGeneratorAction(fDetector));
    SetUserAction(new InjectorEventAction());
    SetUserAction(new InjectorSteppingAction(fDetector));
    SetUserAction(new InjectorTrackingAction());
    SetUserAction(new InjectorRunAction());
}

void InjectorUserActionInitialization::BuildForMaster() const {
}

#endif


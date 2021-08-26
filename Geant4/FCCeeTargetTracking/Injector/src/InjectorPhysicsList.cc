#include "InjectorPhysicsList.hh"

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessVector.hh"
#include "G4RunManager.hh"

#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"

InjectorPhysicsList::InjectorPhysicsList():  G4VModularPhysicsList(){

    fTimeStepMin = 2.E2 * CLHEP::angstrom;
    fTransverseVariationMax = 2.E-2 * CLHEP::angstrom;    

    fEmPhysicsList = new G4EmStandardPhysics();
    fParticleList  = new G4DecayPhysics();

}

InjectorPhysicsList::~InjectorPhysicsList(){
    delete fEmPhysicsList;
    delete fParticleList;
}

void InjectorPhysicsList::ConstructProcess(){
    
    AddTransportation();

    fEmPhysicsList->ConstructProcess();
}

void InjectorPhysicsList::ConstructParticle(){
    fParticleList->ConstructParticle();
}

void InjectorPhysicsList::SetCuts()
{
    // These values are used as the default production thresholds for the world volume.
    SetCutsWithDefault();
}


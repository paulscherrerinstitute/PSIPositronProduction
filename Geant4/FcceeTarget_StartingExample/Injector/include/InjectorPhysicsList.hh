#ifndef InjectorPhysicsList_h
#define InjectorPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class InjectorPhysicsList: public G4VModularPhysicsList
{

private:

    G4VPhysicsConstructor*  fEmPhysicsList;
    G4VPhysicsConstructor*  fParticleList;
    
public:

    G4double GetTransverseVariationMax() {return fTransverseVariationMax;};
    void     SetTransverseVariationMax(G4double aDouble) {fTransverseVariationMax = aDouble;};

    G4double GetTimeStepMin() {return fTimeStepMin;};
    void     SetTimeStepMin(G4double aDouble) {fTimeStepMin = aDouble;};
    
public:
    
    InjectorPhysicsList();
    ~InjectorPhysicsList();
    
    // Construct particle and physics
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();
    
    G4double fTimeStepMin;
    G4double fTransverseVariationMax;
};

#endif

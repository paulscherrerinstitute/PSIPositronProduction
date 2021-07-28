#ifndef InjectorPrimaryGeneratorActionMessenger_h
#define InjectorPrimaryGeneratorActionMessenger_h 1

#include "G4UImessenger.hh"

class InjectorPrimaryGeneratorAction;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

class InjectorPrimaryGeneratorActionMessenger: public G4UImessenger{

  public:
    InjectorPrimaryGeneratorActionMessenger(InjectorPrimaryGeneratorAction*);
    virtual ~InjectorPrimaryGeneratorActionMessenger();
    virtual void SetNewValue(G4UIcommand*, G4String);

  private:
    InjectorPrimaryGeneratorAction* fPrimaryGeneratorAction;
    //G4UIcmdWithADoubleAndUnit* fXtalThickCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimaryEnergyCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaECmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaXYCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaZCmd; 
    //G4UIcmdWithADoubleAndUnit* fPrimarySigmaTCmd; 
    G4UIcmdWithADouble*        fPrimarySigmaTCmd; 
    G4UIcmdWithADoubleAndUnit* fPrimarySigmaPxyCmd; 
};

#endif


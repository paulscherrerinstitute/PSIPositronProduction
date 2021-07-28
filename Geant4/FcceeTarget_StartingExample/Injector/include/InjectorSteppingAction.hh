#ifndef InjectorSteppingAction_h
#define InjectorSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class InjectorDetectorConstruction;

class InjectorTrackingAction;

class G4VPhysicalVolume;

class InjectorSteppingAction : public G4UserSteppingAction
{
  public:
    InjectorSteppingAction(InjectorDetectorConstruction*);
   ~InjectorSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  private:
    
    InjectorDetectorConstruction* fDetector;
    InjectorTrackingAction*       fTrackAction;        
    
    G4bool                first; 

    G4VPhysicalVolume*    fWorld;
    G4VPhysicalVolume*    fXtal;
    G4VPhysicalVolume*    fAmorphous;
    
    G4ThreeVector         fPosition;        
    G4ThreeVector         fMomentum;        
    
    G4double 		  energy;
    G4double              time;
};

#endif

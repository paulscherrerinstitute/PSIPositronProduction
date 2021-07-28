#ifndef InjectorTrackingAction_h
#define InjectorTrackingAction_h 1

#include "G4UserTrackingAction.hh"

class InjectorTrackingAction : public G4UserTrackingAction{
    
public:
    
    void PreUserTrackingAction(const G4Track*);
   
    void PostUserTrackingAction(const G4Track*);
};

#endif


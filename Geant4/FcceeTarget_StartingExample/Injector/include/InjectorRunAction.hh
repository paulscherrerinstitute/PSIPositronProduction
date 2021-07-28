#ifndef InjectorRunAction_h
#define InjectorRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class InjectorRunActionMessenger;

class InjectorDetectorConstruction;

class InjectorRunAction : public G4UserRunAction
{

public:

    InjectorRunAction();
    virtual ~InjectorRunAction();

    InjectorRunActionMessenger* fRunActionMessenger;

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    void SetFileName(G4String);

private:

    G4String fFileName;

};

#endif


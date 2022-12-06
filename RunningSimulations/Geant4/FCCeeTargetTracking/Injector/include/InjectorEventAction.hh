#ifndef InjectorEventAction_h
#define InjectorEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>

class InjectorEventActionMessenger;

class InjectorEventAction : public G4UserEventAction
{
  public:
    InjectorEventAction();
    virtual ~InjectorEventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private:
    //D/G4int fSD_ID;
    G4int fVerboseLevel;

  public:
    inline void SetVerbose(G4int val) { fVerboseLevel = val; }
    inline G4int GetVerbose() const { return fVerboseLevel; }
};

#endif

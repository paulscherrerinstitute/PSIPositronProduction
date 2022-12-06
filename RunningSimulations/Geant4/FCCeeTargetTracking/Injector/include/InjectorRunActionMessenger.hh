#ifndef InjectorRunActionMessenger_h
#define InjectorRunActionMessenger_h 1

#include "G4UImessenger.hh"

class InjectorRunAction;
class G4UIcmdWithAString;

class InjectorRunActionMessenger: public G4UImessenger{

  public:

    InjectorRunActionMessenger(InjectorRunAction*);

    virtual ~InjectorRunActionMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    InjectorRunAction*  fRunAction;

    G4UIcmdWithAString* fFileNameCmd; 
};

#endif


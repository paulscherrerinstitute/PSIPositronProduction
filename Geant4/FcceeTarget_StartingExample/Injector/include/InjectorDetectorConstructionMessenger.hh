#ifndef InjectorDetectorConstructionMessenger_h
#define InjectorDetectorConstructionMessenger_h 1

class InjectorDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAIntAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

#include "G4UImessenger.hh"
#include "globals.hh"

class InjectorDetectorConstructionMessenger: public G4UImessenger
{
  public:
    InjectorDetectorConstructionMessenger(
                    InjectorDetectorConstruction* mpga);
    ~InjectorDetectorConstructionMessenger();

    virtual void SetNewValue(G4UIcommand * command,G4String newValues);
    //virtual G4String GetCurrentValue(G4UIcommand * command);

  private:
    InjectorDetectorConstruction * fDet;

    G4UIdirectory* fMyXtalDirectory;
    //G4UIcmdWithAString*  fXtalMaterialCmd;
    G4UIcmdWith3VectorAndUnit* fXtalSizeCmd;
    G4UIcmdWith3VectorAndUnit* fAmorphousSizeCmd;
    G4UIcmdWith3VectorAndUnit* fAmorphousDistanceCmd;
    G4UIcmdWith3VectorAndUnit* fDipoleFieldCmd;
    G4UIcmdWith3VectorAndUnit* fAMDFrontGapSizeCmd;
    G4UIcmdWithAString*        fAMDOptionCmd;
    G4UIcmdWithADoubleAndUnit* fAMDR1minCmd;
    G4UIcmdWithADoubleAndUnit* fAMDR1maxCmd;
    G4UIcmdWithADoubleAndUnit* fAMDR2minCmd;
    G4UIcmdWithADoubleAndUnit* fAMDR2maxCmd;
    G4UIcmdWithADoubleAndUnit* fAMDLengthCmd;
    G4UIcmdWithABool*          fAMDLinFringeActCmd;
    G4UIcmdWithADoubleAndUnit* fAMDLinFringeB0Cmd;
    G4UIcmdWithADouble* fAMDLinFringeKCmd;
    G4UIcmdWithADoubleAndUnit* fAMDLinFringeB0ZPosCmd;
    G4UIcmdWithADoubleAndUnit* fAMDLinFringeTargZPosCmd;
};

#endif


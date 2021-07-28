#ifndef InjectorDetectorConstruction_h
#define InjectorDetectorConstruction_h 1

#include <TH2F.h>
//#include <TGraph2D.h>

#include "G4VUserDetectorConstruction.hh"
#include "InjectorDetectorConstructionMessenger.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "G4MagneticField.hh"
#include "G4ElectricField.hh"

#include "globals.hh"

class TH2F;
//class TGraph2D;

class G4UniformMagField;

class InjectorRunAction;

class InjectorDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    
    InjectorDetectorConstruction();
    ~InjectorDetectorConstruction();
    
    void DefineMaterials();
    G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

    TH2F* GetHistFieldMapR() {return histFieldMapR;}
    TH2F* GetHistFieldMapZ() {return histFieldMapZ;}
    //TGraph2D* GetGraphFieldMapR() {return graphFieldMapR;}
    //TGraph2D* GetGraphFieldMapZ() {return graphFieldMapZ;}

private:
    InjectorDetectorConstructionMessenger* fMessenger;
    InjectorRunAction* run;

private:
    G4ThreeVector fWorldSize;
    G4Box* fWorldSolid;
    G4LogicalVolume* fWorldLogic;
    G4VPhysicalVolume* fWorldPhysical;
    G4Material* fWorldMaterial;

    TH2F* histFieldMapR;
    TH2F* histFieldMapZ;
    //TGraph2D* graphFieldMapR;
    //TGraph2D* graphFieldMapZ;

public:

    G4ThreeVector GetWorldSize(){return fWorldSize;}

    void SetXtalMaterial(const G4String& name);
    G4String GetXtalMaterial();
    void SetXtalSize(G4ThreeVector);
    G4ThreeVector GetXtalSize() {return fXtalSize;};
    G4ThreeVector GetXtalLocation(){return fXtalLocation;}
    
    void SetAmorphousSize(G4ThreeVector);
    G4ThreeVector GetAmorphousSize() {return fAmorphousSize;};
    G4ThreeVector GetAmorphousLocation() {return fAmorphousLocation;};
    void SetAmorphousDistance(G4ThreeVector);
    G4ThreeVector GetAmorphousDistance() {return fAmorphousDistance;};
    void SetDipoleFieldValue(G4ThreeVector);
    G4ThreeVector GetDipoleFieldValue() {return fDipoleFieldValue;};

    void SetAMDOption(const G4String& name);
    void SetAMDFrontGapSize(G4ThreeVector);
    void SetAMDR1min(G4double);
    void SetAMDR1max(G4double);
    void SetAMDR2min(G4double);
    void SetAMDR2max(G4double);
    void SetAMDLength(G4double);
    void SetAMDLinFringeAct(G4bool);
    void SetAMDLinFringeB0(G4double);
    void SetAMDLinFringeK(G4double);
    void SetAMDLinFringeB0ZPos(G4double);
    void SetAMDLinFringeTargZPos(G4double);
    G4String GetAMDOption();
    G4ThreeVector GetAMDFrontGapSize() {return fAMDFrontGapSize;};
    G4ThreeVector GetAMDFrontGapLocation() {return fAMDFrontGapLocation;};
    G4double GetAMDR1min() {return fAMDR1min;};
    G4double GetAMDR1max() {return fAMDR1max;};
    G4double GetAMDR2min() {return fAMDR2min;};
    G4double GetAMDR2max() {return fAMDR2max;};
    G4double GetAMDLength() {return fAMDLength;};
    G4bool   GetAMDLinFringeAct() {return fAMDLinFringeAct;};
    G4double GetAMDLinFringeB0() {return fAMDLinFringeB0;};
    G4double GetAMDLinFringeK() {return fAMDLinFringeK;};
    G4double GetAMDLinFringeB0ZPos() {return fAMDLinFringeB0ZPos;};
    G4double GetAMDLinFringeTargZPos() {return fAMDLinFringeTargZPos;};
    G4ThreeVector GetAMDLocation() {return fAMDLocation;};

    G4VPhysicalVolume* GetWorld() {return fWorldPhysical;}
    G4VPhysicalVolume* GetXtal() {return fXtalPhysical;}
    G4VPhysicalVolume* GetAmorphous() {return fAmorphousPhysical;}
    G4VPhysicalVolume* GetAMDFrontGap() {return fAMDFrontGapPhysical;}
    G4VPhysicalVolume* GetAMD() {return fAMDPhysical;}

private:
    void ConstructXtalTarget();
    void ConstructAmorphousTarget();
    void ConstructAMDFrontGap();
    void ConstructAMD();

    void FillMagFieldHist();
    //void FillMagFieldGraph();

    G4bool fXtal;
    G4bool fAmorphous;
    G4bool fAMDFrontGap;
    G4bool fAMD;
    
    G4Material* fXtalMaterial;
    G4ThreeVector fXtalSize;
    G4ThreeVector fXtalLocation; 
    G4VSolid* fXtalSolid;
    G4LogicalVolume* fXtalLogic;
    G4VPhysicalVolume* fXtalPhysical;

    G4ThreeVector fAmorphousSize;
    G4ThreeVector fAmorphousLocation;
    G4ThreeVector fAmorphousDistance;
    G4ThreeVector fDipoleFieldValue;
    G4VSolid* fAmorphousSolid;
    G4LogicalVolume* fAmorphousLogic;
    G4VPhysicalVolume* fAmorphousPhysical;

    G4UniformMagField*    fDipoleField;
    G4FieldManager* fFieldMgr;

    G4VSolid* fMagSolid;
    G4LogicalVolume* fMagLogic;
    G4VPhysicalVolume* fMagPhysical;

    G4ThreeVector fAMDFrontGapSize;
    G4ThreeVector fAMDFrontGapLocation;
    G4VSolid* fAMDFrontGapSolid;
    G4LogicalVolume* fAMDFrontGapLogic;
    G4VPhysicalVolume* fAMDFrontGapPhysical;

    G4String fAMDOption;
    G4double fAMDR1min;
    G4double fAMDR1max;
    G4double fAMDR2min;
    G4double fAMDR2max;
    G4double fAMDLength;
    G4bool   fAMDLinFringeAct;
    G4double fAMDLinFringeB0;
    G4double fAMDLinFringeK;
    G4double fAMDLinFringeB0ZPos;
    G4double fAMDLinFringeTargZPos;
    G4ThreeVector fAMDLocation;
    G4VSolid* fAMDSolid;
    G4LogicalVolume* fAMDLogic;
    G4VPhysicalVolume* fAMDPhysical;
    G4VSolid* fAMDVacSolid;
    G4LogicalVolume* fAMDVacLogic;
    G4VPhysicalVolume* fAMDVacPhysical;
};

#endif

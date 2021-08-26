#ifndef InjectorPrimaryGeneratorAction_h
#define InjectorPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleGun.hh"

class G4ParticleGun;
class G4Event;
class InjectorDetectorConstruction;
class InjectorPrimaryGeneratorActionMessenger;

class InjectorPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    //InjectorPrimaryGeneratorAction();
    InjectorPrimaryGeneratorAction(InjectorDetectorConstruction* det);
    virtual ~InjectorPrimaryGeneratorAction();
    
    void GeneratePrimaries(G4Event*);

    //void SetXtalThick(G4double);
    void SetPrimaryEnergy(G4double);
    void SetSigmaE(G4double);
    void SetSigmaXY(G4double);
    void SetSigmaZ(G4double);
    void SetSigmaT(G4double);
    void SetSigmaPxy(G4double);

    static G4double GetVelocity(G4double, G4double);

    InjectorDetectorConstruction* fDet;

private:
    G4ParticleGun* fParticleGun;
    InjectorPrimaryGeneratorActionMessenger* fGunMessenger;

    G4double XtalThick;
    G4double primaryEnergy;
    G4double sigmaUserE;
    G4double sigmaUserXY;
    G4double sigmaUserZ;
    G4double sigmaUserT;
    G4double sigmaUserPxy;
};

#endif


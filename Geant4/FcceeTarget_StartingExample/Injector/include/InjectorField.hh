#ifndef InjectorField_H
#define InjectorField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

//#include <TFile.h>
//#include <TH2F.h>
//#include <TGraph2D.h>
//class TFile;
//class TH2F;
//class TGraph2D;

class InjectorDetectorConstruction;

class InjectorField : public G4MagneticField
{
  public:
    //InjectorField();
    InjectorField(InjectorDetectorConstruction* det);
    virtual ~InjectorField();

    virtual void GetFieldValue(const double Point[3], double *Bfield) const;

    //static void FindField(double r, double z, double &br, double &bz);

    InjectorDetectorConstruction* fDet;

  private:

    G4bool   fAMDLinFringeAct;
    G4double fAMDLinFringeB0;
    G4double fAMDLinFringeK;
    G4double fAMDLinFringeB0ZPos;
    G4double fAMDLinFringeTargZPos;

    //G4double fBz;
    //G4double fRmax_sq;
    //G4double fZmax;
    //TFile* fileFieldMap
    //void FindField(double &br, double &bz, double r, double z);
    //void FindField(double , double , double &, double &);
    //TH2F* histFieldMapR;
    //TH2F* histFieldMapZ;
    double b_dipole_x; // [T]
    double b_dipole_y;
    double b_dipole_z;
    double r_dipole_max;
    double z_dipole_min;
    double z_dipole_max;
    double r_mag_max;
    double z_mag_min;
    double z_mag_max;
    double l_amd_frontgap;
    double z_amd_entrance;
    double z_fringe_max;
    double z_fringe_l;
    double z_fringe_min;
};

#endif


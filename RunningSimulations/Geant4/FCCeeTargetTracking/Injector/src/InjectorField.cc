#include <TH2F.h>
//#include <TGraph2D.h>

#include "InjectorField.hh"
#include "G4SystemOfUnits.hh"
#include "InjectorDetectorConstruction.hh"
#include "G4ThreeVector.hh"

class TH2F;
//class TGraph2D;

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

InjectorField::InjectorField(InjectorDetectorConstruction* det) : fDet(det)
{

  fAMDLinFringeAct      = fDet->GetAMDLinFringeAct(); 		// bool
  fAMDLinFringeB0       = fDet->GetAMDLinFringeB0() / tesla; 	// [T]
  fAMDLinFringeK        = fDet->GetAMDLinFringeK(); 		// [T/mm]
  fAMDLinFringeB0ZPos   = fDet->GetAMDLinFringeB0ZPos() / mm; 	// [mm]
  fAMDLinFringeTargZPos = fDet->GetAMDLinFringeTargZPos() / mm; // [mm]

  G4ThreeVector v3_xtal_size = fDet->GetXtalSize();
  G4ThreeVector v3_xtal_loca = fDet->GetXtalLocation();
  G4ThreeVector v3_amor_size = fDet->GetAmorphousSize();
  G4ThreeVector v3_amor_loca = fDet->GetAmorphousLocation();
  G4ThreeVector v3_amor_dist = fDet->GetAmorphousDistance();
  G4ThreeVector v3_amd_frontgap_size = fDet->GetAMDFrontGapSize();
  G4ThreeVector v3_dipoleField = fDet->GetDipoleFieldValue();
  G4double rmax_amd = fDet->GetAMDR1max();
  G4double l_amd = fDet->GetAMDLength();

  double xtal_H = v3_xtal_size.x();
  double xtal_W = v3_xtal_size.z();
  double amor_H = v3_amor_size.x();
  double amor_W = v3_amor_size.z();
  double amor_D = v3_amor_dist.z();
  double amd_frontgap_H = v3_amd_frontgap_size.x();
  double amd_frontgap_W = v3_amd_frontgap_size.z();

  b_dipole_x = v3_dipoleField.x() / tesla; // [T]
  b_dipole_y = v3_dipoleField.y() / tesla;
  b_dipole_z = v3_dipoleField.z() / tesla;

  r_mag_max = xtal_H/2.0;
  if(r_mag_max < amor_H/2.0) r_mag_max = amor_H/2.0;
  if(r_mag_max < amd_frontgap_H/2.0) r_mag_max = amd_frontgap_H/2.0;
  z_mag_min = 0;
  z_mag_max = xtal_W + amor_D + amor_W + amd_frontgap_W;
  l_amd_frontgap = amd_frontgap_W;

  z_amd_entrance = z_mag_max;

  if ( l_amd > 0 ){
    if(r_mag_max < rmax_amd) r_mag_max = rmax_amd;
    z_mag_max += l_amd;
  }

  z_dipole_min = v3_xtal_loca.z() + v3_xtal_size.z()/2.0;
  z_dipole_max = z_dipole_min + v3_amor_dist.z();
  r_dipole_max = r_mag_max;

  z_fringe_max = v3_amor_loca.z() + v3_amor_size.z()/2.0;
  z_fringe_l = fAMDLinFringeB0 / fAMDLinFringeK - (fAMDLinFringeB0ZPos - fAMDLinFringeTargZPos);
  z_fringe_min = z_fringe_max - z_fringe_l;

}

InjectorField::~InjectorField()
{
}

void InjectorField::GetFieldValue(const double Point[3],double *Bfield) const
{

  double x = Point[0]; // [mm]
  double y = Point[1];
  double z = Point[2];
  double r = sqrt(pow(x,2)+pow(y,2));

  double bx = 0; // [T]
  double by = 0;
  double bz = 0;
  double br = 0;

  // Dipole field between target
  if ( z_dipole_min != z_dipole_max && z >= z_dipole_min && z <= z_dipole_max && r>=0 && r<= r_dipole_max) {
    bx += b_dipole_x; // [T]
    by += b_dipole_y;
    bz += b_dipole_z;
    //std::cout<<"DEBUG:: b_dipole_y = "<<b_dipole_y<<" [T]"<<std::endl;
  }

  // AMD linear fringe field
  if (fAMDLinFringeAct && z_fringe_l >= 0 && z >= z_fringe_min && z <= z_fringe_max){
    bx += -0.5 * fAMDLinFringeK * x;
    by += -0.5 * fAMDLinFringeK * y;
    bz += fAMDLinFringeK * (z-z_fringe_min);
  }

  /*
  bool use_linear_fringe = 0;
  bool use_real_fringe = 0;

  if ( z >= z_mag_min && z <= z_mag_max && r>=0 && r<= r_mag_max) {

    if (use_linear_fringe){
      double a = 0.34; // [T/mm]
      double B1 = 3.0; // [T]
      double z0 = z_mag_max - B1 / a;
      if(z >= z0 && z <= z_amd_entrance){
        bx = -0.5 * a * x;
        by = -0.5 * a * y;
        bz = a * (z - z_amd_entrance) + B1;
      }
      //else if(z > z_amd_entrance){
      //  TH2F* hFieldMapR = (TH2F*) fDet->GetHistFieldMapR();
      //  TH2F* hFieldMapZ = (TH2F*) fDet->GetHistFieldMapZ();
      //  double z_in_map = z - z_amd_entrance;
      //  int jbin = hFieldMapR->GetYaxis()->FindBin(r);
      //  int ibin = hFieldMapR->GetXaxis()->FindBin(z_in_map);
      //  br = hFieldMapR->GetBinContent(ibin,jbin);
      //  bz = hFieldMapZ->GetBinContent(ibin,jbin);
      //  bx = r==0?0:br*x/r;
      //  by = r==0?0:br*y/r;
      //}
    }

    if (use_real_fringe){

      double z_in_map = z - z_amd_entrance;

      TH2F* hFieldMapR = (TH2F*) fDet->GetHistFieldMapR();
      TH2F* hFieldMapZ = (TH2F*) fDet->GetHistFieldMapZ();
      int jbin = hFieldMapR->GetYaxis()->FindBin(r);
      int ibin = hFieldMapR->GetXaxis()->FindBin(z_in_map);
      br = hFieldMapR->GetBinContent(ibin,jbin);
      bz = hFieldMapZ->GetBinContent(ibin,jbin);

      //TGraph2D* gFieldMapR = (TGraph2D*) fDet->GetGraphFieldMapR();
      //TGraph2D* gFieldMapZ = (TGraph2D*) fDet->GetGraphFieldMapZ();
      //br = gFieldMapR->Interpolate(r,z_in_map);
      //bz = gFieldMapZ->Interpolate(r,z_in_map);

      bx = r==0?0:br*x/r;
      by = r==0?0:br*y/r;
    }

  }
  */

  Bfield[0] = bx*tesla;
  Bfield[1] = by*tesla;
  Bfield[2] = bz*tesla;
}


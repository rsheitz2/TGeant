#include "common.hxx"
#include "functions.h"

void boost_CS(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& spin_upStream, TLorentzVector& spin_downStream, 
	      TLorentzVector& virtual_photon, TLorentzVector& lepton1, TLorentzVector& lepton2){
  
  //3 componets of boost vector
  TVector3 photon_boost_z (virtual_photon.BoostVector() );
  
  //Boost along z-axis beam direction by longitudinal photon velocity
  TVector3 boost_z (0, 0, photon_boost_z.Z() );
  beam.Boost(-boost_z);
  target.Boost(-boost_z);
  spin_upStream.Boost(-boost_z);
  spin_downStream.Boost(-boost_z);
  virtual_photon.Boost(-boost_z);
  lepton1.Boost(-boost_z);
  lepton2.Boost(-boost_z);

  //Boost along x-axis transverse photon direction by transverse photon velocity in new frame
  TVector3 photon_boost_x (virtual_photon.BoostVector() );
  TVector3 boost_x (photon_boost_x.X(), 0, 0);
  beam.Boost(-boost_x);
  target.Boost(-boost_x);
  spin_upStream.Boost(-boost_x);
  spin_downStream.Boost(-boost_x);
  virtual_photon.Boost(-boost_x);
  lepton1.Boost(-boost_x);
  lepton2.Boost(-boost_x);

  //I wonder if boost can be along transverse then longitudinal direction...
}//boost_CS

/*void boost_CS(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& Spin, 
	      TLorentzVector& lepton1, TLorentzVector& lepton2){

  //3 componets of boost vector
  TLorentzVector virtual_photon = lepton1 + lepton2;
  TVector3 photon_boost_z (virtual_photon.BoostVector() );
  
  //Boost along z-axis beam direction by longitudinal photon velocity
  TVector3 boost_z (0, 0, photon_boost_z.Z() );
  beam.Boost(-boost_z);
  target.Boost(-boost_z);
  Spin.Boost(-boost_z);
  virtual_photon.Boost(-boost_z);
  lepton1.Boost(-boost_z);
  lepton2.Boost(-boost_z);

  //Boost along x-axis transverse photon direction by transverse photon velocity in new frame
  TVector3 photon_boost_x (virtual_photon.BoostVector() );
  TVector3 boost_x (photon_boost_x.X(), 0, 0);
  beam.Boost(-boost_x);
  target.Boost(-boost_x);
  Spin.Boost(-boost_x);
  lepton1.Boost(-boost_x);
  lepton2.Boost(-boost_x);

  //I wonder if boost can be along transverse then longitudinal direction...
  }//boost_CS//*/

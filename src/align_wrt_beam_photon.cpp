#include "common.hxx"
#include "functions.h"

void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& spin1, TLorentzVector& spin2, 
			   TLorentzVector& virtual_photon, TLorentzVector& lepton1, TLorentzVector& lepton2){
  //spin1 is upstream target or true target spin
  //spin2 is downstream target or simple (always polarized up) spin
  //lepton1 is mu minus
  //lepton2 is mu plus
  TVector3 beam3= beam.Vect();
  TVector3 target3= target.Vect();
  TVector3 spin1_3= spin1.Vect();
  TVector3 spin2_3= spin2.Vect();
  TVector3 vPhoton3 = virtual_photon.Vect();
  TVector3 lepton1_3= lepton1.Vect();
  TVector3 lepton2_3= lepton2.Vect();
  
  TVector3 khat = beam3.Unit();
  TVector3 xDir = vPhoton3 - vPhoton3.Dot(khat)*khat;
  TVector3 ihat = xDir.Unit();
  TVector3 jhat = khat.Cross(ihat);

  beam.SetX(beam3.Dot(ihat) );
  beam.SetY(beam3.Dot(jhat) );
  beam.SetZ(beam3.Dot(khat) );

  target.SetX(target3.Dot(ihat) );
  target.SetY(target3.Dot(jhat) );
  target.SetZ(target3.Dot(khat) );

  spin1.SetX(spin1_3.Dot(ihat) );
  spin1.SetY(spin1_3.Dot(jhat) );
  spin1.SetZ(spin1_3.Dot(khat) );

  spin2.SetX(spin2_3.Dot(ihat) );
  spin2.SetY(spin2_3.Dot(jhat) );
  spin2.SetZ(spin2_3.Dot(khat) );
  
  virtual_photon.SetX(vPhoton3.Dot(ihat) );
  virtual_photon.SetY(vPhoton3.Dot(jhat) );
  virtual_photon.SetZ(vPhoton3.Dot(khat) );

  lepton1.SetX(lepton1_3.Dot(ihat) );
  lepton1.SetY(lepton1_3.Dot(jhat) );
  lepton1.SetZ(lepton1_3.Dot(khat) );

  lepton2.SetX(lepton2_3.Dot(ihat) );
  lepton2.SetY(lepton2_3.Dot(jhat) );
  lepton2.SetZ(lepton2_3.Dot(khat) );

}//align_wrt_beam_photon//*/

void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& Spin,
			   TLorentzVector& lepton1, TLorentzVector& lepton2){

  TLorentzVector virtual_photon = lepton1 + lepton2;
  //Beam angles of rotation
  Double_t beam_azimuthal = beam.Phi();
  Double_t beam_polar = beam.Theta();
  
  //Rotate about azimuthal beam angle (only needed to ensure next rotation is about the Yaxis)
  //Another approach would be to find angle perpendicular to z-axis and beam momentum 
  //       and rotate by beam_polar about this angle
  beam.RotateZ(-beam_azimuthal);
  target.RotateZ(-beam_azimuthal);
  Spin.RotateZ(-beam_azimuthal);
  virtual_photon.RotateZ(-beam_azimuthal);
  lepton1.RotateZ(-beam_azimuthal);
  lepton2.RotateZ(-beam_azimuthal);

  //Rotate about polar angle
  beam.RotateY(-beam_polar);
  target.RotateY(-beam_polar);
  Spin.RotateY(-beam_polar);
  virtual_photon.RotateY(-beam_polar);
  lepton1.RotateY(-beam_polar);
  lepton2.RotateY(-beam_polar);

  //Virtual photon azimuthal angle when z-axis is aligned with the beam
  Double_t vPhoton_azimuthal = virtual_photon.Phi();

  //Rotate about azimuthal angle of virtual photon
  beam.RotateZ(-vPhoton_azimuthal);
  target.RotateZ(-vPhoton_azimuthal);
  Spin.RotateZ(-vPhoton_azimuthal);
  lepton1.RotateZ(-vPhoton_azimuthal);
  lepton2.RotateZ(-vPhoton_azimuthal);

}//align_wrt_beam_photon

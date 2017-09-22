#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__

inline Int_t DYtrigger(Int_t t){
  Int_t MT_LT = 0b010000000000000001;
  Int_t OT_LT = 0b010000000000000100;
  Int_t LT_LT = 0b010000000100000000;

  if (t == MT_LT || t == OT_LT || t == LT_LT) return 1;
  return 0;
}

inline Int_t DYtriggerMaster(Int_t t){
  Int_t master_MT_LT = 0b000000000000000001;
  Int_t master_OT_LT = 0b000000000000000100;
  Int_t master_LT_LT = 0b000000000100000000;

  if (t == master_MT_LT || t == master_OT_LT || t == master_LT_LT) return 1;
  return 0;
}

//Align all input vectors so z-axis is along the beam direction and x-axis is along transverse momentum of the virtual photon
/*void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& spin_upStream, TLorentzVector& spin_downStream, 
  TLorentzVector& virtual_photon, TLorentzVector& lepton1, TLorentzVector& lepton2);//*/

//Align all input vectors so z-axis is along the beam direction and x-axis is along transverse momentum of the virtual photon
void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& Spin,
			   TLorentzVector& lepton1, TLorentzVector& lepton2);

void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& spin_upStream, TLorentzVector& spin_downStream, 
			   TLorentzVector& virtual_photon, TLorentzVector& lepton1, TLorentzVector& lepton2);
  
//Boost from aligned with respect to beam to Collins Soper's center of momentum frame
//This is a boost along beam direction by longitudinal virtual photon momentum
//Followed by a boost along x-axis by the transverse momentum of the virtual photon
void boost_CS(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& spin_upStream, TLorentzVector& spin_downStream, 
	      TLorentzVector& virtual_photon, TLorentzVector& lepton1, TLorentzVector& lepton2);

//Boost from aligned with respect to beam to Collins Soper's center of momentum frame
//This is a boost along beam direction by longitudinal virtual photon momentum
//Followed by a boost along x-axis by the transverse momentum of the virtual photon
void boost_CS(TLorentzVector& beam, TLorentzVector& target, TLorentzVector& Spin, 
	      TLorentzVector& lepton1, TLorentzVector& lepton2);

inline Double_t correctQuadrant(const TLorentzVector& vec, Double_t angle)
{
  //returns 0 < angle < 2*pi rad
  //if (vec.X() <= 0 && vec.Y() > 0) angle = TMath::Pi() - angle; //2nd quadrant 
  //else if (vec.X() < 0 && vec.Y() <= 0) angle += TMath::Pi(); //3rd quadrant
  //else if (vec.X() >= 0 && vec.Y() < 0) angle = 2*TMath::Pi() - angle;//4th quadrant

  //returns -pi < angle < pi rad
  if (vec.X() <= 0 && vec.Y() > 0) angle = TMath::Pi() - angle; //2nd quadrant 
  else if (vec.X() < 0 && vec.Y() <= 0) angle -= TMath::Pi(); //3rd quadrant
  else if (vec.X() >= 0 && vec.Y() < 0) angle *= -1;//4th quadrant

  return angle;
};

inline Double_t azimuthalAngle(const TLorentzVector &vec)
{
  Double_t XYmag = sqrt(pow(vec.X(), 2) + pow(vec.Y(), 2) );
  Double_t angle = TMath::ACos(TMath::Abs(vec.X() )/XYmag);
  return correctQuadrant(vec, angle);
};

void BinData(TH1D** h1, Double_t valueCheck, Double_t valueFill, Double_t* bins, Int_t ih_first);

Bool_t BinAvg(Double_t *Avg, Int_t *count, Double_t binVal,
	      Double_t *binValBound, Double_t avgVal);

#endif

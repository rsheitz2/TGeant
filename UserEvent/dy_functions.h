#ifndef DY_FUNCTIONS_H
#define DY_FUNCTIONS_H

///////////////////////////////////////////////
//Common functions for Drell-Yan analysis
//Functions must be defined inline so their linkage is not external
///////////////////////////////////////////////

inline Bool_t Hit(const PaDetect& det, const PaTPar& p){
  PaTPar extrapolated;
  
  //Bool_t reached = p.Extrapolate(det.Z(), extrapolated, true); //True = Multiple scattering is used! (should be used)
  Bool_t reached = p.Extrapolate(det.Z(), extrapolated, false); //False = Multiple scattering NOT used!
  if (reached && det.InActive(extrapolated.X(), extrapolated.Y()) ) return true;

  return false;
}

inline Bool_t Hit(const PaDetect& det, const PaTrack& p){
  PaTPar extrapolated;
  
  Bool_t reached = p.Extrapolate(det.Z(), extrapolated);
  //Bool_t reached = p.Extrap(det.Z(), extrapolated); 
  if (reached && det.InActive(extrapolated.X(), extrapolated.Y()) ) return true;

  return false;
}

inline Bool_t Hit(const PaDetect& detfront, const PaDetect& detback1,
		  const PaDetect& detback2, const PaTrack& p){
  PaTPar extrap_front;
  Bool_t reached_front = p.Extrapolate(detfront.Z(), extrap_front);
  if (reached_front && detfront.InActive(extrap_front.X(),
					 extrap_front.Y()) ){
    PaTPar extrap_back1;
    Bool_t reached_back1 = extrap_front.Extrapolate(detback1.Z(), extrap_back1);
    if (reached_back1
	&& detback1.InActive(extrap_back1.X(), extrap_back1.Y()) ) return true;
    PaTPar extrap_back2;
    Bool_t reached_back2 = extrap_back1.Extrapolate(detback2.Z(), extrap_back2);
    if (reached_back2
	&& detback2.InActive(extrap_back2.X(), extrap_back2.Y()) ) return true;
  }
  
  return false;
}

inline Bool_t pointsToLAS(const PaTPar& p){
  const PaDetect& H1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG01Y1__") );
  const PaDetect& H2_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG02Y1__") );
  const PaDetect& H2_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG02Y2__") );

  return Hit(H1, p) && (Hit(H2_1, p) || Hit(H2_2, p) );
}

inline Bool_t pointsToLAS(const PaTrack& p){
  const PaDetect& H1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG01Y1__") );
  const PaDetect& H2_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG02Y1__") );
  const PaDetect& H2_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HG02Y2__") );

  //return Hit(H1, p) && (Hit(H2_1, p) || Hit(H2_2, p) );//first
  return (Hit(H1, p) && Hit(H2_1, p) ) || (Hit(H1, p) && Hit(H2_2, p) );//second
  //return Hit(H1, H2_2, H2_1, p);//third
}

inline Bool_t pointsToOuter(const PaTPar& p){
  const PaDetect& HO03Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO03Y1") );
  const PaDetect& HO04Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y1") );
  const PaDetect& HO04Y2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y2") );
  
  return Hit(HO03Y1, p) && (Hit(HO04Y1, p) || Hit(HO04Y2, p) );
}

inline Bool_t pointsToOuter(const PaTrack& p){
  //const PaDetect& HO03Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO03Y1") );//First
  //const PaDetect& HO04Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y1") );
  //const PaDetect& HO04Y2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y2") );

  const PaDetect& HO03Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO03Y1_m") );//second
  const PaDetect& HO04Y1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y1_m") );
  const PaDetect& HO04Y2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HO04Y2_m") );
  
  //return Hit(HO03Y1, p) && (Hit(HO04Y1, p) || Hit(HO04Y2, p) );//first
  return (Hit(HO03Y1, p) && Hit(HO04Y1, p) ) || (Hit(HO03Y1, p) && Hit(HO04Y2, p) );//second
  //return Hit(HO03Y1, HO04Y1, HO04Y2, p);//third
}

inline Bool_t pointsToMiddle(const PaTPar& p){
  //const PaDetect& HM04Y = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1") );
  //const PaDetect& HM05Y = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1") );

  const PaDetect& HM04Y_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1_u") );
  const PaDetect& HM04Y_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1_d") );
  const PaDetect& HM05Y_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1_u") );
  const PaDetect& HM05Y_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1_d") );

  //return Hit(HM04Y, p) && Hit(HM05Y, p);
  return (Hit(HM04Y_1, p) || Hit(HM04Y_2, p) ) && (Hit(HM05Y_1, p) || Hit(HM05Y_2, p) );
}

inline Bool_t pointsToMiddle(const PaTrack& p){
  //const PaDetect& HM04Y = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1") );
  //const PaDetect& HM05Y = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1") );

  const PaDetect& HM04Y_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1_u") );
  const PaDetect& HM04Y_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM04Y1_d") );
  const PaDetect& HM05Y_1 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1_u") );
  const PaDetect& HM05Y_2 = PaSetup::Ref().Detector(PaSetup::Ref().iDetFirst("HM05Y1_d") );
  
  //return Hit(HM04Y, p) && Hit(HM05Y, p);
  return (Hit(HM04Y_1, p) || Hit(HM04Y_2, p) ) && (Hit(HM05Y_1, p) || Hit(HM05Y_2, p) );
}

//inline void common2parVx(vector<Int_t>& vec, Int_t* vx_p1, Int_t size_p1, Int_t* vx_p2, Int_t size_p2){//common2parVx 
inline void common2parVx(vector<Int_t>& vec, Int_t* vx_p1, Int_t size_p1, Int_t* vx_p2, Int_t size_p2){//common2parVx
  for (Int_t i=0; i<size_p1; i++){//p1
    for (Int_t j=0; j<size_p2; j++){//p2
      //cout << "p1vx: " << vx_p1[i]  << "  vx_p2: " << vx_p2[j] << "  i: " << i  << "  j: " << j << endl;
      if (vx_p1[i] == vx_p2[j] ) vec.push_back(vx_p1[i] );
    }//p2
  }//p1/
}//common2parVx

inline void common2parVx(vector<Int_t>& vec, Int_t* vx_p1, Int_t size_p1, Int_t* vx_p2, Int_t size_p2, PaEvent& e){//common2parVx 
  for (Int_t i=0; i<size_p1; i++){//p1
    for (Int_t j=0; j<size_p2; j++){//p2
      //if (vx_p1[i] == vx_p2[j] ) vec.push_back(vx_p1[i] );
      if (vx_p1[i] == vx_p2[j] ) {
	vec.push_back(vx_p1[i] );
	const PaVertex& test_vx_1 = e.vVertex(vx_p1[i]);
	const PaVertex& test_vx_2 = e.vVertex(vx_p2[j]);
	if (test_vx_1.Z() != test_vx_2.Z() ) { cout << "Bad Vx!!!!!!!!!!!!!!!!!!!!!!!1" << endl;}
      }
    }//p2
  }//p1/
}//common2parVx

inline Int_t IsTrigValidation(Int_t trigMask, const PaTrack& tr_p1, const PaTrack& tr_p2){
  Int_t trigValidation = 0;
  if ( ((trigMask >> 8) & 1) && pointsToLAS(tr_p1) && pointsToLAS(tr_p2)) {trigValidation = 1;} //Last_Last
  if ( (trigMask >> 2) & 1){//Last_Outer
    if (pointsToOuter(tr_p1) ){
      if (pointsToLAS(tr_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(tr_p1) ){
      if (pointsToOuter(tr_p2) ) trigValidation = 1;
    }
  }
  //Last Middle trigger
  /*if (trigMask & 1){//Last_Middle
    if (pointsToMiddle(tr_p1) ){
    if (pointsToLAS(tr_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(tr_p1) ){
    if (pointsToMiddle(tr_p2) ) trigValidation = 1;
    }
    }//*/

  return trigValidation;
}

inline Int_t IsTrigValidation(Int_t trigMask, const PaTPar& traj_p1, const PaTPar& traj_p2){
  Int_t trigValidation = 0;
  if ( ((trigMask >> 8) & 1) && pointsToLAS(traj_p1) && pointsToLAS(traj_p2)) {trigValidation = 1;} //Last_Last
  if ( (trigMask >> 2) & 1){//Last_Outer
    if (pointsToOuter(traj_p1) ){
      if (pointsToLAS(traj_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(traj_p1) ){
      if (pointsToOuter(traj_p2) ) trigValidation = 1;
    }
  }
  //Last Middle trigger
  /*if (trigMask & 1){//Last_Middle
    if (pointsToMiddle(traj_p1) ){
    if (pointsToLAS(traj_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(traj_p1) ){
    if (pointsToMiddle(traj_p2) ) trigValidation = 1;
    }
    }//*/

  return trigValidation;
}

inline Int_t IsTrigValidation(Int_t trigMask, PaTPar& traj_p1, PaTPar& traj_p2){
  Int_t trigValidation = 0;
  if ( ((trigMask >> 8) & 1) && pointsToLAS(traj_p1) && pointsToLAS(traj_p2)) {trigValidation = 1;} //Last_Last
  if ( (trigMask >> 2) & 1){//Last_Outer
    if (pointsToOuter(traj_p1) ){
      if (pointsToLAS(traj_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(traj_p1) ){
      if (pointsToOuter(traj_p2) ) trigValidation = 1;
    }
  }
  //Last Middle trigger
  /*if (trigMask & 1){//Last_Middle
    if (pointsToMiddle(traj_p1) ){
    if (pointsToLAS(traj_p2)) trigValidation = 1;
    }
    else if (pointsToLAS(traj_p1) ){
    if (pointsToMiddle(traj_p2) ) trigValidation = 1;
    }
    }//*/

  return trigValidation;
}

inline Int_t IsImageCut(PaTPar traj_p1, PaTPar traj_p2, Int_t trigMask){
  traj_p1(5) = -traj_p1(5);
  traj_p2(5) = -traj_p2(5);

  return IsTrigValidation(trigMask, traj_p1, traj_p2);
}

inline void align_wrt_beam_photon(TLorentzVector& beam, TLorentzVector& target,
				  TLorentzVector& spin1, TLorentzVector& spin2, 
				  TLorentzVector& virtual_photon,
				  TLorentzVector& lepton1,
				  TLorentzVector& lepton2){
  //spin1 is upstream target or true target spin
  //spin2 is downstream target or simple (always polarized up) spin
  //lepton1 is mu minus or mu plus
  //lepton2 is mu plus or mu minus
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

#endif

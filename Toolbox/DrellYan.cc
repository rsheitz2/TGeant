#include "toolbox.hh"

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


void boost_CS(TLorentzVector& beam, TLorentzVector& target,
	      TLorentzVector& spin_upStream, TLorentzVector& spin_downStream, 
	      TLorentzVector& virtual_photon, TLorentzVector& lepton1,
	      TLorentzVector& lepton2){
  
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

  //Boost along x-axis transverse photon direction by transverse photon
  //velocity in new frame
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


/*! \class DrellYan
 *  \brief A template for your own implementation
 *  You can just create a copy of this file, replace 'ToolboxTemplate' with your own classname and start programming!
 */
class DrellYan : ToolboxPlugin
{
public:
  DrellYan(void);
  ~DrellYan(void) {}

  std::string getDescription(void);
  bool processEvent(T4Event* event);
  void endOfEvents(void);
  void beginOfEvents(void);

private:
  // place to declare some histograms
  TH1D *h_Cuts;

  //Before cuts
  TH1D *h_Vx[3];

  //Not understood
  TH1D *h_nBeam;
  TH1D *h_nTraj;

  //Event Tree
  TTree *Event;
  //Positively charged outgoing muon trajectory parameters at vertex
  Double_t phi_muP, theta_muP;
  Double_t qP_muP;
  Double_t muP_X, muP_Y, muP_Z, muP_E;
  //Negatively charged outgoing muon trajectory parameters at vertex
  Double_t phi_muM, theta_muM;
  Double_t qP_muM;
  Double_t muM_X, muM_Y, muM_Z, muM_E;
  //Vertex specific
  Double_t vx_z, vx_y, vx_x;
  //Virtual photon specific/Dimuon
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  Double_t vPhoton_M;
  Double_t vOpenAngle;
  //Beam pion trajectory parameters at vertex
  Double_t phi_pIn, theta_pIn;
  Double_t qP_pIn;
  Double_t pIn_X, pIn_Y, pIn_Z, pIn_E;
  //Drell-Yan Angles
  Double_t PhiS, PhiS_simple, Phi_CS, Theta_CS;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse;

  //Testing
  Int_t nPion, nProton, nNeutron, nMuP, nMuM;

};

static DrellYan* drellYann = new DrellYan();

DrellYan::DrellYan(void)
{
  myName = "DrellYan";
  pluginList::getInstance()->activated_classes.push_back(this);
}

std::string DrellYan::getDescription(void)
{
  std::string description = "Looking at generated z-vertex position and z-vertex position difference.";
  return description;
}

void DrellYan::beginOfEvents(void)
{
  // place to initialize the histograms
  h_Cuts = new TH1D("h_Cuts", "h_Cuts", 130, 0, 130);

  //Before cuts
  h_Vx[0] = new TH1D("h_Vx_z", "Generated Z vertex position", 1000, -500, 200);
  h_Vx[1] = new TH1D("h_Vx_x", "Generated X vertex position", 100, -5, 5);
  h_Vx[2] = new TH1D("h_Vx_y", "Generated Y vertex position", 100, -5, 5);

  //Not understood
  h_nBeam = new TH1D("h_nBeam", "Number of beam particles", 100, 0, 100);
  h_nTraj = new TH1D("h_nTraj", "Number of traj parameters", 100, 0, 100);

  //Event Tree
  Event = new TTree("Event", "Event");
  //Positively charged outgoing muon trajectory parameters at vertex
  Event->Branch("phi_muP", &phi_muP, "phi_muP/D");
  Event->Branch("theta_muP", &theta_muP, "theta_muP/D");
  Event->Branch("qP_muP", &qP_muP, "qP_muP/D");
  Event->Branch("muP_X", &muP_X, "muP_X/D");
  Event->Branch("muP_Y", &muP_Y, "muP_Y/D");
  Event->Branch("muP_Z", &muP_Z, "muP_Z/D");
  Event->Branch("muP_E", &muP_E, "muP_E/D");
  //Negatively charged outgoing muon trajectory parameters at vertex
  Event->Branch("phi_muM", &phi_muM, "phi_muM/D");
  Event->Branch("theta_muM", &theta_muM, "theta_muM/D");
  Event->Branch("qP_muM", &qP_muM, "qP_muM/D");
  Event->Branch("muM_X", &muM_X, "muM_X/D");
  Event->Branch("muM_Y", &muM_Y, "muM_Y/D");
  Event->Branch("muM_Z", &muM_Z, "muM_Z/D");
  Event->Branch("muM_E", &muM_E, "muM_E/D");
  //Vertex specific
  Event->Branch("vx_z", &vx_z, "vx_z/D");
  Event->Branch("vx_x", &vx_x, "vx_x/D");
  Event->Branch("vx_y", &vx_y, "vx_y/D");
  //Virtual photon specific/Dimuon
  Event->Branch("vPhoton_X", &vPhoton_X, "vPhoton_X/D");
  Event->Branch("vPhoton_Y", &vPhoton_Y, "vPhoton_Y/D");
  Event->Branch("vPhoton_Z", &vPhoton_Z, "vPhoton_Z/D");
  Event->Branch("vPhoton_E", &vPhoton_E, "vPhoton_E/D");
  Event->Branch("vPhoton_M", &vPhoton_M, "vPhoton_M/D");
  Event->Branch("vOpenAngle", &vOpenAngle, "vOpenAngle/D");
  //Beam pion trajectory parameters at vertex
  Event->Branch("phi_pIn", &phi_pIn, "phi_pIn/D");
  Event->Branch("theta_pIn", &theta_pIn, "theta_pIn/D");
  Event->Branch("qP_pIn", &qP_pIn, "qP_pIn/D");
  Event->Branch("pIn_X", &pIn_X, "pIn_X/D");
  Event->Branch("pIn_Y", &pIn_Y, "pIn_Y/D");
  Event->Branch("pIn_Z", &pIn_Z, "pIn_Z/D");
  Event->Branch("pIn_E", &pIn_E, "pIn_E/D");
  //DY-variables
  Event->Branch("x_beam", &x_beam, "x_beam/D");
  Event->Branch("x_target", &x_target, "x_target/D");
  Event->Branch("x_feynman", &x_feynman, "x_feynman/D");
  Event->Branch("q_transverse", &q_transverse, "q_transverse/D");

  //Testing
  Event->Branch("nPion", &nPion, "nPion/I");
  Event->Branch("nProton", &nProton, "nProton/I");
  Event->Branch("nNeutron", &nNeutron, "nNeutron/I");
  Event->Branch("nMuP", &nMuP, "nMuP/I");
  Event->Branch("nMuM", &nMuM, "nMuM/I");
  
}

bool DrellYan::processEvent(T4Event* event)
{
  // this function is called for each event
  // return true if this event should be saved, false if not
  // (only if -o is activated)
  Double_t M_proton = 0.938272;

  //cout << event->printEventInfo() << endl;
  //event->printEventInfo();
  //event->printBeamData();

  vx_x = event->beamData.vertexPosition[0]/10.0;
  vx_y = event->beamData.vertexPosition[1]/10.0;
  vx_z = event->beamData.vertexPosition[2]/10.0;

  //Before cuts
  h_Vx[0]->Fill(vx_z);
  h_Vx[1]->Fill(vx_x);
  h_Vx[2]->Fill(vx_y);

  //Not understood
  h_nBeam->Fill(event->beamData.nBeamParticle );
  h_nTraj->Fill(event->beamData.nTrajectories );

  vector<T4BeamParticle> vBeam = event->beamData.beamParticles;
  TLorentzVector muP;
  TLorentzVector muM;
  TLorentzVector pIn;
  TLorentzVector target;
  Bool_t hasMuP=false, hasMuM=false, hasMuBeam=false, hasMuTar=false;
  nPion = 0; nProton = 0; nNeutron = 0; nMuP = 0; nMuM = 0;
  for (vector<T4BeamParticle>::iterator it=vBeam.begin(); it!=vBeam.end(); it++){
    if (it->k[1] == 13 ) {
      muP.SetX(it->p[0]);
      muP.SetY(it->p[1]);
      muP.SetZ(it->p[2]);
      muP.SetE(it->p[3]);
      hasMuP = true;

      nMuP++;
    }
    else if (it->k[1] == -13) {
      muM.SetX(it->p[0]);
      muM.SetY(it->p[1]);
      muM.SetZ(it->p[2]);
      muM.SetE(it->p[3]);
      hasMuM = true;

      nMuM++;
    }
    else if (it->k[1] == 2212 && it->k[0] == -12) {//target
      target.SetX(it->p[0]);
      target.SetY(it->p[1]);
      target.SetZ(it->p[2]);
      target.SetE(it->p[3]);
      hasMuTar = true;

      nProton++;
    }
    else if (it->k[1] == 2112 && it->k[0] == -12 ){//neutron
      target.SetX(it->p[0]);
      target.SetY(it->p[1]);
      target.SetZ(it->p[2]);
      target.SetE(it->p[3]);
      hasMuTar = true;
      
      nNeutron++;
    }
    else if (it->k[1] == -211 && it->k[0] == -12) {//pion
      pIn.SetX(it->p[0]);
      pIn.SetY(it->p[1]);
      pIn.SetZ(it->p[2]);
      pIn.SetE(it->p[3]);
      hasMuBeam = true;

      nPion++;
    }
  }
  
  TLorentzVector diMu;
  if (hasMuP && hasMuM){
    TLorentzVector q = muP + muM;
    diMu.SetX(q.X() );
    diMu.SetY(q.Y() );
    diMu.SetZ(q.Z() );
    diMu.SetE(q.E() );
    
  }

  //Cuts
  Int_t cut_bin = 1, cut_space = 10;
  h_Cuts->Fill(cut_bin); cut_bin += cut_space; //All Data

  if ( vx_z < -320 || vx_z > -140) return false;//NH3 targets
  h_Cuts->Fill(cut_bin); cut_bin += cut_space;
    
  if(TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >= TMath::Power(3, 2)
     ) return false;//NH3 targets
  h_Cuts->Fill(cut_bin); cut_bin += cut_space;

  //Setup Vectors in different coordinate systems
  //Compass frame:
  TLorentzVector lv_Spin(0, 1.0, 0, 0);
  TLorentzVector lv_Spin_simple(0, 1.0, 0, 0);

  //Target frame
  TLorentzVector lv_beam_TF(pIn);
  TLorentzVector lv_target_TF(target);
  TLorentzVector lv_Spin_TF(lv_Spin);
  TLorentzVector lv_Spin_simple_TF(lv_Spin_simple);
  TLorentzVector lv_muPlus_TF(muP);
  TLorentzVector lv_muMinus_TF(muM);
  TLorentzVector lv_virtualPhoton_TF(diMu);
  align_wrt_beam_photon(lv_beam_TF, lv_target_TF, lv_Spin_TF,
			lv_Spin_simple_TF, lv_virtualPhoton_TF, lv_muPlus_TF,
			lv_muMinus_TF);

  //Boost from TF to CS
  TLorentzVector lv_beam_CS(lv_beam_TF);
  TLorentzVector lv_target_CS(lv_target_TF);
  TLorentzVector lv_Spin_CS(lv_Spin_TF);
  TLorentzVector lv_Spin_simple_CS(lv_Spin_simple_TF);
  TLorentzVector lv_muMinus_CS(lv_muMinus_TF);
  TLorentzVector lv_muPlus_CS(lv_muPlus_TF);
  TLorentzVector lv_virtualPhoton_CS(lv_virtualPhoton_TF);
  boost_CS(lv_beam_CS, lv_target_CS, lv_Spin_CS, lv_Spin_simple_CS,
	   lv_virtualPhoton_CS, lv_muPlus_CS, lv_muMinus_CS);

  
  //Event Tree  
  //Positively charged outgoing muon trajectory parameters at vertex
  phi_muP = muP.Phi();
  theta_muP = muP.Theta();
  qP_muP = muP.Vect().Mag();
  muP_X = muP.X();
  muP_Y = muP.Y();
  muP_Z = muP.Z();
  muP_E = muP.E();
  //Negatively charged outgoing muon trajectory parameters at vertex
  phi_muM = muM.Phi();
  theta_muM = muM.Theta();
  qP_muM = muM.Vect().Mag();
  muM_X = muM.X();
  muM_Y = muM.Y();
  muM_Z = muM.Z();
  muM_E = muM.E();
  //Vertex specific  //Defined previously
  
  //Virtual photon specific/Dimuon
  vPhoton_X = diMu.X();
  vPhoton_Y = diMu.Y();
  vPhoton_Z = diMu.Z();
  vPhoton_E = diMu.E();
  vPhoton_M = diMu.M();
  vOpenAngle = muP.Vect().Angle(muM.Vect() );
  //Beam pion trajectory parameters at vertex
  phi_pIn = pIn.Phi();
  theta_pIn = pIn.Theta();
  qP_pIn = pIn.Vect().Mag();
  pIn_X = pIn.X();
  pIn_Y = pIn.Y();
  pIn_Z = pIn.Z();
  pIn_E = pIn.E();
  //Drell-Yan Angles
  PhiS = lv_Spin_TF.Phi();
  PhiS_simple = lv_Spin_simple_TF.Phi();
  Phi_CS = lv_muMinus_CS.Phi();
  Theta_CS = lv_muMinus_CS.Theta();
  
  //DY-variables
  TLorentzVector lv_target(0, 0, 0, M_proton);
  x_beam = diMu.Mag2()/(2*diMu.Dot(pIn) );
  x_target = diMu.Mag2()/(2*diMu.Dot(lv_target) );
  x_feynman = x_beam - x_target;
  q_transverse = diMu.Vect().Cross(pIn.Vect() ).Mag()
    /(pIn.Vect().Mag() );

  Event->Fill();
  
 return true;
}

void DrellYan::endOfEvents(void)
{

  Int_t cut_bin = 1, cut_space = 10;
  h_Cuts->GetXaxis()->SetBinLabel(cut_bin, "AllData"); cut_bin+=cut_space;
  h_Cuts->GetXaxis()->SetBinLabel(cut_bin, "NH3 vx_z"); cut_bin+=cut_space;
  h_Cuts->GetXaxis()->SetBinLabel(cut_bin, "NH3 radius"); cut_bin+=cut_space;

  // place to write the histograms
  h_Cuts->Write();
  h_Vx[0]->Write();
  h_Vx[1]->Write();
  h_Vx[2]->Write();

  //Not understood
  h_nBeam->Write();
  h_nTraj->Write();

  //Event Tree  
  Event->Write();

}

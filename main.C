#include "common.hxx"
#include "functions.h"
#include <bitset>

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "" << endl;
    cout << "To be used with UserEvent410 for Monte Carlo Data" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "Enter data period i.e." << endl;
    cout << "./main W07" << endl;
    cout << "" << endl;
    cout << "Enter an additional argument to write output file i.e." << endl;
    cout << "./main W07 1" << endl;
    cout << "output file is named \"mainW##.root\"" << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  cout << "" << endl;
  TApplication theApp("tapp", &argc, argv); 
    
  TChain* T1 = new TChain("UserEvent410/Particles");
  TString period(argv[1]);
  if (period == "W07"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W07/\
Merged.root");
    cout << "Data from period:   W07" << endl;
  }
  else if (period == "W08"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W08/\
Merged.root");
    cout << "Data from period:   W08" << endl;
  }
  else if (period == "W09"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W09/\
Merged.root");
    cout << "Data from period:   W09" << endl;
  }
  else if (period == "W10"){
    //T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/DY_W10_5000.root");
    //T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/W10_bw_5000_riccardo.root");
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/W10_10x5000.root");
    /*T1->Add("/afs/cern.ch/user/r/rheitz/work/Analysis/TGeant/Output/	\
      dy2015_W10_2200/PhastOutput/DY_W10_2200.root");*/
    cout << "Data from period:   W10" << endl;
  }
  else if (period == "W11"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W11/\
Merged.root");
    cout << "Data from period:   W11" << endl;
  }
  else if (period == "W12"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W12/\
Merged.root");
    cout << "Data from period:   W12" << endl;
  }
  else if (period == "W13"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W13/\
Merged.root");
    cout << "Data from period:   W13" << endl;
  }
  else if (period == "W14"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W14/\
Merged.root");
    cout << "Data from period:   W14" << endl;
  }
  else if(period == "W15"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/W15/\
Merged.root");
    cout << "Data from period:   W15" << endl;
  }
  else if (period == "all"){
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/\
AllPeriods/Merged.root");
    cout << "Data from period:   All Periods Combined" << endl;
  }
  else {
    T1->Add("/afs/cern.ch/user/r/rheitz/work/phast.7.156/Robert/XCheck/\
AllPeriods/Merged.root");
    cout << "Data from period:   All Periods Combined" << endl;
  }
  cout << "" << endl;

  //Internal variables and binning
  Double_t M_proton = 0.938272;

  //Binnings
  Double_t xN_bounds[] = {0.00, 0.13, 0.19, 1.00};
  Double_t xPi_bounds[] = {0.00, 0.40, 0.56, 1.00};
  Double_t xF_bounds[] = {-1.0, 0.21, 0.41, 1.00};
  Double_t M_bounds[] = {4.30, 4.75, 5.50, 8.50};
  
  TVectorD tv_xN_bounds; tv_xN_bounds.Use(4, xN_bounds);
  TVectorD tv_xPi_bounds; tv_xPi_bounds.Use(4, xPi_bounds);
  TVectorD tv_xF_bounds; tv_xF_bounds.Use(4, xF_bounds);
  TVectorD tv_M_bounds; tv_M_bounds.Use(4, M_bounds);
  
  //////UserEvent Information
  /////////////////
  //Positively charged outgoing muon
  Int_t isBeam_p1, numVx_p1, numOutVx_p1, numVxpri_p1;
  //Positively charged outgoing muon trajectory parameters at vertex
  Double_t phi_traj1, theta_traj1;
  Double_t qP_traj1;
  Double_t vP1_X, vP1_Y, vP1_Z, vP1_E;
  //Negatively charged outgoing muon
  Int_t isBeam_p2, numVx_p2, numOutVx_p2, numVxpri_p2;
  //Negatively charged outgoing muon trajectory parameters at vertex
  Double_t phi_traj2, theta_traj2;
  Double_t qP_traj2;
  Double_t vP2_X, vP2_Y, vP2_Z, vP2_E;
  //Vertex specific
  Double_t vx_z, vx_y, vx_x; 
  Double_t vx_xVar, vx_yVar, vx_zVar;
  Double_t vx_Chi, vx_Chi_ndf;
  Int_t vx_ndf, vx_NOutParticles, vx_IsPrimary, vx_IsBestPrimary;
  //Positively charged outgoing muon track parameters
  Double_t Zfirst_tr1, Zlast_tr1;//Zmax_tr1, Zmin_tr1 same as first/last
  Int_t NHits_tr1;
  Double_t Chi2tot_tr1, Chi2tot_Ndf_tr1;
  Int_t Ndf_tr1;
  Double_t meanT_tr1, sigT_tr1;
  Double_t XX0_tr1;
  //Negatively charged outgoing muon track parameters
  Double_t Zfirst_tr2, Zlast_tr2;//Zmax_tr2, Zmin_tr2 same as first/last
  Int_t NHits_tr2;
  Double_t Chi2tot_tr2, Chi2tot_Ndf_tr2;
  Int_t Ndf_tr2;
  Double_t meanT_tr2, sigT_tr2;
  Double_t XX0_tr2;
  //Virtual photon specific/Dimuon
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  Double_t vDiMuon_invM;
  Double_t vOpenAngle;

  ////Beam particle specific
  //Beam pion particle specific
  Int_t isBeam_pIn, numVx_pIn, numOutVx_pIn;
  //Beam pion trajectory parameters at vertex
  Double_t phi_trajPIn, theta_trajPIn;
  Double_t qP_trajPIn;
  Double_t beam_X, beam_Y, beam_Z, beam_E;
  //Beam pion track specific
  Double_t Zlast_trPIn;//first/last are the same for beam particles
  //(first is closest det to target)
  Int_t NHits_trPIn;
  Double_t Chi2tot_trPIn, Chi2tot_Ndf_trPIn;
  Int_t Ndf_trPIn;
  Double_t meanT_trPIn, sigT_trPIn;
  Double_t XX0_trPIn;

  //Event specific
  Int_t NParticle, NTrack, NVertex; 
  Int_t trigMask;
  Long64_t event;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse;

  //Positively charged outgoing muon
  T1->SetBranchAddress("isBeam_p1", &isBeam_p1);
  T1->SetBranchAddress("numVx_p1", &numVx_p1);
  T1->SetBranchAddress("numOutVx_p1", &numOutVx_p1);
  T1->SetBranchAddress("numVxpri_p1", &numVxpri_p1);
  //Positively charnged outgoing muon trajectory parameters at vertex
  T1->SetBranchAddress("phi_traj1", &phi_traj1);
  T1->SetBranchAddress("theta_traj1", &theta_traj1);
  T1->SetBranchAddress("qP_traj1", &qP_traj1);
  T1->SetBranchAddress("vP1_X", &vP1_X);
  T1->SetBranchAddress("vP1_Y", &vP1_Y);
  T1->SetBranchAddress("vP1_Z", &vP1_Z);
  T1->SetBranchAddress("vP1_E", &vP1_E);
  //Negatively charged outgoing muon
  T1->SetBranchAddress("isBeam_p2", &isBeam_p2);
  T1->SetBranchAddress("numVx_p2", &numVx_p2);
  T1->SetBranchAddress("numOutVx_p2", &numOutVx_p2);
  T1->SetBranchAddress("numVxpri_p2", &numVxpri_p2);
  //Negatively charnged outgoing muon trajectory parameters at vertex
  T1->SetBranchAddress("phi_traj2", &phi_traj2);
  T1->SetBranchAddress("theta_traj2", &theta_traj2);
  T1->SetBranchAddress("qP_traj2", &qP_traj2);
  T1->SetBranchAddress("vP2_X", &vP2_X);
  T1->SetBranchAddress("vP2_Y", &vP2_Y);
  T1->SetBranchAddress("vP2_Z", &vP2_Z);
  T1->SetBranchAddress("vP2_E", &vP2_E);
  //Vertex specific
  T1->SetBranchAddress("vx_z", &vx_z);
  T1->SetBranchAddress("vx_x", &vx_x);
  T1->SetBranchAddress("vx_y", &vx_y);
  T1->SetBranchAddress("vx_zVar", &vx_zVar);
  T1->SetBranchAddress("vx_xVar", &vx_xVar);
  T1->SetBranchAddress("vx_yVar", &vx_yVar);
  T1->SetBranchAddress("vx_Chi", &vx_Chi);
  T1->SetBranchAddress("vx_Chi_ndf", &vx_Chi_ndf);
  T1->SetBranchAddress("vx_ndf", &vx_ndf);
  T1->SetBranchAddress("vx_NOutParticles", &vx_NOutParticles);
  T1->SetBranchAddress("vx_IsPrimary", &vx_IsPrimary);
  T1->SetBranchAddress("vx_IsBestPrimary", &vx_IsBestPrimary);
  //Positively charnged outgoing muon track parameters
  T1->SetBranchAddress("Zfirst_tr1", &Zfirst_tr1);
  T1->SetBranchAddress("Zlast_tr1", &Zlast_tr1);
  T1->SetBranchAddress("NHits_tr1", &NHits_tr1);
  T1->SetBranchAddress("Chi2tot_tr1", &Chi2tot_tr1);
  T1->SetBranchAddress("Chi2tot_Ndf_tr1", &Chi2tot_Ndf_tr1);
  T1->SetBranchAddress("Ndf_tr1", &Ndf_tr1);
  T1->SetBranchAddress("meanT_tr1", &meanT_tr1);
  T1->SetBranchAddress("sigT_tr1", &sigT_tr1);
  T1->SetBranchAddress("XX0_tr1", &XX0_tr1);
  //Negatively charnged outgoing muon track parameters
  T1->SetBranchAddress("Zfirst_tr2", &Zfirst_tr2);
  T1->SetBranchAddress("Zlast_tr2", &Zlast_tr2);
  T1->SetBranchAddress("NHits_tr2", &NHits_tr2);
  T1->SetBranchAddress("Chi2tot_tr2", &Chi2tot_tr2);
  T1->SetBranchAddress("Chi2tot_Ndf_tr2", &Chi2tot_Ndf_tr2);
  T1->SetBranchAddress("Ndf_tr2", &Ndf_tr2);
  T1->SetBranchAddress("meanT_tr2", &meanT_tr2);
  T1->SetBranchAddress("sigT_tr2", &sigT_tr2);
  T1->SetBranchAddress("XX0_tr2", &XX0_tr2);
  //Virtual photon specific/Dimuon
  T1->SetBranchAddress("vPhoton_X", &vPhoton_X);
  T1->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  T1->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  T1->SetBranchAddress("vPhoton_E", &vPhoton_E);
  T1->SetBranchAddress("vDiMuon_invM", &vDiMuon_invM);
  T1->SetBranchAddress("vOpenAngle", &vOpenAngle);

  ////Beam particle specific
  //Beam pion particle specific
  T1->SetBranchAddress("isBeam_pIn", &isBeam_pIn);
  T1->SetBranchAddress("numVx_pIn", &numVx_pIn);
  T1->SetBranchAddress("numOutVx_pIn", &numOutVx_pIn);
  //Beam pion trajectory parameters at vertex
  T1->SetBranchAddress("phi_trajPIn", &phi_trajPIn);
  T1->SetBranchAddress("theta_trajPIn", &theta_trajPIn);
  T1->SetBranchAddress("qP_trajPIn", &qP_trajPIn);
  T1->SetBranchAddress("beam_X", &beam_X);
  T1->SetBranchAddress("beam_Y", &beam_Y);
  T1->SetBranchAddress("beam_Z", &beam_Z);
  T1->SetBranchAddress("beam_E", &beam_E);
  //Beam pion track specific
  T1->SetBranchAddress("Zlast_trPIn", &Zlast_trPIn);
  T1->SetBranchAddress("NHits_trPIn", &NHits_trPIn);
  T1->SetBranchAddress("Chi2tot_trPIn", &Chi2tot_trPIn);
  T1->SetBranchAddress("Chi2tot_Ndf_trPIn", &Chi2tot_Ndf_trPIn);
  T1->SetBranchAddress("Ndf_trPIn", &Ndf_trPIn);
  T1->SetBranchAddress("meanT_trPIn", &meanT_trPIn);
  T1->SetBranchAddress("sigT_trPIn", &sigT_trPIn);
  T1->SetBranchAddress("XX0_trPIn", &XX0_trPIn);
  
  //Event specific
  T1->SetBranchAddress("NParticle", &NParticle);
  T1->SetBranchAddress("NTrack", &NTrack);
  T1->SetBranchAddress("NVertex", &NVertex);
  T1->SetBranchAddress("trigMask", &trigMask);
  T1->SetBranchAddress("event", &event);
  //DY-variables
  T1->SetBranchAddress("x_beam", &x_beam);
  T1->SetBranchAddress("x_target", &x_target);
  T1->SetBranchAddress("x_feynman", &x_feynman);
  T1->SetBranchAddress("q_transverse", &q_transverse);

  
  TH1D* hCuts = new TH1D("hCuts", "hCuts", 200, 0, 200);
  Int_t cut_bin = 1, cut_space = 10;

  TTree *tree = new TTree("pT_Weighted", "pT_Weighted");
  Double_t PhiS, PhiS_simple, Phi_CS, Theta_CS;
  Int_t targetPosition;
  Double_t xTarg;
  Double_t Spin;
  tree->Branch("PhiS", &PhiS, "PhiS/D");
  tree->Branch("PhiS_simple", &PhiS_simple, "PhiS_simple/D");
  tree->Branch("Phi_CS", &Phi_CS, "Phi_CS/D");
  tree->Branch("Theta_CS", &Theta_CS, "Theta_CS/D");
  tree->Branch("x_beam", &x_beam, "x_beam/D");
  tree->Branch("x_target", &x_target, "x_target/D");
  tree->Branch("xTarg", &xTarg, "xTarg/D");
  tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
  tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
  tree->Branch("Mmumu", &vDiMuon_invM, "Mmumu/D");
  tree->Branch("targetPosition", &targetPosition, "targetPosition/I");
  tree->Branch("theta_traj1", &theta_traj1, "theta_traj1/D");
  tree->Branch("qP_traj1", &qP_traj1, "qP_traj1/D");
  tree->Branch("theta_traj2", &theta_traj2, "theta_traj2/D");
  tree->Branch("qP_traj2", &qP_traj2, "qP_traj2/D");
  tree->Branch("vx_z", &vx_z, "vx_z/D");
  tree->Branch("vx_x", &vx_x, "vx_x/D");
  tree->Branch("vx_y", &vx_y, "vx_y/D");

  
  Int_t tree_entries = T1->GetEntries();
  //Int_t tree_entries = 1000;//Debug
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  for (Int_t i=0; i<tree_entries; i++){
    T1->GetEntry(i, 0);
    
    //Cuts
    cut_bin = 1;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//All Data

    TLorentzVector lv_p1_Mu(vP1_X, vP1_Y, vP1_Z, vP1_E);
    TLorentzVector lv_p2_Mu(vP2_X, vP2_Y, vP2_Z, vP2_E);
    TLorentzVector lv_diMu = lv_p1_Mu + lv_p2_Mu;
    TLorentzVector lv_target_1 (0, 0, 0, M_proton);
    xTarg = lv_diMu.Mag2()/(2*lv_diMu.Dot(lv_target_1) );

    if (x_beam < 0.0 || x_beam > 1.0) continue;
    if (x_target < 0.0 || x_target > 1.0) continue;
    if (x_feynman < -1.0 || x_feynman > 1.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Physical Kinematics

    if (q_transverse < 0.4 || q_transverse > 5.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//qT cuts
    
    /*if ( (vx_z < -294.5 || vx_z > -239.3) && (vx_z < -219.5 || vx_z > -164.3)
	 ) continue;//NH3 targets
	 hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target z-cut
    
    if(TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >= TMath::Power(1.9, 2)
       ) continue;//NH3 targets
       hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target radial cut*/
        
    ////All data after cuts
    //////////////
    ///////////////General useful quantities
    Double_t beam[] = {beam_X, beam_Y, beam_Z, beam_E};
    Double_t target[] = {M_proton};
    Double_t muPlus[] = {vP1_X, vP1_Y, vP1_Z, vP1_E};
    Double_t muMinus[] = {vP2_X, vP2_Y, vP2_Z, vP2_E};
    Double_t virtualPhoton[] = {vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E};
    TLorentzVector lv_photon_main(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);
    TLorentzVector lv_beam_main(beam_X, beam_Y, beam_Z, beam_E);
    //Int_t period = -1;//period 1 defined as upstream up, downstream down

    //if (vx_z >= -294.5 && vx_z <= -239.3){//Up stream NH3
    if (vx_z <= -230){//No target cuts
      Spin = 1.0;
      targetPosition = 0;

    }//Up stream
    //else if (vx_z >= -219.5 && vx_z <= -164.3){//Down stream NH3
    else if (vx_z >= -230){//No target cuts
      Spin = -1.0;
      targetPosition = 1;

    }
    else{
      Spin = -10.0;
      targetPosition = -10;
      cout << "Target position is wrong after cuts" << endl;
      exit(EXIT_FAILURE);
    }

    //Setup Vectors in different coordinate systems
    //Compass frame:
    TLorentzVector lv_beam(beam[0], beam[1], beam[2], beam[3]);
    TLorentzVector lv_target(0, 0, 0, target[0]);
    TLorentzVector lv_Spin(0, Spin, 0, 0);
    TLorentzVector lv_Spin_simple(0, 1, 0, 0);
    TLorentzVector lv_muPlus(muPlus[0], muPlus[1], muPlus[2], muPlus[3]);
    TLorentzVector lv_muMinus(muMinus[0], muMinus[1], muMinus[2], muMinus[3]);
    TLorentzVector lv_virtualPhoton(virtualPhoton[0], virtualPhoton[1],
				    virtualPhoton[2], virtualPhoton[3]);

    //Target frame:
    TLorentzVector lv_beam_TF(lv_beam);
    TLorentzVector lv_target_TF(lv_target);
    TLorentzVector lv_Spin_TF(lv_Spin);
    TLorentzVector lv_Spin_simple_TF(lv_Spin_simple);
    TLorentzVector lv_muPlus_TF(lv_muPlus);
    TLorentzVector lv_muMinus_TF(lv_muMinus);
    TLorentzVector lv_virtualPhoton_TF(lv_virtualPhoton);
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

    Double_t PhiS_lab = lv_Spin.Phi() - lv_virtualPhoton.Phi();
    if(PhiS_lab > TMath::Pi()) PhiS_lab = -2*TMath::Pi() + PhiS_lab;
    else if (PhiS_lab < -1.0*TMath::Pi() ) PhiS_lab = 2*TMath::Pi() + PhiS_lab;
    PhiS = lv_Spin_TF.Phi();
    PhiS_simple = lv_Spin_simple_TF.Phi();
    Phi_CS = lv_muMinus_CS.Phi();
    Theta_CS = lv_muMinus_CS.Theta();

    tree->Fill();
  }//tree entries    

  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  //Cuts histogram
  const Int_t nCuts = 5;
  TString cutNames[nCuts] = {"All Data", "xPion,xN,xF", "0.4<qT<5",
			     "Target z-cut", "Target radius"};
  for (Int_t i=0, j=1; i<nCuts; i++, j+=cut_space){
    Int_t bin_index = hCuts->GetXaxis()->FindBin(j);
    hCuts->GetXaxis()->SetBinLabel(bin_index, cutNames[i]);
  }
  //hCuts->Draw();

  TString outFile="main";
  if (argc < 3) cout << "No file output" << endl;
  else {
    outFile += argv[1]; outFile += ".root";
    TFile *myFile = new TFile(outFile, "RECREATE");
    hCuts->Write();
    tree->Write();
  
    tv_xN_bounds.Write("tv_xN_bounds");
    tv_xPi_bounds.Write("tv_xPi_bounds");
    tv_xF_bounds.Write("tv_xF_bounds");
    tv_M_bounds.Write("tv_M_bounds");

    cout << myFile->GetName() << " was written" << endl;
    myFile->Close();
  }
    
  theApp.Run();//Needed to make root graphics work on C++
}//main

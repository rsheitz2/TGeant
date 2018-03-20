#ifndef DY_VARIABLES_H
#define DY_VARIABLES_H

//Constant ([GeV])
const static Double_t M_mu = G3partMass[5]; //mass muon+ & muon-
const static Double_t M_pi = G3partMass[8]; //mass pion+ & pion-
const static Double_t M_proton = G3partMass[14]; //mass proton
//const static Int_t maxCuts = 17;//Should be the same as nCuts in JobEnd.cc
//const static Int_t nCutHist = 21;
const static Int_t maxCuts = 15;//Should be the same as nCuts in JobEnd.cc
const static Int_t nCutHist = 12;
const static Int_t n2D_cutHist = 3;

//Positively charged outgoing muon
static Int_t isBeam_p1, numVx_p1, numOutVx_p1, numVxpri_p1;
//Positively charged outgoing muon trajectory parameters at vertex
static Double_t phi_traj1, theta_traj1;
static Double_t qP_traj1;
static Double_t vP1_X, vP1_Y, vP1_Z, vP1_E;
//Negatively charged outgoing muon
static Int_t isBeam_p2, numVx_p2, numOutVx_p2, numVxpri_p2;
//Negatively charged outgoing muon trajectory parameters at vertex
static Double_t phi_traj2, theta_traj2;
static Double_t qP_traj2;
static Double_t vP2_X, vP2_Y, vP2_Z, vP2_E;
//Vertex specific
static Double_t vx_z, vx_y, vx_x; 
static Double_t vx_xVar, vx_xyVar, vx_yVar, vx_xzVar, vx_yzVar, vx_zVar;
static Double_t vx_Chi, vx_Chi_ndf;
static Int_t vx_ndf, vx_NOutParticles, vx_IsPrimary, vx_IsBestPrimary;
//Positively charged outgoing muon track parameters
static Double_t Zfirst_tr1, Zlast_tr1;//Zmax_tr1, Zmin_tr1 same as first/last
static Int_t NHits_tr1;
static Double_t Chi2tot_tr1, Chi2tot_Ndf_tr1;
static Int_t Ndf_tr1;
static Double_t meanT_tr1, sigT_tr1;
static Double_t XX0_tr1;
//Negatively charged outgoing muon track parameters
static Double_t Zfirst_tr2, Zlast_tr2;//Zmax_tr2, Zmin_tr2 same as first/last
static Int_t NHits_tr2;
static Double_t Chi2tot_tr2, Chi2tot_Ndf_tr2;
static Int_t Ndf_tr2;
static Double_t meanT_tr2, sigT_tr2;
static Double_t XX0_tr2;
//Virtual photon specific/Dimuon
static Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
static Double_t vDiMuon_invM;
static Double_t vOpenAngle;

////Beam particle specific
//Beam pion particle specific
static Int_t isBeam_pIn, numVx_pIn, numOutVx_pIn;
//Beam pion trajectory parameters at vertex
static Double_t phi_trajPIn, theta_trajPIn;
static Double_t qP_trajPIn;
static Double_t beam_X, beam_Y, beam_Z, beam_E;
//Beam pion track specific
static Double_t Zlast_trPIn;//first/last are the same for beam particles
//(first is closest det to target)
static Int_t NHits_trPIn;
static Double_t Chi2tot_trPIn, Chi2tot_Ndf_trPIn;
static Int_t Ndf_trPIn;
static Double_t meanT_trPIn, sigT_trPIn;
static Double_t XX0_trPIn;

//Event specific
static Int_t NParticle, NTrack, NVertex; 
static Int_t trigMask, MasterTrigMask;
static Long64_t event, RunNum, SpillNum;
//DY-variables
static Double_t x_beam, x_target, x_feynman, q_transverse;
//Target Polarization
static Double_t avgUpStream, avgDownStream;
static Double_t N14_UpStream, N14_DownStream;
static Double_t upStreamCoil1, upStreamCoil2, upStreamCoil3, upStreamCoil4;
static Double_t upStreamCoil5;
static Double_t downStreamCoil6, downStreamCoil7, downStreamCoil8;
static Double_t downStreamCoil9, downStreamCoil10;
static Double_t Polarization, dilutionFactor, error_dilutionFactor;

//Monte Carlo specific
//MC Positively charged outgoing muon track parameters
static Int_t pid_MCtr1, NMCHits_tr1;
static Double_t Pinv_MCtr1, theta_MCtr1, phi_MCtr1;
static Double_t vMCtr1_X, vMCtr1_Y, vMCtr1_Z, vMCtr1_E;
//MC Negatively charged outgoing muon track parameters
static Int_t pid_MCtr2, NMCHits_tr2;
static Double_t Pinv_MCtr2, theta_MCtr2, phi_MCtr2;
static Double_t vMCtr2_X, vMCtr2_Y, vMCtr2_Z, vMCtr2_E;
//MC Beam muon track parameters
static Int_t pid_MCtrIn, NMCHits_trIn, IsBeam_MCtrIn;
static Double_t Pinv_MCtrIn, theta_MCtrIn, phi_MCtrIn;
static Double_t vMCtrIn_X, vMCtrIn_Y, vMCtrIn_Z, vMCtrIn_E;
//MC DY-variables
static Double_t MC_x_beam, MC_x_target, MC_x_feynman, MC_q_transverse;

inline void SetupMCTreeVariables(TTree *tree){
  //Positively charged outgoing muon
  tree->Branch("isBeam_p1", &isBeam_p1, "isBeam_p1/I");
  tree->Branch("numVx_p1", &numVx_p1, "numVx_p1/I");
  tree->Branch("numOutVx_p1", &numOutVx_p1, "numOutVx_p1/I");
  tree->Branch("numVxpri_p1", &numVxpri_p1, "numVxpri_p1/I");
  tree->Branch("phi_traj1", &phi_traj1, "phi_traj1/D");
  tree->Branch("theta_traj1", &theta_traj1, "theta_traj1/D");
  tree->Branch("qP_traj1", &qP_traj1, "qP_traj1/D");
  tree->Branch("vP1_X", &vP1_X, "vP1_X/D");
  tree->Branch("vP1_Y", &vP1_Y, "vP1_Y/D");
  tree->Branch("vP1_Z", &vP1_Z, "vP1_Z/D");
  tree->Branch("vP1_E", &vP1_E, "vP1_E/D");
  //Negatively charged outgoing muon
  tree->Branch("isBeam_p2", &isBeam_p2, "isBeam_p2/I");
  tree->Branch("numVx_p2", &numVx_p2, "numVx_p2/I");
  tree->Branch("numOutVx_p2", &numOutVx_p2, "numOutVx_p2/I");
  tree->Branch("numVxpri_p2", &numVxpri_p2, "numVxpri_p2/I");
  //Negatively charged outgoing muon trajectory parameters at vertex
  tree->Branch("phi_traj2", &phi_traj2, "phi_traj2/D");
  tree->Branch("theta_traj2", &theta_traj2, "theta_traj2/D");
  tree->Branch("qP_traj2", &qP_traj2, "qP_traj2/D");
  tree->Branch("vP2_X", &vP2_X, "vP2_X/D");
  tree->Branch("vP2_Y", &vP2_Y, "vP2_Y/D");
  tree->Branch("vP2_Z", &vP2_Z, "vP2_Z/D");
  tree->Branch("vP2_E", &vP2_E, "vP2_E/D");
  //Vertex specific
  tree->Branch("vx_z", &vx_z, "vx_z/D");
  tree->Branch("vx_x", &vx_x, "vx_x/D");
  tree->Branch("vx_y", &vx_y, "vx_y/D");
  tree->Branch("vx_zVar", &vx_zVar, "vx_zVar/D");
  tree->Branch("vx_xVar", &vx_xVar, "vx_xVar/D");
  tree->Branch("vx_xyVar", &vx_xyVar, "vx_xyVar/D");
  tree->Branch("vx_yVar", &vx_yVar, "vx_yVar/D");
  tree->Branch("vx_xzVar", &vx_xzVar, "vx_xzVar/D");
  tree->Branch("vx_yzVar", &vx_yzVar, "vx_yzVar/D");
  tree->Branch("vx_Chi", &vx_Chi, "vx_Chi/D");
  tree->Branch("vx_Chi_ndf", &vx_Chi_ndf, "vx_Chi_ndf/D");
  tree->Branch("vx_ndf", &vx_ndf, "vx_ndf/I");
  tree->Branch("vx_NOutParticles", &vx_NOutParticles, "vx_NOutParticles/I");
  tree->Branch("vx_IsPrimary", &vx_IsPrimary, "vx_IsPrimary/I");
  tree->Branch("vx_IsBestPrimary", &vx_IsBestPrimary, "vx_IsBestPrimary/I");
  //Positively charged outgoing muon track parameters
  tree->Branch("Zfirst_tr1", &Zfirst_tr1, "Zfirst_tr1/D");
  tree->Branch("Zlast_tr1", &Zlast_tr1, "Zlast_tr1/D");
  tree->Branch("NHits_tr1", &NHits_tr1, "NHits_tr1/I");
  tree->Branch("Chi2tot_tr1", &Chi2tot_tr1, "Chi2tot_tr1/D");
  tree->Branch("Chi2tot_Ndf_tr1", &Chi2tot_Ndf_tr1, "Chi2tot_Ndf_tr1/D");
  tree->Branch("Ndf_tr1", &Ndf_tr1, "Ndf_tr1/I");
  tree->Branch("meanT_tr1", &meanT_tr1, "meanT_tr1/D");
  tree->Branch("sigT_tr1", &sigT_tr1, "sigT_tr1/D");
  tree->Branch("XX0_tr1", &XX0_tr1, "XX0_tr1/D");
  //Negatively charged outgoing muon track parameters
  tree->Branch("Zfirst_tr2", &Zfirst_tr2, "Zfirst_tr2/D");
  tree->Branch("Zlast_tr2", &Zlast_tr2, "Zlast_tr2/D");
  tree->Branch("NHits_tr2", &NHits_tr2, "NHits_tr2/I");
  tree->Branch("Chi2tot_tr2", &Chi2tot_tr2, "Chi2tot_tr2/D");
  tree->Branch("Chi2tot_Ndf_tr2", &Chi2tot_Ndf_tr2, "Chi2tot_Ndf_tr2/D");
  tree->Branch("Ndf_tr2", &Ndf_tr2, "Ndf_tr2/I");
  tree->Branch("meanT_tr2", &meanT_tr2, "meanT_tr2/D");
  tree->Branch("sigT_tr2", &sigT_tr2, "sigT_tr2/D");
  tree->Branch("XX0_tr2", &XX0_tr2, "XX0_tr2/D");
  //Virtual photon specific/Dimuon
  tree->Branch("vPhoton_X", &vPhoton_X, "vPhoton_X/D");
  tree->Branch("vPhoton_Y", &vPhoton_Y, "vPhoton_Y/D");
  tree->Branch("vPhoton_Z", &vPhoton_Z, "vPhoton_Z/D");
  tree->Branch("vPhoton_E", &vPhoton_E, "vPhoton_E/D");
  tree->Branch("vDiMuon_invM", &vDiMuon_invM, "vDiMuon_invM/D");
  tree->Branch("vOpenAngle", &vOpenAngle, "vOpenAngle/D");

  ////Beam particle specific
  //Beam pion particle specific
  tree->Branch("isBeam_pIn", &isBeam_pIn, "isBeam_pIn/I");
  tree->Branch("numVx_pIn", &numVx_pIn, "numVx_pIn/I");
  tree->Branch("numOutVx_pIn", &numOutVx_pIn, "numOutVx_pIn/I");
  //Beam pion trajectory parameters at vertex
  tree->Branch("phi_trajPIn", &phi_trajPIn, "phi_trajPIn/D");
  tree->Branch("theta_trajPIn", &theta_trajPIn, "theta_trajPIn/D");
  tree->Branch("qP_trajPIn", &qP_trajPIn, "qP_trajPIn/D");
  tree->Branch("beam_X", &beam_X, "beam_X/D");
  tree->Branch("beam_Y", &beam_Y, "beam_Y/D");
  tree->Branch("beam_Z", &beam_Z, "beam_Z/D");
  tree->Branch("beam_E", &beam_E, "beam_E/D");
  //Beam pion track specific
  tree->Branch("Zlast_trPIn", &Zlast_trPIn, "Zlast_trPIn/D");
  tree->Branch("NHits_trPIn", &NHits_trPIn, "NHits_trPIn/I");
  tree->Branch("Chi2tot_trPIn", &Chi2tot_trPIn, "Chi2tot_trPIn/D");
  tree->Branch("Chi2tot_Ndf_trPIn", &Chi2tot_Ndf_trPIn, "Chi2tot_Ndf_trPIn/D");
  tree->Branch("Ndf_trPIn", &Ndf_trPIn, "Ndf_trPIn/I");
  tree->Branch("meanT_trPIn", &meanT_trPIn, "meanT_trPIn/D");
  tree->Branch("sigT_trPIn", &sigT_trPIn, "sigT_trPIn/D");
  tree->Branch("XX0_trPIn", &XX0_trPIn, "XX0_trPIn/D");

  //Event specific
  tree->Branch("NParticle", &NParticle, "NParticle/I");
  tree->Branch("NTrack", &NTrack, "NTrack/I");
  tree->Branch("NVertex", &NVertex, "NVertex/I");
  tree->Branch("trigMask", &trigMask, "trigMask/I");
  tree->Branch("event", &event, "event/L");
  //DY-variables
  tree->Branch("x_beam", &x_beam, "x_beam/D");
  tree->Branch("x_target", &x_target, "x_target/D");
  tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
  tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
    
  //Monte Carlo specific
  //MC Positively charged outgoing muon track parameters
  tree->Branch("pid_MCtr1", &pid_MCtr1, "pid_MCtr1/I");
  tree->Branch("NMCHits_tr1", &NMCHits_tr1, "NMCHits_tr1/I");
  tree->Branch("Pinv_MCtr1", &Pinv_MCtr1, "Pinv_MCtr1/D");
  tree->Branch("theta_MCtr1", &theta_MCtr1, "theta_MCtr1/D");
  tree->Branch("phi_MCtr1", &phi_MCtr1, "phi_MCtr1/D");
  tree->Branch("vMCtr1_X", &vMCtr1_X, "vMCtr1_X/D");
  tree->Branch("vMCtr1_Y", &vMCtr1_Y, "vMCtr1_Y/D");
  tree->Branch("vMCtr1_Z", &vMCtr1_Z, "vMCtr1_Z/D");
  tree->Branch("vMCtr1_E", &vMCtr1_E, "vMCtr1_E/D");
  //MC Negatively charged outgoing muon track parameters
  tree->Branch("pid_MCtr2", &pid_MCtr2, "pid_MCtr2/I");
  tree->Branch("NMCHits_tr2", &NMCHits_tr2, "NMCHits_tr2/I");
  tree->Branch("Pinv_MCtr2", &Pinv_MCtr2, "Pinv_MCtr2/D");
  tree->Branch("theta_MCtr2", &theta_MCtr2, "theta_MCtr2/D");
  tree->Branch("phi_MCtr2", &phi_MCtr2, "phi_MCtr2/D");
  tree->Branch("vMCtr2_X", &vMCtr2_X, "vMCtr2_X/D");
  tree->Branch("vMCtr2_Y", &vMCtr2_Y, "vMCtr2_Y/D");
  tree->Branch("vMCtr2_Z", &vMCtr2_Z, "vMCtr2_Z/D");
  tree->Branch("vMCtr2_E", &vMCtr2_E, "vMCtr2_E/D");
  //MC Beam muon track parameters
  tree->Branch("pid_MCtrIn", &pid_MCtrIn, "pid_MCtrIn/I");
  tree->Branch("NMCHits_trIn", &NMCHits_trIn, "NMCHits_trIn/I");
  tree->Branch("IsBeam_MCtrIn", &IsBeam_MCtrIn, "IsBeam_MCtrIn/I");
  tree->Branch("Pinv_MCtrIn", &Pinv_MCtrIn, "Pinv_MCtrIn/D");
  tree->Branch("theta_MCtrIn", &theta_MCtrIn, "theta_MCtrIn/D");
  tree->Branch("phi_MCtrIn", &phi_MCtrIn, "phi_MCtrIn/D");
  tree->Branch("vMCtrIn_X", &vMCtrIn_X, "vMCtrIn_X/D");
  tree->Branch("vMCtrIn_Y", &vMCtrIn_Y, "vMCtrIn_Y/D");
  tree->Branch("vMCtrIn_Z", &vMCtrIn_Z, "vMCtrIn_Z/D");
  tree->Branch("vMCtrIn_E", &vMCtrIn_E, "vMCtrIn_E/D");
  //MC DY-variables
  tree->Branch("MC_x_beam", &MC_x_beam, "MC_x_beam/D");
  tree->Branch("MC_x_target", &MC_x_target, "MC_x_target/D");
  tree->Branch("MC_x_feynman", &MC_x_feynman, "MC_x_feynman/D");
  tree->Branch("MC_q_transverse", &MC_q_transverse, "MC_q_transverse/D");
}//SetupMCTreeVariables


inline void SetupRealTreeVariables(TTree *tree){
//Positively charged outgoing muon
    tree->Branch("isBeam_p1", &isBeam_p1, "isBeam_p1/I");
    tree->Branch("numVx_p1", &numVx_p1, "numVx_p1/I");
    tree->Branch("numOutVx_p1", &numOutVx_p1, "numOutVx_p1/I");
    tree->Branch("numVxpri_p1", &numVxpri_p1, "numVxpri_p1/I");
    //Positively charnged outgoing muon trajectory parameters at vertex
    tree->Branch("phi_traj1", &phi_traj1, "phi_traj1/D");
    tree->Branch("theta_traj1", &theta_traj1, "theta_traj1/D");
    tree->Branch("qP_traj1", &qP_traj1, "qP_traj1/D");
    tree->Branch("vP1_X", &vP1_X, "vP1_X/D");
    tree->Branch("vP1_Y", &vP1_Y, "vP1_Y/D");
    tree->Branch("vP1_Z", &vP1_Z, "vP1_Z/D");
    tree->Branch("vP1_E", &vP1_E, "vP1_E/D");
    //Negatively charged outgoing muon
    tree->Branch("isBeam_p2", &isBeam_p2, "isBeam_p2/I");
    tree->Branch("numVx_p2", &numVx_p2, "numVx_p2/I");
    tree->Branch("numOutVx_p2", &numOutVx_p2, "numOutVx_p2/I");
    tree->Branch("numVxpri_p2", &numVxpri_p2, "numVxpri_p2/I");
    //Negatively charnged outgoing muon trajectory parameters at vertex
    tree->Branch("phi_traj2", &phi_traj2, "phi_traj2/D");
    tree->Branch("theta_traj2", &theta_traj2, "theta_traj2/D");
    tree->Branch("qP_traj2", &qP_traj2, "qP_traj2/D");
    tree->Branch("vP2_X", &vP2_X, "vP2_X/D");
    tree->Branch("vP2_Y", &vP2_Y, "vP2_Y/D");
    tree->Branch("vP2_Z", &vP2_Z, "vP2_Z/D");
    tree->Branch("vP2_E", &vP2_E, "vP2_E/D");
    //Vertex specific
    tree->Branch("vx_z", &vx_z, "vx_z/D");
    tree->Branch("vx_x", &vx_x, "vx_x/D");
    tree->Branch("vx_y", &vx_y, "vx_y/D");
    tree->Branch("vx_zVar", &vx_zVar, "vx_zVar/D");
    tree->Branch("vx_xVar", &vx_xVar, "vx_xVar/D");
    tree->Branch("vx_xyVar", &vx_xyVar, "vx_xyVar/D");
    tree->Branch("vx_yVar", &vx_yVar, "vx_yVar/D");
    tree->Branch("vx_xzVar", &vx_xzVar, "vx_xzVar/D");
    tree->Branch("vx_yzVar", &vx_yzVar, "vx_yzVar/D");
    tree->Branch("vx_Chi", &vx_Chi, "vx_Chi/D");
    tree->Branch("vx_Chi_ndf", &vx_Chi_ndf, "vx_Chi_ndf/D");
    tree->Branch("vx_ndf", &vx_ndf, "vx_ndf/I");
    tree->Branch("vx_NOutParticles", &vx_NOutParticles, "vx_NOutParticles/I");
    tree->Branch("vx_IsPrimary", &vx_IsPrimary, "vx_IsPrimary/I");
    tree->Branch("vx_IsBestPrimary", &vx_IsBestPrimary, "vx_IsBestPrimary/I");
    //Positively charnged outgoing muon track parameters
    tree->Branch("Zfirst_tr1", &Zfirst_tr1, "Zfirst_tr1/D");
    tree->Branch("Zlast_tr1", &Zlast_tr1, "Zlast_tr1/D");
    tree->Branch("NHits_tr1", &NHits_tr1, "NHits_tr1/I");
    tree->Branch("Chi2tot_tr1", &Chi2tot_tr1, "Chi2tot_tr1/D");
    tree->Branch("Chi2tot_Ndf_tr1", &Chi2tot_Ndf_tr1, "Chi2tot_Ndf_tr1/D");
    tree->Branch("Ndf_tr1", &Ndf_tr1, "Ndf_tr1/I");
    tree->Branch("meanT_tr1", &meanT_tr1, "meanT_tr1/D");
    tree->Branch("sigT_tr1", &sigT_tr1, "sigT_tr1/D");
    tree->Branch("XX0_tr1", &XX0_tr1, "XX0_tr1/D");
    //Negatively charnged outgoing muon track parameters
    tree->Branch("Zfirst_tr2", &Zfirst_tr2, "Zfirst_tr2/D");
    tree->Branch("Zlast_tr2", &Zlast_tr2, "Zlast_tr2/D");
    tree->Branch("NHits_tr2", &NHits_tr2, "NHits_tr2/I");
    tree->Branch("Chi2tot_tr2", &Chi2tot_tr2, "Chi2tot_tr2/D");
    tree->Branch("Chi2tot_Ndf_tr2", &Chi2tot_Ndf_tr2, "Chi2tot_Ndf_tr2/D");
    tree->Branch("Ndf_tr2", &Ndf_tr2, "Ndf_tr2/I");
    tree->Branch("meanT_tr2", &meanT_tr2, "meanT_tr2/D");
    tree->Branch("sigT_tr2", &sigT_tr2, "sigT_tr2/D");
    tree->Branch("XX0_tr2", &XX0_tr2, "XX0_tr2/D");
    //Virtual photon specific/Dimuon
    tree->Branch("vPhoton_X", &vPhoton_X, "vPhoton_X/D");
    tree->Branch("vPhoton_Y", &vPhoton_Y, "vPhoton_Y/D");
    tree->Branch("vPhoton_Z", &vPhoton_Z, "vPhoton_Z/D");
    tree->Branch("vPhoton_E", &vPhoton_E, "vPhoton_E/D");
    tree->Branch("vDiMuon_invM", &vDiMuon_invM, "vDiMuon_invM/D");
    tree->Branch("vOpenAngle", &vOpenAngle, "vOpenAngle/D");

    ////Beam particle specific
    //Beam pion particle specific
    tree->Branch("isBeam_pIn", &isBeam_pIn, "isBeam_pIn/I");
    tree->Branch("numVx_pIn", &numVx_pIn, "numVx_pIn/I");
    tree->Branch("numOutVx_pIn", &numOutVx_pIn, "numOutVx_pIn/I");
    //Beam pion trajectory parameters at vertex
    tree->Branch("phi_trajPIn", &phi_trajPIn, "phi_trajPIn/D");
    tree->Branch("theta_trajPIn", &theta_trajPIn, "theta_trajPIn/D");
    tree->Branch("qP_trajPIn", &qP_trajPIn, "qP_trajPIn/D");
    tree->Branch("beam_X", &beam_X, "beam_X/D");
    tree->Branch("beam_Y", &beam_Y, "beam_Y/D");
    tree->Branch("beam_Z", &beam_Z, "beam_Z/D");
    tree->Branch("beam_E", &beam_E, "beam_E/D");
    //Beam pion track specific
    tree->Branch("Zlast_trPIn", &Zlast_trPIn, "Zlast_trPIn/D");
    tree->Branch("NHits_trPIn", &NHits_trPIn, "NHits_trPIn/I");
    tree->Branch("Chi2tot_trPIn", &Chi2tot_trPIn, "Chi2tot_trPIn/D");
    tree->Branch("Chi2tot_Ndf_trPIn", &Chi2tot_Ndf_trPIn,
		 "Chi2tot_Ndf_trPIn/D");
    tree->Branch("Ndf_trPIn", &Ndf_trPIn, "Ndf_trPIn/I");
    tree->Branch("meanT_trPIn", &meanT_trPIn, "meanT_trPIn/D");
    tree->Branch("sigT_trPIn", &sigT_trPIn, "sigT_trPIn/D");
    tree->Branch("XX0_trPIn", &XX0_trPIn, "XX0_trPIn/D");

    //Event specific
    tree->Branch("NParticle", &NParticle, "NParticle/I");
    tree->Branch("NTrack", &NTrack, "NTrack/I");
    tree->Branch("NVertex", &NVertex, "NVertex/I");
    tree->Branch("trigMask", &trigMask, "trigMask/I");
    tree->Branch("MasterTrigMask", &MasterTrigMask, "MasterTrigMask/I");
    tree->Branch("RunNum", &RunNum, "RunNum/L");
    tree->Branch("SpillNum", &SpillNum, "SpillNum/L");
    tree->Branch("event", &event, "event/L");
    //DY-variables
    tree->Branch("x_beam", &x_beam, "x_beam/D");
    tree->Branch("x_target", &x_target, "x_target/D");
    tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
    tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
    //Target Polarization
    tree->Branch("avgUpStream", &avgUpStream, "avgUpStream/D");
    tree->Branch("avgDownStream", &avgDownStream, "avgDownStream/D");
    tree->Branch("N14_UpStream", &N14_UpStream, "N14_UpStream/D");
    tree->Branch("N14_DownStream", &N14_DownStream, "N14_DownStream/D");
    tree->Branch("upStreamCoil1", &upStreamCoil1, "upStreamCoil1/D");
    tree->Branch("upStreamCoil2", &upStreamCoil2, "upStreamCoil2/D");
    tree->Branch("upStreamCoil3", &upStreamCoil3, "upStreamCoil3/D");
    tree->Branch("upStreamCoil4", &upStreamCoil4, "upStreamCoil4/D");
    tree->Branch("upStreamCoil5", &upStreamCoil5, "upStreamCoil5/D");
    tree->Branch("downStreamCoil6", &downStreamCoil6, "downStreamCoil6/D");
    tree->Branch("downStreamCoil7", &downStreamCoil7, "downStreamCoil7/D");
    tree->Branch("downStreamCoil8", &downStreamCoil8, "downStreamCoil8/D");
    tree->Branch("downStreamCoil9", &downStreamCoil9, "downStreamCoil9/D");
    tree->Branch("downStreamCoil10", &downStreamCoil10,
		      "downStreamCoil10/D");
    tree->Branch("Polarization", &Polarization, "Polarization/D");
    tree->Branch("dilutionFactor", &dilutionFactor, "dilutionFactor/D");
    tree->Branch("error_dilutionFactor", &error_dilutionFactor,
		      "error_dilutionFactor/D");  
}//SetupRealTreeVariables


inline void AssignPhastParVariables(PaEvent &e, const PaParticle &p1,
				    const PaParticle &p2,
				    const PaParticle &pIn){

  //Event specific
  NParticle = e.NParticle();
  NTrack = e.NTrack();
  NVertex = e.NVertex();
  trigMask = e.TrigMask();
  MasterTrigMask = e.MasterTriggerMask();
  event = e.UniqueEvNum();
  RunNum = e.RunNum();
  SpillNum = e.SpillNum();

  //Positively charged outgoing muon
  isBeam_p1 = p1.IsBeam();
  numVx_p1 = p1.NVertex();
  numOutVx_p1 = p1.NOutVertex();

  //Negatively charged outgoing muon
  isBeam_p2 = p2.IsBeam();
  numVx_p2 = p2.NVertex();
  numOutVx_p2 = p2.NOutVertex();

  //Beam pion particle specific
  isBeam_pIn = pIn.IsBeam();
  numVx_pIn = pIn.NVertex();
  numOutVx_pIn = pIn.NOutVertex();
}//AssignPhastParVariables


inline void AssignPhastTrackVariables(const PaVertex &vx, const PaTrack &tr1,
				     const PaTrack &tr2, const PaTrack &trIn){
  //Vertex specific
  vx_z = vx.Z();
  vx_x = vx.X();
  vx_y = vx.Y();
  vx_xVar = vx.Cov(0);
  vx_xyVar = vx.Cov(1);
  vx_yVar = vx.Cov(2);
  vx_xzVar = vx.Cov(3);
  vx_yzVar = vx.Cov(4);
  vx_zVar = vx.Cov(5);
  vx_Chi = vx.Chi2();
  vx_ndf = vx.Ndf();
  vx_Chi_ndf = vx_Chi/(1.0*vx_ndf);
  vx_NOutParticles = vx.NOutParticles();
  vx_IsPrimary = vx.IsPrimary();	
  vx_IsBestPrimary = vx.IsBestPrimary();

  //Positively charged outgoing muon track parameters
  Zfirst_tr1 = tr1.ZFirst();
  Zlast_tr1 = tr1.ZLast();
  NHits_tr1 = tr1.NHits();
  Chi2tot_tr1 = tr1.Chi2tot();
  Ndf_tr1 = tr1.Ndf();
  Chi2tot_Ndf_tr1 = Chi2tot_tr1/(1.0*Ndf_tr1);
  meanT_tr1 = tr1.MeanTime();
  sigT_tr1 = tr1.SigmaTime();
  XX0_tr1 = tr1.XX0();
  
  //Negatively charged outgoing muon track parameters
  Zfirst_tr2 = tr2.ZFirst();
  Zlast_tr2 = tr2.ZLast();
  NHits_tr2 = tr2.NHits();
  Chi2tot_tr2 = tr2.Chi2tot();
  Ndf_tr2 = tr2.Ndf();
  Chi2tot_Ndf_tr2 = Chi2tot_tr2/(2.0*Ndf_tr2);
  meanT_tr2 = tr2.MeanTime();
  sigT_tr2 = tr2.SigmaTime();
  XX0_tr2 = tr2.XX0();

  //Beam pion track specific
  Zlast_trPIn = trIn.ZLast();
  NHits_trPIn = trIn.NHits();
  Chi2tot_trPIn = trIn.Chi2tot();
  Ndf_trPIn = trIn.Ndf();
  Chi2tot_Ndf_trPIn = Chi2tot_trPIn/(1.0*Ndf_trPIn);
  meanT_trPIn = trIn.MeanTime();
  sigT_trPIn = trIn.SigmaTime();
  XX0_trPIn = trIn.XX0();
}//AssignPhastTrackVaribles


inline void AssignPhastTrajVariables(const PaTPar &traj_p1,
				    const PaTPar &traj_p2,
				    const PaTPar &traj_pIn){
  
  //Positively charged outgoing muon trajectory parameters at vertex
  phi_traj1 = traj_p1.Phi();
  theta_traj1 = traj_p1.Theta();
  qP_traj1 = traj_p1.qP();

  //Negatively charged outgoing muon trajectory parameters at vertex
  phi_traj2 = traj_p2.Phi();
  theta_traj2 = traj_p2.Theta();
  qP_traj2 = traj_p2.qP();

  //Beam pion trajectory parameters at vertex
  phi_trajPIn = traj_pIn.Phi();
  theta_trajPIn = traj_pIn.Theta();
  qP_trajPIn = traj_pIn.qP();
}//AssignPhastTrajVaribles


inline void AssignMCtrack(const PaMCtrack &MCtr1, const PaMCtrack &MCtr2,
			  const PaMCtrack &MCtrIn){
  pid_MCtr1 = MCtr1.Pid();
  NMCHits_tr1 = MCtr1.NMCHits();
  Pinv_MCtr1 = MCtr1.Pinv();

  pid_MCtr2 = MCtr2.Pid();
  NMCHits_tr2 = MCtr2.NMCHits();
  Pinv_MCtr2 = MCtr2.Pinv();

  pid_MCtrIn = MCtrIn.Pid();
  NMCHits_trIn = MCtrIn.NMCHits();
  Pinv_MCtrIn = MCtrIn.Pinv();
}//AssignMCtrack


inline void AssignMCtrackVoid(){
  pid_MCtr1 = -999.0;
  NMCHits_tr1 = -999.0;
  Pinv_MCtr1 = -999.0;
  theta_MCtr1 = -999.0;
  phi_MCtr1 = -999.0;
  vMCtr1_X = -999.0;
  vMCtr1_Y = -999.0;
  vMCtr1_Z = -999.0;

  pid_MCtr2 = -999.0;
  NMCHits_tr2 = -999.0;
  Pinv_MCtr2 = -999.0;
  theta_MCtr2 = -999.0;
  phi_MCtr2 = -999.0;
  vMCtr2_X = -999.0;
  vMCtr2_Y = -999.0;
  vMCtr2_Z = -999.0;
	
  pid_MCtrIn = -999.0;
  NMCHits_trIn = -999.0;
  Pinv_MCtrIn = -999.0;
  theta_MCtrIn = -999.0;
  phi_MCtrIn = -999.0;
  vMCtrIn_X = -999.0;
  vMCtrIn_Y = -999.0;
  vMCtrIn_Z = -999.0;

  MC_x_beam = -999.0;
  MC_x_target = -999.0;
  MC_x_feynman = -999.0;
  MC_q_transverse = -999.0;
}//AssignMCtrackVoid


inline void AssignTarget(PaEvent &e){
  //Target Polarization
  vector<float> TargetPol = PaMetaDB::Ref().Polarizations(e.RunNum() );
  Int_t it = 0;
  avgUpStream = TargetPol.at(it); it++;
  avgDownStream = TargetPol.at(it); it++;
  N14_UpStream =  TargetPol.at(it); it++;
  N14_DownStream =  TargetPol.at(it); it = 8;
  upStreamCoil1 = TargetPol.at(it); it++;
  upStreamCoil2 = TargetPol.at(it); it++;
  upStreamCoil3 = TargetPol.at(it); it++;
  upStreamCoil4 = TargetPol.at(it); it++;
  upStreamCoil5 = TargetPol.at(it); it++;
  downStreamCoil6 = TargetPol.at(it); it++;
  downStreamCoil7 = TargetPol.at(it); it++;
  downStreamCoil8 = TargetPol.at(it); it++;
  downStreamCoil9 = TargetPol.at(it); it++;
  downStreamCoil10 = TargetPol.at(it);

  Polarization = PaAlgo::GetDYtargetPolarization(e, vx_z);
  if (vx_z > -294.5 && vx_z < -239.3){
    PaAlgo::GetDYdilutionFactor(e, x_target, q_transverse,
				vDiMuon_invM*vDiMuon_invM, 'U',
				dilutionFactor, error_dilutionFactor);
  }
  else if (vx_z > -219.5 && vx_z < -164.3){
    PaAlgo::GetDYdilutionFactor(e, x_target, q_transverse,
				vDiMuon_invM*vDiMuon_invM, 'D',
				dilutionFactor, error_dilutionFactor);
  }
  else {
    dilutionFactor = 0.0;
    error_dilutionFactor = 0.0;
  }
}//AssignTarget


inline void HistArraySetup(TH1D **hist, TH1D *hCut[][maxCuts], Int_t nBins,
			   Double_t min,
			   Double_t max, Int_t iCut, TString hName){
  
  Phast::Ref().HistFileDir(Form("%s_CutImpact", hName.Data() ) );
  for (Int_t i=0; i<maxCuts; i++) {
    hist[i] = new TH1D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins, min, max);

    hCut[iCut][i] = hist[i];
  }
}//HistArraySetup


inline void Hist2D_ArraySetup(TH2D **hist, TH2D *hCut[][maxCuts],
			      Int_t nBins_x, Double_t min_x, Double_t max_x,
			      Int_t nBins_y, Double_t min_y, Double_t max_y,
			      Int_t iCut, TString hName){
  
  Phast::Ref().HistFileDir(Form("%s_CutImpact", hName.Data() ) );
  for (Int_t i=0; i<maxCuts; i++) {
    hist[i] = new TH2D(Form("hCut%i_%s", i, hName.Data() ),
		       Form("hCut%i_%s", i, hName.Data() ),
		       nBins_x, min_x, max_x,
		       nBins_y, min_y, max_y);

    hCut[iCut][i] = hist[i];
  }
}//Hist2D_ArraySetup


inline void FillCuts(TH1D *hist[][maxCuts], Double_t *variables,
		     Int_t cut){
  for (Int_t i=0; i<nCutHist; i++) hist[i][cut]->Fill(variables[i] );
}//FillCuts


inline void Fill2D_Cuts(TH2D *hist[][maxCuts], Double_t *variables,
			Int_t cut){
  for (Int_t i=0; i<n2D_cutHist; i++) {
    hist[i][cut]->Fill(variables[2*i], variables[2*i+1]);
  }
}//Fill2D_Cuts


#endif


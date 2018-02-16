#include "common.hxx"
#include "functions.h"
#include "setup.h"
#include <bitset>

using namespace std;

int main(int argc, char **argv){
    if(argc < 2){
    cout << "" << endl;
    cout << "To be used with Real Data" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-Pperiod] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
    cout << "Option:  -P period         (which period to take bad spills from)"
	 << endl;
    cout << "(i.e W07, W08...  Can also enter \"WAll\" for all periods)";
    cout << "" << endl;
    cout << "Option:  -u ##		(new UserEvent number, default==420)"
	 << endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  cout << "" << endl;
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  Int_t uflag=0, wflag=0, Qflag=0, fflag=0, Pflag=0;
  Int_t c;
  TString userNum = "", fname = "", outFile = "", period = "";
  
  while ((c = getopt (argc, argv, "wu:f:Q:P:")) != -1) {
    switch (c) {
    case 'u':
      uflag = 1;
      userNum += optarg;
      break;
    case 'w':
      wflag = 1;
      break;
    case 'Q':
      Qflag = 1;
      outFile += optarg;
      break;
    case 'f':
      fflag = 1;
      fname += optarg;
      break;
    case 'P':
      Pflag = 1;
      period += optarg;
      break;
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'P')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'Q')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint (optopt))
	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
	fprintf (stderr,
		 "Unknown option character `\\x%x'.\n",
		 optopt);
      return 1;
    default:
      abort ();
    }
  }

  if(wflag && Qflag){
    cout << "Please only enter -w or -Qoutfile.  Not both" << endl;
    exit(EXIT_FAILURE);
  }

  TString userEvent = "UserEvent";
  if (!uflag) {
    userEvent += "420/Particles";
    cout << "Default UserEvent420 used" << endl;
  }
  else userEvent += userNum + "/Particles";
  TChain* T1 = new TChain(userEvent);
  
  TString BadSpillPath = "src/BadSpills/t3/";//Get badspill information
  map <Long64_t, vector<Long64_t> > BadMap;
  Long64_t BadRun, BadSpill;
  if (!Pflag){
    cout << "Please enter a period for bad spills list" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    BadSpillPath += period;
    ifstream fBadSpills(BadSpillPath+"_BadSpills.txt");
    while(fBadSpills >> BadRun >> BadSpill) BadMap[BadRun].push_back(BadSpill);

    cout << "Data from period:   " << period << endl;
  }

  if (!fflag) {
    cout << "Please enter an input file" << endl;
    exit(EXIT_FAILURE);
  }
  else{
    TFile *f1 = TFile::Open(fname);
    if(!f1){
      cout << fname << " does not exist" <<endl;
      exit(EXIT_FAILURE);
    }
    f1->Close();
  }
  T1->Add(fname);
  cout << "" << endl;

  if (BadMap.size() == 0){
    cout << "Bad spills file did not open" << endl;
    exit(EXIT_FAILURE);    
  }

  //Internal variables and binning
  Double_t M_proton = 0.938272;

  //Binnings
  Double_t xN_bounds[] = {0.00, 0.13, 0.19, 1.00};
  Double_t xPi_bounds[] = {0.00, 0.40, 0.56, 1.00};
  Double_t xF_bounds[] = {-1.0, 0.21, 0.41, 1.00};
  Double_t pT_bounds[] = {0.4, 0.9, 1.4, 5.0};
  //Double_t M_bounds[] = {4.30, 4.75, 5.50, 8.50};
  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "JPsi mass bounds" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << " " << endl;
  Double_t M_bounds[] = {2.5, 3.2, 3.8, 4.30};//JPsi
  
  TVectorD tv_xN_bounds; tv_xN_bounds.Use(4, xN_bounds);
  TVectorD tv_xPi_bounds; tv_xPi_bounds.Use(4, xPi_bounds);
  TVectorD tv_xF_bounds; tv_xF_bounds.Use(4, xF_bounds);
  TVectorD tv_pT_bounds; tv_pT_bounds.Use(4, pT_bounds);
  TVectorD tv_M_bounds; tv_M_bounds.Use(4, M_bounds);

  //Averages
  Double_t AvgPolarization=0.0, AvgDilution=0.0, AvgDilution_corrected=0.0;
  Int_t AvgPolarization_count=0, AvgDilution_count=0;
  
  Double_t AvgPol_xN[3]={0.0}, AvgDil_xN[3]={0.0};
  Int_t AvgPol_xN_count[3]={0}, AvgDil_xN_count[3]={0};
  Double_t AvgPol_xN_UpStream[3]={0.0}, AvgDil_xN_UpStream[3]={0.0};
  Int_t AvgPol_xN_count_UpStream[3]={0}, AvgDil_xN_count_UpStream[3]={0};
  Double_t AvgPol_xN_DownStream[3]={0.0}, AvgDil_xN_DownStream[3]={0.0};
  Int_t AvgPol_xN_count_DownStream[3]={0}, AvgDil_xN_count_DownStream[3]={0};
  
  Double_t AvgPol_xPi[3]={0.0}, AvgDil_xPi[3]={0.0};
  Int_t AvgPol_xPi_count[3]={0}, AvgDil_xPi_count[3]={0};
  Double_t AvgPol_xPi_UpStream[3]={0.0}, AvgDil_xPi_UpStream[3]={0.0};
  Int_t AvgPol_xPi_count_UpStream[3]={0}, AvgDil_xPi_count_UpStream[3]={0};
  Double_t AvgPol_xPi_DownStream[3]={0.0}, AvgDil_xPi_DownStream[3]={0.0};
  Int_t AvgPol_xPi_count_DownStream[3]={0}, AvgDil_xPi_count_DownStream[3]={0};
  
  Double_t AvgPol_xF[3]={0.0}, AvgDil_xF[3]={0.0};
  Int_t AvgPol_xF_count[3]={0}, AvgDil_xF_count[3]={0};
  Double_t AvgPol_xF_UpStream[3]={0.0}, AvgDil_xF_UpStream[3]={0.0};
  Int_t AvgPol_xF_count_UpStream[3]={0}, AvgDil_xF_count_UpStream[3]={0};
  Double_t AvgPol_xF_DownStream[3]={0.0}, AvgDil_xF_DownStream[3]={0.0};
  Int_t AvgPol_xF_count_DownStream[3]={0}, AvgDil_xF_count_DownStream[3]={0};

  Double_t AvgPol_pT[3]={0.0}, AvgDil_pT[3]={0.0};
  Int_t AvgPol_pT_count[3]={0}, AvgDil_pT_count[3]={0};
  Double_t AvgPol_pT_UpStream[3]={0.0}, AvgDil_pT_UpStream[3]={0.0};
  Int_t AvgPol_pT_count_UpStream[3]={0}, AvgDil_pT_count_UpStream[3]={0};
  Double_t AvgPol_pT_DownStream[3]={0.0}, AvgDil_pT_DownStream[3]={0.0};
  Int_t AvgPol_pT_count_DownStream[3]={0}, AvgDil_pT_count_DownStream[3]={0};

  Double_t AvgPol_M[3]={0.0}, AvgDil_M[3]={0.0};
  Int_t AvgPol_M_count[3]={0}, AvgDil_M_count[3]={0};
  Double_t AvgPol_M_UpStream[3]={0.0}, AvgDil_M_UpStream[3]={0.0};
  Int_t AvgPol_M_count_UpStream[3]={0}, AvgDil_M_count_UpStream[3]={0};
  Double_t AvgPol_M_DownStream[3]={0.0}, AvgDil_M_DownStream[3]={0.0};
  Int_t AvgPol_M_count_DownStream[3]={0}, AvgDil_M_count_DownStream[3]={0};
  
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
  Int_t trigMask, MasterTrigMask;
  Long64_t event, RunNum, SpillNum;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse;
  //Target Polarization
  Double_t avgUpStream, avgDownStream;
  Double_t N14_UpStream, N14_DownStream;
  Double_t upStreamCoil1, upStreamCoil2, upStreamCoil3, upStreamCoil4;
  Double_t upStreamCoil5;
  Double_t downStreamCoil6, downStreamCoil7, downStreamCoil8;
  Double_t downStreamCoil9, downStreamCoil10;
  Double_t Polarization, dilutionFactor, error_dilutionFactor;
  //Positions
  Double_t SM1_p1x, SM1_p1y, SM1_p2x, SM1_p2y;
  Double_t SM2_p1x, SM2_p1y, SM2_p2x, SM2_p2y;
  Double_t HG01_p1x, HG01_p1y, HG02_y1_p1x, HG02_y1_p1y;
  Double_t HG02_y2_p1x, HG02_y2_p1y;
  Double_t HG01_p2x, HG01_p2y, HG02_y1_p2x, HG02_y1_p2y;
  Double_t HG02_y2_p2x, HG02_y2_p2y;

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
  T1->SetBranchAddress("MasterTrigMask", &MasterTrigMask);
  T1->SetBranchAddress("RunNum", &RunNum);
  T1->SetBranchAddress("SpillNum", &SpillNum);
  T1->SetBranchAddress("event", &event);
  //DY-variables
  T1->SetBranchAddress("x_beam", &x_beam);
  T1->SetBranchAddress("x_target", &x_target);
  T1->SetBranchAddress("x_feynman", &x_feynman);
  T1->SetBranchAddress("q_transverse", &q_transverse);
  //Target Polarization
  T1->SetBranchAddress("avgUpStream", &avgUpStream);
  T1->SetBranchAddress("avgDownStream", &avgDownStream);
  T1->SetBranchAddress("N14_UpStream", &N14_UpStream);
  T1->SetBranchAddress("N14_DownStream", &N14_DownStream);
  T1->SetBranchAddress("upStreamCoil1", &upStreamCoil1);
  T1->SetBranchAddress("upStreamCoil2", &upStreamCoil2);
  T1->SetBranchAddress("upStreamCoil3", &upStreamCoil3);
  T1->SetBranchAddress("upStreamCoil4", &upStreamCoil4);
  T1->SetBranchAddress("upStreamCoil5", &upStreamCoil5);
  T1->SetBranchAddress("downStreamCoil6", &downStreamCoil6);
  T1->SetBranchAddress("downStreamCoil7", &downStreamCoil7);
  T1->SetBranchAddress("downStreamCoil8", &downStreamCoil8);
  T1->SetBranchAddress("downStreamCoil9", &downStreamCoil9);
  T1->SetBranchAddress("downStreamCoil10", &downStreamCoil10);
  T1->SetBranchAddress("Polarization", &Polarization);
  T1->SetBranchAddress("dilutionFactor", &dilutionFactor);
  T1->SetBranchAddress("error_dilutionFactor", &error_dilutionFactor);
  //Positions
  T1->SetBranchAddress("SM1_p1x", &SM1_p1x);
  T1->SetBranchAddress("SM1_p1y", &SM1_p1y);
  T1->SetBranchAddress("SM1_p2x", &SM1_p2x);
  T1->SetBranchAddress("SM1_p2y", &SM1_p2y);
  T1->SetBranchAddress("SM2_p1x", &SM2_p1x);
  T1->SetBranchAddress("SM2_p1y", &SM2_p1y);
  T1->SetBranchAddress("SM2_p2x", &SM2_p2x);
  T1->SetBranchAddress("SM2_p2y", &SM2_p2y);
  T1->SetBranchAddress("HG01_p1x", &HG01_p1x);
  T1->SetBranchAddress("HG01_p1y", &HG01_p1y);
  T1->SetBranchAddress("HG01_p2x", &HG01_p2x);
  T1->SetBranchAddress("HG01_p2y", &HG01_p2y);
  T1->SetBranchAddress("HG02_y1_p1x", &HG02_y1_p1x);
  T1->SetBranchAddress("HG02_y1_p1y", &HG02_y1_p1y);
  T1->SetBranchAddress("HG02_y1_p2x", &HG02_y1_p2x);
  T1->SetBranchAddress("HG02_y1_p2y", &HG02_y1_p2y);
  T1->SetBranchAddress("HG02_y2_p1x", &HG02_y2_p1x);
  T1->SetBranchAddress("HG02_y2_p1y", &HG02_y2_p1y);
  T1->SetBranchAddress("HG02_y2_p2x", &HG02_y2_p2x);
  T1->SetBranchAddress("HG02_y2_p2y", &HG02_y2_p2y);


  //Cut histograms
  const Int_t nCutHist = 11;//Number of impact cut hist
  TH1D* hCuts = new TH1D("hCuts", "hCuts", 200, 0, 200);
  Int_t cut_bin = 1, cut_space = 10;
  
  TH1D *hCut_VxZ[nRealCuts];
  TH1D *hCut_MuPTheta[nRealCuts],*hCut_MuPPhi[nRealCuts],*hCut_MuPqP[nRealCuts];
  TH1D *hCut_MuMTheta[nRealCuts],*hCut_MuMPhi[nRealCuts],*hCut_MuMqP[nRealCuts];
  TH1D *hCut_xN[nRealCuts], *hCut_xPi[nRealCuts], *hCut_xF[nRealCuts];
  TH1D *hCut_qT[nRealCuts];

  TH1D *hImpactCuts[nCutHist][nRealCuts];
  
  Int_t ih = 0;
  HistArraySetupReal(hCut_VxZ, hImpactCuts, 500, -500, 100, ih, "VxZ"); ih++;
  HistArraySetupReal(hCut_MuPTheta, hImpactCuts, 100, 0, 0.3, ih, "MuPTheta");
  ih++;
  HistArraySetupReal(hCut_MuPPhi, hImpactCuts, 100, -TMath::Pi(), TMath::Pi(),
		     ih, "MuPPhi"); ih++;
  HistArraySetupReal(hCut_MuPqP, hImpactCuts, 200, 0, 200, ih, "MuPqP"); ih++;
  HistArraySetupReal(hCut_MuMTheta, hImpactCuts, 100, 0, 0.3, ih, "MuMTheta");
  ih++;
  HistArraySetupReal(hCut_MuMPhi, hImpactCuts, 100, -TMath::Pi(), TMath::Pi(),
		     ih, "MuMPhi"); ih++;
  HistArraySetupReal(hCut_MuMqP, hImpactCuts, 200, -200, 0, ih, "MuMqP"); ih++;
  HistArraySetupReal(hCut_xN, hImpactCuts, 100, 0, 1, ih, "xN"); ih++;
  HistArraySetupReal(hCut_xPi, hImpactCuts, 100, 0, 1, ih, "xPi"); ih++;
  HistArraySetupReal(hCut_xF, hImpactCuts, 200, -1, 1, ih, "xF"); ih++;
  HistArraySetupReal(hCut_qT, hImpactCuts, 200, 0, 5, ih, "qT"); ih++;

  TTree *tree = new TTree("pT_Weighted", "pT_Weighted");
  Double_t PhiS, PhiS_simple, Phi_CS, Theta_CS;
  Int_t targetPosition;
  Double_t Spin[7];
  tree->Branch("PhiS", &PhiS, "PhiS/D");
  tree->Branch("PhiS_simple", &PhiS_simple, "PhiS_simple/D");
  tree->Branch("Phi_CS", &Phi_CS, "Phi_CS/D");
  tree->Branch("Theta_CS", &Theta_CS, "Theta_CS/D");
  tree->Branch("x_beam", &x_beam, "x_beam/D");
  tree->Branch("x_target", &x_target, "x_target/D");
  tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
  tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
  tree->Branch("Mmumu", &vDiMuon_invM, "Mmumu/D");
  for (Int_t i=0; i<7; i++) {
    tree->Branch(Form("Spin_%i", i), &Spin[i], Form("Spin_%i/D", i));
  }
  tree->Branch("targetPosition", &targetPosition, "targetPosition/I");
  tree->Branch("theta_traj1", &theta_traj1, "theta_traj1/D");
  tree->Branch("phi_traj1", &phi_traj1, "phi_traj1/D");
  tree->Branch("qP_traj1", &qP_traj1, "qP_traj1/D");
  tree->Branch("theta_traj2", &theta_traj2, "theta_traj2/D");
  tree->Branch("phi_traj2", &phi_traj2, "phi_traj2/D");
  tree->Branch("qP_traj2", &qP_traj2, "qP_traj2/D");
  tree->Branch("theta_trajPIn", &theta_trajPIn, "theta_trajPIn/D");
  tree->Branch("phi_trajPIn", &phi_trajPIn, "phi_trajPIn/D");
  tree->Branch("qP_trajPIn", &qP_trajPIn, "qP_trajPIn/D");
  tree->Branch("vx_z", &vx_z, "vx_z/D");
  tree->Branch("vx_x", &vx_x, "vx_x/D");
  tree->Branch("vx_y", &vx_y, "vx_y/D");
  tree->Branch("SM1_p1x", &SM1_p1x, "SM1_p1x/D");
  tree->Branch("SM1_p1y", &SM1_p1y, "SM1_p1y/D");
  tree->Branch("SM1_p2x", &SM1_p2x, "SM1_p2x/D");
  tree->Branch("SM1_p2y", &SM1_p2y, "SM1_p2y/D");
  tree->Branch("SM2_p1x", &SM2_p1x, "SM2_p1x/D");
  tree->Branch("SM2_p1y", &SM2_p1y, "SM2_p1y/D");
  tree->Branch("SM2_p2x", &SM2_p2x, "SM2_p2x/D");
  tree->Branch("SM2_p2y", &SM2_p2y, "SM2_p2y/D");
  tree->Branch("HG01_p1x", &HG01_p1x, "HG01_p1x/D");
  tree->Branch("HG01_p1y", &HG01_p1y, "HG01_p1y/D");
  tree->Branch("HG01_p2x", &HG01_p2x, "HG01_p2x/D");
  tree->Branch("HG01_p2y", &HG01_p2y, "HG01_p2y/D");
  tree->Branch("HG02_y1_p1x", &HG02_y1_p1x, "HG02_y1_p1x/D");
  tree->Branch("HG02_y1_p1y", &HG02_y1_p1y, "HG02_y1_p1y/D");
  tree->Branch("HG02_y1_p2x", &HG02_y1_p2x, "HG02_y1_p2x/D");
  tree->Branch("HG02_y1_p2y", &HG02_y1_p2y, "HG02_y1_p2y/D");
  tree->Branch("HG02_y2_p1x", &HG02_y2_p1x, "HG02_y2_p1x/D");
  tree->Branch("HG02_y2_p1y", &HG02_y2_p1y, "HG02_y2_p1y/D");
  tree->Branch("HG02_y2_p2x", &HG02_y2_p2x, "HG02_y2_p2x/D");
  tree->Branch("HG02_y2_p2y", &HG02_y2_p2y, "HG02_y2_p2y/D");  

  Int_t tree_entries = T1->GetEntries();
  //Int_t tree_entries = 1000;//Debug
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  for (Int_t i=0; i<tree_entries; i++){
    T1->GetEntry(i, 0);
    
    //Cuts
    cut_bin = 1;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//All Data
    
    Double_t cut_variables[nCutHist] = {vx_z, theta_traj1,
					phi_traj1, qP_traj1, theta_traj2,
					phi_traj2,
					qP_traj2, x_beam, x_target,
					x_feynman, q_transverse};
    Int_t icut = 0;
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    map<Long64_t, vector<Long64_t> >::iterator it = BadMap.find(RunNum);
    if (it != BadMap.end() ){
      if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(), -310)
	  != BadMap[it->first].end() ) continue;
      else if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(),
			 SpillNum) != BadMap[it->first].end() ) continue;
    }
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Bad spills
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    TLorentzVector lv_p1_Mu(vP1_X, vP1_Y, vP1_Z, vP1_E);
    TLorentzVector lv_p2_Mu(vP2_X, vP2_Y, vP2_Z, vP2_E);
    TLorentzVector lv_diMu = lv_p1_Mu + lv_p2_Mu;
    TLorentzVector lv_target_1 (0, 0, 0, M_proton);

    if (x_beam < 0.0 || x_beam > 1.0) continue;
    if (x_target < 0.0 || x_target > 1.0) continue;
    if (x_feynman < -1.0 || x_feynman > 1.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Physical Kinematics
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if (q_transverse < 0.4 || q_transverse > 5.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//qT cuts
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if ( (vx_z < -294.5 || vx_z > -239.3) && (vx_z < -219.5 || vx_z > -164.3)
	 ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target z-cut
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if(TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >= TMath::Power(1.9, 2)
       ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target radial cut
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;
        
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

    AvgDilution += TMath::Abs(dilutionFactor);
    AvgDilution_corrected += 0.95*TMath::Abs(dilutionFactor);
    AvgDilution_count++;

    Double_t correct_dil = 0.95*TMath::Abs(dilutionFactor);
    BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
    BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
    BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
    BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
    BinAvg(AvgDil_M, AvgDil_M_count, vDiMuon_invM, M_bounds, correct_dil);
    if (vx_z >= -294.5 && vx_z <= -239.3){//Up stream NH3
      Spin[0] = -1.0*avgUpStream/(TMath::Abs(avgUpStream) );
      Spin[1] = avgUpStream;
      Spin[2] = upStreamCoil1;
      Spin[3] = upStreamCoil2;
      Spin[4] = upStreamCoil3;
      Spin[5] = upStreamCoil4;
      Spin[6] = upStreamCoil5;
      targetPosition = 0;

      AvgDilution += TMath::Abs(dilutionFactor);
      AvgDilution_corrected += 0.95*TMath::Abs(dilutionFactor);
      AvgDilution_count++;

      Double_t correct_dil = 0.95*TMath::Abs(dilutionFactor);
      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, vDiMuon_invM, M_bounds, correct_dil);

      BinAvg(AvgDil_xN_UpStream, AvgDil_xN_count_UpStream, x_target,
	     xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi_UpStream, AvgDil_xPi_count_UpStream, x_beam,
	     xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF_UpStream, AvgDil_xF_count_UpStream, x_feynman,
	     xF_bounds, correct_dil);
      BinAvg(AvgDil_pT_UpStream, AvgDil_pT_count_UpStream, q_transverse,
	     pT_bounds, correct_dil);
      BinAvg(AvgDil_M_UpStream, AvgDil_M_count_UpStream, vDiMuon_invM,
	     M_bounds, correct_dil);

      Double_t pol = TMath::Abs(Polarization);
      BinAvg(AvgPol_xN_UpStream, AvgPol_xN_count_UpStream, x_target,
	     xN_bounds, pol);
      BinAvg(AvgPol_xPi_UpStream, AvgPol_xPi_count_UpStream, x_beam,
	     xPi_bounds, pol);
      BinAvg(AvgPol_xF_UpStream, AvgPol_xF_count_UpStream, x_feynman,
	     xF_bounds, pol);
      BinAvg(AvgPol_pT_UpStream, AvgPol_pT_count_UpStream, q_transverse,
	     pT_bounds, pol);
      BinAvg(AvgPol_M_UpStream, AvgPol_M_count_UpStream, vDiMuon_invM,
	     M_bounds, pol);
    }//Up stream
    else if (vx_z >= -219.5 && vx_z <= -164.3){//Down stream NH3
      Spin[0] = -1.0*avgDownStream/(TMath::Abs(avgDownStream) );
      Spin[1] = avgDownStream;
      Spin[2] = downStreamCoil6;
      Spin[3] = downStreamCoil7;
      Spin[4] = downStreamCoil8;
      Spin[5] = downStreamCoil9;
      Spin[6] = downStreamCoil10;
      targetPosition = 1;

      AvgDilution += TMath::Abs(dilutionFactor);
      AvgDilution_corrected += 0.91*TMath::Abs(dilutionFactor);
      AvgDilution_count++;

      Double_t correct_dil = 0.91*TMath::Abs(dilutionFactor);
      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, vDiMuon_invM, M_bounds, correct_dil);

      BinAvg(AvgDil_xN_DownStream, AvgDil_xN_count_DownStream, x_target,
	     xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi_DownStream, AvgDil_xPi_count_DownStream, x_beam,
	     xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF_DownStream, AvgDil_xF_count_DownStream, x_feynman,
	     xF_bounds, correct_dil);
      BinAvg(AvgDil_pT_DownStream, AvgDil_pT_count_DownStream, q_transverse,
	     pT_bounds, correct_dil);
      BinAvg(AvgDil_M_DownStream, AvgDil_M_count_DownStream, vDiMuon_invM,
	     M_bounds, correct_dil);

      Double_t pol = TMath::Abs(Polarization);
      BinAvg(AvgPol_xN_DownStream, AvgPol_xN_count_DownStream, x_target,
	     xN_bounds, pol);
      BinAvg(AvgPol_xPi_DownStream, AvgPol_xPi_count_DownStream, x_beam,
	     xPi_bounds, pol);
      BinAvg(AvgPol_xF_DownStream, AvgPol_xF_count_DownStream, x_feynman,
	     xF_bounds, pol);
      BinAvg(AvgPol_pT_DownStream, AvgPol_pT_count_DownStream, q_transverse,
	     pT_bounds, pol);
      BinAvg(AvgPol_M_DownStream, AvgPol_M_count_DownStream, vDiMuon_invM,
	     M_bounds, pol);
    }//Down Stream

    //Setup Vectors in different coordinate systems
    //Compass frame:
    TLorentzVector lv_beam(beam[0], beam[1], beam[2], beam[3]);
    TLorentzVector lv_target(0, 0, 0, target[0]);
    //TLorentzVector lv_Spin(0, Spin[0], 0, 0);
    TLorentzVector lv_Spin(0, Spin[0], 0, 0);
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

    //Average
    AvgPolarization += TMath::Abs(Polarization);
    AvgPolarization_count++;
    if(Polarization == 0.0 || dilutionFactor == 0.0){
      cout << "Problems with Polarization or dilution value" << endl;
    }
    
    Double_t pol = TMath::Abs(Polarization);
    BinAvg(AvgPol_xN, AvgPol_xN_count, x_target, xN_bounds, pol);
    BinAvg(AvgPol_xPi, AvgPol_xPi_count, x_beam, xPi_bounds, pol);
    BinAvg(AvgPol_xF, AvgPol_xF_count, x_feynman, xF_bounds, pol);
    BinAvg(AvgPol_pT, AvgPol_pT_count, q_transverse, pT_bounds, pol);
    BinAvg(AvgPol_M, AvgPol_M_count, vDiMuon_invM, M_bounds, pol);

    tree->Fill();
  }//tree entries

  TVectorD Dil_xN(3), Dil_xN_UpStream(3), Dil_xN_DownStream(3);
  TVectorD Pol_xN(3), Pol_xN_UpStream(3), Pol_xN_DownStream(3);
  TVectorD Dil_xPi(3), Dil_xPi_UpStream(3), Dil_xPi_DownStream(3);
  TVectorD Pol_xPi(3), Pol_xPi_UpStream(3), Pol_xPi_DownStream(3);
  TVectorD Dil_xF(3), Dil_xF_UpStream(3), Dil_xF_DownStream(3);
  TVectorD Pol_xF(3), Pol_xF_UpStream(3), Pol_xF_DownStream(3);
  TVectorD Dil_pT(3), Dil_pT_UpStream(3), Dil_pT_DownStream(3);
  TVectorD Pol_pT(3), Pol_pT_UpStream(3), Pol_pT_DownStream(3);
  TVectorD Dil_M(3), Dil_M_UpStream(3), Dil_M_DownStream(3);
  TVectorD Pol_M(3), Pol_M_UpStream(3), Pol_M_DownStream(3);
  for (Int_t i=0; i<3; i++) {
    //Dilution
    ///////////////
    Dil_xN[i] = AvgDil_xN[i]/AvgDil_xN_count[i];
    Dil_xPi[i] = AvgDil_xPi[i]/AvgDil_xPi_count[i];
    Dil_xF[i] = AvgDil_xF[i]/AvgDil_xF_count[i];
    Dil_pT[i] = AvgDil_pT[i]/AvgDil_pT_count[i];
    Dil_M[i] = AvgDil_M[i]/AvgDil_M_count[i];

    //UpStream
    Dil_xN_UpStream[i] = AvgDil_xN_UpStream[i]/AvgDil_xN_count_UpStream[i];
    Dil_xPi_UpStream[i] = AvgDil_xPi_UpStream[i]/AvgDil_xPi_count_UpStream[i];
    Dil_xF_UpStream[i] = AvgDil_xF_UpStream[i]/AvgDil_xF_count_UpStream[i];
    Dil_pT_UpStream[i] = AvgDil_pT_UpStream[i]/AvgDil_pT_count_UpStream[i];
    Dil_M_UpStream[i] = AvgDil_M_UpStream[i]/AvgDil_M_count_UpStream[i];

    //DownStream
    Dil_xN_DownStream[i] =AvgDil_xN_DownStream[i]/AvgDil_xN_count_DownStream[i];
    Dil_xPi_DownStream[i]=
      AvgDil_xPi_DownStream[i]/AvgDil_xPi_count_DownStream[i];
    Dil_xF_DownStream[i] =AvgDil_xF_DownStream[i]/AvgDil_xF_count_DownStream[i];
    Dil_pT_DownStream[i] =AvgDil_pT_DownStream[i]/AvgDil_pT_count_DownStream[i];
    Dil_M_DownStream[i] = AvgDil_M_DownStream[i]/AvgDil_M_count_DownStream[i];

    //Polarization
    ////////////////
    Pol_xN[i] = AvgPol_xN[i]/AvgPol_xN_count[i];
    Pol_xPi[i] = AvgPol_xPi[i]/AvgPol_xPi_count[i];
    Pol_xF[i] = AvgPol_xF[i]/AvgPol_xF_count[i];
    Pol_pT[i] = AvgPol_pT[i]/AvgPol_pT_count[i];
    Pol_M[i] = AvgPol_M[i]/AvgPol_M_count[i];

    //UpStream
    Pol_xN_UpStream[i] = AvgPol_xN_UpStream[i]/AvgPol_xN_count_UpStream[i];
    Pol_xPi_UpStream[i] = AvgPol_xPi_UpStream[i]/AvgPol_xPi_count_UpStream[i];
    Pol_xF_UpStream[i] = AvgPol_xF_UpStream[i]/AvgPol_xF_count_UpStream[i];
    Pol_pT_UpStream[i] = AvgPol_pT_UpStream[i]/AvgPol_pT_count_UpStream[i];
    Pol_M_UpStream[i] = AvgPol_M_UpStream[i]/AvgPol_M_count_UpStream[i];

    //DownStream
    Pol_xN_DownStream[i] =AvgPol_xN_DownStream[i]/AvgPol_xN_count_DownStream[i];
    Pol_xPi_DownStream[i]=
      AvgPol_xPi_DownStream[i]/AvgPol_xPi_count_DownStream[i];
    Pol_xF_DownStream[i] =AvgPol_xF_DownStream[i]/AvgPol_xF_count_DownStream[i];
    Pol_pT_DownStream[i] =AvgPol_pT_DownStream[i]/AvgPol_pT_count_DownStream[i];
    Pol_M_DownStream[i] = AvgPol_M_DownStream[i]/AvgPol_M_count_DownStream[i];
  }
  TVectorD Dil_int(1), Pol_int(1);
  Dil_int[0] = AvgDilution_corrected/AvgDilution_count;
  Pol_int[0] = AvgPolarization/AvgPolarization_count;
  
  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  //Cuts histogram
  TString cutNames[nRealCuts] = {"AllData", "GoodSpills", "xPion,xN,xF",
				 "0.4<qT<5",
				 "TargetZ-cut", "TargetRadius"};
  for (Int_t i=0, j=1; i<nRealCuts; i++, j+=cut_space){
    Int_t bin_index = hCuts->GetXaxis()->FindBin(j);
    hCuts->GetXaxis()->SetBinLabel(bin_index, cutNames[i]);
  }

  if (!wflag && !Qflag) cout << "No file output" << endl;
  else{
    if (wflag) outFile += "RealData.root";
    TFile *myFile = new TFile(outFile, "RECREATE");
    hCuts->Write();
    tree->Write();

    tv_xN_bounds.Write("tv_xN_bounds");
    tv_xPi_bounds.Write("tv_xPi_bounds");
    tv_xF_bounds.Write("tv_xF_bounds");
    tv_pT_bounds.Write("tv_pT_bounds");
    tv_M_bounds.Write("tv_M_bounds");

    Dil_int.Write("Dil_int");
    Dil_xN.Write("Dil_xN");
    Dil_xPi.Write("Dil_xPi");
    Dil_xF.Write("Dil_xF");
    Dil_pT.Write("Dil_pT");
    Dil_M.Write("Dil_M");

    Dil_xN_UpStream.Write("Dil_xN_UpStream");
    Dil_xPi_UpStream.Write("Dil_xPi_UpStream");
    Dil_xF_UpStream.Write("Dil_xF_UpStream");
    Dil_pT_UpStream.Write("Dil_pT_UpStream");
    Dil_M_UpStream.Write("Dil_M_UpStream");

    Dil_xN_DownStream.Write("Dil_xN_DownStream");
    Dil_xPi_DownStream.Write("Dil_xPi_DownStream");
    Dil_xF_DownStream.Write("Dil_xF_DownStream");
    Dil_pT_DownStream.Write("Dil_pT_DownStream");
    Dil_M_DownStream.Write("Dil_M_DownStream");

    Pol_int.Write("Pol_int");
    Pol_xN.Write("Pol_xN");
    Pol_xPi.Write("Pol_xPi");
    Pol_xF.Write("Pol_xF");
    Pol_pT.Write("Pol_pT");
    Pol_M.Write("Pol_M");

    Pol_xN_UpStream.Write("Pol_xN_UpStream");
    Pol_xPi_UpStream.Write("Pol_xPi_UpStream");
    Pol_xF_UpStream.Write("Pol_xF_UpStream");
    Pol_pT_UpStream.Write("Pol_pT_UpStream");
    Pol_M_UpStream.Write("Pol_M_UpStream");

    Pol_xN_DownStream.Write("Pol_xN_DownStream");
    Pol_xPi_DownStream.Write("Pol_xPi_DownStream");
    Pol_xF_DownStream.Write("Pol_xF_DownStream");
    Pol_pT_DownStream.Write("Pol_pT_DownStream");
    Pol_M_DownStream.Write("Pol_M_DownStream");
    
    TDirectory *VxZ_CutImpact = myFile->mkdir("VxZ_CutImpact");
    TDirectory *MuPTheta_CutImpact = myFile->mkdir("MuPTheta_CutImpact");
    TDirectory *MuPPhi_CutImpact = myFile->mkdir("MuPPhi_CutImpact");
    TDirectory *MuPqP_CutImpact = myFile->mkdir("MuPqP_CutImpact");
    TDirectory *MuMTheta_CutImpact = myFile->mkdir("MuMTheta_CutImpact");
    TDirectory *MuMPhi_CutImpact = myFile->mkdir("MuMPhi_CutImpact");
    TDirectory *MuMqP_CutImpact = myFile->mkdir("MuMqP_CutImpact");
    TDirectory *xN_CutImpact = myFile->mkdir("xN_CutImpact");
    TDirectory *xPi_CutImpact = myFile->mkdir("xPi_CutImpact");
    TDirectory *xF_CutImpact = myFile->mkdir("xF_CutImpact");
    TDirectory *qT_CutImpact = myFile->mkdir("qT_CutImpact");
    for (Int_t i=0; i<nRealCuts; i++) {
      VxZ_CutImpact->cd();
      hCut_VxZ[i]->Write(cutNames[i]);

      MuPTheta_CutImpact->cd();
      hCut_MuPTheta[i]->Write(cutNames[i]);
      MuPPhi_CutImpact->cd();
      hCut_MuPPhi[i]->Write(cutNames[i]);
      MuPqP_CutImpact->cd();
      hCut_MuPqP[i]->Write(cutNames[i]);

      MuMTheta_CutImpact->cd();
      hCut_MuMTheta[i]->Write(cutNames[i]);
      MuMPhi_CutImpact->cd();
      hCut_MuMPhi[i]->Write(cutNames[i]);
      MuMqP_CutImpact->cd();
      hCut_MuMqP[i]->Write(cutNames[i]);

      xN_CutImpact->cd();
      hCut_xN[i]->Write(cutNames[i]);
      xPi_CutImpact->cd();
      hCut_xPi[i]->Write(cutNames[i]);
      xF_CutImpact->cd();
      hCut_xF[i]->Write(cutNames[i]);
      qT_CutImpact->cd();
      hCut_qT[i]->Write(cutNames[i]);
    }

    cout << myFile->GetName() << " was written" << endl;
    myFile->Close();
  }
    
  theApp.Run();//Needed to make root graphics work on C++
}//main

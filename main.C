#include "common.hxx"
#include "functions.h"
#include "setup.h"
#include <bitset>

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "" << endl;
    cout << "To be used with Monte Carlo Data" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
    cout << "Option:  -u ##		(new UserEvent number, default==420)"
	 << endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "Option:  -b textfile with binning information	";
    cout << "(textfile should be made from Macro/Binning/avgBinBounds.C)"
	 << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t uflag=0, wflag=0, Qflag=0, fflag=0, binFlag=0;
  Int_t c;
  TString userNum = "", fname = "", outFile = "", binFile="";
  
  while ((c = getopt (argc, argv, "wb:u:f:Q:")) != -1) {
    switch (c) {
    case 'u':
      uflag = 1;
      userNum += optarg;
      break;
    case 'b':
      binFlag = 1;
      binFile += optarg;
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
      cout << fname << endl;
      break;
    case '?':
      if (optopt == 'u')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'Q')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'b')
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

  //Basic Checks
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

  vector<Double_t> xN_bounds, xN_xval; 
  vector<Double_t> xPi_bounds, xPi_xval;
  vector<Double_t> xF_bounds, xF_xval; 
  vector<Double_t> pT_bounds, pT_xval;
  vector<Double_t> M_bounds, M_xval;
  vector<Double_t> rad_bounds, rad_xval;
  vector<Double_t> vxZ_upstream_bounds, vxZ_upstream_xval;
  vector<Double_t> vxZ_downstream_bounds, vxZ_downstream_xval;
  xN_bounds.push_back(0.0);
  xPi_bounds.push_back(0.0);
  xF_bounds.push_back(-1.0);
  pT_bounds.push_back(0.4);
  M_bounds.push_back(4.3);
  rad_bounds.push_back(0.0);
  vxZ_upstream_bounds.push_back(-294.5);
  vxZ_downstream_bounds.push_back(-219.5);
  if (binFlag) {
    string line;
    TString dy_type = "";
    Int_t xval = 1;
    ifstream f_bins(binFile);
    if(!f_bins.is_open() ) {
      cout << " " << endl;
      cout << "binFile: " << binFile << " did not open" << endl;
      exit(EXIT_FAILURE); }
    while (!f_bins.eof()) {
      getline(f_bins,line);

      if (line[1] == 'N') {
	if (dy_type == "xN") xval = 1;
	else {
	  dy_type = "xN";
	  xval = 0;
	}			
      }
      else if (line[1] == 'P') {
	if (dy_type == "xPi") xval = 1;
	else {
	  dy_type = "xPi";
	  xval = 0;
	}			
      }
      else if (line[1] == 'F') {
	if (dy_type == "xF") xval = 1;
	else {
	  dy_type = "xF";
	  xval = 0;
	}			
      }
      else if (line[1] == 'T') {
	if (dy_type == "pT") xval = 1;
	else {
	  dy_type = "pT";
	  xval = 0;
	}			
      }
      else if (line[2] == 's') {
	if (dy_type == "M") xval = 1;
	else {
	  dy_type = "M";
	  xval = 0;
	}			
      }
      else if (line[0] == 'r') {
	if (dy_type == "rad") xval = 1;
	else {
	  dy_type = "rad";
	  xval = 0;
	}			
      }
      else if (line[4] == 'u') {
	if (dy_type == "vxZ_upstream") xval = 1;
	else {
	  dy_type = "vxZ_upstream";
	  xval = 0;
	}			
      }
      else if (line[4] == 'd') {
	if (dy_type == "vxZ_downstream") xval = 1;
	else {
	  dy_type = "vxZ_downstream";
	  xval = 0;
	}			
      }

      //Don't read title lines
      if (line[0] == 'x' || line[0] == 'p' || line[1] == 'a' || line[0] == 'v'){
	continue;}

      if (dy_type == "xN"){
	if (xval == 0) xN_bounds.push_back(atof(line.c_str() ) );
	else xN_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "xPi"){
	if (xval == 0) xPi_bounds.push_back(atof(line.c_str() ) );
	else xPi_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "xF"){
	if (xval == 0) xF_bounds.push_back(atof(line.c_str() ) );
	else xF_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "pT"){
	if (xval == 0) pT_bounds.push_back(atof(line.c_str() ) );
	else pT_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "M"){
	if (xval == 0) M_bounds.push_back(atof(line.c_str() ) );
	else M_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "rad"){
	if (xval == 0) rad_bounds.push_back(atof(line.c_str() ) );
	else rad_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "vxZ_upstream"){
	if (xval == 0) vxZ_upstream_bounds.push_back(atof(line.c_str() ) );
	else vxZ_upstream_xval.push_back(atof(line.c_str() ) );
      }
      else if (dy_type == "vxZ_downstream"){
	if (xval == 0) vxZ_downstream_bounds.push_back(atof(line.c_str() ) );
	else vxZ_downstream_xval.push_back(atof(line.c_str() ) );
      }
    }//end file loop
  }//end binFlag
  else {//HM DY binning
    xN_bounds.push_back(0.13);
    xN_bounds.push_back(0.19);
    xPi_bounds.push_back(0.40);
    xPi_bounds.push_back(0.56);
    xF_bounds.push_back(0.21);
    xF_bounds.push_back(0.41);
    pT_bounds.push_back(0.9);
    pT_bounds.push_back(1.4);
    rad_bounds.push_back(0.719511);
    rad_bounds.push_back(1.14036);
    vxZ_upstream_bounds.push_back(-275.021);
    vxZ_upstream_bounds.push_back(-258.531);
    vxZ_downstream_bounds.push_back(-199.956);
    vxZ_downstream_bounds.push_back(-183.598);

    M_bounds.push_back(4.75);
    M_bounds.push_back(5.50);
  }
  xN_bounds.push_back(1.0);
  xPi_bounds.push_back(1.0);
  xF_bounds.push_back(1.0);
  pT_bounds.push_back(5.0);
  rad_bounds.push_back(1.9);
  vxZ_upstream_bounds.push_back(-239.3);
  vxZ_downstream_bounds.push_back(-164.3);
  
  M_bounds.push_back(8.5);
  cout << " " << endl;
  cout << "Warning!!!!!!!" << endl;
  cout << "High Mass" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << " " << endl;

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
  // }}}
  
  //Internal variables and binning
  ////////////////
  // {{{
  Double_t M_proton = 0.938272;
  Int_t nBounds = xN_bounds.size();

  TVectorD tv_xN_bounds(nBounds);
  TVectorD tv_xPi_bounds(nBounds);
  TVectorD tv_xF_bounds(nBounds);
  TVectorD tv_pT_bounds(nBounds);
  TVectorD tv_M_bounds(nBounds);
  TVectorD tv_rad_bounds(nBounds);
  TVectorD tv_vxZ_upstream_bounds(nBounds);
  TVectorD tv_vxZ_downstream_bounds(nBounds);
  for (UInt_t i=0; i<xN_bounds.size(); i++) {
    tv_xN_bounds[i] = xN_bounds.at(i);
    tv_xPi_bounds[i] = xPi_bounds.at(i);
    tv_xF_bounds[i] = xF_bounds.at(i);
    tv_pT_bounds[i] = pT_bounds.at(i);
    tv_M_bounds[i] = M_bounds.at(i);
    tv_rad_bounds[i] = rad_bounds.at(i);
    tv_vxZ_upstream_bounds[i] = vxZ_upstream_bounds.at(i);
    tv_vxZ_downstream_bounds[i] = vxZ_downstream_bounds.at(i);
  }
  Int_t nBins = xN_xval.size();
  if (nBins == 0){
    cout << " " << endl;
    cout << "Error: nBins = 0" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  TVectorD tv_xN_xval(nBins);
  TVectorD tv_xPi_xval(nBins);
  TVectorD tv_xF_xval(nBins);
  TVectorD tv_pT_xval(nBins);
  TVectorD tv_M_xval(nBins);
  TVectorD tv_rad_xval(nBins);
  TVectorD tv_vxZ_upstream_xval(nBins);
  TVectorD tv_vxZ_downstream_xval(nBins);
  for (UInt_t i=0; i<xN_xval.size(); i++) {
    tv_xN_xval[i] = xN_xval.at(i);
    tv_xPi_xval[i] = xPi_xval.at(i);
    tv_xF_xval[i] = xF_xval.at(i);
    tv_pT_xval[i] = pT_xval.at(i);
    tv_M_xval[i] = M_xval.at(i);
    tv_rad_xval[i] = rad_xval.at(i);
    tv_vxZ_upstream_xval[i] = vxZ_upstream_xval.at(i);
    tv_vxZ_downstream_xval[i] = vxZ_downstream_xval.at(i);
  }
  

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

  //Monte Carlo specific
  //MC Positively charged outgoing muon track parameters
  Int_t pid_MCtr1, NMCHits_tr1;
  Double_t Pinv_MCtr1, theta_MCtr1, phi_MCtr1;
  Double_t vMCtr1_X, vMCtr1_Y, vMCtr1_Z, vMCtr1_E;
  //MC Negatively charged outgoing muon track parameters
  Int_t pid_MCtr2, NMCHits_tr2;
  Double_t Pinv_MCtr2, theta_MCtr2, phi_MCtr2;
  Double_t vMCtr2_X, vMCtr2_Y, vMCtr2_Z, vMCtr2_E;
  //MC Beam muon track parameters
  Int_t pid_MCtrIn, NMCHits_trIn, IsBeam_MCtrIn;
  Double_t Pinv_MCtrIn, theta_MCtrIn, phi_MCtrIn;
  Double_t vMCtrIn_X, vMCtrIn_Y, vMCtrIn_Z, vMCtrIn_E;
  //MC DY-variables
  Double_t MC_x_beam, MC_x_target, MC_x_feynman, MC_q_transverse;
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
  T1->SetBranchAddress("event", &event);
  //DY-variables
  T1->SetBranchAddress("x_beam", &x_beam);
  T1->SetBranchAddress("x_target", &x_target);
  T1->SetBranchAddress("x_feynman", &x_feynman);
  T1->SetBranchAddress("q_transverse", &q_transverse);

  //Monte Carlo specific
  //MC Positively charged outgoing muon track parameters
  T1->SetBranchAddress("pid_MCtr1", &pid_MCtr1);
  T1->SetBranchAddress("NMCHits_tr1", &NMCHits_tr1);
  T1->SetBranchAddress("Pinv_MCtr1", &Pinv_MCtr1);
  T1->SetBranchAddress("theta_MCtr1", &theta_MCtr1);
  T1->SetBranchAddress("phi_MCtr1", &phi_MCtr1);
  T1->SetBranchAddress("vMCtr1_X", &vMCtr1_X);
  T1->SetBranchAddress("vMCtr1_Y", &vMCtr1_Y);
  T1->SetBranchAddress("vMCtr1_Z", &vMCtr1_Z);
  T1->SetBranchAddress("vMCtr1_E", &vMCtr1_E);
  //MC Negatively charged outgoing muon track parameters
  T1->SetBranchAddress("pid_MCtr2", &pid_MCtr2);
  T1->SetBranchAddress("NMCHits_tr2", &NMCHits_tr2);
  T1->SetBranchAddress("Pinv_MCtr2", &Pinv_MCtr2);
  T1->SetBranchAddress("theta_MCtr2", &theta_MCtr2);
  T1->SetBranchAddress("phi_MCtr2", &phi_MCtr2);
  T1->SetBranchAddress("vMCtr2_X", &vMCtr2_X);
  T1->SetBranchAddress("vMCtr2_Y", &vMCtr2_Y);
  T1->SetBranchAddress("vMCtr2_Z", &vMCtr2_Z);
  T1->SetBranchAddress("vMCtr2_E", &vMCtr2_E);
  //MC Beam muon track parameters
  T1->SetBranchAddress("pid_MCtrIn", &pid_MCtrIn);
  T1->SetBranchAddress("NMCHits_trIn", &NMCHits_trIn);
  T1->SetBranchAddress("IsBeam_MCtrIn", &IsBeam_MCtrIn);
  T1->SetBranchAddress("Pinv_MCtrIn", &Pinv_MCtrIn);
  T1->SetBranchAddress("theta_MCtrIn", &theta_MCtrIn);
  T1->SetBranchAddress("phi_MCtrIn", &phi_MCtrIn);
  T1->SetBranchAddress("vMCtrIn_X", &vMCtrIn_X);
  T1->SetBranchAddress("vMCtrIn_Y", &vMCtrIn_Y);
  T1->SetBranchAddress("vMCtrIn_Z", &vMCtrIn_Z);
  T1->SetBranchAddress("vMCtrIn_E", &vMCtrIn_E);
  //MC DY-variables
  T1->SetBranchAddress("MC_x_beam", &MC_x_beam);
  T1->SetBranchAddress("MC_x_target", &MC_x_target);
  T1->SetBranchAddress("MC_x_feynman", &MC_x_feynman);
  T1->SetBranchAddress("MC_q_transverse", &MC_q_transverse);
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
  // }}}
  
  //Cut histograms
  ///////////////
  // {{{
  const Int_t nCutHist = 13;//Number of impact cut hist
  const Int_t nCutHist_inNH3 = 2;//Spin dependent cuts
  //const Int_t nMCCuts = 6//src/setup.h
  const Int_t n2D_cutHist = 3;
  TH1D* hCuts = new TH1D("hCuts", "hCuts", 200, 0, 200);
  Int_t cut_bin = 1, cut_space = 10;

  TH1D *hCut_VxZ[nMCCuts];
  TH1D *hCut_MuPTheta[nMCCuts],*hCut_MuPPhi[nMCCuts],*hCut_MuPqP[nMCCuts];
  TH1D *hCut_MuMTheta[nMCCuts],*hCut_MuMPhi[nMCCuts],*hCut_MuMqP[nMCCuts];
  TH1D *hCut_xN[nMCCuts], *hCut_xPi[nMCCuts], *hCut_xF[nMCCuts];
  TH1D *hCut_qT[nMCCuts];
  TH1D *hCut_PhiPhoton[nMCCuts], *hCut_PhiPhoton_gen[nMCCuts];
  TH1D *hCut_PhiS_simple[nMCCuts], *hCut_PhiS_simple_gen[nMCCuts];

  TH2D *hCut_MuP_PxPy[nMCCuts], *hCut_MuM_PxPy[nMCCuts];
  TH2D *hCut_Beam_PxPy[nMCCuts];

  TH1D *hImpactCuts[nCutHist+nCutHist_inNH3][nMCCuts];
  TH2D *h2D_ImpactCuts[n2D_cutHist][nMCCuts];
 
  Int_t ih = 0;
  HistArraySetupMC(hCut_VxZ, hImpactCuts, 500, -500, 100, ih, "VxZ"); ih++;
  HistArraySetupMC(hCut_MuPTheta, hImpactCuts, 100, 0, 0.3, ih, "MuPTheta");
  ih++;
  HistArraySetupMC(hCut_MuPPhi, hImpactCuts, 100, -TMath::Pi(), TMath::Pi(),
		   ih, "MuPPhi"); ih++;
  HistArraySetupMC(hCut_MuPqP, hImpactCuts, 200, 0, 200, ih, "MuPqP"); ih++;
  HistArraySetupMC(hCut_MuMTheta, hImpactCuts, 100, 0, 0.3, ih, "MuMTheta");
  ih++;
  HistArraySetupMC(hCut_MuMPhi, hImpactCuts, 100, -TMath::Pi(), TMath::Pi(),
		   ih, "MuMPhi"); ih++;
  HistArraySetupMC(hCut_MuMqP, hImpactCuts, 200, -200, 0, ih, "MuMqP"); ih++;
  HistArraySetupMC(hCut_xN, hImpactCuts, 100, 0, 1, ih, "xN"); ih++;
  HistArraySetupMC(hCut_xPi, hImpactCuts, 100, 0, 1, ih, "xPi"); ih++;
  HistArraySetupMC(hCut_xF, hImpactCuts, 200, -1, 1, ih, "xF"); ih++;
  HistArraySetupMC(hCut_qT, hImpactCuts, 200, 0, 5, ih, "qT"); ih++;
  HistArraySetupMC(hCut_PhiPhoton, hImpactCuts, 200, -TMath::Pi(), TMath::Pi(),
		   ih,"PhiPhoton");ih++;
  HistArraySetupMC(hCut_PhiPhoton_gen, hImpactCuts,
		   200, -TMath::Pi(), TMath::Pi(),
		   ih,"PhiPhoton_gen");ih++;
  HistArraySetupMC(hCut_PhiS_simple, hImpactCuts, 200, -TMath::Pi(),TMath::Pi(),
		   ih, "PhiS_simple"); ih++;
  HistArraySetupMC(hCut_PhiS_simple_gen, hImpactCuts,
		   200, -TMath::Pi(),TMath::Pi(),
		   ih, "PhiS_simple_gen"); ih++;

  ih=0;
  Hist2D_ArraySetupMC(hCut_MuP_PxPy, h2D_ImpactCuts, 100, -5, 5, 100, -5,5,ih,
		    "MuP_PxPy"); ih++;
  Hist2D_ArraySetupMC(hCut_MuM_PxPy, h2D_ImpactCuts, 100, -5, 5, 100, -5,5,ih,
		    "MuM_PxPy"); ih++;
  Hist2D_ArraySetupMC(hCut_Beam_PxPy, h2D_ImpactCuts, 100, -1, 1, 100,-1,1,ih,
		    "Beam_PxPy"); ih++;
  // }}}

  //pT_Weighted tree
  ///////////////
  // {{{
  TTree *tree = new TTree("pT_Weighted", "pT_Weighted");
  Double_t PhiS, PhiS_simple, Phi_CS, Theta_CS, rapidity;
  Double_t Gen_PhiS, Gen_PhiS_simple, Gen_Phi_CS, Gen_Theta_CS, Gen_rapidity;
  Double_t gen_vPhoton_X, gen_vPhoton_Y, gen_vPhoton_Z, gen_vPhoton_E;
  Int_t targetPosition;
  Double_t Spin;
  tree->Branch("PhiS", &PhiS, "PhiS/D");
  tree->Branch("PhiS_simple", &PhiS_simple, "PhiS_simple/D");
  tree->Branch("Phi_CS", &Phi_CS, "Phi_CS/D");
  tree->Branch("Theta_CS", &Theta_CS, "Theta_CS/D");
  tree->Branch("rapidity", &rapidity, "rapidity/D");
  tree->Branch("Gen_PhiS", &Gen_PhiS, "Gen_PhiS/D");
  tree->Branch("Gen_PhiS_simple", &Gen_PhiS_simple, "Gen_PhiS_simple/D");
  tree->Branch("Gen_Phi_CS", &Gen_Phi_CS, "Gen_Phi_CS/D");
  tree->Branch("Gen_Theta_CS", &Gen_Theta_CS, "Gen_Theta_CS/D");
  tree->Branch("Gen_rapidity", &Gen_rapidity, "Gen_rapidity/D");
  tree->Branch("trigMask", &trigMask, "trigMask/I");
  tree->Branch("vPhoton_X", &vPhoton_X, "vPhoton_X/D");
  tree->Branch("vPhoton_Y", &vPhoton_Y, "vPhoton_Y/D");
  tree->Branch("vPhoton_Z", &vPhoton_Z, "vPhoton_Z/D");
  tree->Branch("vPhoton_E", &vPhoton_E, "vPhoton_E/D");
  tree->Branch("gen_vPhoton_X", &gen_vPhoton_X, "gen_vPhoton_X/D");
  tree->Branch("gen_vPhoton_Y", &gen_vPhoton_Y, "gen_vPhoton_Y/D");
  tree->Branch("gen_vPhoton_Z", &gen_vPhoton_Z, "gen_vPhoton_Z/D");
  tree->Branch("gen_vPhoton_E", &gen_vPhoton_E, "gen_vPhoton_E/D");
  tree->Branch("x_beam", &x_beam, "x_beam/D");
  tree->Branch("x_target", &x_target, "x_target/D");
  tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
  tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
  tree->Branch("Mmumu", &vDiMuon_invM, "Mmumu/D");
  tree->Branch("targetPosition", &targetPosition, "targetPosition/I");
  tree->Branch("Spin", &Spin, "Spin/D");
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
  tree->Branch("vx_zVar", &vx_zVar, "vx_zVar/D");
  tree->Branch("vx_xVar", &vx_xVar, "vx_xVar/D");
  tree->Branch("vx_yVar", &vx_yVar, "vx_yVar/D");
  tree->Branch("vOpenAngle", &vOpenAngle, "vOpenAngle/D");
  tree->Branch("theta_MCtr1", &theta_MCtr1, "theta_MCtr1/D");
  tree->Branch("phi_MCtr1", &phi_MCtr1, "phi_MCtr1/D");
  tree->Branch("Pinv_MCtr1", &Pinv_MCtr1, "Pinv_MCtr1/D");
  tree->Branch("theta_MCtr2", &theta_MCtr2, "theta_MCtr2/D");
  tree->Branch("phi_MCtr2", &phi_MCtr2, "phi_MCtr2/D");
  tree->Branch("Pinv_MCtr2", &Pinv_MCtr2, "Pinv_MCtr2/D");
  tree->Branch("theta_MCtrIn", &theta_MCtrIn, "theta_MCtrIn/D");
  tree->Branch("phi_MCtrIn", &phi_MCtrIn, "phi_MCtrIn/D");
  tree->Branch("Pinv_MCtrIn", &Pinv_MCtrIn, "Pinv_MCtrIn/D");
  tree->Branch("MC_x_beam", &MC_x_beam, "MC_x_beam/D");
  tree->Branch("MC_x_target", &MC_x_target, "MC_x_target/D");
  tree->Branch("MC_x_feynman", &MC_x_feynman, "MC_x_feynman/D");
  tree->Branch("MC_q_transverse", &MC_q_transverse, "MC_q_transverse/D");
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
  // }}}

  //Cuts (Turn on/off)
  // {{{
  Bool_t physKinematics=true, qTcut=true, vxZ_NH3=true, rad_NH3=true;
  Bool_t gen_trIn_Z=true, gen_tr1_tr2=true;
  // }}}

  Int_t tree_entries = T1->GetEntries();
  //Int_t tree_entries = 10000; cout << "Debugging" << endl;//Debug
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  Bool_t first = true;
  for (Int_t i=0; i<tree_entries; i++){
    T1->GetEntry(i, 0);

    //Settings
    if (first || i==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!" << endl;
      cout << "Physical kinematics	=    " << physKinematics << endl;
      cout << "qTcut			=    " << qTcut <<  endl;
      cout << "vxZ_NH3			=    " << vxZ_NH3 << endl;
      cout << "rad_NH3			=    " << rad_NH3 << endl;
      cout << "Positive vMCtrIn_Z	=    " << gen_trIn_Z << endl;
      cout << "Physical gen_tr1_tr2	=    " << gen_tr1_tr2 << endl;
      cout << " " << endl;

      first = false;
    }
    
    //Cuts
    cut_bin = 1;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//All Data

    ///////////////General useful quantities and Compass/TF frame setups
    // {{{
    //Compass frame:
    TLorentzVector lv_p1_Mu(vP1_X, vP1_Y, vP1_Z, vP1_E);//MuPlus
    TLorentzVector lv_p2_Mu(vP2_X, vP2_Y, vP2_Z, vP2_E);//MuMinus
    TLorentzVector lv_beam(beam_X, beam_Y, beam_Z, beam_E);
    TLorentzVector lv_target(0, 0, 0, M_proton);
    TLorentzVector lv_Spin(0, 1.0, 0, 0);//same as Spin_simple now
    TLorentzVector lv_Spin_simple(0, 1.0, 0, 0);
    TLorentzVector lv_diMu = lv_p1_Mu + lv_p2_Mu;

    TLorentzVector lv_Gen_muPlus(vMCtr1_X, vMCtr1_Y, vMCtr1_Z, vMCtr1_E);
    TLorentzVector lv_Gen_muMinus(vMCtr2_X, vMCtr2_Y, vMCtr2_Z, vMCtr2_E);
    TLorentzVector lv_Gen_beam(vMCtrIn_X, vMCtrIn_Y, vMCtrIn_Z,
			       vMCtrIn_E);
    TLorentzVector lv_Gen_virtualPhoton = lv_Gen_muPlus + lv_Gen_muMinus;

    //Target frame:
    TLorentzVector lv_beam_TF(lv_beam);
    TLorentzVector lv_target_TF(lv_target);
    TLorentzVector lv_Spin_TF(lv_Spin);
    TLorentzVector lv_Spin_simple_TF(lv_Spin_simple);
    TLorentzVector lv_muPlus_TF(lv_p1_Mu);
    TLorentzVector lv_muMinus_TF(lv_p2_Mu);
    TLorentzVector lv_virtualPhoton_TF(lv_diMu);

    TLorentzVector lv_Gen_beam_TF(lv_Gen_beam);
    TLorentzVector lv_Gen_target_TF(lv_target);
    TLorentzVector lv_Gen_Spin_TF(lv_Spin);
    TLorentzVector lv_Gen_Spin_simple_TF(lv_Spin_simple);
    TLorentzVector lv_Gen_muPlus_TF(lv_Gen_muPlus);
    TLorentzVector lv_Gen_muMinus_TF(lv_Gen_muMinus);
    TLorentzVector lv_Gen_virtualPhoton_TF(lv_Gen_virtualPhoton);

    Double_t cut_variables[nCutHist] = {vx_z, theta_traj1,
					phi_traj1, qP_traj1, theta_traj2,
					phi_traj2,
					qP_traj2, x_beam, x_target,
					x_feynman, q_transverse, lv_diMu.Phi(),
					lv_Gen_virtualPhoton.Phi() };

    Double_t cut2D_variables[2*n2D_cutHist] = {lv_p1_Mu.X(), lv_p1_Mu.Y(),
					       lv_p2_Mu.X(), lv_p2_Mu.Y(),
					       lv_beam.X(), lv_beam.Y()};

    Bool_t inNH3 = ( (vx_z>-294.5 && vx_z<-239.3) ||
		     (vx_z>-219.5 && vx_z<-164.3) ) ? true : false;
    
    // }}}
    
    //Perform Cuts
    // {{{
    Int_t icut = 0;
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); 
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); icut++;
    
    if (physKinematics && (x_beam < 0.0 || x_beam > 1.0) ) continue;
    if (physKinematics && (x_target < 0.0 || x_target > 1.0) ) continue;
    if (physKinematics && (x_feynman < -1.0 || x_feynman > 1.0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Physical Kinematics
    if (inNH3){
      align_wrt_beam_photon(lv_beam_TF, lv_target_TF, lv_Spin_TF,
			  lv_Spin_simple_TF, lv_virtualPhoton_TF, lv_muPlus_TF,
			  lv_muMinus_TF);
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );

      align_wrt_beam_photon(lv_Gen_beam_TF, lv_Gen_target_TF, lv_Gen_Spin_TF,
			  lv_Gen_Spin_simple_TF, lv_Gen_virtualPhoton_TF,
			  lv_Gen_muPlus_TF, lv_Gen_muMinus_TF);
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );
    }
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if (qTcut && (q_transverse < 0.4 || q_transverse > 5.0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//qT cuts
    if (inNH3){
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );}
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if (vxZ_NH3 && (vx_z < -294.5 || vx_z > -239.3) &&
	(vx_z < -219.5 || vx_z > -164.3)
	) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target z-cut
    if (inNH3){
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );}
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    
    if(rad_NH3 && (TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >=
		   TMath::Power(1.9, 2) ) ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target radial cut
    if (inNH3){
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );}
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (gen_trIn_Z && (vMCtrIn_Z < 0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Positive Pz Pion
    if (inNH3){
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );
    }
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (gen_tr1_tr2 &&
	(vMCtr1_X < -200 || vMCtr2_X < -200 || vMCtrIn_X < -200) ) continue;
    if (gen_tr1_tr2 &&
	(vMCtr1_Y < -200 || vMCtr2_Y < -200 || vMCtrIn_Y < -200) ) continue;
    if (gen_tr1_tr2 &&
	(vMCtr1_Z < -200 || vMCtr2_Z < -200 || vMCtrIn_Z < -200) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//tr1/tr2 physical momentum
    if (inNH3){
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() );
      hCut_PhiS_simple_gen[icut]->Fill(lv_Gen_Spin_simple_TF.Phi() );}
    Fill2D_CutsMC(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist);
    FillCutsMC(hImpactCuts, cut_variables, icut, nCutHist); icut++;
    // }}}
        
    ////All data after cuts
    //////////////
    if (vx_z >= -294.5 && vx_z <= -239.3){//Up stream NH3
      //if (vx_z <= -230){//No target cuts
      Spin = 1.0;
      targetPosition = 0;

    }//Up stream
    else if (vx_z >= -219.5 && vx_z <= -164.3){//Down stream NH3
      //else if (vx_z >= -230){//No target cuts
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
    ///////////////
    // {{{
    //Compass frame:
    gen_vPhoton_X = lv_Gen_virtualPhoton.X();
    gen_vPhoton_Y = lv_Gen_virtualPhoton.Y();
    gen_vPhoton_Z = lv_Gen_virtualPhoton.Z();
    gen_vPhoton_E = lv_Gen_virtualPhoton.E();
    
    //Target frame:
    //Performed after Physical Kinematics cut
    /*align_wrt_beam_photon(lv_beam_TF, lv_target_TF, lv_Spin_TF,
			  lv_Spin_simple_TF, lv_virtualPhoton_TF, lv_muPlus_TF,
			  lv_muMinus_TF);
    align_wrt_beam_photon(lv_Gen_beam_TF, lv_Gen_target_TF, lv_Gen_Spin_TF,
			  lv_Gen_Spin_simple_TF, lv_Gen_virtualPhoton_TF,
			  lv_Gen_muPlus_TF, lv_Gen_muMinus_TF);*/
    
    
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

    TLorentzVector lv_Gen_beam_CS(lv_Gen_beam_TF);
    TLorentzVector lv_Gen_target_CS(lv_Gen_target_TF);
    TLorentzVector lv_Gen_Spin_CS(lv_Gen_Spin_TF);
    TLorentzVector lv_Gen_Spin_simple_CS(lv_Gen_Spin_simple_TF);
    TLorentzVector lv_Gen_muMinus_CS(lv_Gen_muMinus_TF);
    TLorentzVector lv_Gen_muPlus_CS(lv_Gen_muPlus_TF);
    TLorentzVector lv_Gen_virtualPhoton_CS(lv_Gen_virtualPhoton_TF);
    boost_CS(lv_Gen_beam_CS, lv_Gen_target_CS, lv_Gen_Spin_CS,
	     lv_Gen_Spin_simple_CS, lv_Gen_virtualPhoton_CS, lv_Gen_muPlus_CS,
	     lv_Gen_muMinus_CS);
    // }}}

    
    Double_t PhiS_lab = lv_Spin.Phi() - lv_diMu.Phi();
    if(PhiS_lab > TMath::Pi()) PhiS_lab = -2*TMath::Pi() + PhiS_lab;
    else if (PhiS_lab < -1.0*TMath::Pi() ) PhiS_lab = 2*TMath::Pi() + PhiS_lab;
    PhiS = lv_Spin_TF.Phi();//same as PhiS_simple
    PhiS_simple = lv_Spin_simple_TF.Phi();
    Phi_CS = lv_muMinus_CS.Phi();
    Theta_CS = lv_muMinus_CS.Theta();
    rapidity = 0.5*TMath::Log(x_beam/x_target);

    Gen_PhiS = lv_Gen_Spin_TF.Phi();
    Gen_PhiS_simple = lv_Gen_Spin_simple_TF.Phi();
    Gen_Phi_CS = lv_Gen_muMinus_CS.Phi();
    Gen_Theta_CS = lv_Gen_muMinus_CS.Theta();
    Gen_rapidity = 0.5*TMath::Log(MC_x_beam/MC_x_target);

    tree->Fill();
  }//tree entries    

  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  //Cuts histogram
  TString cutNames[nMCCuts] = {"AllData", "xPion,xN,xF", "0.4<qT<5",
			       "TargetZ-cut", "TargetRadius", "Positive_Pz_In",
			       "Physical_Pxyz_tr1_tr2"};
  for (Int_t i=0, j=1; i<nMCCuts; i++, j+=cut_space){
    Int_t bin_index = hCuts->GetXaxis()->FindBin(j);
    hCuts->GetXaxis()->SetBinLabel(bin_index, cutNames[i]);
  }

  //Write Output
  ////////////////
  // {{{
  if (!wflag && !Qflag) cout << "No file output" << endl;
  else{ 
    if (wflag) outFile += "RealData.root";
    TFile *myFile = new TFile(outFile, "RECREATE");
    hCuts->Write();
    tree->Write();

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
    TDirectory *PhiPhoton_CutImpact = myFile->mkdir("PhiPhoton_CutImpact");
    TDirectory *PhiPhoton_gen_CutImpact=myFile->mkdir("PhiPhoton_gen_CutImpact");
    TDirectory *PhiS_simple_CutImpact = myFile->mkdir("PhiS_simple_CutImpact");
    TDirectory *PhiS_simple_gen_CutImpact=myFile->mkdir("PhiS_simple_gen_CutImpact");
    
    TDirectory *MuP_PxPy_CutImpact = myFile->mkdir("MuP_PxPy_CutImpact");
    TDirectory *MuM_PxPy_CutImpact = myFile->mkdir("MuM_PxPy_CutImpact");
    TDirectory *Beam_PxPy_CutImpact = myFile->mkdir("Beam_PxPy_CutImpact");
    for (Int_t i=0; i<nMCCuts; i++) {
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
      
      PhiPhoton_CutImpact->cd();
      hCut_PhiPhoton[i]->Write(cutNames[i]);
      PhiPhoton_gen_CutImpact->cd();
      hCut_PhiPhoton_gen[i]->Write(cutNames[i]);
      PhiS_simple_CutImpact->cd();
      hCut_PhiS_simple[i]->Write(cutNames[i]);
      PhiS_simple_gen_CutImpact->cd();
      hCut_PhiS_simple_gen[i]->Write(cutNames[i]);

      MuP_PxPy_CutImpact->cd();
      hCut_MuP_PxPy[i]->Write(cutNames[i]);
      MuM_PxPy_CutImpact->cd();
      hCut_MuM_PxPy[i]->Write(cutNames[i]);
      Beam_PxPy_CutImpact->cd();
      hCut_Beam_PxPy[i]->Write(cutNames[i]);
    }

    myFile->cd();
    
    tv_xN_bounds.Write("tv_xN_bounds");
    tv_xPi_bounds.Write("tv_xPi_bounds");
    tv_xF_bounds.Write("tv_xF_bounds");
    tv_pT_bounds.Write("tv_pT_bounds");
    tv_M_bounds.Write("tv_M_bounds");

    tv_xN_xval.Write("tv_xN_xval");
    tv_xPi_xval.Write("tv_xPi_xval");
    tv_xF_xval.Write("tv_xF_xval");
    tv_pT_xval.Write("tv_pT_xval");
    tv_M_xval.Write("tv_M_xval");

    
    cout << myFile->GetName() << " was written" << endl;
    myFile->Close();
  }
  // }}}


  theApp.Run();//Needed to make root graphics work on C++
}//main

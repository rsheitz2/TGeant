#include "common.hxx"
#include "functions.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "To be used with Real Data or MC Data" << endl;
    cout << "To determine the bin bounds and average bin values";
    cout << " for the given number of Bins" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-Pperiod] [-nBins] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
	cout << "---Needed Options---" << endl;
    cout << "Option:  -P period         (which period to take bad spills from)"
	 << endl;
    cout << "         (i.e W07, W08...  Can also enter \"WAll\" for ALL periods or" << 
      " \"MC\" for NO periods)" << endl;
	cout << "Option:  -n bins           (How many bins to make)" << endl;
    cout << "" << endl;
	cout << "---Write Option---" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "    Otherwise bin values are output and over written to ";
    cout << "\"binValues.txt\"" << endl;
    cout << "" << endl;
	cout << "---Additional Cut Options---" << endl;
    cout << "Option:  -i minMass (to specify a minimum mass cut)";
	cout << "          (default value is 0.0)" << endl;
    cout << "Option:  -a maxMass (to specify a maximum mass cut)";
	cout << "          (default value is 16.0)" << endl;
	cout << "" << endl;
	cout << "---Additional Options---" << endl;
    cout << "Option:  -u ##		(new UserEvent number, default==420)" << endl;
    cout << "" << endl;
	
    exit(EXIT_FAILURE);
  }
  cout << "" << endl;
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  Int_t uflag=0, Qflag=0, fflag=0, Pflag=0, binflag=0, iflag=0, aflag=0;
  Int_t c;
  TString userNum = "", fname = "", outFile = "", period = "";
  Double_t M_min = 0.0, M_max=16.0;
  Int_t nBins;
  
  while ((c = getopt (argc, argv, "n:u:f:Q:P:i:a:")) != -1) {
    switch (c) {
    case 'n':
      binflag = 1;
      nBins = atoi(optarg);
      break;
    case 'u':
      uflag = 1;
      userNum += optarg;
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
    case 'i':
      iflag = 1;
	  M_min = stof(optarg);
      break;
    case 'a':
      aflag = 1;
	  M_max = stof(optarg);
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
      else if (optopt == 'i')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'a')
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
  if(!binflag){
    cout << "number of bins must be specified: -nBins" << endl;
    exit(EXIT_FAILURE);
  }

  TString userEvent = "UserEvent";
  if (!uflag) {
    userEvent += "420/Particles";
    cout << "Default UserEvent420 used" << endl;
  }
  else userEvent += userNum + "/Particles";
  TChain* T1 = new TChain(userEvent);
  
  TString BadSpillPath = "Src/BadSpills/t3/";//Get badspill information
  map <Long64_t, vector<Long64_t> > BadMap;
  Long64_t BadRun, BadSpill;
  if (!Pflag){
    cout << "Please enter a period for bad spills list" << endl;
    exit(EXIT_FAILURE);
  }
  else if (period == "MC"){}
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

  if (period == "MC"){}
  else if (BadMap.size() == 0){
    cout << "\"MC\" option not specified and Bad spills file did not open" << endl;
    exit(EXIT_FAILURE);    
  }

  //Internal variables and binning
  Double_t M_proton = 0.938272;

  //Averages
  Double_t AvgPolarization=0.0, AvgDilution=0.0, AvgDilution_corrected=0.0;
  Int_t AvgPolarization_count=0, AvgDilution_count=0;

  //Vectors for sorting
  vector<Double_t> sort_xTarg;
  vector<Double_t> sort_xBeam;
  vector<Double_t> sort_xF;
  vector<Double_t> sort_pT;
  vector<Double_t> sort_mass;
  vector<Double_t> sort_rad;
  vector<Double_t> sort_vxZ_upstream;
  vector<Double_t> sort_vxZ_downstream;

  
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
    

  Int_t tree_entries = T1->GetEntries();
  //Int_t tree_entries = 1000;//Debug
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  Bool_t first = true;
  for (Int_t i=0; i<tree_entries; i++){
    T1->GetEntry(i, 0);
    
    //Settings
    if (first || i==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!" << endl;
	  cout << "Normal DY cuts" << endl;

	  if (iflag || aflag) {
		  cout << "Additional Mass cut" << endl;
		  cout << "    Mass range " << M_min << " - " << M_max << endl;
	  }

      first = false;
    }

	//Additional cuts
	if ( (vDiMuon_invM < M_min) || (vDiMuon_invM > M_max) ) continue;

    //Cuts
    cut_bin = 1;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//All Data
    map<Long64_t, vector<Long64_t> >::iterator it = BadMap.find(RunNum);
    if (it != BadMap.end() ){
      if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(), -310)
	  != BadMap[it->first].end() ) continue;
      else if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(),
			 SpillNum) != BadMap[it->first].end() ) continue;
    }
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Bad spills
  
    TLorentzVector lv_p1_Mu(vP1_X, vP1_Y, vP1_Z, vP1_E);
    TLorentzVector lv_p2_Mu(vP2_X, vP2_Y, vP2_Z, vP2_E);
    TLorentzVector lv_diMu = lv_p1_Mu + lv_p2_Mu;
    TLorentzVector lv_target_1 (0, 0, 0, M_proton);

    if (x_beam < 0.0 || x_beam > 1.0) continue;
    if (x_target < 0.0 || x_target > 1.0) continue;
    if (x_feynman < -1.0 || x_feynman > 1.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Physical Kinematics
      
    if (q_transverse < 0.4 || q_transverse > 5.0) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//qT cuts
      
    if ( (vx_z < -294.5 || vx_z > -239.3) && (vx_z < -219.5 || vx_z > -164.3)
	 ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target z-cut
      
    if(TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >= TMath::Power(1.9, 2)
       ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target radial cut

    
    ////All data after cuts
    //////////////
    Double_t radius = TMath::Sqrt(vx_x*vx_x+ vx_y*vx_y);


    if (vx_z >= -294.5 && vx_z <= -239.3){//Up stream NH3
      if (period != "MC") {
	AvgDilution += TMath::Abs(dilutionFactor);
	AvgDilution_corrected += 0.95*TMath::Abs(dilutionFactor);
	AvgDilution_count++;

	AvgPolarization += TMath::Abs(Polarization);
	AvgPolarization_count++;
      }

      sort_vxZ_upstream.push_back(vx_z);
    }//Up stream
    else if (vx_z >= -219.5 && vx_z <= -164.3){//Down stream NH3
      if (period != "MC") {
	AvgDilution += TMath::Abs(dilutionFactor);
	AvgDilution_corrected += 0.91*TMath::Abs(dilutionFactor);
	AvgDilution_count++;

	AvgPolarization += TMath::Abs(Polarization);
	AvgPolarization_count++;
      }

      sort_vxZ_downstream.push_back(vx_z);
    }//Down Stream

    //Both targets
    sort_xTarg.push_back(x_target);
    sort_xBeam.push_back(x_beam);
    sort_xF.push_back(x_feynman);
    sort_pT.push_back(q_transverse);
    sort_mass.push_back(vDiMuon_invM);
    sort_rad.push_back(radius);

  }//tree entries
  

  //Print bin boundaries and average values
  if (sort_xF.size() != 0){
    if (!Qflag) outFile += "binValues.txt";
    ofstream binValueFile(outFile, std::ofstream::trunc);
    
    std::sort(sort_xTarg.begin(), sort_xTarg.end() );
    PrintBin(binValueFile, sort_xTarg, nBins, "xN");
    std::sort(sort_xBeam.begin(), sort_xBeam.end() );
    PrintBin(binValueFile, sort_xBeam, nBins, "xPi");
    std::sort(sort_xF.begin(), sort_xF.end() );
    PrintBin(binValueFile, sort_xF, nBins, "xF");
    std::sort(sort_pT.begin(), sort_pT.end() );
    PrintBin(binValueFile, sort_pT, nBins, "pT");
    std::sort(sort_mass.begin(), sort_mass.end() );
    PrintBin(binValueFile, sort_mass, nBins, "mass");
    std::sort(sort_rad.begin(), sort_rad.end() );
    PrintBin(binValueFile, sort_rad, nBins, "rad");
    std::sort(sort_vxZ_upstream.begin(), sort_vxZ_upstream.end() );
    PrintBin(binValueFile, sort_vxZ_upstream, nBins, "vxZ_upstream");
    std::sort(sort_vxZ_downstream.begin(), sort_vxZ_downstream.end() );
    PrintBin(binValueFile, sort_vxZ_downstream, nBins, "vxZ_downstream");
    
    binValueFile.close();

    cout << "File: " << outFile << " was written" << endl;
  }

  if (period != "MC") {
    cout << " " << endl;
    cout << "Average dilution factor " << AvgDilution_corrected/AvgDilution_count;
    cout << endl;
    cout << "Average polarization " << AvgPolarization/AvgPolarization_count;
    cout << endl;;
  }
  
  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  //Cuts histogram
  Int_t nRealCuts = 6;
  TString cutNames[nRealCuts] = {"AllData", "GoodSpills", "xPion,xN,xF",
				 "0.4<qT<5",
				 "TargetZ-cut", "TargetRadius"};
  for (Int_t i=0, j=1; i<nRealCuts; i++, j+=cut_space){
    Int_t bin_index = hCuts->GetXaxis()->FindBin(j);
    hCuts->GetXaxis()->SetBinLabel(bin_index, cutNames[i]);
  }

  theApp.Run();//Needed to make root graphics work on C++
}//main

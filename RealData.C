#include "common.hxx"
#include "functions.h"
#include "setup.h"
#include <bitset>

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "To be used with Real Data" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-Pperiod] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
	cout << "---Needed Options---" << endl;
    cout << "Option:  -P period         (which period to take bad spills from)"
	 << endl;
    cout << "(i.e W07, W08...  Can also enter \"WAll\" for all periods)" <<endl;
    cout << "" << endl;
	cout << "---Write Options---" << endl;
    cout << "Option:  -w		(write output to file)" << endl;
    cout << "        default output file is named \"Output.root\"" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "" << endl;
	cout << "---Binning Options---" << endl;
    cout << "Option:  -b textfile with binning information	";
    cout << "(textfile should be made from Macro/Binning/avgBinBounds.C)"<<endl;
    cout << "Option:  -M (\"HM\", \"JPsi\", \"AMDY\") to specify which mass ";
    cout << "range to use for \"binning information\" " << endl;
    cout << "   (default mass range is high mass)" << endl;
    cout << "" << endl;
	cout << "---Additional Cut Options---" << endl;
    cout << "Option:  -i minMass (to specify a minimum mass cut)";
	cout << "          (default value is 0.0)" << endl;
    cout << "Option:  -a maxMass (to specify a maximum mass cut)";
	cout << "          (default value is 16.0)" << endl;
	cout << "" << endl;
	cout << "---Additional Options---" << endl;
    cout << "Option:  -u ##		(new UserEvent number, default==420)"
	 << endl;
    exit(EXIT_FAILURE);
  }
  cout << "" << endl;
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  ///////////////
  // {{{
  Int_t uflag=0, wflag=0, Qflag=0, fflag=0, Pflag=0, binFlag=0, iflag=0, aflag=0;
  Int_t Mflag=0;
  Int_t c;
  TString userNum = "", fname = "", outFile = "", period = "", binFile = "";
  TString massRange="";
  Double_t M_min=0.0, M_max=16.0;

  while ((c = getopt (argc, argv, "wM:b:u:f:Q:P:")) != -1) {
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
    case 'M':
      Mflag = 1;
	  massRange += optarg;
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
      else if (optopt == 'b')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'i')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'a')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'M')
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
  rad_bounds.push_back(0.0);
  vxZ_upstream_bounds.push_back(-294.5);
  vxZ_downstream_bounds.push_back(-219.5);

  if (!Mflag || massRange=="HM") M_bounds.push_back(4.3);//High mass
  else if (massRange=="JPsi")M_bounds.push_back(2.5);//JPsi mass
  else if (massRange=="AMDY")M_bounds.push_back(0.0);//All Mass DY
  else {
    cout << "Invalid mass range specified" << endl;
    exit(EXIT_FAILURE);
  }

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
	continue;
      }
      
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
  else {
    cout << " " << endl;
    cout << "No bin flag file given" << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }
  xN_bounds.push_back(1.0);
  xPi_bounds.push_back(1.0);
  xF_bounds.push_back(1.0);
  pT_bounds.push_back(5.0);
  rad_bounds.push_back(1.9);
  vxZ_upstream_bounds.push_back(-239.3);
  vxZ_downstream_bounds.push_back(-164.3);
  if(xN_xval.size()==0 || xPi_xval.size()==0 || xF_xval.size()==0 ||
     pT_xval.size()==0 || M_xval.size()==0 || rad_xval.size()==0 ||
     vxZ_upstream_xval.size()==0 || vxZ_downstream_xval.size()==0){
    cout << "Error:" << endl;
    cout << "Modern xval values not specifed in " << binFile << endl;
    cout << " " << endl;
    exit(EXIT_FAILURE);
  }

  if (!Mflag || massRange=="HM") M_bounds.push_back(8.5);//High mass
  else if (massRange=="JPsi")M_bounds.push_back(4.3);//JPsi mass
  else if (massRange=="AMDY")M_bounds.push_back(16.0);//All Mass DY
  cout << " " << endl;
  cout << "Mass range set to:" << endl;
  (!Mflag) ? cout << "HM" << endl : cout << massRange << endl;
  cout << " " << endl;

  if (M_bounds.at(0) > M_bounds.at(1) || M_bounds.back() < M_xval.front() ){
	  cout << "Mass Range not setup correct" << endl;
    exit(EXIT_FAILURE);
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
  // }}}

  //Internal variables and binning
  ///////////////
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


  //Averages
  Double_t AvgPolarization=0.0, AvgDilution=0.0, AvgDilution_corrected=0.0;
  Int_t AvgPolarization_count=0, AvgDilution_count=0;

  vector<Double_t> AvgPol_xN(nBins, 0.0), AvgDil_xN(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count(nBins, 0), AvgDil_xN_count(nBins, 0);
  vector<Double_t> AvgPol_xN_UpStream(nBins, 0.0), AvgDil_xN_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream(nBins, 0), AvgDil_xN_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream(nBins, 0.0), AvgDil_xN_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream(nBins, 0), AvgDil_xN_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xN_UpStream_Up(nBins, 0.0), AvgDil_xN_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream_Up(nBins, 0), AvgDil_xN_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream_Up(nBins, 0.0), AvgDil_xN_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream_Up(nBins, 0), AvgDil_xN_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xN_UpStream_Down(nBins, 0.0), AvgDil_xN_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_UpStream_Down(nBins, 0), AvgDil_xN_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xN_DownStream_Down(nBins, 0.0), AvgDil_xN_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xN_count_DownStream_Down(nBins, 0), AvgDil_xN_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_xPi(nBins, 0.0), AvgDil_xPi(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count(nBins, 0), AvgDil_xPi_count(nBins, 0);
  vector<Double_t> AvgPol_xPi_UpStream(nBins, 0.0), AvgDil_xPi_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream(nBins, 0), AvgDil_xPi_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream(nBins, 0.0), AvgDil_xPi_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream(nBins, 0), AvgDil_xPi_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xPi_UpStream_Up(nBins, 0.0), AvgDil_xPi_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream_Up(nBins, 0), AvgDil_xPi_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream_Up(nBins, 0.0), AvgDil_xPi_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream_Up(nBins, 0), AvgDil_xPi_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xPi_UpStream_Down(nBins, 0.0), AvgDil_xPi_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_UpStream_Down(nBins, 0), AvgDil_xPi_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xPi_DownStream_Down(nBins, 0.0), AvgDil_xPi_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xPi_count_DownStream_Down(nBins, 0), AvgDil_xPi_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_xF(nBins, 0.0), AvgDil_xF(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count(nBins, 0), AvgDil_xF_count(nBins, 0);
  vector<Double_t> AvgPol_xF_UpStream(nBins, 0.0), AvgDil_xF_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream(nBins, 0), AvgDil_xF_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream(nBins, 0.0), AvgDil_xF_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream(nBins, 0), AvgDil_xF_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_xF_UpStream_Up(nBins, 0.0), AvgDil_xF_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream_Up(nBins, 0), AvgDil_xF_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream_Up(nBins, 0.0), AvgDil_xF_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream_Up(nBins, 0), AvgDil_xF_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_xF_UpStream_Down(nBins, 0.0), AvgDil_xF_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_UpStream_Down(nBins, 0), AvgDil_xF_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_xF_DownStream_Down(nBins, 0.0), AvgDil_xF_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_xF_count_DownStream_Down(nBins, 0), AvgDil_xF_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_pT(nBins, 0.0), AvgDil_pT(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count(nBins, 0), AvgDil_pT_count(nBins, 0);
  vector<Double_t> AvgPol_pT_UpStream(nBins, 0.0), AvgDil_pT_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream(nBins, 0), AvgDil_pT_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream(nBins, 0.0), AvgDil_pT_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream(nBins, 0), AvgDil_pT_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_pT_UpStream_Up(nBins, 0.0), AvgDil_pT_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream_Up(nBins, 0), AvgDil_pT_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream_Up(nBins, 0.0), AvgDil_pT_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream_Up(nBins, 0), AvgDil_pT_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_pT_UpStream_Down(nBins, 0.0), AvgDil_pT_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_UpStream_Down(nBins, 0), AvgDil_pT_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_pT_DownStream_Down(nBins, 0.0), AvgDil_pT_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_pT_count_DownStream_Down(nBins, 0), AvgDil_pT_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_M(nBins, 0.0), AvgDil_M(nBins, 0.0);
  vector<Int_t> AvgPol_M_count(nBins, 0), AvgDil_M_count(nBins, 0);
  vector<Double_t> AvgPol_M_UpStream(nBins, 0.0), AvgDil_M_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream(nBins, 0), AvgDil_M_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream(nBins, 0.0), AvgDil_M_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream(nBins, 0), AvgDil_M_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_M_UpStream_Up(nBins, 0.0), AvgDil_M_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream_Up(nBins, 0), AvgDil_M_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream_Up(nBins, 0.0), AvgDil_M_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream_Up(nBins, 0), AvgDil_M_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_M_UpStream_Down(nBins, 0.0), AvgDil_M_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_UpStream_Down(nBins, 0), AvgDil_M_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_M_DownStream_Down(nBins, 0.0), AvgDil_M_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_M_count_DownStream_Down(nBins, 0), AvgDil_M_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_rad(nBins, 0.0), AvgDil_rad(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count(nBins, 0), AvgDil_rad_count(nBins, 0);
  vector<Double_t> AvgPol_rad_UpStream(nBins, 0.0), AvgDil_rad_UpStream(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream(nBins, 0), AvgDil_rad_count_UpStream(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream(nBins, 0.0), AvgDil_rad_DownStream(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream(nBins, 0), AvgDil_rad_count_DownStream(nBins, 0);
  //Polarization by target
  vector<Double_t> AvgPol_rad_UpStream_Up(nBins, 0.0), AvgDil_rad_UpStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream_Up(nBins, 0), AvgDil_rad_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream_Up(nBins, 0.0), AvgDil_rad_DownStream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream_Up(nBins, 0), AvgDil_rad_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_rad_UpStream_Down(nBins, 0.0), AvgDil_rad_UpStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_UpStream_Down(nBins, 0), AvgDil_rad_count_UpStream_Down(nBins, 0);
  vector<Double_t> AvgPol_rad_DownStream_Down(nBins, 0.0), AvgDil_rad_DownStream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_rad_count_DownStream_Down(nBins, 0), AvgDil_rad_count_DownStream_Down(nBins, 0);

  vector<Double_t> AvgPol_vxZ_upstream(nBins, 0.0), AvgDil_vxZ_upstream(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream(nBins, 0), AvgDil_vxZ_count_UpStream(nBins, 0);

  //Polarization by target
  vector<Double_t> AvgPol_vxZ_upstream_Up(nBins, 0.0), AvgDil_vxZ_upstream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream_Up(nBins, 0), AvgDil_vxZ_count_UpStream_Up(nBins, 0);
  vector<Double_t> AvgPol_vxZ_upstream_Down(nBins, 0.0), AvgDil_vxZ_upstream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_UpStream_Down(nBins, 0), AvgDil_vxZ_count_UpStream_Down(nBins, 0);

  vector<Double_t> AvgPol_vxZ_downstream(nBins, 0.0), AvgDil_vxZ_downstream(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream(nBins, 0), AvgDil_vxZ_count_DownStream(nBins, 0);

  //Polarization by target
  vector<Double_t> AvgPol_vxZ_downstream_Up(nBins, 0.0), AvgDil_vxZ_downstream_Up(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream_Up(nBins, 0), AvgDil_vxZ_count_DownStream_Up(nBins, 0);
  vector<Double_t> AvgPol_vxZ_downstream_Down(nBins, 0.0), AvgDil_vxZ_downstream_Down(nBins, 0.0);
  vector<Int_t> AvgPol_vxZ_count_DownStream_Down(nBins, 0), AvgDil_vxZ_count_DownStream_Down(nBins, 0);

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
  // }}}

  //Cut histograms
  ///////////////
  // {{{

  const Int_t nCutHist = 12;//Number of impact cut hist
  const Int_t nCutHist_inNH3 = 2;//Spin dependent cuts
  //const Int_t nRealCuts = 7;//Number of cuts made (setup.h)
  const Int_t n2D_cutHist = 3;
  TH1D* hCuts = new TH1D("hCuts", "hCuts", 200, 0, 200);
  Int_t cut_bin = 1, cut_space = 10;

  TH1D *hCut_VxZ[nRealCuts];
  TH1D *hCut_MuPTheta[nRealCuts],*hCut_MuPPhi[nRealCuts],*hCut_MuPqP[nRealCuts];
  TH1D *hCut_MuMTheta[nRealCuts],*hCut_MuMPhi[nRealCuts],*hCut_MuMqP[nRealCuts];
  TH1D *hCut_xN[nRealCuts], *hCut_xPi[nRealCuts], *hCut_xF[nRealCuts];
  TH1D *hCut_qT[nRealCuts], *hCut_PhiPhoton[nRealCuts];
  TH1D *hCut_PhiS[nRealCuts], *hCut_PhiS_simple[nRealCuts];

  TH2D *hCut_MuP_PxPy[nRealCuts], *hCut_MuM_PxPy[nRealCuts];
  TH2D *hCut_Beam_PxPy[nRealCuts];

  TH1D *hImpactCuts[nCutHist+nCutHist_inNH3][nRealCuts];
  TH2D *h2D_ImpactCuts[n2D_cutHist][nRealCuts];
  
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
  HistArraySetupMC(hCut_PhiPhoton, hImpactCuts, 200, -TMath::Pi(), TMath::Pi(),
		   ih,"PhiPhoton");ih++;
  HistArraySetupMC(hCut_PhiS, hImpactCuts, 200, -TMath::Pi(), TMath::Pi(),
		   ih,"PhiS");ih++;
  HistArraySetupMC(hCut_PhiS_simple, hImpactCuts,
		   200, -TMath::Pi(), TMath::Pi(),
		   ih,"PhiS_simple");ih++;

  ih=0;
  Hist2D_ArraySetupReal(hCut_MuP_PxPy, h2D_ImpactCuts, 100, -5, 5, 100, -5,5,ih,
		    "MuP_PxPy"); ih++;
  Hist2D_ArraySetupReal(hCut_MuM_PxPy, h2D_ImpactCuts, 100, -5, 5, 100, -5,5,ih,
		    "MuM_PxPy"); ih++;
  Hist2D_ArraySetupReal(hCut_Beam_PxPy, h2D_ImpactCuts, 100, -5, 5, 100,-5,5,ih,
		    "Beam_PxPy"); ih++;
  // }}}

  //pT_Weighted tree
  ///////////////
  // {{{
  TTree *tree = new TTree("pT_Weighted", "pT_Weighted");
  Double_t PhiS, PhiS_simple, Phi_CS, Theta_CS, rapidity;
  Int_t targetPosition;
  Double_t Spin[7];
  tree->Branch("PhiS", &PhiS, "PhiS/D");
  tree->Branch("PhiS_simple", &PhiS_simple, "PhiS_simple/D");
  tree->Branch("Phi_CS", &Phi_CS, "Phi_CS/D");
  tree->Branch("Theta_CS", &Theta_CS, "Theta_CS/D");
  tree->Branch("rapidity", &rapidity, "rapidity/D");
  tree->Branch("MasterTrigMask", &MasterTrigMask, "MasterTrigMask/I");
  tree->Branch("trigMask", &trigMask, "trigMask/I");
  tree->Branch("vPhoton_X", &vPhoton_X, "vPhoton_X/D");
  tree->Branch("vPhoton_Y", &vPhoton_Y, "vPhoton_Y/D");
  tree->Branch("vPhoton_Z", &vPhoton_Z, "vPhoton_Z/D");
  tree->Branch("vPhoton_E", &vPhoton_E, "vPhoton_E/D");
  tree->Branch("x_beam", &x_beam, "x_beam/D");
  tree->Branch("x_target", &x_target, "x_target/D");
  tree->Branch("x_feynman", &x_feynman, "x_feynman/D");
  tree->Branch("q_transverse", &q_transverse, "q_transverse/D");
  tree->Branch("Mmumu", &vDiMuon_invM, "Mmumu/D");
  for (Int_t i=0; i<7; i++) {
    tree->Branch(Form("Spin_%i", i), &Spin[i], Form("Spin_%i/D", i));
  }
  tree->Branch("dilutionFactor", &dilutionFactor, "dilutionFactor/D");
  tree->Branch("Polarization", &Polarization, "Polarization/D");
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
  tree->Branch("vx_zVar", &vx_zVar, "vx_zVar/D");
  tree->Branch("vx_xVar", &vx_xVar, "vx_xVar/D");
  tree->Branch("vx_yVar", &vx_yVar, "vx_yVar/D");
  tree->Branch("vOpenAngle", &vOpenAngle, "vOpenAngle/D");
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
  Bool_t badSpillCut=true, physKinematics=true, qTcut=true, dilCut=true;
  Bool_t vxZ_NH3=true, rad_NH3=true;
  // }}}

  Int_t tree_entries = T1->GetEntries();
  //Int_t tree_entries = 1000;//Debug
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  Bool_t first = true;
  for (Int_t ev=0; ev<tree_entries; ev++){
    T1->GetEntry(ev, 0);

    //Settings
    if (first || ev==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!" << endl;
      cout << "Bad spills		=    " << badSpillCut << endl;
      cout << "Physical kinematics	=    " << physKinematics << endl;
      cout << "qTcut			=    " << qTcut <<  endl;
      cout << "Nonzero dilution factor  =    " << dilCut << endl;
      cout << "vxZ_NH3			=    " << vxZ_NH3 << endl;
      cout << "rad_NH3			=    " << rad_NH3 << endl;
      cout << " " << endl;
	  if (iflag || aflag) {
		  cout << "Additional Mass cut" << endl;
		  cout << "    Mass range " << M_min << " - " << M_max << endl;
	  }
      cout << " " << endl;

      first = false;
    }

	//Additional cuts
	if ( (vDiMuon_invM < M_min) || (vDiMuon_invM > M_max) ) continue;

    //Perform Cuts
    // {{{
    //Setup Vectors in different coordinate systems
    //Compass frame:
    TLorentzVector lv_beam(beam_X, beam_Y, beam_Z, beam_E);
    TLorentzVector lv_target(0, 0, 0, M_proton);
    TLorentzVector lv_Spin(0, Spin[0], 0, 0);
    TLorentzVector lv_Spin_simple(0, 1, 0, 0);
    TLorentzVector lv_p1_Mu(vP1_X, vP1_Y, vP1_Z, vP1_E);//MuPlus
    TLorentzVector lv_p2_Mu(vP2_X, vP2_Y, vP2_Z, vP2_E);//MuMinus
    TLorentzVector lv_diMu(vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E);

    //Target frame:
    TLorentzVector lv_beam_TF(lv_beam);
    TLorentzVector lv_target_TF(lv_target);
    TLorentzVector lv_Spin_TF(lv_Spin);
    TLorentzVector lv_Spin_simple_TF(lv_Spin_simple);
    TLorentzVector lv_muPlus_TF(lv_p1_Mu);
    TLorentzVector lv_muMinus_TF(lv_p2_Mu);
    TLorentzVector lv_virtualPhoton_TF(lv_diMu);
    
    cut_bin = 1;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//All Data

    Double_t cut_variables[nCutHist] = {vx_z, theta_traj1,
					phi_traj1, qP_traj1, theta_traj2,
					phi_traj2,
					qP_traj2, x_beam, x_target,
					x_feynman, q_transverse, lv_diMu.Phi()};

    Double_t cut2D_variables[2*n2D_cutHist] = {lv_p1_Mu.X(), lv_p1_Mu.Y(),
					       lv_p2_Mu.X(), lv_p2_Mu.Y(),
					       lv_beam.X(), lv_beam.Y()};

    Bool_t inNH3 = ( (vx_z>-294.5 && vx_z<-239.3) ||
		     (vx_z>-219.5 && vx_z<-164.3) ) ? true : false;
    
    Int_t icut = 0;
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist);
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    if (inNH3){
      align_wrt_beam_photon(lv_beam_TF, lv_target_TF, lv_Spin_TF,
			  lv_Spin_simple_TF, lv_virtualPhoton_TF, lv_muPlus_TF,
			  lv_muMinus_TF);
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); } icut++;

    map<Long64_t, vector<Long64_t> >::iterator it = BadMap.find(RunNum);
    if (badSpillCut && it != BadMap.end() ){
      if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(), -310)
	  != BadMap[it->first].end() ) continue;
      else if (std::find(BadMap[it->first].begin(), BadMap[it->first].end(),
			 SpillNum) != BadMap[it->first].end() ) continue;
    }
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Bad spills
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (physKinematics && (x_beam < 0.0 || x_beam > 1.0) ) continue;
    if (physKinematics && (x_target < 0.0 || x_target > 1.0) ) continue;
    if (physKinematics && (x_feynman < -1.0 || x_feynman > 1.0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Physical Kinematics
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (qTcut && (q_transverse < 0.4 || q_transverse > 5.0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//qT cuts
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() );
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (dilCut && (dilutionFactor == 0.0) ) continue;
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//dilution Factor
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if (vxZ_NH3 && (vx_z < -294.5 || vx_z > -239.3) &&
	(vx_z < -219.5 || vx_z > -164.3)
	 ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target z-cut
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    if(rad_NH3 &&
       (TMath::Power(vx_x, 2) + TMath::Power(vx_y, 2) >= TMath::Power(1.9, 2) )
       ) continue;//NH3 targets
    hCuts->Fill(cut_bin-1); cut_bin += cut_space;//Target radial cut
    if (inNH3) {
      hCut_PhiS[icut]->Fill(lv_Spin_TF.Phi() ); 
      hCut_PhiS_simple[icut]->Fill(lv_Spin_simple_TF.Phi() ); }
    Fill2D_CutsReal(h2D_ImpactCuts, cut2D_variables, icut, n2D_cutHist); 
    FillCutsReal(hImpactCuts, cut_variables, icut, nCutHist); icut++;

    // }}}

    ////All data after cuts
    //////////////
    ///////////////General useful quantities
    Double_t radius = TMath::Sqrt(vx_x*vx_x+vx_y*vx_y);

    if (vx_z >= -294.5 && vx_z <= -239.3){//Up stream NH3
      Spin[0] = -1.0*avgUpStream/(TMath::Abs(avgUpStream) );
      Spin[1] = avgUpStream;
      Spin[2] = upStreamCoil1;
      Spin[3] = upStreamCoil2;
      Spin[4] = upStreamCoil3;
      Spin[5] = upStreamCoil4;
      Spin[6] = upStreamCoil5;
      targetPosition = 0;

      //Dilution
      // {{{
      AvgDilution += TMath::Abs(dilutionFactor);
      AvgDilution_corrected += 0.95*TMath::Abs(dilutionFactor);
      AvgDilution_count++;

      Double_t correct_dil = 0.95*TMath::Abs(dilutionFactor);
      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, vDiMuon_invM, M_bounds, correct_dil);
      BinAvg(AvgDil_rad, AvgDil_rad_count, radius, rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_upstream, AvgDil_vxZ_count_UpStream, vx_z,
	     vxZ_upstream_bounds, correct_dil);

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
      BinAvg(AvgDil_rad_UpStream, AvgDil_rad_count_UpStream, radius, rad_bounds,
	     correct_dil);
      // }}}

      //Polarization
      // {{{
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
      BinAvg(AvgPol_rad_UpStream, AvgPol_rad_count_UpStream, radius,
	     rad_bounds, pol);
      BinAvg(AvgPol_vxZ_upstream,AvgPol_vxZ_count_UpStream, vx_z,
	     vxZ_upstream_bounds, pol);
      // }}}

      if (Spin[0] > 0) {//Polarized Up
	//Dilution/Polarization
	// {{{
	BinAvg(AvgDil_xN_UpStream_Up, AvgDil_xN_count_UpStream_Up, x_target,
	       xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_UpStream_Up, AvgDil_xPi_count_UpStream_Up, x_beam,
	       xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_UpStream_Up, AvgDil_xF_count_UpStream_Up, x_feynman,
	       xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_UpStream_Up,AvgDil_pT_count_UpStream_Up,q_transverse,
	       pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_UpStream_Up, AvgDil_M_count_UpStream_Up, vDiMuon_invM,
	       M_bounds, correct_dil);
	BinAvg(AvgDil_rad_UpStream_Up, AvgDil_rad_count_UpStream_Up, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_upstream_Up, AvgDil_vxZ_count_UpStream_Up, vx_z,
	       vxZ_upstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_UpStream_Up, AvgPol_xN_count_UpStream_Up, x_target,
	       xN_bounds, pol);
	BinAvg(AvgPol_xPi_UpStream_Up, AvgPol_xPi_count_UpStream_Up, x_beam,
	       xPi_bounds, pol);
	BinAvg(AvgPol_xF_UpStream_Up, AvgPol_xF_count_UpStream_Up, x_feynman,
	       xF_bounds, pol);
	BinAvg(AvgPol_pT_UpStream_Up,AvgPol_pT_count_UpStream_Up,q_transverse,
	       pT_bounds, pol);  
	BinAvg(AvgPol_M_UpStream_Up, AvgPol_M_count_UpStream_Up, vDiMuon_invM,
	       M_bounds, pol);
	BinAvg(AvgPol_rad_UpStream_Up, AvgPol_rad_count_UpStream_Up, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_upstream_Up, AvgPol_vxZ_count_UpStream_Up, vx_z,
	       vxZ_upstream_bounds, pol);
	// }}}
      }
      else if (Spin[0] < 0){//Polarized Down
	//Dilution/Polarization
	// {{{
	BinAvg(AvgDil_xN_UpStream_Down, AvgDil_xN_count_UpStream_Down, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_UpStream_Down, AvgDil_xPi_count_UpStream_Down, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_UpStream_Down, AvgDil_xF_count_UpStream_Down, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_UpStream_Down,AvgDil_pT_count_UpStream_Down,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_UpStream_Down, AvgDil_M_count_UpStream_Down,
	       vDiMuon_invM, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_UpStream_Down, AvgDil_rad_count_UpStream_Down, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_upstream_Down, AvgDil_vxZ_count_UpStream_Down, vx_z,
	       vxZ_upstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_UpStream_Down, AvgPol_xN_count_UpStream_Down, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_UpStream_Down, AvgPol_xPi_count_UpStream_Down, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_UpStream_Down, AvgPol_xF_count_UpStream_Down, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_UpStream_Down,AvgPol_pT_count_UpStream_Down,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_UpStream_Down, AvgPol_M_count_UpStream_Down,
	       vDiMuon_invM, M_bounds, pol);
	BinAvg(AvgPol_rad_UpStream_Down, AvgPol_rad_count_UpStream_Down, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_upstream_Down, AvgPol_vxZ_count_UpStream_Down, vx_z,
	       vxZ_upstream_bounds, pol);
	// }}}
      }
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

      //Dilution
      // {{{
      AvgDilution += TMath::Abs(dilutionFactor);
      AvgDilution_corrected += 0.91*TMath::Abs(dilutionFactor);
      AvgDilution_count++;

      Double_t correct_dil = 0.91*TMath::Abs(dilutionFactor);
      BinAvg(AvgDil_xN, AvgDil_xN_count, x_target, xN_bounds, correct_dil);
      BinAvg(AvgDil_xPi, AvgDil_xPi_count, x_beam, xPi_bounds, correct_dil);
      BinAvg(AvgDil_xF, AvgDil_xF_count, x_feynman, xF_bounds, correct_dil);
      BinAvg(AvgDil_pT, AvgDil_pT_count, q_transverse, pT_bounds, correct_dil);
      BinAvg(AvgDil_M, AvgDil_M_count, vDiMuon_invM, M_bounds, correct_dil);
      BinAvg(AvgDil_rad, AvgDil_rad_count, radius, rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_downstream, AvgDil_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, correct_dil);

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
      BinAvg(AvgDil_rad_DownStream, AvgDil_rad_count_DownStream, radius,
	     rad_bounds, correct_dil);
      BinAvg(AvgDil_vxZ_downstream,AvgDil_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, correct_dil);
      // }}}

      //Polarization
      // {{{
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
      BinAvg(AvgPol_rad_DownStream, AvgPol_rad_count_DownStream, radius,
	     rad_bounds, pol);
      BinAvg(AvgPol_vxZ_downstream,AvgPol_vxZ_count_DownStream, vx_z,
	     vxZ_downstream_bounds, pol);
      // }}}

      if (Spin[0] > 0) {//Polarized Up
	//Dilution/Polarization
	// {{{
	BinAvg(AvgDil_xN_DownStream_Up, AvgDil_xN_count_DownStream_Up, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_DownStream_Up, AvgDil_xPi_count_DownStream_Up, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_DownStream_Up, AvgDil_xF_count_DownStream_Up, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_DownStream_Up,AvgDil_pT_count_DownStream_Up,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_DownStream_Up, AvgDil_M_count_DownStream_Up,
	       vDiMuon_invM, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_DownStream_Up, AvgDil_rad_count_DownStream_Up, radius,
	       rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_downstream_Up, AvgDil_vxZ_count_DownStream_Up, vx_z,
	       vxZ_downstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_DownStream_Up, AvgPol_xN_count_DownStream_Up, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_DownStream_Up, AvgPol_xPi_count_DownStream_Up, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_DownStream_Up, AvgPol_xF_count_DownStream_Up, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_DownStream_Up,AvgPol_pT_count_DownStream_Up,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_DownStream_Up, AvgPol_M_count_DownStream_Up,
	       vDiMuon_invM, M_bounds, pol);
	BinAvg(AvgPol_rad_DownStream_Up, AvgPol_rad_count_DownStream_Up, radius,
	       rad_bounds, pol);
	BinAvg(AvgPol_vxZ_downstream_Up, AvgPol_vxZ_count_DownStream_Up, vx_z,
	       vxZ_downstream_bounds, pol);
	// }}}
      }
      else if (Spin[0] < 0){//Polarized Down
	//Dilution/Polarization
	// {{{
	BinAvg(AvgDil_xN_DownStream_Down, AvgDil_xN_count_DownStream_Down, 
	       x_target, xN_bounds, correct_dil);
	BinAvg(AvgDil_xPi_DownStream_Down, AvgDil_xPi_count_DownStream_Down, 
	       x_beam, xPi_bounds, correct_dil);
	BinAvg(AvgDil_xF_DownStream_Down, AvgDil_xF_count_DownStream_Down, 
	       x_feynman, xF_bounds, correct_dil);
	BinAvg(AvgDil_pT_DownStream_Down,AvgDil_pT_count_DownStream_Down,
	       q_transverse, pT_bounds, correct_dil);  
	BinAvg(AvgDil_M_DownStream_Down, AvgDil_M_count_DownStream_Down,
	       vDiMuon_invM, M_bounds, correct_dil);
	BinAvg(AvgDil_rad_DownStream_Down, AvgDil_rad_count_DownStream_Down,
	       radius, rad_bounds, correct_dil);
	BinAvg(AvgDil_vxZ_downstream_Down, AvgDil_vxZ_count_DownStream_Down,
	       vx_z, vxZ_downstream_bounds, correct_dil);

	BinAvg(AvgPol_xN_DownStream_Down, AvgPol_xN_count_DownStream_Down, 
	       x_target, xN_bounds, pol);
	BinAvg(AvgPol_xPi_DownStream_Down, AvgPol_xPi_count_DownStream_Down, 
	       x_beam, xPi_bounds, pol);
	BinAvg(AvgPol_xF_DownStream_Down, AvgPol_xF_count_DownStream_Down, 
	       x_feynman, xF_bounds, pol);
	BinAvg(AvgPol_pT_DownStream_Down,AvgPol_pT_count_DownStream_Down,
	       q_transverse, pT_bounds, pol);  
	BinAvg(AvgPol_M_DownStream_Down, AvgPol_M_count_DownStream_Down,
	       vDiMuon_invM, M_bounds, pol);
	BinAvg(AvgPol_rad_DownStream_Down, AvgPol_rad_count_DownStream_Down,
	       radius, rad_bounds, pol);
	BinAvg(AvgPol_vxZ_downstream_Down, AvgPol_vxZ_count_DownStream_Down,
	       vx_z, vxZ_downstream_bounds, pol);
	// }}}
      }
    }//Down Stream

    //Setup Vectors in different coordinate systems
    // {{{
    //Performed after bad spills cut
    /*align_wrt_beam_photon(lv_beam_TF, lv_target_TF, lv_Spin_TF,
			  lv_Spin_simple_TF, lv_virtualPhoton_TF, lv_muPlus_TF,
			  lv_muMinus_TF);*/


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
    // }}}

    Double_t PhiS_lab = lv_Spin.Phi() - lv_diMu.Phi();
    if(PhiS_lab > TMath::Pi()) PhiS_lab = -2*TMath::Pi() + PhiS_lab;
    else if (PhiS_lab < -1.0*TMath::Pi() ) PhiS_lab = 2*TMath::Pi() + PhiS_lab;
    PhiS = lv_Spin_TF.Phi();
    PhiS_simple = lv_Spin_simple_TF.Phi();
    Phi_CS = lv_muMinus_CS.Phi();
    Theta_CS = lv_muMinus_CS.Theta();
    rapidity = 0.5*TMath::Log(x_beam/x_target);

    //Average
    AvgPolarization += TMath::Abs(Polarization);
    AvgPolarization_count++;
    if(Polarization == 0.0 || dilutionFactor == 0.0){
      cout << "Problems with Polarization or dilution value" << endl;
      cout << "Polarization = " << Polarization << endl;
      cout << "Dilution = " << dilutionFactor << endl;
    }

    Double_t pol = TMath::Abs(Polarization);
    BinAvg(AvgPol_xN, AvgPol_xN_count, x_target, xN_bounds, pol);
    BinAvg(AvgPol_xPi, AvgPol_xPi_count, x_beam, xPi_bounds, pol);
    BinAvg(AvgPol_xF, AvgPol_xF_count, x_feynman, xF_bounds, pol);
    BinAvg(AvgPol_pT, AvgPol_pT_count, q_transverse, pT_bounds, pol);
    BinAvg(AvgPol_M, AvgPol_M_count, vDiMuon_invM, M_bounds, pol);
    BinAvg(AvgPol_rad, AvgPol_rad_count, radius, rad_bounds, pol);
    
    tree->Fill();
  }//tree entries


  //Dilution and Polarization TVectorD
  // {{{
  TVectorD Dil_xN(nBins), Dil_xN_UpStream(nBins), Dil_xN_DownStream(nBins);
  TVectorD Dil_xN_UpStream_Up(nBins), Dil_xN_DownStream_Up(nBins);
  TVectorD Dil_xN_UpStream_Down(nBins), Dil_xN_DownStream_Down(nBins);
  TVectorD Pol_xN(nBins), Pol_xN_UpStream(nBins), Pol_xN_DownStream(nBins);
  TVectorD Pol_xN_UpStream_Up(nBins), Pol_xN_DownStream_Up(nBins);
  TVectorD Pol_xN_UpStream_Down(nBins), Pol_xN_DownStream_Down(nBins);

  TVectorD Dil_xPi(nBins), Dil_xPi_UpStream(nBins), Dil_xPi_DownStream(nBins);
  TVectorD Dil_xPi_UpStream_Up(nBins), Dil_xPi_DownStream_Up(nBins);
  TVectorD Dil_xPi_UpStream_Down(nBins), Dil_xPi_DownStream_Down(nBins);
  TVectorD Pol_xPi(nBins), Pol_xPi_UpStream(nBins), Pol_xPi_DownStream(nBins);
  TVectorD Pol_xPi_UpStream_Up(nBins), Pol_xPi_DownStream_Up(nBins);
  TVectorD Pol_xPi_UpStream_Down(nBins), Pol_xPi_DownStream_Down(nBins);
  
  TVectorD Dil_xF(nBins), Dil_xF_UpStream(nBins), Dil_xF_DownStream(nBins);
  TVectorD Dil_xF_UpStream_Up(nBins), Dil_xF_DownStream_Up(nBins);
  TVectorD Dil_xF_UpStream_Down(nBins), Dil_xF_DownStream_Down(nBins);
  TVectorD Pol_xF(nBins), Pol_xF_UpStream(nBins), Pol_xF_DownStream(nBins);
  TVectorD Pol_xF_UpStream_Up(nBins), Pol_xF_DownStream_Up(nBins);
  TVectorD Pol_xF_UpStream_Down(nBins), Pol_xF_DownStream_Down(nBins);

  TVectorD Dil_pT(nBins), Dil_pT_UpStream(nBins), Dil_pT_DownStream(nBins);
  TVectorD Dil_pT_UpStream_Up(nBins), Dil_pT_DownStream_Up(nBins);
  TVectorD Dil_pT_UpStream_Down(nBins), Dil_pT_DownStream_Down(nBins);
  TVectorD Pol_pT(nBins), Pol_pT_UpStream(nBins), Pol_pT_DownStream(nBins);
  TVectorD Pol_pT_UpStream_Up(nBins), Pol_pT_DownStream_Up(nBins);
  TVectorD Pol_pT_UpStream_Down(nBins), Pol_pT_DownStream_Down(nBins);

  TVectorD Dil_M(nBins), Dil_M_UpStream(nBins), Dil_M_DownStream(nBins);
  TVectorD Dil_M_UpStream_Up(nBins), Dil_M_DownStream_Up(nBins);
  TVectorD Dil_M_UpStream_Down(nBins), Dil_M_DownStream_Down(nBins);
  TVectorD Pol_M(nBins), Pol_M_UpStream(nBins), Pol_M_DownStream(nBins);
  TVectorD Pol_M_UpStream_Up(nBins), Pol_M_DownStream_Up(nBins);
  TVectorD Pol_M_UpStream_Down(nBins), Pol_M_DownStream_Down(nBins);

  TVectorD Dil_rad(nBins), Dil_rad_UpStream(nBins), Dil_rad_DownStream(nBins);
  TVectorD Dil_rad_UpStream_Up(nBins), Dil_rad_DownStream_Up(nBins);
  TVectorD Dil_rad_UpStream_Down(nBins), Dil_rad_DownStream_Down(nBins);
  TVectorD Pol_rad(nBins), Pol_rad_UpStream(nBins), Pol_rad_DownStream(nBins);
  TVectorD Pol_rad_UpStream_Up(nBins), Pol_rad_DownStream_Up(nBins);
  TVectorD Pol_rad_UpStream_Down(nBins), Pol_rad_DownStream_Down(nBins);

  TVectorD Dil_vxZ_upstream(nBins);
  TVectorD Dil_vxZ_upstream_Up(nBins);
  TVectorD Dil_vxZ_upstream_Down(nBins);
  TVectorD Pol_vxZ_upstream(nBins);
  TVectorD Pol_vxZ_upstream_Up(nBins);
  TVectorD Pol_vxZ_upstream_Down(nBins);

  TVectorD Dil_vxZ_downstream(nBins);
  TVectorD Dil_vxZ_downstream_Up(nBins);
  TVectorD Dil_vxZ_downstream_Down(nBins);
  TVectorD Pol_vxZ_downstream(nBins);
  TVectorD Pol_vxZ_downstream_Up(nBins);
  TVectorD Pol_vxZ_downstream_Down(nBins);

  for (Int_t i=0; i<nBins; i++) {
    //Dilution
    ///////////////
    Dil_xN[i] = AvgDil_xN[i]/AvgDil_xN_count[i];
    Dil_xPi[i] = AvgDil_xPi[i]/AvgDil_xPi_count[i];
    Dil_xF[i] = AvgDil_xF[i]/AvgDil_xF_count[i];
    Dil_pT[i] = AvgDil_pT[i]/AvgDil_pT_count[i];
    Dil_M[i] = AvgDil_M[i]/AvgDil_M_count[i];
    Dil_rad[i] = AvgDil_rad[i]/AvgDil_rad_count[i];
    Dil_vxZ_upstream[i] = AvgDil_vxZ_upstream[i]/AvgDil_vxZ_count_UpStream[i];
    Dil_vxZ_downstream[i] =
      AvgDil_vxZ_downstream[i]/AvgDil_vxZ_count_DownStream[i];

    //UpStream
    Dil_xN_UpStream[i] = AvgDil_xN_UpStream[i]/AvgDil_xN_count_UpStream[i];
    Dil_xPi_UpStream[i] = AvgDil_xPi_UpStream[i]/AvgDil_xPi_count_UpStream[i];
    Dil_xF_UpStream[i] = AvgDil_xF_UpStream[i]/AvgDil_xF_count_UpStream[i];
    Dil_pT_UpStream[i] = AvgDil_pT_UpStream[i]/AvgDil_pT_count_UpStream[i];
    Dil_M_UpStream[i] = AvgDil_M_UpStream[i]/AvgDil_M_count_UpStream[i];
    Dil_rad_UpStream[i] = AvgDil_rad_UpStream[i]/AvgDil_rad_count_UpStream[i];

    //UpStream Polarized Up
    Dil_xN_UpStream_Up[i] = 
      AvgDil_xN_UpStream_Up[i]/AvgDil_xN_count_UpStream_Up[i];
    Dil_xPi_UpStream_Up[i] = 
      AvgDil_xPi_UpStream_Up[i]/AvgDil_xPi_count_UpStream_Up[i];
    Dil_xF_UpStream_Up[i] = 
      AvgDil_xF_UpStream_Up[i]/AvgDil_xF_count_UpStream_Up[i];
    Dil_pT_UpStream_Up[i] = 
      AvgDil_pT_UpStream_Up[i]/AvgDil_pT_count_UpStream_Up[i];
    Dil_M_UpStream_Up[i] = 
      AvgDil_M_UpStream_Up[i]/AvgDil_M_count_UpStream_Up[i];
    Dil_rad_UpStream_Up[i] = 
      AvgDil_rad_UpStream_Up[i]/AvgDil_rad_count_UpStream_Up[i];
    Dil_vxZ_upstream_Up[i] = 
      AvgDil_vxZ_upstream_Up[i]/AvgDil_vxZ_count_UpStream_Up[i];

    //UpStream Polarized Down
    Dil_xN_UpStream_Down[i] = 
      AvgDil_xN_UpStream_Down[i]/AvgDil_xN_count_UpStream_Down[i];
    Dil_xPi_UpStream_Down[i] = 
      AvgDil_xPi_UpStream_Down[i]/AvgDil_xPi_count_UpStream_Down[i];
    Dil_xF_UpStream_Down[i] = 
      AvgDil_xF_UpStream_Down[i]/AvgDil_xF_count_UpStream_Down[i];
    Dil_pT_UpStream_Down[i] =
      AvgDil_pT_UpStream_Down[i]/AvgDil_pT_count_UpStream_Down[i];
    Dil_M_UpStream_Down[i] =
      AvgDil_M_UpStream_Down[i]/AvgDil_M_count_UpStream_Down[i];
    Dil_rad_UpStream_Down[i] = 
      AvgDil_rad_UpStream_Down[i]/AvgDil_rad_count_UpStream_Down[i];
    Dil_vxZ_upstream_Down[i] = 
      AvgDil_vxZ_upstream_Down[i]/AvgDil_vxZ_count_UpStream_Down[i];

    //DownStream
    Dil_xN_DownStream[i] =AvgDil_xN_DownStream[i]/AvgDil_xN_count_DownStream[i];
    Dil_xPi_DownStream[i]=
      AvgDil_xPi_DownStream[i]/AvgDil_xPi_count_DownStream[i];
    Dil_xF_DownStream[i] =AvgDil_xF_DownStream[i]/AvgDil_xF_count_DownStream[i];
    Dil_pT_DownStream[i] =AvgDil_pT_DownStream[i]/AvgDil_pT_count_DownStream[i];
    Dil_M_DownStream[i] = AvgDil_M_DownStream[i]/AvgDil_M_count_DownStream[i];
    Dil_rad_DownStream[i] =
      AvgDil_rad_DownStream[i]/AvgDil_rad_count_DownStream[i];

    //DownStream Polarized Up
    Dil_xN_DownStream_Up[i] =
      AvgDil_xN_DownStream_Up[i]/AvgDil_xN_count_DownStream_Up[i];
    Dil_xPi_DownStream_Up[i]=
      AvgDil_xPi_DownStream_Up[i]/AvgDil_xPi_count_DownStream_Up[i];
    Dil_xF_DownStream_Up[i] =
      AvgDil_xF_DownStream_Up[i]/AvgDil_xF_count_DownStream_Up[i];
    Dil_pT_DownStream_Up[i] =
      AvgDil_pT_DownStream_Up[i]/AvgDil_pT_count_DownStream_Up[i];
    Dil_M_DownStream_Up[i] =
      AvgDil_M_DownStream_Up[i]/AvgDil_M_count_DownStream_Up[i];
    Dil_rad_DownStream_Up[i] = 
      AvgDil_rad_DownStream_Up[i]/AvgDil_rad_count_DownStream_Up[i];
    Dil_vxZ_downstream_Up[i] = 
      AvgDil_vxZ_downstream_Up[i]/AvgDil_vxZ_count_DownStream_Up[i];
	
    //DownStream Polarized Down
    Dil_xN_DownStream_Down[i] =
      AvgDil_xN_DownStream_Down[i]/AvgDil_xN_count_DownStream_Down[i];
    Dil_xPi_DownStream_Down[i]=
      AvgDil_xPi_DownStream_Down[i]/AvgDil_xPi_count_DownStream_Down[i];
    Dil_xF_DownStream_Down[i] =
      AvgDil_xF_DownStream_Down[i]/AvgDil_xF_count_DownStream_Down[i];
    Dil_pT_DownStream_Down[i] =
      AvgDil_pT_DownStream_Down[i]/AvgDil_pT_count_DownStream_Down[i];
    Dil_M_DownStream_Down[i] =
      AvgDil_M_DownStream_Down[i]/AvgDil_M_count_DownStream_Down[i];
    Dil_rad_DownStream_Down[i] = 
      AvgDil_rad_DownStream_Down[i]/AvgDil_rad_count_DownStream_Down[i];
    Dil_vxZ_downstream_Down[i] = 
      AvgDil_vxZ_downstream_Down[i]/AvgDil_vxZ_count_DownStream_Down[i];


    //Polarization
    ////////////////
    Pol_xN[i] = AvgPol_xN[i]/AvgPol_xN_count[i];
    Pol_xPi[i] = AvgPol_xPi[i]/AvgPol_xPi_count[i];
    Pol_xF[i] = AvgPol_xF[i]/AvgPol_xF_count[i];
    Pol_pT[i] = AvgPol_pT[i]/AvgPol_pT_count[i];
    Pol_M[i] = AvgPol_M[i]/AvgPol_M_count[i];
    Pol_rad[i] = AvgPol_rad[i]/AvgPol_rad_count[i];
    Pol_vxZ_upstream[i] = AvgPol_vxZ_upstream[i]/AvgPol_vxZ_count_UpStream[i];
    Pol_vxZ_downstream[i] =
      AvgPol_vxZ_downstream[i]/AvgPol_vxZ_count_DownStream[i];

    //UpStream////////////
    Pol_xN_UpStream[i] = AvgPol_xN_UpStream[i]/AvgPol_xN_count_UpStream[i];
    Pol_xPi_UpStream[i] = AvgPol_xPi_UpStream[i]/AvgPol_xPi_count_UpStream[i];
    Pol_xF_UpStream[i] = AvgPol_xF_UpStream[i]/AvgPol_xF_count_UpStream[i];
    Pol_pT_UpStream[i] = AvgPol_pT_UpStream[i]/AvgPol_pT_count_UpStream[i];
    Pol_M_UpStream[i] = AvgPol_M_UpStream[i]/AvgPol_M_count_UpStream[i];
    Pol_rad_UpStream[i] = AvgPol_rad_UpStream[i]/AvgPol_rad_count_UpStream[i];

    //UpStream Polarized Up
    Pol_xN_UpStream_Up[i] = 
      AvgPol_xN_UpStream_Up[i]/AvgPol_xN_count_UpStream_Up[i];
    Pol_xPi_UpStream_Up[i] = 
      AvgPol_xPi_UpStream_Up[i]/AvgPol_xPi_count_UpStream_Up[i];
    Pol_xF_UpStream_Up[i] = 
      AvgPol_xF_UpStream_Up[i]/AvgPol_xF_count_UpStream_Up[i];
    Pol_pT_UpStream_Up[i] = 
      AvgPol_pT_UpStream_Up[i]/AvgPol_pT_count_UpStream_Up[i];
    Pol_M_UpStream_Up[i] = 
      AvgPol_M_UpStream_Up[i]/AvgPol_M_count_UpStream_Up[i];
    Pol_rad_UpStream_Up[i] = 
      AvgPol_rad_UpStream_Up[i]/AvgPol_rad_count_UpStream_Up[i];
    Pol_vxZ_upstream_Up[i] = 
      AvgPol_vxZ_upstream_Up[i]/AvgPol_vxZ_count_UpStream_Up[i];

    //UpStream Polarized Down
    Pol_xN_UpStream_Down[i] = 
      AvgPol_xN_UpStream_Down[i]/AvgPol_xN_count_UpStream_Down[i];
    Pol_xPi_UpStream_Down[i] = 
      AvgPol_xPi_UpStream_Down[i]/AvgPol_xPi_count_UpStream_Down[i];
    Pol_xF_UpStream_Down[i] = 
      AvgPol_xF_UpStream_Down[i]/AvgPol_xF_count_UpStream_Down[i];
    Pol_pT_UpStream_Down[i] =
      AvgPol_pT_UpStream_Down[i]/AvgPol_pT_count_UpStream_Down[i];
    Pol_M_UpStream_Down[i] =
      AvgPol_M_UpStream_Down[i]/AvgPol_M_count_UpStream_Down[i];
    Pol_rad_UpStream_Down[i] = 
      AvgPol_rad_UpStream_Down[i]/AvgPol_rad_count_UpStream_Down[i];
    Pol_vxZ_upstream_Down[i] = 
      AvgPol_vxZ_upstream_Down[i]/AvgPol_vxZ_count_UpStream_Down[i];
	
    //DownStream////////////
    Pol_xN_DownStream[i] =AvgPol_xN_DownStream[i]/AvgPol_xN_count_DownStream[i];
    Pol_xPi_DownStream[i]=
      AvgPol_xPi_DownStream[i]/AvgPol_xPi_count_DownStream[i];
    Pol_xF_DownStream[i] =AvgPol_xF_DownStream[i]/AvgPol_xF_count_DownStream[i];
    Pol_pT_DownStream[i] =AvgPol_pT_DownStream[i]/AvgPol_pT_count_DownStream[i];
    Pol_M_DownStream[i] = AvgPol_M_DownStream[i]/AvgPol_M_count_DownStream[i];
    Pol_rad_DownStream[i] =
      AvgPol_rad_DownStream[i]/AvgPol_rad_count_DownStream[i];

    //DownStream Polarized Up
    Pol_xN_DownStream_Up[i] =
      AvgPol_xN_DownStream_Up[i]/AvgPol_xN_count_DownStream_Up[i];
    Pol_xPi_DownStream_Up[i]=
      AvgPol_xPi_DownStream_Up[i]/AvgPol_xPi_count_DownStream_Up[i];
    Pol_xF_DownStream_Up[i] =
      AvgPol_xF_DownStream_Up[i]/AvgPol_xF_count_DownStream_Up[i];
    Pol_pT_DownStream_Up[i] =
      AvgPol_pT_DownStream_Up[i]/AvgPol_pT_count_DownStream_Up[i];
    Pol_M_DownStream_Up[i] =
      AvgPol_M_DownStream_Up[i]/AvgPol_M_count_DownStream_Up[i];
    Pol_rad_DownStream_Up[i] = 
      AvgPol_rad_DownStream_Up[i]/AvgPol_rad_count_DownStream_Up[i];
    Pol_vxZ_downstream_Up[i] = 
      AvgPol_vxZ_downstream_Up[i]/AvgPol_vxZ_count_DownStream_Up[i];
    
	
    //DownStream Polarized Down
    Pol_xN_DownStream_Down[i] =
      AvgPol_xN_DownStream_Down[i]/AvgPol_xN_count_DownStream_Down[i];
    Pol_xPi_DownStream_Down[i]=
      AvgPol_xPi_DownStream_Down[i]/AvgPol_xPi_count_DownStream_Down[i];
    Pol_xF_DownStream_Down[i] =
      AvgPol_xF_DownStream_Down[i]/AvgPol_xF_count_DownStream_Down[i];
    Pol_pT_DownStream_Down[i] =
      AvgPol_pT_DownStream_Down[i]/AvgPol_pT_count_DownStream_Down[i];
    Pol_M_DownStream_Down[i] =
      AvgPol_M_DownStream_Down[i]/AvgPol_M_count_DownStream_Down[i];
    Pol_rad_DownStream_Down[i] = 
      AvgPol_rad_DownStream_Down[i]/AvgPol_rad_count_DownStream_Down[i];
    Pol_vxZ_downstream_Down[i] = 
      AvgPol_vxZ_downstream_Down[i]/AvgPol_vxZ_count_DownStream_Down[i];
	
  }//Dilution and polariztion setup loop

  TVectorD Dil_int(1), Pol_int(1);
  Dil_int[0] = AvgDilution_corrected/AvgDilution_count;
  Pol_int[0] = AvgPolarization/AvgPolarization_count;
  // }}}

  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  //Cuts histogram
  TString cutNames[nRealCuts] = {"AllData", "GoodSpills", "xPion,xN,xF",
				 "0.4<qT<5", "dilutionFactor",
				 "TargetZ-cut", "TargetRadius"};
  for (Int_t i=0, j=1; i<nRealCuts; i++, j+=cut_space){
    Int_t bin_index = hCuts->GetXaxis()->FindBin(j);
    hCuts->GetXaxis()->SetBinLabel(bin_index, cutNames[i]);
  }

  //Write Output
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
    TDirectory *PhiS_CutImpact = myFile->mkdir("PhiS_CutImpact");
    TDirectory *PhiS_simple_CutImpact = myFile->mkdir("PhiS_simple_CutImpact");
    
    TDirectory *MuP_PxPy_CutImpact = myFile->mkdir("MuP_PxPy_CutImpact");
    TDirectory *MuM_PxPy_CutImpact = myFile->mkdir("MuM_PxPy_CutImpact");
    TDirectory *Beam_PxPy_CutImpact = myFile->mkdir("Beam_PxPy_CutImpact");
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

      PhiPhoton_CutImpact->cd();
      hCut_PhiPhoton[i]->Write(cutNames[i]);
      PhiS_CutImpact->cd();
      hCut_PhiS[i]->Write(cutNames[i]);
      PhiS_simple_CutImpact->cd();
      hCut_PhiS_simple[i]->Write(cutNames[i]);

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
    tv_rad_bounds.Write("tv_rad_bounds");
    tv_vxZ_upstream_bounds.Write("tv_vxZ_upstream_bounds");
    tv_vxZ_downstream_bounds.Write("tv_vxZ_downstream_bounds");

    tv_xN_xval.Write("tv_xN_xval");
    tv_xPi_xval.Write("tv_xPi_xval");
    tv_xF_xval.Write("tv_xF_xval");
    tv_pT_xval.Write("tv_pT_xval");
    tv_M_xval.Write("tv_M_xval");
    tv_rad_xval.Write("tv_rad_xval");
    tv_vxZ_upstream_xval.Write("tv_vxZ_upstream_xval");
    tv_vxZ_downstream_xval.Write("tv_vxZ_downstream_xval");

    Dil_int.Write("Dil_int");
    Dil_xN.Write("Dil_xN");
    Dil_xPi.Write("Dil_xPi");
    Dil_xF.Write("Dil_xF");
    Dil_pT.Write("Dil_pT");
    Dil_M.Write("Dil_M");
    Dil_rad.Write("Dil_rad");
    Dil_vxZ_upstream.Write("Dil_vxZ_upstream");
    Dil_vxZ_downstream.Write("Dil_vxZ_downstream");
    
    Dil_xN_UpStream.Write("Dil_xN_UpStream");
    Dil_xPi_UpStream.Write("Dil_xPi_UpStream");
    Dil_xF_UpStream.Write("Dil_xF_UpStream");
    Dil_pT_UpStream.Write("Dil_pT_UpStream");
    Dil_M_UpStream.Write("Dil_M_UpStream");
    Dil_rad_UpStream.Write("Dil_rad_UpStream");

    Dil_xN_UpStream_Up.Write("Dil_xN_UpStream_Up");
    Dil_xPi_UpStream_Up.Write("Dil_xPi_UpStream_Up");
    Dil_xF_UpStream_Up.Write("Dil_xF_UpStream_Up");
    Dil_pT_UpStream_Up.Write("Dil_pT_UpStream_Up");
    Dil_M_UpStream_Up.Write("Dil_M_UpStream_Up");
    Dil_rad_UpStream_Up.Write("Dil_rad_UpStream_Up");
    Dil_vxZ_upstream_Up.Write("Dil_vxZ_upstream_Up");
	
    Dil_xN_UpStream_Down.Write("Dil_xN_UpStream_Down");
    Dil_xPi_UpStream_Down.Write("Dil_xPi_UpStream_Down");
    Dil_xF_UpStream_Down.Write("Dil_xF_UpStream_Down");
    Dil_pT_UpStream_Down.Write("Dil_pT_UpStream_Down");
    Dil_M_UpStream_Down.Write("Dil_M_UpStream_Down");
    Dil_rad_UpStream_Down.Write("Dil_rad_UpStream_Down");
    Dil_vxZ_upstream_Down.Write("Dil_vxZ_upstream_Down");

    Dil_xN_DownStream.Write("Dil_xN_DownStream");
    Dil_xPi_DownStream.Write("Dil_xPi_DownStream");
    Dil_xF_DownStream.Write("Dil_xF_DownStream");
    Dil_pT_DownStream.Write("Dil_pT_DownStream");
    Dil_M_DownStream.Write("Dil_M_DownStream");
    Dil_rad_DownStream.Write("Dil_rad_DownStream");

    Dil_xN_DownStream_Up.Write("Dil_xN_DownStream_Up");
    Dil_xPi_DownStream_Up.Write("Dil_xPi_DownStream_Up");
    Dil_xF_DownStream_Up.Write("Dil_xF_DownStream_Up");
    Dil_pT_DownStream_Up.Write("Dil_pT_DownStream_Up");
    Dil_M_DownStream_Up.Write("Dil_M_DownStream_Up");
    Dil_rad_DownStream_Up.Write("Dil_rad_DownStream_Up");
    Dil_vxZ_downstream_Up.Write("Dil_vxZ_downstream_Up");
	
    Dil_xN_DownStream_Down.Write("Dil_xN_DownStream_Down");
    Dil_xPi_DownStream_Down.Write("Dil_xPi_DownStream_Down");
    Dil_xF_DownStream_Down.Write("Dil_xF_DownStream_Down");
    Dil_pT_DownStream_Down.Write("Dil_pT_DownStream_Down");
    Dil_M_DownStream_Down.Write("Dil_M_DownStream_Down");
    Dil_rad_DownStream_Down.Write("Dil_rad_DownStream_Down");
    Dil_vxZ_downstream_Down.Write("Dil_vxZ_downstream_Down");

    Pol_int.Write("Pol_int");
    Pol_xN.Write("Pol_xN");
    Pol_xPi.Write("Pol_xPi");
    Pol_xF.Write("Pol_xF");
    Pol_pT.Write("Pol_pT");
    Pol_M.Write("Pol_M");
    Pol_rad.Write("Pol_rad");
    Pol_vxZ_upstream.Write("Pol_vxZ_upstream");
    Pol_vxZ_downstream.Write("Pol_vxZ_downstream");

    Pol_xN_UpStream.Write("Pol_xN_UpStream");
    Pol_xPi_UpStream.Write("Pol_xPi_UpStream");
    Pol_xF_UpStream.Write("Pol_xF_UpStream");
    Pol_pT_UpStream.Write("Pol_pT_UpStream");
    Pol_M_UpStream.Write("Pol_M_UpStream");
    Pol_rad_UpStream.Write("Pol_rad_UpStream");

    Pol_xN_UpStream_Up.Write("Pol_xN_UpStream_Up");
    Pol_xPi_UpStream_Up.Write("Pol_xPi_UpStream_Up");
    Pol_xF_UpStream_Up.Write("Pol_xF_UpStream_Up");
    Pol_pT_UpStream_Up.Write("Pol_pT_UpStream_Up");
    Pol_M_UpStream_Up.Write("Pol_M_UpStream_Up");
    Pol_rad_UpStream_Up.Write("Pol_rad_UpStream_Up");
    Pol_vxZ_upstream_Up.Write("Pol_vxZ_upstream_Up");    
	
    Pol_xN_UpStream_Down.Write("Pol_xN_UpStream_Down");
    Pol_xPi_UpStream_Down.Write("Pol_xPi_UpStream_Down");
    Pol_xF_UpStream_Down.Write("Pol_xF_UpStream_Down");
    Pol_pT_UpStream_Down.Write("Pol_pT_UpStream_Down");
    Pol_M_UpStream_Down.Write("Pol_M_UpStream_Down");
    Pol_rad_UpStream_Down.Write("Pol_rad_UpStream_Down");
    Pol_vxZ_upstream_Down.Write("Pol_vxZ_upstream_Down");    

    Pol_xN_DownStream.Write("Pol_xN_DownStream");
    Pol_xPi_DownStream.Write("Pol_xPi_DownStream");
    Pol_xF_DownStream.Write("Pol_xF_DownStream");
    Pol_pT_DownStream.Write("Pol_pT_DownStream");
    Pol_M_DownStream.Write("Pol_M_DownStream");
    Pol_rad_DownStream.Write("Pol_rad_DownStream");

    Pol_xN_DownStream_Up.Write("Pol_xN_DownStream_Up");
    Pol_xPi_DownStream_Up.Write("Pol_xPi_DownStream_Up");
    Pol_xF_DownStream_Up.Write("Pol_xF_DownStream_Up");
    Pol_pT_DownStream_Up.Write("Pol_pT_DownStream_Up");
    Pol_M_DownStream_Up.Write("Pol_M_DownStream_Up");
    Pol_rad_DownStream_Up.Write("Pol_rad_DownStream_Up");
    Pol_vxZ_downstream_Up.Write("Pol_vxZ_downstream_Up");    
	
    Pol_xN_DownStream_Down.Write("Pol_xN_DownStream_Down");
    Pol_xPi_DownStream_Down.Write("Pol_xPi_DownStream_Down");
    Pol_xF_DownStream_Down.Write("Pol_xF_DownStream_Down");
    Pol_pT_DownStream_Down.Write("Pol_pT_DownStream_Down");
    Pol_M_DownStream_Down.Write("Pol_M_DownStream_Down");
    Pol_rad_DownStream_Down.Write("Pol_rad_DownStream_Down");
    Pol_vxZ_downstream_Down.Write("Pol_vxZ_downstream_Down");    

    cout << myFile->GetName() << " was written" << endl;
    myFile->Close();
  }
  // }}}

  theApp.Run();//Needed to make root graphics work on C++
}//main

#include "common.hxx"
#include "functions.h"

using namespace std;

int main(int argc, char **argv){
  if(argc < 2){
    cout << "To be used with Generated MC data only!!!!" << endl;
    cout << "To determine the bin bounds and average bin values";
    cout << " for the given number of Bins" << endl;
    cout << "" << endl;
    cout << "Usage:" << endl;
    cout << "./main [options] [-nBins] [-ffilename]" << endl;
    cout << "filename should be the full path name" << endl;
    cout << "" << endl;
    cout << "---Needed Options---" << endl;
    cout << "Option:  -n bins           (How many bins to make)" << endl;
    cout << "" << endl;
    cout << "---Write Option---" << endl;
    cout << "Option:  -Q outName	(write output to file to outName)"
	 << endl;
    cout << "    Otherwise bin values are output and over written to ";
    cout << "\"genBinValues.txt\"" << endl;
    cout << "" << endl;
    cout << "---Additional Options---" << endl;
    cout << "Option:  -D        (debug mode, only loop over 10000 events)"<<endl;
    cout << "Option:  -R        (reducde mode, only loop over 100000 events)"<<endl;
	cout << "             (To be used with large data sets)" << endl;
	cout << "      (-R and -D options cannot be used at the same time)" << endl;
    cout << "" << endl;

    exit(EXIT_FAILURE);
  }
  cout << "" << endl;
  TApplication theApp("tapp", &argc, argv);

  //Read input arguments
  Int_t Qflag=0, fflag=0, binflag=0;
  Int_t Dflag=0, Rflag=0;
  Int_t c;
  TString userNum = "", fname = "", outFile = "";
  Int_t nBins;

  while ((c = getopt (argc, argv, "n:f:Q:DR")) != -1) {
    switch (c) {
    case 'n':
      binflag = 1;
      nBins = atoi(optarg);
      break;
    case 'Q':
      Qflag = 1;
      outFile += optarg;
      break;
    case 'f':
      fflag = 1;
      fname += optarg;
      break;
    case 'D':
      Dflag = 1;
      break;
    case 'R':
      Rflag = 1;
      break;
    case '?':
      if (optopt == 'n')
	fprintf (stderr, "Option -%c requires an argument.\n", optopt);
      else if (optopt == 'f')
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

  //Basic Checks
  if(!binflag){
    cout << "number of bins must be specified: -nBins" << endl;
    exit(EXIT_FAILURE);
  }

  TChain* T1 = new TChain("Events");
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
    TList *li = f1->GetListOfKeys();
    TIter iter( li->MakeIterator() );
    while(TObject* obj = iter()){
      TKey* theKey = (TKey*)obj;
      if (strncmp (theKey->GetClassName(),"TTree",4) == 0){
	T1->Add( fname+"/"+obj->GetName()+";"+Form("%i",theKey->GetCycle())  );
      }
    }
    f1->Close();
  }
  cout << "" << endl;

  if(Rflag && Dflag){
	  cout << "-R and -D options cannot be used at the same time!" << endl;
	  exit(EXIT_FAILURE);
  }
  else if(Dflag){
    cout << "Debug mode ====> Only 10000 events considered in tree" << endl;
  }
  else if(Rflag){
    cout << "Reduced mode ====> Only 100000 events considered in tree" << endl;
  }


  //Internal variables and binning
  Double_t M_proton = 0.938272;

  //Vectors for sorting
  vector<Double_t> sort_xTarg;
  vector<Double_t> sort_xBeam;
  vector<Double_t> sort_xF;
  vector<Double_t> sort_pT;
  vector<Double_t> sort_mass;
  vector<Double_t> sort_rad;
  vector<Double_t> sort_vxZ_upstream;
  vector<Double_t> sort_vxZ_downstream;


  //TTree Variables
  //Vertex specific
  Double_t vx_z, vx_x, vx_y;
  Int_t targetPosition;
  //Drell-Yan Angles
  Double_t PhiS_simple, Theta_CS;
  //Virtual Photon
  Double_t muM_X, muM_Y, muM_Z, muM_E;
  Double_t muP_X, muP_Y, muP_Z, muP_E;
  Double_t vPhoton_X, vPhoton_Y, vPhoton_Z, vPhoton_E;
  Double_t vOpenAngle;
  //DY-variables
  Double_t x_beam, x_target, x_feynman, q_transverse, Mmumu;

//Vertex specific
  T1->SetBranchAddress("vx_z", &vx_z);
  T1->SetBranchAddress("vx_x", &vx_x);
  T1->SetBranchAddress("vx_y", &vx_y);
  T1->SetBranchAddress("targetPosition", &targetPosition);
  //Drell-Yan Angles
  T1->SetBranchAddress("PhiS_simple", &PhiS_simple);
  T1->SetBranchAddress("Theta_CS", &Theta_CS);
  //Virtual Photon
  T1->SetBranchAddress("muM_X", &muM_X);
  T1->SetBranchAddress("muM_Y", &muM_Y);
  T1->SetBranchAddress("muM_Z", &muM_Z);
  T1->SetBranchAddress("muM_E", &muM_E);
  T1->SetBranchAddress("muP_X", &muP_X);
  T1->SetBranchAddress("muP_Y", &muP_Y);
  T1->SetBranchAddress("muP_Z", &muP_Z);
  T1->SetBranchAddress("muP_E", &muP_E);
  T1->SetBranchAddress("vPhoton_X", &vPhoton_X);
  T1->SetBranchAddress("vPhoton_Y", &vPhoton_Y);
  T1->SetBranchAddress("vPhoton_Z", &vPhoton_Z);
  T1->SetBranchAddress("vPhoton_E", &vPhoton_E);
  T1->SetBranchAddress("vOpenAngle", &vOpenAngle);
  //DY-variables
  T1->SetBranchAddress("x_beam", &x_beam);
  T1->SetBranchAddress("x_target", &x_target);
  T1->SetBranchAddress("x_feynman", &x_feynman);
  T1->SetBranchAddress("q_transverse", &q_transverse);
  T1->SetBranchAddress("Mmumu", &Mmumu);


  Int_t tree_entries;
  if (Dflag) tree_entries = 10000;
  else if (Rflag) tree_entries = 100000;
  else tree_entries = T1->GetEntries();
  cout << "Entries in tree = " << T1->GetEntries() << endl;
  cout << "Entries considered = " << tree_entries << endl;
  Bool_t first = true;
  for (Int_t ev=0; ev<tree_entries; ev++){
    T1->GetEntry(ev, 0);

    //Settings
    if (first || ev==tree_entries-1){
      cout << " " << endl;
      cout << "Setup!!!!!!!!!" << endl;
      cout << "No new cuts" << endl;

      first = false;
    }

    ////All data after cuts
    //////////////
    Double_t radius = TMath::Sqrt(vx_x*vx_x+ vx_y*vx_y);

    ////All data after cuts
    //////////////
    if (vx_z < -229.4) sort_vxZ_upstream.push_back(vx_z);//Up Stream
    else sort_vxZ_downstream.push_back(vx_z);//Down Stream

    //Both targets
    sort_xTarg.push_back(x_target);
    sort_xBeam.push_back(x_beam);
    sort_xF.push_back(x_feynman);
    sort_pT.push_back(q_transverse);
    sort_mass.push_back(Mmumu);
    sort_rad.push_back(radius);

  }//tree entries


  //Print bin boundaries and average values
  if (sort_xF.size() != 0){
    if (!Qflag) outFile += "genBinValues.txt";
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
  else{
    cout << "" << endl;
    cout << "" << endl;
    cout <<"Error:"<< endl;
    cout << "sorting vectors are 0 size" << endl;
    cout << "No file written" << endl;
    cout << "" << endl;
    cout << "" << endl;
  }

  cout << "!!!!!!!!!!!!!!!" << endl;
  cout << "Code Finished" << endl;
  cout << "!!!!!!!!!!!!!!!" << endl;

  theApp.Run();//Needed to make root graphics work on C++
}//main

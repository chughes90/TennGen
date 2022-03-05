//This file creates the interface b/w root and pythia libraries for macros you want to run in root that use
//pythia libraries

void maindriver_PreGen(Int_t nEvents = 100000000, Int_t jobID = 7301 , Int_t HF = 0, Int_t cent_bin = 0, Int_t Num_BKGD_Files = 1, Bool_t RT_Stats = kFALSE , Bool_t GRID = kFALSE ){


if( !GRID ){

}
else if( GRID){

  Num_BKGD_Files = 1; // can only run over 1 at a time on the grid - in here as a catch statement

}

  gRandom->SetSeed(jobID); //needed for random number generatrion DO NOT FORGET

//__________________________________________________________________________________________________//

//compile macro

//YOU MUST EXPORT THESE ENVIRONMENTAL VARIABLES TO RUN MY JET CODE

  gSystem->Load("libPhysics");
  gROOT->ProcessLine(".L BkgrLoad.cxx+");
  gROOT->ProcessLine(".L TennGen_PreGen_Macro.C++");

  TennGen_PreGen_Macro( nEvents, jobID , HF, cent_bin, Num_BKGD_Files, RT_Stats , GRID );

}	

//This file creates the interface b/w root and pythia libraries for macros you want to run in root that use
//pythia libraries

void maindriver(Int_t nEvents = 20, Int_t jobID =2072, Int_t HF = 0 , Int_t cent_bin = 0){

  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

  gSystem->Load("libPhysics");
  gRandom->SetSeed(jobID);
  gROOT->ProcessLine(".L TennGen.cxx+");
  //gROOT->ProcessLine(".L BkgrLoad.cxx+"); separate macro to use pre-run background events
  gROOT->ProcessLine(".L TennGen_Test_Macro.C++");

  TennGen_Test_Macro( nEvents, jobID , HF , cent_bin ) ;
}	

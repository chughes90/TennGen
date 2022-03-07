//// Author: Charles Philip Hughes <chughe26@vols.utk.edu>
// Update: 2019-09-20 14:40:27-0500
// Copyright (c) 2017; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// University of Tennessee ALICE collaboration
      

#include "TApplication.h"
#include "TFile.h"
#include "TError.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TSystemDirectory.h"
#include <cstdlib>
#include "TMath.h"
#include "TMCParticle.h"
#include "TParticle.h"
#include "TTimeStamp.h"
#include "TPythia6.h"
#include <iostream>
#include <cstring>
#include <utility>
#include "TTree.h"
#include "TVector3.h"
//#include "TFilePtr.h"
#include <algorithm>
#include "BkgrLoad.h"

using namespace std;


ClassImp(BkgrLoad)
//____________________________
//____________________________________________________________________________________________________________________//
//_________________Load the current TTree into memory (must already have called pass in settings !!!!!)________________//
//____________________________________________________________________________________________________________________//
void BkgrLoad::Load()
{

  fileEventPos = 0; //keep track of running through all the events in the root file	
  activeFile = 0;   //keep track of how many files you have run through 
  eventsRemaining = 0;
  
  cent = centralityFlag;
  harm = harmonicsFlag;
  TString filePath = datapath;


  TFile *rootFile;
  TList *fileList;
  TSystemDirectory *dir = new TSystemDirectory(filePath.Data(),filePath.Data());

  if (!filePath.EndsWith("/")) filePath = filePath + '/'; //check to see if you end with forward slash

    fileList = dir->GetListOfFiles(); //get the list of files in the directory
      if (fileList) 
      { 
	  TSystemFile *file; //create pointer to an as of yet unspecified file
          TString fname;
          TIter next(fileList); //get the next file in the list
	  while ((file=(TSystemFile*)next()))  //while there is a next file
          {
	    fname = file->GetName(); //get the file name
	    if (!file->IsDirectory() && fname.EndsWith("HF_0_cent_bin_0.root")) //make sure you are not talking about a director
	    { 
	      fileNames.push_back(fname); //add name to the list of filenames
	      TString fullPath = filePath + fname; //get the FULL path
	      //cout << fullPath << " path"<< endl;
	      UInt_t tharm, tcent, tevents, tversion;
	      TFile *rootFile = new TFile(fullPath.Data(),"READ"); //get the current root file in question
	      rootFile->cd(); //get inito the foot file
	      if (rootFile)
	      {
	        TTree *tt = (TTree*)rootFile->Get("HughesMC"); //get the header TTree
	        if (tt){
		  //cout << "Found \"HughesMC\"\n";
	          tt->SetBranchAddress("centrality",&tcent); //make the apppropriate branches
		  tt->SetBranchAddress("harmonics",&tharm);
		  tt->SetBranchAddress("events",&tevents);
		  tt->SetBranchAddress("version",&tversion);
		  tt->GetEntry(0); //get the zeroth header entry
		    if ((tcent == cent) && (tharm==harm) && (tversion==version)) //check to see if its right
		    {
                      eventsperfile = tevents;
                      //cout<<"\nEvents = "<<tevents<<endl;
                      TVector3 *eve  = new TVector3(tevents,0,0);
                      //cout<<"\nEvents = "<<eve->X()<<endl;
                      //cout<<"\nIS SEGFAULT HERE @ LINE 95\n"<<endl;
                      valid_Files->Add(rootFile);
                      //cout<<"\nIS SEGFAULT HERE @ LINE 97\n"<<endl;
		      eventsRemaining += tevents;
                      events_arr -> Add(eve);
                      //cout<<"\nIS SEGFAULT HERE @ LINE 100\n"<<endl;
                      delete eve;
		    }
		  else rootFile->Close();
		}
	      else rootFile->Close();
            }
	  }
        }
		
	//randomize order
        valid_Files->Randomize();
        numFiles = valid_Files->GetEntries();
        totalEvents = eventsRemaining;
        cout << "BGLoad found " << numFiles << " matching data files with " << eventsRemaining << " events\n";
		
      }
    else cout << "Invalid directory";
    delete dir;

}
//____________________________________________________________________________________________________________________//
//_________________Fill a TClonesArray and return it to the user (must already have called Load !!!!!)________________//
//____________________________________________________________________________________________________________________//
TClonesArray *BkgrLoad::GetEventBG()
{
  //cout<<"\nIS SEGFAULT HERE @ LINE 126\n"<<endl;
  //for(Int_t b = 0 ; b < 30; b++ ){
    TClonesArray *BGparticles = new TClonesArray("TMCParticle" , 1200);
  //cout<<"\nIS SEGFAULT HERE @ LINE 129\n"<<endl;	
    // Check if Events remain
    if (activeFile >= (numFiles)) {
      //cout << "No more events available,b = "<<b<<"\n"<<endl;
      cout << "No more events available"<<endl;
      //continue;
      return 0;
    }
  //cout<<"\nIS SEGFAULT HERE @ LINE 135\n"<<endl;
  //valid_Files->Write();
  //cout<<"\n\nvalid_Files_array length = "<<valid_Files->GetEntries()<<endl;
  //cout<<"\nIS SEGFAULT HERE @ LINE 138\n"<<endl;
  TFile *curFile = (TFile *)valid_Files->At(activeFile);
  //cout<<"\n\ncurFile name = "<<curFile->GetName()<<endl;
  TVector3 *cureve = (TVector3 *)events_arr->At(activeFile);

  curFile->cd();
	
  // if = 0 load reference indices, shuffle order
  if (fileEventPos==0) {
		
    fileEventOrder.resize(eventsperfile);
    Int_t ev = fileEventOrder.size();
    for (Int_t y=0; y < ev; y++) fileEventOrder[y] = y;
      std::random_shuffle(fileEventOrder.begin(),fileEventOrder.end());
    }

    TTree *tt = (TTree*)curFile->Get(TString::Itoa(fileEventOrder[fileEventPos],10));
    Float_t px, py, pz;//, m , E;
    Short_t kf;

    if (tt){
	//cout << "Event " << fileEventOrder[fileEventPos] << endl;
		
      eventsRemaining--;
      tt->SetBranchAddress("px",&px);
      tt->SetBranchAddress("py",&py);
      tt->SetBranchAddress("pz",&pz);
      tt->SetBranchAddress("kf",&kf);
      //tt->SetBranchAddress("m", &m );
      //tt->SetBranchAddress("E", &E );
		
      UInt_t nEntries = tt->GetEntries();
      cout<<"\n\n\n\nNum Particles/Event = "<<nEntries<<"\n\n\n\n"<<endl;
      TMCParticle* mcp;
      for (UInt_t j = 0; j < nEntries; j++){
        // get particles
        tt->GetEntry(j);
        mcp = (TMCParticle*) BGparticles->ConstructedAt(j);
        mcp->SetPx(px);
        mcp->SetPy(py);
        mcp->SetPz(pz);
        //mcp->SetMass(m);
	mcp->SetKF(kf);
	//mcp->SetEnergy(E);
        //cout<<"| px = "<<px<<" | py = "<<py<<" | pz = "<<pz<<" | m = "<<m<<" | E = "<<E<<" | KF = "<<kf<<" |"<<endl;
      }
		
    }

    fileEventPos++;
    cout<<"fileEventPos is now = "<<fileEventPos<<endl;
    if (fileEventPos==eventsperfile){
      fileEventPos = 0;
      activeFile++;
      cout<<"activeFile is now = "<<activeFile<<endl;
    }
    //cout<<"\nIS SEGFAULT HERE @ LINE 194\n"<<endl;
 // }	
    return BGparticles;
}
//____________________________________________________________________________________________________________________//
//____________________________________________________________________________________________________________________//
//_________________Fill a TClonesArray and return it to the user (must already have called Load !!!!!)________________//
//____________________________________________________________________________________________________________________//
void BkgrLoad::PrintStatus()
{
  cout << "*********************** BGLoad *****************************\n";
  cout << "Active Dir: " << filePath << endl;
  cout << "Harmonics: " << harm << " | Centrality bin: " << cent << endl;
  cout << "Found " << totalEvents << " events of which " << eventsRemaining << " remain\n";
  if (kTRUE) {
    cout << "List of applicable data files:\n";
    for (Int_t x = 0; x < numFiles; x++) cout << '\t' << fileNames[x] << endl;
  }

}
//____________________________________________________________________________________________________________________//
//______________________________________________________________________________
BkgrLoad::BkgrLoad()//constructor
{

}

//______________________________________________________________________________
BkgrLoad::~BkgrLoad()//destructor
{
   //for (Int_t x = 0; x < numFiles; x++){
     //TFile *newFile = (TFile *)valid_Files->At(x);
     //newFile->Close();
    // delete newFile;
  // }
   /// close all files
  //delete valid_Files;
  //delete events_arr;

}
//______________________________________________________________________________
void BkgrLoad::CloseFiles()//Close File Function
{
  for (Int_t x = 0; x < numFiles; x++){
    TFile *newFile = (TFile *)valid_Files->At(x);
    newFile->Close();
    delete newFile;
  }
   /// close all files
   //  //delete valid_Files;
   //    //delete events_arr;
}
//______________________________________________________________________________

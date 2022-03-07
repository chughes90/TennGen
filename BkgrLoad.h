// Author: Christian Holm Christensen <cholm@hilux15.nbi.dk>
// Update: 2002-08-16 16:40:27+0200
// Copyright: 2002 (C) Christian Holm Christensen
// Copyright (C) 2006, Rene Brun and Fons Rademakers.
// For the licensing terms see $ROOTSYS/LICENSE.
//
// 2nd Author: Richard Corke, Copyright 2013
// Roughly follows analysis of:
// * T. Aaltonen et al. [CDF Collaboration],
//
// 3rd Author: Joel Mazer <jmazer@utk.edu> (UTK collaboration)
//
// MODIFIED BY CHARLES HUGHES <chughe26@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2016-01-15 11:45 EST
//
// FURTHER MODIFIED BY JAMES NEUHAUS <jneuhau1@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2016-09-23 12:05 EST
//
// FURTHER MODIFIED BY CHARLES HUGHES <chughe26@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2019-09-15 16:04 EST

#ifndef ROOT_BkgrLoad
#define ROOT_BkgrLoad

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// BkgrLoad                                                             //
//                                                                      //
// Description of the public objects                                    //
//                                                                      // 
// void   PassInSparse(THnSparseD spars){sparsforcorr = spars;}         //
// -- pass in the sparse of correlations that we use to                 //
//    fill the pT arr to unfold                                         //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
// TH2D *GetArray();  -- get detector level pT jet vs pT assoc          //
//                        for unfolding                                 //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

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

class BkgrLoad  : public TObject {

private:
	TString filePath;
	UInt_t cent;
	UInt_t harm;
	Int_t activeFile;
	Int_t numFiles;
        Int_t eventsperfile;
	Int_t eventsRemaining;	
	std::vector<Int_t> fileEventOrder;
	Int_t fileEventPos;
	Int_t totalEvents;
	std::vector<TString> fileNames;

        TObjArray *valid_Files = new TObjArray();
        TObjArray *events_arr = new TObjArray();

        char* datapath;
        Int_t centralityFlag;
        Int_t harmonicsFlag;
        UInt_t version;
        Int_t seed;

public:
        BkgrLoad();
        virtual ~BkgrLoad();
        void PassInSettings(char* filepath,  Int_t cent_bin , Int_t HF , Int_t ver , Int_t seeed ){datapath = filepath; centralityFlag = cent_bin; harmonicsFlag = HF; version = ver; seed = seeed;};//inline function - nothing needed in .cxx, declare your settings to run through
        void Load(); //Load the trees from the files into memory
        void PrintStatus(); //Print the current status
        void CloseFiles(); //close all the damn files

        TClonesArray *GetEventBG(); //Fill the TClonesArray with TMC particles and actually return the TClonesArray to work with
        Int_t EventsRemaining() {return eventsRemaining;}; //get the number of events that remain, inline definition nothing neeeded in .cxx

        Int_t GetActiveFile() {return activeFile;}; //get the current active file, inline definition nothing needed in .cxx
        Int_t GetNumFiles() {return numFiles;}; //get the total number of files in wd, inline definition nothing needed in .cxx

        ClassDef(BkgrLoad,1)  //Event structure
}; //end class structure

#endif



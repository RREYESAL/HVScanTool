#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBrowser.h"
#include "TH1.h"
#include "TBranch.h"
#include "TChain.h"
#include "TUnixSystem.h"
#include "TCanvas.h"
#include "TGraph.h"


std::vector< std::pair<std::string,std::string> > dictionary(const char*);// function prototype



// To do a map of the rolls. Every roll has an id match a name   
std::vector< std::pair<std::string,std::string> > dictionary(const char* subd){
  ifstream f;
  Double_t id;
  std::string name;
  if (strcmp(subd,"endcap")==0)f.open("../data/detIdEndCap.txt");
  else f.open("../data/detIdBarrel.txt");
  std::vector< std::pair<std::string,std::string> > dictionary_;
  while (1) {
        f >> name
          >> id;
        if (f.eof()) break;
        stringstream s; std::string stid;
        s << UInt_t(id);
        stid=s.str();
        std::pair<std::string,std::string > pair_map;
        pair_map.first = name ;
        pair_map.second = stid;
        dictionary_.push_back(pair_map);
   }
   f.close();
   return dictionary_ ;
}


void MakeASummary1(const char* subdetect, bool BL=false){

  std::ifstream fResEff, fResCls;
  std::ofstream rollclean, rollremoved; 
  std::ofstream rollData;
  
  std::vector< std::pair<std::string,std::string> > map;
  map = dictionary(subdetect);
  std::cout << map.size()<< std::endl;
   
  bool match = false;
  vector<std::string> blacklistrollv;
  if (BL){
  std::string blacklistroll;
  double id;
  std::ifstream blacklist_file;
  blacklist_file.open("../data/BlacklistExclusive_2016.txt");
  //blacklist_file.open("../data/BlacklistInclusive_2016_2017.txt");
  while(blacklist_file.good()){
     blacklist_file >> blacklistroll;
    
    if (blacklist_file.eof())break; 
     blacklistrollv.push_back(blacklistroll);
     }
  blacklist_file.close();
  }
  
  Char_t   RollName[38];
  Double_t WorkingPoint;
  Double_t slope50, emax, hv50, chi2, EffWP, clsWP, chi2cls;
  Double_t emaxErr, slopeErr, hv50Err;
  Double_t slope;
  Int_t nevents = 0;
   
/*  TString outputFile;
  if (strcmp(subdetect,"barrel" ) ==0 && BL==true )outputFile  = "/afs/cern.ch/work/r/rreyesal/public/RPC/HVScanSW/ComparisonAge/summary/barrel_summary_Inc2017BlackList.root";
  else if (strcmp(subdetect,"endcap")==0 && BL==true ) outputFile = "/afs/cern.ch/work/r/rreyesal/public/RPC/HVScanSW/ComparisonAge/summary/endcap_summary_Inc2017BlackList.root";
  else if (strcmp(subdetect,"barrel")==0 && BL==false )outputFile  = "/afs/cern.ch/work/r/rreyesal/public/RPC/HVScanSW/ComparisonAge/summary/barrel_summary_Inc2017.root";
  else if (strcmp(subdetect,"endcap")==0 && BL==false )outputFile  = "/afs/cern.ch/work/r/rreyesal/public/RPC/HVScanSW/ComparisonAge/summary/endcap_summary_Inc2017.root";
  else{}; 

  TFile *file = new TFile(outputFile,"RECREATE");
  TTree *filtered = new TTree("filtered","summary");
  TTree *removed = new TTree("removed","summary");
    

  //To save the filtered 
  filtered->Branch("RollName",&RollName,"RollName/C");
  filtered->Branch("WorkingPoint",&WorkingPoint,"WorkingPoint/D");
  filtered->Branch("emax",&emax,"emax/D");
  filtered->Branch("slope",&slope,"slope/D");
  filtered->Branch("hv50",&hv50,"hv50/D");
  filtered->Branch("chi2",&chi2,"chi2/D");
  filtered->Branch("slope50",&slope50,"slope50/D");
  filtered->Branch("EffWP",&EffWP,"EffWP/D");
  filtered->Branch("clsWP",&clsWP,"clsWP/D");
  filtered->Branch("chi2cls",&chi2cls,"chi2cls/D");
  filtered->Branch("emaxErr",&emaxErr,"emaxErr/D");
  filtered->Branch("slopeErr",&slopeErr,"slopeErr/D");
  filtered->Branch("hv50Err",&hv50Err,"hv50Err/D");

   // Only if the Black list is applied, save only the blacklisted rolls 
  if (BL){
  removed->Branch("RollName",&RollName,"RollName/C");
  removed->Branch("WorkingPoint",&WorkingPoint,"WorkingPoint/D");
  removed->Branch("emax",&emax,"emax/D");
  removed->Branch("slope",&slope,"slope/D");
  removed->Branch("hv50",&hv50,"hv50/D");
  removed->Branch("chi2",&chi2,"chi2/D");
  removed->Branch("slope50",&slope50,"slope50/D");
  removed->Branch("EffWP",&EffWP,"EffWP/D");
  removed->Branch("clsWP",&clsWP,"clsWP/D");
  removed->Branch("chi2cls",&chi2cls,"chi2cls/D");
  removed->Branch("emaxErr",&emaxErr,"emaxErr/D");
  removed->Branch("slopeErr",&slopeErr,"slopeErr/D");
  removed->Branch("hv50Err",&hv50Err,"hv50Err/D");
  }else ; 
*/

  std::string name_;
  std::string id_;
  Double_t FitResEff[10];
  Double_t FitResCls[7];
 
  if (BL){
	if (strcmp(subdetect,"barrel")==0)rollremoved.open("../summary/Exclusiveremoved_barrel_rolls_2017on2016.txt");
	else {rollremoved.open("../summary/Exclusive_removed_endcap_rolls_2017on2016.txt");}
	}
  if (strcmp(subdetect,"barrel")==0)rollData.open("../summary/Exclusive_rollDataSummary_barrel_2017on2016.txt");
  else rollData.open("../summary/Exclusive_rollDataSummary_endcap_2017on216.txt");
  for (vector<std::pair<std::string,std::string> >::const_iterator itmap = map.begin() ;itmap != map.end(); itmap++  ){
  	 name_ =itmap->first;
  	 id_ =itmap->second;
         fResEff.open(("../results/"+name_+"/fitData.txt").c_str());
         while (1){
          fResEff >>FitResEff[0]//wp
                  >>FitResEff[1]//slope50
                  >>FitResEff[2]//emax
                  >>FitResEff[3]//hv50
                  >>FitResEff[4]//chi2
                  >>FitResEff[5]//effwp
                  >>FitResEff[6]//emaxerr
                  >>FitResEff[7]//slopeerr
                  >>FitResEff[8]//hv50err
                  >>FitResEff[9];//slope
                  if (fResEff.eof())break;
                  }
          fResEff.close();
          fResCls.open(("../results/"+name_+"/fitDataCls.txt").c_str());
          while (1){
          fResCls >>FitResCls[0]//a
                  >>FitResCls[1]//b
                  >>FitResCls[2]//c
                  >>FitResCls[3]//d
                  >>FitResCls[4]//chi2cls
                  >>FitResCls[5]//wp
                  >>FitResCls[6];//clswp
                  if (fResEff.eof())break;
          }
          fResCls.close();

          if (BL){
		  match=false;
                  for (int i=0; i <int(blacklistrollv.size()); i++)if(name_ == blacklistrollv.at(i))match=true; 
		  }
          if (match){
		  strcpy(RollName, name_.c_str());
		  WorkingPoint 	= FitResEff[0]; 
		  slope50 	= FitResEff[1]; 
		  emax 		= FitResEff[2]; 
		  hv50 		= FitResEff[3]; 
		  chi2 		= FitResEff[4]; 
		  EffWP 	= FitResEff[5]; 
		  emaxErr 	= FitResEff[6]; 
		  slopeErr 	= FitResEff[7];
		  hv50Err 	= FitResEff[8];
		  slope 	= FitResEff[9];
		  clsWP 	= FitResCls[6];
		  chi2cls	= FitResCls[4]; 
//		  removed->Fill();	
		  rollremoved << RollName <<  " " 
			 <<   id_         << " "  
			 <<  FitResEff[0] << " "  
 			 <<  FitResEff[1] << " "  
 			 <<  FitResEff[2] << " "  
			 <<  FitResEff[3] << " "  
			 <<  FitResEff[4] << " "  
			 <<  FitResEff[5] << " "  
			 <<  FitResEff[6] << " "  
			 <<  FitResEff[7] << " " 
			 <<  FitResEff[8] << " " 
			 <<  FitResEff[9] << " " 
		  	 <<  FitResCls[6] << " " 
			 <<  FitResCls[4] << std::endl; 
  	}
	  else {
		  
	          if (BL) rollclean << RollName << std::endl; 
		  strcpy(RollName, name_.c_str());
		  WorkingPoint 	= FitResEff[0]; 
		  slope50 	= FitResEff[1]; 
		  emax 		= FitResEff[2]; 
		  hv50 		= FitResEff[3]; 
		  chi2 		= FitResEff[4]; 
		  EffWP 	= FitResEff[5]; 
		  emaxErr 	= FitResEff[6]; 
		  slopeErr 	= FitResEff[7];
		  hv50Err 	= FitResEff[8];
		  slope 	= FitResEff[9];
		  clsWP 	= FitResCls[6];
		  chi2cls	= FitResCls[4]; 
//		  filtered->Fill();
                  
	          rollData << RollName <<  " " 
			 <<   id_         << " "  
			 <<  FitResEff[0] << " "  
 			 <<  FitResEff[1] << " "  
 			 <<  FitResEff[2] << " "  
			 <<  FitResEff[3] << " "  
			 <<  FitResEff[4] << " "  
			 <<  FitResEff[5] << " "  
			 <<  FitResEff[6] << " "  
			 <<  FitResEff[7] << " " 
			 <<  FitResEff[8] << " " 
			 <<  FitResEff[9] << " " 
		  	 <<  FitResCls[6] << " " 
			 <<  FitResCls[4] << std::endl; 
			  nevents++; 
		  } 
	}
       std::cout << "rolls->  " <<  nevents << std::endl; 
	if (BL)rollremoved.close();
	rollData.close();
//	filtered->Write("",TObject::kOverwrite);
//	if (BL)removed->Write("",TObject::kOverwrite);
//	file->Write("",TObject::kOverwrite);
//	file->Close();
}





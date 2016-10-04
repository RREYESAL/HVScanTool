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
  if (subd =="endcap")f.open("../data/detIdEndCap.txt");
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


void MakeASummary(const char* subdetect, bool BL=false){

  std::ifstream fResEff, fResCls;
  std::ofstream rollclean;
  std::vector< std::pair<std::string,std::string> > map;
  map = dictionary(subdetect);
  std::cout << map.size()<< std::endl;
   
  bool match = false;
  vector<std::string> blacklistrollv;
  if (BL){
  std::string blacklistroll;
  std::ifstream blacklist_file;
  blacklist_file.open("../data/blacklist.txt");
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
   
  TString outputFile;
  if (subdetect=="barrel" && BL==true )outputFile  = "../summary/barrel_summary_2016BlackList.root";
  else if (subdetect=="endcap" && BL==true ) outputFile = "../summary/endcap_summary_2016BlackList.root";
  else if (subdetect=="barrel" && BL==false )outputFile  = "../summary/barrel_summary_2016.root";
  else if (subdetect=="endcap" && BL==false )outputFile  = "../summary/endcap_summary_2016.root";
  else{}; 

  Int_t nevents = 0;
  TFile *file = new TFile(outputFile,"RECREATE");
  TTree *T = new TTree("T","summary");

  T->Branch("RollName",&RollName,"RollName/C");
  T->Branch("WorkingPoint",&WorkingPoint,"WorkingPoint/D");
  T->Branch("emax",&emax,"emax/D");
  T->Branch("slope",&slope,"slope/D");
  T->Branch("hv50",&hv50,"hv50/D");
  T->Branch("chi2",&chi2,"chi2/D");
  T->Branch("slope50",&slope50,"slope50/D");
  T->Branch("EffWP",&EffWP,"EffWP/D");
  T->Branch("clsWP",&clsWP,"clsWP/D");
  T->Branch("chi2cls",&chi2cls,"chi2cls/D");
  T->Branch("emaxErr",&emaxErr,"emaxErr/D");
  T->Branch("slopeErr",&slopeErr,"slopeErr/D");
  T->Branch("hv50Err",&hv50Err,"hv50Err/D");

  std::string id_;
  Double_t FitResEff[10];
  Double_t FitResCls[7];
 
  if (BL)rollclean.open("../data/cleanedrolls.txt");
  for (vector<std::pair<std::string,std::string> >::const_iterator itmap = map.begin() ;itmap != map.end(); itmap++  ){
  	 id_ =itmap->first;
         fResEff.open(("../results/"+id_+"/fitData.txt").c_str());
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
          fResCls.open(("../results/"+id_+"/fitDataCls.txt").c_str());
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
                  for (int i=0; i <int(blacklistrollv.size()); i++)if(id_ == blacklistrollv.at(i))match=true; 
		  }
          if (match)continue;
	  else {
		  
	          if (BL) rollclean << RollName << std::endl; 
		  strcpy(RollName, id_.c_str());
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
		  nevents++; 
		  T->Fill();
		  } 
	}
        std::cout << "rolls->  " <<  nevents << std::endl; 
	if (BL)rollclean.close();
	T->Write("",TObject::kOverwrite);
	file->Write("",TObject::kOverwrite);
	file->Close();
}





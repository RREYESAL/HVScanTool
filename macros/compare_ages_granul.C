
#include <vector>
#include "TCanvas.h"
#include "TAxis.h"
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <TChain.h>

#include "Data.C"

using namespace TMath;
using namespace std;


//# ifndef __CINT__
//int main()
void compare_ages_granul(TString& eraA, TString &eraB, TString& module)
{
    std::vector<int> samples;
    if (module == "barrel"){
    samples.push_back(1);
    samples.push_back(2);
     }   
    else if (module=="endcap"){
    samples.push_back(-1);
    samples.push_back(-2);
    }
    else cout<<"choose barrel or endcap";       
    
    // -------- D E F I N E   T H E    H I S T O G R A M S  
    Double_t MaxDiffEff=10., MinDiffEff=-10.;
    Double_t MaxDiffWp=80., MinDiffWp=-80.;
    
    TH2F *effwpZRB= new TH2F("effwpZRB"," ",5,0,5,4,0,4); 
    effwpZRB->SetTitle("#Delta EffatWP("+eraB+"-"+eraA+") Z/R ");
    TH2F *effwpRPhiB= new TH2F("effwpRPhiB"," ",12,0,12,10,0,10); 
    effwpRPhiB->SetTitle("#Delta EffatWP("+eraB+"-"+eraA+") #varphi/R ");   
    TH2F *effwpZRlB= new TH2F("effwpZRlB"," ",5,0,5,10,0,10); 
    effwpZRlB->SetTitle("#Delta effwp("+eraB+"-"+eraA+") Z/R ");
    TH2F *effwpZPhiB= new TH2F("effwpZPhiB"," ",12,0,12,5,0,5); 
    effwpZPhiB->SetTitle("#Delta Effatwp("+eraB+"-"+eraA+") Z/#varphi ");
    
   
    TH2F *effwpRPhiEC= new TH2F("effwpRPhiEC"," ",36,0,36,6,0,6); 
    effwpRPhiEC->SetTitle("#Delta EffatWP("+eraB+"-"+eraA+") #varphi/R ");
    TH2F *effwpZREC= new TH2F("effwpZREC"," ",6,-3.,3.,2,0,2); 
    effwpZREC->SetTitle("#Delta EffatWP("+eraB+"-"+eraA+") Z/R ");
 
    TH2F *effwpZPhiEC= new TH2F("effwpZPhiEC"," ",36,0,36,6,0,6);
    effwpZPhiEC->SetTitle("#Delta EffatWP("+eraB+"-"+eraA+") Z/#varphi ");

   
    std::vector<double> wpBARa,  wpBARb, wpENDa,  wpENDb; 
    std::vector<double> hv50BARa, hv50BARb, hv50ENDa, hv50ENDb ; 
    std::vector<double> ewBARa, ewBARb, ewENDa, ewENDb; 
    vector<string> RollBARa, RollENDa; 
    vector<string> RollBARb, RollENDb; 
  
  
   //+++++++++++ D e f i n e   a n d   c r e a t e   y o u r   m a p s   
    vector <string> sector, layer; 
    vector <string> channel, eta_part; 
    vector<string> disk, ring;
    vector<string> wheel, stations;
    
    //++++ E n d   C a p

    for(int i=1;i<=36;i++){
         ostringstream s;
	 string name;
	 s << i; 
	 if (i<10) name = "0"+s.str();
         else      name = s.str();
	channel.push_back(name);
    }
    eta_part.push_back("C");eta_part.push_back("B");eta_part.push_back("A");
    std::vector< std::pair<std::string,std::string >  > EndCap_CH;
    for(int i=0;i<channel.size();i++){
	    for(int j=0;j<eta_part.size();j++){
		    std::pair<std::string,std::string > pair_map;
		    pair_map.first  = channel.at(i);
	            pair_map.second = eta_part.at(j);
		    EndCap_CH.push_back(pair_map);
	    }
    }
 

    std::vector< std::pair<std::string,std::string >  > EndCap_Map; 
    for(int k= -3;k <= 3; k++){
	    if (k!=0){
		    ostringstream s;  string st;
	            s << k; 
		    st=s.str();
		    if (k<0)disk.push_back("RE"+st);
	            else disk.push_back("RE+"+st);
	    }
    }
    ring.push_back("R2"); ring.push_back("R3");
    for (int i=0; i < disk.size(); i++){
      for (int j=0; j <ring.size(); j++){
        std::pair<std::string,std::string > pair_map;
        pair_map.first = disk.at(i);
        pair_map.second = ring.at(j);
        EndCap_Map.push_back(pair_map);  
      }
    }
    std::vector< std::pair<std::string,std::string >  > EndCap_ZPhi;
    for (int i=0; i < channel.size(); i++){
      for (int j=0; j <disk.size(); j++){
        std::pair<std::string,std::string > pair_map;
        pair_map.first = channel.at(i);
        pair_map.second = disk.at(j);
        EndCap_ZPhi.push_back(pair_map);
      }
    }


/*   for (vector <std::pair<std::string,std::string>>::const_iterator it = EndCap_CH.begin() ;it != EndCap_CH.end(); it++  )
	   cout << it->first << " " << it->second << endl;

return;
*/
    //+++++ B a r r e l 
    //
    wheel.push_back("W-2");
    wheel.push_back("W-1");
    wheel.push_back("W+0");
    wheel.push_back("W+1");
    wheel.push_back("W+2");
    stations.push_back("RB1");
    stations.push_back("RB2");
    stations.push_back("RB3");
    stations.push_back("RB4");
    layer.push_back("RB1in");
    layer.push_back("RB1out");
    layer.push_back("RB2in");
    layer.push_back("RB2out");
    layer.push_back("RB3-");
    layer.push_back("RB3+");
    layer.push_back("RB4-");
    layer.push_back("RB4+");
    layer.push_back("RB4--");
    layer.push_back("RB4++");

    for(int i=1;i<=12;i++){
         ostringstream s;
	 string name;
	 s << i; 
	 if (i<10) name = "S0"+s.str();
         else      name = "S"+s.str();
	 sector.push_back(name);
    }

    std::vector< std::pair<std::string,std::string >  > Barrel_S;
    for(int i=0;i<sector.size();i++){
	    for(int j=0;j<layer.size();j++){
		    std::pair<std::string,std::string > pair_map;
		    pair_map.first  = sector.at(i);
	            pair_map.second = layer.at(j);
		    Barrel_S.push_back(pair_map);
	    }
    }
   
    std::vector< std::pair<std::string,std::string >  > Barrel_Map;
    for (int i=0; i < wheel.size(); i++){
      for (int j=0; j <stations.size(); j++){
        std::pair<std::string,std::string > pair_map;
        pair_map.first = wheel.at(i);
        pair_map.second = stations.at(j);
        Barrel_Map.push_back(pair_map);  
      }
    }
    std::vector< std::pair<std::string,std::string >  > Barrel_lW;
    for (int i=0; i < wheel.size(); i++){
      for (int j=0; j <layer.size(); j++){
        std::pair<std::string,std::string > pair_map;
        pair_map.first = wheel.at(i);
        pair_map.second = layer.at(j);
        Barrel_lW.push_back(pair_map);  
      }
    }
    std::vector< std::pair<std::string,std::string >  > Barrel_WS;
    for (int i=0; i < sector.size(); i++){
      for (int j=0; j <wheel.size(); j++){
        std::pair<std::string,std::string > pair_map;
        pair_map.first = sector.at(i);
        pair_map.second = wheel.at(j);
        Barrel_WS.push_back(pair_map);  
      }
    }
    
     //--- G e t   t h e   d a t a -------------------   
    
    gROOT->ProcessLine(".L Data.C+");
    for (std::vector<int>::const_iterator Sample = samples.begin(); Sample != samples.end(); Sample++) {
        //std::cout << (*sample) << std::endl;
        TChain *ch = new TChain("T");
   
        if ( (*Sample)==-1 ) { ch->Add("endcap_summary_"+eraA+".root/T"); }
        if ( (*Sample)==-2 ) { ch->Add("endcap_summary_"+eraB+".root/T"); }
        
        if ( (*Sample)==1 ) { ch->Add("barrel_summary_"+eraA+".root/T"); }
        if ( (*Sample)==2 ) { ch->Add("barrel_summary_"+eraB+".root/T"); }

        TTree *tree = (TTree*) ch;
        Data t(tree);
        Long64_t nentries = t.fChain->GetEntries();
        cout <<" Entries : " << nentries <<endl;
        //Int_t nTen = nentries;
        Int_t nbytes = 0;
        Int_t nb = 0;
        for(Long64_t jentry=0; jentry<nentries;jentry++)
        {
            Long64_t ientry = t.LoadTree(jentry);
            if (ientry < 0)break;
            nb = t.fChain->GetEntry(jentry); nbytes += nb;
            //if(jentry%nTen==0) cout<<10*(jentry/nTen)<<"%-"<<flush;
            if(jentry==nentries-1) cout<<endl;
            
            if ( (*Sample)==-1 ) {
        			   RollENDa.push_back(t.RollName);
                                   wpENDa.push_back(1000*t.WorkingPoint);
                                   hv50ENDa.push_back(1000*t.HVfifty);  
                                   ewENDa.push_back(t.EffWP);  
         			 }

            if ( (*Sample)==-2 ) { 
 				   RollENDb.push_back(t.RollName);
                                   wpENDb.push_back(1000*t.WorkingPoint);
                                   hv50ENDb.push_back(1000*t.HVfifty);  
                                   ewENDb.push_back(t.EffWP);  
                                  } 
            if ( (*Sample)==1 ) { 
 				  RollBARa.push_back(t.RollName);
                                  wpBARa.push_back(1000*t.WorkingPoint);
				  hv50BARa.push_back(1000*t.HVfifty);
                                  ewBARa.push_back(t.EffWP);  
                                 }
            if ( (*Sample)==2 ) { 
                                  RollBARb.push_back(t.RollName);
                                  wpBARb.push_back(1000*t.WorkingPoint);
				  hv50BARb.push_back(1000*t.HVfifty);
                                  ewBARb.push_back(t.EffWP);  
                                  }
        }
    }
    
   //--- Do the operations whatever you want roll by roll
   //
   //
 Int_t y=0;
 nentries = 0;
  Int_t bin=1;
   if(module=="barrel"){ 
   for (vector <std::pair<std::string,std::string>>::const_iterator it = Barrel_Map.begin() ;it != Barrel_Map.end(); it++  ){
       TString nameEW = "hEffWP_"+it->first+"_"+it->second;   
       TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
       for (int i = 0; i<RollBARb.size(); i++) {
         for (int j = 0; j<RollBARa.size(); j++){   
           	if (RollBARb.at(i)!=RollBARa.at(j))continue; 
           	if((RollBARa.at(j)).find(it->first) != std::string::npos )
           		if((RollBARa.at(j)).find(it->second) != std::string::npos ){
            		h1->Fill(ewBARb.at(i) - ewBARa.at(j));    
            		}     
		}
      }
   y++;
   effwpZRB->SetBinContent(bin,y,h1->GetMean()); 
   if (y == 4 ){y=0;   bin++; }  
   }

   y=0;
   bin=1;  
   for (vector <std::pair<std::string,std::string>>::const_iterator it = Barrel_S.begin() ;it != Barrel_S.end(); it++  ){
       TString nameEW = "hEffWP_"+it->first+"_"+it->second;   
       TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
       for (int i = 0; i<RollBARb.size(); i++) {
         for (int j = 0; j<RollBARa.size(); j++){   
           if (RollBARb.at(i)!=RollBARa.at(j))continue; 
           if((RollBARa.at(j)).find(it->first) != std::string::npos )
           if((RollBARa.at(j)).find(it->second) != std::string::npos ){
           h1->Fill(ewBARb.at(i) - ewBARa.at(j));    
          }
         }
      }
   y++;
  
   effwpRPhiB->SetBinContent(bin,y,h1->GetMean()); 
   if( (it->second == "RB4++" || it->second == "RB4--" ) && !(it->first=="S04") ){
   //cout << it->first<< " " << it->second << endl; 
   effwpRPhiB->SetBinContent(bin,y,MinDiffEff-100); 
   }
   if( (it->second == "RB4+") && (it->first=="S09" || it->first=="S11") ){
 //  cout << it->first<< " " << it->second << endl; 
     effwpRPhiB->SetBinContent(bin,y,MinDiffEff-100); 
   } 
   if (y == 10 ){y=0;bin++;}  
  }

  y=0;
  bin=1;  
  for (vector <std::pair<std::string,std::string>>::const_iterator it = Barrel_lW.begin() ;it != Barrel_lW.end(); it++  ){
	  TString nameEW = "hEffWP_"+it->first+"_"+it->second;   
          TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
          for (int i = 0; i<RollBARb.size(); i++) {
		  for (int j = 0; j<RollBARa.size(); j++){   
			  if (RollBARb.at(i)!=RollBARa.at(j))continue; 
           		  if((RollBARa.at(j)).find(it->first) != std::string::npos )
           		  if((RollBARa.at(j)).find(it->second) != std::string::npos ){
				  h1->Fill(ewBARb.at(i) - ewBARa.at(j));    
         	//		  cout << it->first<< " " << it->second <<" "<<RollBARa.at(j)<<" "<< ewBARb.at(i)  <<" "<< ewBARa.at(j)<<" " <<ewBARb.at(i) - ewBARa.at(j)<<  endl; 
            		  }     
		  }
          } 
	  y++;
          //cout << it->first<< " " << it->second << endl; 
	  effwpZRlB->SetBinContent(bin,y,h1->GetMean()); 
	  if (y == 10 ){y=0; bin++;}  
   }
  y=0;
  bin=1;  
  for (vector <std::pair<std::string,std::string>>::const_iterator it = Barrel_WS.begin() ;it != Barrel_WS.end(); it++  ){
	  TString nameEW = "hEffWP_"+it->first+"_"+it->second;   
          TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
          for (int i = 0; i<RollBARb.size(); i++) {
		  for (int j = 0; j<RollBARa.size(); j++){   
			  if (RollBARb.at(i)!=RollBARa.at(j))continue; 
           		  if((RollBARa.at(j)).find(it->first) != std::string::npos )
           		  if((RollBARa.at(j)).find(it->second) != std::string::npos ){
				  h1->Fill(ewBARb.at(i) - ewBARa.at(j));    
            		  }     
		  }
          } 
	  y++;
       nentries += h1->GetEntries(); 
        //  cout << it->first<< " " << it->second << endl; 
	  effwpZPhiB->SetBinContent(bin,y,h1->GetMean()); 
	  if (y == 5 ){y=0; bin++;}  
   }
  
   cout<< "entries" << nentries << endl;


  
 }
//return;
 
 else {  
	 y=0;
	 bin=1;
         nentries=0;  
	 for (vector <std::pair<std::string,std::string>>::const_iterator it = EndCap_Map.begin() ;it != EndCap_Map.end(); it++  ){
         TString nameEW = "hEffWP_"+it->first+"_"+it->second;
         TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
         for (int i = 0; i<RollENDb.size(); i++) {
		  for (int j = 0; j<RollENDa.size(); j++){
			  if (RollENDb.at(i) != RollENDa.at(j))continue;   
	      		  if((RollENDa.at(j)).find(it->first) != std::string::npos )
				  if((RollENDa.at(j)).find(it->second) != std::string::npos ){
					  h1->Fill(ewENDb.at(i) - ewENDa.at(j));
					  }
		}
          }
	  y++; 
          nentries += h1->GetEntries(); 
   	  effwpZREC->SetBinContent(bin,y,h1->GetMean()); 
   	  if (y == 2 ){y=0;bin++;}
       }

       y=0;bin=1;
       for (vector <std::pair<std::string,std::string>>::const_iterator it = EndCap_CH.begin() ;it != EndCap_CH.end(); it++  ){
       TString name2 = "hWP_R2"+it->first+"_"+it->second;
       TH1F *h1 = new TH1F(name2," ",100,-100,100);
       TString name3 = "hEffWP_R3_"+it->first+"_"+it->second;
       TH1F *h2 = new TH1F(name3," ",100,-100,100);
      
        for (int i = 0; i<RollENDb.size(); i++) {
          for (int j = 0; j<RollENDa.size(); j++){
             if (RollENDb.at(i) != RollENDa.at(j))continue;   
	      if((RollENDa.at(j)).find("R2_CH"+it->first+"_"+it->second) != std::string::npos ){
                   h1->Fill(ewENDb.at(i) - ewENDa.at(j));
  //    	   cout << it->first<< " " << it->second <<" "<<RollENDa.at(j)<<" "<< ewENDb.at(i)  <<" "<< ewENDa.at(j)<<" " <<ewENDb.at(i) - ewENDa.at(j)<<  endl; 
	      }
	      if((RollENDa.at(j)).find("R3_CH"+it->first+"_"+it->second) != std::string::npos ){
                   h2->Fill(ewENDb.at(i) - ewENDa.at(j));
    //               cout << it->first<< " " << it->second <<" "<<RollENDa.at(j)<<" "<< ewENDb.at(i)  <<" "<< ewENDa.at(j)<<" " <<ewENDb.at(i) - ewENDa.at(j)<<  endl; 
	       }

	   }
       
      }
      y++; 
     // cout << it->first<< " " << it->second << endl; 
      effwpRPhiEC->SetBinContent(bin,y,h1->GetMean()); 
      effwpRPhiEC->SetBinContent(bin,y+3,h2->GetMean()); 
      if (y == 3 ){y=0;bin++;}
      }
	 
      y=0;
	 bin=1;  
	 for (vector <std::pair<std::string,std::string>>::const_iterator it = EndCap_ZPhi.begin() ;it != EndCap_ZPhi.end(); it++  ){
         TString nameEW = "hEffWP_"+it->first+"_"+it->second;
         TH1F *h1 = new TH1F(nameEW," ",100,-100,100);
         for (int i = 0; i<RollENDb.size(); i++) {
		  for (int j = 0; j<RollENDa.size(); j++){
			  if (RollENDb.at(i) != RollENDa.at(j))continue;   
	      		  if((RollENDa.at(j)).find(it->first) != std::string::npos )
				  if((RollENDa.at(j)).find(it->second) != std::string::npos ){
					  h1->Fill(ewENDb.at(i) - ewENDa.at(j));
					  }
		}
          }
	  y++; 
          //cout << it->first<< " " << it->second << endl;    
	  effwpZPhiEC->SetBinContent(bin,y,h1->GetMean()); 
   	  if (y == 6 ){y=0;bin++;}
       }
     
       
          cout<< "entries" << nentries << endl;
  
 }

   //========== P r i n t i n g   t h e   p l o t s
  
  if(module=="barrel"){
   
   cout<< wheel.size()<<endl; 
  //---- Put the label in axis   
    for (i=0;i<wheel.size();i++){ 
	   char* label= wheel.at(i).c_str();
	   effwpZRB->GetXaxis()->SetBinLabel(i+1,label);
           effwpZRlB->GetXaxis()->SetBinLabel(i+1,label);
           effwpZPhiB->GetYaxis()->SetBinLabel(i+1,label);

        }

    for (i=0;i<stations.size();i++){ 
	    const char* label= stations.at(i).c_str();
      	    effwpZRB->GetYaxis()->SetBinLabel(i+1,label);
     }
    for (i=0;i<sector.size();i++){
	     const char* label= sector.at(i).c_str();
	     effwpRPhiB->GetXaxis()->SetBinLabel(i+1,label);
	     effwpZPhiB->GetXaxis()->SetBinLabel(i+1,label);
     }
    for (i=0;i<layer.size();i++){ 
      const char* label= layer.at(i).c_str();
      effwpRPhiB->GetYaxis()->SetBinLabel(i+1,label);
      effwpZRlB->GetYaxis()->SetBinLabel(i+1,label);
     }

  }
  
  else {
  //---- Put the label in axis   
    for (i=0;i<disk.size();i++){ 
      const char* label= disk.at(i).c_str();
      effwpZREC->GetXaxis()->SetBinLabel(i+1,label);
      effwpZPhiEC->GetYaxis()->SetBinLabel(i+1,label);
     }
    for (i=0;i<ring.size();i++){ 
      const char* label=  ring.at(i).c_str();
      effwpZREC->GetYaxis()->SetBinLabel(i+1,label);
     }
    for (i=0;i<channel.size();i++){ 
      const  char* label= channel.at(i).c_str();
      effwpRPhiEC->GetXaxis()->SetBinLabel(i+1,label);
      effwpZPhiEC->GetXaxis()->SetBinLabel(i+1,label);
     }
    for (i=0;i<eta_part.size();i++){ 
      const char* label= eta_part.at(i).c_str();
      if(i==1)const char*  label =("Ring2 "+eta_part.at(i)).c_str(); 
      effwpRPhiEC->GetYaxis()->SetBinLabel(i+1,label);
      if(i==1)const char*  label =("Ring3 "+eta_part.at(i)).c_str(); 
      effwpRPhiEC->GetYaxis()->SetBinLabel(i+4,label);
      }
  }

 
  // Draw the Map in a histogram         

  gStyle->SetOptStat(0); 
  gStyle->SetPalette(1); 
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1); 
  gStyle->SetPalette(1); 
  TCanvas* c = new TCanvas("c", "myCanvas",1000,600);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetRightMargin(0.155);
   c->SetFrameBorderMode(0);
   c->SetFrameBorderMode(0);
   c->SetGrid();
  TCanvas* c1 = new TCanvas("c1", " ",800,800);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetRightMargin(0.136);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   c1->SetGrid();
  TCanvas* c2 = new TCanvas("c2", "myCanvas",1000,600);
   c2->SetFillColor(0);
   c2->SetBorderMode(0);
   c2->SetBorderSize(2);
   c2->SetRightMargin(0.155);
   c2->SetFrameBorderMode(0);
   c2->SetFrameBorderMode(0);
   c2->SetGrid();
  TCanvas* c3 = new TCanvas("c3", " ",1000,600);
   c3->SetFillColor(0);
   c3->SetBorderMode(0);
   c3->SetBorderSize(2);
   c3->SetRightMargin(0.136);
   c3->SetFrameBorderMode(0);
   c3->SetFrameBorderMode(0);
   c3->SetGrid();


  TString axis1="[#varphi]";
  TString axis2="[R]";
  TString axis3="[Z]";
  TString axis4="[%]";
  
  if(module=="barrel"){
  
  effwpRPhiB->GetZaxis()->SetTitle(axis4);
  effwpRPhiB->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpRPhiB->GetXaxis()->SetTitle(axis1);
  effwpRPhiB->GetYaxis()->SetTitle(axis2);
  effwpRPhiB->GetXaxis()->SetLabelSize(.05);
  effwpRPhiB->GetYaxis()->SetLabelSize(.05);  
  effwpRPhiB->GetYaxis()->SetTitleOffset(1.3);  
  effwpRPhiB->GetXaxis()->SetTitleOffset(1.2);  
  effwpRPhiB->GetZaxis()->SetTitleOffset(1.);  
  
  effwpZRB->GetZaxis()->SetTitle(axis4);
  effwpZRB->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpZRB->GetXaxis()->SetTitle(axis3);
  effwpZRB->GetYaxis()->SetTitle(axis2);
  effwpZRB->GetXaxis()->SetLabelSize(.05);
  effwpZRB->GetYaxis()->SetLabelSize(.05);  
  effwpZRB->GetYaxis()->SetTitleOffset(1.3);  
  effwpZRB->GetXaxis()->SetTitleOffset(1.2);  
  effwpZRB->GetZaxis()->SetTitleOffset(1.);  
  
  c->cd();
  effwpRPhiB->Draw("COLZ"); 
  c->SaveAs("test2/RPhi_"+eraB+"_"+eraA+"_barrel.png");
  c->Update(); 
  c->Modified(); 
 
  c1->cd();
  effwpZRB->Draw("COLZ");
  c1->SaveAs("test2/ZR_"+eraB+"_"+eraA+"_barrel.png");
  c1->Update(); 
  c1->Modified(); 
  
  effwpZPhiB->GetZaxis()->SetTitle("[%]");
  effwpZPhiB->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpZPhiB->GetXaxis()->SetTitle(axis1);
  effwpZPhiB->GetYaxis()->SetTitle(axis3);
  effwpZPhiB->GetXaxis()->SetLabelSize(.05);
  effwpZPhiB->GetYaxis()->SetLabelSize(.05);  
  effwpZPhiB->GetYaxis()->SetTitleOffset(1.3);  
  effwpZPhiB->GetXaxis()->SetTitleOffset(1.2);  
  effwpZPhiB->GetZaxis()->SetTitleOffset(1.);  
  
  effwpZRlB->GetZaxis()->SetTitle("[%]");
  effwpZRlB->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpZRlB->GetXaxis()->SetTitle(axis3);
  effwpZRlB->GetYaxis()->SetTitle(axis2);
  effwpZRlB->GetXaxis()->SetLabelSize(.05);
  effwpZRlB->GetYaxis()->SetLabelSize(.05);  
  effwpZRlB->GetYaxis()->SetTitleOffset(1.3);  
  effwpZRlB->GetXaxis()->SetTitleOffset(1.2);  
  effwpZRlB->GetZaxis()->SetTitleOffset(1.);  
  
  c2->cd();
  effwpZPhiB->Draw("COLZ"); 
  c2->SaveAs("test2/ZPhi_"+eraB+"_"+eraA+"_barrel.png");
  c2->Update(); 
  c2->Modified(); 
 
  c3->cd();
  effwpZRlB->Draw("COLZ");
  c3->SaveAs("test2/ZRl_"+eraB+"_"+eraA+"_barrel.png");
  c3->Update(); 
  c3->Modified(); 
  
  }
 else {
 
  effwpRPhiEC->GetZaxis()->SetTitle("[%]");
  effwpRPhiEC->GetXaxis()->SetTitle(axis1);
  effwpRPhiEC->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpRPhiEC->GetYaxis()->SetTitle(axis2);
  effwpRPhiEC->GetXaxis()->SetLabelSize(.04);
  effwpRPhiEC->GetYaxis()->SetLabelSize(.04);  
  c->cd();
  effwpRPhiEC->Draw("COLZ");
  c->SaveAs("test2/RPhi_"+eraB+"_"+eraA+"_EndCap.png");
  c->Update(); 
  c->Modified(); 
  effwpZREC->GetZaxis()->SetTitle("[%]");
  effwpZREC->GetXaxis()->SetTitle(axis3);
  effwpZREC->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpZREC->GetYaxis()->SetTitle(axis2);
  effwpZREC->GetXaxis()->SetLabelSize(.05);
  effwpZREC->GetYaxis()->SetLabelSize(.05);  
  c1->cd();
  effwpZREC->Draw("COLZ");
  c1->SaveAs("test2/ZR_"+eraB+"_"+eraA+"_EndCap.png");
  c1->Update(); 
  c1->Modified();
  effwpZPhiEC->GetZaxis()->SetTitle("[%]");
  effwpZPhiEC->GetZaxis()->SetRangeUser(MinDiffEff,MaxDiffEff);
  effwpZPhiEC->GetXaxis()->SetTitle(axis1);
  effwpZPhiEC->GetYaxis()->SetTitle(axis3);
  effwpZPhiEC->GetXaxis()->SetLabelSize(.05);
  effwpZPhiEC->GetYaxis()->SetLabelSize(.05);
  effwpZPhiEC->GetYaxis()->SetTitleOffset(1.3);
  effwpZPhiEC->GetXaxis()->SetTitleOffset(1.2);
  effwpZPhiEC->GetZaxis()->SetTitleOffset(1.);
  c2->cd();
  effwpZPhiEC->Draw("COLZ");
  c2->SaveAs("test2/ZPhi_"+eraB+"_"+eraA+"_EndCap.png");
  c2->Update();
  c2->Modified();
 }
 wheel.clear();
 layer.clear();
 sector.clear();
 stations.clear();
 disk.clear();
 ring.clear();
 eta_part.clear();
 EndCap_CH.clear();  
 EndCap_Map.clear();  
 Barrel_S.clear();  
 Barrel_WS.clear();  
 Barrel_lW.clear();  
 Barrel_Map.clear();  
 
 return;
}


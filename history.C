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
#include <cstdio>
#include <cstdlib>
#include <TChain.h>

#include "Data.C"

using namespace TMath;
using namespace std;


//# ifndef __CINT__
//int main()
void history(TString& module)
{
   
    // Define the arrays to graph 
    const int YEARS = 10;
    const NUMPOINTS = 4;
    const char *year[YEARS]  = {"2011","2012a","2012b","2013a", "2013b","2014a","2014b", "2015a","2015b","2016"};
    
    Double_t x[NUMPOINTS]= {0.5,1.5,2.5,9.5};
    Double_t xerr[NUMPOINTS] = {.0,.0,.0,.0};
    
    Double_t yEff[NUMPOINTS]; 
    Double_t yEfferr[NUMPOINTS]; 
    Double_t yWP[NUMPOINTS];
    Double_t yWPerr[NUMPOINTS]; 
   
    // Define the histograms to xtract information  
    TH1F* heffwp1 = new TH1F("effwp_2011", "eff wp", 101, 50., 110);
    TH1F* heffwp2I = new TH1F("effwp_2012I", "eff wp", 101, 50, 110);
    TH1F* heffwp2II = new TH1F("effwp_2012II", "eff wp", 101, 50, 110);
    TH1F* heffwp6  = new TH1F("effwp_2016", "eff wp", 101, 50, 110);
    TH1F* hwp1 = new TH1F("hwp1", "wp", 60, 8900, 10300);
    TH1F* hwp2I = new TH1F("hwp2I", "wp ", 60, 8900, 10300);
    TH1F* hwp2II = new TH1F("hwp2II", "wp", 60, 8900, 10300);
    TH1F* hwp6  = new TH1F("hwp6", "wp", 60, 8900, 10300);
    
   //Define the vector to fill the histograms  
    vector<double>  effwp1, effwp2II, effwp2I, effwp5, effwp6, effwp5b;
    vector<double>  wp1, wp2II, wp2I, wp5, wp6, wp5b;
    vector<string> Roll1, Roll2I; 
    vector<string> Roll2II,Roll6; 
 
    // The Data are come from barrel or endcap 
    std::vector<int> samples;
    if (module=="barrel"){
    samples.push_back(1);
    samples.push_back(2);
    samples.push_back(3);
    samples.push_back(4);
     }   
   else if (module=="endcap"){
    samples.push_back(-1);
    samples.push_back(-2);
    samples.push_back(-3);
    samples.push_back(-4);
    }
    
  else cout<<"choose barrel or endcap";       
    
 
    //--- Get the data  

    gROOT->ProcessLine(".L Data.C+");
    
    for (std::vector<int>::const_iterator Sample = samples.begin(); Sample != samples.end(); Sample++) {
        //std::cout << (*sample) << std::endl;
        TChain *ch = new TChain("T");
   
        if ( (*Sample)==-1 ) { ch->Add("endcap_summary_2011.root/T"); }
        if ( (*Sample)==-2 ) { ch->Add("endcap_summary_2012Iv2.root/T"); }
        if ( (*Sample)==-3 ) { ch->Add("endcap_summary_2012IIv2.root/T"); }
        if ( (*Sample)==-4 ) { ch->Add("endcap_summary_2016BLv2.root/T"); }
        
        if ( (*Sample)==1 ) { ch->Add("barrel_summary_2011.root/T"); }
        if ( (*Sample)==2 ) { ch->Add("barrel_summary_2012Iv2.root/T"); }
        if ( (*Sample)==3 ) { ch->Add("barrel_summary_2012IIv2.root/T"); }
        if ( (*Sample)==4 ) { ch->Add("barrel_summary_2016BLv2.root/T"); }


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
     
        // Fill the vectors and histograms        
            if ( (*Sample)==-1 ) { 
 				  Roll1.push_back(t.RollName);
       				  effwp1.push_back(t.EffWP);
       				  wp1.push_back(1000*t.WorkingPoint);
        			  heffwp1->Fill(t.EffWP);
        			  hwp1->Fill(1000*t.WorkingPoint);
				  }
            if ( (*Sample)==-2 ) { 
                                  Roll2I.push_back(t.RollName);
                                  effwp2I.push_back(t.EffWP);
                                  wp2I.push_back(1000*t.WorkingPoint);
                                 // heffwp2I->Fill(t.EffWP);
                                 }
            if ( (*Sample)==-3 ) { 
                                  Roll2II.push_back(t.RollName);
                                  effwp2II.push_back(t.EffWP);
                                  wp2II.push_back(1000*t.WorkingPoint);
 				  heffwp2II->Fill(t.EffWP);
 				  hwp2II->Fill(1000*t.WorkingPoint);
                                 }
            if ( (*Sample)==-4 ) { 
                                  Roll6.push_back(t.RollName);
                                  effwp6.push_back(t.EffWP);//heffwp6->Fill(t.EffWP);
                                  wp6.push_back(1000*t.WorkingPoint);//heffwp6->Fill(t.EffWP);
                                  }

           if ( (*Sample)==1 ) { 
 				Roll1.push_back(t.RollName);
       				effwp1.push_back(t.EffWP);
       				wp1.push_back(1000*t.WorkingPoint);
 				heffwp1->Fill(t.EffWP);
				hwp1->Fill(1000*t.WorkingPoint);
                                 }
            if ( (*Sample)==2 ) { 
                                 Roll2I.push_back(t.RollName);
       				 wp2I.push_back(1000*t.WorkingPoint);
                                 effwp2I.push_back(t.EffWP);//heffwp2I->Fill(t.EffWP);
                                 }
            if ( (*Sample)==3 ) { 
                                 Roll2II.push_back(t.RollName);
                                 effwp2II.push_back(t.EffWP);
                                 wp2II.push_back(1000*t.WorkingPoint);
                                 heffwp2II->Fill(t.EffWP);
                                 hwp2II->Fill(1000*t.WorkingPoint);
                                 }
            if ( (*Sample)==4 ) { 
                                  Roll6.push_back(t.RollName);
                                  effwp6.push_back(t.EffWP);
                                  wp6.push_back(1000*t.WorkingPoint);
				  heffwp6->Fill(t.EffWP);
				  hwp6->Fill(1000*t.WorkingPoint);
                                  }
        }
    }
   //---Filling the histograms  
   if(module=="endcap"){
     for(int i=0; i<Roll1.size();i++){
         for(int j=0; j<Roll6.size();j++){
 	   if(Roll1.at(i) != Roll6.at(j))continue;
              //cout << "2016 "<< Roll6.at(j)<< " 2011 " << Roll1.at(i) << std::endl;
	      heffwp6->Fill(effwp6.at(j));
	      hwp6->Fill(wp6.at(j));
	}
     }
     for(int i=0; i<Roll2II.size();i++){
         for(int j=0; j<Roll2I.size();j++){
 	   if(Roll2II.at(i) != Roll2I.at(j))continue;
              //cout << "2016 "<< Roll6.at(j)<< " 2011 " << Roll1.at(i) << std::endl;
	      heffwp2I->Fill(effwp2I.at(j));
	      hwp2I->Fill(wp2I.at(j));
	}
     }
      
  }
  else {
  
     for(int i=0; i<Roll2II.size();i++){
         for(int j=0; j<Roll2I.size();j++){
 	   if(Roll2II.at(i) != Roll2I.at(j))continue;
              //cout << "2016 "<< Roll6.at(j)<< " 2011 " << Roll1.at(i) << std::endl;
	      heffwp2I->Fill(effwp2I.at(j));
	      hwp2I->Fill(wp2I.at(j));
	}
     }
  
   
  }
 // Filling the arrays 

     *(yEff)= heffwp1->GetMean();
     *(yEfferr)= (heffwp1->GetStdDev())/sqrt(heffwp1->GetEntries());
     *(yEff+1)= heffwp2I->GetMean();
     *(yEfferr+1)= (heffwp2I->GetStdDev())/sqrt(heffwp2I->GetEntries());
     *(yEff+2)= heffwp2II->GetMean();
     *(yEfferr+2)= (heffwp2II->GetStdDev())/sqrt(heffwp2II->GetEntries());
     *(yEff+3)= heffwp6->GetMean();
     *(yEfferr+3)= (heffwp6->GetStdDev())/sqrt(heffwp6->GetEntries());
     
     *(yWP)= hwp1->GetMean();
     *(yWPerr)= (hwp1->GetStdDev())/sqrt(hwp1->GetEntries());
     *(yWP+1)= hwp2I->GetMean();
     *(yWPerr+1)= (hwp2I->GetStdDev())/sqrt(hwp2I->GetEntries());
     *(yWP+2)= hwp2II->GetMean();
     *(yWPerr+2)= (hwp2II->GetStdDev())/sqrt(hwp2II->GetEntries());
     *(yWP+3)= hwp6->GetMean();
     *(yWPerr+3)= (hwp6->GetStdDev())/sqrt(hwp6->GetEntries());
  
  //heffwp2II->Draw();
 
  //heffwp->Set
 // +++ D R A W  T H E  G R A P H S  
  Double_t xmin=0, xmax=YEARS, ymin, ymax; 
  //pad->GetUxmin();
  if (module =="endcap" ){ymin = 9400;  ymax = 10000;}
  else { ymin = 9200; ymax = 9800;}
  //create a TCanvas
  TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
  TPad *pad = new TPad("pad","",0,0,1,1);  
  pad->Draw();
  pad->cd();
  // draw a frame to define the range
  TH1F *hr = pad->DrawFrame(0.0,30,5.,110);
  TH1F *h = new TH1F("h"," ",YEARS,0,YEARS);
  // PUT THE LABEL AXIS 
  for (i=1;i<=YEARS;i++)
	  h->GetXaxis()->SetBinLabel(i,year[i-1]);
	      
  h->SetStats(0);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetTickLength(0.02);
  h->SetMaximum(98);
  h->SetMinimum(90);
  h->Draw("text");
  //  DEFINE THE GRAPHS 
  TGraphErrors *gr1 = new TGraphErrors(NUMPOINTS,x,yEff,xerr,yEfferr);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(21);
  gr1->SetName("gr1");
  gr1->Draw("LP");
  TExec *ex = new TExec("ex","drawtext();");
  gr1->GetListOfFunctions()->Add(ex);

  TGraphErrors *gr2 = new TGraphErrors(NUMPOINTS,x,yWP,xerr,yWPerr);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(20);
  gr2->SetName("gr2");
  TExec *ex1 = new TExec("ex1","drawtextWP();");
  gr2->GetListOfFunctions()->Add(ex1);
  c1->cd();
 
  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000); 
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();
  TH1F *hframe = overlay->DrawFrame(xmin,ymin,xmax,ymax);
  hframe->SetMarkerColor(2);
  hframe->GetYaxis()->SetTitle("eff[%]");
  hframe->GetYaxis()->SetTitleSize(0.04);
  hframe->GetXaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetTitleOffset(1.0);
  hframe->GetYaxis()->SetTickLength(0);
  hframe->GetXaxis()->SetTickLength(0);
  gr2->Draw("LP");
  //Draw an axis on the right side
   TGaxis *axis = new TGaxis(xmax,ymin,xmax, ymax,ymin,ymax,510,"+L");
   axis->SetLineColor(1);
   axis->SetLabelColor(1);
   axis->SetLabelSize(0.035);
   axis->SetTitleSize(0.045);
   axis->SetTitleOffset(1.08);
   axis->SetTickLength(0.02);
   axis->SetTitle("WP [V]");
   axis->Draw();
  
  Leg = new TLegend(0.2, 0.2, 0.3, 0.3);
  Leg->AddEntry(gr1,"Eff at WP", "lpl");
  Leg->AddEntry(gr2, "WP", "lpl");
  //Leg->SetTextAlign(22);
  //Leg->SetTextSize(TextSize);
  Leg->Draw();  
  if(module=="endcap")c1->SaveAs("historyPlot_EffWP_EndCap.C");
  else c1->SaveAs("historyPlot_EffWP_Barrel.C");
}

void drawtext()
{  
   Int_t i,n;
   Double_t x,y;
   TLatex *l;
   
   TH1F *h = (TH1F*)gPad->GetListOfPrimitives()->FindObject("h");
   TGraphErrors *g = (TGraphErrors*)gPad->GetListOfPrimitives()->FindObject("gr1");
   n = g->GetN();
   //cout << n  << endl; 
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      l = new TLatex(x,y+.5,Form("%2.2f",y));
      l->SetTextSize(0.03);
      l->SetTextFont(42);
      l->SetTextColor(4);
      l->SetTextAlign(21);
      l->Paint();
   }
}
void drawtextWP()
{  
   Int_t i,n;
   Double_t x,y;
   TLatex *l;
   
   TH1F *h = (TH1F*)gPad->GetListOfPrimitives()->FindObject("hframe");
   TGraphErrors *g = (TGraphErrors*)gPad->GetListOfPrimitives()->FindObject("gr2");
   n = g->GetN();
   //cout << n  << endl; 
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      l = new TLatex(x,y+12,Form("%4.0f",y));
      l->SetTextSize(0.03);
      l->SetTextFont(42);
      l->SetTextColor(2);
      l->SetTextAlign(21);
      l->Paint();
   }
}


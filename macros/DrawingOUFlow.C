#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
// ROOT 
#include "TPad.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLegend.h"
#include "TText.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TMath.h"
#include "TTree.h"
#include "TROOT.h"
#include "TFile.h"
#include "TBrowser.h"
#include "TH2.h"
#include "TBranch.h"
#include "TChain.h"
#include "TUnixSystem.h"

void Title(TH1F* ,const char* ,const char* ,const char* , float);
void Format(TH1F* , int , int , int , int);
void Legend2(float, float, float, float,const char*, TH1F*,const char*, TH1F*,const char*, float);
void Legend3(float, float, float, float, const char*, TH1F*, const char*, TH1F*, const char*, TH1F*, const  char*, float);
void DrawWithOFUF(TH1F*, bool, int);
void DrawDiffWithOFUF(TH1F* , bool, bool);

using namespace std; 

void DrawingOUFlow()
{
  
  const int NUM = 9; 
  TCanvas* MyCanvas[NUM];  
  //TPad* MyPad[NUM];  
  TH1F* HistosB[NUM];
  TH1F* HistosEC[NUM];
  TH1F* HistosRE4[NUM];
  std::string TitlesHistos[9]; 
  std::string TitlesPNG[9]; 
  // Histos
  HistosB[0] 	= new TH1F("wp_Barrel", "Working Point", 70, 8.6, 10);
  HistosEC[0] 	= new TH1F("wp_Endcap", "Working Point", 70, 8.6, 10);
  HistosRE4[0]	= new TH1F("wp_RE4", "Working Point ", 70, 8.6, 10);
  TitlesHistos[0] = "Working Point [V]";   
  TitlesPNG[0] = "wp";   


  HistosB[1] 	= new TH1F("hv50_Barrel", "HV50", 70, 8.6, 10);
  HistosEC[1] 	= new TH1F("hv50_Endcap", "HV50", 70, 8.6, 10); 
  HistosRE4[1] 	= new TH1F("hv50_RE4", "HV50 ", 70, 8.6, 10);
  TitlesHistos[1] = "HV50 [V]";
  TitlesPNG[1] = "hv50";


  HistosB[2] = new TH1F("knee_Barrel", "Knee", 70, 8.6, 10);
  HistosEC[2] = new TH1F("knee_Endcap", "Knee", 70, 8.6, 10);
  HistosRE4[2] = new TH1F("knee_RE4", "Knee ", 70, 8.6, 10);
  TitlesHistos[2] = "Knee [V]";
  TitlesPNG[2] = "knee";


  HistosB[3] = new TH1F("emax_Barrel", "emax", 25, 70, 120); 
  HistosEC[3] = new TH1F("emax_Endcap", "emax", 25, 70, 120);
  HistosRE4[3] = new TH1F("emax_RE4", "emax ", 25, 70, 120);
  TitlesHistos[3] = "Emax [%]";
  TitlesPNG[3] = "emax";

  HistosB[4] = new TH1F("slope50_Barrel", "slope50", 50, 0, 400);
  HistosEC[4] = new TH1F("slope50_Endcap", "slope50", 50, 0, 400);
  HistosRE4[4] = new TH1F("slope50_RE4", "slope50 ", 50, 0, 400);
  TitlesHistos[4] = "Solpe50 [%]";
  TitlesPNG[4] = "slope50";

  HistosB[5] = new TH1F("Chi2_Barrel", "chi2", 44, -2, 20);
  HistosEC[5] = new TH1F("Chi2_Endcap", "chi2", 44, -2, 20);
  HistosRE4[5] = new TH1F("Chi2_RE4", "chi2 ", 44, -2, 20);
  TitlesHistos[5] = "Chi2";
  TitlesPNG[5] = "chi2";
  
  HistosB[6] = new TH1F("effwp_Barrel", "effwp", 30, 55, 115);
  HistosEC[6] = new TH1F("effwp_Endcap", "effwp", 30, 55, 115);
  HistosRE4[6] = new TH1F("effwp_RE4", "effwp ", 30, 55, 115);
  TitlesHistos[6] = "Efficiency at WP [%]";
  TitlesPNG[6] = "effwp";

  HistosB[7] = new TH1F("cls_Barrel", "Cluster Size", 48, -2, 6);
  HistosEC[7] = new TH1F("cls_Endcap", "Cluster Size", 48, -2, 6);
  HistosRE4[7] = new TH1F("cls_RE4", "Cluster Size", 48, -2, 6);
  TitlesHistos[7] = "Cluster Size";       
  TitlesPNG[7] = "clswp";       


  HistosB[8] = new TH1F("Chi2Cls_Barrel", "effwp", 44, -2, 20);
  HistosEC[8] = new TH1F("Chi2Cls_Endcap", "effwp", 44, -2, 20);
  HistosRE4[8] = new TH1F("Chi2Cls_RE4", "effwp ", 44, -2, 20);
  TitlesHistos[8] = "Chi2 ClusterSize";  
  TitlesPNG[8] = "chi2clswp";  

  vector <double> wp, effmax, HV50, slope;
  vector <double> effwp, chi2eff, clswp, chi2clswp;
  vector <std::string> roll;
  Char_t   RollName[38];
  Double_t WorkingPoint;
  Double_t slope50, emax, hv50, chi2, EffWP, clsWP, chi2cls;
  TChain *ch = new TChain("ch","");
  ch->Add("../summary/barrel_summary_2016.root/T");
  ch->Add("../summary/endcap_summary_2016.root/T");

  TTree *T = (TTree*)ch;

  cout<<" Entries : "<<T->GetEntries()<<endl;
  Long64_t nentries = T->GetEntries();
  T->SetBranchAddress("RollName", &RollName);
  T->SetBranchAddress("WorkingPoint", &WorkingPoint);
  T->SetBranchAddress("emax", &emax);
  T->SetBranchAddress("hv50", &hv50);
  T->SetBranchAddress("chi2", &chi2);
  T->SetBranchAddress("slope50", &slope50);
  T->SetBranchAddress("EffWP", &EffWP);
  T->SetBranchAddress("clsWP", &clsWP);
  T->SetBranchAddress("chi2cls", &chi2cls);

 for(Long64_t i=0;i<nentries;i++){
           T->GetEntry(i);
           std::string chamber = RollName; 
           roll.push_back(RollName);
           wp.push_back(WorkingPoint);
           effmax.push_back(emax);
           HV50.push_back(hv50);
           slope.push_back(slope50);  
           effwp.push_back(EffWP);  
           clswp.push_back(clsWP);  
           chi2eff.push_back(chi2);  
           chi2clswp.push_back(chi2cls);  

         }
  gROOT->Reset(); 
  for (int k=0;k<int(roll.size()); k++){      
     if((roll.at(k)).find("W")!= std::string::npos){
              HistosB[0]->Fill(wp.at(k));
              HistosB[1]->Fill(HV50.at(k));
              HistosB[2]->Fill(wp.at(k) - .100);//Knee 
              HistosB[3]->Fill(effmax.at(k));
              HistosB[4]->Fill(slope.at(k));
              HistosB[5]->Fill(chi2eff.at(k));  
              HistosB[6]->Fill(effwp.at(k));  
              HistosB[7]->Fill(clswp.at(k)); 
              HistosB[8]->Fill(chi2clswp.at(k)); 
            }
            else if ((roll.at(k)).find("RE+4")!= std::string::npos || (roll.at(k)).find("RE-4")!= std::string::npos ){
              HistosRE4[0]->Fill(wp.at(k));
              HistosRE4[1]->Fill(HV50.at(k));
              HistosRE4[2]->Fill(wp.at(k) - .120);//Knee 
              HistosRE4[3]->Fill(effmax.at(k));
              HistosRE4[4]->Fill(slope.at(k));
              HistosRE4[5]->Fill(chi2eff.at(k));  
              HistosRE4[6]->Fill(effwp.at(k));  
              HistosRE4[7]->Fill(clswp.at(k)); 
              HistosRE4[8]->Fill(chi2clswp.at(k)); 
                   }
             else {
              HistosEC[0]->Fill(wp.at(k));
              HistosEC[1]->Fill(HV50.at(k));
              HistosEC[2]->Fill(wp.at(k) - .120);//Knee 
              HistosEC[3]->Fill(effmax.at(k));
              HistosEC[4]->Fill(slope.at(k));
              HistosEC[5]->Fill(chi2eff.at(k));  
              HistosEC[6]->Fill(effwp.at(k));  
              HistosEC[7]->Fill(clswp.at(k)); 
              HistosEC[8]->Fill(chi2clswp.at(k)); 
                 }
           
      }
     
     gStyle->SetOptStat(1111);
      // create your canvas and the pads 
     for(int j=0; j< NUM; j++){
     	std::stringstream ss;
        ss << j;
        std::string canvasstr = "canvas"+ss.str(); 
        std::string padstr = "pad"+ss.str(); 
        const char* name1 = canvasstr.c_str();
        const char* name2 = padstr.c_str();
        MyCanvas[j] = new TCanvas(name1," ", 1024, 768);  
       
        TPad *pad = new TPad(name2," ",0.0, 0.0, 1.0, 1.0, 0);  
        MyCanvas[j]->cd();
        MyCanvas[j]->Draw(); 
        pad->Draw();
        pad->SetTicks(1, 1);
        pad->cd();
 
        Title(HistosB[j], TString(TitlesHistos[j]), TString(TitlesHistos[j]), "Number of Rolls", 1.2);
        Format(HistosB[j], 2, 2, 0, 0);
        DrawWithOFUF(HistosB[j], 0, 0);
       Format(HistosEC[j], 2, 4, 0, 0);
        DrawWithOFUF(HistosEC[j], 0, 1);
        Format(HistosRE4[j], 2, 3, 0, 0);
        DrawWithOFUF(HistosRE4[j], 0, 2);
        // Legend
        Legend3(0.12, 0.70, 0.30, 0.88,
                   TString(TitlesHistos[j]),
                   HistosB[j], "Barrel",
                   HistosEC[j], "Endcap",
                   HistosRE4[j], "RE4",
                   0.028);
        std::string pngname= TitlesPNG[j];
        pad->Modified();
        MyCanvas[j]->Update();
        MyCanvas[j]->SaveAs("../summary/"+TString(pngname)+".png");
        
       }
}



void Title(TH1F* Hist,const char* Title,const char* XTitle,const char* YTitle, float TitleOffset)
{
	Hist->SetTitle(Title);
	Hist->GetXaxis()->SetTitle(XTitle);
	Hist->GetXaxis()->CenterTitle();	
	// Hist->GetXaxis()->SetTitleOffset(TitleOffset);
	Hist->GetYaxis()->SetTitle(YTitle);
	Hist->GetYaxis()->CenterTitle();
	Hist->GetYaxis()->SetTitleOffset(TitleOffset);
}

void Format(TH1F* Hist, int LineWidth, int LineColor, int FillStyle, int FillColor)
{
	Hist->SetLineWidth(LineWidth);
	Hist->SetLineColor(LineColor);
	Hist->SetFillStyle(FillStyle);
	Hist->SetFillColor(FillColor);
}

void Legend2(float x1, float y1, float x2, float y2,const char* Header, TH1F* Entry1,const char* Desc1, TH1F* Entry2,const char* Desc2, float TextSize)
{
	TLegend *Leg = new TLegend(x1, y1, x2, y2);
	Leg->AddEntry(Entry1, Desc1, "L");
	Leg->AddEntry(Entry2, Desc2, "L");
	//Leg->SetTextAlign(22);
	Leg->SetTextSize(TextSize);
	Leg->Draw();
}

void Legend3(float x1, float y1, float x2, float y2, const char* Header, TH1F* Entry1, const char* Desc1, TH1F* Entry2, const char* Desc2, TH1F* Entry3,const  char* Desc3, float TextSize)
{
	TLegend *Leg = new TLegend(x1, y1, x2, y2);
	Leg->AddEntry(Entry1, Desc1, "L");
	Leg->AddEntry(Entry2, Desc2, "L");	
	Leg->AddEntry(Entry3, Desc3, "L");
	//Leg->SetTextAlign(22);
	Leg->SetTextSize(TextSize);
	Leg->Draw();
}

void DrawWithOFUF(TH1F* Hist, bool Norm, int Same)
{
	const char* Name = Hist->GetName();
	const char* Title = Hist->GetTitle();

	Int_t Nx    = Hist->GetNbinsX()+2;
	Double_t BW = Hist->GetBinWidth(0);
	Double_t x1 = Hist->GetBinLowEdge(1)-BW;
	BW = Hist->GetBinWidth(Nx-1);
	Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW+BW;
	// Book a temporary histogram having extra bins_overflows and underflows
	///TH1F *HTmp = new TH1F(Name, Title, Nx, x1, x2);
        TH1F *HTmp = (TH1F*)Hist->Clone("hTmp");
        HTmp->SetName(Name);
	HTmp->SetBins(Nx, x1, x2);
        // Fill the new hitogram including the extra bin_overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());

	// Make title and format same as original
	//Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
	Format(HTmp, HTmp->GetLineWidth(), HTmp->GetLineColor(), HTmp->GetFillStyle(), HTmp->GetFillColor());

	if(Norm)
		HTmp->Scale(1/HTmp->Integral());
	HTmp->GetYaxis()->SetRangeUser(HTmp->GetMinimum()*1.2, HTmp->GetMaximum()*1.2);

	// Draw the temporary histogram
	if(Same==2)
	{
		HTmp->Draw("sames");
		// StatBox
		gPad->Update();
		TPaveStats* St3 = (TPaveStats*) HTmp->FindObject("stats");
		St3->SetX1NDC(0.74);
		St3->SetX2NDC(0.88);
		St3->SetY1NDC(0.45);
		St3->SetY2NDC(0.58);
	}
	else if(Same==1)
	{
		HTmp->Draw("sames");
		// StatBox
		gPad->Update();
		TPaveStats* St2 = (TPaveStats*) HTmp->FindObject("stats");
		St2->SetX1NDC(0.74);
		St2->SetX2NDC(0.88);
		St2->SetY1NDC(0.60);
		St2->SetY2NDC(0.73);
	}
	else if (Same==0)
	{
		HTmp->Draw();
		// Overflow Text
		TText *Tex = new TText(x2-BW/2, Hist->GetBinContent(Nx)+20, "Overflow");
		Tex->SetTextAngle(90);
		Tex->SetTextAlign(12);
		Tex->SetTextSize(0.03);
		Tex->SetTextColor(1);
		Tex->Draw();
		// Underflow Text
		TText *Tex2 = new TText(x1+BW/2, Hist->GetBinContent(0)+20, "Underflow");
		Tex2->SetTextAngle(90);
		Tex2->SetTextAlign(12);
		Tex2->SetTextSize(0.03);
		Tex2->SetTextColor(1);
		Tex2->Draw();
		// StatBox
		gPad->Update();
		TPaveStats* St1 = (TPaveStats*) HTmp->FindObject("stats");
		St1->SetX1NDC(0.74);
		St1->SetX2NDC(0.88);
		St1->SetY1NDC(0.75);
		St1->SetY2NDC(0.88);
	}


}

void DrawDiffWithOFUF(TH1F* Hist, bool Norm, bool Same)
{
	const char* Name = Hist->GetName();
	const char* Title = Hist->GetTitle();

	Int_t Nx    = Hist->GetNbinsX()+2;
	Double_t BW = Hist->GetBinWidth(0);
	Double_t x1 = Hist->GetBinLowEdge(1)-BW;
	BW = Hist->GetBinWidth(Nx-1);
	Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW+BW;
        Double_t ymax = Hist->GetMaximum;
	// Book a temporary histogram having extra bins_overflows and underflows
	TH1F *HTmp = new TH1F(Name, Title, Nx, x1, x2);

	// Fill the new hitogram including the extra bin_overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());
        
	// Make title and format same as original
	//Title(HTmp, HTmp->GetTitle(), HTmp->GetXaxis()->GetTitle(), HTmp->GetYaxis()->GetTitle(), HTmp->GetYaxis()->GetTitleOffset());
	Format(HTmp, Hist->GetLineWidth(), Hist->GetLineColor(), Hist->GetFillStyle(), Hist->GetFillColor());

	if(Norm)
		HTmp->Scale(1/HTmp->Integral());
	HTmp->GetYaxis()->SetRangeUser(HTmp->GetMinimum()*1.2, HTmp->GetMaximum()*1.2);

	// Draw the temporary histogram
	if(Same)
		HTmp->Draw("same");
	else
		HTmp->Draw();
}

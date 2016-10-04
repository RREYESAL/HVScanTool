#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <math.h>
#include <sstream>

#include "TH1.h"
#include "TMath.h"
#include "TGraph.h"

using namespace std; 

Double_t EffCompute(double hv,double S,double emax,double hv50 ){
	return emax / (1.0 + TMath::Exp( S *( hv - hv50 ) ) );}

void EffAtKneeDist(bool plots)
{
   gROOT->Reset();
        
   const int NUMPOINTS = 40; 
   Double_t x[NUMPOINTS];
   Double_t xerr[NUMPOINTS]={0};
   
   Double_t yB[NUMPOINTS];
   Double_t yBerr[NUMPOINTS];
   Double_t yEC[NUMPOINTS];
   Double_t yECerr[NUMPOINTS];
   Double_t yECre4[NUMPOINTS];
   Double_t yECre4err[NUMPOINTS];

   TH1F* HistosB[NUMPOINTS];
   TH1F* HistosEC[NUMPOINTS];
   TH1F* HistosRE4[NUMPOINTS];

   Double_t kneeEC, kneeB;        
   Double_t effatkneeEC, effatkneeB;        
   Double_t WorkingPoint, emax, HVfifty, S, Knee; 
   Char_t RollName[40]; 
   vector <double> wp, effmax, hv50, slope; 
   vector <std::string> roll;  

   TChain *ch = new TChain("ch","");
   ch->Add("../summary/barrel_summary_2016.root/T");
   ch->Add("../summary/endcap_summary_2016.root/T");

   TTree *T = (TTree*)ch;

   cout<<" Entries : "<<T->GetEntries()<<endl;
   Long64_t nentries = T->GetEntries();
   T->SetBranchAddress("RollName", &RollName);
   T->SetBranchAddress("WorkingPoint", &WorkingPoint);
   T->SetBranchAddress("emax", &emax);
   T->SetBranchAddress("HVfifty", &HVfifty);
   T->SetBranchAddress("S", &S);
   
   for (Long64_t i=0;i<nentries;i++){
	   T->GetEntry(i);
	   roll.push_back(RollName); 
	   wp.push_back(WorkingPoint); 
	   effmax.push_back(emax); 
	   hv50.push_back(HVfifty);
	   slope.push_back(S);
   }
    
 
   for (int j=0; j < NUMPOINTS; j++){
	   std::stringstream ss;
	   ss << j; 	   
	   TString name1 = "hEC"+ss.str();
	   TString name2 = "hBa"+ss.str();
	   TString name3 = "hEC4"+ss.str();
	   TH1F *hEC = new TH1F(name1,"EndCap",60, 60.,120.);  
	   TH1F *hBa = new TH1F(name2,"Barrel",60, 60.,120.);  
	   TH1F *hEC4 = new TH1F(name3,"re4",60, 60.,120.);
	   for (int k=0;k<roll.size(); k++ ){

		   if((roll.at(k)).find("W")!= std::string::npos){
                         kneeB=wp.at(k)-0.100; 
			 effatkneeB=EffCompute(kneeB+((j)*(0.02)),slope.at(k),effmax.at(k),hv50.at(k));
                         hBa->Fill(effatkneeB);
		   }
		   else if ((roll.at(k)).find("RE+4")!= std::string::npos || (roll.at(k)).find("RE-4")!= std::string::npos ){
           		 kneeEC=wp.at(k)-0.120; 
			 effatkneeEC=EffCompute(kneeEC+((j)*(0.02)),slope.at(k),effmax.at(k),hv50.at(k));
			 hEC4->Fill(effatkneeEC);
		   }
		   else {
                      kneeEC=wp.at(k)-0.120;
                      effatkneeEC=EffCompute(kneeEC+((j)*(0.02)),slope.at(k),effmax.at(k),hv50.at(k));
                      hEC->Fill(effatkneeEC);
		 }
	   }
	   *(x+j)    = j*(20.);
	   *(xerr+j) = 0;
           *(yB+j)    = hBa->GetMean(); 
	   *(yBerr+j) = (hBa->GetStdDev())/sqrt(hBa->GetEntries());      	 
	   *(yEC+j)    = hEC->GetMean(); 
	   *(yECerr+j) = (hEC->GetStdDev())/sqrt(hEC->GetEntries());
	   *(yECre4+j)    = hEC4->GetMean(); 
	   *(yECre4err+j) = (hEC4->GetStdDev())/sqrt(hEC4->GetEntries());
           *(HistosB+j) = hBa;
           *(HistosEC+j) = hEC;
           *(HistosRE4+j) = hEC4;

  } 
 gStyle->SetOptStat(1111); 
  int max=0; 
   
  // Canvas

  TCanvas *c1 = new TCanvas("c1","gerrors2",200,10,700,500);
  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->SetTicks();
  pad->Draw();
  pad->cd();
                  
  TH1F *h = pad->DrawFrame(0.0,0,850,105);
  h->SetStats(0);
  h->GetXaxis()->SetTitleOffset(1.1);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitle("effwp [%]");
  h->GetXaxis()->SetTitle("Added Voltage to the knee [V]");
  h->SetMaximum(100);
  h->SetMinimum(92);
  //h->Draw("text");
 
 
  TGraphErrors *gr1 = new TGraphErrors(NUMPOINTS,x,yB,xerr,yBerr);
  gr1->SetMarkerSize(0.9);
  gr1->SetMarkerColor(4);
  gr1->SetMarkerStyle(20);
  gr1->SetName("gr1");
  //gr1->Draw();
  TGraphErrors *gr2 = new TGraphErrors(NUMPOINTS,x,yEC,xerr,yECerr);
  gr2->SetMarkerSize(0.9);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(23);
  gr2->SetName("gr2");
  
  TGraphErrors *gr3 = new TGraphErrors(NUMPOINTS,x,yECre4,xerr,yECre4err);
  gr3->SetMarkerSize(.9);
  gr3->SetMarkerColor(3);
  gr3->SetMarkerStyle(21);
  gr3->SetName("gr3");

  TMultiGraph * mg = new TMultiGraph("mg","mg");
  mg->Add( gr1 );
  mg->Add( gr2 );
  mg->Add( gr3 );
  mg->Draw("LP");
 
  
  Leg = new TLegend(0.2, 0.7, 0.3, 0.85);
  Leg->AddEntry(gr1,"Barrel", "lpl");
  Leg->AddEntry(gr2,"EndCap", "lpl");
  Leg->AddEntry(gr3,"RE4", "lpl");
  Leg->Draw();
  if (plots){
  drawHistos(HistosB, NUMPOINTS, "Barrel");
  drawHistos(HistosEC, NUMPOINTS, "EndCap");
  drawHistos(HistosRE4, NUMPOINTS, "RE4");}
  return;

}

void drawtext()
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;

   TH1F *h = (TH1F*)gPad->GetListOfPrimitives()->FindObject("hframe");
   TGraphErrors *g = (TGraphErrors*)gPad->GetListOfPrimitives()->FindObject("gr1");
   n = g->GetN();
   cout << n  << endl;
   for (i=1; i<n; i++) {
      g->GetPoint(i,x,y);
      l = new TLatex(x,y+.2,Form("%2.1f",y));
      l->SetTextSize(0.025);
      l->SetTextFont(42);
      l->SetTextAlign(21);
      l->Paint();
   }
}
void Title(TH1F* Hist, char* Title, char* XTitle, char* YTitle, float TitleOffset)
{       
        Hist->SetTitle(Title);
        Hist->GetXaxis()->SetTitle(XTitle);
        //Hist->GetXaxis()->CenterTitle();        
        Hist->GetYaxis()->SetTitle(YTitle);
        //Hist->GetYaxis()->CenterTitle();
        Hist->GetYaxis()->SetTitleOffset(TitleOffset);
}

void Format(TH1F* Hist, int LineWidth, int LineColor, int FillStyle, int FillColor)
{       
        Hist->SetLineWidth(LineWidth);
        Hist->SetLineColor(LineColor);
        Hist->SetFillStyle(FillStyle);
        Hist->SetFillColor(FillColor);
}

void Legend2(float x1, float y1, float x2, float y2, char* Header, TH1F* Entry1, char* Desc1, TH1F* Entry2, char* Desc2, float TextSize)
{       
        Leg = new TLegend(x1, y1, x2, y2);
        Leg->AddEntry(Entry1, Desc1, "L");
        Leg->AddEntry(Entry2, Desc2, "L");
        //Leg->SetTextAlign(22);
        Leg->SetTextSize(TextSize);
        Leg->Draw();
}

void Legend3(float x1, float y1, float x2, float y2, char* Header, TH1F* Entry1, char* Desc1, TH1F* Entry2, char* Desc2, TH1F* Entry3, char* Desc3, float TextSize)
{       
        Leg = new TLegend(x1, y1, x2, y2);
        Leg->AddEntry(Entry1, Desc1, "L");
        Leg->AddEntry(Entry2, Desc2, "L");
        Leg->AddEntry(Entry3, Desc3, "L");
        //Leg->SetTextAlign(22);
        Leg->SetTextSize(TextSize);
        Leg->Draw();
}

void Legend5(float x1, float y1, float x2, float y2, char* Header, TH1F* Entry1, char* Desc1, TH1F* Entry2, char* Desc2, TH1F* Entry3, char* Desc3, TH1F* Entry4, char* Desc4, TH1F* Entry5, char* Desc5, float TextSize)
{       
        Leg = new TLegend(x1, y1, x2, y2);
        Leg->AddEntry(Entry1, Desc1, "L");
        Leg->AddEntry(Entry2, Desc2, "L");
        Leg->AddEntry(Entry3, Desc3, "L");
        Leg->AddEntry(Entry4, Desc4, "L");
        Leg->AddEntry(Entry5, Desc5, "L");
        //Leg->SetTextAlign(22);
        Leg->SetTextSize(TextSize);
        Leg->Draw();
}

void DrawWithOFUF(TH1F* Hist, bool Norm, int Same)
{
        char* Name = Hist->GetName();
        char* Title = Hist->GetTitle();

        Int_t Nx    = Hist->GetNbinsX()+2;
        Double_t BW = Hist->GetBinWidth(0);
        Double_t x1 = Hist->GetBinLowEdge(1);
        BW = Hist->GetBinWidth(Nx-1);
        Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW;

        // Book a temporary histogram having extra bins for overflows and underflows
        TH1F *HTmp = (TH1F*)Hist->Clone("hTmp"); 
        HTmp->SetName(Name);

        // Fill the new hitogram including the extra bin for overflows
        for (Int_t i=1; i<=Nx; i++)
                HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

        // Restore the number of entries
        HTmp->SetEntries(Hist->GetEntries());

        // Make title and format same as original
        Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
        Format(HTmp, Hist->GetLineWidth(), Hist->GetLineColor(), Hist->GetFillStyle(), Hist->GetFillColor());

        //HTmpIntgrl = HTmp->Integral();
        if(Norm)
                HTmp->Scale(1/HTmp->Integral());
        HTmp->GetYaxis()->SetRangeUser(HTmp->GetMinimum()*1.2, HTmp->GetMaximum()*1.2);

        // Draw the temporary histogram
        if(Same==4)
        {
                HTmp->Draw("sames");
                // StatBox
                //gPad->SetLogy();
                gPad->Update();
                TPaveStats* St5 = (TPaveStats*) HTmp->FindObject("stats");
                St5->SetX1NDC(0.74);
                St5->SetX2NDC(0.88);
                St5->SetY1NDC(0.15);
                St5->SetY2NDC(0.28);
        }
        else if(Same==3)
        {
                HTmp->Draw("sames");
                // StatBox
                gPad->Update();
                TPaveStats* St4 = (TPaveStats*) HTmp->FindObject("stats");
                St4->SetX1NDC(0.74);
                St4->SetX2NDC(0.88);
                St4->SetY1NDC(0.30);
                St4->SetY2NDC(0.43);
        }
        else if(Same==2)
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
                TText *Tex2 = new TText(x1+BW/2, Hist->GetBinContent(1)+20, "Underflow");
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
void drawHistos(TH1F** histo, int num, string module){

  for (int j = 0; j < num ;j++){
     std::stringstream ss;			                  
     ss << j*(20);
     string namestr = module+"_knee+"+ss.str();
     const char* name = namestr.c_str();							                  
     histo[j]->SetName(name);
     //histo[j]->GetXaxis()->SetRangeUser(60.,120.);
     histo[j]->Print();
          
   } 
 
  for (int j = 0; j < num ;j+=5){
        std::stringstream ss, ss1, ss2,ss3, ss4, ss5;
        ss << j;
        ss1 << (j+0)*20;
        ss2 << (j+1)*20;
        ss3 << (j+2)*20;
        ss4 << (j+3)*20;
        ss5 << (j+4)*20;
        string namestr;
        namestr =module+"_knee"+ss.str();
        const char* name = namestr.c_str();
        TCanvas* c = new TCanvas(name," ", 1024, 768);
        //MyCanvas->SetGrid(1, 1);
        c->SetTicks(1, 1);
        TPad *Pad1 = new TPad("Pad1", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
        Pad1->Draw();
        Pad1->SetTicks(1, 1);


        Pad1->cd();
        gStyle->SetOptStat(1111);
        if(module == "Barrel")Title(histo[j], "Barrel Efficiency", "Efficiency[%]", "Number of Rolls", 1.2);
        else if(module == "EndCap")Title(histo[j], "EndCap Efficiency", "Efficiency[%]", "Number of Rolls", 1.2);
        else if(module == "RE4")Title(histo[j], "RE4 Efficiency", "Efficiency[%]", "Number of Rolls", 1.2);
        Format(histo[j], 2, 2, 0, 0);
        DrawWithOFUF(histo[j], 0, 0);
        Format(histo[j+1], 2, 3, 0, 0);
        DrawWithOFUF(histo[j+1], 0, 1);
        Format(histo[j+2], 2, 4, 0, 0);
        DrawWithOFUF(histo[j+2], 0, 2);
        Format(histo[j+3], 2, 5, 0, 0);
        DrawWithOFUF(histo[j+3], 0, 3);
        Format(histo[j+4], 2, 6, 0, 0);
        DrawWithOFUF(histo[j+4], 0, 4);
        // Legend
        Legend5(0.12, 0.70, 0.30, 0.88,
                "effwp",
                histo[j+0], ("knee+"+ss1.str()+" V").c_str(),
                histo[j+1], ("knee+"+ss2.str()+" V").c_str(),
                histo[j+2], ("knee+"+ss3.str()+" V").c_str(),
                histo[j+3], ("knee+"+ss4.str()+" V").c_str(),
                histo[j+4], ("knee+"+ss5.str()+" V").c_str(),
                0.030);
       
       if(module == "Barrel") c->SaveAs("Barrel_EffKnee"+TString(ss.str())+".png");       
       else if(module == "EndCap") c->SaveAs("EndCap_EffKnee"+TString(ss.str())+".png");       
       else if(module == "RE4") c->SaveAs("RE4_EffKnee"+TString(ss.str())+".png");       
       c->Modified();
       c->Update();

 } 


}





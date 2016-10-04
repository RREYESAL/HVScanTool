#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include "TCanvas.h"
#include "TAxis.h"
#include <TCanvas.h>
#include <TAxis.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <TChain.h>

#include "Data.C"

using namespace TMath;
using namespace std;

//int cpu = 2;

void WPChannel()
{
 
    //---Define the Roll Names and variables      
    vector<string> rollBAR, rollEND; 
    vector<Float_t> wpBAR, wpEND, wp_chBAR, wp_chEND;
    Float_t wp_chEC, wp_chB;
    Float_t MAX=0.0; 
    Float_t MIN=1000000.0;
    Float_t MEAN,sum;
    int nrolls;   
 
   TH1F* hmaxEC = new TH1F ("WPmaxEC","#Delta_{wp}(MAX_{ROLL} - WP_{CH})",100,0.,600.);
   TH1F* hmaxB = new TH1F ("WPmaxB","#Delta_{wp}(MAX_{ROLL} - WP_{CH})",100,0.,600.);
   TH1F* hminEC = new TH1F ("WPminEC","#Delta_{wp}(WP_{CH} - MIN_{ROLL})",100,0,500.);  
   TH1F* hminB = new TH1F ("WPminB","#Delta_{wp}(WP_{CH} - MIN_{ROLL})",100,0,600.);
   TH1F* hwpEC = new TH1F ("WPEC","(WP(CH)",100,8800,10200.);  
   TH1F* hwpB = new TH1F ("WPB","WP(CH)",50,8800,10200.);
   TH1F* hwpRE4 = new TH1F ("WPRE4","WP(CH)",50,8800,10200.);
   
   ofstream wp_ch_EC, wp_ch_B; 
   
   //--- Get the Name of the Chambers
    ifstream map_channel; 
    //map_channel.open("ChambersName_EndCap_noRE4.txt");  
    map_channel.open("ChambersName_EndCap.txt");  
    string CH0, CH1;
    vector<string> ch_roll0,ch_roll1;      
    while (map_channel.good()) {
        map_channel >> CH0
                    >> CH1;
        if(map_channel.eof()) break; 
        ch_roll0.push_back(CH0+"_"); ch_roll1.push_back(CH1+"_"); 
    }      
    map_channel.close();
     
   
    map_channel.open("ChambersName_Barrel.txt");
    string CHB;
    vector<string> chB_roll;
    while (map_channel.good()) {
        map_channel >> CHB; 
  	if (map_channel.eof())break;        
        chB_roll.push_back(CHB+"_"); 
    }
    map_channel.close();    
      

    //--- Get the data summary   
    std::vector<int> samples;
    samples.push_back(1);
    samples.push_back(-1);
    gROOT->ProcessLine(".L Data.C+");

    for (std::vector<int>::const_iterator Sample = samples.begin(); Sample != samples.end(); Sample++) {
        TChain *ch = new TChain("T");
        
 	if ( (*Sample)==-1 ) { ch->Add("../BlackList/endcap_summary_2016.root/T"); }
        if ( (*Sample)==1 ) { ch->Add("../BlackList/barrel_summary_2016.root/T"); }


        TTree *tree = (TTree*) ch;
        Data t(tree);
        Long64_t nentries = t.fChain->GetEntries();
        cout <<" Entries : " << nentries <<endl;
        Int_t nbytes = 0;
        Int_t nb = 0;
        for(Long64_t jentry=0; jentry<nentries;jentry++)
        {
            Long64_t ientry = t.LoadTree(jentry);
            if (ientry < 0)break;
            nb = t.fChain->GetEntry(jentry); nbytes += nb;
            if(jentry==nentries-1) cout<<endl;

            if ( (*Sample)==-1 ) {
                                   rollEND.push_back(t.RollName);
                                   wpEND.push_back(1000*t.WorkingPoint);
                                 }

            if ( (*Sample)==1 ) {
                                  rollBAR.push_back(t.RollName);
                                  wpBAR.push_back(1000*t.WorkingPoint);
				}
       
        }



   }
 
   //--- Define the HV applied to each chamber 
 
   string stroll1, stroll2;   
   //+++ For the End Cap 

     std::vector< std::pair<std::string,std::string >> MapWPCHEC;
     std::vector< std::pair<std::string,std::string >> MapRollCHEC;

   string::size_type match1, matchmaxEC;
   string::size_type matchREp4, matchREm4;
   Int_t MAXLESS_100, MAXGREATER_100; 
   Int_t MINLESS_100, MINGREATER_100; 
   for (int i=0; i<ch_roll0.size(); i++){
     	stroll1="c_";
     	stroll2="c_";
        MIN=100000.0; 
     	MAX=0.0;
     	sum = 0.0;
     	nrolls = 0;
        std::pair<std::string,std::string > pair_mapWP;
        std::pair<std::string,std::string > pair_mapROLL;
 	for (int j=0; j<rollEND.size(); j++){ 
       	   match1= (rollEND.at(j)).find(ch_roll0.at(i)); 
           matchmaxEC= (rollEND.at(j)).find(ch_roll1.at(i));
          if(match1!=std::string::npos || matchmaxEC!=std::string::npos) {
		   std::stringstream ss2;
           	   ss2<<wpEND.at(j);  
		  
		  if (match1!=std::string::npos){
			  if ((rollEND.at(j)).find("_A")!=std::string::npos)stroll2 +="_wpCH1A_"+ss2.str();
		          else if ((rollEND.at(j)).find("_B")!=std::string::npos)stroll2 +="_wpCH1B_"+ss2.str();
		          else stroll2 +="_wpCH1C_"+ss2.str();
		  }
		  else{

			  if ((rollEND.at(j)).find("_A")!=std::string::npos)stroll2 +="_wpCH2A_"+ss2.str();
		          else if ((rollEND.at(j)).find("_B")!=std::string::npos)stroll2 +="_wpCH2B_"+ss2.str();
		          else stroll2 +="_wpCH2C_"+ss2.str();

		  }
	   	   stroll1 +=rollEND.at(j)+"_";
	   	   if(wpEND.at(j)<MIN)MIN=wpEND.at(j);
           	   if(wpEND.at(j)>MAX)MAX=wpEND.at(j);
	   	   sum +=wpEND.at(j);
           	   ++nrolls;      
	   }
      
        }

      	if (nrolls == 0 )MEAN = 0;
      	else MEAN = sum/nrolls;
        std::stringstream ss0, ss1;
        ss0 << nrolls; 
          
    	
        wp_chEC = (MAX-MIN >= 100.) ? MIN+100.0: MEAN ;  
        //++++++++++++ 
        pair_mapROLL.first =ch_roll0.at(i)+ch_roll1.at(i);
        pair_mapROLL.second =stroll2;
	MapRollCHEC.push_back(pair_mapROLL);
	
	pair_mapWP.first =ch_roll0.at(i)+ch_roll1.at(i);
	ss1<<wp_chEC;  
        pair_mapWP.second =ss0.str()+"_wp"+ss1.str();  
	MapWPCHEC.push_back(pair_mapWP);

	matchREp4= (ch_roll0.at(i)).find("RE+4");
        matchREm4= (ch_roll0.at(i)).find("RE-4");   	
        if(nrolls !=0){ 
              if (matchREm4!=std::string::npos || matchREp4!=std::string::npos  )hwpRE4->Fill(wp_chEC);
  	      else hwpEC->Fill(wp_chEC); 
              
              hmaxEC->Fill(MAX-wp_chEC);
     	      hmaxEC->Fill(wp_chEC-MIN);
      	}
        
        if(MAX-wp_chEC<100){
        	++MAXLESS_100;}
    	else ++MAXGREATER_100;
        if(wp_chEC-MIN<100)++MINLESS_100;
        else ++MINGREATER_100;
          
   
   }
   //+++ For the Barrel  
   std::vector< std::pair<std::string,std::string >> MapWPCHB;
   std::vector< std::pair<std::string,std::string >> MapRollCHB;
   string::size_type match;
   Int_t MAXLESSB_100, MAXGREATERB_100;
   Int_t MINLESSB_100, MINGREATERB_100;
   for (int i=0; i<chB_roll.size(); i++){
        stroll1="c_";
	stroll2="c_";
	MIN=100000.0;
        MAX=0.0;
        sum = 0.0;
        nrolls = 0;
        std::pair<std::string,std::string > pair_mapWPB;
        std::pair<std::string,std::Int_t > pair_mapNumB;
        std::pair<std::string,std::string > pair_mapROLLB;
        for (int j=0; j<rollBAR.size(); j++){
           std::stringstream ss2;
           ss2<<wpBAR.at(j);  
           match= (rollBAR.at(j)).find(chB_roll.at(i));
           if(match!=std::string::npos){
            if ((rollBAR.at(j)).find("_Forward")!=std::string::npos)stroll2 +="_wpFor_"+ss2.str();
            else if ((rollBAR.at(j)).find("_Backward")!=std::string::npos)stroll2 +="_wpBac_"+ss2.str();
            else stroll2 +="_wpMid_"+ss2.str();


           if(wpBAR.at(j)<MIN)MIN=wpBAR.at(j);
           if(wpBAR.at(j)>MAX)MAX=wpBAR.at(j);
           sum +=wpBAR.at(j);
           ++nrolls;
           }

        }

        if (nrolls == 0 )MEAN = 0.;
        else MEAN = sum/nrolls;
        
       wp_chB = (MAX-MIN >= 100.) ? MIN+100.0: MEAN !=0 ? MEAN : 0 ;
       
       if(nrolls!=0){ 
       hwpB->Fill(wp_chB);   
       hmaxB->Fill(MAX-wp_chB);
       hmaxB->Fill(wp_chB-MIN);
       }
 
        pair_mapROLLB.first =chB_roll.at(i);
        pair_mapROLLB.second =stroll2;
	MapRollCHB.push_back(pair_mapROLLB);
	
        std::stringstream ss0, ss1;
        ss0 << nrolls; 
	pair_mapWPB.first =chB_roll.at(i);
	ss1<<wp_chB;  
        pair_mapWPB.second =ss0.str()+"_wp"+ss1.str();  
	MapWPCHB.push_back(pair_mapWPB);
       
       
       if(MAX-wp_chB<=100)++MAXLESSB_100;
       else  ++MAXGREATERB_100;
       hminB->Fill(wp_chB-MIN);
       if(wp_chB-MIN<100)++MINLESSB_100;
       else ++MINGREATERB_100;
       //wp_ch_B << chB_roll.at(i) << "   wpCH=" << wp_chB<< endl;        

  }
  // wp_ch_B.close();
  //ofstream WPCHout;
  //WPCHout.open("WPpCH_out.txt");
  for (vector <std::pair<std::string,std::string>>::const_iterator it =  MapRollCHEC.begin() ;it != MapRollCHEC.end(); it++ ){
       //cout <<it->first<<" "<<it->second<<"  "<< endl;
      for (vector <std::pair<std::string,std::string>>::const_iterator itr =  MapWPCHEC.begin() ;itr != MapWPCHEC.end(); itr++  ){
     		if(it->first != itr->first)continue;
	      	   //WPCHout <<it->first<<"  "<<itr->second <<" "<<it->second<<endl;  
	      	   cout  <<it->first<<"  "<<itr->second <<" "<<it->second<<endl;  
      }
  }
 /* 
  //WPCHout.seekp(0L,ios::end); 
  for (vector <std::pair<std::string,std::string>>::const_iterator it =  MapRollCHB.begin() ;it != MapRollCHB.end(); it++ ){
       //cout <<it->first<<" "<<it->second<<"  "<< endl;
      for (vector <std::pair<std::string,std::string>>::const_iterator itr =  MapWPCHB.begin() ;itr != MapWPCHB.end(); itr++  ){
     		if(it->first != itr->first)continue;
	      	   cout <<it->first<<"  "<<itr->second <<" "<<it->second<<endl;  
      }
  }*/
  /*
  WPCHout.close();
  system("sed -e 's:_wp: :g' -e 's:_A_:_A :g' -e 's:_B_:_B :g' -e 's:_C_:_C :g' WPpCH_out.txt > WPpCH_out_.txt");
  system("mv WPpCH_out_.txt WPpCH_out.txt");
 */
  
  
  
  // ---- Plot your results 
  
  gStyle->SetOptStat(111111);
  TCanvas* c = new TCanvas("c","c", 1024, 768);
  c->SetTicks(1, 1);
  TPad *Pad1 = new TPad("Pad1", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
  Pad1->Draw();
  Pad1->SetTicks(1, 1);

  Pad1->cd();

  Title(hwpB,"WP_{CH} Black-list "," WP_{ch}[V] ", "Events", 1.2);
  Format(hwpB, 2, 2, 3001, 2);
  DrawWithOFUF(hwpB, 0, 0); 
  Format(hwpEC, 2, 4, 3001, 4);
  DrawWithOFUF(hwpEC, 0, 1);
  Format(hwpRE4, 2, 3, 3001, 3);
  DrawWithOFUF(hwpRE4, 0, 2);

  Legend3(0.12, 0.75, 0.28, 0.88,
       "hwp",
        hwpB, "WP_{CH} Barrel",
        hwpEC,"WP_{CH} EnCap",
        hwpRE4,"WP_{CH} RE4",
        0.028);  

 c->Update();
 Pad1->cd();
 Title(hmaxEC,"Diff MAX MIN wrt WP_{CH}"," WP_{ch}[V] ", "Events", 1.2);
 Format(hmaxEC, 2, 4, 3001, 4);
 DrawWithOFUF(hmaxEC, 0, 0);
 Format(hmaxB, 2, 2, 3001, 4);
 DrawWithOFUF(hmaxB, 0, 1);

 Legend2(0.12, 0.75, 0.28, 0.88,
       "Diff MAX MIN wrt WP_{CH}",
        hmaxB, "Barrel",
        hmaxEC,"EnCap",
        0.028);

			 
 return;

}


void Title(TH1F* Hist, TString &Title, char* XTitle, char* YTitle, float TitleOffset)
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
        Leg->SetTextSize(TextSize);
        Leg->Draw();
}

void DrawWithOFUF(TH1F* Hist, bool Norm, int Same)
{
	char* Name = Hist->GetName();
	char* Title = Hist->GetTitle();

	Int_t Nx    = Hist->GetNbinsX()+2;
	Double_t BW = Hist->GetBinWidth(0);
	Double_t x1 = Hist->GetBinLowEdge(1)-BW;
	BW = Hist->GetBinWidth(Nx-1);
	Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW+BW;

	// Book a temporary histogram having extra bins for overflows and underflows
        TH1F *HTmp = (TH1F*)Hist->Clone("hTmp");
        HTmp->SetName(Name);
        HTmp->SetBins(Nx, x1, x2);
	// Fill the new hitogram including the extra bin for overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());

	// Make title and format same as original
	Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
	Format(HTmp, Hist->GetLineWidth(), Hist->GetLineColor(), Hist->GetFillStyle(), Hist->GetFillColor());

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
	char* Name = Hist->GetName();
	char* Title = Hist->GetTitle();

	Int_t Nx    = Hist->GetNbinsX()+2;
	Double_t BW = Hist->GetBinWidth(0);
	Double_t x1 = Hist->GetBinLowEdge(1)-BW;
	BW = Hist->GetBinWidth(Nx-1);
	Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW+BW;

	// Book a temporary histogram having extra bins for overflows and underflows
        TH1F *HTmp = (TH1F*)Hist->Clone("hTmp");
        HTmp->SetName(Name);
        HTmp->SetBins(Nx, x1, x2);
	// Fill the new hitogram including the extra bin for overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());

	// Make title and format same as original
	Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
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


#include <iostream>
#include "TPad.h"
#include "TCanvas.h"
#include "TH1F.h"
//#include "gROOT.h"

using namespace std; 

void overflow()
{

//gROOT->LoadMacro("tdrstyle.C");
//	setTDRStyle();

	ifstream barrel, endcap, re4;
        
	barrel.open("results_summary2/barrel_summary.txt");
  	endcap.open("results_summary2/endcap_summary.txt");
  	re4.open("results_summary2/re4_summary.txt");

          cout<<"Open succesful"<<endl;

	// The Variables
	char RollName[25];
	string clswp, clsWpDef;
	double wp, slope50, emax, hv50, chi2, effwp, wpDef, effWpDef, effmaxerror, Serror, hv50error, S, parameter1;

	// Histos
	TH1F* wpBAR = new TH1F("wp for Barrel", "Working Point", 70, 8.6, 10);
	TH1F* wpEND = new TH1F("wp for Endcap", "Working Point", 70, 8.6, 10);
	TH1F* wpRE4 = new TH1F("wp for RE4", "Working Point ", 70, 8.6, 10);

	TH1F* hv50BAR = new TH1F("hv50 for Barrel", "HV50", 70, 8.6, 10);
	TH1F* hv50END = new TH1F("hv50 for Endcap", "HV50", 70, 8.6, 10);
	TH1F* hv50RE4 = new TH1F("hv50 for RE4", "HV50 ", 70, 8.6, 10);

	TH1F* kneeBAR = new TH1F("knee for Barrel", "Knee", 70, 8.6, 10);
	TH1F* kneeEND = new TH1F("knee for Endcap", "Knee", 70, 8.6, 10);
	TH1F* kneeRE4 = new TH1F("knee for RE4", "Knee ", 70, 8.6, 10);

	TH1F* emaxBAR = new TH1F("emax for Barrel", "emax", 25, 70, 120);
	TH1F* emaxEND = new TH1F("emax for Endcap", "emax", 25, 70, 120);
	TH1F* emaxRE4 = new TH1F("emax for RE4", "emax ", 25, 70, 120);

	TH1F* slope50BAR = new TH1F("slope50 for Barrel", "slope50", 50, 0, 400);
	TH1F* slope50END = new TH1F("slope50 for Endcap", "slope50", 50, 0, 400);
	TH1F* slope50RE4 = new TH1F("slope50 for RE4", "slope50 ", 50, 0, 400);

	TH1F* chi2BAR = new TH1F("chi2 for Barrel", "chi2", 44, -2, 20);
	TH1F* chi2END = new TH1F("chi2 for Endcap", "chi2", 44, -2, 20);
	TH1F* chi2RE4 = new TH1F("chi2 for RE4", "chi2 ", 44, -2, 20);
       
	TH1F* effwpBAR = new TH1F("effwp for Barrel", "effwp", 30, 55, 115);
	TH1F* effwpEND = new TH1F("effwp for Endcap", "effwp", 30, 55, 115);
	TH1F* effwpRE4 = new TH1F("effwp for RE4", "effwp ", 30, 55, 115);


	// Read the data
	while(1)
	{
		barrel>>RollName>>wp>>slope50>>emax>>hv50>>chi2>>clswp>>effwp>>wpDef>>effWpDef>>clsWpDef>>effmaxerror>>Serror>>hv50error>>S>>parameter1;
		if (!barrel.good()) break;
	
			wpBAR->Fill(wp);
			hv50BAR->Fill(hv50);
			kneeBAR->Fill(wp-0.1);
			emaxBAR->Fill(emax);
			slope50BAR->Fill(slope50);
			chi2BAR->Fill(chi2);
			effwpBAR->Fill(effwp);

	}

	while(1)
	{
		endcap>>RollName>>wp>>slope50>>emax>>hv50>>chi2>>clswp>>effwp>>wpDef>>effWpDef>>clsWpDef>>effmaxerror>>Serror>>hv50error>>S>>parameter1;
		if (!endcap.good()) break;
	
			wpEND->Fill(wp);
			hv50END->Fill(hv50);
			kneeEND->Fill(wp-0.12);
			emaxEND->Fill(emax);
			slope50END->Fill(slope50);
			chi2END->Fill(chi2);
			effwpEND->Fill(effwp);

	}

	while(1)
	{
		re4>>RollName>>wp>>slope50>>emax>>hv50>>chi2>>clswp>>effwp>>wpDef>>effWpDef>>clsWpDef>>effmaxerror>>Serror>>hv50error>>S>>parameter1;
		if (!re4.good()) break;
	
			wpRE4->Fill(wp);
			hv50RE4->Fill(hv50);
			kneeRE4->Fill(wp-0.12);
			emaxRE4->Fill(emax);
			slope50RE4->Fill(slope50);
			chi2RE4->Fill(chi2);
			effwpRE4->Fill(effwp);

	}
	

	// Canvas
	TCanvas* c1 = new TCanvas("c1","c", 1024, 768);
	c1->SetGrid(1, 1);
	c1->SetTicks(1, 1);
	TPad *Pad1 = new TPad("Pad1", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad1->Draw();
	Pad1->SetGrid(1, 1);
	Pad1->SetTicks(1, 1);
	TCanvas* c2 = new TCanvas("c2","c", 1024, 768);
	c2->SetGrid(1, 1);
	c2->SetTicks(1, 1);
	TPad *Pad2 = new TPad("Pad2", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad2->Draw();
	Pad2->SetGrid(1, 1);
	Pad2->SetTicks(1, 1);
	TCanvas* c3 = new TCanvas("c3","c", 1024, 768);
	c3->SetGrid(1, 1);
	c3->SetTicks(1, 1);
	TPad *Pad3= new TPad("Pad3", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad3->Draw();
	Pad3->SetGrid(1, 1);
	Pad3->SetTicks(1, 1);
	TCanvas* c4 = new TCanvas("c4","c", 1024, 768);
	c4->SetGrid(1, 1);
	c4->SetTicks(1, 1);
	TPad *Pad4 = new TPad("Pad4", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad4->Draw();
	Pad4->SetGrid(1, 1);
	Pad4->SetTicks(1, 1);
	TCanvas* c5 = new TCanvas("c5","c", 1024, 768);
	c5->SetGrid(1, 1);
	c5->SetTicks(1, 1);
	TPad *Pad5 = new TPad("Pad5", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad5->Draw();
	Pad5->SetGrid(1, 1);
	Pad5->SetTicks(1, 1);
	TCanvas* c6 = new TCanvas("c6","c", 1024, 768);
	c6->SetGrid(1, 1);
	c6->SetTicks(1, 1);
	TPad *Pad6 = new TPad("Pad6", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad6->Draw();
	Pad6->SetGrid(1, 1);
	Pad6->SetTicks(1, 1);
	TCanvas* c7 = new TCanvas("c7","c", 1024, 768);
	c7->SetGrid(1, 1);
	c7->SetTicks(1, 1);
	TPad *Pad7 = new TPad("Pad7", "The pad 100% of the height", 0.0, 0.0, 1.0, 1.0, 0);
	Pad7->Draw();
	Pad7->SetGrid(1, 1);
	Pad7->SetTicks(1, 1);


	// Draw the histos
	// wp
	c1->cd();
	Pad1->cd();
	gStyle->SetOptStat(1111);
	Title(wpEND, "Working Point", "WP", "Events", 1.2);
	Format(wpEND, 2, 4, 0, 0);
	DrawWithOFUF(wpEND, 0, 0);
	Format(wpBAR, 2, 2, 0, 0);
	DrawWithOFUF(wpBAR, 0, 1);
	Format(wpRE4, 2, 3, 0, 0);
	DrawWithOFUF(wpRE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.28, 0.88,
		"wp",
		wpBAR, "wp_Barrel",
		wpEND, "wp_Endcap",
		wpRE4, "wp_RE4",
		0.030);

	c1->Update();
	c1->SaveAs("histograms/wp.png");

	// hv50
	c2->cd();
	Pad2->cd();
	gStyle->SetOptStat(1111);
	Title(hv50END, "HV50", "HV50", "Events", 1.2);
	Format(hv50END, 2, 4, 0, 0);
	DrawWithOFUF(hv50END, 0, 0);
	Format(hv50BAR, 2, 2, 0, 0);
	DrawWithOFUF(hv50BAR, 0, 1);
	Format(hv50RE4, 2, 3, 0, 0);
	DrawWithOFUF(hv50RE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.28, 0.88,
		"hv50",
		hv50BAR, "hv50_Barrel",
		hv50END, "hv50_Endcap",
		hv50RE4, "hv50_RE4",
		0.028);

	c2->SaveAs("histograms/hv50.png");

	// knee
	c3->cd();
	Pad3->cd();
	gStyle->SetOptStat(1111);
	Title(kneeEND, "Knee", "Knee", "Events", 1.2);
	Format(kneeEND, 2, 4, 0, 0);
	DrawWithOFUF(kneeEND, 0, 0);
	Format(kneeBAR, 2, 2, 0, 0);
	DrawWithOFUF(kneeBAR, 0, 1);
	Format(kneeRE4, 2, 3, 0, 0);
	DrawWithOFUF(kneeRE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.28, 0.88,
		"knee",
		kneeBAR, "knee_Barrel",
		kneeEND, "knee_Endcap",
		kneeRE4, "knee_RE4",
		0.028);

	c3->SaveAs("histograms/knee.png");


	// emax
	c4->cd();
	Pad4->cd();
	gStyle->SetOptStat(1111);
	Title(emaxBAR, "emax", "emax", "Events", 1.2);
	Format(emaxBAR, 2, 2, 0, 0);
	DrawWithOFUF(emaxBAR, 0, 0);
	Format(emaxEND, 2, 4, 0, 0);
	DrawWithOFUF(emaxEND, 0, 1);
	Format(emaxRE4, 2, 3, 0, 0);
	DrawWithOFUF(emaxRE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.30, 0.88,
		"emax",
		emaxBAR, "emax_Barrel",
		emaxEND, "emax_Endcap",
		emaxRE4, "emax_RE4",
		0.030);

	c4->SaveAs("histograms/emax.png");

	// slope50
	c5->cd();
	Pad5->cd();
	gStyle->SetOptStat(1111);
	Title(slope50BAR, "Slope50", "Slope50", "Events", 1.2);
	Format(slope50BAR, 2, 2, 0, 0);
	DrawWithOFUF(slope50BAR, 0, 0);
	Format(slope50END, 2, 4, 0, 0);
	DrawWithOFUF(slope50END, 0, 1);
	Format(slope50RE4, 2, 3, 0, 0);
	DrawWithOFUF(slope50RE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.30, 0.88,
		"slope50",
		slope50BAR, "slope50_Barrel",
		slope50END, "slope50_Endcap",
		slope50RE4, "slope50_RE4",
		0.028);

	c5->SaveAs("histograms/slope50.png");

	// chi2
	c6->cd();
	Pad6->cd();
	gStyle->SetOptStat(1111);
	Title(chi2BAR, "Chi2", "Chi2", "Events", 1.2);
	Format(chi2END, 2, 4, 0, 0);
	DrawWithOFUF(chi2END, 0, 0);
    Format(chi2BAR, 2, 2, 0, 0);
    DrawWithOFUF(chi2BAR, 0, 1);
	Format(chi2RE4, 2, 3, 0, 0);
	DrawWithOFUF(chi2RE4, 0, 2);
	// Legend
	Legend3(0.32, 0.75, 0.48, 0.88,
		"chi2",
		chi2BAR, "chi2_Barrel",
		chi2END, "chi2_Endcap",
		chi2RE4, "chi2_RE4",
		0.030);

	c6->SaveAs("histograms/chi2.png");

	// effwp
	c7->cd();
	Pad7->cd();
	gStyle->SetOptStat(1111);
	Title(effwpBAR, "Eff at WP", "Eff at WP", "Events", 1.2);
	Format(effwpEND, 2, 4, 0, 0);
	DrawWithOFUF(effwpEND, 0, 0);
    Format(effwpBAR, 2, 2, 0, 0);
    DrawWithOFUF(effwpBAR, 0, 1);
	Format(effwpRE4, 2, 3, 0, 0);
	DrawWithOFUF(effwpRE4, 0, 2);
	// Legend
	Legend3(0.12, 0.75, 0.30, 0.88,
		"effwp",
		effwpBAR, "effwp_Barrel",
		effwpEND, "effwp_Endcap",
		effwpRE4, "effwp_RE4",
		0.030);

	c7->SaveAs("histograms/effwp.png");


	endcap.close();
	barrel.close();
	re4.close();


}

void Title(TH1F* Hist, char* Title, char* XTitle, char* YTitle, float TitleOffset)
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
	Double_t x1 = Hist->GetBinLowEdge(1)-BW;
	BW = Hist->GetBinWidth(Nx-1);
	Double_t x2 = Hist->GetBinLowEdge(Nx-2)+BW+BW;

	// Book a temporary histogram having extra bins for overflows and underflows
	TH1F *HTmp = new TH1F(Name, Title, Nx, x1, x2);

	// Fill the new hitogram including the extra bin for overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());

	// Make title and format same as original
	Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
	Format(HTmp, Hist->GetLineWidth(), Hist->GetLineColor(), Hist->GetFillStyle(), Hist->GetFillColor());

	HTmpIntgrl = HTmp->Integral();
	if(Norm)
		HTmp->Scale(1/HTmpIntgrl);
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
	TH1F *HTmp = new TH1F(Name, Title, Nx, x1, x2);

	// Fill the new hitogram including the extra bin for overflows
	for (Int_t i=1; i<=Nx; i++)
		HTmp->Fill(HTmp->GetBinCenter(i), Hist->GetBinContent(i-1));

	// Restore the number of entries
	HTmp->SetEntries(Hist->GetEntries());

	// Make title and format same as original
	Title(HTmp, Hist->GetTitle(), Hist->GetXaxis()->GetTitle(), Hist->GetYaxis()->GetTitle(), Hist->GetYaxis()->GetTitleOffset());
	Format(HTmp, Hist->GetLineWidth(), Hist->GetLineColor(), Hist->GetFillStyle(), Hist->GetFillColor());

	HTmpIntgrl = HTmp->Integral();
	if(Norm)
		HTmp->Scale(1/HTmpIntgrl);
	HTmp->GetYaxis()->SetRangeUser(HTmp->GetMinimum()*1.2, HTmp->GetMaximum()*1.2);

	// Draw the temporary histogram
	if(Same)
		HTmp->Draw("same");
	else
		HTmp->Draw();
}

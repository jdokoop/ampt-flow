#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TText.h"

using namespace std;

//------------------------------------------------
// Global Variables
//------------------------------------------------

const int NPT = 10;
const float dPhi = 0.0785398;
const float dEta = 0.05;

string collSyst = "dAu";

//The pT bins
float ptl[NPT] = {0.20, 0.34, 0.48, 0.62, 0.76, 0.90, 1.04, 1.18, 1.32, 1.46};
float pth[NPT] = {0.34, 0.48, 0.62, 0.76, 0.90, 1.04, 1.18, 1.32, 1.46, 1.60};

TFile *inputFile;
TFile *inputFile_pp;

//Histograms from input file
TH3F *hdEta_dPhi_pT;
TH3F *hdEta_dPhi_pT_pp;

TH2F *hEta_pT;
TH2F *hEta_pT_pp;

//Two-particle correlations
TH2D *hdEta_dPhi[NPT];
TH2D *hdEta_dPhi_pp[NPT];

//Correlations as a function of pT
TH1F *hCorrelation_He3Au[NPT];
TH1F *hCorrelation_pp[NPT];

TH1F *hCorrelation_NoZYAM_He3Au[NPT];
TH1F *hCorrelation_NoZYAM_pp[NPT];

TH1F *hCorrelation_NoZYAM_NoFrag_He3Au[NPT];

//Fit functions for ZYAM
float bZYAM[NPT];
float bZYAM_pp[NPT];
TF1 *ZYAMFit[NPT];
TF1 *ZYAMFit_pp[NPT];

//Run configuration
float sigma = 1.5;
int etaLim = 2;
float deltaEta_lo = 2;
float deltaEta_hi = 3;

//Fourier coefficients
double an_coeffs[NPT][5];
double vn_coeffs[NPT][5];
double vn_err_coeffs[NPT][5];

TGraphErrors *gv2;
TGraphErrors *gv3;
TGraphErrors *gv4;

//Fourier reconstructions of correlation functions
TF1 *reconstruction[NPT];

void readHistogramsFromFile()
{
	//Raw data file names
	string filename = "/Users/jdok/Documents/CUBoulder/JNLab/Spring15/He3_AMPT/RawData/ampt_"+collSyst+"_200GeV.root";
	inputFile = new TFile(filename.c_str());
	inputFile_pp = new TFile("/Users/jdok/Documents/CUBoulder/JNLab/Spring15/He3_AMPT/RawData/ampt_pp_200GeV.root");

	//Get histograms from file
	hdEta_dPhi_pT = (TH3F*) inputFile->Get("hdEta_dPhi_pT");
	hdEta_dPhi_pT->SetName("hdEta_dPhi_pT");
	hdEta_dPhi_pT_pp = (TH3F*) inputFile_pp->Get("hdEta_dPhi_pT");
	hdEta_dPhi_pT_pp->SetName("hdEta_dPhi_pT_pp");

	hEta_pT = (TH2F*) inputFile->Get("hEta_pT");
	hEta_pT->SetName("hEta_pT");
	hEta_pT_pp = (TH2F*) inputFile_pp->Get("hEta_pT");
	hEta_pT_pp->SetName("hEta_pT_pp");
}

void prepareInputHistograms()
{
	int binLow = 0;
	int binHigh = 0;

	//TCanvas *c = new TCanvas("c","c",1000,700);
	//c->Divide(2,7);

	for(int i=0; i<NPT; i++)
	{
		//For He+Au histograms extract 2D correlation function for each pT bin
		binLow = hdEta_dPhi_pT->GetZaxis()->FindBin(ptl[i]);
		binHigh = hdEta_dPhi_pT->GetZaxis()->FindBin(pth[i]);

		hdEta_dPhi_pT->GetZaxis()->SetRange(binLow,binHigh);

		hdEta_dPhi[i] = (TH2D*) hdEta_dPhi_pT->Project3D("yx");
		hdEta_dPhi[i]->SetName(Form("hdEta_dPhi_%i",i));
	}

	binLow = 0;
	binHigh = 0;

	for(int i=0; i<NPT; i++)
	{
		//c->cd(i+1);

		//For p+p histograms extract 2D correlation function for each pT bin
		binLow = hdEta_dPhi_pT_pp->GetZaxis()->FindBin(ptl[i]);
		binHigh = hdEta_dPhi_pT_pp->GetZaxis()->FindBin(pth[i]);

		hdEta_dPhi_pT_pp->GetZaxis()->SetRange(binLow,binHigh);

		hdEta_dPhi_pp[i] = (TH2D*) hdEta_dPhi_pT_pp->Project3D("yx");
		hdEta_dPhi_pp[i]->SetName(Form("hdEta_dPhi_pp_%i",i));
		//hdEta_dPhi_pp[i]->Draw("COLZ");
	}
}

void normalizeInputHistograms()
{
	int binLow = 0;
	int binHigh = 0;

	int binEtaLow = 0;
	int binEtaHigh = 0;

	float nTrig = 0;

	for(int i=0; i<NPT; i++)
	{
		binLow = hEta_pT->GetYaxis()->FindBin(ptl[i]);
		binHigh = hEta_pT->GetYaxis()->FindBin(pth[i]);

		binEtaLow = hEta_pT->GetXaxis()->FindBin(-1*etaLim);
		binEtaHigh = hEta_pT->GetXaxis()->FindBin(etaLim);

		TH1F *hEta = (TH1F*) hEta_pT->ProjectionX(Form("hEta_%i",i),binLow,binHigh);
		nTrig = hEta->Integral(binEtaLow,binEtaHigh);

		hdEta_dPhi[i]->Scale(1/nTrig);
		hdEta_dPhi[i]->Scale(1/dEta);
		hdEta_dPhi[i]->Scale(1/dPhi);
	}

	binLow = 0;
	binHigh = 0;
	nTrig = 0;

	for(int i=0; i<NPT; i++)
	{
		binLow = hEta_pT_pp->GetYaxis()->FindBin(ptl[i]);
		binHigh = hEta_pT_pp->GetYaxis()->FindBin(pth[i]);

		TH1F *hEta = (TH1F*) hEta_pT_pp->ProjectionX(Form("hEta_pp_%i",i),binLow,binHigh);
		nTrig = hEta->Integral();

		hdEta_dPhi_pp[i]->Scale(1/nTrig);
		hdEta_dPhi_pp[i]->Scale(1/dEta);
		hdEta_dPhi_pp[i]->Scale(1/dPhi);

	}
}

void integrateCorrelation()
{
	//Integrate correlation functions to extract 1D long-range azimuthal per-trigger yields
	int xBinLo = 0;
	int xBinHi = 0;

	for(int i=0; i<NPT; i++)
	{
		xBinLo = hdEta_dPhi[i]->GetXaxis()->FindBin(deltaEta_lo);
		xBinHi = hdEta_dPhi[i]->GetXaxis()->FindBin(deltaEta_hi);

		hCorrelation_He3Au[i] = (TH1F*) hdEta_dPhi[i]->ProjectionY(Form("hCorrelation_He3Au_%i",i),xBinLo,xBinHi);
		hCorrelation_He3Au[i]->SetBinContent(1,hCorrelation_He3Au[i]->GetBinContent(2));
	}

	for(int i=0; i<NPT; i++)
	{
		xBinLo = hdEta_dPhi_pp[i]->GetXaxis()->FindBin(deltaEta_lo);
		xBinHi = hdEta_dPhi_pp[i]->GetXaxis()->FindBin(deltaEta_hi);

		hCorrelation_pp[i] = (TH1F*) hdEta_dPhi_pp[i]->ProjectionY(Form("hCorrelation_pp_%i",i),xBinLo,xBinHi);
		hCorrelation_pp[i]->SetBinContent(1,hCorrelation_pp[i]->GetBinContent(2));
	}
}

void drawIntegratedCorrelations()
{
	gStyle->SetOptStat(0);
	TCanvas *integratedCorr = new TCanvas("integratedCorr","Raw Correlation Functions He3+Au",1200,600);
	integratedCorr->Divide(5,2);

	for(int i=0; i<NPT; i++)
	{
		integratedCorr->cd(i+1);
		TF1 *f = (TF1*) hCorrelation_He3Au[i]->GetFunction(Form("ZYAMFit_%i",i));
		hCorrelation_He3Au[i]->Draw();
		f->Draw("same");
	}

	TCanvas *integratedCorr_pp = new TCanvas("integratedCorr_pp","Raw Correlation Functions p+p",1200,600);
	integratedCorr_pp->Divide(5,2);

	for(int i=0; i<NPT; i++)
	{
		integratedCorr_pp->cd(i+1);
		hCorrelation_pp[i]->SetLineColor(kBlack);
		hCorrelation_pp[i]->Draw();
	}

	TCanvas *integratedCorr_noZYAM = new TCanvas("integratedCorr_noZYAM","integratedCorr_noZYAM",1200,600);
	integratedCorr_noZYAM->Divide(5,2);

	for(int i=0; i<NPT; i++)
	{
		integratedCorr_noZYAM->cd(i+1);
		hCorrelation_NoZYAM_He3Au[i]->Draw();
		hCorrelation_NoZYAM_pp[i]->SetLineColor(kRed);
		hCorrelation_NoZYAM_pp[i]->Draw("same");
	}

	TCanvas *integratedCorr_noFrag = new TCanvas("integratedCorr_noFrag","Correlation Function - pp Background Removed",1200,600);
	integratedCorr_noFrag->Divide(5,2);

	for(int i=0; i<NPT; i++)
	{
		integratedCorr_noFrag->cd(i+1);
		hCorrelation_NoZYAM_NoFrag_He3Au[i]->SetLineColor(kOrange+7);
		hCorrelation_NoZYAM_NoFrag_He3Au[i]->Draw();

		//Draw Fourier reconstructions
		//reconstruction[i]->Draw("same");
	}
}

void computeBZYAM()
{
	//Compute bZYAM by fitting with a polynomial
	for(int i=0; i<NPT; i++)
	{
		ZYAMFit[i] = new TF1(Form("ZYAMFit_%i",i),"pol8",-1*TMath::Pi()/2,1.5*TMath::Pi());
		//ZYAMFit[i] = new TF1(Form("ZYAMFit_%i",i),"pol2",0.5,2);
		if(hCorrelation_He3Au[i]->GetEntries() > 0)
		{
			hCorrelation_He3Au[i]->Fit(Form("ZYAMFit_%i",i),"Q0R");
			TF1 *f = (TF1*) hCorrelation_He3Au[i]->GetFunction(Form("ZYAMFit_%i",i));
			bZYAM[i] = ZYAMFit[i]->GetMinimum(0,TMath::Pi());
		}
		else
		{
			bZYAM[i] = 0;
		}

		ZYAMFit_pp[i] = new TF1(Form("ZYAMFit_pp_%i",i),"pol8",-1*TMath::Pi()/2,1.5*TMath::Pi());
		if(hCorrelation_pp[i]->GetEntries() > 0)
		{
			hCorrelation_pp[i]->Fit(Form("ZYAMFit_pp_%i",i),"Q0R");
			TF1 *g = (TF1*) hCorrelation_pp[i]->GetFunction(Form("ZYAMFit_pp_%i",i));
			bZYAM_pp[i] = ZYAMFit_pp[i]->GetMinimum(-1*TMath::Pi()/2,1.5*TMath::Pi());
			cout <<"BZYAM PP "<< bZYAM_pp[i] << endl;
		}
		else
		{
			bZYAM_pp[i] = 0;
		}
	}
}

void removeBZYAM()
{
	int nbins = 0;
	float cont = 0;
	for(int i=0; i<NPT; i++)
	{
		hCorrelation_NoZYAM_He3Au[i] = (TH1F*) hCorrelation_He3Au[i]->Clone(Form("hCorrelation_NoZYAM_He3Au_%i",i));
		hCorrelation_NoZYAM_He3Au[i]->Reset();

		nbins = hCorrelation_He3Au[i]->GetNbinsX();
		for(int j=1; j<=nbins; j++)
		{
			cont = hCorrelation_He3Au[i]->GetBinContent(j);
			hCorrelation_NoZYAM_He3Au[i]->SetBinContent(j,cont-bZYAM[i]);
		}
	}

	for(int i=0; i<NPT; i++)
		{
			hCorrelation_NoZYAM_pp[i] = (TH1F*) hCorrelation_pp[i]->Clone(Form("hCorrelation_NoZYAM_pp_%i",i));
			hCorrelation_NoZYAM_pp[i]->Reset();

			nbins = hCorrelation_pp[i]->GetNbinsX();
			for(int j=1; j<=nbins; j++)
			{
				cont = hCorrelation_pp[i]->GetBinContent(j);
				hCorrelation_NoZYAM_pp[i]->SetBinContent(j,cont-bZYAM_pp[i]);
			}
		}
}

void removeFragBrackground()
{
	for(int i=0; i<NPT; i++)
	{
		hCorrelation_NoZYAM_NoFrag_He3Au[i] = (TH1F*) hCorrelation_NoZYAM_He3Au[i]->Clone(Form("hCorrelation_NoZYAM_NoFrag_He3Au_%i",i));
		hCorrelation_NoZYAM_NoFrag_He3Au[i]->Add(hCorrelation_NoZYAM_pp[i],-1);
	}
}

void An(int n, const TH1* h, double &an, double &an_err, TString opt)
{
  // Extract fourier coefficient n directly from h. If opt is "hi" or
  // "lo", the coefficient is calculated for the points + or - their
  // errors. "hi" and "lo" are meant for an h with systematic (not
  // statistical) errors, so an_err is not calculated.
  // Code from Andrew Adare
  // https://www.phenix.bnl.gov/viewvc/viewvc.cgi/phenix/offline/analysis/sickles_dAu/macros/DphiDFTAnalysis.C?revision=1.9&view=co
  an = an_err = 0.;

  double period = h->GetNbinsX();

  // Calculate and set coeff. an
  for (int ib=1; ib<=h->GetNbinsX(); ib++) {
    double deltaphi = h->GetBinCenter(ib);
    double weight   = h->GetBinContent(ib);

    if (opt.Contains("hi"))
      weight += h->GetBinError(ib);
    else if (opt.Contains("lo"))
      weight -= h->GetBinError(ib);

    an += cos(n * deltaphi) * weight / period;
  }

  if (opt.Contains("hi") || opt.Contains("lo"))
    return;

  // Statistical uncertainty an_err
  for (int ib=1; ib<=h->GetNbinsX(); ib++) {
    double dphi = h->GetBinCenter(ib);
    double dw   = h->GetBinError(ib);
    double e    = dw * cos(n * dphi ) / period;
    an_err += e*e;
  }
  an_err = TMath::Sqrt(an_err);

  return;
}

void extractFourierCoefficients(TH1D *hProjection, float bZYAM, int pTBin)
{
	double a0 = hProjection->Integral()/hProjection->GetNbinsX();
	double an, an_err,cn,cn_err,vn,vn_err;

	const int n = 4; //Number of orders of anisotropy to compute
	int ccColor[n] = {kBlack, kRed, kBlue, kGreen};

	TF1 *f_an[n];

	an_coeffs[pTBin][0] = a0;

	for(int i=0; i<n; i++)
	{
		an = 0;
		an_err = 0;
		cn = 0;
		cn_err = 0;

		An(i+1,hProjection,an,an_err,"");

		cn = an/(bZYAM + a0);
		cn_err = an_err/(bZYAM + a0);

		if(cn >= 0)
		{
			vn = TMath::Sqrt(cn);
			vn_err = 0.5*cn_err;
		}
		else
		{
			vn = -9999;
			vn_err = -9999;
		}

		an_coeffs[pTBin][i+1] = an;
		vn_coeffs[pTBin][i+1] = vn;
		vn_err_coeffs[pTBin][i+1] = vn_err;

		cout << endl;
		cout << "pT Bin = " << pTBin << endl;
		cout << "a_0 = " << a0 << endl;
		cout << "a_" << i+1 << " = " << an << " pm " << an_err << endl;
		cout << "c_" << i+1 << " = " << cn << " pm " << cn_err << endl;
		cout << "v_" << i+1 << " = " << vn << " pm " << vn_err << endl;
	    cout << endl;
	}
}

void drawFlowMoments()
{
	//Arrays for TGraphErrors
	double epsilon = 0.007;

	double pT[NPT] = {0.27,0.41,0.55,0.69,0.83,0.97,1.11,1.25,1.39,1.53};
	double pT_plus_epsilon[NPT] = {0.27+epsilon,0.41+epsilon,0.55+epsilon,0.69+epsilon,0.83+epsilon,0.97+epsilon,1.11+epsilon,1.25+epsilon,1.39+epsilon,1.53+epsilon};//,1.6625+epsilon,1.8875+epsilon};
	double pT_minus_epsilon[NPT] = {0.27-epsilon,0.41-epsilon,0.55-epsilon,0.69-epsilon,0.83-epsilon,0.97-epsilon,1.11-epsilon,1.25-epsilon,1.39-epsilon,1.53-epsilon};//,1.6625-epsilon,1.8875-epsilon};

	double ex[NPT] = {0};

	double v2[NPT] = {0};
	double v3[NPT] = {0};
	double v4[NPT] = {0};

	double e_v2[NPT] = {0};
	double e_v3[NPT] = {0};
	double e_v4[NPT] = {0};

	for(int i=0; i<NPT; i++)
	{
		v2[i] = vn_coeffs[i][2];
		v3[i] = vn_coeffs[i][3];
		v4[i] = vn_coeffs[i][4];

		e_v2[i] = vn_err_coeffs[i][2];
		e_v3[i] = vn_err_coeffs[i][3];
		e_v4[i] = vn_err_coeffs[i][4];

		cout << "Drawing - pT Bin " << i << endl;
		cout << "v2 " << v2[i] << endl;
		cout << "v3 " << v3[i] << endl << endl;
	}

	//TGraphs
	gv2 = new TGraphErrors(NPT,pT,v2,ex,e_v2);
	gv2->SetMarkerStyle(20);
	gv2->SetLineColor(kBlue);

	gv3 = new TGraphErrors(NPT,pT,v3,ex,e_v3);
	gv3->SetMarkerStyle(20);
	gv3->SetLineColor(kRed);

	gv4 = new TGraphErrors(NPT,pT,v4,ex,e_v4);
	gv4->SetMarkerStyle(20);
	gv4->SetLineColor(kOrange+7);

	TCanvas *cFourier = new TCanvas("cFourier","cFourier",700,600);
	gv2->Draw("ALP");
	gv3->Draw("LP,same");
	//gv4->Draw("LP,same");
}

void writeToFile()
{
	string outname = "/Users/jdok/Documents/CUBoulder/JNLab/Spring15/He3_AMPT/FlowData/AMPT_"+collSyst+"_sigma75_vn.root";
	TFile *fOut = new TFile(outname.c_str(),"RECREATE");
	gv2->SetName("gv2");
	gv2->Write();
	gv3->SetName("gv3");
	gv3->Write();
}

void drawSelectedCorrelations()
{
	gStyle->SetOptStat(0);
	int selectedBin = 7;
	TCanvas *cCorrSelected = new TCanvas("cCorrSelected","cCorrSelected",500,500);
	hCorrelation_NoZYAM_He3Au[selectedBin]->SetYTitle("#frac{1}{N_{Trig}} #frac{dN}{d#Delta#phi} - b_{ZYAM}");
	hCorrelation_NoZYAM_He3Au[selectedBin]->SetXTitle("#Delta#phi [rad]");
	hCorrelation_NoZYAM_He3Au[selectedBin]->Draw();
	hCorrelation_NoZYAM_pp[selectedBin]->SetLineColor(kRed);
	hCorrelation_NoZYAM_pp[selectedBin]->Draw("same");

	TLegend *tl = new TLegend(0.3,0.6,0.5,0.8);
	tl->AddEntry(hCorrelation_NoZYAM_He3Au[selectedBin],"^{3}He+Au Central","L");
	tl->AddEntry(hCorrelation_NoZYAM_pp[selectedBin],"p+p","L");
	//tl->Draw("same");

	TText *tPanel1 = new TText(0.03,0.03,"(a)");
	//tPanel1->Draw("same");

	TCanvas *cCorrSelected2 = new TCanvas("cCorrSelected2","cCorrSelected2",500,500);
	hCorrelation_NoZYAM_NoFrag_He3Au[selectedBin]->SetYTitle("#frac{1}{N_{Trig}} #frac{dN}{d#Delta#phi} - b_{ZYAM}");
	hCorrelation_NoZYAM_NoFrag_He3Au[selectedBin]->SetXTitle("#Delta#phi [rad]");
	hCorrelation_NoZYAM_NoFrag_He3Au[selectedBin]->Draw();

	TText *tPanel2 = new TText(0.03,0.03,"(b)");
	//tPanel2->Draw("same");
}

void computeFourierReconstruction()
{
	for(int i=0; i<NPT; i++)
	{
		float a0 = an_coeffs[i][0];
		float a1 = an_coeffs[i][1];
		float a2 = an_coeffs[i][2];
		float a3 = an_coeffs[i][3];
		float a4 = an_coeffs[i][4];

		reconstruction[i] = new TF1(Form("reconstruction_%i",i), "[0] + 2*[1]*TMath::Cos(x) + 2*[2]*TMath::Cos(2*x) + 2*[3]*TMath::Cos(3*x) + 2*[4]*TMath::Cos(4*x)",-TMath::PiOver2(),3*TMath::PiOver2());
		reconstruction[i]->SetParameters(a0,a1,a2,a3,a4);
	}
}

void FlowCoefficientsAnalyzer()
{
	//Read histograms from file
	readHistogramsFromFile();

	//Prepare histograms for analysis
	prepareInputHistograms();

	//Normalize by the number of trigger particles
	normalizeInputHistograms();

	//Integrate correlations over eta range
	integrateCorrelation();

	//Compute bZYAM
	computeBZYAM();

	//Remove bZYAM
	removeBZYAM();

	//Remove background from jet fragmentation
	removeFragBrackground();

	//Compute Fourier coefficients
	for(int i=0; i<NPT; i++)
	{
		extractFourierCoefficients((TH1D*) hCorrelation_NoZYAM_NoFrag_He3Au[i],bZYAM[i],i);
	}

	//Check Fourier reconstruction
	computeFourierReconstruction();

	//Draw flow moments
	drawFlowMoments();

	//Draw correlation functions
	drawIntegratedCorrelations();

	//Draw correlations for selected pT bin
	//drawSelectedCorrelations();

	//Calculate invariant yield
	//computeInvariantYield();

	//Charged particle multiplicity
	//computeChargedPartMult();

	//Write to file
	//writeToFile();
}

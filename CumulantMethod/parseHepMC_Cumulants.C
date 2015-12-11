//------------------------------------------
// Code to test the implementation of
// two- and four-particle cumulants on
// output from the AMPT generator
//------------------------------------------

#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TFile.h"
#include "TMath.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TProfile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

using namespace std;

//-------------------------------------------
// Variables
//-------------------------------------------
//Number of bins in histograms
const int NOB = 100;

//Data structure to contain particles
struct particle
{
  int id;
  float px;
  float py;
  float pz;
  float x;
  float y;
  float z;
  float phi;
  float rsquare;
  float pT;
  float eta;
};

vector<particle> finalparticles;

//Event characterization variables
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

//Event-wise reference flow
//Must be cleared after every event
float q2x  = 0;
float q2y  = 0;
float q4x  = 0;
float q4y  = 0;
float sb_2 = 0;
float sb_4 = 0;

TH1F *h_sb_2prime;
TH1F *h_sb_4prime;

TH1F *hp2x;
TH1F *hp2y;
TH1F *hp4x;
TH1F *hp4y;

TH1F *hn;

//Average over events
TProfile *h_db_2;
TProfile *h_db_4;
TProfile *h_db_2prime;
TProfile *h_db_4prime;

TProfile *h_db_2_etagap;
TProfile *h_db_2prime_etagap;

//Diagnostic histograms
TH1F *hEta;
TH1F *hpT;
TH1F *hpT_EtaCut;
TH1F *hPhi;

int evtnumber = 0;

float etaGap = 0.75;

const int NBIN = 50;

//-------------------------------------------
// Functions
//-------------------------------------------

float calc4_event(float Xn, float Yn, float X2n, float Y2n, float M)
{

  float Qn2 = Xn * Xn + Yn * Yn;
  float Qn2d = Xn * Xn - Yn * Yn;

  float one   = Qn2 * Qn2;
  float two   = X2n * X2n + Y2n * Y2n;
  float three = (2 * (X2n * Qn2d + 2 * Y2n * Xn * Yn));
  float four  = 2 * (2 * (M - 2) * Qn2);
  float five  = 2 * M * (M - 3);

  float numerator = one + two - three - four + five;
  float denominator = M * (M - 1) * (M - 2) * (M - 3);

  return numerator / denominator;

}

float calc4_track(float xn, float yn, float x2n, float y2n, float Xn, float Yn, float X2n, float Y2n, float M)
{

  float one   = (xn * Xn + yn * Yn) * (Xn * Xn + Yn * Yn);
  float two   = x2n * Xn * Xn - x2n * Yn * Yn + 2 * y2n * Xn * Yn;
  float three = xn * Xn * X2n + xn * Yn * Y2n - yn * (X2n * Yn - Xn * Y2n);
  float four  = 2 * M * (xn * Xn + yn * Yn);
  float five  = 2 * (Xn * Xn + Yn * Yn);
  float six   = 7 * (xn * Xn + yn * Yn);
  float seven = xn * Xn + yn * Yn;
  float eight = x2n * X2n + y2n * Y2n;
  float nine = 2 * (xn * Xn + yn * Yn);

  float numerator = one - two - three - four - five + six - seven + eight + nine + 2 * M - 6;
  float denominator = (M - 1) * (M - 2) * (M - 3);

  return numerator / denominator;

}

void computeFlow()
{
  float c2_2        = h_db_2->GetBinContent(1);
  float c2_2_etagap = h_db_2_etagap->GetBinContent(1);
  float c2_4        = h_db_4->GetBinContent(1) - 2 * pow(h_db_2->GetBinContent(1), 2);

  cout << "c2_2    " << c2_2 << endl;
  cout << "c2_2 eg " << c2_2_etagap << endl;
  cout << "c2_4    " << c2_4 << endl;

  TProfile *hd2_2 = (TProfile*) h_db_2prime->Clone("hd2_2");
  hd2_2->SetTitle("Two Particle Cumulants");
  hd2_2->Scale(1.0 / sqrt(c2_2));

  TProfile *hd2_2_etagap = (TProfile*) h_db_2prime_etagap->Clone("hd2_2_etagap");
  hd2_2_etagap->SetTitle("Two Particle Cumulanta with Eta Gap");
  hd2_2_etagap->Scale(1.0/sqrt(c2_2_etagap));

  TH1D *h_proj_db_4prime = h_db_4prime->ProjectionX("h_proj_db_4prime", "E");
  TH1D *h_proj_db_2prime = h_db_2prime->ProjectionX("h_proj_db_2prime", "E");

  h_proj_db_2prime->Scale(2 * h_db_2->GetBinContent(1));
  TH1D *hd2_4 = (TH1D*) h_proj_db_4prime->Clone("hd2_4");
  hd2_4->Add(h_proj_db_2prime, -1);

  hd2_4->Scale(-1.0 / pow(-1 * c2_4, 0.75));

  hd2_2_etagap->Draw();
  hd2_2->SetLineColor(kRed);
  hd2_2->Draw("same");
}

void writeHistosToFile(char *outFileName = "")
{
  TFile* fout = new TFile(outFileName, "RECREATE");
  h_db_2prime->Write();
  h_db_4prime->Write();
  h_db_2->Write();
  h_db_4->Write();
  fout->Close();
}

void processEvent()
{
  int mult = finalparticles.size();
  int mult_A = 0;
  int mult_B = 0;

  // --- first track loop, q-vectors
  float Q2x = 0;
  float Q2y = 0;
  float Q4x = 0;
  float Q4y = 0;

  float Q2x_A = 0;
  float Q2y_A = 0;
  float Q4x_A = 0;
  float Q4y_A = 0;

  float Q2x_B = 0;
  float Q2y_B = 0;
  float Q4x_B = 0;
  float Q4y_B = 0;

  for (int itrk = 0; itrk < mult; itrk++)
  {
    float phi = finalparticles[itrk].phi;
    float eta = finalparticles[itrk].eta;

    Q2x += cos(2 * phi);
    Q2y += sin(2 * phi);
    Q4x += cos(4 * phi);
    Q4y += sin(4 * phi);

    if (eta < -1 * etaGap)
    {
      Q2x_A += cos(2 * phi);
      Q2y_A += sin(2 * phi);
      Q4x_A += cos(4 * phi);
      Q4y_A += sin(4 * phi);

      mult_A++;
    }
    else if (eta > etaGap)
    {
      Q2x_B += cos(2 * phi);
      Q2y_B += sin(2 * phi);
      Q4x_B += cos(4 * phi);
      Q4y_B += sin(4 * phi);

      mult_B++;
    }

  } // End of track loop

  float two = ( Q2x * Q2x + Q2y * Q2y - mult) / (mult * mult - mult);
  float two_etagap = (Q2x_A * Q2x_B + Q2y_A * Q2y_B) / (mult_A * mult_B);
  float four = calc4_event(Q2x, Q2y, Q4x, Q4y, mult);

  h_db_2->Fill(0.5, two);
  h_db_4->Fill(0.5, four);
  h_db_2_etagap->Fill(0.5, two_etagap);

  for (int itrk = 0; itrk < mult; itrk++)
  {
    particle p = finalparticles[itrk];
    float phi = p.phi;
    float pt = p.pT;
    float eta = p.eta;

    float u2x = cos(2 * phi);
    float u2y = sin(2 * phi);

    float twoprime = ( u2x * Q2x + u2y * Q2y - 1) / (mult - 1);
    h_db_2prime->Fill(pt, twoprime);

    float u4x = cos(4 * phi);
    float u4y = sin(4 * phi);

    float fourprime = calc4_track(u2x, u2y, u4x, u4y, Q2x, Q2y, Q4x, Q4y, mult);
    h_db_4prime->Fill(pt, fourprime);

    if (eta < -1 * etaGap)
    {
      float u2x_A = cos(2 * phi);
      float u2y_A = sin(2 * phi);

      float twoprime_etagap = (u2x_A * Q2x_B + u2y_A * Q2y_B) / mult_B;
      h_db_2prime_etagap->Fill(pt, twoprime_etagap);
    }

  } // End of track loop
}

void parseHepMC_Cumulants(char *filename = "/direct/phenix+hhj2/jdok/AMPT_PhobosGlauber/ana/ampt.dat", char* outFileName = "ampt_correlation_out.root")
{
  h_db_2      = new TProfile("h_db_2", "h_db_2", 1, 0, 1, -500, 500);
  h_db_4      = new TProfile("h_db_4", "h_db_4", 1, 0, 1, -500, 500);
  h_db_2prime = new TProfile("h_db_2prime", "h_db_2prime", NBIN, 0, 5, -500, 500);
  h_db_4prime = new TProfile("h_db_4prime", "h_db_4prime", NBIN, 0, 5, -500, 500);

  h_db_2_etagap = new TProfile("h_db_2_etagap", "h_db_2_etagap", 1, 0, 1, -500, 500);
  h_db_2prime_etagap = new TProfile("h_db_2prime_etagap", "h_db_2prime_etagap", NBIN, 0, 5, -500, 500);

  h_sb_2prime = new TH1F("h_sb_2prime", "h_sb_2prime", NBIN, 0, 5);
  h_sb_4prime = new TH1F("h_sb_4prime", "h_sb_4prime", NBIN, 0, 5);
  hn          = new TH1F("hn", "hn", NBIN, 0, 5);
  hp2x        = new TH1F("hp2x", "hp2x", NBIN, 0, 5);
  hp2y        = new TH1F("hp2y", "hp2y", NBIN, 0, 5);
  hp4x        = new TH1F("hp4x", "hp4x", NBIN, 0, 5);
  hp4y        = new TH1F("hp4y", "hp4y", NBIN, 0, 5);

  hEta         = new TH1F("hEta", "hEta;#eta;Counts", NOB, -6, 6);
  hpT          = new TH1F("hpT", "hpT;p_{T};Counts", NOB, 0, 10);
  hpT_EtaCut   = new TH1F("hpT_EtaCut", "hpT_EtaCut;p_{T};Counts", NOB, 0, 20);
  hPhi         = new TH1F("hPhi", "hPhi;#phi;Counts", NOB, -TMath::Pi(), TMath::Pi());

  for (int fnum = 0; fnum < 15; fnum++)
  {
    cout << "OPEN NEW FILE " << fnum << endl;
    ifstream myFile;
    myFile.open(Form("/direct/phenix+hhj2/jdok/AMPT_dAu_200_5M/ampt_%i.dat", fnum));
    //myFile.open(filename);

    if (!myFile)
    {
      printf("Input file does not exist!\n");
      continue;
    }
    else
    {
      cout << "--> Successfully opened file " << filename << endl << endl;
      evtnumber = 0;
    }

    //--------------------------------------
    // Run over events
    //--------------------------------------

    while (myFile)
    {
      if (evtnumber % 500 == 0)
      {
        cout << "-----> // // Reading Event No. " << evtnumber << endl;
      }

      //Read the event header
      int testnum;
      int nlist; //Number of finalparticles in output list
      double impactpar;
      int npartproj;
      int nparttarg;
      int npartprojelas;
      int npartprojinelas;
      int nparttargelas;
      int nparttarginelas;
      double junk;

      myFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

      //Avoid double reading the last entry
      if (!myFile) break;

      //Loop over each particle in the event
      for (int i = 0; i < nlist; i++)
      {
        TMCParticle *part = new TMCParticle();

        int partid;
        float pvec[3];
        float mass;
        double spacetime[4];

        myFile >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

        float energy =  sqrt (pow(pvec[0], 2) + pow(pvec[1], 2) + pow(pvec[2], 2) + pow(mass, 2));

        part->SetEnergy(energy);
        part->SetKF(partid);
        part->SetMass(mass);
        part->SetPx(pvec[0]);
        part->SetPy(pvec[1]);
        part->SetPz(pvec[2]);

        TLorentzVector ev(pvec[0], pvec[1], pvec[2], energy);

        if (ev.Pt() < 0.001)
        {
          continue;
        }

        int KF = TMath::Abs(partid);

        if (!(KF == 211 || KF == 321 || KF == 2212)) //Pions, kaons and protons
        {
          continue;
        }

        if (TMath::Abs(ev.PseudoRapidity()) > 2) continue;

        particle p;
        p.px = pvec[0];
        p.py = pvec[1];
        p.pz = pvec[2];
        p.pT = ev.Pt();
        p.phi = ev.Phi();
        p.eta = ev.PseudoRapidity();

        finalparticles.push_back(p);

      } // End loop over finalparticles

      //Do cumulant analysis here only with high multiplicity events
      if (finalparticles.size() > 50)
      {
        processEvent();
      }

      finalparticles.clear();
    }

    myFile.close();
  }

  computeFlow();
  //writeHistosToFile(outFileName);
}

  //------------------------------------------
  // Code to test the implementation of 
  // two- and four-particle cumulants on
  // output from the AMPT generator
  //------------------------------------------

  #include "TLorentzVector.h"
  #include "TMCParticle.h"
  #include "TFile.h"
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
vector<TMCParticle*> finalparticles;

  //Complex number 
struct cmpx
{
  float re;
  float im;
};

  //Number of bins in histograms
const int NOB = 100;

  //Event characterization variables
int npart          = 0;
int npartsum       = 0;
int nspectator     = 0;
int numevent       = 0;
int process_number = 0;

  //Variables for reference four-particle cumulant analysis
float q2x  = 0;
float q2y  = 0;
float q4x  = 0;
float q4y  = 0;
int n      = 0;
float cn_2 = 0;
float cn_4 = 0;
float db_2 = 0;
float db_4 = 0;

  //Variables for differential four-particle cumulant analysis
TH1F *hpT_Template;
TProfile *hdb_2prime;
TProfile *hdb_4prime;
TProfile *hdb_2;
TProfile *hdb_4;
TProfile *hd2_2;
TH1F *hd2_4;
TH1F *hv2_4;

int evtnumber = 0;

  //-------------------------------------------
  // Functions
  //-------------------------------------------

void writeHistosToFile(char *outFileName="")
{
  TFile* fout = new TFile(outFileName, "RECREATE");
  //hv2_pT->Write();
  //hv3_pT->Write();
  //hNcoll_Yield->Write();
}

cmpx multiplyComplex(float x, float y, float u, float v)
{
  //z1 = x + iy
  //z2 = u + iv
  cmpx c;
  c.re = x*u -y*v;
  c.im = x*v + y*u;

  return c;
}

float getNormSq(float x, float y)
{
  //z = x + iy
  return x*x + y*y;
}

void computeCumulants()
{
  db_2 = hdb_2->GetBinContent(1);
  db_4 = hdb_4->GetBinContent(1);

  cn_2 = db_2;
  cn_4 = db_4 - 2*pow(db_2,2);

  for(int i=1; i<=hpT_Template->GetNbinsX(); i++)
  {
    float cont = (hdb_4prime->GetBinContent(i)) - 2*(hdb_2prime->GetBinContent(i))*db_2;
    hd2_4->SetBinContent(i, cont);
  }

  cout << "<<2>> = " << db_2 << endl;
  cout << "<<4>> = " << db_4 << endl;

  //Two-particle cumulant result
  hd2_2 = (TProfile*) hdb_2prime->Clone("hd2_2");
  hd2_2->Scale(1.0/sqrt(cn_2));
  //hd2_2->Draw();

  //Four-particle cumulant result
  hv2_4 = (TH1F*)hd2_4->Clone("hv2_4");
  float scaleFactor = -1*pow(-1*cn_4,3.0/4.0);
  hv2_4->Scale(1.0/scaleFactor);
  hv2_4->Draw();
}

void processEventCumulants()
{
  //Initialize variables
  n   = 0;
  q2x = 0;
  q2y = 0;
  q4x = 0;
  q4y = 0;

    //Compute reference flow for the event at hand
  for(int i=0; i<finalparticles.size(); i++)
  {
    TMCParticle *part = finalparticles[i];
    float px = part->GetPx();
    float py = part->GetPy();
    float pz = part->GetPz();
    float energy = part->GetEnergy();
    TLorentzVector ev(px,py,pz,energy);
    double pT = ev.Pt();
    double phi = ev.Phi();
    double eta = ev.Eta();

      //Only consider particles within -2 < eta < 2 and 0 < pT [GeV/c] < 5
    if(TMath::Abs(eta) > 2.0 || pT > 5.0) continue;

    q2x = q2x + TMath::Cos(2*phi);
    q2y = q2y + TMath::Sin(2*phi);

    q4x = q4x + TMath::Cos(4*phi);
    q4y = q4y + TMath::Sin(4*phi);

    n++;
  }

  float sb_2 = (getNormSq(q2x,q2y) - n)/(n*(n-1));

  cmpx c1 = multiplyComplex(q2x, -1*q2y, q2x, -1*q2y); 
  cmpx c2 = multiplyComplex(q4x, q4y, c1.re, c1.im);
  float sb_4 = (getNormSq(q2x,q2y)*getNormSq(q2x,q2y) + getNormSq(q4x,q4y) - 2*c2.re - 4*(n-2)*getNormSq(q2x,q2y) + 2*n*(n-3))/(n*(n-1)*(n-2)*(n-3));

  //cout << "<2> = " << sb_2 << endl;
  //cout << "<4> = " << sb_4 << endl << endl;

  hdb_2->Fill(0.5,sb_2,1);
  hdb_4->Fill(0.5,sb_4,1);

  //Compute differential flow for the event at hand
  TH1F *hp2x         = (TH1F*) hpT_Template->Clone("hp2x");
  TH1F *hp2y         = (TH1F*) hpT_Template->Clone("hp2y");
  TH1F *hp4x         = (TH1F*) hpT_Template->Clone("hp4x");
  TH1F *hp4y         = (TH1F*) hpT_Template->Clone("hp4y");
  TH1F *hsb_2prime   = (TH1F*) hpT_Template->Clone("hsb_2prime");
  TH1F *hsb_4prime   = (TH1F*) hpT_Template->Clone("hsb_4prime");

  float m[NOB] = {0};

  for(int i=0; i<finalparticles.size(); i++)
  {
    TMCParticle *part = finalparticles[i];
    float px = part->GetPx();
    float py = part->GetPy();
    float pz = part->GetPz();
    float energy = part->GetEnergy();
    TLorentzVector ev(px,py,pz,energy);
    double pT = ev.Pt();
    double phi = ev.Phi();
    double eta = ev.Eta();

    if(TMath::Abs(eta) > 2.0 || pT > 5.0) continue;

    int bin = hpT_Template->FindBin(pT);
    hp2x->SetBinContent(bin, hp2x->GetBinContent(bin) + TMath::Cos(2*phi));
    hp2y->SetBinContent(bin, hp2y->GetBinContent(bin) + TMath::Sin(2*phi));
    hp4x->SetBinContent(bin, hp4x->GetBinContent(bin) + TMath::Cos(4*phi));
    hp4y->SetBinContent(bin, hp4y->GetBinContent(bin) + TMath::Sin(4*phi));

    m[bin-1]++;
  }

  for(int i=1; i<= hsb_2prime->GetNbinsX(); i++)
  {
    cmpx aux1 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
    float aux2 = m[i-1]*(n-1); 
    hsb_2prime->SetBinContent(i, (aux1.re-m[i-1])/aux2);
  }

  for(int i=1; i<hdb_4prime->GetNbinsX(); i++)
  {
    cmpx aux1_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, q2y);
    cmpx aux1_1 = multiplyComplex(q2x, -1*q2y, q2x, -1*q2y);
    cmpx aux1_2 = multiplyComplex(aux1_0.re, aux1_0.im, aux1_1.re, aux1_1.im);
    float aux1  = aux1_2.re;

    cmpx aux2_0 = multiplyComplex(hp4x->GetBinContent(i), hp4y->GetBinContent(i), q2x, -1*q2y);
    cmpx aux2_1 = multiplyComplex(aux2_0.re, aux2_0.im, q2x, -1*q2y);
    float aux2  = aux2_1.re;

    cmpx aux3_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, q2y);
    cmpx aux3_1 = multiplyComplex(aux3_0.re, aux3_0.im, q4x, -1*q4y);
    float aux3  = aux2_1.re;

    cmpx aux4_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2y, -1*q2y);
    float aux4  = 2*n*aux4_0.re;

    float aux5  = 2*m[i]*getNormSq(q2x,q2y);

    cmpx aux6_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
    float aux6  = 7*aux6_0.re;

    cmpx aux7_0 = multiplyComplex(q2x, q2y, hp2x->GetBinContent(i), -1*hp2y->GetBinContent(i));
    float aux7  = aux7_0.re;

    cmpx aux8_0 = multiplyComplex(hp4x->GetBinContent(i), hp4y->GetBinContent(i), q4x, -1*q4y);
    float aux8  = aux8_0.re;

    cmpx aux9_0 = multiplyComplex(hp2x->GetBinContent(i), hp2y->GetBinContent(i), q2x, -1*q2y);
    float aux9  = 2*aux9_0.re;

    float aux10 = 2*m[i]*n;
    float aux11 = 6*m[i];
    float aux12 = (m[i]*n-3*m[i])*(n-1)*(n-2);

    hsb_4prime->SetBinContent(i, (aux1-aux2-aux3-aux4-aux5+aux6-aux7+aux8+aux9+aux10-aux11)/(aux12));
  }

  for(int i=1; i<=hpT_Template->GetNbinsX(); i++)
  {
    hdb_2prime->Fill(hsb_2prime->GetBinCenter(i), hsb_2prime->GetBinContent(i),1);
    hdb_4prime->Fill(hsb_4prime->GetBinCenter(i), hsb_4prime->GetBinContent(i),1);
  }
}

void parseHepMC_Cumulants(char *filename = "ampt_dAu_200.dat", char* outFileName = "ampt_correlation_out.root")
{
  //Initialize histograms
  hpT_Template = new TH1F("hpT_Template","hpT_Template",NOB,0,5);
  hdb_2prime = new TProfile("hdb_2prime","hdb_2prime",NOB,0,5,-500,500);
  hdb_4prime = new TProfile("hdb_4prime","hdb_4prime",NOB,0,5,-500,500);
  hd2_4       = new TH1F("hd2_4","hd2_4",NOB,0,5);
  hdb_2      = new TProfile("hdb_2","hdb_2",1,0,1,-100,100);
  hdb_4      = new TProfile("hdb_4","hdb_4",1,0,1,-100,100);

  //Open File
  ifstream myFile;
  myFile.open(filename);

  if (!myFile) 
  {
    printf("Input file does not exist!\n");
    return;
  }
  else
  {
    cout << "--> Successfully opened file " << filename << endl << endl;
  }

  //--------------------------------------
  // Run over events
  //--------------------------------------

  while (myFile) 
  {
    if(evtnumber % 100 == 0)
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
      for (int i=0;i<nlist;i++) 
      {
        TMCParticle *part = new TMCParticle();

        int partid;
        float pvec[3];
        float mass;
        double spacetime[4];

        myFile >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

        float energy =  sqrt (pow(pvec[0],2)+pow(pvec[1],2)+pow(pvec[2],2)+pow(mass,2));

        part->SetEnergy(energy);
        part->SetKF(partid);
        part->SetMass(mass);
        part->SetPx(pvec[0]);
        part->SetPy(pvec[1]);
        part->SetPz(pvec[2]);

        TLorentzVector ev(pvec[0],pvec[1],pvec[2],energy);

        if(ev.Pt() < 0.001)
        {
          continue;
        }

        int KF = TMath::Abs(partid);

    if(!(KF == 211 || KF == 321 || KF == 2212)) //Pions, kaons and protons
    {
      continue;
    }

    finalparticles.push_back(part);  
    
  } // End loop over finalparticles

      //Do cumulant analysis here
  processEventCumulants();
  finalparticles.clear();
}
  computeCumulants();
  //writeHistosToFile(outFileName);
}

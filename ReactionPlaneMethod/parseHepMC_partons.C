#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TFile.h"
#include "TNtuple.h"
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
vector<TMCParticle*> particles;
vector<float> psi_2;
vector<float> psi_3;
vector<float> ncoll;

//Flow coefficients at midrapidity
TH2F *hv2_pT;
TH2F *hv3_pT;

//For study of flow as a function of eta
TH3F *hv2_pT_eta;
TH3F *hv3_pT_eta;

//Charged particle yield
TH1F *hNch_eta;

TH2F *hNcoll_Yield;

int evtnumber = 0;

//-------------------------------------------
// Functions
//-------------------------------------------

void loadNcoll()
{
  ifstream ncoll_file;
  ncoll_file.open("ncoll_file.txt");

  float ncoll_value = 0;
  while(ncoll_file)
    {
      ncoll_file >> ncoll_value;
      ncoll.push_back(ncoll_value);
    }

  ncoll_file.close();
}

void loadEventPlane()
{
  //Psi_2
  ifstream psi_2_file;
  psi_2_file.open("psi_2.txt");

  float psi_2_value = 0;
  while(psi_2_file)
    {
      psi_2_file >> psi_2_value;
      psi_2.push_back(psi_2_value);
    }

  psi_2_file.close();

  //Psi_3
  ifstream psi_3_file;
  psi_3_file.open("psi_3.txt");

  float psi_3_value = 0;
  while(psi_3_file)
    {
      psi_3_file >> psi_3_value;
      psi_3.push_back(psi_3_value);
    }

  psi_3_file.close();
}

void writeHistosToFile(char *outFileName="")
{
  TFile* fout = new TFile(outFileName, "RECREATE");
  hv2_pT->Write();
  hv3_pT->Write();
  //hv2_pT_eta->Write();
  //hv3_pT_eta->Write();
  //hNch_eta->Write();
  hNcoll_Yield->Write();
}

void loopParticles()
{
  //Get event plane info for event at hand
  float psi_2_angle = psi_2[evtnumber-1];
  float psi_3_angle = psi_3[evtnumber-1];
  float event_ncoll = ncoll[evtnumber-1];
  
  for(int ipart = 0; ipart < particles.size(); ipart++)
    {
      //Get the particle information
      TMCParticle *part = particles[ipart];
      float px = part->GetPx();
      float py = part->GetPy();
      float pz = part->GetPz();
      float energy = part->GetEnergy();
      TLorentzVector ev(px,py,pz,energy);
      double pT = ev.Pt();
      double phi = ev.Phi();
      double eta = ev.Eta();

      double v2 = TMath::Cos(2*(phi-psi_2_angle));
      double v3 = TMath::Cos(3*(phi-psi_3_angle));

      hNch_eta->Fill(eta);

      hv2_pT_eta->Fill(v2,pT,eta);
      hv3_pT_eta->Fill(v3,pT,eta);
      
      if(abs(eta) > 2.0)
	{
	  continue;
	}

      hv2_pT->Fill(v2,pT);
      hv3_pT->Fill(v3,pT);
    }

  hNcoll_Yield->Fill(event_ncoll,particles.size());
}

void parseHepMC_partons(char *filename = "zpc.dat", char* outFileName = "ampt_correlation_out.root")
{
  //Initialize histograms
  hv2_pT = new TH2F("hv2_pT","hv2_pT;v_{2};p_{T}",1200,-1.2,1.2,50,0,5);
  hv3_pT = new TH2F("hv3_pT","hv3_pT;v_{3};p_{T}",1200,-1.2,1.2,50,0,5);

  hNcoll_Yield = new TH2F("hNcoll_Yield","hNcoll_Yield;N_{coll};Yield",40,0,40,700,0,700);

  hv2_pT_eta = new TH3F("hv2_pT_eta","hv2_pT_eta",1200,-1.2,1.2,50,0,5,300,-6,6);
  hv3_pT_eta = new TH3F("hv3_pT_eta","hv3_pT_eta",1200,-1.2,1.2,50,0,5,300,-6,6);

  hNch_eta = new TH1F("hNch_eta","hNch_eta",300,-6,6);

  //Load event plane angles for each event
  loadEventPlane();
  //Load ncoll information for each event
  loadNcoll();

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
  //07-13-15: Modified to parse zpc.dat, with parton information at freezeout
  //          Order of columns is changed with respect to ampt.dat  
   
  while (myFile) 
    {
      if(evtnumber % 100 == 0)
	{
	  cout << "-----> // // Reading Event No. " << evtnumber << endl;
	}
      evtnumber++;
           
      //Read the event header
      /*
      int testnum;
      int nlist; //Number of particles in output list
      double impactpar;
      int npartproj;
      int nparttarg;
      int npartprojelas;
      int npartprojinelas;
      int nparttargelas;
      int nparttarginelas;
      double junk; 
      */

      //Read the event header for partons 
      int npartons;
      double impactpar;
      int npartprojelas;
      int npartprojinelas;
      int nparttargelas;
      int nparttarginelas;
      double junk;

      //myFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

      myFile >> junk >> npartons >> impactpar >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas;

      //Avoid double reading the last entry
      if (!myFile) break;

      //Loop over each particle in the event
      for (int i=0;i<npartons;i++) 
	{
	  TMCParticle *part = new TMCParticle();
	  
	  int partid;
	  float pvec[3];
	  float mass;
	  double spacetime[4];
	  
	  //myFile >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

	  myFile >> pvec[0] >> pvec[1] >> pvec[2] >> partid >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];
	 
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

	  /*
	  if(!(KF == 211 || KF == 321 || KF == 2212)) //Pions, kaons and protons
	    {
	      continue;
	    }
	  */

	  particles.push_back(part);  
	  
	} // End loop over particles

      loopParticles();
      particles.clear();
    }
  
  writeHistosToFile(outFileName);
}

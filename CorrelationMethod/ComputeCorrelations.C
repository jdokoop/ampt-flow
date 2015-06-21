//---------------------------------------------------------
// Code to compute two-particle correlations from an
// input file in HEPMC format.
// In this particular version of the code, the trigger
// and associated pT bins are selected to be identical.
//---------------------------------------------------------

#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

//--------------------------------------
// Variables
//--------------------------------------

vector<TMCParticle*> particles;

//Number of pT bins
const int NPT = 10;

//Simulated detector acceptance |eta| < 2
const float etaCut = 2.0;

//3D Histograms
TH3F* hdEta_dPhi_pT;

//2D Histograms
TH2F* hEta_pT;

//1D Histograms
TH1F* hdN_dEta;
TH1F* hNEvents;
TH1F* hNpart;

//Associated and trigger pT bins
float ptl[NPT] = {0.20, 0.34, 0.48, 0.62, 0.76, 0.90, 1.04, 1.18, 1.32, 1.46};
float pth[NPT] = {0.34, 0.48, 0.62, 0.76, 0.90, 1.04, 1.18, 1.32, 1.46, 1.60};
float pTBins[NPT+1] = {0.20, 0.34, 0.48, 0.62, 0.76, 0.90, 1.04, 1.18, 1.32, 1.46, 1.60};

//--------------------------------------
// Functions
//--------------------------------------

/*
 * Function to get the pT bin of a given particle given its transverse momentum
 */
void getPtBin(double val)
{
  for(int i=0; i<NPT; i++)
    {
      if(ptl[i] <= val && val < pth[i])
	{
	  return i;
	}
    }

  return -1;
}

/*
 * Function that loops through particle list to compute two-particle correlations
 */
void computeCorrelation()
{
  for(int ipart1 = 0; ipart1 < particles.size(); ipart1++)
    {
      //Get the particle information
      TMCParticle *part1 = particles[ipart1];
      float px1 = part1->GetPx();
      float py1 = part1->GetPy();
      float pz1 = part1->GetPz();
      float energy1 = part1->GetEnergy();
      TLorentzVector ev1(px1,py1,pz1,energy1);
      float y = ev1.Rapidity();
      double eta1 = ev1.Eta();
      double phi1 = ev1.Phi();
      double pT1 = ev1.Pt();
      int pTBin1 = getPtBin(pT1);

      //Acceptance cut
      if(fabs(eta1) >= etaCut)
	{
	  continue;
	}

      //Reject particles outside of pT range of interest
      if(pTBin1 == -1)
	{
	  continue;
	}

      //Fill pseudorapidity distribution of trigger particles                                                                                      
      hEta_pT->Fill(eta1,pT1);

      //Fill charged particle distribution
      hdN_dEta->Fill(eta1);
   
      for(int ipart2 = 0; ipart2 < particles.size(); ipart2++)
	{
	  //Avoid self-correlation
	  if(ipart1 == ipart2)
	    {
	      continue;
	    }

	  //Get the particle information for the second particle
	  TMCParticle *part2 = particles[ipart2];
	  float px2 = part2->GetPx();
	  float py2 = part2->GetPy();
	  float pz2 = part2->GetPz();
	  float energy2 = part2->GetEnergy();
	  TLorentzVector ev2(px2,py2,pz2,energy2);
	  double eta2 = ev2.Eta();
	  double phi2 = ev2.Phi();
	  double pT2 = ev2.Pt();
	  int pTBin2 = getPtBin(pT2);

	  //Reject particle outside of pT range of interest
	  if(pTBin2 == -1)
	    {
	      continue;
	    }

	  //Acceptance cut
	  if(fabs(eta2) >= etaCut)
	    {
	      continue;
	    }

	  //Only correlate particles in the same bin
	  if(pTBin1 != pTBin2)
	    {
	      continue;
	    }

	  double dPhi = 0;
	  double dEta = eta1-eta2;

	  //Account for wrap-around effects
	  if(phi1-phi2 > 1.5*TMath::Pi())
	    {
	      dPhi = phi1-phi2-2*TMath::Pi();
	    }
	  else if(phi1-phi2 < -0.5*TMath::Pi())
	    {
	      dPhi = phi1-phi2+2*TMath::Pi();
	    }
	  else
	    {
	      dPhi = phi1-phi2;
	    }

	  //Fill histogram
	  hdEta_dPhi_pT->Fill(dEta,dPhi,pT2);	  
	}
    }
}

/*
 * Write histograms to file
 */
void writeHistosToFile(char *outFileName="")
{
  TFile* fout = new TFile(outFileName, "RECREATE");

  hdEta_dPhi_pT->Write();
  hEta_pT->Write();
  hdN_dEta->Write();
  hNEvents->Write();
  hNpart->Write();
}

void ComputeCorrelations(char *filename = "ampt.dat", char* outFileName = "ampt_correlation2_out.root")
{ 
  //Initialize Variables
 
  //3D Histograms
  hdEta_dPhi_pT = new TH3F("hdEta_dPhi_pT","hdEta_dPhi_pT",200,-5.025,4.975,80,-41*TMath::Pi()/80,119*TMath::Pi()/80,30,0.2,3.2);

  //2D Histograms 
  hEta_pT = new TH2F("hEta_pT","hEta_pT",160,-4.025,3.975,15,0.2,3.2);
 
  //1D Histograms
  hdN_dEta = new TH1F("hdN_dEta","hdN_dEta",100,-2.222,2.178);
  hNEvents = new TH1F("hNEvents","hNEvents",1,0,2000);
  hNpart = new TH1F("hNpart","hNpart",200,0,50);

  //Open File
  ifstream myFile;
  myFile.open(filename);

  //Check to see if file was found
  if (!myFile) 
    {
      printf("Input file does not exist!\n");
      return;
    }
  else
    {
      cout << endl;
      cout << "--> Successfully opened file " << filename << endl;
      cout << endl;
    }

  //--------------------------------------
  // Run over events
  //--------------------------------------

  int evtnumber = 0;
  
  while (myFile) 
    {
      if(evtnumber % 100 == 0)
	{
	  cout << "-----> // // Reading Event No. " << evtnumber << endl;
	}
           
      //Read the event header
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

      myFile >> evtnumber >> testnum >> nlist >> impactpar >> npartproj >> nparttarg >> npartprojelas >> npartprojinelas >> nparttargelas >> nparttarginelas >> junk;

      //Avoid double reading the last entry
      if (!myFile) break;

      hNEvents->Fill(nlist);
      hNpart->Fill(npartprojinelas+nparttarginelas);

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

	  particles.push_back(part);  
	  
	} // End loop over particles

      //Compute 2-particle correlations
      computeCorrelation();
      particles.clear();
      evtnumber++;
    } 

  writeHistosToFile(outFileName);
}




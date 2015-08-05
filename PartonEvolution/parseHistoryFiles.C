//----------------------------------------------------
// Code to parse the partonic history files from AMPT
// and trace the scattering evolution of quarks prior
// to hadronization.
//
// Author: J. Orjuela-Koop
// 07-29-2015
//----------------------------------------------------
#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TRegexp.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TMath.h"
#include "TVector3.h"

#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>

using namespace std;

//-----------------------------------
// Parton Data Structure
//-----------------------------------
struct parton
{
  //Three-momentum
  float px;
  float py;
  float pz;

  //Four-position
  float x;
  float y;
  float z;
  float t;

  //Mass
  float m;

  //Angle by which it scatters in the CoM frame
  float angleCMS;
};

//-----------------------------------
// Variables
//-----------------------------------

//Categories for the number of scattering events. NSCAT = 0, 1, 2, >=3
const int NSCATT = 4;

//Vector of vectors to contain parton evolution for a given event
vector<vector<parton> > eventParticles;

//Vector to contain the participant plane angle for every event
vector<float> psi2;

//Output histograms
TProfile *hNscatt_pT;
TProfile *hv2_pT[NSCATT];
TH1F *hdN_dpT[NSCATT];
TH1F *hDeltaRapidity[NSCATT];
TH1F *hRapidity[NSCATT];
TH1F *hAngle[NSCATT];
TH2F *hEndStatev2[NSCATT];

//Counter for the number of events processed
int evtnumber = 0;

//-----------------------------------
// Functions
//-----------------------------------

void writeHistograms(char *outFileName="")
{
  TFile* fout = new TFile(outFileName, "RECREATE");
  for(int i=0; i<NSCATT; i++)
    {
      hv2_pT[i]->Write();
      hdN_dpT[i]->Write();
      hDeltaRapidity[i]->Write();
      hRapidity[i]->Write();
      hEndStatev2[i]->Write();
    }
  hNscatt_pT->Write();
}

/*
 * Load previously calculated participant plane angles for each event into a vector
 */
void loadParticipantPlanes()
{
  ifstream myFilePsi2;
  myFilePsi2.open("psi_2.txt");

  if (!myFilePsi2)
    {
      cout << "File does not exist!" << endl;
      return;
    }
  else
    {
      cout << "--> Successfully opened participant plane file" << endl;
    }

  float angle = 0;
  while(myFilePsi2)
    {
      myFilePsi2 >> angle;
      psi2.push_back(angle);
    }
}

/*
 * Draw spectra and v_2 for different number of scattering events
 */
void draw()
{
  gStyle->SetOptStat(0);

  TCanvas *cDeltaRapidity = new TCanvas("cDeltaRapidity","cDeltaRapidity",500,500);
  hDeltaRapidity[0]->Rebin(2);
  hDeltaRapidity[1]->Rebin(2);
  hDeltaRapidity[2]->Rebin(2);
  hDeltaRapidity[3]->Rebin(2);

  hDeltaRapidity[0]->SetLineColor(kViolet-3);
  hDeltaRapidity[1]->SetLineColor(kAzure-3);
  hDeltaRapidity[2]->SetLineColor(kSpring-6);
  hDeltaRapidity[3]->SetLineColor(kOrange-3);

  hDeltaRapidity[0]->Scale(1/(hDeltaRapidity[0]->GetMaximum()));
  hDeltaRapidity[1]->Scale(1/(hDeltaRapidity[1]->GetMaximum()));
  hDeltaRapidity[2]->Scale(1/(hDeltaRapidity[2]->GetMaximum()));
  hDeltaRapidity[3]->Scale(1/(hDeltaRapidity[3]->GetMaximum()));

  hDeltaRapidity[0]->SetTitle("");
  hDeltaRapidity[0]->Draw();
  hDeltaRapidity[1]->Draw("same");
  hDeltaRapidity[2]->Draw("same");
  hDeltaRapidity[3]->Draw("same");

  TCanvas *cRapidity = new TCanvas("cRapidity","cRapidity",500,500);
  hRapidity[0]->Rebin(2);
  hRapidity[1]->Rebin(2);
  hRapidity[2]->Rebin(2);
  hRapidity[3]->Rebin(2);

  hRapidity[0]->SetLineColor(kViolet-3);
  hRapidity[1]->SetLineColor(kAzure-3);
  hRapidity[2]->SetLineColor(kSpring-6);
  hRapidity[3]->SetLineColor(kOrange-3);

  hRapidity[0]->SetTitle("");
  hRapidity[0]->Draw();
  hRapidity[1]->Draw("same");
  hRapidity[2]->Draw("same");
  hRapidity[3]->Draw("same");

  TCanvas *cdNdpT = new TCanvas("cdNdpT","cdNdpT",500,500);
  hdN_dpT[0]->Rebin(2);
  hdN_dpT[1]->Rebin(2);
  hdN_dpT[2]->Rebin(2);
  hdN_dpT[3]->Rebin(2);

  hdN_dpT[0]->SetLineColor(kViolet-3);
  hdN_dpT[1]->SetLineColor(kAzure-3);
  hdN_dpT[2]->SetLineColor(kSpring-6);
  hdN_dpT[3]->SetLineColor(kOrange-3);

  hdN_dpT[0]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[1]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[2]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));
  hdN_dpT[3]->Scale(1.0/(hdN_dpT[0]->GetBinWidth(1)));

  hdN_dpT[0]->Scale(1.0/(float)evtnumber);
  hdN_dpT[1]->Scale(1.0/(float)evtnumber);
  hdN_dpT[2]->Scale(1.0/(float)evtnumber);
  hdN_dpT[3]->Scale(1.0/(float)evtnumber);

  hdN_dpT[0]->Draw();
  hdN_dpT[1]->Draw("same");
  hdN_dpT[2]->Draw("same");
  hdN_dpT[3]->Draw("same");

  TCanvas *cv2 = new TCanvas("cv2","cv2",500,500);
  hv2_pT[0]->SetLineColor(kViolet-3);
  hv2_pT[1]->SetLineColor(kAzure-3);
  hv2_pT[2]->SetLineColor(kSpring-6);
  hv2_pT[3]->SetLineColor(kOrange-3);

  hv2_pT[0]->Draw();
  hv2_pT[1]->Draw("same");
  hv2_pT[2]->Draw("same");
  hv2_pT[3]->Draw("same");

  TCanvas *cNscatt_pT = new TCanvas("cNscatt_pT","cNscatt_pT",500,500);
  hNscatt_pT->Rebin(2);
  hNscatt_pT->Draw();

  TCanvas *cScatteringAngle = new TCanvas("cScatteringAngle","cScatteringAngle",500,500);
  hAngle[0]->SetLineColor(kViolet-3);
  hAngle[1]->SetLineColor(kAzure-3);
  hAngle[2]->SetLineColor(kSpring-6);
  hAngle[3]->SetLineColor(kOrange-3);

  hAngle[1]->Scale(1.0/(hAngle[1]->GetMaximum()));
  hAngle[2]->Scale(1.0/(hAngle[2]->GetMaximum()));
  hAngle[3]->Scale(1.0/(hAngle[3]->GetMaximum()));

  hAngle[1]->Draw();
  hAngle[2]->Draw("same");
  hAngle[3]->Draw("same");
}

/*
 * Boost collision between to partons to center of mass frame
 * Partons 1 and 2 constitute the initial state
 * Partons 3 and 4 constitute the final state
 */
void boostCMS(parton p1, parton p2, parton p3, parton p4, float *angle_1_4, float *angle_2_3)
{
  //cout << "-> " << p1.px << ", " << p1.py << ", " << p1.pz << endl;
  //cout << "-> " << p2.px << ", " << p2.py << ", " << p2.pz << endl;
  //cout << "-> " << p3.px << ", " << p3.py << ", " << p3.pz << endl;
  //cout << "-> " << p4.px << ", " << p4.py << ", " << p4.pz << endl;

  float px1 = p1.px;
  float py1 = p1.py;
  float pz1 = p1.pz;
  float mass1 = p1.m;
  float energy1 =  sqrt (pow(px1,2)+pow(py1,2)+pow(pz1,2)+pow(mass1,2));
  TLorentzVector ev1(px1,py1,pz1,energy1);

  float px2 = p2.px;
  float py2 = p2.py;
  float pz2 = p2.pz;
  float mass2 = p2.m;
  float energy2 =  sqrt (pow(px2,2)+pow(py2,2)+pow(pz2,2)+pow(mass2,2));
  TLorentzVector ev2(px2,py2,pz2,energy2);

  float px3 = p3.px;
  float py3 = p3.py;
  float pz3 = p3.pz;
  float mass3 = p3.m;
  float energy3 =  sqrt (pow(px3,2)+pow(py3,2)+pow(pz3,2)+pow(mass3,2));
  TLorentzVector ev3(px3,py3,pz3,energy3);

  float px4 = p4.px;
  float py4 = p4.py;
  float pz4 = p4.pz;
  float mass4 = p4.m;
  float energy4 = sqrt (pow(px4,2)+pow(py4,2)+pow(pz4,2)+pow(mass4,2));
  TLorentzVector ev4(px4,py4,pz4,energy4);

  //Define boost vector to CoM
  TLorentzVector cms = ev1+ev2;
  ev1.Boost(-cms.BoostVector());
  ev2.Boost(-cms.BoostVector());
  ev3.Boost(-cms.BoostVector());
  ev4.Boost(-cms.BoostVector());
  cms.Boost(-cms.BoostVector());

  //Total momentum in CoM should be identically zero
  //cout << "--Total momentum in CoM = " << cms.P() << endl << endl;

  //cout << ev1.Px() << ", " << ev1.Py() << "; " << ev1.Pz() << endl;
  //cout << ev2.Px() << ", " << ev2.Py() << "; " << ev2.Pz() << endl;
  //cout << ev3.Px() << ", " << ev3.Py() << "; " << ev3.Pz() << endl;
  //cout << ev4.Px() << ", " << ev4.Py() << "; " << ev4.Pz() << endl;

  float angleLeft = ev1.Angle(ev4.Vect());
  float angleRight = ev2.Angle(ev3.Vect());

  //Wrap around to have angle between 0 and 2Pi
  if(angleLeft < 0)
    {
      angleLeft = angleLeft + 2*TMath::Pi();
    }
  else if(angleLeft > 2*TMath::Pi())
    {
      angleLeft = angleLeft - 2*TMath::Pi();
    }

  if(angleRight < 0)
    {
      angleRight = angleRight + 2*TMath::Pi();
    }
  else if(angleRight > 2*TMath::Pi())
    {
      angleRight = angleRight - 2*TMath::Pi();
    }

  //hAngle->Fill(angle);
  (*angle_1_4) = angleLeft;
  (*angle_2_3) = angleRight;
}

/*
 * Compute the rapidity of a parton given its mass and momentum
 */
float computeRapidity(parton p)
{
  float px = p.px;
  float py = p.py;
  float pz = p.pz;
  float mass = p.m;

  float energy =  sqrt (pow(px,2)+pow(py,2)+pow(pz,2)+pow(mass,2));
  TLorentzVector ev(px,py,pz,energy);

  return ev.Rapidity();
}

/*
 * Compute azimuthal angle of a given parton
 */
float computePhi(parton p)
{
  float px = p.px;
  float py = p.py;
  float pz = p.pz;
  float mass = p.m;

  float energy =  sqrt (pow(px,2)+pow(py,2)+pow(pz,2)+pow(mass,2));
  TLorentzVector ev(px,py,pz,energy);

  return ev.Phi();
}

/*
 * Loop over partons in each event, filling desired histograms
 */
void processEvent()
{
  //Get participant plane angle for the event at hand
  float psi2_angle = psi2[evtnumber-1];
  
  //Each entry in the 'eventParticles' vector corresponds to a single parton and its history
  for(int i=0; i<eventParticles.size(); i++)
    {   
      //Determine number of scattering events undergone by the parton
      vector<parton> v = eventParticles[i];
      int numstages = v.size();
      int numscatterings = v.size()-1;

      if(numscatterings > 2) numscatterings = 3;

      //Fill histograms
      float y = computeRapidity(v[numstages-1]);
      float delta_y = y - computeRapidity(v[0]);
      float pT = sqrt(pow(v[numstages-1].px,2) + pow(v[numstages-1].py,2));
      float phi = computePhi(v[numstages-1]);
      float v2 = TMath::Cos(2*(phi-psi2_angle));

      //Fill histograms for scattering angle in CoM for different number of scatterings
      for(int j=0; j<v.size(); j++)
	{
	  float scattAngle = v[j].angleCMS;
	  if(j<=2)
	    {
	      hAngle[j]->Fill(scattAngle);
	    }
	  else
	    {
	      hAngle[3]->Fill(scattAngle);
	    }
	}

      hRapidity[numscatterings]->Fill(y);
      hDeltaRapidity[numscatterings]->Fill(delta_y);
      hdN_dpT[numscatterings]->Fill(pT);
      hv2_pT[numscatterings]->Fill(pT,v2,1);
      hEndStatev2[numscatterings]->Fill(v2,pT);
      hNscatt_pT->Fill(pT,numscatterings,1);
    }
}

void parseHistoryFiles(char *initialInfoFile = "/direct/phenix+hhj/jdok/He3Au_Correlation/ampt_vn_pT/data/parton-initial-afterPropagation.dat", char *evolInfoFile = "/direct/phenix+hhj/jdok/He3Au_Correlation/ampt_vn_pT/data/parton-collisionsHistory.dat", char *outputFile = "evolution_out.root")
{
  //Load participant plane angles from file
  loadParticipantPlanes();

  //Initialize histograms
  for(int i=0; i<NSCATT; i++)
    {
      hDeltaRapidity[i] = new TH1F(Form("hDeltaRapidity_%i",i),Form("hDeltaRapidity_%i;#Delta y",i),500,-10,10);
      hRapidity[i] = new TH1F(Form("hRapidity_%i",i),Form("hRapidity_%i;y",i),500,-5,5);
      hdN_dpT[i] = new TH1F(Form("dN_dpT_%i",i),Form("dN_dpT_%i;p_{T};dN/dp_{T}",i),100,0,4);
      hv2_pT[i] = new TProfile(Form("hv2_pT_%i",i),Form("hv2_pT_%i;p_{T};v_{2}",i),50,0,5,-2,2);
      hAngle[i] = new TH1F(Form("hAngle_%i",i),Form("hAngle_%i;#theta",i),500,0,2*TMath::Pi());
      hEndStatev2[i] = new TH2F(Form("hEndStatev2_%i",i),Form("hEndStatev2_%i",i),1200,-1.2,1.2,50,0,5);
    }

  hNscatt_pT = new TProfile("hNscatt_pT","Profile of Nscatt vs pT",100,0,5,0,20);

  //Read initial parton information file
  ifstream myFileInitialInfo;
  myFileInitialInfo.open(initialInfoFile);

  //Read parton evolution from file
  ifstream myFileEvolInfo;
  myFileEvolInfo.open(evolInfoFile);

  if (!myFileInitialInfo)
    {
      cout << "File does not exist!" << endl;
      return;
    }
  else
    {
      cout << "--> Successfully opened file " << initialInfoFile << endl;
    }

  while (myFileInitialInfo)
    {
      if(evtnumber % 100 == 0)
        {
          cout << "-----> // // Reading Event No. " << evtnumber << endl;
        }
      
      //Read the event header                                                                                                                                                
      int nlist;
      int iterindex;
      int nbaryons;
      int nmesons;
      int nparticles;
      int nparticleszpc;
      
      myFileInitialInfo >> evtnumber >> iterindex >> nlist >> nbaryons >> nmesons >> nparticles >> nparticleszpc;
      
      //Avoid double reading the last entry                                                                                                                                  
      if (!myFileInitialInfo) break;
      
      //Loop over each parton in the event                                                                                                                                   
      for (int i=0;i<nlist;i++)
        {
          int partid;
          float pvec[3];
          float mass;
          double spacetime[4];
	  
          myFileInitialInfo >> partid >> pvec[0] >> pvec[1] >> pvec[2] >> mass >> spacetime[0] >> spacetime[1] >> spacetime[2] >> spacetime[3];

	  parton part;
	  part.px = pvec[0];
	  part.py = pvec[1];
	  part.pz = pvec[2];
	  part.x = spacetime[0];
	  part.y = spacetime[1];
	  part.z = spacetime[2];
	  part.t = spacetime[3];
	  part.m = mass;
	  part.angleCMS = -9999;

	  //Add the parton as the first stage of its own scattering evolution
	  vector<parton> aux;
	  aux.push_back(part);
	  eventParticles.push_back(aux);
	}
      
      //Get evolution information
      string line;
      TString tline;
      string description;
      int junk1;
      int junk2;
      int partonindex1;
      int partonindex2;
      int partid_init_1;
      float pvec_init_1[3];
      float mass_init_1;
      float spacetime_init_1[4];
      int partid_final_1;
      float pvec_final_1[3];
      float mass_final_1;
      float spacetime_final_1[4];
      int partid_init_2;
      float pvec_init_2[3];
      float mass_init_2;
      float spacetime_init_2[4];
      int partid_final_2;
      float pvec_final_2[3];
      float mass_final_2;
      float spacetime_final_2[4];
       
      string regexp = Form(" event,miss,iscat,jscat= %i [0-9]+ [0-9]+ [0-9]+", evtnumber);
      string regexp_break = Form("%i [0-9]+", evtnumber+1);
      while(myFileEvolInfo)
	{
	  TRegexp e(regexp);
	  TRegexp b(regexp_break);
	  std::getline(myFileEvolInfo,line);
	  tline = line;

	  if(tline.Contains(b))
	    {
	      break;
	    }

	  if(!tline.Contains(e))
	    {
	      continue;
	    }

	  stringstream headerstream(line);
	  headerstream >> description >> junk1 >> junk2 >> partonindex1 >> partonindex2;

	  myFileEvolInfo >> partid_init_1 >> pvec_init_1[0] >> pvec_init_1[1] >> pvec_init_1[2] >> mass_init_1 >> spacetime_init_1[0] >> spacetime_init_1[1] >> spacetime_init_1[2] >> spacetime_init_1[3];	
	  
	  myFileEvolInfo >> partid_init_2 >> pvec_init_2[0] >> pvec_init_2[1] >> pvec_init_2[2] >> mass_init_2 >> spacetime_init_2[0] >> spacetime_init_2[1] >> spacetime_init_2[2] >> spacetime_init_2[3];	 
	  
	  myFileEvolInfo >> partid_final_1 >> pvec_final_1[0] >> pvec_final_1[1] >> pvec_final_1[2] >> mass_final_1 >> spacetime_final_1[0] >> spacetime_final_1[1] >> spacetime_final_1[2] >> spacetime_final_1[3];
	  
	  myFileEvolInfo >> partid_final_2 >> pvec_final_2[0] >> pvec_final_2[1] >> pvec_final_2[2] >> mass_final_2 >> spacetime_final_2[0] >> spacetime_final_2[1] >> spacetime_final_2[2] >> spacetime_final_2[3];	 
	    
	  //Create scattered parton object from file information
	  parton part1;
	  part1.px = pvec_final_1[0];
	  part1.py = pvec_final_1[1];
	  part1.pz = pvec_final_1[2];
	  part1.x = spacetime_final_1[0];
	  part1.y = spacetime_final_1[1];
	  part1.z = spacetime_final_1[2];
	  part1.t = spacetime_final_1[3];
	  part1.m = mass_final_1;

	  parton part2;
	  part2.px = pvec_final_2[0];
	  part2.py = pvec_final_2[1];
	  part2.pz = pvec_final_2[2];
	  part2.x = spacetime_final_2[0];
	  part2.y = spacetime_final_2[1];
	  part2.z = spacetime_final_2[2];
	  part2.t = spacetime_final_2[3];
	  part2.m = mass_final_2;
	  
	  parton part3;
	  part3.px = pvec_init_1[0];
	  part3.py = pvec_init_1[1];
	  part3.pz = pvec_init_1[2];
	  part3.x = spacetime_init_1[0];
	  part3.y = spacetime_init_1[1];
	  part3.z = spacetime_init_1[2];
	  part3.t = spacetime_init_1[3];
	  part3.m = mass_init_1;

	  parton part4;
	  part4.px = pvec_init_2[0];
	  part4.py = pvec_init_2[1];
	  part4.pz = pvec_init_2[2];
	  part4.x = spacetime_init_2[0];
	  part4.y = spacetime_init_2[1];
	  part4.z = spacetime_init_2[2];
	  part4.t = spacetime_init_2[3];
	  part4.m = mass_init_2;
 
	  //Find the angular distribution of scattered particles in the CoM frame
	  float angle_1_4 = -9999.0;
	  float angle_2_3 = -9999.0;
	  boostCMS(part3, part4, part1, part2, &angle_1_4, &angle_2_3);	   

	  part1.angleCMS = angle_2_3;
	  part2.angleCMS = angle_1_4;

	  eventParticles[partonindex1-1].push_back(part1);
	  eventParticles[partonindex2-1].push_back(part2);
	}

      //Run analysis on particles on an event-by-event basis
      processEvent();
      eventParticles.clear();
    }
  writeHistograms(outputFile);
  draw();
}

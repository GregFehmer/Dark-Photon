#include <TH1D.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TChain.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <TGraphAsymmErrors.h>
#include <TVector3.h>
#include <TRandom.h>
#include <TMath.h>

using std::string;

double MUON_MASS = 105.6583745e-03; //GeV
double Electron_Mass = 0.51099895e-03; // GeV
int nM = 0;
int nE = 0;
int NM = 0;
int NE = 0;

bool sameValHighPrecision(double a, double b)
{
   return fabs(a - b) < 1.000e-08;
}

TLorentzVector fillTLorentzVector(double pT, double eta, double phi, double mass)
{
  TLorentzVector object_p4;
  object_p4.SetPtEtaPhiM(pT, eta, phi, mass);
  return object_p4;
}

typedef struct
{
  float pt; 
  float phi;
  float eta;
  TLorentzVector lep_lv;
  Bool_t id;
  int charge;
  float deltaR;
  
} leptonInfo;

int Histogram(std::string infile, std::string treeStr)
{
  std::string inputfilename=(infile+".root").c_str();
  TFile *inputFile = new TFile((inputfilename).c_str());
  TChain *tree=new TChain(treeStr.c_str());
  tree->Add(inputfilename.c_str());

  Float_t         Muon_pt[15]; 
  Float_t         Muon_phi[15];
  Float_t         Muon_eta[15];  
  UInt_t          nMuon; 
  Bool_t          Muon_tightId[15];
  Int_t           Muon_charge[15];
  Float_t 	Electron_pt[15];
  Float_t 	Electron_phi[15];
  Float_t 	Electron_eta[15];
  UInt_t 	nElectron;
  Bool_t 	Electron_mvaFall17V2noIso_WPL[15];
  Int_t 	Electron_charge[15];
  //mvaFall17V2noIso
  
  
  tree->SetBranchAddress("Muon_pt", &Muon_pt); 
  tree->SetBranchAddress("Muon_phi", &Muon_phi);
  tree->SetBranchAddress("Muon_eta", &Muon_eta);
  tree->SetBranchAddress("Muon_tightId", &Muon_tightId);
  tree->SetBranchAddress("nMuon", &nMuon);
  tree->SetBranchAddress("Muon_charge", &Muon_charge);

  tree->SetBranchAddress("Electron_pt", &Electron_pt);
  tree->SetBranchAddress("Electron_phi", &Electron_phi);
  tree->SetBranchAddress("Electron_eta", &Electron_eta);
  tree->SetBranchAddress("Electron_mvaFall17V2noIso_WPL", &Electron_mvaFall17V2noIso_WPL);
  tree->SetBranchAddress("nElectron", &nElectron);
  tree->SetBranchAddress("Electron_charge", &Electron_charge);
  


  unsigned int nEvents=tree->GetEntries();
  std::cout << "Reading events = " << nEvents << std::endl;
  TH1D *h_m1 = new TH1D("h_m1", "Dimuon mass; M_{dimuon} [GeV]; Events/GeV", 100, 0, 5);h_m1->Sumw2();
  TH1D *h_m2 = new TH1D("h_m2", "Dielectron mass; M_{dielectron} [GeV]; Events/GeV", 200, 0, 20);h_m2->Sumw2();
  TH1D *h_nMuons = new TH1D("h_nMuons", "Number of Muons per Event; nMuons; nEvents", 10, 0, 10);
  TH1D *h_nElectrons = new TH1D("h_nElectrons", "Number of Electrons per Event; nElectrons; nEvents", 10, 0, 10);
  TH1D *h_deltaRMuons = new TH1D("h_deltaRMuons", "DeltaR Muons; DeltaR; nEvents", 100, 0, 6);
  TH1D *h_deltaRElectrons = new TH1D("h_deltaRElectrons", "DeltaR Electrons; DeltaR; nEvents", 100, 0, 6);
  
  TH1D *h_PtMuons = new TH1D("h_PtMuons", "Pt Electrons; Pt; nEvents", 1000, 0, 1000);
  TH1D *h_PtElectrons = new TH1D("h_PtElectrons", "Pt Electrons; P; nEvents", 1000, 0, 1000);
  TH1D *h_PhiMuons = new TH1D("h_PhiMuons", "Phi Electrons; Phi; nEvents", 100, -3.5, 3.5);
  TH1D *h_PhiElectrons = new TH1D("h_PhiElectrons", "Phi Electrons; Phi; nEvents", 100, -3.5, 3.5);
  TH1D *h_EtaMuons = new TH1D("h_EtaMuons", "Eta Electrons; Eta; nEvents", 100, -3, 3);
  TH1D *h_EtaElectrons = new TH1D("h_EtaElectrons", "Eta Electrons; Eta; nEvents", 100, -3, 3);



  for (unsigned int i=0; i<nEvents; ++i)
  {

    tree->GetEvent(i); 
    std::vector<leptonInfo> v_muons;
    for (unsigned int i=0; i<nMuon; i++)
    {
      leptonInfo lepton;
      lepton.pt = Muon_pt[i];
      lepton.phi = Muon_phi[i];
      lepton.eta = Muon_eta[i];
      lepton.id = Muon_tightId[i];
      lepton.charge = Muon_charge[i];
      lepton.lep_lv.SetPtEtaPhiM(Muon_pt[i], Muon_eta[i], Muon_phi[i], MUON_MASS);

      if(lepton.id==1) 
      {
        v_muons.push_back(lepton);
      }
    }

    std::vector<leptonInfo> v_Electrons;
    for (unsigned int i = 0; i < nElectron; i++)
      {
       leptonInfo lepton;
       lepton.pt = Electron_pt[i];
       lepton.phi = Electron_phi[i];
       lepton.eta = Electron_eta[i];
       lepton.id = Electron_mvaFall17V2noIso_WPL[i];
       lepton.charge = Electron_charge[i];
       lepton.lep_lv.SetPtEtaPhiM(Electron_pt[i], Electron_eta[i], Electron_phi[i], Electron_Mass);
       

      if (lepton.id == 1)
      {
        v_Electrons.push_back(lepton);
      }
    }
      
   if(v_muons.size() >= 2)
    {
      NM++;
      TLorentzVector l1 = v_muons.at(0).lep_lv;
      TLorentzVector l2 = v_muons.at(1).lep_lv;
      Float_t deltaR_M = l1.DeltaR(l2);

      h_m1->Fill((l1+l2).M());
      h_deltaRMuons->Fill(deltaR_M);
      
      h_PtMuons->Fill(v_muons.at(0).pt);
      h_PhiMuons->Fill(v_muons.at(0).phi);
      h_EtaMuons->Fill(v_muons.at(0).eta);
      
      h_PtMuons->Fill(v_muons.at(1).pt);
      h_PhiMuons->Fill(v_muons.at(1).phi);
      h_EtaMuons->Fill(v_muons.at(1).eta);
      
    }
      
   if(v_Electrons.size() >= 2)
    {
      NE++;
      TLorentzVector l3 = v_Electrons.at(0).lep_lv;
      TLorentzVector l4 = v_Electrons.at(1).lep_lv;
      Float_t deltaR_E = l3.DeltaR(l4);

      h_m2->Fill((l3+l4).M());
      h_deltaRElectrons->Fill(deltaR_E);
      
      h_PtElectrons->Fill(v_Electrons.at(0).pt);
      h_PhiElectrons->Fill(v_Electrons.at(0).phi);
      h_EtaElectrons->Fill(v_Electrons.at(0).eta);
      
      h_PtElectrons->Fill(v_Electrons.at(1).pt);
      h_PhiElectrons->Fill(v_Electrons.at(1).phi);
      h_EtaElectrons->Fill(v_Electrons.at(1).eta);
         
    }
    
    h_nMuons->Fill(v_muons.size());
    h_nElectrons->Fill(v_Electrons.size());
    
   
    }

  std::string histfilename=("histograms_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_m1->Write();
  h_m2->Write();
  h_nMuons->Write(); 
  h_nElectrons->Write();
  h_deltaRMuons->Write();
  h_deltaRElectrons->Write();
  
  h_PtMuons->Write();
  h_PtElectrons->Write();
  h_PhiMuons->Write();
  h_PhiElectrons->Write();
  h_EtaMuons->Write();
  h_EtaElectrons->Write();
  
  tFile->Close();
  inputFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
} 

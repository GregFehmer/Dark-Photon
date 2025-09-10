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
//double Muon-Mass2023=105.6583755*10e-03;
double Electron_Mass = 0.51099895e-03; // GeV
int nM = 0;
int nE = 0;
int NM = 0;
int NE = 0;
int ngen = 0;
int npho = 0;

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
  
  float genmass;
  float genpt;
  float genphi;
  float geneta;
  int genID;
  int gensize;
  
  float phomass;
  float phopt;
  float phophi;
  float phoeta;
  //int genID;
  //int gensize;
} leptonInfo;

int ReadDPFile(std::string infile, std::string treeStr)
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
  
  Float_t GenDressedLepton_mass[15];
  Float_t GenDressedLepton_pt[15];
  Float_t GenDressedLepton_phi[15];
  Float_t GenDressedLepton_eta[15];
  Int_t GenDressedLepton_pdgId[15];
  UInt_t nGenDressedLepton;
  
  Float_t Photon_pt[15];
  Float_t Photon_phi[15];
  Float_t Photon_eta[15];
  
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
  
  tree->SetBranchAddress("GenDressedLepton_mass", &GenDressedLepton_mass);
  tree->SetBranchAddress("GenDressedLepton_pt", &GenDressedLepton_pt);
  tree->SetBranchAddress("GenDressedLepton_phi", &GenDressedLepton_phi);
  tree->SetBranchAddress("GenDressedLepton_eta", &GenDressedLepton_eta);
  tree->SetBranchAddress("GenDressedLepton_pdgId", &GenDressedLepton_pdgId);
  tree->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton);
  
  tree->SetBranchAddress("Photon_pt", &Photon_pt);
  tree->SetBranchAddress("Photon_phi", &Photon_phi);
  tree->SetBranchAddress("Photon_eta", &Photon_eta);

  unsigned int nEvents=tree->GetEntries();
  std::cout << "Reading events = " << nEvents << std::endl;
  TH1D *h_m1 = new TH1D("h_m1", "Dimuon mass; M_{dimuon} [GeV]; Events/GeV", 100, 0, 5);h_m1->Sumw2();
  TH1D *h_m2 = new TH1D("h_m2", "Dielectron mass; M_{dielectron} [GeV]; Events/GeV", 200, 0, 20);h_m2->Sumw2();
  TH1D *h_m3 = new TH1D("h_m3", "Dilepton mass; M_{e+Muon+} [GeV]; Events/GeV", 200, 0, 20);h_m3->Sumw2();
  TH1D *h_nMuons = new TH1D("h_nMuons", "Number of Muons per Event; nMuons; nEvents", 10, 0, 10);
  TH1D *h_nElectrons = new TH1D("h_nElectrons", "Number of Electrons per Event; nElectrons; nEvents", 10, 0, 10);
  TH1D *h_deltaRMuons = new TH1D("h_deltaRMuons", "DeltaR Muons; DeltaR; nEvents", 100, 0, 6);
  TH1D *h_deltaRElectrons = new TH1D("h_deltaRElectrons", "DeltaR Electrons; DeltaR; nEvents", 100, 0, 6);
  TH1D *h_deltaPhi = new TH1D("h_deltaPhiElectrons", "DeltaPhi^2 Electrons; DeltaPhi; nEvents", 100, 0, 6);
  TH1D *h_deltaEta = new TH1D("h_deltaEtaElectrons", "DeltaEta^2 Electrons; DeltaEta; nEvents", 100, 0, 6);
  TH1D *h_Phi = new TH1D("h_PhiElectrons", "Phi Electrons; Phi; nEvents", 100, -3.2, 3.2);
  
  TH1D *h_GenMass = new TH1D("h_GenMass", "Gen Mass", 100, 0, 0.2);
  TH1D *h_GenPhi = new TH1D("h_GenPhi", "Gen Phi", 100, 3.5, 3.5);
  TH1D *h_GenEta = new TH1D("h_GenEta", "Gen Eta", 100, 3.5, 3.5);


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
       
       lepton.phopt = Photon_pt[i];
       lepton.phophi = Photon_phi[i];
       lepton.phoeta = Photon_eta[i];
       
       
      if (lepton.id == 1)
      {
        v_Electrons.push_back(lepton);
      }
    }
      
    std::vector<leptonInfo> v_genElectrons;
    for (unsigned int i = 0; i < nGenDressedLepton; i++)
      {
       leptonInfo lepton;
       lepton.genmass = GenDressedLepton_mass[i];
       lepton.genpt = GenDressedLepton_pt[i];
       lepton.genphi = GenDressedLepton_phi[i];
       lepton.geneta = GenDressedLepton_eta[i];
       lepton.genID = GenDressedLepton_pdgId[i];
       
       v_genElectrons.push_back(lepton);
    }
    
    if(v_muons.size() >= 2 and v_muons.at(0).charge*v_muons.at(1).charge == -1 and v_muons.at(0).pt > 5 and v_muons.at(1).pt > 5)
    {
      NM++;
      TLorentzVector l1 = v_muons.at(0).lep_lv;
      TLorentzVector l2 = v_muons.at(1).lep_lv;
      Float_t deltaR_M = 5;
      
      for (unsigned int i=1; i<v_muons.size() ; ++i)
      {
        TLorentzVector ln = v_muons.at(i).lep_lv;
        Float_t deltaR_n = l1.DeltaR(ln);
        
        if(deltaR_n < deltaR_M) 
        {
          l2 = ln;
          deltaR_M = l1.DeltaR(l2);
        }
      }

      h_m1->Fill((l1+l2).M());
      h_deltaRMuons->Fill(deltaR_M);
      
     if (deltaR_M <= 0.1)
     {
       nM++;
     }
      
    }
      
   if(v_Electrons.size() >= 2 and v_Electrons.at(0).charge*v_Electrons.at(1).charge == -1 and v_Electrons.at(0).pt > 5 and v_Electrons.at(1).pt > 5)
    {
      NE++;
      TLorentzVector l3 = v_Electrons.at(0).lep_lv;
      TLorentzVector l4 = v_Electrons.at(1).lep_lv;
      Float_t deltaR_E = 5;
      
      for (unsigned int i=1; i<v_Electrons.size() ; ++i)
      {
        TLorentzVector ln = v_Electrons.at(i).lep_lv;
        Float_t deltaR_n = l3.DeltaR(ln);
        
        if(deltaR_n < deltaR_E) 
        {
          l4 = ln;
          deltaR_E = l3.DeltaR(l4);
        }
      }
      
      h_m2->Fill((l3+l4).M());
      h_deltaRElectrons->Fill(deltaR_E);
      
      //h_deltaPhi->Fill(dPhi_E);
      //h_deltaEta->Fill(dEta_E);
      
      h_Phi->Fill(v_Electrons.at(0).phi);
      h_Phi->Fill(v_Electrons.at(1).phi);
      
     if (deltaR_E <= 0.1)
     {
       nE++;
     }
     
     npho++;
     
 /*    if (npho>=50 and npho<=60)
     {
       int size = v_Electrons.size();
       int size2 = v_genElectrons.size();
       
       std::cout<<" "<<std::endl;
       std::cout<<"Electron size = "<<size<<std::endl;
       
       for (unsigned int i=0; i<size; ++i)
       {
         std::cout<<"Electron at index " <<i<<  " has pt = " <<v_Electrons.at(i).pt<<", phi = " <<v_Electrons.at(i).phi<<", eta = "<<v_Electrons.at(i).eta<<", and charge = "<<v_Electrons.at(i).charge<<std::endl;
         
        }
        
       for (unsigned int i=0; i<size; ++i)
       {
         std::cout<<"Photon at index " <<i<< " has pt = " <<v_Electrons.at(i).phopt<<", phi = " <<v_Electrons.at(i).phophi<<", and eta = "<<v_Electrons.at(i).phoeta<<std::endl;
         
        }
        
       for (unsigned int i=0; i<size2; ++i)
       {
         std::cout<<"Gen Electron at index " <<i<<  " has pt = " <<v_genElectrons.at(i).genpt<<", phi = " <<v_genElectrons.at(i).genphi<<", eta = "<<v_genElectrons.at(i).geneta<<", and ID = "<<v_genElectrons.at(i).genID<<std::endl;
         
        }
      } */
    } 
    
    h_nMuons->Fill(v_muons.size());
    h_nElectrons->Fill(v_Electrons.size());
   
  //if(v_genElectrons.size() >= 4 and v_genElectrons.at(0).genpt > 5 and v_genElectrons.at(1).genpt > 5 )
  /*
   if(v_Electrons.size() >= 2 and v_Electrons.at(0).charge*v_Electrons.at(1).charge == -1 and v_Electrons.at(0).pt > 5 and v_Electrons.at(1).pt > 5)
  
     {
     ngen++;
     
     if (ngen>=100 and ngen<=110)
     {
       int size = v_genElectrons.size();
     
       std::cout<<"Electron size = "<<size<<std::endl;
       
       for (unsigned int i=0; i<size; ++i)
       {
         std::cout<<"GenDressedLepton at index " <<i<<  " has ID = "<<v_genElectrons.at(i).genID<<", phi = " <<v_genElectrons.at(i).genphi<<", and eta = "<<v_genElectrons.at(i).geneta<<std::endl;
         std::cout<<" "<<std::endl;
       }
     }
     
    
     h_GenMass->Fill(v_Electrons.at(0).genmass);
     h_GenPhi->Fill(v_Electrons.at(0).genphi);
     h_GenEta->Fill(v_Electrons.at(0).geneta);
     
     h_GenMass->Fill(v_Electrons.at(1).genmass);
     h_GenPhi->Fill(v_Electrons.at(1).genphi);
     h_GenEta->Fill(v_Electrons.at(1).geneta);
     
     
    }   
*/
    
   if(v_muons.size()==1 && v_Electrons.size() == 1 and v_muons.at(0).charge*v_Electrons.at(0).charge==-1)
    {
      TLorentzVector l5 = v_muons.at(0).lep_lv;
      TLorentzVector l6 = v_Electrons.at(0).lep_lv;
      //std::cout << "(l1+l2).M() = " << (l1+l2).M() << std::endl;
      h_m3->Fill((l5+l6).M());
    }

  }
  
  std::cout << "Total muons counted = " << NM << std::endl;
  std::cout << "Total electrons counted = " << NE << std::endl;
  std::cout << "Muons that pass efficiency = " << nM << std::endl;
  std::cout << "Electrons that pass efficiency = " << nE << std::endl;

  std::string histfilename=("output_"+infile+".root").c_str();
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");
  h_m1->Write();
  h_m2->Write();
  h_m3->Write();
  h_nMuons->Write(); 
  h_nElectrons->Write();
  h_deltaRMuons->Write();
  h_deltaRElectrons->Write();
  h_deltaPhi->Write();
  h_deltaEta->Write();
  h_deltaPhi->Write();
  
  h_GenMass->Write();
  h_GenPhi->Write();
  h_GenEta->Write();
  tFile->Close();
  inputFile->Close();
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  return 0;
} 

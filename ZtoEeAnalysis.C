#define ZtoEeAnalysis_cxx
#include "ZtoEeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TF1.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TParticleClassPDG.h>

////////////////////////////
// Breit-Wigner Functions //
////////////////////////////

//                2              Γ^2 * M^2
//B(m;M,Γ) = N * ---  -------------------------------       ----> Relativistic Breit-Wigner Distribution
//                π    (m^2 - M^2)^2 + m^4(Γ^2 / M^2)

//Γ = Gamma
//π = Pi
//x[0] = m Gamma=par[1]  M=par[2]
Double_t mybwMass(Double_t *x, Double_t *par)
{
  Double_t arg1 = 14.0 / 22.0;                                                               // 2 over pi
  Double_t arg2 = par[1] * par[1] * par[2] * par[2];                                         // Γ^2 * M^2
  Double_t arg3 = ((x[0] * x[0]) - (par[2] * par[2])) * ((x[0] * x[0]) - (par[2] * par[2])); // (m^2 - M^2)^2
  Double_t arg4 = x[0] * x[0] * x[0] * x[0] * ((par[1] * par[1]) / (par[2] * par[2]));       // m^4(Γ^2 / M^2)
  return par[0] * arg1 * arg2 / (arg3 + arg4);                                               // return
}
void ZtoEeAnalysis::Loop() // t is getter
{
  //   In a ROOT session, you can do:
  //      root> .L ZtoEeAnalysis.C
  //      root> ZtoEeAnalysis t
  //      root> t.GetEntry(12); // Fill t data members with entry number 12
  //      root> t.Show();       // Show values of entry 12
  //      root> t.Show(16);     // Read and show values of entry 16
  //      root> t.Loop();       // Loop on all entries
//     This is the loop skeleton where:
  //    jentry is the global entry number in the chain
  //    ientry is the entry number in the current Tree
  //  Note that the argument to GetEntry must be:
  //    jentry for TChain::GetEntry
  //    ientry for TTree::GetEntry and TBranch::GetEntry
  //
  //       To read only selected branches, Insert statements like:
  // METHOD1:
  //    fChain->SetBranchStatus("*",0);  // disable all branches
  //    fChain->SetBranchStatus("branchname",1);  // activate branchname
  // METHOD2: replace line
  //    fChain->GetEntry(jentry);       //read all branches
  //by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0)
    return;

  // configs
  int rapMin = -20;
  int rapMax =  20;

  /// TPave Settings for Statistics and Fitting Box
  gStyle->SetOptStat("KSiouRMen");
  gStyle->SetOptFit(11112);
  //output file
  TFile *output = new TFile("ZtoEeAnalysisoutput.root", "RECREATE");

  //define histos
  TH1F *h_eventnumber = new TH1F("h_eventnumber", "NumberOfEvents", 5, 0, 5);
 
  // For Z Boson
  TH1F *histMassPairZ = new TH1F("histMassPairZ", ";m_{ll}(GeV/c^{2});Entries", 5. / 0.02, 0., 150.);
  histMassPairZ->SetTitle("MassPairZ Histogram");
  TH1F *histPtZ = new TH1F("histPtZ", ";p_{T}(GeV/c);Entries", 100, 0., 15.);
  histPtZ->SetTitle("PtZ Histogram");
  TH1F *histRapidityZ = new TH1F("histRapidityZ", ";y;Entries", 100, rapMin, rapMax);
  histRapidityZ->SetTitle("Rapidity Z Histogram");
  TH1F *histEtaZ = new TH1F("histEtaZ", ";#eta;Entries", 100, rapMin, rapMax);
  histEtaZ->SetTitle("EtaZ Histogram");
  TH1F *histPhiZ = new TH1F("histPhiZ", ";#phi;Entries", 100, rapMin, rapMax);
  histPhiZ->SetTitle("PhiZ Histogram");

  // For Electron
  TH1F *histMassPairElectron = new TH1F("histMassPairElectron", ";m_{ll}(GeV/c^{2});Entries", 100, -1, 1);
  histMassPairElectron->SetTitle("MassPairElectron Histogram");
  TH1F *histPtElectron = new TH1F("histPtElectron", ";p_{T}(GeV/c);Entries", 100, 0., 100.);
  histPtElectron->SetTitle("PtElectron Histogram");
  TH1F *histRapidityElectron = new TH1F("histRapidityElectron", ";y;Entries", 100, rapMin, rapMax);
  histRapidityElectron->SetTitle("RapidityElectron Histogram");
  TH1F *histEtaElectron = new TH1F("histEtaElectron", ";#eta;Entries", 100, rapMin, rapMax);
  histEtaElectron->SetTitle("EtaElectron Histogram");
  TH1F *histPhiElectron = new TH1F("histPhiElectron", ";#phi;Entries", 100, rapMin, rapMax);
  histPhiElectron->SetTitle("PhiElectron Histogram");

  // For AntiElectron
  TH1F *histMassPairAntiElectron = new TH1F("histMassPairAntiElectron", ";m_{ll}(GeV/c^{2});Entries", 5. / 0.02, -1., 1.);
  histMassPairAntiElectron->SetTitle("MassPairElectron Histogram");
  TH1F *histPtAntiElectron = new TH1F("histPtAntiElectron", ";p_{T}(GeV/c);Entries", 100, 0., 100.);
  histPtAntiElectron->SetTitle("PtAntiElectron Histogram");
  TH1F *histRapidityAntiElectron = new TH1F("histRapidityAntiElectron", ";y;Entries", 100, rapMin, rapMax);
  histRapidityAntiElectron->SetTitle("RapidityAntiElectron Histogram");
  TH1F *histEtaAntiElectron = new TH1F("histEtaAntiElectron", ";#eta;Entries", 100, rapMin, rapMax);
  histEtaAntiElectron->SetTitle("EtaAntiElectron Histogram");
  TH1F *histPhiAntiElectron = new TH1F("histPhiAntiElectron", ";#phi;Entries", 100, rapMin, rapMax);
  histPhiAntiElectron->SetTitle("PhiAntiElectron Histogram");

  //get event number
  int eventnumber = 0;
Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  { // start main event loop
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    eventnumber++;
    //cout << "event number = " <<  eventnumber << endl;
    h_eventnumber->Fill(1);
    // Define Variables as 4 Vectors
    TLorentzVector m1;
    TLorentzVector m2;
    TLorentzVector kZMass;

    // Declarations
    m1.SetPtEtaPhiM(Particle_PT[2], Particle_Eta[2], Particle_Phi[2], Particle_M[2]);
    m2.SetPtEtaPhiM(Particle_PT[3], Particle_Eta[3], Particle_Phi[3], Particle_M[3]);

    // Mass Calculation for Z Boson
    kZMass = m1 + m2;

    ////////////////////////////
    // Filling the Histograms //
    ////////////////////////////

    // Electron
    histMassPairElectron->Fill(m1.M());
    histPtElectron->Fill(m1.Pt());
    histRapidityElectron->Fill(m1.Rapidity());
    histEtaElectron->Fill(m1.Eta());
    histPhiElectron->Fill(m1.Phi());

    // AntiElectron
    histMassPairAntiElectron->Fill(m2.M());
    histPtAntiElectron->Fill(m2.Pt());
    histRapidityElectron->Fill(m1.Rapidity());
    histEtaAntiElectron->Fill(m2.Eta());
    histPhiAntiElectron->Fill(m2.Phi());

// Z
    histMassPairZ->Fill(kZMass.M());
    histPtZ->Fill(kZMass.Pt());
    histRapidityZ->Fill(kZMass.Rapidity());
    histEtaZ->Fill(kZMass.Eta());
    histPhiZ->Fill(kZMass.Phi());

  } // end main event loop

  ////////////////////////////////////////////////
  // Declare some bin parameters for BW Fitting //
  ////////////////////////////////////////////////

  //for Z Boson
  int divisionMass = histMassPairZ->GetNbinsX();
  float massMIN = histMassPairZ->GetBinLowEdge(1);
  float massMAX = histMassPairZ->GetBinLowEdge(divisionMass + 1);
  //float BIN_SIZE_MASS = histMassPairZ->GetBinWidth(1);

  TF1 *bwFuncMass = new TF1("mybwMass", mybwMass, massMIN, massMAX, 3);
  bwFuncMass->SetParameter(0, 1.0);
  bwFuncMass->SetParName(0, "const");
  bwFuncMass->SetParameter(1, 91);
  bwFuncMass->SetParName(1, "sigma");
  bwFuncMass->SetParameter(2, 2.365);
  bwFuncMass->SetParName(2, "mean");
 histMassPairZ->Fit("mybwMass", "QR");
  //TF1 *fithistMassPairZ = histMassPairZ->GetFunction("mybwMass");

  output->cd();

  h_eventnumber->Write();
  histMassPairElectron->Write();
  histPtElectron->Write();
  histRapidityElectron->Write();
  histEtaElectron->Write();
  histPhiElectron->Write();

  histMassPairAntiElectron->Write();
  histPtAntiElectron->Write();
  histRapidityAntiElectron->Write();
  histEtaAntiElectron->Write();
  histPhiAntiElectron->Write();

  histMassPairZ->Write();
  histPtZ->Write();
  histRapidityZ->Write();
  histEtaZ->Write();
  histPhiZ->Write();

  output->Close();
  cout << "Finished!!!" << endl;
}
 



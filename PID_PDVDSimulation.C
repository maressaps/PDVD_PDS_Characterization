/*
 * =============================================================================
 * Photon Distribution and Muon Selection Analysis
 * =============================================================================
 *
 * Purpose: compare the photon distributions on the two ProtoDUNE-VD
 * cathode sides for simulated electrons, muons, and pions.
 *
 * For each particle sample, it:
 *   - reads photoncounter/PhotonsTree;
 *   - sums photons associated with cathode 1 (even channels) and cathode 2 (odd channels);
 *   - creates 1D photon distributions for each cathode;
 *   - creates a 2D cathode-1 vs cathode-2 photon scatter plot;
 *   - defines an elliptical muon-selection region;
 *   - scans ellipse parameters to maximize efficiency x purity;
 *   - calculates the final efficiency, purity, and confusion matrix;
 *   - stores diagnostic optimization heatmaps.
 *
 * Cathode channel grouping
 * ------------------------
 *   Cathode 1: channels 4, 6, 8, 10
 *   Cathode 2: channels 5, 7, 9, 11
 *
 * Event grouping
 * --------------
 * The macro assumes 20 consecutive PhotonsTree entries per physical event.
 * After every 20 entries, the accumulated cathode photon totals are used to
 * fill the event-level observables and are then reset.
 *
 * Muon selection
 * --------------
 * Muons are treated as signal and electrons+pions as background.
 *
 *   efficiency = TP / (TP + FN)
 *   purity     = TP / (TP + FP)
 *   score      = efficiency * purity
 *
 * The selection region is a rotated ellipse in the plane
 *   x = photons on cathode 1
 *   y = photons on cathode 2.
 *
 * Inputs
 * ------
 *   inputfileELectrons : electron ROOT sample
 *   inputfileMuons     : muon ROOT sample
 *   inputfilePions     : pion ROOT sample
 *   energy             : beam energy [GeV]
 *
 * Output
 * ------
 *   PhotonsPerCathode_<energy>GeV.root
 *
 * Contact
 * -------
  * Written by Maressa Sampaio and Laura Paulucci
 * maressap@ifi.unicamp.br
 * =============================================================================
 */
 #include <cstdio>
#include <iostream>
#include <string.h>
#include <RtypesCore.h>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>
#include <TDirectoryFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TClonesArray.h>
#include <THn.h>
#include <THnSparse.h>
#include <TH1D.h>
#include <TROOT.h>
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

struct Point {
    double x;
    double y;
};

bool inEllipse(double x, double y,
               double mux, double muy,
               double a, double b,
               double theta) {

  double dx = x - mux;
  double dy = y - muy;

  double xr =  cos(theta)*dx + sin(theta)*dy;
  double yr = -sin(theta)*dx + cos(theta)*dy;

  double val = (xr*xr)/(a*a) + (yr*yr)/(b*b);

  return val < 1.0;
}



void compute_metrics(double mux, double muy,
                     double a, double b, double theta,
                     double &eff, double &pur,
                     const vector<Point>& mu,
                     const vector<Point>& e,
                     const vector<Point>& pi) {

        int TP = 0;
        int FP = 0;

        for (auto &p : mu) {
          if (inEllipse(p.x, p.y, mux, muy, a, b, theta))
            TP++;
        }

        for (auto &p : e) {
          if (inEllipse(p.x, p.y, mux, muy, a, b, theta))
            FP++;
        }

        for (auto &p : pi) {
          if (inEllipse(p.x, p.y, mux, muy, a, b, theta))
            FP++;
        }

        eff = (double)TP / mu.size();
        pur = TP + FP > 0 ? (double)TP / (TP + FP) : 0;
}


vector<double> linspace(double start, double stop, int N) {
    vector<double> values;

    if (N == 1) {
        values.push_back(start);
        return values;
    }

    double step = (stop - start) / (N - 1);

    for (int i = 0; i < N; i++) {
        values.push_back(start + i * step);
    }

    return values;
}

struct EllipseObjective {
    vector<Point> *X_mu;
    vector<Point> *X_e;
    vector<Point> *X_pi;
    double theta;

    double operator()(const double *par) {
        double mux = par[0];
        double muy = par[1];
        double a   = par[2];
        double b   = par[3];


        double eff, pur;
        compute_metrics(mux, muy, a, b, theta,
                        eff, pur, *X_mu, *X_e, *X_pi);

        return -(eff * pur);
    }
};



void PerCathodeplots(const TString inputfileELectrons, const TString inputfileMuons, const TString inputfilePions, const int energy){
    
    TFile* out = new TFile(Form("PhotonsPerCathode_%dGeV.root", energy), "RECREATE");

    gStyle->SetGridColor(kGray);
    gStyle->SetGridStyle(3);   // dashed
    gStyle->SetGridWidth(1);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    gROOT->ForceStyle();
    TGaxis::SetMaxDigits(3);
    
    // 1D photon distributions for cathode 1 and cathode 2, separated by particle type.
    TH1D *histElectrons1 = new TH1D(" ", " ", 80, 0, 1000000);
    histElectrons1->GetXaxis()->SetTitle("Average photons on cathode 1");
    histElectrons1->GetYaxis()->SetTitle("Counts");
    histElectrons1->SetTitle(Form("%d GeV", energy));
    histElectrons1->SetTitleSize(0.04);
    histElectrons1->SetLineWidth(2);
    histElectrons1->SetFillStyle(0);
    histElectrons1->SetLineStyle(1);
    histElectrons1->SetLineColor(kBlue+1);
    TH1D *histElectrons2 = new TH1D(" ", " ", 80, 0, 1000000);
    histElectrons2->GetXaxis()->SetTitle("Average photons on cathode 2");
    histElectrons2->GetYaxis()->SetTitle("Counts");
    histElectrons2->SetTitle(Form("%d GeV", energy));
    histElectrons2->SetTitleSize(0.04);
    histElectrons2->SetLineWidth(2);
    histElectrons2->SetFillStyle(0);
    histElectrons2->SetLineStyle(1);
    histElectrons2->SetLineColor(kBlue+1);
    
    TH1D *histMuons1 = new TH1D(" ", " ", 80, 0, 1000000);
    histMuons1->SetLineWidth(2);
    histMuons1->SetFillStyle(0);
    histMuons1->SetLineStyle(1);
    histMuons1->SetLineColor(kOrange+7);
    TH1D *histMuons2 = new TH1D(" ", " ", 80, 0, 1000000);
    histMuons2->GetXaxis()->SetTitle("Average photons on cathode 2");
    histMuons2->SetLineWidth(2);
    histMuons2->SetFillStyle(0);
    histMuons2->SetLineStyle(1);
    histMuons2->SetLineColor(kOrange+7);
    
    TH1D *histPions1 = new TH1D(" ", " ", 80, 0, 1000000);
    histPions1->SetLineWidth(2);
    histPions1->SetFillStyle(0);
    histPions1->SetLineStyle(1);
    histPions1->SetLineColor(kGreen+2);
    TH1D *histPions2 = new TH1D(" ", " ", 80, 0, 1000000);
    histPions2->SetLineWidth(2);
    histPions2->SetFillStyle(0);
    histPions2->SetLineStyle(1);
    histPions2->SetLineColor(kGreen+2);
    
    TH1D *histAll1 = new TH1D(" ", " ", 80, 0, 1000000);
    histAll1->SetLineWidth(2);
    histAll1->SetFillStyle(0);
    histAll1->SetLineStyle(2);
    histAll1->SetLineColor(kBlack);
    
    TH1D *histAll2 = new TH1D(" ", " ", 80, 0, 1000000);
    histAll2->SetLineWidth(2);
    histAll2->SetFillStyle(0);
    histAll2->SetLineStyle(2);
    histAll2->SetLineColor(kBlack);
    
    // ============================ ELECTRONS ============================   
    TFile *inputElectrons = new TFile(inputfileELectrons, "READ");
    TTree *photonTreeElectrons = (TTree*)inputElectrons->Get("photoncounter/PhotonsTree");
    // Select the variables we need
    int nphotonsElectrons; photonTreeElectrons->SetBranchAddress("number_photons", &nphotonsElectrons);
    int chElectrons;       photonTreeElectrons->SetBranchAddress("ch", &chElectrons);
    int eventElectrons;    photonTreeElectrons->SetBranchAddress("event", &eventElectrons);
    TGraph *gr1 = new TGraph();
    gr1->SetMarkerStyle(24);
    gr1->SetMarkerColor(kBlue+1);
    // Number of events
    int nentriesElectrons = photonTreeElectrons->GetEntries();
    //Variables to store the total number of photons landing on PDs in each cathode
    double ncat1Electrons=0, ncat2Electrons=0;
    // Loop through the events and fill histogram
    for (int i = 0; i < nentriesElectrons; i++) {
        photonTreeElectrons->GetEntry(i);
        if(chElectrons==4 || chElectrons==6 || chElectrons==8 || chElectrons==10){
            ncat1Electrons+=nphotonsElectrons;
        }
        if(chElectrons==5 || chElectrons==7 || chElectrons==9 || chElectrons==11){
            ncat2Electrons+=nphotonsElectrons;
        }
     if((i+1)%20==0&&i!=0){
       gr1->SetPoint((i+1)/20, ncat1Electrons, ncat2Electrons);
       histElectrons1->Fill(ncat1Electrons/4.0);
       histAll1->Fill(ncat1Electrons/4.0);
       histElectrons2->Fill(ncat2Electrons/4.0);
       histAll2->Fill(ncat2Electrons/4.0);
       ncat1Electrons=0; ncat2Electrons=0;
       } //Starting a new event. There are 20 entries per event
     }
     
    // ============================== MUONS ==============================   
    TFile *inputMuons = new TFile(inputfileMuons, "READ");
    TTree *photonTreeMuons = (TTree*)inputMuons->Get("photoncounter/PhotonsTree");
    // Select the variables we need
    int nphotonsMuons; photonTreeMuons->SetBranchAddress("number_photons", &nphotonsMuons);
    int chMuons;       photonTreeMuons->SetBranchAddress("ch", &chMuons);
    int eventMuons;    photonTreeMuons->SetBranchAddress("event", &eventMuons);
    TGraph *gr2 = new TGraph();
    gr2->SetMarkerStyle(25);
    gr2->SetMarkerColor(kOrange+7);
    // Number of events
    int nentriesMuons = photonTreeMuons->GetEntries();
    // Variables to store the total number of photons landing on PDs in each cathode
    double ncat1Muons=0, ncat2Muons=0;
    
    // Loop through the events and fill histogram
    for (int i = 0; i < nentriesMuons; i++) {
        photonTreeMuons->GetEntry(i);
        if(chMuons==4 || chMuons==6 || chMuons==8 || chMuons==10){
            ncat1Muons+=nphotonsMuons;
        }
        if(chMuons==5 || chMuons==7 || chMuons==9 || chMuons==11){
            ncat2Muons+=nphotonsMuons;
        }
     if((i+1)%20==0&&i!=0){
       gr2->SetPoint((i+1)/20, ncat1Muons, ncat2Muons);
       histMuons1->Fill(ncat1Muons/4.0);
       histAll1->Fill(ncat1Muons/4.0);
       histMuons2->Fill(ncat2Muons/4.0);
       histAll2->Fill(ncat2Muons/4.0);
       ncat1Muons=0; ncat2Muons=0;
       } //Starting a new event. There are 20 entries per event
     }
    
     
     
    // ============================== PIONS ==============================   
    // Load the hist file we produced
    TFile *inputPions = new TFile(inputfilePions, "READ");
    TTree *photonTreePions = (TTree*)inputPions->Get("photoncounter/PhotonsTree");
    // Select the variables we need
    int nphotonsPions; photonTreePions->SetBranchAddress("number_photons", &nphotonsPions);
    int chPions;       photonTreePions->SetBranchAddress("ch", &chPions);
    int eventPions;    photonTreePions->SetBranchAddress("event", &eventPions);
    TGraph *gr3 = new TGraph();
    gr3->SetMarkerStyle(26);
    gr3->SetMarkerColor(kGreen+2);
    // Number of events
    int nentriesPions = photonTreePions->GetEntries();
    // Variables to store the total number of photons landing on PDs in each cathode
    double ncat1Pions=0, ncat2Pions=0;
    
    // Loop through the events and fill histogram
    for (int i = 0; i < nentriesPions; i++) {
        photonTreePions->GetEntry(i);
        if(chPions==4 || chPions==6 || chPions==8 || chPions==10){
            ncat1Pions+=nphotonsPions;
        }
        if(chPions==5 || chPions==7 || chPions==9 || chPions==11){
            ncat2Pions+=nphotonsPions;
        }
     if((i+1)%20==0&&i!=0){
       gr3->SetPoint((i+1)/20, ncat1Pions, ncat2Pions);
       histPions1->Fill(ncat1Pions/4.0);
       histAll1->Fill(ncat1Pions/4.0);
       histPions2->Fill(ncat2Pions/4.0);
       histAll2->Fill(ncat2Pions/4.0);
       ncat1Pions=0; ncat2Pions=0;
       } //Starting a new event. There are 20 entries per event
     }
     // ================= MUON SELECTION OPTIMIZATION =================
     //CONFUSION MATRIX:
    int nMuonTotal = gr2->GetN();
    int nElectronTotal = gr1->GetN();
    int nPionTotal = gr3->GetN();
    int nNotMuonTotal = nElectronTotal + nPionTotal;

    int TP = 0;
    int FP = 0;
    int FN = 0;
    int TN = 0;

    double best_score = 0;
    double best_mux, best_muy, best_a, best_b, best_theta, best_eff, best_pur;
    vector<Point> X_mu, X_e, X_pi;
    
    // Copy muon TGraph points into X_mu
    for (int i = 0; i < gr2->GetN(); i++) {
        double x,y;
        gr2->GetPoint(i,x,y);
        X_mu.push_back({x,y});
    }
    
    // Copy electron TGraph points into X_e
    for (int i = 0; i < gr1->GetN(); i++) {
        double x,y;
        gr1->GetPoint(i,x,y);
        X_e.push_back({x,y});
    }
    
    // Copy pion TGraph points into X_pi
    for (int i = 0; i < gr3->GetN(); i++) {
        double x,y;
        gr3->GetPoint(i,x,y);
        X_pi.push_back({x,y});
    }
    
    double xmin = 1e30, xmax = -1;
    double ymin = 1e30, ymax = -1;

    for (int i = 0; i < gr2->GetN(); i++) {
        double x,y;
        gr2->GetPoint(i, x, y);

        if (x < xmin) xmin = x;
        if (x > xmax) xmax = x;

        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }
    
    // Coarse scan grids for ellipse semi-axes, center coordinates, and rotation angle
    vector<double> a_vals  = linspace(1e4, 8e5, 20);
    vector<double> b_vals  = linspace(1e4, 8e5, 20);
    vector<double> mux_vals = linspace(xmin, xmax, 20);
    vector<double> muy_vals = linspace(ymin, ymax, 20);
    vector<double> theta_vals = linspace(-90*TMath::DegToRad(), 90*TMath::DegToRad(), 10);
    
    auto [a_min, a_max] = std::minmax_element(a_vals.begin(), a_vals.end());
    auto [b_min, b_max] = std::minmax_element(b_vals.begin(), b_vals.end());

    TH2F *h_ab = new TH2F("h_ab",";a;b", 20, *a_min, *a_max, 20, *b_min, *b_max);

    auto [mux_min, mux_max] = std::minmax_element(mux_vals.begin(), mux_vals.end());
    auto [muy_min, muy_max] = std::minmax_element(muy_vals.begin(), muy_vals.end());

    TH2F *h_mu = new TH2F("h_mu",";x_center; y_center", 20, *mux_min, *mux_max, 20, *muy_min, *muy_max);
        
    for (double mux : mux_vals) {
      for (double muy : muy_vals) {
        for (double a : a_vals) {
          for (double b : b_vals) {
            for (double theta : theta_vals) {
            
            // Evaluate efficiency and purity for this candidate ellipse
            double eff, pur;
            compute_metrics(mux, muy, a, b, theta,
                            eff, pur,
                            X_mu, X_e, X_pi);
            
            // Optimization score used throughout the scan
            double score = eff * pur;
            
            std::cout << "score=" << score 
           << "  mux=" << mux 
           << "  muy=" << muy
           << "  a=" << a 
           << "  b=" << b
           << "  theta=" << theta << std::endl;
            
            int bin_ab_x = h_ab->GetXaxis()->FindBin(a);
            int bin_ab_y = h_ab->GetYaxis()->FindBin(b);

            int bin_mu_x = h_mu->GetXaxis()->FindBin(mux);
            int bin_mu_y = h_mu->GetYaxis()->FindBin(muy);

            // keep the BEST score per bin
            double old_ab = h_ab->GetBinContent(bin_ab_x, bin_ab_y);
            if (score > old_ab)
                h_ab->SetBinContent(bin_ab_x, bin_ab_y, score);

            double old_mu = h_mu->GetBinContent(bin_mu_x, bin_mu_y);
            if (score > old_mu)
                h_mu->SetBinContent(bin_mu_x, bin_mu_y, score);
            
            // Update the globally best ellipse when a better score is found
            if (score > best_score) {
              best_score = score;
              best_mux = mux;
              best_muy = muy;
              best_a = a;
              best_b = b;
              best_theta=theta;
              best_eff=eff;
              best_pur=pur;
            }
          }
        }
      }
    }
    }
    
    std::cout << "best score=" << best_score 
         << "  best_mux=" << best_mux 
         << "  best_muy=" << best_muy
         << "  best_a=" << best_a 
         << "  best_b=" << best_b
         << "  best_theta=" << best_theta 
         << "  best_eff=" << best_eff 
         << "  best_pur=" << best_pur << std::endl;
    
    EllipseObjective obj;
    obj.X_mu = &X_mu;
    obj.X_e  = &X_e;
    obj.X_pi = &X_pi;
    obj.theta = best_theta;

    ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer("Minuit2", "Simplex");

    ROOT::Math::Functor f(obj, 5);

    min->SetFunction(f);

    min->SetVariable(0, "mux", best_mux, 1e3);
    min->SetVariable(1, "muy", best_muy, 1e3);
    min->SetVariable(2, "a",   best_a,   1e3);
    min->SetVariable(3, "b",   best_b,   1e3);
    min->SetVariable(4, "theta", best_theta, 1e-3);
    
    min->SetVariableLimits(0, 1e3, 5e5);   // a
    min->SetVariableLimits(1, 1e3, 5e5);
    min->SetVariableLimits(2, 1e3, 5e5);   // a
    min->SetVariableLimits(3, 1e3, 5e5);   // b
    min->SetVariableLimits(4, -90*TMath::DegToRad(), 90*TMath::DegToRad()); // theta

    min->Minimize();
    
    // Classify muons with the selected ellipse to build the final confusion matrix
    for (auto &p : X_mu) {
        if (inEllipse(p.x, p.y, best_mux, best_muy, best_a, best_b, best_theta))
            TP++;
        else
            FN++;
    }
    
    // Electrons selected by the ellipse are false positives
    for (auto &p : X_e) {
        if (inEllipse(p.x, p.y, best_mux, best_muy, best_a, best_b, best_theta))
            FP++;
        else
            TN++;
    }
    
    // Pions selected by the ellipse are also false positives
    for (auto &p : X_pi) {
        if (inEllipse(p.x, p.y, best_mux, best_muy, best_a, best_b, best_theta))
            FP++;
        else
            TN++;
    }
    
    TEfficiency eff_teff("eff","",1,0,1);
    eff_teff.SetTotalEvents(1, TP + FN);
    eff_teff.SetPassedEvents(1, TP);

    double eff     = eff_teff.GetEfficiency(1);
    double eff_lo  = eff_teff.GetEfficiencyErrorLow(1);
    double eff_hi  = eff_teff.GetEfficiencyErrorUp(1);
    
    TEfficiency pur_teff("pur","",1,0,1);
    pur_teff.SetTotalEvents(1, TP + FP);
    pur_teff.SetPassedEvents(1, TP);

    double pur     = pur_teff.GetEfficiency(1);
    double pur_lo  = pur_teff.GetEfficiencyErrorLow(1);
    double pur_hi  = pur_teff.GetEfficiencyErrorUp(1);
    
    std::cout << "Efficiency = " << eff
          << " +" << eff_hi
          << " -" << eff_lo << std::endl;

    std::cout << "Purity = " << pur
            << " +" << pur_hi
            << " -" << pur_lo << std::endl;

    out->cd();
    
    
    TCanvas *c1 = new TCanvas("c1", "Number of photons 5GeV");
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.08);
    c1->SetTopMargin(0.08);
    c1->SetGrid();
    histElectrons1->GetXaxis()->SetTitleSize(0.045);
    histElectrons1->GetYaxis()->SetTitleSize(0.045);
    histElectrons1->GetXaxis()->SetLabelSize(0.045);
    histElectrons1->GetYaxis()->SetLabelSize(0.045);
    histElectrons1->GetYaxis()->SetTitleOffset(1.2);
    histElectrons1->GetXaxis()->SetTitleOffset(1.2);
    
    histElectrons1->Draw("HIST");
    histMuons1->Draw("HIST SAME");
    histPions1->Draw("HIST SAME");
    histAll1->Draw("HIST SAME");
    histElectrons1->GetYaxis()->SetRangeUser(0, 500);
    auto legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend->AddEntry(histElectrons1, "Electron", "l");
    legend->AddEntry(histMuons1, "Muon", "l");
    legend->AddEntry(histPions1, "Pion", "l");
    legend->AddEntry(histAll1, "All particles", "l");
    legend->SetBorderSize(0);
    legend->SetFillStyle(1001);
    legend->SetTextSize(0.04);
    legend->Draw();
    c1->Write(Form("HistCathode1%dGeV.pdf", energy));
    
  
    
    TCanvas *c2 = new TCanvas("c2", "Number of photons 5GeV");
    c2->SetLeftMargin(0.15);
    c2->SetBottomMargin(0.15);
    c2->SetRightMargin(0.08);
    c2->SetTopMargin(0.08);
    c2->SetGrid();
    histMuons2->GetXaxis()->SetTitleSize(0.045);
    histMuons2->GetYaxis()->SetTitleSize(0.045);
    histMuons2->GetXaxis()->SetLabelSize(0.045);
    histMuons2->GetYaxis()->SetLabelSize(0.045);
    histMuons2->GetYaxis()->SetTitleOffset(1.2);
    histMuons2->GetXaxis()->SetTitleOffset(1.2);
    histMuons2->Draw("HIST");
    histElectrons2->Draw("HIST SAME");
    histPions2->Draw("HIST SAME");
    histAll2->Draw("HIST SAME");
    histMuons2->GetYaxis()->SetRangeUser(0, 500);
    auto legend2 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend2->AddEntry(histElectrons2, "Electron", "l");
    legend2->AddEntry(histMuons2, "Muon", "l");
    legend2->AddEntry(histPions2, "Pion", "l");
    legend2->AddEntry(histAll2, "All particles", "l");
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(1001);
    legend2->SetTextSize(0.04);
    legend2->Draw();
    c2->Write(Form("HistCathode2%dGeV.pdf", energy));
    
    
    TCanvas *c3 = new TCanvas("c3", "Number of photons 5GeV");
    c3->SetLeftMargin(0.15);
    c3->SetBottomMargin(0.15);
    c3->SetRightMargin(0.15);
    c3->SetTopMargin(0.08);
    c3->SetCanvasSize(500,500);
    
    
     gr2->GetXaxis()->SetTitle("Photons on cathode 1");
     gr2->GetYaxis()->SetTitle("Photons on cathode 2");
     gr2->SetTitle(Form("%d GeV", energy));
     gr2->GetXaxis()->SetTitleSize(0.045); 
     gr2->GetYaxis()->SetTitleSize(0.045);             
     gr2->GetXaxis()->SetLabelSize(0.045);
     gr2->GetYaxis()->SetLabelSize(0.045);
     gr2->GetYaxis()->SetTitleOffset(1.5);
     gr2->GetXaxis()->SetTitleOffset(1.2);
     gr2->GetXaxis()->SetLimits(0, 3000000);
     gr2->GetYaxis()->SetRangeUser(0, 3000000);
     gr2->Draw("AP");   // A = axes, P = points
     gr3->Draw("P SAME");
     gr1->Draw("P SAME");
     gr2->GetXaxis()->SetLimits(0, 3e6);
     gr2->SetMinimum(0);
     gr2->SetMaximum(3e6);

      // NOW enforce aspect ratio
      gPad->SetFixedAspectRatio();
      gPad->Update();
 
      TGraph *ellipse = new TGraph(200);
      for (int i = 0; i < 200; i++) {
          double t = 2*TMath::Pi()*i/200.;

          double x_local = best_a * cos(t);
          double y_local = best_b * sin(t);

          double x = best_mux + cos(best_theta)*x_local - sin(best_theta)*y_local;
          double y = best_muy + sin(best_theta)*x_local + cos(best_theta)*y_local;

          ellipse->SetPoint(i,x,y);
      }
      ellipse->SetLineColor(kBlack);
      ellipse->SetLineWidth(2);
      ellipse->Draw("L SAME");
      
      auto legend3 = new TLegend(0.6, 0.7, 0.8, 0.9);
       legend3->AddEntry(gr1, "Electron", "p");
       legend3->AddEntry(gr2, "Muon", "p");
       legend3->AddEntry(gr3, "Pion", "p");
       legend3->SetBorderSize(1);
       legend3->SetFillStyle(0);
       legend3->SetTextSize(0.03);
       legend3->Draw();

        c3->cd();
        c3->Modified();
        c3->Update();

        out->cd();
        c3->Write(Form("APA1vsAPA2%dGeV.pdf", energy)); 
     
     
    gStyle->SetNumberContours(10); // smoother gradient
    
    
    TCanvas *c4 = new TCanvas("c4", "");
    c4->SetLeftMargin(0.15);
    c4->SetBottomMargin(0.15);
    c4->SetRightMargin(0.15);
    c4->SetTopMargin(0.08);
    
    TH2F *hConf = new TH2F("hConf","",
                      2,0,2,2,0,2);
    
    int total_muon = TP + FN;
    int total_not  = FP + TN;

    hConf->SetBinContent(1,2, 100.0*TP/total_muon); // correctly selected muons
    hConf->SetBinContent(2,2, 100.0*FN/total_muon); // missed muons

    hConf->SetBinContent(1,1, 100.0*FP/total_not);  // false positives
    hConf->SetBinContent(2,1, 100.0*TN/total_not);  // correctly rejected
    hConf->GetXaxis()->SetBinLabel(1,"Muon");
    hConf->GetXaxis()->SetBinLabel(2,"Not Muon");
    hConf->GetYaxis()->SetBinLabel(1,"Rejected");
    hConf->GetYaxis()->SetBinLabel(2,"Selected");
    hConf->GetXaxis()->SetTitleSize(0.055);
    hConf->GetYaxis()->SetTitleSize(0.055);
    hConf->GetZaxis()->SetTitleSize(0.055);
    hConf->GetXaxis()->SetLabelSize(0.055);
    hConf->GetYaxis()->SetLabelSize(0.055);
    hConf->GetZaxis()->SetLabelSize(0.055);
    hConf->GetYaxis()->SetTitleOffset(1.2);
    hConf->GetXaxis()->SetTitleOffset(1.2);
    
    hConf->SetMinimum(0);
    hConf->SetMaximum(50); 
    
    gStyle->SetPalette(kViridis);
    gStyle->SetPaintTextFormat("4.1f");

    hConf->GetXaxis()->CenterLabels();
    hConf->GetYaxis()->CenterLabels();

    hConf->SetMarkerSize(3.0);

    hConf->Draw("COLZ TEXT");
    c4->Write(Form("ConfusionMatrix%dGeV.pdf", energy)); 
    
    
    TCanvas *c_ab = new TCanvas("c_ab","",600,600);
    c_ab->SetLeftMargin(0.15);
    c_ab->SetBottomMargin(0.15);
    c_ab->SetRightMargin(0.2);
    c_ab->SetTopMargin(0.08);
    
    h_ab->GetXaxis()->SetTitleSize(0.045);
    h_ab->GetYaxis()->SetTitleSize(0.045);
    h_ab->GetXaxis()->SetLabelSize(0.045);
    h_ab->GetYaxis()->SetLabelSize(0.045);
    h_ab->GetYaxis()->SetTitleOffset(1.2);
    h_ab->GetXaxis()->SetTitleOffset(1.2);
    h_ab->Draw("COLZ");
    c_ab->Write(Form("Heatmap_ab%dGeV.pdf", energy)); 
    
    
    
    TCanvas *c_mu = new TCanvas("c_mu","",600,600);
    c_mu->SetLeftMargin(0.15);
    c_mu->SetBottomMargin(0.15);
    c_mu->SetRightMargin(0.2);
    c_mu->SetTopMargin(0.08);
    h_mu->GetXaxis()->SetTitleSize(0.045);
    h_mu->GetYaxis()->SetTitleSize(0.045);
    h_mu->GetXaxis()->SetLabelSize(0.045);
    h_mu->GetYaxis()->SetLabelSize(0.045);
    h_mu->GetYaxis()->SetTitleOffset(1.2);
    h_mu->GetXaxis()->SetTitleOffset(1.2);
    h_mu->Draw("COLZ");
    c_mu->Write(Form("Heatmap_center%dGeV.pdf", energy));    
                      
     inputElectrons->Close();
     inputMuons->Close();
     inputPions->Close();
     
 }

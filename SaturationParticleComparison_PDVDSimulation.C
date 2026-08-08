/*
 * =============================================================================
 * Saturation Summary Macro
 * =============================================================================
 *
 * Purpose: read the output of the saturation analyzer for three
 * particle samples at the same beam energy and compares:
 *
 *   1. the number of saturated waveforms per event;
 *   2. the total number of waveforms per event;
 *   3. the number of saturated bins per saturated waveform.
 *
 *  -v5-sat.root files for the same energy but different particles
 *   inputfile1 -> electrons
 *   inputfile2 -> pions
 *   inputfile3 -> muons
 *
 * For each input file, the macro reads:
 *
 *   saturationcounter/CountWaveforms
 *
 * and uses the branches:
 *
 *   EventID          : event number
 *   nwaveformsSat    : number of saturated waveforms in the event
 *   nwaveformsTotal  : total number of waveforms in the event
 *   binsatVector     : number of saturated ADC bins for each saturated waveform
 *
 * Three histograms are produced for each particle sample:
 *
 *   hist[f]       : distribution of saturated waveforms per event
 *   histtotal[f]  : distribution of total waveforms per event
 *   histBins[f]   : distribution of saturated-bin counts per waveform
 *
 * Written by Maressa Sampaio and Laura Paulucci
 * maressap@ifi.unicamp.br
 *
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
#include <vector>

void SatData(const TString inputfile1, const TString inputfile2, const TString inputfile3, const int energy){
  
    TString files[3] = {inputfile1, inputfile2, inputfile3};
    TH1D *hist[3];        // saturated
    TH1D *histtotal[3];  // total
    TH1D *histBins[3]; 
    
    gStyle->SetTitleFontSize(0.045);
    gStyle->SetOptStat(0);
    gStyle->SetCanvasPreferGL(kTRUE);
    
    
    for (int f = 0; f < 3; f++) {

    TFile *input = new TFile(files[f], "READ");
    TTree *CountWaveforms = (TTree*)input->Get("saturationcounter/CountWaveforms");
    // Select the variables we need
    int Event;                CountWaveforms->SetBranchAddress("EventID", &Event);
    int nwaveformsSat;        CountWaveforms->SetBranchAddress("nwaveformsSat", &nwaveformsSat);
    int nwaveformsTotal;      CountWaveforms->SetBranchAddress("nwaveformsTotal", &nwaveformsTotal);
    float BinsatVector[10000]; CountWaveforms->SetBranchAddress("binsatVector", &BinsatVector);
    int nentries = CountWaveforms->GetEntries();

    hist[f] = new TH1D(Form("hist_sat_%d", f), "", 150, 0, 160);
    histtotal[f] = new TH1D(Form("hist_tot_%d", f), "", 150, 0, 160);
    histBins[f] = new TH1D(Form("hist_bins_%d", f), "", 150, 0, 600);



    for (int i = 0; i < nentries; i++) {
        CountWaveforms->GetEntry(i);

        hist[f]->Fill(nwaveformsSat);
        histtotal[f]->Fill(nwaveformsTotal);

        for (int j = 0; j < nwaveformsSat; j++) {
            histBins[f]->Fill(BinsatVector[j]);
        }
    }
}

    int colorshist[3] = {kOrange+3, kOrange+7, kOrange+0};
    int colorstotal[3] = {kBlue+3, kAzure, kCyan};
    int colorsbins[3] = {kGreen+3, kTeal-7, kGreen+0};

    for (int f = 0; f < 3; f++) {
        hist[f]->SetFillStyle(1001);
        hist[f]->SetLineWidth(2);
        hist[f]->SetLineStyle(1);
        hist[f]->SetFillColorAlpha(colorshist[f], 0.5);
        hist[f]->SetLineColor(colorshist[f]);
        hist[f]->GetXaxis()->SetTitle("Number of waveforms per event");
        hist[f]->GetYaxis()->SetTitle("Number of events");
        hist[f]->GetXaxis()->SetTitleSize(0.045);
        hist[f]->GetYaxis()->SetTitleSize(0.045);
        hist[f]->GetXaxis()->SetLabelSize(0.045);
        hist[f]->GetYaxis()->SetLabelSize(0.045);
        hist[f]->GetYaxis()->SetTitleOffset(1.2);
        hist[f]->GetXaxis()->SetTitleOffset(1.2);
        
        histtotal[f]->SetFillStyle(1001);
        histtotal[f]->SetFillColorAlpha(colorstotal[f], 0.5);
        histtotal[f]->SetLineColor(colorstotal[f]);
        histtotal[f]->SetLineStyle(1); 
        histtotal[f]->SetLineWidth(2);
        histtotal[f]->GetXaxis()->SetTitle("Number of waveforms per event");
        histtotal[f]->GetYaxis()->SetTitle("Number of events");
        histtotal[f]->GetXaxis()->SetTitleSize(0.045);
        histtotal[f]->GetYaxis()->SetTitleSize(0.045);
        histtotal[f]->GetXaxis()->SetLabelSize(0.045);
        histtotal[f]->GetYaxis()->SetLabelSize(0.045);
        histtotal[f]->GetYaxis()->SetTitleOffset(1.2);
        histtotal[f]->GetXaxis()->SetTitleOffset(1.2);

        histBins[f]->SetLineWidth(2);
        histBins[f]->SetFillStyle(1001);
        histBins[f]->SetFillColorAlpha(colorsbins[f], 0.5);
        histBins[f]->SetLineColor(colorsbins[f]);
        histBins[f]->SetLineStyle(1); 
        histBins[f]->GetXaxis()->SetTitle("Number of saturated bins per waveform");
        histBins[f]->GetYaxis()->SetTitle("Number of waveforms");
        histBins[f]->GetXaxis()->SetTitleSize(0.045);
        histBins[f]->GetYaxis()->SetTitleSize(0.045);
        histBins[f]->GetXaxis()->SetLabelSize(0.045);
        histBins[f]->GetYaxis()->SetLabelSize(0.045);
        histBins[f]->GetYaxis()->SetTitleOffset(1.2);
        histBins[f]->GetXaxis()->SetTitleOffset(1.2);
    }
    
        TCanvas *c1 = new TCanvas("c1", "Number of saturated waveforms", 800, 600);
        c1->SetLeftMargin(0.15);
        c1->SetBottomMargin(0.15);
        c1->SetRightMargin(0.08);

        hist[0]->Draw("HIST");
        hist[1]->Draw("HIST SAME");
        hist[2]->Draw("HIST SAME");
        
        hist[0]->SetLineWidth(2);
        hist[1]->SetLineWidth(2);
        hist[2]->SetLineWidth(2);
        auto legend = new TLegend(0.8, 0.75, 0.95, 0.9);
        legend->AddEntry(hist[0], "Electron", "l");
        legend->AddEntry(hist[1], "Pion", "l");
        legend->AddEntry(hist[2], "Muon", "l");
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.04);        
        legend->Draw();
        c1->Print(Form("SatWaveforms%dGeV.pdf", energy));
        
        TCanvas *c3 = new TCanvas("c3", "Number of saturated waveforms", 800, 600);
        c3->SetLeftMargin(0.15);
        c3->SetBottomMargin(0.15);
        c3->SetRightMargin(0.08);

        histtotal[0]->Draw("HIST");
        histtotal[1]->Draw("HIST SAME");
        histtotal[2]->Draw("HIST SAME");
        
        auto legend3 = new TLegend(0.8, 0.75, 0.95, 0.9);
        legend3->AddEntry(histtotal[0], "Electron", "l");
        legend3->AddEntry(histtotal[1], "Pion", "l");
        legend3->AddEntry(histtotal[2], "Muon", "l");
        legend3->SetBorderSize(0);
        legend3->SetFillStyle(0);
        legend3->SetTextSize(0.04);

        hist[0]->SetLineWidth(2);
        hist[1]->SetLineWidth(2);
        hist[2]->SetLineWidth(2);
        legend3->Draw();
        c3->Print(Form("TotalWaveforms%dGeV.pdf", energy));

        TCanvas *c2 = new TCanvas("c2", "Bins saturated", 800, 600);
        c2->SetLeftMargin(0.15);
        c2->SetBottomMargin(0.15);  
        c2->SetRightMargin(0.08);

        histBins[2]->Draw("HIST");
        histBins[1]->Draw("HIST SAME");
        histBins[0]->Draw("HIST SAME");
        
        auto legend2 = new TLegend(0.8, 0.75, 0.95, 0.9);
        legend2->AddEntry(histBins[0], "Electron", "l");
        legend2->AddEntry(histBins[1], "Pion", "l");
        legend2->AddEntry(histBins[2], "Muon", "l");
        legend2->SetBorderSize(0);
        legend2->SetFillStyle(0);
        legend2->SetTextSize(0.04);

        histBins[0]->SetLineWidth(2);
        histBins[1]->SetLineWidth(2);
        histBins[2]->SetLineWidth(2);
        legend2->Draw();
        c2->Print(Form("BinsSat%dGeV.pdf", energy));

}

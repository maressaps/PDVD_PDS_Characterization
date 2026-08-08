/*
=============================================================================
 * ProtoDUNE-VD simulation light yield analysis
 * =============================================================================
 *
 * Purpose: Read photon-visibility maps separately for argon and xenon cases
 * and converts the optical visibility into expected light yield (LY).
 *
 * The analysis:
 *   - selects detector slices at fixed x, y, or z;
 *   - sums the optical visibility of the selected PD channels;
 *   - weights each channel by its efficiency;
 *   - converts visibility into light yield using the assumed scintillation
 *     photon yield in the active and buffer regions;
 *   - produces 2D LY maps for each detector slice;
 *   - averages the slices to obtain X, Y, and Z projected LY maps;
 *   - combines the argon and xenon contributions with the chosen mixture
 *     fractions;
 *   - calculates mean, minimum, and maximum LY values, including values
 *     restricted to the active volume projection.
 *
 * Main inputs
 * -----------
 *   input_fileAr : ROOT visibility file (with x, y, z and the visibility in each PD) for the argon component.
 *   input_fileXe : ROOT visibility file for the xenon component.
 *   pd           : photon-detector selection:
 *                    0 = all photon detectors
 *                    1 = X-Arapucas only
 *                    2 = PMTs only
 *
 * Main output
 * -----------
 *   LightYieldMaps_Ar&XeMixture.root, with X-, Y-, and Z-projection subdirectories.
 *
 * Important analysis parameters defined below include:
 *   - active region photon yield = 24000 photons/MeV;
 *   - buffer region photon yield = 40000 photons/MeV;
 *   - channel efficiencies;
 *   - detector/field-cage boundaries;
 *   - Ar and Xe mixture weights.
 *
 * Written by Maressa Sampaio and Laura Paulucci
 * maressap@ifi.unicamp.br
 *
 * =============================================================================
 */

#include <cstdio>
#include <iostream>
#include <algorithm>
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

/*
 * Helper function: make_LYplanes: fills one 2D light-yield histogram for a thin slice of the detector.
 *
 * plane  : coordinate held fixed (0=x, 1=y, 2=z)
 * plane2 : coordinate placed on the histogram y axis
 * plane3 : coordinate placed on the histogram x axis
 * position: center of the selected slice [cm]
 * pd     : 0=all PDs, 1=X-Arapucas only, 2=PMTs only
 * element: 0=argon component, 1=xenon component
 *
 * For each visibility point inside the selected slice, the function sums
 * the direct visibility over all PD channels after applying the corresponding
 * detection efficiency. The result is multiplied by the scintillation photon
 * yield and filled into the requested 2D projection.
 */
 
void make_LYplanes(TH2F* &hist, TTree* tAct, TTree* tBuf, int plane, int plane2, int plane3, double position, int n_pt_act, int n_pt_buf, int numPDs, int pd, int element){
  
  // Accumulator for the total efficiency-weighted visibility
  double total_visDirect;
  // Scintillation photon yields used to convert optical visibility into LY
  int photonYield=24000, photonYieldBuff=40000;
  // Define a thin +/-0.1 cm interval around the requested slice position
  double lowerlimit=position-0.1;
  double upperlimit=position+0.1;
  
  double coords[3];   tAct->SetBranchAddress("coords", coords);
  double coordsBuf[3];   tBuf->SetBranchAddress("coords", coordsBuf);
  double partial_visDirect[numPDs]; tAct->SetBranchAddress("opDet_visDirect", partial_visDirect);
  double partial_visBuf[numPDs]; tBuf->SetBranchAddress("opDet_visDirectBuff", partial_visBuf);
  
  // The efficiency vector has one value per channel.
  std::vector<float> eff;
  std::vector<double> xmin_tpc, xmax_tpc;
  // Nominal internal FC limits [cm].
  xmin_tpc={-350.0, -350.0, 0}, xmax_tpc={350.0, 350.0, 300.0}; // Field cage internal limits
  
  // Efficiency for argon
  if (element==0){
    if(pd==0){ //all PDs
      eff = { 0.03,  0.03,  0.03, 0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.03, 
              0.03,  0.03,  0.12, 0.036, 0.12,  0.12,  0.0,   0.03,  0.12,  0.036, 
              0.12,  0.12,  0.03, 0.03,  0.036, 0.036, 0.036, 0.036, 0.036, 0.0, 
              0.036, 0.036, 0.0,  0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.0 };
    }
      else if(pd==1){//just arapucas
        eff = { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 
                0.03, 0.03, 0.0,  0.0,  0.0,  0.0,  0.0,  0.03, 0.0,  0.0, 
                0.0,  0.0,  0.03, 0.03, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
                0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };
        }
        else{//just pmts
          eff = { 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
                  0.0,   0.0,   0.12, 0.036, 0.12,  0.12,  0.0,   0.0,   0.12,  0.036, 
                  0.12,  0.12,  0.0,  0.0,   0.036, 0.036, 0.036, 0.036, 0.036, 0.0, 
                  0.036, 0.036, 0.0,  0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.0 };
      
    }}
  
  // Efficiency for xenon
  if(element==1){
    if(pd==0){ //all PDs
      eff = { 0.03,  0.03,  0.03, 0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.03, 
              0.03,  0.03,  0.12, 0.036, 0.12,  0.12,  0.03,   0.03,  0.12,  0.036, 
              0.12,  0.12,  0.03, 0.03,  0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 
              0.036, 0.036, 0.0,  0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036 };
    }
      else if(pd==1){//just arapucas
        eff = { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 
                0.03, 0.03, 0.0,  0.0,  0.0,  0.0,  0.03,  0.03, 0.0,  0.0, 
                0.0,  0.0,  0.03, 0.03, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0, 
                0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0 };
        }
        else{//just pmts
          eff = { 0.0,   0.0,   0.0,  0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0, 
                  0.0,   0.0,   0.12, 0.036, 0.12,  0.12,  0.0,   0.0,   0.12,  0.036, 
                  0.12,  0.12,  0.0,  0.0,   0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 
                  0.036, 0.036, 0.0,  0.036, 0.036, 0.036, 0.036, 0.036, 0.036, 0.036 };
      
    }}

  
  // Loop over every point in the active visibility map and keep only those that lie inside the requested detector slice.
  for (int ievent = 0; ievent < n_pt_act; ievent++) {
      	tAct->GetEntry(ievent); 
        total_visDirect=0;
        if (coords[plane]>lowerlimit && coords[plane]<upperlimit){
    	  // for each accepted entry we need to sum the channels visibility
    	  for (int i_channel=0; i_channel < numPDs; i_channel++){
    	        total_visDirect += partial_visDirect[i_channel]*eff[i_channel];
    	        
          }
          total_visDirect= total_visDirect*photonYield;
          hist->Fill(coords[plane3], coords[plane2], total_visDirect );
          
        }}
   
   for (int ievent = 0; ievent < n_pt_buf; ievent++) {
      	tBuf->GetEntry(ievent); 
        total_visDirect=0;
        if (coordsBuf[plane]>lowerlimit && coordsBuf[plane]<upperlimit){
    	  for (int i_channel=0; i_channel < numPDs; i_channel++){
    	        total_visDirect += partial_visBuf[i_channel]*eff[i_channel]; 
          }
          total_visDirect= total_visDirect*photonYieldBuff;
          hist->Fill( coordsBuf[plane3], coordsBuf[plane2], total_visDirect );
          
        }}
   
               	        
}

/*
 * Main analysis function: protodunevd_vis
 * ---------------------------------------
 * input_fileAr : ROOT file containing the argon visibility map
 * input_fileXe : ROOT file containing the xenon visibility map
 * pd           : 0=all PDs, 1=X-Arapucas only, 2=PMTs only
 *
 * The function generates slice-by-slice LY maps for the three coordinate
 * directions, computes projected averages, combines the Ar and Xe components,
 * and evaluates summary statistics.
 */
void protodunevd_vis(const TString input_fileAr, const TString input_fileXe, int pd){

	
	// Open the two input ROOT files containing the Ar and Xe visibility maps
        TFile* finAr = new TFile(input_fileAr,"UPDATE");
        TFile* finXe = new TFile(input_fileXe,"UPDATE");
        
        // Load the necessary Trees
        // visibility for buffer region
        TTree* visbufAr =  (TTree*)finAr->Get("vismap/photoVisMapBuffer");
        TTree* visbufXe =  (TTree*)finXe->Get("vismap/photoVisMapBuffer");
        // visibility for active volume
        TTree* visactiveAr =  (TTree*)finAr->Get("vismap/photoVisMap");
        TTree* visactiveXe =  (TTree*)finXe->Get("vismap/photoVisMap");
        // Optical detectors parameters
        TTree* opdetAr =  (TTree*)finAr->Get("vismap/opDetMap");
        TTree* opdetXe =  (TTree*)finXe->Get("vismap/opDetMap");

        // -------------------------------------------------------------------
        // Spatial boundaries and binning
        // -------------------------------------------------------------------
        // xmin/xmax describe the full mapped cryostat region.
        // dx defines the nominal projection width in each direction [cm].
        std::vector<double> xmin, xmax, dx;
        xmin={-375.0, -427.4, -277.75}, xmax={415.0, 427.4, 577.05}, dx={10,10,10};
        xmin_tpc={-350.0, -350.0, 0}, xmax_tpc={350.0, 350.0, 300.0};
        cout << "dx= " << dx[0] << " " << dx[1] <<" " << dx[2] << endl;
        
        // Read the ROOT grid histograms used to obtain map boundaries and binning information.
        TH1D* h_tmp[3] = {0}; 
        h_tmp[0] = (TH1D*)fin->Get("vismap/hgrid0"); //X axis
        h_tmp[1] = (TH1D*)fin->Get("vismap/hgrid1"); //Y axis
        h_tmp[2] = (TH1D*)fin->Get("vismap/hgrid2"); //Z axis
        std::vector<int> nbins;
        nbins= {h_tmp[0]->GetNbinsX(), h_tmp[1]->GetNbinsX(), h_tmp[2]->GetNbinsX()};
        std::vector<double> xmin_tpc, xmax_tpc;
        xmin_tpc ={
          h_tmp[0]->GetXaxis()->GetXmin(), 
          h_tmp[1]->GetXaxis()->GetXmin(), 
          h_tmp[2]->GetXaxis()->GetXmin()};
        xmax_tpc ={
          h_tmp[0]->GetXaxis()->GetXmax(), 
          h_tmp[1]->GetXaxis()->GetXmax(), 
          h_tmp[2]->GetXaxis()->GetXmax()}; 
        std::vector<int> nbins;
        nbins={79, 85, 85};
        nbins[0] = (int)((xmax[0]-xmin[0])/dx[0]);
        nbins[1] = (int)((xmax[1]-xmin[1])/dx[1]);
        nbins[2] = (int)((xmax[2]-xmin[2])/dx[2]);
        cout << "nbins= " << nbins[0] << " " << nbins[1] <<" " << nbins[2] << endl;
        
        // Number of points/entries in each active and buffer visibility tree
        const int n_pt_Ar = visactiveAr->GetEntries();
        const int n_pt_Xe = visactiveXe->GetEntries();
        const int n_pt_Ar_buf = visbufAr->GetEntries();
        const int n_pt_Xe_buf = visbufXe->GetEntries();
    	int numPDs = 40;
        
        // Output ROOT file that will contain all maps and projections
        TFile* out = new TFile("LightYieldMaps_Ar&XeMixture.root", "RECREATE");
      
        // Relative Ar and Xe contributions used when constructing the final mixed-light maps. Concentration is equivalent to the proportion of light production considering a 17% reduction due to absorption of the argon component
        double concentration0 = 0.30;
        double concentration1 = 0.53;
    
  	// Create subdirectories
        TDirectory *ProjX = out->mkdir("ProjX");
        TDirectory *ProjY = out->mkdir("ProjY");
        TDirectory *ProjZ = out->mkdir("ProjZ");

        int plane=0, plane2=1, plane3=2;// LY(y,z) for fixed x
        double position; //position of the projection 	
    	ProjX->cd(); 
    	
   	// Histograms used for individual Ar/Xe slices and their accumulated projection averages
   	TCanvas *c1 =new TCanvas();
   	
   	TH2F *histXAr= new TH2F("histXAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histXAr->GetXaxis()->SetTitle("z (cm)");
        histXAr->GetYaxis()->SetTitle("y (cm)");
        
        TH2F *histXXe= new TH2F("histXXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histXXe->GetXaxis()->SetTitle("z (cm)");
        histXXe->GetYaxis()->SetTitle("y (cm)"); 
        
        TH2F *histAveXAr= new TH2F("histAveXAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveXAr->GetXaxis()->SetTitle("z (cm)");
        histAveXAr->GetYaxis()->SetTitle("y (cm)");   
        
        TH2F *histAveXXe= new TH2F("histAveXXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveXXe->GetXaxis()->SetTitle("z (cm)");
        histAveXXe->GetYaxis()->SetTitle("y (cm)"); 

   	TString fname;
        cout << "ProjX" << endl;
  	   	
	for(int i=0; i<nbins[plane];i++){
    	    position=xmin[plane]+(double(i)+0.5)*dx[plane];
            make_LYplanes(histXAr, visactiveAr, visbufAr, plane, plane2, plane3, position, n_pt_Ar, n_pt_Ar_buf, numPDs, pd, 0);
            fname.Form("LY_Ar_x=%g",position);
            histXAr->SetName(fname);
            histXAr->Draw("colz");
            histXAr->Write();
            histAveXAr->Add(histXAr);
            histXAr->Reset();
            

            make_LYplanes(histXXe, visactiveXe, visbufXe, plane, plane2, plane3, position, n_pt_Xe, n_pt_Xe_buf, numPDs, pd, 1); 
            fname.Form("LY_Xe_x=%g",position);
            histXXe->SetName(fname);
            histXXe->Draw("colz");
            histXXe->Write();
            histAveXXe->Add(histXXe);
            histXXe->Reset();
            }
        
        // Divide the accumulated argon map by the number of x slices to obtain the projected mean LY(y,z).
	histAveXAr->Scale(1./nbins[plane]);
	histAveXAr->Write();
	
	histAveXXe->Scale(1./nbins[plane]);
	histAveXXe->Write();
	
	
	TH2F *mixAveX= new TH2F("mixAveX", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
	mixAveX->GetXaxis()->SetTitle("z (cm)");
        mixAveX->GetYaxis()->SetTitle("y (cm)"); 
        
        TH2F* hArX_scaled = (TH2F*)histAveXAr;
        TH2F* hXeX_scaled = (TH2F*)histAveXXe;
        
        // Apply the Ar and Xe mixture weights.
        hArX_scaled->Scale(concentration0);
        hXeX_scaled->Scale(concentration1);

        mixAveX->Add(hArX_scaled, hXeX_scaled); 
        mixAveX->Write(); 
        
        int nBinsX = mixAveX->GetNbinsX();
        int nBinsY = mixAveX->GetNbinsY();

        double sum = 0; //active volume only
        double sum2 =0; //total volume
        double max = -100000;
        double min = 100000;
        int total = 0;
        int total2 = 0;
        double binX;
        double binY;

        for (int i = 1; i <= nBinsX; ++i) {
            for (int j = 1; j <= nBinsY; ++j) {
            binX=mixAveX->GetXaxis()->GetBinCenter(i);
            binY=mixAveX->GetYaxis()->GetBinCenter(j);
            sum2 += mixAveX->GetBinContent(i, j);
            total2++;
            if(binX > xmin_tpc[2]  && binX < xmax_tpc[2] && binY > xmin_tpc[1]  && binY < xmax_tpc[1]){
                    sum += mixAveX->GetBinContent(i, j);
                    total++;
                    if (mixAveX->GetBinContent(i, j) < min){min=mixAveX->GetBinContent(i, j);}
                    if (mixAveX->GetBinContent(i, j) > max){max=mixAveX->GetBinContent(i, j);}
                }
            }}

        double mean= sum / total;
        double mean2= sum2 / total2;
        cout << "Mean LY value for the X projection active volume considering the Ar/Xe mix: " << mean << endl;
        cout << "Mean LY value for the X projection total volume considering the Ar/Xe mix: " << mean2 << endl;
        cout << "Maximum LY value for the X projection considering the Ar/Xe mix: " << max << endl;
        cout << "Minimum LY value for the X projection considering the Ar/Xe mix: " << min << endl;
        
        
        plane=1, plane2=0, plane3=2; //LY(x,z) for fixed y
        ProjY->cd(); 
   	
   	TH2F *histYAr= new TH2F("histYAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histYAr->GetXaxis()->SetTitle("z (cm)");
        histYAr->GetYaxis()->SetTitle("x (cm)");   
        
        TH2F *histYXe= new TH2F("histYXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histYXe->GetXaxis()->SetTitle("z (cm)");
        histYXe->GetYaxis()->SetTitle("x (cm)"); 
        
        
        TH2F *histAveYAr= new TH2F("histAveYAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveYAr->GetXaxis()->SetTitle("z (cm)");
        histAveYAr->GetYaxis()->SetTitle("x (cm)"); 
        
        TH2F *histAveYXe= new TH2F("histAveYXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveYXe->GetXaxis()->SetTitle("z (cm)");
        histAveYXe->GetYaxis()->SetTitle("x (cm)"); 

        cout << "ProjY" << endl;
	for(int i=0; i<nbins[plane];i++){
    	    position=xmin[plane]+(double(i)+0.5)*dx[plane];
    	    std::cout<<  "entrou1" <<endl;
            make_LYplanes(histYAr, visactiveAr, visbufAr, plane, plane2, plane3, position, n_pt_Ar, n_pt_Ar_buf, numPDs, pd, 0);
            fname.Form("LY_Ar_y=%g",position);
            histYAr->SetName(fname);
            histYAr->Draw("colz");
            histAveYAr->Add(histYAr);
            histYAr->Write();
            histYAr->Reset();
            
            make_LYplanes(histYXe, visactiveXe, visbufXe, plane, plane2, plane3, position, n_pt_Xe, n_pt_Xe_buf, numPDs, pd, 1); 
            
            fname.Form("LY_Xe_y=%g",position);
            histYXe->SetName(fname);
            histYXe->Draw("colz");
            histAveYXe->Add(histYXe);
            histYXe->Write();
            histYXe->Reset();
        }
	histAveYAr->Scale(1./nbins[plane]);
	histAveYAr->Write();
	
	histAveYXe->Scale(1./nbins[plane]);
	histAveYXe->Write();
	
	TH2F *mixAveY= new TH2F("mixAveY", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
	mixAveY->GetXaxis()->SetTitle("z (cm)");
        mixAveY->GetYaxis()->SetTitle("x (cm)"); 
        
        TH2F* hArY_scaled = (TH2F*)histAveYAr;
        TH2F* hXeY_scaled = (TH2F*)histAveYXe;
        
        hArY_scaled->Scale(concentration0);
        hXeY_scaled->Scale(concentration1);

        mixAveY->Add(hArY_scaled, hXeY_scaled); 
        mixAveY->Write();
        
        nBinsX = mixAveY->GetNbinsX();
        nBinsY = mixAveY->GetNbinsY();

        sum = 0;
        sum2 = 0;
        max = -100000;
        min = 100000;
        total = 0;
        total2 = 0;
        
       
        for (int i = 1; i <= nBinsX; ++i) {
            for (int j = 1; j <= nBinsY; ++j) {
                binX=mixAveY->GetXaxis()->GetBinCenter(i);
                binY=mixAveY->GetYaxis()->GetBinCenter(j);
                sum2 += mixAveY->GetBinContent(i, j);
                total2++;
                if(binX > xmin_tpc[2]  && binX < xmax_tpc[2] && binY > xmin_tpc[0]  && binY < xmax_tpc[0]){
                    sum += mixAveY->GetBinContent(i, j);
                    total++;
                    if (mixAveY->GetBinContent(i, j) < min){min=mixAveY->GetBinContent(i, j);}
                    if (mixAveY->GetBinContent(i, j) > max){max=mixAveY->GetBinContent(i, j);}
            }}}
        
        mean=0;
        mean2=0;
        mean= sum / total;
        mean2= sum2 / total2;
        cout << "Mean LY value for the Y projection active volume considering the Ar/Xe mix: " << mean << endl;
        cout << "Mean LY value for the Y projection total volume considering the Ar/Xe mix: " << mean2 << endl;
        cout << "Maximum LY value for the Y projection considering the Ar/Xe mix: " << max << endl;
        cout << "Minimum LY value for the Y projection considering the Ar/Xe mix: " << min << endl;
           

        plane=2, plane2=0, plane3=1; // LY(x,y) for fixed z
        ProjZ->cd(); 
   	
   	TH2F *histZAr= new TH2F("histZAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
    	histZAr->GetXaxis()->SetTitle("y (cm)");
    	histZAr->GetYaxis()->SetTitle("x (cm)");
    	
    	TH2F *histZXe= new TH2F("histZXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
    	histZXe->GetXaxis()->SetTitle("y (cm)");
    	histZXe->GetYaxis()->SetTitle("x (cm)");
        
        TH2F *histAveZAr= new TH2F("histAveZAr", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveZAr->GetXaxis()->SetTitle("y (cm)");
        histAveZAr->GetYaxis()->SetTitle("x (cm)");
        
        TH2F *histAveZXe= new TH2F("histAveZXe", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        histAveZXe->GetXaxis()->SetTitle("y (cm)");
        histAveZXe->GetYaxis()->SetTitle("x (cm)");
        
       
        cout << "ProjZ" << endl;
        for(int i=0; i<nbins[plane];i++){
    	    position=xmin[plane]+(double(i)+0.5)*dx[plane];
            make_LYplanes(histYAr, visactiveAr, visbufAr, plane, plane2, plane3, position, n_pt_Ar, n_pt_Ar_buf, numPDs, pd, 0);
            fname.Form("LY_Ar_z=%g",position);
            histZAr->SetName(fname);
            histZAr->Draw("colz");
            histAveZAr->Add(histZAr);
            histZAr->Write();
            histZAr->Reset();
            
            make_LYplanes(histZXe, visactiveXe, visbufXe, plane, plane2, plane3, position, n_pt_Xe, n_pt_Xe_buf, numPDs, pd, 1); 
            
            fname.Form("LY_Xe_z=%g",position);
            histZXe->SetName(fname);
            histZXe->Draw("colz");
            histAveZXe->Add(histZXe);
            histZXe->Write();
            histZXe->Reset();
        }
	histAveZAr->Scale(1./nbins[plane]);
	histAveZAr->Write();
	
	histAveZXe->Scale(1./nbins[plane]);
	histAveZXe->Write();
	
        
        TH2F *mixAveZ= new TH2F("mixAveZ", " ", nbins[plane3], xmin[plane3], xmax[plane3], nbins[plane2], xmin[plane2], xmax[plane2]);
        mixAveZ->GetXaxis()->SetTitle("y (cm)");
        mixAveZ->GetYaxis()->SetTitle("x (cm)"); 
              
        TH2F* hArZ_scaled = (TH2F*)histAveZAr;
        TH2F* hXeZ_scaled = (TH2F*)histAveZXe;
        
        hArZ_scaled->Scale(concentration0);
        hXeZ_scaled->Scale(concentration1);

        mixAveZ->Add(hArZ_scaled, hXeZ_scaled); 

        mixAveZ->Write();
        
        nBinsX = mixAveZ->GetNbinsX();
        nBinsY = mixAveZ->GetNbinsY();

        sum = 0;
        sum2 = 0;
        max = -100000;
        min = 100000;
        total = 0;
        total2 = 0;
        
        for (int i = 1; i <= nBinsX; ++i) {
            for (int j = 1; j <= nBinsY; ++j) {
                binX=mixAveZ->GetXaxis()->GetBinCenter(i);
                binY=mixAveZ->GetYaxis()->GetBinCenter(j);
                sum2 += mixAveZ->GetBinContent(i, j);
                total2++;
                //check if in active volume
                if(binX > xmin_tpc[1]  && binX < xmax_tpc[1] && binY > xmin_tpc[0]  && binY < xmax_tpc[0]){
                    sum += mixAveZ->GetBinContent(i, j);
                    total++;
                    if (mixAveZ->GetBinContent(i, j) < min){min=mixAveZ->GetBinContent(i, j);}
                    if (mixAveZ->GetBinContent(i, j) > max){max=mixAveZ->GetBinContent(i, j);}
            }}}
        
        mean =0;
        mean2 =0;
        mean= sum / total;
        mean2= sum2 / total2;
        cout << "Mean LY value for the Z projection active volume considering the Ar/Xe mix: " << mean << endl;
        cout << "Mean LY value for the Z projection total volume considering the Ar/Xe mix: " << mean2 << endl;
        cout << "Maximum LY value for the Z projection considering the Ar/Xe mix: " << max << endl;
        cout << "Minimum LY value for the Z projection considering the Ar/Xe mix: " << min << endl;
    	}



  	

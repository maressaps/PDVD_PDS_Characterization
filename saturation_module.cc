/*
 * =============================================================================
 * Channel Saturation Analyzer
 * =============================================================================
 *
 * Purpose: This ART/LArSoft analyzer studies saturation in waveforms.
 * For every event, it reads raw::OpDetWaveform objects, determines whether each
 * waveform reaches the channel ADC saturation threshold, counts the
 * number of saturated waveforms, and measures the number of consecutive
 * saturated ADC samples around the waveform maximum.
 *
 * The analyzer writes one entry per event to the TTree "CountWaveforms" with:
 *
 *   EventID         : ART event number
 *   nwaveformsSat   : number of saturated optical waveforms in the event
 *   nwaveformsTotal : total number of optical waveforms in the event
 *   binsatVector    : number of saturated ADC bins for each saturated waveform
 *
 * Saturation threshold
 * --------------------
 * The nominal threshold is calculated from the configured dynamic range:
 *
 *   max = 2^(DynamicRange)
 *
 * The threshold is shifted using the measured channel baseline relative to the reference baseline of
 * 1854 ADC:
 *
 *   max = 2^(DynamicRange) - (baseline[channel] - 1854)
 *
 * A waveform is classified as saturated when its maximum ADC value is greater
 * than or equal to this threshold.
 *
 * Saturated-bin counting
 * ----------------------
 * Starting from the waveform maximum, the analyzer counts consecutive samples
 * at or above the saturation threshold both before and after the peak. The
 * resulting width of the saturated plateau is stored in binsatVector.
 *
 * Configuration parameters
 * ------------------------
 *   OpDetWaveformLabel : input label for the raw::OpDetWaveform collection
 *   DynamicRange       : ADC dynamic range in bits
 *
 * Written by Maressa Sampaio and Laura Paulucci
 * maressap@ifi.unicamp.br
 *
 * =============================================================================
 */

#ifndef saturation_H
#define saturation_H 1

// ROOT includes
#include "TH1.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/PhotonBackTrackerService.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// ART includes.
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"

namespace opdet {

  class saturation : public art::EDAnalyzer{
  public:

    // Standard constructor and destructor for an ART module.
    saturation(const fhicl::ParameterSet&);
    virtual ~saturation();

    void beginJob();

    // Main analyzer routine, called once for each ART event.
    void analyze (const art::Event& evt);

  private:

    Int_t fEventID;
    std::string fOpDetWaveformLabel;
    float fDynamicRange;
    TTree * fCountTree;
    Int_t fnwaveformsSat;
    Int_t fnwaveformsTotal;
    Int_t fnbinsat;
    float fBinsatVector[10000];
    std::map<int, float> baselineMap;
  };

}

#endif // saturation_H

namespace opdet {

  //-----------------------------------------------------------------------
  // Constructor
  saturation::saturation(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Read configuration parameters from the FHiCL file.
    fOpDetWaveformLabel = pset.get<std::string>("OpDetWaveformLabel","");
    fDynamicRange       = pset.get<float>("DynamicRange"); // maximum ADC value of the detector in bits
    // The original simulation was done with a baseline of 1864 ADC. To be more comparable to the data, we implement channel dependent baselines.
    baselineMap.insert(std::make_pair(1010, 2940.9));
    baselineMap.insert(std::make_pair(1011, 3035.2));
    baselineMap.insert(std::make_pair(1020, 2462.0));
    baselineMap.insert(std::make_pair(1021, 3399.9));
    baselineMap.insert(std::make_pair(1030, 2633.1));
    baselineMap.insert(std::make_pair(1031, 2827.8));
    baselineMap.insert(std::make_pair(1040, 2985.3));
    baselineMap.insert(std::make_pair(1041, 2684.8));
    baselineMap.insert(std::make_pair(1050, 4477.2));
    baselineMap.insert(std::make_pair(1051, 4386.8));
    baselineMap.insert(std::make_pair(1060, 1963.9));
    baselineMap.insert(std::make_pair(1061, 2047.8));
    baselineMap.insert(std::make_pair(1070, 2048.4));
    baselineMap.insert(std::make_pair(1071, 2068.7));
    baselineMap.insert(std::make_pair(1080, 2250.3));
    baselineMap.insert(std::make_pair(1081, 1854.0));

    art::ServiceHandle< art::TFileService > tfs;

    fCountTree = tfs->make<TTree>("CountWaveforms","CountWaveforms");
    fCountTree->Branch("EventID",       &fEventID,      "EventID/I");
    fCountTree->Branch("nwaveformsSat", &fnwaveformsSat, "nwaveformsSat/I");
    fCountTree->Branch("nwaveformsTotal", &fnwaveformsTotal, "nwaveformsTotal/I");
    fCountTree->Branch("binsatVector", &fBinsatVector, "fBinsatVector[nwaveformsSat]/F");
  }

  //-----------------------------------------------------------------------
  // Destructor
  saturation::~saturation()
  {}

  //-----------------------------------------------------------------------
  void saturation::beginJob()
  {}

  //-----------------------------------------------------------------------
  void saturation::analyze(const art::Event& evt)
  {

    art::ServiceHandle< art::TFileService > tfs;

    // Record the event ID
    fEventID = evt.id().event();

    /////////////////////
    // Count waveforms //
    ////////////////////

    fnwaveformsSat = 0;
    fnwaveformsTotal = 0;
    // Stores the previous number of saturated waveforms
    Int_t nwvfSatI=0; 

    auto wfHandle = evt.getHandle< std::vector< raw::OpDetWaveform > >(fOpDetWaveformLabel);
    if (wfHandle) {
      fnwaveformsTotal = wfHandle->size();
      short int max;

      for (auto wf: *wfHandle) {
        // Nominal ADC saturation threshold from the configured dynamic range
        max = pow(2,fDynamicRange);
        if(wf.ChannelNumber()>1000&&wf.ChannelNumber()<2000) {max=max-(baselineMap[wf.ChannelNumber()]-1854);}
        std::cout << "channel " << wf.ChannelNumber() << " max " << max << std::endl;
        fnbinsat=0;
        // Locate the maximum ADC sample in the waveform
        auto peak = max_element(std::begin(wf), std::end(wf));
        // A waveform is considered saturated if its maximum reaches or exceeds
        // the effective saturation threshold.
        if ( *peak >= max) { 
          ++fnwaveformsSat;
          std::cout << *peak << std::endl;
        }
        while(*peak>=max){
          std::cout << wf.ChannelNumber() << " " << fnbinsat << " " << fnwaveformsSat << " " << *peak << std::endl;
          peak--;
          ++fnbinsat;
        }
        peak = max_element(std::begin(wf), std::end(wf));
        peak++;
        // Count consecutive saturated samples after the maximum.
        while(*peak>=max){
          std::cout << fnwaveformsSat << " " << fnbinsat << " " << *peak << std::endl;
          peak++;
          ++fnbinsat;
        }
        if((fnwaveformsSat-nwvfSatI)!=0){ 
          fBinsatVector[fnwaveformsSat-1]=fnbinsat;
          nwvfSatI=fnwaveformsSat;
        }
      } //end for wvfHandle
      fCountTree->Fill();
        memset(fBinsatVector, 0, sizeof(fBinsatVector)); // Sets all bytes to 0
    }//end if wvfHandle

  }//end analyze
} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(saturation)
}

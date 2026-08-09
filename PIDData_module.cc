/*
 * =============================================================================
 * Data Particle Identification
 * =============================================================================
 *
 * Purpose: this ART/LArSoft analyzer extracts beam information together
 * with the PD response in a time window around the largest
 * optical hit in each event.
 *
 * The resulting TTree can be used to study the spatial/light pattern observed
 * by the PDS and to investigate particle identification.
 *
 * For each ART event, the analyzer:
 *   1. Stores run, subrun, and event identifiers.
 *   2. Reads beam momentum, TOF, and Cherenkov statuses.
 *   3. Reads reconstructed optical hits (recob::OpHit).
 *   4. Finds the optical hit with the largest PE value in the event.
 *   5. Uses that hit time as the reference time.
 *   6. Defines an asymmetric integration window around the reference hit:
 *          start = peakTime - 0.04 * Window
 *          end   = peakTime + 0.96 * Window
 *   7. Selects all OpHits inside that interval.
 *   8. Sums PE independently for each optical channel.
 *   9. Stores the channel numbers and corresponding PE sums in pdsTree.
 *
 * Output tree: pdsTree
 *   run, subrun, event : ART identifiers
 *   channel            : optical channels with accepted PE
 *   pesum              : PE sum for each corresponding channel
 *   momentum           : reconstructed beam momentum
 *   tof                : beam time of flight
 *   xcetl              : low-pressure Cherenkov status
 *   xceth              : high-pressure Cherenkov status
 *
 *
 * Adapted by M. Sampaio (maressap@ifi.unicamp.br) and L. Paulucci
 * =============================================================================
 */

#ifndef AnaWidthData_H
#define AnaWidthData_H 1

// ROOT includes
#include "TH1.h"
#include "TEfficiency.h"
#include "TTree.h"

// C++ includes
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include "math.h"
#include <climits>

// LArSoft includes
#include "larcore/Geometry/WireReadout.h"
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "lardataobj/RawData/OpDetWaveform.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "duneopdet/OpticalDetector/OpFlashSort.h"
#include "lardataobj/Simulation/SimChannel.h"
//#include "duneobj/ProtoDUNE/ProtoDUNEBeamEvent.h"
#include "dunecore/DuneObj/ProtoDUNEBeamEvent.h"

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

#include <iomanip>

namespace opdet {

  class AnaPIDData : public art::EDAnalyzer{
  public:

    AnaPIDData(const fhicl::ParameterSet&);
    virtual ~AnaPIDData();
    void beginJob();

    void analyze (const art::Event& evt);


  private:

  std::string fOpHitLabel;
  std::string fBeamModuleLabel;

  double fWindow;        // in ticks
  

    // output tree that will contain the sum of PEs in the desired integrated time (which is defined as a fcl parameter input)
    TTree* fPIDTree;
    int fRun;
    int fSubRun;
    int fEventID;
    double fMomentum;
    double fTOF;
    int    fXCETL;
    int    fXCETH;

    std::vector<int>    fOpChannel;
    std::vector<double> fPEsum;
    
  };
}//end namespace

#endif // AnaPIDData_H

namespace opdet {

  AnaPIDData::AnaPIDData(fhicl::ParameterSet const& pset)
    : EDAnalyzer(pset)
  {

    // Reconstructed hit input label
    fOpHitLabel        = pset.get<std::string>("OpHitLabel");
    // Total integration window width in optical ticks.
    // Default: 256 ticks.
    fWindow      = pset.get<double>("WindowTicks", 256.0); 
    // ProtoDUNE beam-event input label
    fBeamModuleLabel = pset.get<std::string>("BeamModuleLabel");

    art::ServiceHandle< art::TFileService > tfs;


    if (!fOpHitLabel.empty()) {

      fPIDTree = tfs->make<TTree>("pdsTree","PDS Beam Window Tree");
      fPIDTree->Branch("run",    &fRun);
      fPIDTree->Branch("subrun", &fSubRun);
      fPIDTree->Branch("event",  &fEventID);
      fPIDTree->Branch("channel", &fOpChannel);
      fPIDTree->Branch("pesum",   &fPEsum);
      fPIDTree->Branch("momentum", &fMomentum);
      fPIDTree->Branch("tof",      &fTOF);
      fPIDTree->Branch("xcetl",    &fXCETL);
      fPIDTree->Branch("xceth",    &fXCETH);
    }

  }
  
  //-----------------------------------------------------------------------
  // Destructor
  AnaPIDData::~AnaPIDData()
  {}

  //-----------------------------------------------------------------------
  void AnaPIDData::beginJob()
  {}

  //-----------------------------------------------------------------------
  void AnaPIDData::analyze(const art::Event& evt)
  {
    // Get the required services
    art::ServiceHandle< geo::Geometry > geom;
    art::ServiceHandle< art::TFileService > tfs;  
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);
    double tickPeriod = clockData.OpticalClock().TickPeriod();

    // Record the event ID
    fRun    = evt.run();
    fSubRun = evt.subRun();
    fEventID  = evt.id().event();
    

    fOpChannel.clear();
    fPEsum.clear();
    
    // reset
    fMomentum = -999.;
    fTOF      = -999.;
    fXCETL    = -1;
    fXCETH    = -1;

auto beamHandle = evt.getHandle<std::vector<beam::ProtoDUNEBeamEvent>>(fBeamModuleLabel);
std::cout << "Beam label = "<< fBeamModuleLabel<< std::endl;
if (!beamHandle) {
    std::cout << "NO BEAM HANDLE FOUND" << std::endl;
}

if (beamHandle) {

  auto const& beamEvent = beamHandle->front();
  
  fTOF   = beamEvent.GetTOF();
  fXCETH = beamEvent.GetCKov0Status(); // HIGH PRESSURE
  fXCETL = beamEvent.GetCKov1Status(); // LOW PRESSURE

  std::cout << "\n=== Beam info ===\n";

  std::cout << "TOF = " << beamEvent.GetTOF() << std::endl;
  std::cout << "TOFChan = " << beamEvent.GetTOFChan() << std::endl;

  std::cout << "NTOFs = "
            << beamEvent.GetTOFs().size()
            << std::endl;

  for(size_t i=0;i<beamEvent.GetTOFs().size();++i){
    std::cout
      << "TOF[" << i << "] = "
      << beamEvent.GetTOFs()[i]
      << "  Chan="
      << beamEvent.GetTOFChans()[i]
      << std::endl;
  }

  std::cout << "CKov0 = "
            << beamEvent.GetCKov0Status()
            << std::endl;

  std::cout << "CKov1 = "
            << beamEvent.GetCKov1Status()
            << std::endl;

  std::cout << "NRecoMom = "
            << beamEvent.GetRecoBeamMomenta().size()
            << std::endl;

  for(size_t i=0;
      i<beamEvent.GetRecoBeamMomenta().size();
      ++i){

    std::cout
      << "P[" << i << "] = "
      << beamEvent.GetRecoBeamMomentum(i)
      << std::endl;
    
    fMomentum= beamEvent.GetRecoBeamMomentum(i);
  }
}
       
    if (!fOpHitLabel.empty()) {

    auto ophits = evt.getHandle<std::vector<recob::OpHit>>(fOpHitLabel);

    if (ophits) {

        // Find global maximum PE hit in the event
        double maxPE       = -1.0;
        double peakTimeUS  = -9999.0;
        int    peakChannel = -1;
        

        for (auto const& hit : *ophits) {

            double hitTimeUS = hit.PeakTime() * tickPeriod;   // convert ticks to microseconds
            double pe        = hit.PE();

            if (pe > maxPE) {
                maxPE       = pe;
                peakTimeUS  = hitTimeUS;
                peakChannel = hit.OpChannel();
            }
        }

        double windowUS = fWindow * tickPeriod;
        double startTime = peakTimeUS - (windowUS*0.04);
        double endTime   = peakTimeUS + (windowUS*0.96);

        
        std::map<int,double> peMap;

        int acceptedHits = 0;

        for (auto const& hit : *ophits) {

            int ch = hit.OpChannel();

            double hitTimeUS = hit.PeakTime() * tickPeriod;
            double pe        = hit.PE();

            if (hitTimeUS < startTime) continue;
            if (hitTimeUS >= endTime) continue;

            peMap[ch] += pe;
            acceptedHits++;
        }

        for (auto const& kv : peMap) {
            fOpChannel.push_back(kv.first);
            fPEsum.push_back(kv.second);
        }

        std::cout << std::fixed << std::setprecision(6);

        std::cout << "\n=====================================================\n";
        std::cout << "Peak-search timing mode\n";
        std::cout << "-----------------------------------------------------\n";
        std::cout << "Run/SubRun/Event   : "
                  << fRun << " / "
                  << fSubRun << " / "
                  << fEventID << "\n";

        std::cout << "Total OpHits       : "
                  << ophits->size() << "\n";

        std::cout << "Maximum PE hit\n";
        std::cout << "   Channel         : " << peakChannel << "\n";
        std::cout << "   Peak PE         : " << maxPE << "\n";
        std::cout << "   Peak time       : " << peakTimeUS << " us\n";

        std::cout << "Integration window : ["
                  << startTime << " , "
                  << endTime   << "] us\n";

        std::cout << "Window width       : "
                  << (endTime - startTime)
                  << " us\n";

        std::cout << "Accepted hits      : "
                  << acceptedHits << "\n";

        std::cout << "Channels with PE   : "
                  << peMap.size() << "\n";

        std::cout << "-----------------------------------------------------\n";
        std::cout << "Channel sums:\n";

        for (auto const& kv : peMap) {
            std::cout << "   Ch "
                      << kv.first
                      << "  --> "
                      << kv.second
                      << " PE\n";
        }

        std::cout << "=====================================================\n";
    }

    fPIDTree->Fill();
}
    
    
  }//void

 

} // namespace opdet

namespace opdet {
  DEFINE_ART_MODULE(AnaPIDData)
}

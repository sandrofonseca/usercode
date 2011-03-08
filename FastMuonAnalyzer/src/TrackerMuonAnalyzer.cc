

#include "FastMuonBremAnalyzers/FastMuonAnalyzer/interface/TrackerMuonAnalyzer.h"

// Collaborating Class Header
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
//Tracker Geometry
// #include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"
//

#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

// Muon Geometry
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "Geometry/DTGeometry/interface/DTTopology.h"

// Digis
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"

// Local and global quantities
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include <Math/VectorUtil.h>
#include <Math/Point3D.h>

// namespaces
using namespace std;
using namespace edm;

// constants, enums and typedefs
// typedef std::vector<L1MuGMTCand> L1MuonCollection;
typedef std::vector<l1extra::L1MuonParticle> L1MuonCollection;


/// Constructor
TrackerMuonAnalyzer::TrackerMuonAnalyzer(const ParameterSet& pset){

  theRootFileName = pset.getUntrackedParameter<string>("rootFileName");
  theSimTrackLabel = pset.getParameter<edm::InputTag>("labelSimTrack");
  theSimVertexLabel = pset.getParameter<edm::InputTag>("labelSimVertex");

  theDTHits = pset.getParameter<edm::InputTag>("dtSimHits");
  theCSCHits = pset.getParameter<edm::InputTag>("cscSimHits");
  theRPCHits = pset.getParameter<edm::InputTag>("rpcSimHits");
 
  theRecoTrackLabelFast = pset.getParameter<edm::InputTag>("labelRecoTrackFast");
  theL1MuonLabelFast = pset.getParameter<edm::InputTag>("labelL1MuonFast");
  theL2MuonLabelFast = pset.getParameter<edm::InputTag>("labelL2MuonFast");
  theL3MuonLabelFast = pset.getParameter<edm::InputTag>("labelL3MuonFast");
  theGlobalMuonLabelFast = pset.getParameter<edm::InputTag>("labelGlobalMuonFast");
  theGlobalMuonTrackLabelFast = pset.getParameter<edm::InputTag>("labelGlobalMuonTrackFast");
  theStandAloneMuonTrackLabelFast = pset.getParameter<edm::InputTag>("labelStandAloneMuonTrackFast");

  // //Tracker
//   theSimHitTOBLabel = pset.getParameter<edm::InputTag>("labelSimHitsTOB");
//   theSimHitTIBLabel = pset.getParameter<edm::InputTag>("labelSimHitsTIB");
  //
  isFastSim = pset.getUntrackedParameter<bool>("FastSim");
  ignoreMissingCollections = pset.getUntrackedParameter<bool>("IgnoreMissingCollections");
  theRootTree = pset.getUntrackedParameter<bool>("RootTree");
  theDebug = pset.getUntrackedParameter<bool>("Debug");

  //  SegmentationBarrel = 2;
  SegmentationBarrel = 12;
  //  SegmentationEndcap = 2;
  SegmentationEndcap = 36;

  readMuonHits = theRootTree && (theDTHits.label() != "MISSING");
}

/// Destructor
TrackerMuonAnalyzer::~TrackerMuonAnalyzer(){
}

void TrackerMuonAnalyzer::beginJob(){

  // Initialize counters
  numberOfEvents=0;
  theFirstEvent = true;

  // Create the root file
  theFile = new TFile(theRootFileName.c_str(), "RECREATE");
  theFile->cd();


  // Parametri dei plots
  Int_t NBinEta = 120;
  Int_t NBinPhi = 100;
  Int_t NBinRes = 100;
  Float_t etamax = 2.4;
  Float_t phibarmax = 360./SegmentationBarrel;
  Float_t phiendmax = 360./SegmentationEndcap;

  const Int_t npt = 137;
  Float_t vecpt[npt] = {
          0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,
         10.,  11.,  12.,  13.,  14.,  15.,  16.,  17.,  18.,  19.,
         20.,  21.,  22.,  23.,  24.,  25.,  26.,  27.,  28.,  29.,
         30.,  31.,  32.,  33.,  34.,  35.,  36.,  37.,  38.,  39.,
         40.,  41.,  42.,  43.,  44.,  45.,  46.,  47.,  48.,  49.,
         50.,  51.,  52.,  53.,  54.,  55.,  56.,  57.,  58.,  59.,
         60.,  61.,  62.,  63.,  64.,  65.,  66.,  67.,  68.,  69.,
         70.,  71.,  72.,  73.,  74.,  75.,  76.,  77.,  78.,  79.,
         80.,  81.,  82.,  83.,  84.,  85.,  86.,  87.,  88.,  89.,
         90.,  91.,  92.,  93.,  94.,  95.,  96.,  97.,  98.,  99.,
        100., 110., 120., 130., 140., 150., 160., 170., 180., 190.,
        200., 210., 220., 230., 240., 250., 260., 270., 280., 290.,
        300., 310., 320., 330., 340., 350., 360., 370., 380., 390.,
        400., 500., 600., 700., 800., 900., 1000. };

  const Int_t nptl1 = 32;
  Float_t vecptl1[nptl1] = {
             0.0,   1.5,   2.0,   2.5,   3.0,   3.5,   4.0,   4.5,
             5.0,   6.0,   7.0,   8.0,  10.0,  12.0,  14.0,  16.0,
            18.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,  50.0,
            60.0,  70.0,  80.0,  90.0, 100.0, 120.0, 140.0, 160.0  };

  //
  // Histogrammi:
  //
  if(test_Track)
 {
 trackerPTIn= new TH1F("trackerPtIn","trackerIn p_{T}",100,999.90,1000.10);
 trackerPTOn= new TH1F("trackerPtOn","trackerOut p_{T}",100,999.90,1000.10);
 }


  // Generated muons
   if (useGen) {
  GenPtHisto = new TH1F("GenPtHisto","Gen. p_{T}",npt-1,vecpt);
  GenPtDetHisto[0] = new TH1F("GenPtBarrelHisto","Gen. p_{T}",npt-1,vecpt);
  GenPtDetHisto[1] = new TH1F("GenPtOverlapHisto","Gen. p_{T}",npt-1,vecpt);
  GenPtDetHisto[2] = new TH1F("GenPtEndcapHisto","Gen. p_{T}",npt-1,vecpt);
  GenPtDetHisto[3] = new TH1F("GenPtExternalHisto","Gen. p_{T}",npt-1,vecpt);
  GenEtaHisto = new TH1F("GenEtaHisto","Gen. #eta",NBinEta,-etamax,etamax);
  GenAbsEtaHisto = new TH1F("GenAbsEtaHisto","Gen. #eta",NBinEta,0,etamax);
  GenPhiHisto[0] = new TH1F("GenPhiBarrelHisto","Gen. #phi",NBinPhi,0,phibarmax);
  GenPhiHisto[1] = new TH1F("GenPhiOverlapHisto","Gen. #phi",NBinPhi,0,phibarmax);
  GenPhiHisto[2] = new TH1F("GenPhiEndcapHisto","Gen. #phi",NBinPhi,0,phiendmax);
  GenPhiHisto[3] = new TH1F("GenPhiExternalHisto","Gen. #phi",NBinPhi,0,phiendmax);
   }

  // Reco Muons (either Fast or Full)
if (useReco) {
  Lv1PtHisto[0] = new TH1F("Lv1PtBarrelHisto","L1 p_{T}",npt-1,vecpt);
  Lv1PtHisto[1] = new TH1F("Lv1PtOverlapHisto","L1 p_{T}",npt-1,vecpt);
  Lv1PtHisto[2] = new TH1F("Lv1PtEndcapHisto","L1 p_{T}",npt-1,vecpt);
  Lv1PtHisto[3] = new TH1F("Lv1PtExternalHisto","L1 p_{T}",npt-1,vecpt);
  Lv1ChaPtHisto[0] = new TH1F("Lv1ChaPtBarrelHisto","L1 wrong charge",npt-1,vecpt);
  Lv1ChaPtHisto[1] = new TH1F("Lv1ChaPtOverlapHisto","L1 wrong charge",npt-1,vecpt);
  Lv1ChaPtHisto[2] = new TH1F("Lv1ChaPtEndcapHisto","L1 wrong charge",npt-1,vecpt);
  Lv1ChaPtHisto[3] = new TH1F("Lv1ChaPtExternalHisto","L1 wrong charge",npt-1,vecpt);
  Lv1EtaHisto = new TH1F("Lv1EtaHisto","L1 #eta",NBinEta,-etamax,etamax);
  Lv1AbsEtaHisto = new TH1F("Lv1AbsEtaHisto","L1 #eta",NBinEta,0,etamax);
  Lv1ChaEtaHisto = new TH1F("Lv1ChaEtaHisto","L1 #eta",NBinEta,-etamax,etamax);
  Lv1ChaAbsEtaHisto = new TH1F("Lv1ChaAbsEtaHisto","L1 #eta",NBinEta,0,etamax);
  Lv1PhiHisto[0] = new TH1F("Lv1PhiBarrelHisto","L1 #phi",NBinPhi,0,phibarmax);
  Lv1PhiHisto[1] = new TH1F("Lv1PhiOverlapHisto","L1 #phi",NBinPhi,0,phibarmax);
  Lv1PhiHisto[2] = new TH1F("Lv1PhiEndcapHisto","L1 #phi",NBinPhi,0,phiendmax);
  Lv1PhiHisto[3] = new TH1F("Lv1PhiExternalHisto","L1 #phi",NBinPhi,0,phiendmax);
  Lv1PtOutHisto[0] = new TH1F("Lv1PtOutBarrelHisto","L1 p_{T}^{out}",nptl1-1,vecptl1);
  Lv1PtOutHisto[1] = new TH1F("Lv1PtOutOverlapHisto","L1 p_{T}^{out}",nptl1-1,vecptl1);
  Lv1PtOutHisto[2] = new TH1F("Lv1PtOutEndcapHisto","L1 p_{T}^{out}",nptl1-1,vecptl1);
  Lv1PtOutHisto[3] = new TH1F("Lv1PtOutExternalHisto","L1 p_{T}^{out}",nptl1-1,vecptl1);

  Lv2PtHisto = new TH1F("Lv2PtHisto","L2 p_{T}",npt-1,vecpt);
  Lv2ChaPtHisto = new TH1F("Lv2ChaPtHisto","L2 wrong charge",npt-1,vecpt);
  Lv2EtaHisto = new TH1F("Lv2EtaHisto","L2 #eta",NBinEta,-etamax,etamax);
  Lv2AbsEtaHisto = new TH1F("Lv2AbsEtaHisto","L2 #eta",NBinEta,0,etamax);
  Lv2ChaEtaHisto = new TH1F("Lv2ChaEtaHisto","L2 #eta",NBinEta,-etamax,etamax);
  Lv2ChaAbsEtaHisto = new TH1F("Lv2ChaAbsEtaHisto","L2 #eta",NBinEta,0,etamax);
  Lv2PhiHisto[0] = new TH1F("Lv2PhiBarrelHisto","L2 #phi",NBinPhi,0,phibarmax);
  Lv2PhiHisto[1] = new TH1F("Lv2PhiOverlapHisto","L2 #phi",NBinPhi,0,phibarmax);
  Lv2PhiHisto[2] = new TH1F("Lv2PhiEndcapHisto","L2 #phi",NBinPhi,0,phiendmax);
  Lv2PhiHisto[3] = new TH1F("Lv2PhiExternalHisto","L2 #phi",NBinPhi,0,phiendmax);
  Lv2PtResTotHisto = new TH1F("Lv2PtResTotHisto","p_{T}(L2) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv2PtResHisto[0] = new TH1F("Lv2PtResBarrelHisto","p_{T}(L2) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv2PtResHisto[1] = new TH1F("Lv2PtResOverlapHisto","p_{T}(L2) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv2PtResHisto[2] = new TH1F("Lv2PtResEndcapHisto","p_{T}(L2) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv2PtResHisto[3] = new TH1F("Lv2PtResExternalHisto","p_{T}(L2) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv2InvPtResTotHisto = new TH1F("Lv2InvPtResTotHisto","1/p_{T}(L2) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv2InvPtResHisto[0] = new TH1F("Lv2InvPtResBarrelHisto","1/p_{T}(L2) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv2InvPtResHisto[1] = new TH1F("Lv2InvPtResOverlapHisto","1/p_{T}(L2) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv2InvPtResHisto[2] = new TH1F("Lv2InvPtResEndcapHisto","1/p_{T}(L2) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv2InvPtResHisto[3] = new TH1F("Lv2InvPtResExternalHisto","1/p_{T}(L2) / 1/p_{T}(gen)",NBinRes,0.50,1.50);

  Lv3PtHisto = new TH1F("Lv3PtHisto","L3 p_{T}",npt-1,vecpt);
  Lv3ChaPtHisto = new TH1F("Lv3ChaPtHisto","L3 wrong charge",npt-1,vecpt);
  Lv3EtaHisto = new TH1F("Lv3EtaHisto","L3 #eta",NBinEta,-etamax,etamax);
  Lv3AbsEtaHisto = new TH1F("Lv3AbsEtaHisto","L3 #eta",NBinEta,0,etamax);
  Lv3ChaEtaHisto = new TH1F("Lv3ChaEtaHisto","L3 #eta",NBinEta,-etamax,etamax);
  Lv3ChaAbsEtaHisto = new TH1F("Lv3ChaAbsEtaHisto","L3 #eta",NBinEta,0,etamax);
  Lv3PhiHisto[0] = new TH1F("Lv3PhiBarrelHisto","L3 #phi",NBinPhi,0,phibarmax);
  Lv3PhiHisto[1] = new TH1F("Lv3PhiOverlapHisto","L3 #phi",NBinPhi,0,phibarmax);
  Lv3PhiHisto[2] = new TH1F("Lv3PhiEndcapHisto","L3 #phi",NBinPhi,0,phiendmax);
  Lv3PhiHisto[3] = new TH1F("Lv3PhiExternalHisto","L3 #phi",NBinPhi,0,phiendmax);
  Lv3PtResTotHisto = new TH1F("Lv3PtResTotHisto","p_{T}(L3) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv3PtResHisto[0] = new TH1F("Lv3PtResBarrelHisto","p_{T}(L3) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv3PtResHisto[1] = new TH1F("Lv3PtResOverlapHisto","p_{T}(L3) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv3PtResHisto[2] = new TH1F("Lv3PtResEndcapHisto","p_{T}(L3) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv3PtResHisto[3] = new TH1F("Lv3PtResExternalHisto","p_{T}(L3) / p_{T}(gen)",2*NBinRes,0.,3.);
  Lv3InvPtResTotHisto = new TH1F("Lv3InvPtResTotHisto","1/p_{T}(L3) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv3InvPtResHisto[0] = new TH1F("Lv3InvPtResBarrelHisto","1/p_{T}(L3) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv3InvPtResHisto[1] = new TH1F("Lv3InvPtResOverlapHisto","1/p_{T}(L3) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv3InvPtResHisto[2] = new TH1F("Lv3InvPtResEndcapHisto","1/p_{T}(L3) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  Lv3InvPtResHisto[3] = new TH1F("Lv3InvPtResExternalHisto","1/p_{T}(L3) / 1/p_{T}(gen)",NBinRes,0.50,1.50);

  GLPtHisto = new TH1F("GLPtHisto","GL p_{T}",npt-1,vecpt);
  GLChaPtHisto = new TH1F("GLChaPtHisto","GL wrong charge",npt-1,vecpt);
  GLEtaHisto = new TH1F("GLEtaHisto","GL #eta",NBinEta,-etamax,etamax);
  GLAbsEtaHisto = new TH1F("GLAbsEtaHisto","GL #eta",NBinEta,0,etamax);
  GLChaEtaHisto = new TH1F("GLChaEtaHisto","GL #eta",NBinEta,-etamax,etamax);
  GLChaAbsEtaHisto = new TH1F("GLChaAbsEtaHisto","GL #eta",NBinEta,0,etamax);
  GLPhiHisto[0] = new TH1F("GLPhiBarrelHisto","GL #phi",NBinPhi,0,phibarmax);
  GLPhiHisto[1] = new TH1F("GLPhiOverlapHisto","GL #phi",NBinPhi,0,phibarmax);
  GLPhiHisto[2] = new TH1F("GLPhiEndcapHisto","GL #phi",NBinPhi,0,phiendmax);
  GLPhiHisto[3] = new TH1F("GLPhiExternalHisto","GL #phi",NBinPhi,0,phiendmax);
  GLPtResTotHisto = new TH1F("GLPtResTotHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLPtResHisto[0] = new TH1F("GLPtResBarrelHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLPtResHisto[1] = new TH1F("GLPtResOverlapHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLPtResHisto[2] = new TH1F("GLPtResEndcapHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLPtResHisto[3] = new TH1F("GLPtResExternalHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLInvPtResTotHisto = new TH1F("GLInvPtResTotHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLInvPtResHisto[0] = new TH1F("GLInvPtResBarrelHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLInvPtResHisto[1] = new TH1F("GLInvPtResOverlapHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLInvPtResHisto[2] = new TH1F("GLInvPtResEndcapHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLInvPtResHisto[3] = new TH1F("GLInvPtResExternalHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);

  GLTrackPtHisto = new TH1F("GLTrackPtHisto","GL p_{T}",npt-1,vecpt);
  GLTrackChaPtHisto = new TH1F("GLTrackChaPtHisto","GL wrong charge",npt-1,vecpt);
  GLTrackEtaHisto = new TH1F("GLTrackEtaHisto","GL #eta",NBinEta,-etamax,etamax);
  GLTrackAbsEtaHisto = new TH1F("GLTrackAbsEtaHisto","GL #eta",NBinEta,0,etamax);
  GLTrackChaEtaHisto = new TH1F("GLTrackChaEtaHisto","GL #eta",NBinEta,-etamax,etamax);
  GLTrackChaAbsEtaHisto = new TH1F("GLTrackChaAbsEtaHisto","GL #eta",NBinEta,0,etamax);
  GLTrackPhiHisto[0] = new TH1F("GLTrackPhiBarrelHisto","GL #phi",NBinPhi,0,phibarmax);
  GLTrackPhiHisto[1] = new TH1F("GLTrackPhiOverlapHisto","GL #phi",NBinPhi,0,phibarmax);
  GLTrackPhiHisto[2] = new TH1F("GLTrackPhiEndcapHisto","GL #phi",NBinPhi,0,phiendmax);
  GLTrackPhiHisto[3] = new TH1F("GLTrackPhiExternalHisto","GL #phi",NBinPhi,0,phiendmax);
  GLTrackPtResTotHisto = new TH1F("GLTrackPtResTotHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLTrackPtResHisto[0] = new TH1F("GLTrackPtResBarrelHisto","p_{T}(GL) / 1/p_{T}(gen)",2*NBinRes,0.,3.);
  GLTrackPtResHisto[1] = new TH1F("GLTrackPtResOverlapHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLTrackPtResHisto[2] = new TH1F("GLTrackPtResEndcapHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLTrackPtResHisto[3] = new TH1F("GLTrackPtResExternalHisto","p_{T}(GL) / p_{T}(gen)",2*NBinRes,0.,3.);
  GLTrackInvPtResTotHisto = new TH1F("GLTrackInvPtResTotHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLTrackInvPtResHisto[0] = new TH1F("GLTrackInvPtResBarrelHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLTrackInvPtResHisto[1] = new TH1F("GLTrackInvPtResOverlapHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLTrackInvPtResHisto[2] = new TH1F("GLTrackInvPtResEndcapHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  GLTrackInvPtResHisto[3] = new TH1F("GLTrackInvPtResExternalHisto","1/p_{T}(GL) / 1/p_{T}(gen)",NBinRes,0.50,1.50);

  STATrackPtHisto = new TH1F("STATrackPtHisto","STA p_{T}",npt-1,vecpt);
  STATrackChaPtHisto = new TH1F("STATrackChaPtHisto","STA wrong charge",npt-1,vecpt);
  STATrackEtaHisto = new TH1F("STATrackEtaHisto","STA #eta",NBinEta,-etamax,etamax);
  STATrackAbsEtaHisto = new TH1F("STATrackAbsEtaHisto","STA #eta",NBinEta,0,etamax);
  STATrackChaEtaHisto = new TH1F("STATrackChaEtaHisto","STA #eta",NBinEta,-etamax,etamax);
  STATrackChaAbsEtaHisto = new TH1F("STATrackChaAbsEtaHisto","STA #eta",NBinEta,0,etamax);
  STATrackPhiHisto[0] = new TH1F("STATrackPhiBarrelHisto","STA #phi",NBinPhi,0,phibarmax);
  STATrackPhiHisto[1] = new TH1F("STATrackPhiOverlapHisto","STA #phi",NBinPhi,0,phibarmax);
  STATrackPhiHisto[2] = new TH1F("STATrackPhiEndcapHisto","STA #phi",NBinPhi,0,phiendmax);
  STATrackPhiHisto[3] = new TH1F("STATrackPhiExternalHisto","STA #phi",NBinPhi,0,phiendmax);
  STATrackPtResTotHisto = new TH1F("STATrackPtResTotHisto","p_{T}(STA)/p_{T}(gen)",2*NBinRes,0.,3.);
  STATrackPtResHisto[0] = new TH1F("STATrackPtResBarrelHisto","p_{T}(STA)/p_{T}(gen)",2*NBinRes,0.,3.);
  STATrackPtResHisto[1] = new TH1F("STATrackPtResOverlapHisto","p_{T}(STA)/p_{T}(gen)",2*NBinRes,0.,3.);
  STATrackPtResHisto[2] = new TH1F("STATrackPtResEndcapHisto","p_{T}(STA)/p_{T}(gen)",2*NBinRes,0.,3.);
  STATrackPtResHisto[3] = new TH1F("STATrackPtResExternalHisto","p_{T}(STA)/p_{T}(gen)",2*NBinRes,0.,3.);
  STATrackInvPtResTotHisto = new TH1F("STATrackInvPtResTotHisto","1/p_{T}(STA) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  STATrackInvPtResHisto[0] = new TH1F("STATrackInvPtResBarrelHisto","1/p_{T}(STA) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  STATrackInvPtResHisto[1] = new TH1F("STATrackInvPtResOverlapHisto","1/p_{T}(STA) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  STATrackInvPtResHisto[2] = new TH1F("STATrackInvPtResEndcapHisto","1/p_{T}(STA) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  STATrackInvPtResHisto[3] = new TH1F("STATrackInvPtResExternalHisto","1/p_{T}(STA) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  STATrackChi2TotHisto = new TH1F("STATrackChi2TotHisto","#chi^2",200,0.,200.);
  STATrackNdofTotHisto = new TH1F("STATrackNdofTotHisto","N_{dof}",160,0.,80.);
  STATrackNormChi2TotHisto = new TH1F("STATrackNormChi2TotHisto","#chi^2",200,0.,20.);
  STATrackChi2Histo[0] = new TH1F("STATrackChi2BarrelHisto","#chi^2",200,0.,200.);
  STATrackChi2Histo[1] = new TH1F("STATrackChi2OverlapHisto","#chi^2",200,0.,200.);
  STATrackChi2Histo[2] = new TH1F("STATrackChi2EndcapHisto","#chi^2",200,0.,200.);
  STATrackChi2Histo[3] = new TH1F("STATrackChi2ExternalHisto","#chi^2",200,0.,200.);
  STATrackNdofHisto[0] = new TH1F("STATrackNdofBarrelHisto","N_{dof}",160,0.,80.);
  STATrackNdofHisto[1] = new TH1F("STATrackNdofOverlapHisto","N_{dof}",160,0.,80.);
  STATrackNdofHisto[2] = new TH1F("STATrackNdofEndcapHisto","N_{dof}",200,0.,80.);
  STATrackNdofHisto[3] = new TH1F("STATrackNdofExternalHisto","N_{dof}",160,0.,80.);
  STATrackNormChi2Histo[0] = new TH1F("STATrackNormChi2BarrelHisto","#chi^2",200,0.,20.);
  STATrackNormChi2Histo[1] = new TH1F("STATrackNormChi2OverlapHisto","#chi^2",200,0.,20.);
  STATrackNormChi2Histo[2] = new TH1F("STATrackNormChi2EndcapHisto","#chi^2",200,0.,20.);
  STATrackNormChi2Histo[3] = new TH1F("STATrackNormChi2ExternalHisto","#chi^2",200,0.,20.);

  TKPtHisto = new TH1F("TKPtHisto","TK p_{T}",npt-1,vecpt);
  TKChaPtHisto = new TH1F("TKChaPtHisto","TK wrong charge",npt-1,vecpt);
  TKEtaHisto = new TH1F("TKEtaHisto","TK #eta",NBinEta,-etamax,etamax);
  TKAbsEtaHisto = new TH1F("TKAbsEtaHisto","TK #eta",NBinEta,0,etamax);
  TKChaEtaHisto = new TH1F("TKChaEtaHisto","TK #eta",NBinEta,-etamax,etamax);
  TKChaAbsEtaHisto = new TH1F("TKChaAbsEtaHisto","TK #eta",NBinEta,0,etamax);
  TKPhiHisto[0] = new TH1F("TKPhiBarrelHisto","STA #phi",NBinPhi,0,phibarmax);
  TKPhiHisto[1] = new TH1F("TKPhiOverlapHisto","STA #phi",NBinPhi,0,phibarmax);
  TKPhiHisto[2] = new TH1F("TKPhiEndcapHisto","STA #phi",NBinPhi,0,phiendmax);
  TKPhiHisto[3] = new TH1F("TKPhiExternalHisto","STA #phi",NBinPhi,0,phiendmax);
  TKPtResTotHisto = new TH1F("TKPtResTotHisto","p_{T}(TK)/p_{T}(gen)",2*NBinRes,0.,3.);
  TKPtResHisto[0] = new TH1F("TKPtResBarrelHisto","p_{T}(TK)/p_{T}(gen)",2*NBinRes,0.,3.);
  TKPtResHisto[1] = new TH1F("TKPtResOverlapHisto","p_{T}(TK)/p_{T}(gen)",2*NBinRes,0.,3.);
  TKPtResHisto[2] = new TH1F("TKPtResEndcapHisto","p_{T}(TK)/p_{T}(gen)",2*NBinRes,0.,3.);
  TKPtResHisto[3] = new TH1F("TKPtResExternalHisto","p_{T}(TK)/p_{T}(gen)",2*NBinRes,0.,3.);
  TKInvPtResTotHisto = new TH1F("TKInvPtResTotHisto","1/p_{T}(TK) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  TKInvPtResHisto[0] = new TH1F("TKInvPtResBarrelHisto","1/p_{T}(TK) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  TKInvPtResHisto[1] = new TH1F("TKInvPtResOverlapHisto","1/p_{T}(TK) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  TKInvPtResHisto[2] = new TH1F("TKInvPtResEndcapHisto","1/p_{T}(TK) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
  TKInvPtResHisto[3] = new TH1F("TKInvPtResExternalHisto","1/p_{T}(TK) / 1/p_{T}(gen)",NBinRes,0.50,1.50);
 }
if (useEff) {
  EffLv1PtHisto[0] = new TH1F("EffLv1PtBarrelHisto","Eff L1 vs p_{T}",npt-1,vecpt);
  EffLv1PtHisto[1] = new TH1F("EffLv1PtOverlapHisto","Eff L1 vs p_{T}",npt-1,vecpt);
  EffLv1PtHisto[2] = new TH1F("EffLv1PtEndcapHisto","Eff L1 vs p_{T}",npt-1,vecpt);
  EffLv1PtHisto[3] = new TH1F("EffLv1PtExternalHisto","Eff L1 vs p_{T}",npt-1,vecpt);
  FracWrgChaLv1PtHisto[0] = new TH1F("FracWrgChaLv1PtBarrelHisto","Fract L1 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  FracWrgChaLv1PtHisto[1] = new TH1F("FracWrgChaLv1PtOverlapHisto","Fract L1 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  FracWrgChaLv1PtHisto[2] = new TH1F("FracWrgChaLv1PtEndcapHisto","Fract L1 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  FracWrgChaLv1PtHisto[3] = new TH1F("FracWrgChaLv1PtExternalHisto","Fract L1 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffLv1EtaHisto = new TH1F("EffLv1EtaHisto","Eff L1 vs #eta",NBinEta,-etamax,etamax);
  EffLv1AbsEtaHisto = new TH1F("EffLv1AbsEtaHisto","Eff L1 vs #eta",NBinEta,0,etamax);
  FracWrgChaLv1EtaHisto = new TH1F("FracWrgChaLv1EtaHisto"," Fract L1 with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaLv1AbsEtaHisto = new TH1F("FracWrgChaLv1AbsEtaHisto"," Fract L1 with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffLv1PhiHisto[0] = new TH1F("EffLv1PhiBarrelHisto","Eff L1 vs #phi",NBinPhi,0,phibarmax);
  EffLv1PhiHisto[1] = new TH1F("EffLv1PhiOverlapHisto","Eff L1 vs #phi",NBinPhi,0,phibarmax);
  EffLv1PhiHisto[2] = new TH1F("EffLv1PhiEndcapHisto","Eff L1 vs #phi",NBinPhi,0,phiendmax);
  EffLv1PhiHisto[3] = new TH1F("EffLv1PhiExternalHisto","Eff L1 vs #phi",NBinPhi,0,phiendmax);

  EffLv2PtHisto = new TH1F("EffLv2PtHisto","Eff L2 vs p_{T}",npt-1,vecpt);
  FracWrgChaLv2PtHisto = new TH1F("FracWrgChaLv2PtHisto","Fract L2 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffLv2EtaHisto = new TH1F("EffLv2EtaHisto","Eff L2 vs #eta",NBinEta,-etamax,etamax);
  EffLv2AbsEtaHisto = new TH1F("EffLv2AbsEtaHisto","Eff L2 vs #eta",NBinEta,0,etamax);
  FracWrgChaLv2EtaHisto = new TH1F("FracWrgChaLv2EtaHisto"," Fract L2 with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaLv2AbsEtaHisto = new TH1F("FracWrgChaLv2AbsEtaHisto"," Fract L2 with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffLv2PhiHisto[0] = new TH1F("EffLv2PhiBarrelHisto","Eff L2 vs #phi",NBinPhi,0,phibarmax);
  EffLv2PhiHisto[1] = new TH1F("EffLv2PhiOverlapHisto","Eff L2 vs #phi",NBinPhi,0,phibarmax);
  EffLv2PhiHisto[2] = new TH1F("EffLv2PhiEndcapHisto","Eff L2 vs #phi",NBinPhi,0,phiendmax);
  EffLv2PhiHisto[3] = new TH1F("EffLv2PhiExternalHisto","Eff L2 vs #phi",NBinPhi,0,phiendmax);

  EffLv3PtHisto = new TH1F("EffLv3PtHisto","Eff L3 vs p_{T}",npt-1,vecpt);
  FracWrgChaLv3PtHisto = new TH1F("FracWrgChaLv3PtHisto","Fract L3 with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffLv3EtaHisto = new TH1F("EffLv3EtaHisto","Eff L3 vs #eta",NBinEta,-etamax,etamax);
  EffLv3AbsEtaHisto = new TH1F("EffLv3AbsEtaHisto","Eff L3 vs #eta",NBinEta,0,etamax);
  FracWrgChaLv3EtaHisto = new TH1F("FracWrgChaLv3EtaHisto"," Fract L3 with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaLv3AbsEtaHisto = new TH1F("FracWrgChaLv3AbsEtaHisto"," Fract L3 with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffLv3PhiHisto[0] = new TH1F("EffLv3PhiBarrelHisto","Eff L3 vs #phi",NBinPhi,0,phibarmax);
  EffLv3PhiHisto[1] = new TH1F("EffLv3PhiOverlapHisto","Eff L3 vs #phi",NBinPhi,0,phibarmax);
  EffLv3PhiHisto[2] = new TH1F("EffLv3PhiEndcapHisto","Eff L3 vs #phi",NBinPhi,0,phiendmax);
  EffLv3PhiHisto[3] = new TH1F("EffLv3PhiExternalHisto","Eff L3 vs #phi",NBinPhi,0,phiendmax);

  EffGLPtHisto = new TH1F("EffGLPtHisto","Eff GL vs p_{T}",npt-1,vecpt);
  FracWrgChaGLPtHisto = new TH1F("FracWrgChaGLPtHisto","Fract GL with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffGLEtaHisto = new TH1F("EffGLEtaHisto","Eff GL vs #eta",NBinEta,-etamax,etamax);
  EffGLAbsEtaHisto = new TH1F("EffGLAbsEtaHisto","Eff GL vs #eta",NBinEta,0,etamax);
  FracWrgChaGLEtaHisto = new TH1F("FracWrgChaGLEtaHisto"," Fract GL with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaGLAbsEtaHisto = new TH1F("FracWrgChaGLAbsEtaHisto"," Fract GL with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffGLPhiHisto[0] = new TH1F("EffGLPhiBarrelHisto","Eff GL vs #phi",NBinPhi,0,phibarmax);
  EffGLPhiHisto[1] = new TH1F("EffGLPhiOverlapHisto","Eff GL vs #phi",NBinPhi,0,phibarmax);
  EffGLPhiHisto[2] = new TH1F("EffGLPhiEndcapHisto","Eff GL vs #phi",NBinPhi,0,phiendmax);
  EffGLPhiHisto[3] = new TH1F("EffGLPhiExternalHisto","Eff GL vs #phi",NBinPhi,0,phiendmax);

  EffGLTrackPtHisto = new TH1F("EffGLTrackPtHisto","Eff GL vs p_{T}",npt-1,vecpt);
  FracWrgChaGLTrackPtHisto = new TH1F("FracWrgChaGLTrackPtHisto","Fract GL with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffGLTrackEtaHisto = new TH1F("EffGLTrackEtaHisto","Eff GL vs #eta",NBinEta,-etamax,etamax);
  EffGLTrackAbsEtaHisto = new TH1F("EffGLTrackAbsEtaHisto","Eff GL vs #eta",NBinEta,0,etamax);
  FracWrgChaGLTrackEtaHisto = new TH1F("FracWrgChaGLTrackEtaHisto"," Fract GL with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaGLTrackAbsEtaHisto = new TH1F("FracWrgChaGLTrackAbsEtaHisto"," Fract GL with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffGLTrackPhiHisto[0] = new TH1F("EffGLTrackPhiBarrelHisto","Eff GL vs #phi",NBinPhi,0,phibarmax);
  EffGLTrackPhiHisto[1] = new TH1F("EffGLTrackPhiOverlapHisto","Eff GL vs #phi",NBinPhi,0,phibarmax);
  EffGLTrackPhiHisto[2] = new TH1F("EffGLTrackPhiEndcapHisto","Eff GL vs #phi",NBinPhi,0,phiendmax);
  EffGLTrackPhiHisto[3] = new TH1F("EffGLTrackPhiExternalHisto","Eff GL vs #phi",NBinPhi,0,phiendmax);

  EffSTATrackPtHisto = new TH1F("EffSTATrackPtHisto","Eff STA vs p_{T}",npt-1,vecpt);
  FracWrgChaSTATrackPtHisto = new TH1F("FracWrgChaSTATrackPtHisto","Fract STA with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffSTATrackEtaHisto = new TH1F("EffSTATrackEtaHisto","Eff STA vs #eta",NBinEta,-etamax,etamax);
  EffSTATrackAbsEtaHisto = new TH1F("EffSTATrackAbsEtaHisto","Eff STA vs #eta",NBinEta,0,etamax);
  FracWrgChaSTATrackEtaHisto = new TH1F("FracWrgChaSTATrackEtaHisto"," Fract STA with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaSTATrackAbsEtaHisto = new TH1F("FracWrgChaSTATrackAbsEtaHisto"," Fract STA with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffSTATrackPhiHisto[0] = new TH1F("EffSTATrackPhiBarrelHisto","Eff STA vs #phi",NBinPhi,0,phibarmax);
  EffSTATrackPhiHisto[1] = new TH1F("EffSTATrackPhiOverlapHisto","Eff STA vs #phi",NBinPhi,0,phibarmax);
  EffSTATrackPhiHisto[2] = new TH1F("EffSTATrackPhiEndcapHisto","Eff STA vs #phi",NBinPhi,0,phiendmax);
  EffSTATrackPhiHisto[3] = new TH1F("EffSTATrackPhiExternalHisto","Eff STA vs #phi",NBinPhi,0,phiendmax);

  EffTKPtHisto = new TH1F("EffTKPtHisto","Eff TK vs p_{T}",npt-1,vecpt);
  FracWrgChaTKPtHisto = new TH1F("FracWrgChaTKPtHisto","Fract TK with wrong charge ass. vs p_{T}",npt-1,vecpt);
  EffTKEtaHisto = new TH1F("EffTKEtaHisto","Eff TK vs #eta",NBinEta,-etamax,etamax);
  EffTKAbsEtaHisto = new TH1F("EffTKAbsEtaHisto","Eff TK vs #eta",NBinEta,0,etamax);
  FracWrgChaTKEtaHisto = new TH1F("FracWrgChaTKEtaHisto"," Fract TK with wrong charge ass. vs #eta",NBinEta,-etamax,etamax);
  FracWrgChaTKAbsEtaHisto = new TH1F("FracWrgChaTKAbsEtaHisto"," Fract TK with wrong charge ass. vs #eta",NBinEta,0,etamax);
  EffTKPhiHisto[0] = new TH1F("EffTKPhiBarrelHisto","Eff TK vs #phi",NBinPhi,0,phibarmax);
  EffTKPhiHisto[1] = new TH1F("EffTKPhiOverlapHisto","Eff TK vs #phi",NBinPhi,0,phibarmax);
  EffTKPhiHisto[2] = new TH1F("EffTKPhiEndcapHisto","Eff TK vs #phi",NBinPhi,0,phiendmax);
  EffTKPhiHisto[3] = new TH1F("EffTKPhiExternalHisto","Eff TK vs #phi",NBinPhi,0,phiendmax);
 }

  if (theRootTree) {
  // I root trees

  theSumTree = new TTree("T0","TrackerMuonAnalyzer summary root tree",0);
  theSumTree->Branch("run",&run,"run/I");
  theSumTree->Branch("evt",&evt,"evt/I");
  theSumTree->Branch("nsimu",&nsimu,"nsimu/i");
  theSumTree->Branch("ntkFast",&ntkFast,"ntkFast/s");
  theSumTree->Branch("nl1Fast",&nl1Fast,"nl1Fast/s");
  theSumTree->Branch("nl2Fast",&nl2Fast,"nl2Fast/s");
  theSumTree->Branch("nl3Fast",&nl3Fast,"nl3Fast/s");
  theSumTree->Branch("nglFast",&nglFast,"nglFast/s");
  theSumTree->Branch("ngltkFast",&ngltkFast,"ngltkFast/s");
  theSumTree->Branch("nstatkFast",&nstatkFast,"nstatkFast/s");


  theTree = new TTree("T1","TrackerMuonAnalyzer root tree",0);
  theTree->Branch("run",&run,"run/I");
  theTree->Branch("evt",&evt,"evt/I");
  theTree->Branch("nSimMu",&nSimMu,"nSimMu/i");

  theTree->Branch("pdgSim",&pdgSim,"pdgSim[nSimMu]/I");
  theTree->Branch("pSim",&pSim,"pSim[nSimMu]/F");
  theTree->Branch("pTrackerOut",&pTrackerOut,"pTrackerOut[nSimMu]/F");
  theTree->Branch("pMuonIn",&pMuonIn,"pMuonIn[nSimMu]/F");
  theTree->Branch("pMuonOut",&pMuonOut,"pMuonOut[nSimMu]/F");
  theTree->Branch("ptSim",&ptSim,"ptSim[nSimMu]/F");
  theTree->Branch("ptTrackerOut",&ptTrackerOut,"ptTrackerOut[nSimMu]/F");
  theTree->Branch("ptMuonIn",&ptMuonIn,"ptMuonIn[nSimMu]/F");
  theTree->Branch("ptMuonOut",&ptMuonOut,"ptMuonOut[nSimMu]/F");
  theTree->Branch("phiTrackerOut",&phiTrackerOut,"phiTrackerOut[nSimMu]/F");
  theTree->Branch("phiMuonIn",&phiMuonIn,"phiMuonIn[nSimMu]/F");
  theTree->Branch("phiMuonOut",&phiMuonOut,"phiMuonOut[nSimMu]/F");
  theTree->Branch("phiSim",&phiSim,"phiSim[nSimMu]/F");
  theTree->Branch("etaSim",&etaSim,"etaSim[nSimMu]/F");
  theTree->Branch("chaSim",&chaSim,"chaSim[nSimMu]/I");
  theTree->Branch("pdgMot",&pdgMot,"pdgMot[nSimMu]/I");
  theTree->Branch("ptMot",&ptMot,"ptMot[nSimMu]/F");
  theTree->Branch("phiMot",&phiMot,"phiMot[nSimMu]/F");
  theTree->Branch("etaMot",&etaMot,"etaMot[nSimMu]/F");
  theTree->Branch("nDTHits",&nDTHits,"nDTHits[nSimMu]/i");
  theTree->Branch("nCSCHits",&nCSCHits,"nCSCHits[nSimMu]/i");
  theTree->Branch("nRPCHits",&nRPCHits,"nRPCHits[nSimMu]/i");

  theTree->Branch("ptL1",&ptL1,"ptL1[nSimMu]/F");
  theTree->Branch("phiL1",&phiL1,"phiL1[nSimMu]/F");
  theTree->Branch("etaL1",&etaL1,"etaL1[nSimMu]/F");
  theTree->Branch("chaL1",&chaL1,"chaL1[nSimMu]/I");

  theTree->Branch("ptL2",&ptL2,"ptL2[nSimMu]/F");
  theTree->Branch("phiL2",&phiL2,"phiL2[nSimMu]/F");
  theTree->Branch("etaL2",&etaL2,"etaL2[nSimMu]/F");
  theTree->Branch("chaL2",&chaL2,"chaL2[nSimMu]/I");

  theTree->Branch("ptL3",&ptL3,"ptL3[nSimMu]/F");
  theTree->Branch("phiL3",&phiL3,"phiL3[nSimMu]/F");
  theTree->Branch("etaL3",&etaL3,"etaL3[nSimMu]/F");
  theTree->Branch("chaL3",&chaL3,"chaL3[nSimMu]/I");

  theTree->Branch("ptGL",&ptGL,"ptGL[nSimMu]/F");
  theTree->Branch("phiGL",&phiGL,"phiGL[nSimMu]/F");
  theTree->Branch("etaGL",&etaGL,"etaGL[nSimMu]/F");
  theTree->Branch("chaGL",&chaGL,"chaGL[nSimMu]/I");

  theTree->Branch("ptTK",&ptTK,"ptTK[nSimMu]/F");
  theTree->Branch("phiTK",&phiTK,"phiTK[nSimMu]/F");
  theTree->Branch("etaTK",&etaTK,"etaTK[nSimMu]/F");
  theTree->Branch("chaTK",&chaTK,"chaTK[nSimMu]/I");

  theTree->Branch("ptSTATK",&ptSTATK,"ptSTATK[nSimMu]/F");
  theTree->Branch("phiSTATK",&phiSTATK,"phiSTATK[nSimMu]/F");
  theTree->Branch("etaSTATK",&etaSTATK,"etaSTATK[nSimMu]/F");
  theTree->Branch("chaSTATK",&chaSTATK,"chaSTATK[nSimMu]/I");
  theTree->Branch("chi2STATK",&chi2STATK,"chi2STATK[nSimMu]/F");
  theTree->Branch("ndofSTATK",&ndofSTATK,"ndofSTATK[nSimMu]/F");
  theTree->Branch("dtHitsSTATK",&dtHitsSTATK,"dtHitsSTATK[nSimMu]/I");
  theTree->Branch("cscHitsSTATK",&cscHitsSTATK,"cscHitsSTATK[nSimMu]/I");
  theTree->Branch("rpcHitsSTATK",&rpcHitsSTATK,"rpcHitsSTATK[nSimMu]/I");

  theTree->Branch("ptGLTK",&ptGLTK,"ptGLTK[nSimMu]/F");
  theTree->Branch("phiGLTK",&phiGLTK,"phiGLTK[nSimMu]/F");
  theTree->Branch("etaGLTK",&etaGLTK,"etaGLTK[nSimMu]/F");
  theTree->Branch("chaGLTK",&chaGLTK,"chaGLTK[nSimMu]/I");
  }

}

void TrackerMuonAnalyzer::endJob(){
  cout << endl;
  cout << " --- Total number of events processed : " << numberOfEvents << endl;
    
  // Calcola le efficienze
  
  if (theL1MuonLabelFast.label() != "MISSING") {
    for (int iDet=0;iDet<4;iDet++) {
      EffLv1PtHisto[iDet]->Divide(Lv1PtHisto[iDet],GenPtDetHisto[iDet]);
      FracWrgChaLv1PtHisto[iDet]->Divide(Lv1ChaPtHisto[iDet],Lv1PtHisto[iDet]);
    }
    EffLv1EtaHisto->Divide(Lv1EtaHisto,GenEtaHisto);
    EffLv1AbsEtaHisto->Divide(Lv1AbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaLv1EtaHisto->Divide(Lv1ChaEtaHisto,Lv1EtaHisto);
    FracWrgChaLv1AbsEtaHisto->Divide(Lv1ChaAbsEtaHisto,Lv1AbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffLv1PhiHisto[iDet]->Divide(Lv1PhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theL2MuonLabelFast.label() != "MISSING") {
    EffLv2PtHisto->Divide(Lv2PtHisto,GenPtHisto);
    FracWrgChaLv2PtHisto->Divide(Lv3ChaPtHisto,Lv2PtHisto);
    EffLv2EtaHisto->Divide(Lv2EtaHisto,GenEtaHisto);
    EffLv2AbsEtaHisto->Divide(Lv2AbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaLv2EtaHisto->Divide(Lv2ChaEtaHisto,Lv2EtaHisto);
    FracWrgChaLv2AbsEtaHisto->Divide(Lv2ChaAbsEtaHisto,Lv2AbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffLv2PhiHisto[iDet]->Divide(Lv2PhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theL3MuonLabelFast.label() != "MISSING") {
    EffLv3PtHisto->Divide(Lv3PtHisto,GenPtHisto);
    FracWrgChaLv3PtHisto->Divide(Lv3ChaPtHisto,Lv3PtHisto);
    EffLv3EtaHisto->Divide(Lv3EtaHisto,GenEtaHisto);
    EffLv3AbsEtaHisto->Divide(Lv3AbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaLv3EtaHisto->Divide(Lv3ChaEtaHisto,Lv3EtaHisto);
    FracWrgChaLv3AbsEtaHisto->Divide(Lv3ChaAbsEtaHisto,Lv3AbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffLv3PhiHisto[iDet]->Divide(Lv3PhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theGlobalMuonLabelFast.label() != "MISSING") {
    EffGLPtHisto->Divide(GLPtHisto,GenPtHisto);
    FracWrgChaGLPtHisto->Divide(GLChaPtHisto,GLPtHisto);
    EffGLEtaHisto->Divide(GLEtaHisto,GenEtaHisto);
    EffGLAbsEtaHisto->Divide(GLAbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaGLEtaHisto->Divide(GLChaEtaHisto,GLEtaHisto);
    FracWrgChaGLAbsEtaHisto->Divide(GLChaAbsEtaHisto,GLAbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffGLPhiHisto[iDet]->Divide(GLPhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theGlobalMuonTrackLabelFast.label() != "MISSING") { 
    EffGLTrackPtHisto->Divide(GLTrackPtHisto,GenPtHisto);
    FracWrgChaGLTrackPtHisto->Divide(GLTrackChaPtHisto,GLTrackPtHisto);
    EffGLTrackEtaHisto->Divide(GLTrackEtaHisto,GenEtaHisto);
    EffGLTrackAbsEtaHisto->Divide(GLTrackAbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaGLTrackEtaHisto->Divide(GLTrackChaEtaHisto,GLTrackEtaHisto);
    FracWrgChaGLTrackAbsEtaHisto->Divide(GLTrackChaAbsEtaHisto,GLTrackAbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffGLTrackPhiHisto[iDet]->Divide(GLTrackPhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theStandAloneMuonTrackLabelFast.label() != "MISSING") {
    EffSTATrackPtHisto->Divide(STATrackPtHisto,GenPtHisto);
    FracWrgChaSTATrackPtHisto->Divide(STATrackChaPtHisto,STATrackPtHisto);
    EffSTATrackEtaHisto->Divide(STATrackEtaHisto,GenEtaHisto);
    EffSTATrackAbsEtaHisto->Divide(STATrackAbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaSTATrackEtaHisto->Divide(STATrackChaEtaHisto,STATrackEtaHisto);
    FracWrgChaSTATrackAbsEtaHisto->Divide(STATrackChaAbsEtaHisto,STATrackAbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffSTATrackPhiHisto[iDet]->Divide(STATrackPhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }
  if (theRecoTrackLabelFast.label() != "MISSING") {
    EffTKPtHisto->Divide(TKPtHisto,GenPtHisto);
    FracWrgChaTKPtHisto->Divide(TKChaPtHisto,TKPtHisto);
    EffTKEtaHisto->Divide(TKEtaHisto,GenEtaHisto);
    EffTKAbsEtaHisto->Divide(TKAbsEtaHisto,GenAbsEtaHisto);
    FracWrgChaTKEtaHisto->Divide(TKChaEtaHisto,TKEtaHisto);
    FracWrgChaTKAbsEtaHisto->Divide(TKChaAbsEtaHisto,TKAbsEtaHisto);
    for (int iDet=0;iDet<4;iDet++) {
      EffTKPhiHisto[iDet]->Divide(TKPhiHisto[iDet],GenPhiHisto[iDet]);
    }
  }


  // Write the histos to file
  theFile->cd();

  GenPtHisto->Write();
  
  trackerPTIn->Write();
  trackerPTOn->Write();


  for (int iDet=0;iDet<4;iDet++) GenPtDetHisto[iDet]->Write();
  GenEtaHisto->Write();
  GenAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GenPhiHisto[iDet]->Write();

  for (int iDet=0;iDet<4;iDet++) Lv1PtHisto[iDet]->Write();
  for (int iDet=0;iDet<4;iDet++) Lv1ChaPtHisto[iDet]->Write();
  Lv1EtaHisto->Write();
  Lv1AbsEtaHisto->Write();
  Lv1ChaEtaHisto->Write();
  Lv1ChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv1PhiHisto[iDet]->Write();
  for (int iDet=0;iDet<4;iDet++) Lv1PtOutHisto[iDet]->Write();

  Lv2PtHisto->Write();
  Lv2ChaPtHisto->Write();
  Lv2EtaHisto->Write();
  Lv2AbsEtaHisto->Write();
  Lv2ChaEtaHisto->Write();
  Lv2ChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv2PhiHisto[iDet]->Write();
  Lv2PtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv2PtResHisto[iDet]->Write();
  Lv2InvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv2InvPtResHisto[iDet]->Write();

  Lv3PtHisto->Write();
  Lv3ChaPtHisto->Write();
  Lv3EtaHisto->Write();
  Lv3AbsEtaHisto->Write();
  Lv3ChaEtaHisto->Write();
  Lv3ChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv3PhiHisto[iDet]->Write();
  Lv3PtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv3PtResHisto[iDet]->Write();
  Lv3InvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) Lv3InvPtResHisto[iDet]->Write();

  GLPtHisto->Write();
  GLChaPtHisto->Write();
  GLEtaHisto->Write();
  GLAbsEtaHisto->Write();
  GLChaEtaHisto->Write();
  GLChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLPhiHisto[iDet]->Write();
  for (int iDet=0;iDet<4;iDet++) GLPhiHisto[iDet]->Write();
  GLPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLPtResHisto[iDet]->Write();
  GLInvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLInvPtResHisto[iDet]->Write();

  GLTrackPtHisto->Write();
  GLTrackChaPtHisto->Write();
  GLTrackEtaHisto->Write();
  GLTrackAbsEtaHisto->Write();
  GLTrackChaEtaHisto->Write();
  GLTrackChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLTrackPhiHisto[iDet]->Write();
  GLTrackPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLTrackPtResHisto[iDet]->Write();
  GLTrackInvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) GLTrackInvPtResHisto[iDet]->Write();

  STATrackPtHisto->Write();
  STATrackChaPtHisto->Write();
  STATrackEtaHisto->Write();
  STATrackAbsEtaHisto->Write();
  STATrackChaEtaHisto->Write();
  STATrackChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) STATrackPhiHisto[iDet]->Write();
  STATrackPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) STATrackPtResHisto[iDet]->Write();
  STATrackInvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) STATrackInvPtResHisto[iDet]->Write();
  STATrackChi2TotHisto->Write();
  STATrackNdofTotHisto->Write();
  STATrackNormChi2TotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) {
    STATrackChi2Histo[iDet]->Write();
    STATrackNdofHisto[iDet]->Write();
    STATrackNormChi2Histo[iDet]->Write();
  }

  TKPtHisto->Write();
  TKChaPtHisto->Write();
  TKEtaHisto->Write();
  TKAbsEtaHisto->Write();
  TKChaEtaHisto->Write();
  TKChaAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) TKPhiHisto[iDet]->Write();
  TKPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) TKPtResHisto[iDet]->Write();
  TKInvPtResTotHisto->Write();
  for (int iDet=0;iDet<4;iDet++) TKInvPtResHisto[iDet]->Write();

  for (int iDet=0;iDet<4;iDet++) EffLv1PtHisto[iDet]->Write();
  for (int iDet=0;iDet<4;iDet++) FracWrgChaLv1PtHisto[iDet]->Write();
  EffLv1EtaHisto->Write();
  EffLv1AbsEtaHisto->Write();
  FracWrgChaLv1EtaHisto->Write();
  FracWrgChaLv1AbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffLv1PhiHisto[iDet]->Write();

  EffLv2PtHisto->Write();
  FracWrgChaLv2PtHisto->Write();
  EffLv2EtaHisto->Write();
  EffLv2AbsEtaHisto->Write();
  FracWrgChaLv2EtaHisto->Write();
  FracWrgChaLv2AbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffLv2PhiHisto[iDet]->Write();

  EffLv3PtHisto->Write();
  FracWrgChaLv3PtHisto->Write();
  EffLv3EtaHisto->Write();
  EffLv3AbsEtaHisto->Write();
  FracWrgChaLv3EtaHisto->Write();
  FracWrgChaLv3AbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffLv3PhiHisto[iDet]->Write();

  EffGLPtHisto->Write();
  FracWrgChaGLPtHisto->Write();
  EffGLEtaHisto->Write();
  EffGLAbsEtaHisto->Write();
  FracWrgChaGLEtaHisto->Write();
  FracWrgChaGLAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffGLPhiHisto[iDet]->Write();

  EffGLTrackPtHisto->Write();
  FracWrgChaGLTrackPtHisto->Write();
  EffGLTrackEtaHisto->Write();
  EffGLTrackAbsEtaHisto->Write();
  FracWrgChaGLTrackEtaHisto->Write();
  FracWrgChaGLTrackAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffGLTrackPhiHisto[iDet]->Write();

  EffSTATrackPtHisto->Write();
  FracWrgChaSTATrackPtHisto->Write();
  EffSTATrackEtaHisto->Write();
  EffSTATrackAbsEtaHisto->Write();
  FracWrgChaSTATrackEtaHisto->Write();
  FracWrgChaSTATrackAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffSTATrackPhiHisto[iDet]->Write();

  EffTKPtHisto->Write();
  FracWrgChaTKPtHisto->Write();
  EffTKEtaHisto->Write();
  EffTKAbsEtaHisto->Write();
  FracWrgChaTKEtaHisto->Write();
  FracWrgChaTKAbsEtaHisto->Write();
  for (int iDet=0;iDet<4;iDet++) EffTKPhiHisto[iDet]->Write();

  if (theRootTree) {
    theSumTree->Write();
    theTree->Write();
  }

  theFile->Close();
}


void TrackerMuonAnalyzer::analyze(const Event & event, const EventSetup& eventSetup){

  run =  event.id().run();
  evt =  event.id().event();

  if (theFirstEvent) {
    theFirstEvent = false;
    // Get the DT Geometry
    eventSetup.get<MuonGeometryRecord>().get(muonDTGeom);
  // Get the Tracker Geometry
    //    eventSetup.get<TrackerDigiGeometryRecord>().get(trackerGeom);

    if (!muonDTGeom.isValid()) {
      std::cout << " Unable to find MuonGeometryRecord for the DTGeometry in event!" << std::endl;
      return;
    }


    


  }

  float theLightestMuonP=1E6;    float theLightestAntiMuonP=1E6;
  float theHeaviestMuonP=0.;     float theHeaviestAntiMuonP=0.;
  float theLightestMuonPhi=1E6;  float theLightestAntiMuonPhi=1E6;
  float theHeaviestMuonPhi=1E6;  float theHeaviestAntiMuonPhi=1E6;
  float theLightestMuonPT=1E6;   float theLightestAntiMuonPT=1E6;
  float theHeaviestMuonPT=0.;    float theHeaviestAntiMuonPT=0.;

  if (theDebug) {
    cout << "Run: " << run << " Event: " << evt << endl;
  }
  numberOfEvents++;
  if ( numberOfEvents%500 == 0 )
    cout << " --- Processed " << numberOfEvents << " events" << endl;
  
  if (theDebug && isFastSim) {
    // Get the whole SimTrack collection in the FastSim events
    Handle<SimTrackContainer> allsimTracks;
    event.getByLabel("famosSimHits",allsimTracks);
    SimTrackContainer::const_iterator allsimTrack;
    cout << " LIST OF ALL SIMTRACKS IN THE EVENT " << endl;
    for (allsimTrack = allsimTracks->begin(); allsimTrack != allsimTracks->end(); ++allsimTrack){
      cout << " (" << (*allsimTrack).type() << ") -> " 
	   << (*allsimTrack).momentum().eta() << " , "
	   << sqrt(pow((*allsimTrack).momentum().px(),2)+pow((*allsimTrack).momentum().py(),2) ) << endl;
    }
    cout << " END OF LIST OF ALL SIMTRACKS " << endl;
  }


  // Get the SimTrack collection from the event
  Handle<SimTrackContainer> simTracks;
  event.getByLabel(theSimTrackLabel,simTracks);
  SimTrackContainer::const_iterator simTrack;
  nsimu =  simTracks->size();



  // Get the SimVertex collection from the event
  Handle<SimVertexContainer> simVertices;
  event.getByLabel(theSimVertexLabel,simVertices);
  SimVertexContainer::const_iterator simVertex;


  UInt_t nDTHitsMuon = 0;
  UInt_t nDTHitsAntiMuon = 0;
  UInt_t nCSCHitsMuon = 0;
  UInt_t nCSCHitsAntiMuon = 0;
  UInt_t nRPCHitsMuon = 0;
  UInt_t nRPCHitsAntiMuon = 0;



  // Look for hits in the muon detectors
  if (readMuonHits) {
    Handle<PSimHitContainer> muonDTHits;
    event.getByLabel(theDTHits,muonDTHits);
    Handle<PSimHitContainer> muonCSCHits;
    event.getByLabel(theCSCHits,muonCSCHits);
    Handle<PSimHitContainer> muonRPCHits;
    event.getByLabel(theRPCHits,muonRPCHits);
    // Check whether there are hits in DT/CSC/RPC,
    // and save the most/less energetic muon/antimuon:
    PSimHitContainer::const_iterator simMuonHit=muonDTHits->begin();
    PSimHitContainer::const_iterator endDTHit=muonDTHits->end();
    PSimHitContainer::const_iterator endCSCHit=muonCSCHits->end();
    PSimHitContainer::const_iterator endRPCHit=muonRPCHits->end();
    for ( ; simMuonHit!=endDTHit; ++simMuonHit) {
      float theP = simMuonHit->pabs();
      if (simMuonHit->particleType() == 13) {
	const DTWireId wireId( (*simMuonHit).detUnitId() );
	const DTLayer* layer = muonDTGeom->layer(wireId.layerId()); 
	LocalVector localMomentumAtEntry = (*simMuonHit).momentumAtEntry();
	GlobalVector globalMomentumAtEntry = layer->toGlobal(localMomentumAtEntry);
	if (theP<theLightestMuonP) {
	  theLightestMuonP = theP;
	  theLightestMuonPT = globalMomentumAtEntry.perp();
	  theLightestMuonPhi = globalMomentumAtEntry.phi();
	}
	if (theP>theHeaviestMuonP) {
	  theHeaviestMuonP = theP;
	  theHeaviestMuonPT = globalMomentumAtEntry.perp();
	  theHeaviestMuonPhi = globalMomentumAtEntry.phi();
	}
	nDTHitsMuon++;
      }
      else if (simMuonHit->particleType() == -13) {
	const DTWireId wireId( (*simMuonHit).detUnitId() );
	const DTLayer* layer = muonDTGeom->layer(wireId.layerId()); 
	LocalVector localMomentumAtEntry = (*simMuonHit).momentumAtEntry();
	GlobalVector globalMomentumAtEntry = layer->toGlobal(localMomentumAtEntry);
	if (theP<theLightestAntiMuonP) {
	  theLightestAntiMuonP = theP;
	  theLightestAntiMuonPT = globalMomentumAtEntry.perp();
	  theLightestAntiMuonPhi = globalMomentumAtEntry.phi();
	}
	if (theP>theHeaviestAntiMuonP) {
	  theHeaviestAntiMuonP = theP;
	  theHeaviestAntiMuonPT = globalMomentumAtEntry.perp();
	  theHeaviestAntiMuonPhi = globalMomentumAtEntry.phi();
	}
	nDTHitsAntiMuon++;
      }
    }
    simMuonHit=muonCSCHits->begin();
    for ( ; simMuonHit!=endCSCHit; ++simMuonHit) {
      float theP = simMuonHit->pabs();
      if (simMuonHit->particleType() == 13) {
	if (theP<theLightestMuonP) theLightestMuonP = theP;
	if (theP>theHeaviestMuonP) theHeaviestMuonP = theP;
	nCSCHitsMuon++;
      }
      else if (simMuonHit->particleType() == -13) {
	if (theP<theLightestAntiMuonP) theLightestAntiMuonP = theP;
	if (theP>theHeaviestAntiMuonP) theHeaviestAntiMuonP = theP;
	nCSCHitsAntiMuon++;
      }
    }
    simMuonHit=muonRPCHits->begin();
    for ( ; simMuonHit!=endRPCHit; ++simMuonHit) {
      float theP = simMuonHit->pabs();
      if (simMuonHit->particleType() == 13) {
	if (theP<theLightestMuonP) theLightestMuonP = theP;
	if (theP>theHeaviestMuonP) theHeaviestMuonP = theP;
	nRPCHitsMuon++;
      }
      else if (simMuonHit->particleType() == -13) {
	if (theP<theLightestAntiMuonP) theLightestAntiMuonP = theP;
	if (theP>theHeaviestAntiMuonP) theHeaviestAntiMuonP = theP;
	nRPCHitsAntiMuon++;
      }
    }
  }


  // Get the reco tracks
  Handle<reco::TrackCollection> theTracksFast;
  reco::TrackCollection::const_iterator trkFast;
  ntkFast=0;
  if (theRecoTrackLabelFast.label() != "MISSING") {
    event.getByLabel(theRecoTrackLabelFast,theTracksFast);
    trkFast = theTracksFast->begin();
    ntkFast =  theTracksFast->size();
  }


  // Get the L1 Muons
  Handle<L1MuonCollection> l1MuonsFast;
  L1MuonCollection::const_iterator l1mFast;
  nl1Fast=0;
  if (theL1MuonLabelFast.label() != "MISSING") {
    if (event.getByLabel(theL1MuonLabelFast,l1MuonsFast) || !ignoreMissingCollections) {
      l1mFast = l1MuonsFast->begin();
      nl1Fast = l1MuonsFast->size();
    }
  }


  // Get the L2 Muons
  Handle<reco::TrackCollection>  l2MuonsFast;
  reco::TrackCollection::const_iterator  l2mFast;
  nl2Fast=0;
  if (theL2MuonLabelFast.label() != "MISSING") {
    if (event.getByLabel(theL2MuonLabelFast,l2MuonsFast) || !ignoreMissingCollections) {
      l2mFast = l2MuonsFast->begin();
      nl2Fast = l2MuonsFast->size();
    }
  }


  // Get the L3 Muons
  Handle<reco::TrackCollection>  l3MuonsFast;
  reco::TrackCollection::const_iterator  l3mFast;
  nl3Fast=0;
  if (theL3MuonLabelFast.label() != "MISSING") {
    if (event.getByLabel(theL3MuonLabelFast,l3MuonsFast) || !ignoreMissingCollections) {
      l3mFast = l3MuonsFast->begin();
      nl3Fast = l3MuonsFast->size();
    }
  }


  // Get the Global Muons
  Handle<reco::MuonCollection> glMuonsFast;
  reco::MuonCollection::const_iterator glmFast;
  nglFast=0;
  if (theGlobalMuonLabelFast.label() != "MISSING") {
    event.getByLabel(theGlobalMuonLabelFast,glMuonsFast);
    glmFast = glMuonsFast->begin();
    nglFast = glMuonsFast->size();
  }

 
  // Get the Global Muons Tracks (from MuonSimHits)
  Handle<reco::TrackCollection> glMuonTracksFast;
  reco::TrackCollection::const_iterator gltkFast;
  ngltkFast=0;
  if (theGlobalMuonTrackLabelFast.label() != "MISSING") {
    event.getByLabel(theGlobalMuonTrackLabelFast,glMuonTracksFast);
    gltkFast = glMuonTracksFast->begin();
    ngltkFast = glMuonTracksFast->size();
  }

 
  // Get the StanAlone Muons Tracks (from MuonSimHits)
  Handle<reco::TrackCollection> staMuonTracksFast;
  reco::TrackCollection::const_iterator statkFast;
  nstatkFast=0;
  if (theStandAloneMuonTrackLabelFast.label() != "MISSING") {
    event.getByLabel(theStandAloneMuonTrackLabelFast,staMuonTracksFast);
    statkFast = staMuonTracksFast->begin();
    nstatkFast = staMuonTracksFast->size();
  }
 
  if (theRootTree) {
    Bool_t someSimu = (nsimu>0);
    Bool_t someFast = (ntkFast>0) || (nl1Fast>0) || (nl2Fast) || (nl3Fast>0) ||
                      (nglFast>0) || (ngltkFast>0) || (nstatkFast>0);
    if (someSimu || someFast) theSumTree->Fill();
  }


  if (theDebug) {
    cout << " -- Simulated tracks: " << nsimu << endl;
    cout << " -- tk , l1 , l2 , l3, gl , gltk, statk : "
	 << ntkFast << " , " << nl1Fast << " , " << nl2Fast << " , " << nl3Fast << " , " 
	 << nglFast << " , " << ngltkFast << " , " << nstatkFast << endl;
  }

  //
  // Start loop on Sim Muons
  //

  unsigned int iMu = 0;
  unsigned int iTk = 0;

  for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
    // Only look for simulated muons in FullSim
    if (!isFastSim && abs((*simTrack).type()) != 13) continue;

    if (theDebug) {
      cout << " SimTrack # " << iTk << " (" << (*simTrack).type() << ") -> " 
	   << (*simTrack).momentum().eta() << " , "
	   << (*simTrack).momentum().phi() << " , "
	   << "PT TrackerIn_FullSim: "<<sqrt(pow((*simTrack).momentum().px(),2)+pow((*simTrack).momentum().py(),2)) << " , "
           << "PT TrackerOut_FullSim: "<<sqrt(pow((*simTrack).trackerSurfaceMomentum().px(),2)+pow((*simTrack).trackerSurfaceMomentum().py(),2)) << endl;
    }
    iTk++;

    if (iMu >= nMuMax) {
      iMu++;
      continue;
    }

    // Look for the mother only in FastSim (at the moment), where just
    // muons and their mothers are saved in the simMuon collection.
    if (isFastSim && abs((*simTrack).type()) != 13) {
      pdgMot[iMu] = (*simTrack).type();      
      ptMot[iMu] =  sqrt(pow((*simTrack).momentum().px(),2) + 
			 pow((*simTrack).momentum().py(),2) );
      phiMot[iMu] = (*simTrack).momentum().phi();
      etaMot[iMu] = (*simTrack).momentum().eta();
      ++simTrack;
      if ( simTrack!=simTracks->end() && abs((*simTrack).type())==13 ) {
	if (theDebug) {
	  cout << " Daughter SimTrack # " << iMu << " (" << (*simTrack).type() << ") -> " 
	       << (*simTrack).momentum().eta() << " , "
	       << sqrt(pow((*simTrack).momentum().px(),2)+pow((*simTrack).momentum().py(),2) ) << endl;
	}
      }
      else {
	ptMot[iMu] =  -1;
	phiMot[iMu] = -999;
	etaMot[iMu] = -999;
	cout << " WARNING: mother and daughter are both non-muons !" << endl;
	continue;
      }
    }
    else {
      // FullSim muons end up here
      pdgMot[iMu] =  0;
      ptMot[iMu] =  -1;
      phiMot[iMu] = -999;
      etaMot[iMu] = -999;
    }

    pSim[iMu] = (*simTrack).momentum().P();
    pTrackerOut[iMu] = (*simTrack).trackerSurfaceMomentum ().P();
    

    ptTrackerIn[iMu] = sqrt(pow((*simTrack).momentum().px(),2) +
                             pow((*simTrack).momentum().py(),2) );

    ptTrackerOut[iMu] = sqrt(pow((*simTrack).trackerSurfaceMomentum().px(),2) + 
			     pow((*simTrack).trackerSurfaceMomentum().py(),2) );


  
    trackerPTIn->Fill(ptTrackerIn[iMu]);
 
    trackerPTOn->Fill(ptTrackerOut[iMu]);
     
    phiTrackerOut[iMu] = (*simTrack).trackerSurfaceMomentum ().phi();
    
    if ((*simTrack).type()==13) {
      pMuonIn[iMu]  = theHeaviestMuonP;
      pMuonOut[iMu] = theLightestMuonP;
      ptMuonIn[iMu]  = theHeaviestMuonPT;
      ptMuonOut[iMu] = theLightestMuonPT;
      phiMuonIn[iMu]  = theHeaviestMuonPhi;
      phiMuonOut[iMu] = theLightestMuonPhi;
      nDTHits[iMu]  = nDTHitsMuon;
      nCSCHits[iMu] = nCSCHitsMuon;
      nRPCHits[iMu] = nRPCHitsMuon;
    }
    else if  ((*simTrack).type()==-13) {
      pMuonIn[iMu]  = theHeaviestAntiMuonP;
      pMuonOut[iMu] = theLightestAntiMuonP;
      ptMuonIn[iMu]  = theHeaviestAntiMuonPT;
      ptMuonOut[iMu] = theLightestAntiMuonPT;
      phiMuonIn[iMu]  = theHeaviestAntiMuonPhi;
      phiMuonOut[iMu] = theLightestAntiMuonPhi;
      nDTHits[iMu]  = nDTHitsAntiMuon;
      nCSCHits[iMu] = nCSCHitsAntiMuon;
      nRPCHits[iMu] = nRPCHitsAntiMuon;
    }

    pdgSim[iMu] = (*simTrack).type();
    ptSim[iMu] =  sqrt(pow((*simTrack).momentum().px(),2) + 
		       pow((*simTrack).momentum().py(),2) );

    //   trackerPTIn->Fill();

    phiSim[iMu] = (*simTrack).momentum().phi();
    etaSim[iMu] = (*simTrack).momentum().eta();
    chaSim[iMu] = (Int_t) ( (*simTrack).charge() * 1.05 );

    Float_t xpt = ptSim[iMu];
    Int_t cha = chaSim[iMu];
    Float_t eta =  etaSim[iMu];
    Float_t abseta = abs(eta);
    Float_t phi = phiSim[iMu]*180./3.14159265;
    if (phi<0) phi += 360.;

    int iDet = -1;
    if (abseta <= 0.9) 	iDet = 0;
    else if (abseta <= 1.3)   iDet = 1;
    else if (abseta <= 2.1)   iDet = 2;
    else if (abseta <= 2.4)   iDet = 3;

    //      bool ptcut = xpt > 10.; => Metti in un param. del cfg
    bool ptcut = xpt > 8.;
    //    bool etacut = abseta<=2.4;
    
    Float_t phiDet=-1.;
    if (iDet>=0) {
      if (iDet==0 || iDet ==1) phiDet = 360./SegmentationBarrel;
      else if (iDet==2 || iDet==3) phiDet = 360./SegmentationEndcap;
      Int_t phiInt = Int_t (phi/phiDet);
      phiDet = phi - phiDet*phiInt;
    }

    
    /*
    // Look for the production vertex of the SimTrack
    if (theDebug) {
      int vtxIndex = (*simTrack).vertIndex();
      int numVtx = 0;
      if (vtxIndex>=0) {
	for (simVertex = simVertices->begin(); simVertex != simVertices->end(); ++simVertex) {
	  double x = (*simVertex).position().x();
	  double y = (*simVertex).position().y();
	  double radius = sqrt (x*x+y*y);
	  cout << vtxIndex << " : " <<  numVtx << " -> " << radius 
	       << " from " << (*simVertex).position() << endl;
	  numVtx++;
	}
      }
    }
    */

    if (iDet>=0) {
      GenPtHisto->Fill(xpt);
      GenPtDetHisto[iDet]->Fill(xpt);
    }
    if (ptcut) {
      GenEtaHisto->Fill(eta);
      GenAbsEtaHisto->Fill(abseta);
      if (iDet>=0) GenPhiHisto[iDet]->Fill(phiDet);
    }
    math::XYZVector mcDir( (*simTrack).momentum() );
    bool recoExists;
    double resDr = 0.5;
    

    if (nl1Fast>0) {
      recoExists = false;
      double minDr = resDr;
      for (l1mFast = l1MuonsFast->begin(); l1mFast != l1MuonsFast->end(); ++l1mFast){
	math::XYZVector recoDir( (*l1mFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptL1[iMu]  = (*l1mFast).pt();
	  phiL1[iMu] = (*l1mFast).phi();
	  etaL1[iMu] = (*l1mFast).eta();
	  chaL1[iMu] = (*l1mFast).charge();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaL1[iMu] != cha;
	if (iDet>=0) {
	  Lv1PtHisto[iDet]->Fill(xpt);
	  Lv1PtOutHisto[iDet]->Fill(ptL1[iMu]);
	  if (wrongCharge) Lv1ChaPtHisto[iDet]->Fill(xpt);
	}
	if (ptcut) {
	  Lv1EtaHisto->Fill(eta);
	  Lv1AbsEtaHisto->Fill(abseta);
	  if (wrongCharge) Lv1ChaEtaHisto->Fill(eta);
	  if (wrongCharge) Lv1ChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) Lv1PhiHisto[iDet]->Fill(phiDet);          
	}
      }
      else{
	ptL1[iMu]  = -1.;
	phiL1[iMu] = -999.;
	etaL1[iMu] = -999.;
	chaL1[iMu] = -999 ;
      }
    }
    

    if (nl2Fast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      for (l2mFast = l2MuonsFast->begin(); l2mFast != l2MuonsFast->end(); ++l2mFast){
	math::XYZVector recoDir( (*l2mFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptL2[iMu]  = (*l2mFast).pt();
	  phiL2[iMu] = (*l2mFast).phi();
	  etaL2[iMu] = (*l2mFast).eta();
	  chaL2[iMu] = (*l2mFast).charge();
          ptres = (*l2mFast).pt()/xpt;
	  invptres = xpt/(*l2mFast).pt();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaL2[iMu] != cha;
	if (iDet>=0) {
	  Lv2PtHisto->Fill(xpt);
	  if (wrongCharge) Lv2ChaPtHisto->Fill(xpt);
	}
	if (ptcut) {
	  Lv2EtaHisto->Fill(eta);
	  Lv2AbsEtaHisto->Fill(abseta);
	  if (wrongCharge) Lv2ChaEtaHisto->Fill(eta);
	  if (wrongCharge) Lv2ChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) Lv2PhiHisto[iDet]->Fill(phiDet);
	}
	Lv2PtResTotHisto->Fill(ptres);
	if (iDet>=0) Lv2PtResHisto[iDet]->Fill(ptres);
	Lv2InvPtResTotHisto->Fill(invptres);
	if (iDet>=0) Lv2InvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptL2[iMu]  = -1.;
	phiL2[iMu] = -999.;
	etaL2[iMu] = -999.;
      }
    }
    
    
    if (nl3Fast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      for (l3mFast = l3MuonsFast->begin(); l3mFast != l3MuonsFast->end(); ++l3mFast){
	math::XYZVector recoDir( (*l3mFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptL3[iMu]  = (*l3mFast).pt();
	  phiL3[iMu] = (*l3mFast).phi();
	  etaL3[iMu] = (*l3mFast).eta();
	  chaL3[iMu] = (*l3mFast).charge();
	  ptres = (*l3mFast).pt()/xpt;
	  invptres = xpt/(*l3mFast).pt();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaL3[iMu] != cha;
	if (iDet>=0) {
	  Lv3PtHisto->Fill(xpt);
	  if (wrongCharge) Lv3ChaPtHisto->Fill(xpt);
	}
	if (ptcut) {
	  Lv3EtaHisto->Fill(eta);
	  Lv3AbsEtaHisto->Fill(abseta);
	  if (wrongCharge) Lv3ChaEtaHisto->Fill(eta);
	  if (wrongCharge) Lv3ChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) Lv3PhiHisto[iDet]->Fill(phiDet);
	}
	Lv3PtResTotHisto->Fill(ptres);
	if (iDet>=0) Lv3PtResHisto[iDet]->Fill(ptres);
	Lv3InvPtResTotHisto->Fill(invptres);
	if (iDet>=0) Lv3InvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptL3[iMu]  = -1.;
	phiL3[iMu] = -999.;
	etaL3[iMu] = -999.;
      }
    }
    
    
    if (nglFast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      for (glmFast = glMuonsFast->begin(); glmFast != glMuonsFast->end(); ++glmFast){
	math::XYZVector recoDir( (*glmFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptGL[iMu] = (*glmFast).pt();
	  phiGL[iMu] = (*glmFast).phi();
	  etaGL[iMu] = (*glmFast).eta();
	  chaGL[iMu] = (*glmFast).charge();
	  ptres = (*glmFast).pt()/xpt;
	  invptres = xpt/ptGL[iMu];
	}
      }
      if (recoExists) {
        bool wrongCharge = chaGL[iMu] != cha;
	if (iDet>=0) {
	  GLPtHisto->Fill(xpt);
	  if (wrongCharge) GLChaPtHisto->Fill(xpt);
	}
	if (ptcut) {
	  GLEtaHisto->Fill(eta);
	  GLAbsEtaHisto->Fill(abseta);
	  if (wrongCharge) GLChaEtaHisto->Fill(eta);
	  if (wrongCharge) GLChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) GLPhiHisto[iDet]->Fill(phiDet);
	}
	GLPtResTotHisto->Fill(ptres);
	if (iDet>=0) GLPtResHisto[iDet]->Fill(ptres);
	GLInvPtResTotHisto->Fill(invptres);
	if (iDet>=0) GLInvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptGL[iMu]  = -1.;
	phiGL[iMu] = -999.;
	etaGL[iMu] = -999.;
      }
    }
    
    
    if (ntkFast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      for (trkFast = theTracksFast->begin(); trkFast != theTracksFast->end(); ++trkFast){
	math::XYZVector recoDir( (*trkFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptTK[iMu]  = (*trkFast).pt();
	  phiTK[iMu] = (*trkFast).phi();
	  etaTK[iMu] = (*trkFast).eta();
	  chaTK[iMu] = (*trkFast).charge();
	  ptres = (*trkFast).pt()/xpt;
	  invptres = xpt/(*trkFast).pt();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaTK[iMu] != cha;
	if (iDet>=0) {
	  TKPtHisto->Fill(xpt);
	  if (wrongCharge) TKChaPtHisto->Fill(xpt);
	}
	TKPtResTotHisto->Fill(ptres);
	if (iDet>=0)  TKPtResHisto[iDet]->Fill(ptres);
	if (ptcut) {
	  TKEtaHisto->Fill(eta);
	  TKAbsEtaHisto->Fill(abseta);
	  if (wrongCharge) TKChaEtaHisto->Fill(eta);
	  if (wrongCharge) TKChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) TKPhiHisto[iDet]->Fill(phiDet);
	}
	TKPtResTotHisto->Fill(ptres);
	if (iDet>=0) TKPtResHisto[iDet]->Fill(ptres);
	TKInvPtResTotHisto->Fill(invptres);
	if (iDet>=0) TKInvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptTK[iMu]  = -1.;
	phiTK[iMu] = -999.;
	etaTK[iMu] = -999.;
      }
    }
    
    
    if (ngltkFast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      for (gltkFast = glMuonTracksFast->begin(); gltkFast != glMuonTracksFast->end(); ++gltkFast){
	math::XYZVector recoDir( (*gltkFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptGLTK[iMu] = (*gltkFast).pt();
	  phiGLTK[iMu] = (*gltkFast).phi();
	  etaGLTK[iMu] = (*gltkFast).eta();
	  chaGLTK[iMu] = (*gltkFast).charge();
	  ptres = (*gltkFast).pt()/xpt;
	  invptres = xpt/(*gltkFast).pt();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaGLTK[iMu] != cha;
	if (iDet>=0) {
	  GLTrackPtHisto->Fill(xpt);
	  if (wrongCharge) GLTrackChaPtHisto->Fill(xpt);
	}
	if (ptcut) {
	  GLTrackEtaHisto->Fill(eta);
	  GLTrackAbsEtaHisto->Fill(abseta);
	  if (wrongCharge) GLTrackChaEtaHisto->Fill(eta);
	  if (wrongCharge) GLTrackChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) GLTrackPhiHisto[iDet]->Fill(phiDet);
	}
	GLTrackPtResTotHisto->Fill(ptres);
	if (iDet>=0) GLTrackPtResHisto[iDet]->Fill(ptres);
	GLTrackInvPtResTotHisto->Fill(invptres);
	if (iDet>=0) GLTrackInvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptGLTK[iMu]  = -1.;
	phiGLTK[iMu] = -999.;
	etaGLTK[iMu] = -999.;
      }
    }
    
    
    if (nstatkFast>0) {
      recoExists = false;
      double minDr = resDr;
      Float_t ptres = 999.;
      Float_t invptres = 999.;
      Float_t chi2 = -999.;
      Float_t ndof = -999.;
      int iSta = 0;
      for (statkFast = staMuonTracksFast->begin(); statkFast != staMuonTracksFast->end(); ++statkFast){
	math::XYZVector recoDir( (*statkFast).momentum() );
	double tmpDr = ROOT::Math::VectorUtil::DeltaR(mcDir, recoDir);
	if (theDebug) {
	  cout << " STATrack # " << iSta << " -> "
	       << (*statkFast).eta() << " , "
	       << (*statkFast).phi() << " , " << (*statkFast).pt() << endl;
	  iSta++;
	}
	if ( tmpDr < minDr ) {
	  recoExists = true;
	  minDr = tmpDr;
	  ptSTATK[iMu] = (*statkFast).pt();
	  phiSTATK[iMu] = (*statkFast).phi();
	  etaSTATK[iMu] = (*statkFast).eta();
	  chaSTATK[iMu] = (*statkFast).charge();
          chi2STATK[iMu] = (*statkFast).chi2();
          ndofSTATK[iMu] = (*statkFast).ndof();
	  dtHitsSTATK[iMu] = (*statkFast).hitPattern().numberOfValidMuonDTHits();
	  cscHitsSTATK[iMu] = (*statkFast).hitPattern().numberOfValidMuonCSCHits();
	  rpcHitsSTATK[iMu] = (*statkFast).hitPattern().numberOfValidMuonRPCHits();
	  ptres = (*statkFast).pt()/xpt;
	  invptres = xpt/(*statkFast).pt();
	  chi2 = (*statkFast).chi2();
	  ndof = (*statkFast).ndof();
	}
      }
      if (recoExists) {
        bool wrongCharge = chaSTATK[iMu] != cha;
	if (iDet>=0) {
	  STATrackPtHisto->Fill(xpt);
	  if (wrongCharge) STATrackChaPtHisto->Fill(xpt);
	}
	if (ptcut) {
	  STATrackEtaHisto->Fill(eta);
	  STATrackAbsEtaHisto->Fill(abseta);
	  if (wrongCharge) STATrackChaEtaHisto->Fill(eta);
	  if (wrongCharge) STATrackChaAbsEtaHisto->Fill(abseta);
	  if (iDet>=0) STATrackPhiHisto[iDet]->Fill(phiDet);
	  STATrackChi2TotHisto->Fill(chi2);
	  STATrackNdofTotHisto->Fill(ndof);
	  if (ndof>0) STATrackNormChi2TotHisto->Fill(chi2/ndof);
	  if (iDet>=0) {
	    STATrackChi2Histo[iDet]->Fill(chi2);
	    STATrackNdofHisto[iDet]->Fill(ndof);
	    if (ndof>0) STATrackNormChi2Histo[iDet]->Fill(chi2/ndof);
	  }
	}
	STATrackPtResTotHisto->Fill(ptres);
	if (iDet>=0) STATrackPtResHisto[iDet]->Fill(ptres);
	STATrackInvPtResTotHisto->Fill(invptres);
	if (iDet>=0) STATrackInvPtResHisto[iDet]->Fill(invptres);
      }
      else {
	ptSTATK[iMu]  = -1.;
	phiSTATK[iMu] = -999.;
	etaSTATK[iMu] = -999.;
	chaSTATK[iMu] = 0.;
	chi2STATK[iMu] = -999.;
	ndofSTATK[iMu] = -999;
	dtHitsSTATK[iMu] = 0;
	cscHitsSTATK[iMu] = 0;
	rpcHitsSTATK[iMu] = 0;
      }
    }
    
    iMu++;
    
  } // End of the loop over this simulated muon in the event
  
  
  nSimMu = iMu;
  if (theRootTree) {
    if (iMu<nMuMax) theTree->Fill();
    else if (iMu>0) cout << " Warning: " << iMu << " simulated muons found ...";
  }
  if (theDebug) cout<<"-----"<<endl;  

}

DEFINE_FWK_MODULE(TrackerMuonAnalyzer);

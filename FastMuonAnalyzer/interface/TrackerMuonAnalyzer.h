#ifndef FastMuonBremAnalyzers_FastMuonAnalyzer_TrackerMuonAnalyzer_H
#define FastMuonBremAnalyzers_FastMuonAnalyzer_TrackerMuonAnalyzer_H

#define nMuMax 6
//#define nMuMax 1


// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
// Root headers
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"

// Geometry
#include "Geometry/DTGeometry/interface/DTGeometry.h"

// Data Formats
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/L1Trigger/interface/L1MuonParticle.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}


class TrackerMuonAnalyzer: public edm::EDAnalyzer {
public:
  /// Constructor
  TrackerMuonAnalyzer(const edm::ParameterSet& pset);

  /// Destructor
  virtual ~TrackerMuonAnalyzer();

  // Operations

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);

  virtual void beginJob() ;
  virtual void endJob() ;
protected:

private:

  edm::ESHandle<DTGeometry> muonDTGeom;
  //  edm::ESHandle<TrackerGeometry> trackerGeom;
 
  std::string theRootFileName;

  edm::InputTag  theSimTrackLabel;
  edm::InputTag  theSimVertexLabel;
  edm::InputTag  theDTHits;
  edm::InputTag  theCSCHits;
  edm::InputTag  theRPCHits;
  edm::InputTag  theRecoTrackLabelFast;
  edm::InputTag  theL1MuonLabelFast;
  edm::InputTag  theL2MuonLabelFast;
  edm::InputTag  theL3MuonLabelFast;
  edm::InputTag  theGlobalMuonLabelFast;
  edm::InputTag  theGlobalMuonTrackLabelFast;
  edm::InputTag  theStandAloneMuonTrackLabelFast;

 /*  edm::InputTag  theSimHitTOBLabel; */
/*   edm::InputTag  theSimHitTIBLabel; */
 

 TFile* theFile;

  bool ignoreMissingCollections;
  bool theRootTree;
  bool theDebug;
  bool isFastSim;
  bool readMuonHits;
  bool useGen;
  bool useReco;
  bool useEff;
  bool test_Track;

 

  /*
  // Histograms
  TH1F * hPtDiff;
  */

  TTree *theTree;
  TTree *theSumTree;

  Int_t run;
  Int_t evt;
  UInt_t nsimu;
  UShort_t ntkFast;
  UShort_t nl1Fast;
  UShort_t nl2Fast;
  UShort_t nl3Fast;
  UShort_t nglFast;
  UShort_t ngltkFast;
  UShort_t nstatkFast;

  UInt_t  nSimMu;
  Int_t   pdgSim[nMuMax];
  Float_t pSim[nMuMax];

  Float_t abspSim[nMuMax];
  Float_t abspTrackerOut[nMuMax];
  Float_t pTrackerOut[nMuMax];
 
 Float_t pMuonIn[nMuMax];
  Float_t pMuonOut[nMuMax];
  Float_t ptSim[nMuMax];
  Float_t ptTrackerIn[nMuMax];
  Float_t ptTrackerOut[nMuMax];
  Float_t ptMuonIn[nMuMax];
  Float_t ptMuonOut[nMuMax];
  Float_t phiTrackerOut[nMuMax];
  Float_t phiMuonIn[nMuMax];
  Float_t phiMuonOut[nMuMax];
  Float_t phiSim[nMuMax];
  Float_t etaSim[nMuMax];
  Int_t   chaSim[nMuMax];
  Int_t   pdgMot[nMuMax];
  Float_t ptMot[nMuMax];
  Float_t phiMot[nMuMax];
  Float_t etaMot[nMuMax];
  UInt_t  nDTHits[nMuMax];
  UInt_t  nCSCHits[nMuMax];
  UInt_t  nRPCHits[nMuMax];
  Float_t ptL1[nMuMax];
  Float_t phiL1[nMuMax];
  Float_t etaL1[nMuMax];
  Int_t   chaL1[nMuMax];
  Float_t ptL2[nMuMax];
  Float_t phiL2[nMuMax];
  Float_t etaL2[nMuMax];
  Int_t   chaL2[nMuMax];
  Float_t ptL3[nMuMax];
  Float_t phiL3[nMuMax];
  Float_t etaL3[nMuMax];
  Int_t   chaL3[nMuMax];
  Float_t ptGL[nMuMax];
  Float_t phiGL[nMuMax];
  Float_t etaGL[nMuMax];
  Int_t   chaGL[nMuMax];
  Float_t ptTK[nMuMax];
  Float_t phiTK[nMuMax];
  Float_t etaTK[nMuMax];
  Int_t   chaTK[nMuMax];
  Float_t ptSTATK[nMuMax];
  Float_t phiSTATK[nMuMax];
  Float_t etaSTATK[nMuMax];
  Int_t   chaSTATK[nMuMax];
  Float_t chi2STATK[nMuMax];
  Float_t ndofSTATK[nMuMax];
  Int_t   dtHitsSTATK[nMuMax];
  Int_t   cscHitsSTATK[nMuMax];
  Int_t   rpcHitsSTATK[nMuMax];
  Float_t ptGLTK[nMuMax];
  Float_t phiGLTK[nMuMax];
  Float_t etaGLTK[nMuMax];
  Int_t   chaGLTK[nMuMax];

// List of Histos to be filled :
////////////////////////////////////
//Tracker plots(S.Fonseca)
  TH1F *trackerPTOn;
  TH1F *trackerPTIn;
  TH1F *DifftrackerPoutPgen;
  TH1F *DifftrackerPoutPgenInvPgen;
  /////////////////////////////////
  TH1F *GenPtHisto;
  TH1F *GenPtDetHisto[4];
  TH1F *GenEtaHisto , *GenAbsEtaHisto;
  TH1F *GenPhiHisto[4];

  TH1F *Lv1PtHisto[4];
  TH1F *Lv1ChaPtHisto[4];
  TH1F *Lv1PtOutHisto[4];
  TH1F *Lv1EtaHisto , *Lv1AbsEtaHisto;
  TH1F *Lv1ChaEtaHisto , *Lv1ChaAbsEtaHisto;
  TH1F *Lv1PhiHisto[4];

  TH1F *Lv2PtHisto;
  TH1F *Lv2ChaPtHisto;
  TH1F *Lv2EtaHisto , *Lv2AbsEtaHisto;
  TH1F *Lv2ChaEtaHisto , *Lv2ChaAbsEtaHisto;
  TH1F *Lv2PhiHisto[4];
  TH1F *Lv2PtResTotHisto;
  TH1F *Lv2PtResHisto[4];
  TH1F *Lv2InvPtResTotHisto;
  TH1F *Lv2InvPtResHisto[4];

  TH1F *Lv3PtHisto;
  TH1F *Lv3ChaPtHisto;
  TH1F *Lv3EtaHisto , *Lv3AbsEtaHisto;
  TH1F *Lv3ChaEtaHisto , *Lv3ChaAbsEtaHisto;
  TH1F *Lv3PhiHisto[4];
  TH1F *Lv3PtResTotHisto;
  TH1F *Lv3PtResHisto[4];
  TH1F *Lv3InvPtResTotHisto;
  TH1F *Lv3InvPtResHisto[4];

  TH1F *GLPtHisto;
  TH1F *GLChaPtHisto;
  TH1F *GLEtaHisto , *GLAbsEtaHisto;
  TH1F *GLChaEtaHisto , *GLChaAbsEtaHisto;
  TH1F *GLPhiHisto[4];
  TH1F *GLPtResTotHisto;
  TH1F *GLPtResHisto[4];
  TH1F *GLInvPtResTotHisto;
  TH1F *GLInvPtResHisto[4];

  TH1F *GLTrackPtHisto;
  TH1F *GLTrackChaPtHisto;
  TH1F *GLTrackEtaHisto , *GLTrackAbsEtaHisto;
  TH1F *GLTrackChaEtaHisto , *GLTrackChaAbsEtaHisto;
  TH1F *GLTrackPhiHisto[4];
  TH1F *GLTrackPtResTotHisto;
  TH1F *GLTrackPtResHisto[4];
  TH1F *GLTrackInvPtResTotHisto;
  TH1F *GLTrackInvPtResHisto[4];

  TH1F *STATrackPtHisto;
  TH1F *STATrackChaPtHisto;
  TH1F *STATrackEtaHisto , *STATrackAbsEtaHisto;
  TH1F *STATrackChaEtaHisto , *STATrackChaAbsEtaHisto;
  TH1F *STATrackPhiHisto[4];
  TH1F *STATrackPtResTotHisto;
  TH1F *STATrackPtResHisto[4];
  TH1F *STATrackInvPtResTotHisto;
  TH1F *STATrackInvPtResHisto[4];
  TH1F *STATrackChi2TotHisto;
  TH1F *STATrackNdofTotHisto;
  TH1F *STATrackNormChi2TotHisto;
  TH1F *STATrackChi2Histo[4];
  TH1F *STATrackNdofHisto[4];
  TH1F *STATrackNormChi2Histo[4];

  TH1F *TKPtHisto;
  TH1F *TKChaPtHisto;
  TH1F *TKEtaHisto , *TKAbsEtaHisto;
  TH1F *TKChaEtaHisto , *TKChaAbsEtaHisto;
  TH1F *TKPhiHisto[4];
  TH1F *TKPtResTotHisto;
  TH1F *TKPtResHisto[4];
  TH1F *TKInvPtResTotHisto;
  TH1F *TKInvPtResHisto[4];
 
  TH1F *EffLv1PtHisto[4];
  TH1F *FracWrgChaLv1PtHisto[4];
  TH1F *FracWrgChaLv1EtaHisto , *FracWrgChaLv1AbsEtaHisto;
  TH1F *EffLv1EtaHisto , *EffLv1AbsEtaHisto;
  TH1F *EffLv1PhiHisto[4];

  TH1F *EffLv2PtHisto;
  TH1F *FracWrgChaLv2PtHisto;
  TH1F *EffLv2EtaHisto , *EffLv2AbsEtaHisto;
  TH1F *FracWrgChaLv2EtaHisto , *FracWrgChaLv2AbsEtaHisto;
  TH1F *EffLv2PhiHisto[4];

  TH1F *EffLv3PtHisto;
  TH1F *FracWrgChaLv3PtHisto;
  TH1F *EffLv3EtaHisto , *EffLv3AbsEtaHisto;
  TH1F *FracWrgChaLv3EtaHisto , *FracWrgChaLv3AbsEtaHisto;
  TH1F *EffLv3PhiHisto[4];

  TH1F *EffGLPtHisto;
  TH1F *FracWrgChaGLPtHisto;
  TH1F *EffGLEtaHisto , *EffGLAbsEtaHisto;
  TH1F *FracWrgChaGLEtaHisto , *FracWrgChaGLAbsEtaHisto;
  TH1F *EffGLPhiHisto[4];

  TH1F *EffGLTrackPtHisto;
  TH1F *FracWrgChaGLTrackPtHisto;
  TH1F *EffGLTrackEtaHisto , *EffGLTrackAbsEtaHisto;
  TH1F *FracWrgChaGLTrackEtaHisto , *FracWrgChaGLTrackAbsEtaHisto;
  TH1F *EffGLTrackPhiHisto[4];

  TH1F *EffSTATrackPtHisto;
  TH1F *FracWrgChaSTATrackPtHisto;
  TH1F *EffSTATrackEtaHisto , *EffSTATrackAbsEtaHisto;
  TH1F *FracWrgChaSTATrackEtaHisto , *FracWrgChaSTATrackAbsEtaHisto;
  TH1F *EffSTATrackPhiHisto[4];

  TH1F *EffTKPtHisto;
  TH1F *FracWrgChaTKPtHisto;
  TH1F *EffTKEtaHisto , *EffTKAbsEtaHisto;
  TH1F *FracWrgChaTKEtaHisto , *FracWrgChaTKAbsEtaHisto;
  TH1F *EffTKPhiHisto[4];


  // Counters
  int numberOfEvents;
  bool theFirstEvent;

  // Analyzer constant data
  int SegmentationBarrel , SegmentationEndcap;
};
#endif


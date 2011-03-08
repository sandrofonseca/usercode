import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )

# Geometry:
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.DTGeometryBuilder.dtGeometry_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_3XY_V21::All"
process.GlobalTag.globaltag = "MC_3XY_V26::All"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#     'file:/afs/cern.ch/user/a/aperrott/scratch0/CMSSW_2_2_0_pre1/src/APAnalyzers/RelVal220pre1SingleMuPt100.root'
#     'rfio:/castor/cern.ch/user/a/aperrott/FastSim/CMSSW_3_5_8_SingleMuPt1000_cfi_py_GEN_FASTSIM_VALIDATION.root'

#
## Single mu, pt = 10 GeV
#
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0017/DA1E5AF4-8E52-DF11-9F84-003048678A7E.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/FE6EB788-6552-DF11-A463-002618943849.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/A04B3587-6552-DF11-A7A1-00261894386B.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/8C8EFC03-6652-DF11-B40D-002618943829.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/7CAA1D87-6552-DF11-AF08-0026189438E3.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/1AF5D003-6652-DF11-AF5D-00304867901A.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/0C1CF085-6552-DF11-81CE-003048679266.root'


#
## Single mu, pt = 100 GeV
#
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0017/B4F83CAD-8F52-DF11-8206-0026189438D9.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/EAF21794-5752-DF11-B741-00261894390C.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/9E4EE303-5852-DF11-A88C-001A92971B82.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/3EC42D95-5752-DF11-8AEC-001A92971B04.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/3C74FB04-5852-DF11-99AD-001A928116EE.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/18E84C91-5752-DF11-887F-001A92971B82.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-DIGI-RECO/MC_3XY_V26_FastSim-v1/0016/10D0C592-5752-DF11-BA8F-001A928116EE.root'


## Single mu, pt = 1000 GeV
#
#     'rfio:/castor/cern.ch/user/a/aperrott/FastSim/CMSSW_3_5_8_SingleMuPt1000_cfi_py_GEN_FASTSIM_VALIDATION.root'
#FastSim standard
#'file:SingleMuPlusPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION.root'    
#FastSim muonBrem
'file:SingleMuPlusPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_MuonBrem.root'
    



    )
)

process.TrackerMuonAnalyzer = cms.EDAnalyzer("TrackerMuonAnalyzer",

    #rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE.root'),
    rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FASTfromFILE_MuonBrem.root'),
   
    FastSim = cms.untracked.bool(True),
    IgnoreMissingCollections = cms.untracked.bool(True),
    RootTree = cms.untracked.bool(False),
    Debug = cms.untracked.bool(True),
    useGen = cms.untracked.bool(False),
    useReco = cms.untracked.bool(False),
    useEff = cms.untracked.bool(False),
    test_Track= cms.untracked.bool(True),

    labelSimVertex = cms.InputTag("famosSimHits"),
    labelSimTrack = cms.InputTag("famosSimHits","MuonSimTracks"),
    dtSimHits = cms.InputTag("MuonSimHits","MuonDTHits"),
    cscSimHits = cms.InputTag("MuonSimHits","MuonCSCHits"),
    rpcSimHits = cms.InputTag("MuonSimHits","MuonRPCHits"),

  ##   #Tracker Labels
##     labelSimHitsTOB = cms.InputTag("famosSimHits","TrackerHits"),
##     labelSimHitsTIB = cms.InputTag("famosSimHits","TrackerHits"),   
                                             

    # VERIFICA SEMPRE TUTTE QUESTE LABELS!
    labelRecoTrackFast = cms.InputTag("generalTracks"),
    labelL1MuonFast = cms.InputTag("l1extraParticles"),
#    labelL1MuonFast = cms.InputTag("MISSING"),
    labelL2MuonFast = cms.InputTag("hltL2Muons"),
    labelL3MuonFast = cms.InputTag("hltL3Muons"),
#    labelL3MuonFast = cms.InputTag("MISSING"),
    labelGlobalMuonFast = cms.InputTag("muons"),
    labelGlobalMuonTrackFast = cms.InputTag("globalMuons"),
    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons"),
#    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons","UpdatedAtVtx"),

)

process.evtContent = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(process.TrackerMuonAnalyzer)

## To write out events (not need: FastSimulation _is_ fast!)
#process.o1 = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("NewMuonAnalyzerOutputFile.root"),
##    outputCommands = cms.untracked.vstring("keep *",
##                                           "drop *_mix_*_*",
##                                          "drop *_*_*_RECO",
##                                           "drop *_*_*_HLT")
#    )
#process.outpath = cms.EndPath(process.o1)

# Keep the logging output to a nice level #
# process.Timing =  cms.Service("Timing")
process.load("FWCore/MessageService/MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'ERROR'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.destinations = cms.untracked.vstring("pyDetailedInfo.txt")

# Make the job crash in case of missing product
process.options = cms.untracked.PSet( Rethrow = cms.untracked.vstring('ProductNotFound') )

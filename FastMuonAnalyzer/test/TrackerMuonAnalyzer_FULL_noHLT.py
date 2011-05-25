import FWCore.ParameterSet.Config as cms

process = cms.Process("PROD")

# The number of events to be processed.
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2) )
    
# Geometry:
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi")
process.load("Geometry.DTGeometryBuilder.dtGeometry_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_3XY_V15::All"
process.GlobalTag.globaltag = "MC_3XY_V25::All"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#     'file:/afs/cern.ch/user/a/aperrott/scratch0/CMSSW_2_2_0_pre1/src/APAnalyzers/RelVal220pre1SingleMuPt100.root'

#
#-- Single muon, pt = 10 GeV
#
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0017/C8C47517-8F52-DF11-83C9-00304866C398.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0016/C6216236-6A52-DF11-8F90-001A92810AD0.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0016/BA0A2D10-6B52-DF11-B4EE-003048678B12.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0016/8A77F08F-6852-DF11-9352-003048679188.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt10/GEN-SIM-RECO/MC_3XY_V26-v1/0016/6853B313-6852-DF11-A1FA-001A92811708.root'

#
#-- Single muon, pt = 100 GeV
#
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-RECO/MC_3XY_V26-v1/0017/E8728F29-9052-DF11-A46F-0030486791BA.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt100/GEN-SIM-RECO/MC_3XY_V26-v1/0016/B6CBC3F4-3C52-DF11-BF95-0026189437FD.root'
        
#
#-- Single muon, pt = 1000 GeV
#
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt1000/GEN-SIM-RECO/MC_3XY_V26-v1/0017/1260FEBD-8E52-DF11-8EE0-003048678A7E.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt1000/GEN-SIM-RECO/MC_3XY_V26-v1/0016/84DFCED7-4E52-DF11-90CC-0018F3D09692.root',
#    '/store/relval/CMSSW_3_5_8/RelValSingleMuPt1000/GEN-SIM-RECO/MC_3XY_V26-v1/0016/100ED752-4D52-DF11-A2B8-002618943907.root'

#        '/store/relval/CMSSW_3_10_0_pre7/RelValSingleMuPt1000/GEN-SIM-RECO/MC_310_V1-v1/0103/9EBB6657-45FD-DF11-9E3C-002618943918.root',
#        '/store/relval/CMSSW_3_10_0_pre7/RelValSingleMuPt1000/GEN-SIM-RECO/MC_310_V1-v1/0101/F2216E7E-E9FC-DF11-93F4-003048678B0E.root',
#        '/store/relval/CMSSW_3_10_0_pre7/RelValSingleMuPt1000/GEN-SIM-RECO/MC_310_V1-v1/0101/64190210-E7FC-DF11-8F95-002618943914.root'


        '/store/relval/CMSSW_4_1_2/RelValSingleMuPt1000/GEN-SIM-RECO/MC_311_V2-v1/0018/5268ABFF-F044-E011-A5B3-00304867903E.root',
        '/store/relval/CMSSW_4_1_2/RelValSingleMuPt1000/GEN-SIM-RECO/MC_311_V2-v1/0021/2CE84A88-9145-E011-9525-0018F3D095FA.root'#,
       # '/store/relval/CMSSW_4_2_0_pre8/RelValSingleMuPt1000/GEN-SIM-RECO/MC_42_V7-v1/0043/44D6DB7D-D556-E011-8706-00261894388F.root'


    
    )
)
#######################################################################################################################################33
process.TFileService = cms.Service("TFileService",
 
    fileName = cms.string('TrackerMuonAnalyzer_FULL_RelValSingleMuPt1000_May03.root')
)



#############################################################################################################
process.TrackerMuonAnalyzer = cms.EDAnalyzer("TrackerMuonAnalyzer",

#    rootFileName = cms.untracked.string('TrackerMuonAnalyzer_FULL_RelValSingleMuPt1000_May03.root'),

    FastSim = cms.untracked.bool(False),
    IgnoreMissingCollections = cms.untracked.bool(True),
    #RootTree = cms.untracked.bool(True),
    #Debug = cms.untracked.bool(False),
    RootTree = cms.untracked.bool(False),
    Debug = cms.untracked.bool(False),
    useGen = cms.untracked.bool(False),
    useReco = cms.untracked.bool(False),
    useEff = cms.untracked.bool(False),
    test_Track= cms.untracked.bool(True),

    labelSimVertex = cms.InputTag("g4SimHits"),
    labelSimTrack = cms.InputTag("g4SimHits"),
    labelSimHitsTOB = cms.InputTag("g4SimHits","TrackerHits"),
    labelSimHitsTIB = cms.InputTag("g4SimHits","TrackerHits"),

    dtSimHits = cms.InputTag("MISSING"),
    cscSimHits = cms.InputTag("MISSING"),
    rpcSimHits = cms.InputTag("MISSING"),

    # VERIFICA SEMPRE TUTTE QUESTE LABELS!

#    labelRecoTrackFast = cms.InputTag("MISSING"),
    labelRecoTrackFast = cms.InputTag("generalTracks"),
#    labelL1MuonFast = cms.InputTag("MISSING"),
    labelL1MuonFast = cms.InputTag("hltL1extraParticles"),
    labelL2MuonFast = cms.InputTag("MISSING"),
#    labelL2MuonFast = cms.InputTag("hltL2Muons"),
    labelL3MuonFast = cms.InputTag("MISSING"),
#    labelL3MuonFast = cms.InputTag("hltL3Muons"),
#    labelGlobalMuonFast = cms.InputTag("MISSING"),
    labelGlobalMuonFast = cms.InputTag("muons"),
#    labelGlobalMuonTrackFast = cms.InputTag("MISSING"),
    labelGlobalMuonTrackFast = cms.InputTag("globalMuons"),
#    labelStandAloneMuonTrackFast = cms.InputTag("MISSING"),
    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons")
#    labelStandAloneMuonTrackFast = cms.InputTag("standAloneMuons","UpdatedAtVtx"),

)

process.evtContent = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(process.TrackerMuonAnalyzer)

## To write out events (not need: FastSimulation _is_ fast!)
#process.o1 = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("NewMuonAnalyzerOutputFile.root")
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

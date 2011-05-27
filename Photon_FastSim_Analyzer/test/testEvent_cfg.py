import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")




process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.RandomServiceInitialization_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator7TeV_cfi')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
#process.load('FastSimulation.Configuration.HLT_8E29_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
#process.load('FastSimulation.Configuration.Validation_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FastSimulation.Configuration.CommonInputs_cff')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.Configuration.Geometries_START_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')






SkipEvent = cms.untracked.vstring('ProductNotFound')
process.source = cms.Source("EmptySource")
# Generate 20 events

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
# input = cms.untracked.int32(1000)

)
 
# Include the RandomNumberGeneratorService definition
#process.load("FastSimulation.Configuration.RandomServiceInitialization_cff")

# For histograms
#process.load("DQMServices.Core.DQM_cfg")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = "MC_3XY_V21::All"
process.GlobalTag.globaltag =  'START42_V9::All'

# Additional output definition
#TRACKER
process.famosSimHits.MaterialEffects.MuonBremsstrahlung = True

#ECAL
process.famosSimHits.MaterialEffectsForMuonsInECAL.MuonBremsstrahlung = True
#ECAL
process.famosSimHits.MaterialEffectsForMuonsInECAL.EnergyLoss = True
#HCAL
process.famosSimHits.MaterialEffectsForMuonsInHCAL.MuonBremsstrahlung = True
#HCAL
process.famosSimHits.MaterialEffectsForMuonsInHCAL.EnergyLoss = True 



process.famosPileUp.PileUpSimulator = process.PileUpSimulatorBlock.PileUpSimulator
process.famosPileUp.PileUpSimulator.averageNumber = 0
process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.famosSimHits.SimulateMuons = True
process.famosSimHits.TrackerSimHits.firstLoop= False
process.simulation = cms.Sequence(process.simulationWithFamos)
#process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)

# Set correct vertex smearing
process.Realistic7TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic7TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic7TeVCollisionVtxSmearingParameters



process.VtxSmeared = cms.EDProducer("GaussEvtVtxGenerator",
    MeanX = cms.double(0.0),
    MeanY = cms.double(0.0),
    MeanZ = cms.double(0.0),
    SigmaY = cms.double(0.0015),
    SigmaX = cms.double(0.0015),
    SigmaZ = cms.double(5.3),
    TimeOffset = cms.double(0.0),
    src = cms.InputTag("generator")
)

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(1000.01),
        MinPt = cms.double(999.99),
        PartID = cms.vint32(-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single mu pt 1000'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)


process.ProductionFilterSequence = cms.Sequence(process.generator)

process.test = cms.EDAnalyzer("testEvent",
    # necessary to access true particles 
    ParticleFilter = cms.PSet(
        # All protons with an energy larger with EProton (Gev) are kept
        EProton = cms.double(6000.0),
        # Particles with |eta| > etaMax (momentum direction at primary vertex) 
        # are not simulated 
        etaMax = cms.double(5.0),
        # Charged particles with pT < pTMin (GeV/c) are not simulated
        pTMin = cms.double(0.0),
        # Particles with energy smaller than EMin (GeV) are not simulated
        EMin = cms.double(0.0)
    ),
    GeantInfo = cms.bool(False),
    theDebug = cms.bool(False),
   # TrackerRadius       = cms.double(115.0) ## in cm  (TRACKER)
#    TrackerRadius       = cms.double(2000.0) ## in cm  
#     MinRadius =  cms.double(115.0), ## in cm - ref:http://cmsdoc.cern.ch/cms/TDR/TRACKER/final/chapters/tdr_ch1.pdf(Page-4)
      MinRadius =  cms.double(240.0), ## in cm //test2  

    # MaxRadius =  cms.double(160.0) ## in cm 


     #rootFileName = cms.untracked.string('file:/tmp/sfonseca/FastSimHitAnalyzerBrem_200GeV.root')                        # rootFileName = cms.untracked.string('file:/raid/sfonseca/FastSimHitAnalyzerBrem_200GeV.root')
#     rootFileName = cms.untracked.string('file:FastSimHitAnalyzerBrem_500GeV.root')
     # rootFileName = cms.untracked.string('file:/tmp/sfonseca/FastSimHitAnalyzerBrem_800GeV.root')
      #rootFileName = cms.untracked.string('file:/tmp/sfonseca/FastSimHitAnalyzerBrem_1000GeV.root')               
)

process.TFileService = cms.Service("TFileService",
 # fileName = cms.string('FastSimHitAnalyzerBrem_200GeV.root')
  #fileName = cms.string('FastSimHitAnalyzerBrem_500GeV.root')
#fileName = cms.string('FastSimHitAnalyzerBrem_800GeV.root')
#fileName = cms.string('FastSimHitAnalyzerBrem_1000GeV_MuonBremFalse_withoutphoton_selection.root')
 fileName = cms.string('FastSimHitAnalyzerBrem_1000GeV_MuonBrem.root')

)

# Path and EndPath definitions
process.generation_step = cms.Path(cms.SequencePlaceholder("randomEngineStateProducer"))
process.reconstruction = cms.Path(process.famosWithCaloHits+process.test)
#process.validation_step = cms.EndPath(process.validation_prod)
process.endjob_step = cms.EndPath(process.endOfProcess)
#process.out_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step)
#process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.reconstruction,process.validation_step,process.endjob_step])#,process.out_step])
process.schedule.extend([process.reconstruction,process.endjob_step])#,process.out_step])

# special treatment in case of production filter sequence  
for path in process.paths: 
    getattr(process,path)._seq = process.ProductionFilterSequence*getattr(process,path)._seq

# Keep the logging output to a nice level #
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['brem_detailedInfo_1000GeV.txt']

#process.p = cms.Path(
#    cms.SequencePlaceholder("randomEngineStateProducer")+
#    process.offlineBeamSpot+
#    process.famosPileUp+
#    process.famosSimHits+
#    process.test
#)


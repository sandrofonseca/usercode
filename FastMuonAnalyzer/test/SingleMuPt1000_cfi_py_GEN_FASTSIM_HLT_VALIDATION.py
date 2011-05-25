# Auto generated configuration file
# using: 
# Revision: 1.293 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: SingleMuPt1000_cfi.py -s GEN,FASTSIM,HLT:GRun,VALIDATION --pileup=NoPileUp --geometry DB --conditions=auto:startup --eventcontent=FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RECO -n 1000000
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FastSimulation.Configuration.EventContent_cff')
process.load('FastSimulation.PileUpProducer.PileUpSimulator_NoPileUp_cff')
process.load('FastSimulation.Configuration.Geometries_START_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('FastSimulation.Configuration.FamosSequences_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedParameters_cfi')
process.load('FastSimulation.Configuration.HLT_GRun_cff')
process.load('FastSimulation.Configuration.Validation_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.293 $'),
    annotation = cms.untracked.string('SingleMuPt1000_cfi.py nevts:1000000'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_default_perrotta_test.root'),
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_default_TrueMuonBrem.root'),
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case1.root'),#
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case2.root'),#
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case3.root'),#
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case4.root'),#
#    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case5.root'),#
##    fileName = cms.untracked.string('SingleMuPt1000_cfi_py_GEN_FASTSIM_HLT_VALIDATION_Case6.root'),#

    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RECO')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
process.famosSimHits.MaterialEffects.MuonBremsstrahlung = True #Case 1, 2 and 3
##process.famosSimHits.MaterialEffects.MuonBremsstrahlung = False #Case 4,5 and 6
# Smallest bremsstrahlung energy fraction (wrt to the electron energy)
#process.famosSimHits.MaterialEffects.bremEnergyFraction = cms.double(0.005) #Default Case 

#process.famosSimHits.MaterialEffects.bremEnergyFraction = cms.double(0.0005) #Case_1 and Case 4 
#process.famosSimHits.MaterialEffects.bremEnergyFraction = cms.double(0.05)  #Case_2 and case 5
#process.famosSimHits.MaterialEffects.bremEnergyFraction = cms.double(0.5)   #Case_3 and Case 6

process.famosSimHits.SimulateCalorimetry = True
process.famosSimHits.SimulateTracking = True
process.simulation = cms.Sequence(process.simulationWithFamos)
process.HLTEndSequence = cms.Sequence(process.reconstructionWithFamos)
process.Realistic7TeVCollisionVtxSmearingParameters.type = cms.string("BetaFunc")
process.famosSimHits.VertexGenerator = process.Realistic7TeVCollisionVtxSmearingParameters
process.famosPileUp.VertexGenerator = process.Realistic7TeVCollisionVtxSmearingParameters
process.mix.playback = True
process.GlobalTag.globaltag = 'START311_V2::All'

process.generator = cms.EDProducer("FlatRandomPtGunProducer",
    PGunParameters = cms.PSet(
        MaxPt = cms.double(1000.01),
        MinPt = cms.double(999.99),
        PartID = cms.vint32(-13),
        MaxEta = cms.double(2.5),
        MaxPhi = cms.double(3.14159265359),
        MinEta = cms.double(-2.5),
        MinPhi = cms.double(-3.14159265359)
    ),
    Verbosity = cms.untracked.int32(0),
    psethack = cms.string('single mu pt 1000'),
    AddAntiParticle = cms.bool(True),
    firstRun = cms.untracked.uint32(1)
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen_genonly)
process.reconstruction = cms.Path(process.reconstructionWithFamos)
process.prevalidation_step = cms.Path(process.prevalidation)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.validation_step = cms.EndPath(process.validation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.reconstruction,process.prevalidation_step,process.validation_step,process.endjob_step,process.FEVTDEBUGHLToutput_step])
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

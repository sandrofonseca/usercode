import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.GeometryExtended_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("SimG4Core.Application.g4SimHits_cfi")
process.load("IOMC.RandomEngine.IOMC_cff")


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('testNutple.root')
)



process.source = cms.Source("EmptySource")

process.generator = cms.EDProducer("FlatRandomEGunProducer",
    PGunParameters = cms.PSet(
        PartID = cms.vint32(13),
        #PartID = cms.vint32(11),
        MaxEta = cms.double(0.1),
        MaxPhi = cms.double(0.1),
        MinEta = cms.double(0.1),
        #MinE = cms.double(200.0),
         MinE = cms.double(500.0),
        # MinE = cms.double(800.0),
         #MinE = cms.double(1000.0),
        MinPhi = cms.double(0.1),
        #MaxE = cms.double(200.0)
        MaxE = cms.double(500.0)
        # MaxE = cms.double(800.0)
         # MaxE = cms.double(1000.0)
    ),
    AddAntiParticle = cms.bool(False),
    Verbosity = cms.untracked.int32(0)
)

process.VtxSmeared = cms.EDProducer("GaussEvtVtxGenerator",
    src   = cms.InputTag("generator"),
    MeanX = cms.double(0.0),
    MeanY = cms.double(0.0),
    MeanZ = cms.double(0.0),
    SigmaX = cms.double(0.0015),
    SigmaY = cms.double(0.0015),
    SigmaZ = cms.double(5.3),
    TimeOffset = cms.double(0.0)
)

process.g4SimHits.UseMagneticField = False
process.g4SimHits.Physics.type = 'SimG4Core/Physics/QGSP'
process.g4SimHits.Physics.DefaultCutValue = 0.07
process.g4SimHits.StackingAction.SaveFirstLevelSecondary = True

process.g4SimHits.Watchers = cms.VPSet(
    cms.PSet(
      ElectronBremAnalysis = cms.PSet(
        Verbosity           = cms.int32(3),
        TrackerRadius       = cms.double(1150.0), ## in mm
        photonCut           = cms.double(100.0) ## in MeV (std Geant4 Unit)
      ),
    type = cms.string('ElectronBremAnalysis')
    )
)


#process.o1 = cms.OutputModule("PoolOutputModule",
#       fileName = cms.untracked.string('file:simg4.root')
 #    fileName = cms.untracked.string('file:/raid/sfonseca/sim_g4_em.root')
#    fileName = cms.untracked.string('file:/raid/sfonseca/sim_g4_em_800GeV.root')
#     fileName = cms.untracked.string('file:/raid/sfonseca/sim_g4_em_1000GeV.root')

#)

process.p1 = cms.Path(process.generator*process.VtxSmeared*process.g4SimHits)
#process.outpath = cms.EndPath(process.o1)



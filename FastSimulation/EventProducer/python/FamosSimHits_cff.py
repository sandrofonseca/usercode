import FWCore.ParameterSet.Config as cms

# Here so that python translator can see the names
# Now beta function vertex smearing 
#from FastSimulation.Event.Early10TeVCollisionVertexGenerator_cfi import *
from FastSimulation.Event.Realistic7TeV2011CollisionVertexGenerator_cfi import *
from FastSimulation.Event.ParticleFilter_cfi import *
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from FastSimulation.MaterialEffects.MaterialEffects_cfi import *
from FastSimulation.TrajectoryManager.ActivateDecays_cfi import *
from FastSimulation.TrajectoryManager.TrackerSimHits_cfi import *
from FastSimulation.Calorimetry.Calorimetry_cff import *
famosSimHits = cms.EDProducer("FamosProducer",
    # FastCalorimetry
    FamosCalorimetryBlock,
    # Conditions to save Tracker SimHits 
    TrackerSimHitsBlock,
    # Material effects to be simulated in the tracker material and associated cuts
    MaterialEffectsBlock,
    # Material effects for muons in ECAL
    MaterialEffectsForMuonsInECALBlock,
    # Material effects for muons in HCAL
    MaterialEffectsForMuonsInHCALBlock,
    # Services
     MuonServiceProxy,
    # (De)activate decays of unstable particles (K0S, etc...)
    ActivateDecaysBlock,
    # Kinematic cuts for the particle filter in the SimEvent
    ParticleFilterBlock,
    # The HepMCProduct source
    SourceLabel = cms.InputTag("generator"),
    # The genParticle source (in case there is no HepMCProduct)
    GenParticleLabel = cms.InputTag("genParticles"),
    # If false, no SimTrack collection for Muons is stored
    SimulateMuons = cms.bool(True),
    # The beam spot source
    BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    # Run Number
    RunNumber = cms.untracked.int32(1001),
    Verbosity = cms.untracked.int32(0),
    # If false, no misalignment can be applied in the tracker
    ApplyAlignment = cms.bool(False),
    # If false, no SimHits are simulated in the tracker
    # (but the tracker material is still here)
    SimulateTracking = cms.bool(True),
    # If false, no PCaloHits are produced
    SimulateCalorimetry = cms.bool(True),
    # If the following is false, no B field map is used, 
    # and B is set to 4 T altogether
    UseMagneticField = cms.bool(True),
    # The primary vertex smearing 
    VertexGenerator = cms.PSet(
        myVertexGenerator
    )
)



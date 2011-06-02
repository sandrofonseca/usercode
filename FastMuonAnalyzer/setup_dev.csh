#!/bin/tcsh

#pushd $CMSSW_BASE/src
# New improvements for HCAL and Muon Brem
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/Calorimetry/src/CalorimetryManager.cc
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/Calorimetry/interface/CalorimetryManager.h

cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MaterialEffects.cc
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MuonBremsstrahlungSimulator.cc
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MaterialEffects.h
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MuonBremsstrahlungSimulator.h
cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/python/MaterialEffects_cfi.py


pushd $CMSSW_BASE/src

#Testing in CMSSW_4_2_0
cvs co  Configuration/Generator  
cvs co -r V02-03-19  FastSimulation/Calorimetry                       
cvs co -r V01-22-02  FastSimulation/Configuration                     
cvs co -r V04-05-09  FastSimulation/MaterialEffects                   
cvs co -r V00-03-11  FastSimulation/MuonSimHitProducer                
cvs co -r V01-04-22  FastSimulation/TrajectoryManager  

cp UserCode/Calo_May_30/FastSimulation/Calorimetry/src/CalorimetryManager.cc FastSimulation/Calorimetry/src
cp UserCode/Calo_May_30/FastSimulation/Calorimetry/interface/CalorimetryManager.h FastSimulation/Calorimetry/interface

cp UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MaterialEffects.cc FastSimulation/MaterialEffects/src
cp UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MuonBremsstrahlungSimulator.cc FastSimulation/MaterialEffects/src
cp UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MaterialEffects.h FastSimulation/MaterialEffects/interface
cp UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MuonBremsstrahlungSimulator.h FastSimulation/MaterialEffects/interface
cp UserCode/Calo_May_30/FastSimulation/MaterialEffects/python/MaterialEffects_cfi.py FastSimulation/MaterialEffects/python


#Analysis code for Tracker

cvs co -r V00-00-02 -d FastMuonBremAnalyzers/FastMuonAnalyzer  UserCode/FastMuonBremAnalyzers
cvs co -r V00-00-01 -d FastMuonBremAnalyzers/Photon_FastSim_Analyzer  UserCode/FastMuonBremAnalyzers/Photon_FastSim_Analyzer
# G4 Observer
#cvs co -r V00-00-01  -d  FastMuonBremAnalyzers/ElectronBremAnalysis UserCode/FastMuonBremAnalyzers/ElectronBremAnalysis
#scram build -c

#CPU timing tests
cvs co -d FastMuonBremAnalyzers/timing_cpu_test UserCode/FastMuonBremAnalyzers/timing_cpu_test/timing.cpp 


# New improvements for HCAL and Muon Brem
#cvs co -r V00-04-00 UserCode/Calo_May_30/Calorimetry/src/CalorimetryManager.cc
#cvs co -r V00-04-00 UserCode/Calo_May_30/Calorimetry/interface/CalorimetryManager.h

#cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MaterialEffects.cc
#cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/src/MuonBremsstrahlungSimulator.cc
#cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MaterialEffects.h
#cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/interface/MuonBremsstrahlungSimulator.h
#cvs co -r V00-04-00 UserCode/Calo_May_30/FastSimulation/MaterialEffects/python/MaterialEffects_cfi.py



scramv1 b -j 4
popd

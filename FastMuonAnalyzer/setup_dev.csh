#!/bin/tcsh

pushd $CMSSW_BASE/src
#Testing in CMSSW_4_2_0
cvs co  Configuration/Generator  
cvs co -r V02-03-19  FastSimulation/Calorimetry                       
cvs co -r V01-22-02  FastSimulation/Configuration                     
cvs co -r V04-05-09  FastSimulation/MaterialEffects                   
cvs co -r V00-03-11  FastSimulation/MuonSimHitProducer                
cvs co -r V01-04-22  FastSimulation/TrajectoryManager  

#Analysis code for Tracker
cvs co -d FastMuonBremAnalyzers/FastMuonAnalyzer  UserCode/FastMuonBremAnalyzers

scramv1 b -j 4
popd

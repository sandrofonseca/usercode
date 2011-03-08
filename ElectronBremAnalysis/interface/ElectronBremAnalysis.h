/*
 *  MuonBremAnalysis.h
 *  
 *
 *  Created by Sandro Fonseca de Souza on 14/06/10.
 *  Copyright 2010 UERJ/CERN. All rights reserved.
 *
 */

#ifndef ElectronBremAnalysis_h
#define ElectronBremAnalysis_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Notification/interface/Observer.h"



#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TUnixSystem.h"
#include "TSystem.h"
#include "TMath.h"
#include "TF1.h"

#include <string>
#include <iostream>
#include <memory>
#include <vector>
//
// class decleration
//

class BeginOfRun;
class BeginOfEvent;
class EndOfEvent;
class EndOfRun;
class BeginOfTrack;
class G4Step;


class ElectronBremAnalysis : public SimWatcher,
public Observer<const BeginOfRun*>,
public Observer<const EndOfRun*>,
public Observer<const BeginOfEvent*>,
public Observer<const BeginOfTrack*>,
public Observer<const EndOfEvent*>,
public Observer<const G4Step*> {
	
public:
	ElectronBremAnalysis(const edm::ParameterSet& p);
	virtual ~ElectronBremAnalysis();
	
private:
	
     ElectronBremAnalysis(const ElectronBremAnalysis&); // stop default
     const ElectronBremAnalysis& operator=(const ElectronBremAnalysis&);
	
	void update(const BeginOfRun* );
	void update(const BeginOfEvent* );
	void update(const EndOfEvent*  );
        void update(const EndOfRun*  );
        void update(const G4Step* step );
	void update(const BeginOfTrack* trk );
        int verbosity;
        double  track_radius;
        double photon_cut; 

private:
       const static int i=1000000;
       std::string PName;
       int pdg[i],Nphotons[i],iev,pdgNeutral,motherNeutral;
       double Kenergy[i],TotalEnergyDep[i],px[i],py[i],R[i],phi[i],eta[i],mother[i],EnergyNeutral,phiNeutral,etaNeutral;
       std::vector<int>       *Nphotons_;
       std::vector<double>    *Kenergy_;
       std::vector<int>       *pdg_;
       std::vector<double>    *TotalEnergyDep_;
       std::vector<double>    *phi_;
       std::vector<double>    *eta_;
       std::vector<int>       *Nphotons_evt;
       std::vector<double>    *Kenergy_evt;
       std::vector<int>       *Niev_;
       std::vector<double>    *px_;
       std::vector<double>    *py_;
       std::vector<double>    *R_;
	

        TFile    *file;
        TTree    *t1;
 

		
};
#endif 


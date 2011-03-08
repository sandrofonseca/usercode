// -*- C++ -*-
//
// Package:    ElectronBremAnalysis
// Class:      ElectronBremAnalysis
// 
/**
 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sandro Fonseca De Souza,32 4-C14,+41227674949,
//         Created:  Fri Jun 18 13:56:05 CEST 2010
// $Id$
//
//


 // system include files
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
//#include <stdio.h>
//#include <math.h>


// user include files

#include "SimG4CMS/ElectronBremAnalysis/interface/ElectronBremAnalysis.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <TH1F.h>
#include <TH1I.h>

//
#include "SimG4Core/Notification/interface/BeginOfRun.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "G4Run.hh"
#include "G4RunManager.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"

////////////////////////////////////////////////////////////////////
ElectronBremAnalysis::ElectronBremAnalysis(const edm::ParameterSet& p) 
  :Nphotons_(0), Kenergy_(0),pdg_(0),EnergyNeutral(0),TotalEnergyDep_(0),phi_(0),eta_(0),Nphotons_evt(0), Kenergy_evt(0),Niev_(0),px_(0),py_(0)
{
//http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/SimG4Core/CheckSecondary/src/CheckSecondary.cc?revision=1.5&view=markup
    edm::ParameterSet m_p = p.getParameter<edm::ParameterSet>("ElectronBremAnalysis"); 
    verbosity                = m_p.getParameter<int>("Verbosity");
    //PName               = m_p.getParameter<std::string>("ProcessName");
    track_radius             = m_p.getParameter<double>("TrackerRadius"); 
    photon_cut               =  m_p.getParameter<double>("photonCut");
  //  std::string saveFile = m_p.getUntrackedParameter<std::string>("SaveInFile", "None");

  if (verbosity > 0){
   std::cout<<std::endl;
   std::cout<<"============================================================================"<<std::endl;
   std::cout << "ElectronBremAnalysis:: Initialized as observer"<< std::endl;
   
   //std::cout <<" Ntuple will be created"<< std::endl;
   //  std::cout <<" Ntuple file: "<<NtFileName<<std::endl;
   }

 // ntuple = new TNtuple("NTMuonBremPhotons","NTMuomBremPhotons","evt:trackid:charge:pdgcode:x:y:z:stepl:stepe:eta:phi:vpx:vpy:vpz");

 
//Building vectors
   Nphotons_               = new std::vector<int>();  
   Kenergy_                = new std::vector<double>();
   pdg_                    = new std::vector<int>();
   TotalEnergyDep_         = new std::vector<double>();
   phi_                    = new std::vector<double>();
   eta_                    = new std::vector<double>();
   Niev_                   = new std::vector<int>();
   px_                     = new std::vector<double>();
   py_                     = new std::vector<double>();
   R_                      = new std::vector<double>();

  Kenergy_evt                = new std::vector<double>();
  Nphotons_evt               = new std::vector<int>();

}
///////////////////////////////////////////////////////////
ElectronBremAnalysis::~ElectronBremAnalysis() {
 //std::cout << "Save the Secondary Tree "
  //         <<" Tree= "<< t1->GetName() <<" File= " << file->GetName() << std::endl;


  //file->cd();
  //t1->Write();
  //file->Close();
  //delete file;
 


    std::cout << std::endl << "End of ElectronBremAnalysis"
	      << std::endl; 

  
 // Deleting vectors
    if (Nphotons_)           delete Nphotons_;
    if (Kenergy_)            delete Kenergy_;
    if (pdg_)                delete pdg_;
    if (TotalEnergyDep_)     delete TotalEnergyDep_;
    if (phi_)                delete phi_;
    if (eta_)                delete eta_;
    if (Niev_)               delete Niev_;
    if (px_)                 delete px_;
    if (py_)                 delete py_;
    if (R_)                  delete R_;

    if (Nphotons_evt)          delete Nphotons_evt;
    if (Kenergy_evt)           delete Kenergy_evt;

}


//================================================================== per RUN
void ElectronBremAnalysis::update(const BeginOfRun * run) {

   edm::Service<TFileService> fs;
 //int Nevents = (*run)()->GetNumberOfEvent();
 //std::cout << "Total_Number_of_Events="<< Nevents << std::endl; 

std::cout << std::endl << "ElectronBremAnalysis: Starting Run"<< std::endl; 
//  if (verbosity > 0 ) { 
    std::cout << "ElectronBremAnalysis: output step root file created"<< std::endl;
    //TString filename = NtFileName;
    //file = new TFile(filename,"RECREATE");
    //file->cd();
    t1 = fs->make<TTree>("T1", "Secondary Particle Information");
    t1->Branch("NumberPhotonsSecondaries", "std::vector<int>",         &Nphotons_);
    t1->Branch("KenergyPhotonsSecondaries", "std::vector<double>",      &Kenergy_);
    t1->Branch("PDGIDSecondaries", "std::vector<int>",         &pdg_);
    t1->Branch("phi","std::vector<double>" ,&phi_);
    t1->Branch("eta","std::vector<double>" ,&eta_);
    t1->Branch("NbOfEvents",&iev ,"iev/I");
    t1->Branch("px","std::vector<double>" ,&px_);
    t1->Branch("py","std::vector<double>" ,&py_);
    t1->Branch("R","std::vector<double>" ,&R_);

 iev = 0;
 }


/////////////////////////////////////////////////////////////
void ElectronBremAnalysis::update(const BeginOfEvent* g4Event){
iev = (*g4Event)()->GetEventID();
 if (verbosity > 0){
 std::cout << "MuonBremAnalysis::=====> Begin of event = " 
	   << iev << std::endl;
     }
 //iev++;
 pdg[iev]= 0;
 Kenergy[iev]= 0.0;
 TotalEnergyDep[iev] = 0.0;
 Nphotons[iev]= 0;
 phi[iev] = 0.0;
 eta[iev] = 0.0;
 EnergyNeutral = phiNeutral = etaNeutral = 0.0;
 motherNeutral = pdgNeutral = 0;
 px[iev] = 0.0;
 py[iev] = 0.0;
 R[iev] = 0.0;

 //Clearing vectors
 (*Nphotons_).clear();
 (*Kenergy_).clear();
 (*pdg_).clear();
 (*TotalEnergyDep_).clear();
 (*phi_).clear();
 (*eta_).clear();
 (*Niev_).clear();

 (*Nphotons_evt).clear();
 (*Kenergy_evt).clear();
 (*px_).clear();
 (*py_).clear();
 (*R_).clear();


}
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
void ElectronBremAnalysis::update(const G4Step* theStep){ }
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
void ElectronBremAnalysis::update(const BeginOfTrack * trk){
const G4Track * theTrack    = (*trk)();
if(theTrack == NULL) return;

 Kenergy[iev]= 0.0;

pdg[iev] = theTrack->GetDefinition()->GetPDGEncoding();
Kenergy[iev] = theTrack->GetKineticEnergy();
phi[iev] = theTrack->GetMomentum().phi();
eta[iev] = theTrack->GetMomentum().eta();
mother[iev] = theTrack->GetParentID();
px[iev] = theTrack->GetPosition().x();
py[iev] = theTrack->GetPosition().y();
R[iev]= sqrt((px[iev]*px[iev]) + (py[iev]*py[iev]));
 
//Tracker radius 1150.0 mm (115.0 cm)
//Muon brem: process name : muBrens
//Electron brem: process name: eBrem
 
//PName ="eBrem";

TrackInformation * trkInfo = (TrackInformation *)(theTrack->GetUserInformation());
  
  if (trkInfo) {
	  if (!(trkInfo->isPrimary()) && pdg[iev] == 22 && theTrack->GetCreatorProcess()->GetProcessName()=="muBrems" && Kenergy[iev] > 0.0 && R[iev] < track_radius ) {
           if(Kenergy[iev] > photon_cut ){  
     if (verbosity > 2){
      
               
      std::cout << "Event_ID= "<< iev << std::endl;
    // particle is secondary
      std::cout << " Particle_Secondary_PDG_Encoding= "<< pdg[iev] << "Process Name:" << theTrack->GetCreatorProcess()->GetProcessName()<< std::endl;
      //kinetici energy
      //std::cout << " Kenergy_of Secondary_Particle= "<< Kenergy[iev] << "eV"  << std::endl;
      std::cout << " Kenergy_of Secondary_Particle= "<< G4BestUnit(Kenergy[iev],"Energy") <<"Photon_cut (in Mev) >  " << photon_cut << std::endl;
     // Phi
       std::cout << " Phi= "<< G4BestUnit(phi[iev],"Angle") << std::endl;
     //eta
       std::cout << " Eta= "<< eta[iev] << "ParentID"<< theTrack->GetParentID() <<std::endl;
      //px ,py and R
       std::cout << " px= "<< px[iev] <<" py= " << py[iev] <<" R= " << G4BestUnit(R[iev],"Length") << std::endl;    
      
       //std::cout <<"R_u= "<< R[iev] << std::endl;


   //volume name 
       std::cout << " Volume: " << theTrack->GetVolume()->GetName() << std::endl;
   //process name 
       // std::cout << " ProcessName: " << theTrack->GetStep()->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  << std::endl;
  // #ifdef DebugLog
  if ( theTrack->GetNextVolume() != 0 )
      std::cout << " NextVolume: " << theTrack->GetNextVolume()->GetName() << std::endl;
  else 
      std::cout << " NextVolume: OutOfWorld" << std::endl;
//    #endif

           }
         


   Nphotons[iev]++;
   EnergyNeutral += Kenergy[iev];
//   pdgNeutral = pdg[iev];
//   phiNeutral = phi[iev];
 //  etaNeutral = eta[iev];
  // motherNeutral = mother[iev];

   (*Nphotons_).push_back(Nphotons[iev]);
   (*Kenergy_).push_back(EnergyNeutral);
   (*pdg_).push_back(pdg[iev]);
   (*phi_).push_back(phi[iev]);
   (*eta_).push_back(eta[iev]);
   (*px_).push_back(px[iev]);
   (*py_).push_back(py[iev]);
   (*R_).push_back(R[iev]);
   

	   }

	  }
       
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ElectronBremAnalysis::update(const EndOfEvent* g4Event){

if (Nphotons[iev]>0){
std::cout << "Number_of_Photons= "<< Nphotons[iev]
                                 // <<" PDG_ID= "<< pdgNeutral
                                  //<<" Kinetic Energy= "<< EnergyNeutral << " eV" 
                                  <<" Kinetic Energy= "<< G4BestUnit(EnergyNeutral,"Energy")<< std::endl;
                                 // << "Phi= "<< G4BestUnit(phiNeutral,"Angle")
                                  //<< "Mother= "<< motherNeutral  
                                  //<< "Eta= "<< etaNeutral << std::endl;
                                  
//(*Nphotons_).push_back(Nphotons[iev]);
//(*Kenergy_).push_back(EnergyNeutral);

}
else{

std::cout << "There are not photons in this event = "<< (*g4Event)()->GetEventID() <<std::endl;

}
 
  std::cout << "Photons_Secondary::EndofEvent =====> Event " 
			     << (*g4Event)()->GetEventID() << std::endl;
 //Filling tree 
 t1->Fill();
 

 
}
////////////////////////////////////////////////////////////
void ElectronBremAnalysis::update(const EndOfRun* run){
std::cout << "====================================::EndofRun ==========================================================="<< std::endl;
 

}


/////////////////////////////////////////////////////////////////
//define this as a plug-in
//DEFINE_FWK_MODULE(ElectronBremAnalysis);

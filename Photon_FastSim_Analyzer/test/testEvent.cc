// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/Event/interface/FSimVertex.h"
#include "FastSimulation/Particle/interface/ParticleTable.h"



#include "FWCore/ServiceRegistry/interface/Service.h"
#include <vector>
#include <string>
#include <cmath>

//////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TNtuple.h"

//namespaces
using namespace std;
using namespace edm;

class testEvent : public edm::EDAnalyzer {
public :
  explicit testEvent(const edm::ParameterSet&);
  ~testEvent();

  virtual void analyze(const edm::Event&, const edm::EventSetup& );
  virtual void beginRun(edm::Run const&, edm::EventSetup const&  );
private:
 
   TFile    *file;
   TTree    *t1;
  const static int i=1000000;
  
  bool isGeant,Debug;
  double max_radius, min_radius;
  edm::ParameterSet particleFilter_;
  std::vector<FSimEvent*> mySimEvent;
  int pdg[i],Nphotons[i],iev,pdgNeutral,motherNeutral;
  double Kenergy[i],TotalEnergyDep[i],px[i],py[i],R[i],phi[i],eta[i],mother[i],EnergyNeutral,phiNeutral,etaNeutral;
  double frac[i],muon_momentum[i],muon_momentum_out[i],deltaPoutPgenMuon[i],deltaPoutPgenInvPgenMuon[i],muon_energy[i];
  double Pin[i],Pout[i];

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
       std::vector<double>    *frac_;
       std::vector<double>    *muon_momentum_;
       std::vector<double>    *muon_momentum_out_;      
       std::vector<double>    *deltaPoutPgenMuon_ ;
       std::vector<double>    *deltaPoutPgenInvPgenMuon_;
       std::vector<double>    *muon_energy_;
       std::vector<double>    *Pin_;
       std::vector<double>    *Pout_;
               
 
//counters
  int numberOfEvents;
  int nevts;
 
 
};
///////////////////////////////////////////////////////////
testEvent::testEvent(const edm::ParameterSet& p) :
  isGeant(true),Debug(true),mySimEvent(2, static_cast<FSimEvent*>(0)),
  Kenergy_(0),Nphotons_(0),pdg_(0),phi_(0),eta_(0),px_(0),py_(0) {
  
  particleFilter_ = p.getParameter<edm::ParameterSet> ( "ParticleFilter" );
  isGeant = p.getParameter<bool>("GeantInfo");
  Debug = p.getParameter<bool>("theDebug");
   min_radius = p.getParameter<double>("MinRadius");
  


  // For the full sim
  if ( isGeant) mySimEvent[0] = new FSimEvent(particleFilter_);
  // For the fast sim
  mySimEvent[1] = new FSimEvent(particleFilter_);
  								
//Building vectors
    
   Kenergy_                = new std::vector<double>();
   Nphotons_               = new std::vector<int>();
   pdg_                    = new std::vector<int>();
   TotalEnergyDep_         = new std::vector<double>();
   phi_                    = new std::vector<double>();
   eta_                    = new std::vector<double>();
   Niev_                   = new std::vector<int>();
   px_                     = new std::vector<double>();
   py_                     = new std::vector<double>();
   R_                      = new std::vector<double>();  
   frac_                   = new std::vector<double>();
  
   muon_momentum_          = new std::vector<double>();
   muon_momentum_out_      = new std::vector<double>();
   deltaPoutPgenMuon_       = new std::vector<double>();
   deltaPoutPgenInvPgenMuon_ = new std::vector<double>();
   muon_energy_              = new std::vector<double>();
   Pin_                      = new std::vector<double>();
   Pout_                     = new std::vector<double>();

  

}
////////////////////////////////////////////////////////////////////
testEvent::~testEvent()
{

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
    if (frac_)               delete frac_;
    if (Nphotons_evt)          delete Nphotons_evt;
    if (Kenergy_evt)           delete Kenergy_evt;
   
    if (muon_momentum_)         delete muon_momentum_;
    if (muon_momentum_out_)     delete muon_momentum_out_;

    if (deltaPoutPgenMuon_)       delete deltaPoutPgenMuon_;
    if (deltaPoutPgenInvPgenMuon_)delete deltaPoutPgenInvPgenMuon_;
    if (muon_energy_)             delete muon_energy_;
    if (Pin_)                     delete Pin_;
    if (Pout_)                    delete Pout_;


}
///////////////////////////////////////////////////////////////
void testEvent::beginRun(edm::Run const&, edm::EventSetup const& es)
{ 
//Initialization counters
  numberOfEvents =0;
  nevts=0;
   edm::Service<TFileService> fs;

 

    t1 = fs->make<TTree>("T1", "Secondary Particle Information");
    t1->Branch("NumberPhotonsSecondaries", "std::vector<int>",         &Nphotons_);
    t1->Branch("KenergyPhotonsSecondaries", "std::vector<double>",      &Kenergy_);
    t1->Branch("PDGIDSecondaries", "std::vector<int>",         &pdg_);
    t1->Branch("phi","std::vector<double>" ,&phi_);
    t1->Branch("eta","std::vector<double>" ,&eta_);
    t1->Branch("NbOfEvents",&numberOfEvents ," numberOfEvents/I");
    t1->Branch("px","std::vector<double>" ,&px_);
    t1->Branch("py","std::vector<double>" ,&py_);
    t1->Branch("R","std::vector<double>" ,&R_);
    t1->Branch("frac","std::vector<double>" ,&frac_);
    t1->Branch("muon_momentum","std::vector<double>" ,&muon_momentum_);
    t1->Branch("muon_momentum_out","std::vector<double>" ,&muon_momentum_out_);
    t1->Branch("deltaPoutPgenMuon","std::vector<double>" ,&deltaPoutPgenMuon_);
    t1->Branch("deltaPoutPgenInvPgenMuon","std::vector<double>" ,&deltaPoutPgenInvPgenMuon_);
    t1->Branch("muon_energy","std::vector<double>" ,&muon_energy_);
    t1->Branch("Pin","std::vector<double>" ,&Pin_);
    t1->Branch("Pout","std::vector<double>" ,&Pout_);

 
  // init Particle data table (from Pythia)
  edm::ESHandle < HepPDT::ParticleDataTable > pdt;
  es.getData(pdt);
  if ( !ParticleTable::instance() ) ParticleTable::instance(&(*pdt));
  if ( isGeant ) mySimEvent[0]->initializePdt(&(*pdt));
  mySimEvent[1]->initializePdt(&(*pdt));

}
////////////////////////////////////////////////////////////////////////////////
void
testEvent::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

int nbrems=0;
//int nbremsmin=0;
int  run =  iEvent.id().run();
int  evt =  iEvent.id().event();
  if (Debug) {
    cout << "Run: " << run << " Event: " << evt << endl;
    }
  nevts++;
 numberOfEvents++;

 if ( numberOfEvents%200 == 0 ) cout << " --- Events Processed " << numberOfEvents << endl;

  if ( isGeant ) { 
    edm::Handle<std::vector<SimTrack> > fullSimTracks;
    iEvent.getByLabel("g4SimHits","",fullSimTracks);
    edm::Handle<std::vector<SimVertex> > fullSimVertices;
    iEvent.getByLabel("g4SimHits","",fullSimVertices);
    mySimEvent[0]->fill( *fullSimTracks, *fullSimVertices );
  }

 

  Handle<std::vector<SimTrack> > fastSimTracks;
  iEvent.getByLabel("famosSimHits","",fastSimTracks);
  Handle<std::vector<SimVertex> > fastSimVertices;
  iEvent.getByLabel("famosSimHits","",fastSimVertices);
  mySimEvent[1]->fill( *fastSimTracks, *fastSimVertices );
 
  
    for ( unsigned ievt=1; ievt<3; ++ievt ) {
     if ( isGeant || ievt == 1 ) { 


       std::vector<int> myGammas;
       int ntracks = mySimEvent[ievt]->nTracks();
      
       for (int fsimi=0;fsimi<ntracks;++fsimi){

       FSimTrack& theTrack = mySimEvent[ievt]->track(fsimi);
      
         pdg[iev]= 0;
         Kenergy[iev]= 0.0;
         TotalEnergyDep[iev] = 0.0;
         Nphotons[iev]= 0;
         phi[iev] = 0.0;
         eta[iev] = 0.0;
         px[iev] = 0.0;
         py[iev] = 0.0;
         R[iev] = 0.0; 
         frac[iev] = 0.0;
         muon_momentum[iev] = 0.0;
         muon_momentum_out[iev] = 0.0;
         deltaPoutPgenMuon[iev] = 0.0;
         deltaPoutPgenInvPgenMuon[iev] = 0.0;
         muon_energy[iev] = 0.0;
         Pin[iev] = 0.0;
         Pout[iev] = 0.0;

 //Clearing vectors
 (*Nphotons_).clear();
 (*Kenergy_).clear();
 (*pdg_).clear();
 (*TotalEnergyDep_).clear();
 (*phi_).clear();
 (*eta_).clear();
 (*Niev_).clear();
 (*px_).clear();
 (*py_).clear();
 (*R_).clear();
 (*frac_).clear();
 (*muon_momentum_).clear();
 (*muon_momentum_out_).clear();
 (*deltaPoutPgenMuon_).clear();
 (*deltaPoutPgenInvPgenMuon_).clear();
 (*muon_energy_).clear();
 (*Pin_).clear();
 (*Pout_).clear();

      //select the primary particle
	 if(abs(theTrack.type()) == 13 && theTrack.vertex().noParent()){
	   
           int firstDaughter = -1;
	   int lastDaughter = -1;
	  
	   double feta=fabs(theTrack.momentum().eta());
	   //	   h1[ievt]->Fill(feta);
   
	   if(Debug) cout<< "Primary PDG_ID: "<< theTrack.type() <<" eta: " <<feta << endl;
              // select the muon's daughters
	      if( theTrack.nDaughters() ){ 
		firstDaughter = theTrack.daughters()[0];
                lastDaughter = theTrack.daughters()[theTrack.nDaughters()-1];                                 
                   }
	      //to define of the Muon and your daughters
                XYZTLorentzVector theMuon=theTrack.momentum();

                double Pin = theTrack.momentum().P();
                double Pout = theTrack.trackerSurfaceMomentum().P();

	// 	if(Debug){
//                 cout << " The starting muon " << " P_inTracker(GeV): "<<theMuon.P() << " " 
//  		     << "MuonMomentum_TrackerOut: (GeV) "  << theTrack.trackerSurfaceMomentum().P()<<" " 
// 		     <<" Vertex First Position: "<< theTrack.vertex().position() << " " 
// 		     <<" Vertex Last Position: " << theTrack.endVertex().position() << " "
// 		     <<" N of daughters: "<<  theTrack.nDaughters() << " " 
// 		     <<" Fist daughters: " << firstDaughter << " " 
// 		     << " Last daughers: "<< lastDaughter << " " 
// 		     << endl;
	
// 		}
		//
                    (*Pin_).push_back(Pin);
                    (*Pout_).push_back(Pout);


		    //  daughter photons
	                if (!(firstDaughter<0 ||lastDaughter<0)){
			for(int igamma=firstDaughter;igamma<=lastDaughter;++igamma){
                            FSimTrack myGamma =  mySimEvent[ievt]->track(igamma);
			    // Check photons
			       if(myGamma.type() == 22 ){
                                  int PDGgamma = myGamma.type();                              
                                  double px = myGamma.vertex().position().x();
                                  double py = myGamma.vertex().position().y();
				  double radius = sqrt((px * px) + (py * py));
                                 
                                     if(radius < min_radius){
		cout <<"#######################################"<< endl;
               cout << "Accessing photons brem  using FamosSimHits" << endl;
               cout << "######################################"<<endl;
				 cout << "Event number " << nevts << std::endl;
                               
                                  XYZTLorentzVector theMother = theMuon;
                                  theMuon=theMuon-myGamma.momentum();
                                  nbrems++;
				 
                                  myGammas.push_back(igamma);
                                  double muon_energy = theMuon.e(); 
				  double muon_momentum = theMuon.P();
                                  double muon_momentum_out =  theTrack.trackerSurfaceMomentum().P();
                                  double deltaPoutPgenMuon = abs(muon_momentum_out-muon_momentum);
                                  double deltaPoutPgenInvPgenMuon = deltaPoutPgenMuon/abs(muon_momentum);
                                  double gamma_energy = myGamma.momentum().e();
                                  double gamma_frac_energy = gamma_energy/muon_energy;
                                  //double radius = myGamma.vertex().position().pt();

                                  cout << "gamma_energy: (MeV)"<< gamma_energy << endl;
                                  double eta_photon = myGamma.momentum().eta();
                                  double phi_photon = myGamma.momentum().phi();
 
				 //  if(Debug) {
                                   cout <<" muon_energy: (GeV) "<< muon_energy << " "
					<<" muon_momentum: (GeV) "<< muon_momentum << " "
                                	<<" muon_momentum_out: (GeV) "<< muon_momentum_out << " "
				      	<<" deltaPoutPgenMuon: "<< deltaPoutPgenMuon << " "
                                        <<" deltaPoutPgenInvPgenMuon: "<< deltaPoutPgenInvPgenMuon << " "
                                        <<" gamma_energy: (MeV) "<< gamma_energy << " "
				        <<" gamma_frac_energy: (MeV) "<< gamma_frac_energy <<" "
				        <<" radius : "<< radius <<" "<<" eta_photon : "<< eta_photon <<" "
					<<" px: "<< px <<" " <<" py: "<< py <<" PDG: "<< PDGgamma << endl;  
				   //	    }
				        
                                       (*muon_momentum_).push_back(muon_momentum);
                                       (*muon_momentum_out_).push_back(muon_momentum_out);
                                       (*deltaPoutPgenMuon_).push_back(deltaPoutPgenMuon);
                                       (*deltaPoutPgenInvPgenMuon_).push_back(deltaPoutPgenInvPgenMuon);
                                       (*muon_energy_).push_back(muon_energy);
                                       (*frac_).push_back(gamma_frac_energy);
				       (*R_).push_back(radius);
				       (*eta_).push_back(eta_photon);
                                       (*phi_).push_back(phi_photon);
				       (*Kenergy_).push_back(gamma_energy);
                                       (*py_).push_back(py);
                                       (*px_).push_back(px); 
                                       (*pdg_).push_back(PDGgamma);
                                       (*Nphotons_).push_back(nbrems);
				       
				
				if(Debug) cout <<" n_brem_photons :"<< nbrems << endl;
				
				     }
			       }
			   }
			}

	t1->Fill(); 

  
		     
	 }
       }
     }


    }
 

}



//define this as a plug-in

DEFINE_FWK_MODULE(testEvent);


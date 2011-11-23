//Framework headers 
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Fast Simulation headers
#include "FastSimulation/Calorimetry/interface/CalorimetryManager.h"
#include "FastSimulation/Event/interface/FSimEvent.h"
#include "FastSimulation/Event/interface/FSimTrack.h"
#include "FastSimulation/ShowerDevelopment/interface/EMECALShowerParametrization.h"
#include "FastSimulation/ShowerDevelopment/interface/EMShower.h"
#include "FastSimulation/ShowerDevelopment/interface/HDShowerParametrization.h"
#include "FastSimulation/ShowerDevelopment/interface/HDShower.h"
#include "FastSimulation/ShowerDevelopment/interface/HFShower.h"
#include "FastSimulation/ShowerDevelopment/interface/HDRShower.h"
#include "FastSimulation/ShowerDevelopment/interface/HSParameters.h"
#include "FastSimulation/CaloGeometryTools/interface/CaloGeometryHelper.h"
#include "FastSimulation/CaloHitMakers/interface/EcalHitMaker.h"
#include "FastSimulation/CaloHitMakers/interface/HcalHitMaker.h"
#include "FastSimulation/CaloHitMakers/interface/PreshowerHitMaker.h"
//#include "FastSimulation/Utilities/interface/Histos.h"
#include "FastSimulation/Utilities/interface/RandomEngine.h"
#include "FastSimulation/Utilities/interface/GammaFunctionGenerator.h"
#include "FastSimulation/Utilities/interface/LandauFluctuationGenerator.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"  
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "FastSimulation/Event/interface/FSimTrackEqual.h"
// New headers for Muon Mip Simulation
#include "FastSimulation/MaterialEffects/interface/MaterialEffects.h"
#include "FastSimulation/MaterialEffects/interface/EnergyLossSimulator.h"
// Muon Brem
#include "FastSimulation/MaterialEffects/interface/MuonBremsstrahlungSimulator.h"

//Gflash Hadronic Model
#include "SimGeneral/GFlash/interface/GflashHadronShowerProfile.h"
#include "SimGeneral/GFlash/interface/GflashPiKShowerProfile.h"
#include "SimGeneral/GFlash/interface/GflashProtonShowerProfile.h"
#include "SimGeneral/GFlash/interface/GflashAntiProtonShowerProfile.h"
#include "SimGeneral/GFlash/interface/GflashTrajectoryPoint.h"
#include "SimGeneral/GFlash/interface/GflashHit.h"
#include "SimGeneral/GFlash/interface/Gflash3Vector.h"
// STL headers 
#include <vector>
#include <iostream>

//CMSSW headers 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
//#include "DataFormats/EcalDetId/interface/EcalDetId.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
//SHP includes (Sandro Fonseca de Souza 22 Jun 2011 Beta test on Calo)
//#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
//#include "TrackingTools/TrajectoryParametrization/interface/GlobalTrajectoryParameters.h"

#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "DataFormats/GeometrySurface/interface/PlaneBuilder.h"
#include "DataFormats/GeometrySurface/interface/TangentPlane.h"
#include "FastSimulation/TrajectoryManager/interface/LocalMagneticField.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"


//photon
#include "FastSimulation/Particle/interface/RawParticle.h"
//#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"



using namespace edm;

typedef math::XYZVector XYZVector;
typedef math::XYZVector XYZPoint;

std::vector<std::pair<int,float> > CalorimetryManager::myZero_ = std::vector<std::pair<int,float> >
(1,std::pair<int,float>(0,0.));

CalorimetryManager::CalorimetryManager() : 
  myCalorimeter_(0),
  myHistos(0),
  random(0),initialized_(false)
{;}

CalorimetryManager::CalorimetryManager(FSimEvent * aSimEvent, 
				       const edm::ParameterSet& fastCalo,
				       const edm::ParameterSet& fastMuECAL,
				       const edm::ParameterSet& fastMuHCAL,
                                       const edm::ParameterSet& parGflash,
				       const RandomEngine* engine)
  : 
  mySimEvent(aSimEvent), 
  random(engine),initialized_(false),
  theMuonEcalEffects(0), theMuonHcalEffects (0),magField_(0), ptr_(0)

{

  aLandauGenerator = new LandauFluctuationGenerator(random);
  aGammaGenerator = new GammaFunctionGenerator(random);

  //Gflash
  theProfile = new GflashHadronShowerProfile(parGflash);
  thePiKProfile = new GflashPiKShowerProfile(parGflash);
  theProtonProfile = new GflashProtonShowerProfile(parGflash);
  theAntiProtonProfile = new GflashAntiProtonShowerProfile(parGflash);

  readParameters(fastCalo);

//  EBMapping_.resize(62000,myZero_);
//  EEMapping_.resize(20000,myZero_);
//  HMapping_.resize(10000,myZero_);
  EBMapping_.resize(62000);
  EEMapping_.resize(20000);
  HMapping_.resize(10000);
  theDetIds_.resize(10000);

  unsigned s=(unfoldedMode_)?5:1;
  for(unsigned ic=0;ic<62000;++ic)
    {
      EBMapping_[ic].reserve(s);
      if(ic<20000)
	EEMapping_[ic].reserve(s);
      if(ic<10000)
	HMapping_[ic].reserve(s);
    }
  


  myHistos = 0; 
#ifdef FAMOSDEBUG
  myHistos = Histos::instance();
  myHistos->book("h10",140,-3.5,3.5,100,-0.5,99.5);
  myHistos->book("h20",150,0,150.,100,-0.5,99.5);
  myHistos->book("h100",140,-3.5,3.5,100,0,0.1);
  myHistos->book("h110",140,-3.5,3.5,100,0,10.);
  myHistos->book("h120",200,-5.,5.,100,0,0.5);

  myHistos->book("h200",300,0,3.,100,0.,35.);
  myHistos->book("h210",720,-M_PI,M_PI,100,0,35.);
  myHistos->book("h212",720,-M_PI,M_PI,100,0,35.);

  myHistos->bookByNumber("h30",0,7,300,-3.,3.,100,0.,35.);
  myHistos->book("h310",75,-3.,3.,"");
  myHistos->book("h400",100,-10.,10.,100,0.,35.);
  myHistos->book("h410",720,-M_PI,M_PI);
#endif



 myCalorimeter_ = 
    new CaloGeometryHelper(fastCalo);
  myHDResponse_ = 
    new HCALResponse(fastCalo.getParameter<edm::ParameterSet>("HCALResponse"),
		     random);
  myHSParameters_ = 
    new HSParameters(fastCalo.getParameter<edm::ParameterSet>("HSParameters"));

  // Material Effects for Muons in ECAL (only EnergyLoss implemented so far)

  if ( fastMuECAL.getParameter<bool>("PairProduction") || 
       fastMuECAL.getParameter<bool>("Bremsstrahlung") ||
       fastMuECAL.getParameter<bool>("MuonBremsstrahlung") ||
       fastMuECAL.getParameter<bool>("EnergyLoss") || 
       fastMuECAL.getParameter<bool>("MultipleScattering") )
    theMuonEcalEffects = new MaterialEffects(fastMuECAL,random);

  // Material Effects for Muons in HCAL (only EnergyLoss implemented so far)

  if ( fastMuHCAL.getParameter<bool>("PairProduction") || 
       fastMuHCAL.getParameter<bool>("Bremsstrahlung") ||
       fastMuHCAL.getParameter<bool>("MuonBremsstrahlung") ||
       fastMuHCAL.getParameter<bool>("EnergyLoss") || 
       fastMuHCAL.getParameter<bool>("MultipleScattering") )
    theMuonHcalEffects = new MaterialEffects(fastMuHCAL,random);

}
//////////////////////////////

void CalorimetryManager::SetMagneticField(const MagneticField* magField){


  
  magField_ = magField;
 



}

////////////////////////////////////////////////////////////////////////

void CalorimetryManager::SetPropagator(const Propagator* ptr){


  
  ptr_ = ptr;
 



}

//////////////////////////////////
void CalorimetryManager::clean()
{
  unsigned size=firedCellsEB_.size();
  for(unsigned ic=0;ic<size;++ic)
    {
      EBMapping_[firedCellsEB_[ic]].clear();
    }
  firedCellsEB_.clear();

  size=firedCellsEE_.size();
  for(unsigned ic=0;ic<size;++ic)
    {
      EEMapping_[firedCellsEE_[ic]].clear();
    }
  firedCellsEE_.clear();
  
  size=firedCellsHCAL_.size();
  for(unsigned ic=0;ic<size;++ic)
    {
      HMapping_[firedCellsHCAL_[ic]].clear();
    }
  firedCellsHCAL_.clear();

  ESMapping_.clear();
  muonSimTracks.clear();
}

CalorimetryManager::~CalorimetryManager()
{
#ifdef FAMOSDEBUG
  myHistos->put("Famos.root");
#endif
  if(myCalorimeter_) delete myCalorimeter_;
  if(myHDResponse_) delete myHDResponse_;

  if ( theMuonEcalEffects ) delete theMuonEcalEffects;
  if ( theMuonHcalEffects ) delete theMuonHcalEffects;

  if ( theProfile ) delete theProfile;
}

void CalorimetryManager::reconstruct()
{
  if(evtsToDebug_.size())
    {
      std::vector<unsigned int>::const_iterator itcheck=find(evtsToDebug_.begin(),evtsToDebug_.end(),mySimEvent->id().event());
      debug_=(itcheck!=evtsToDebug_.end());
      if(debug_)
	mySimEvent->print();
    }
  // Clear the content of the calorimeters 
  if(!initialized_)
    {
      const CaloSubdetectorGeometry* geom=myCalorimeter_->getHcalGeometry();
      for(int subdetn=1;subdetn<=4;++subdetn)
	{
	  const std::vector<DetId>& ids(geom->getValidDetIds(DetId::Hcal,subdetn));  
	  for (std::vector<DetId>::const_iterator i=ids.begin(); i!=ids.end(); i++) 
	    {
	      HcalDetId myDetId(*i);
	      unsigned hi=myDetId.hashed_index();
	      theDetIds_[hi]=myDetId;
	    }
	}
      
      // Check if the preshower is really available
      if(simulatePreshower_ && !myCalorimeter_->preshowerPresent())
	{
	  std::cout << " WARNING " << std::endl;
	  std::cout << " The preshower simulation has been turned on; but no preshower geometry is available " << std::endl;
	  std::cout << " Disabling the preshower simulation " << std::endl;
	  simulatePreshower_ = false;
	}

      initialized_=true;
    }
  clean();

  LogInfo("FastCalorimetry") << "Reconstructing " << (int) mySimEvent->nTracks() << " tracks." << std::endl;
  for( int fsimi=0; fsimi < (int) mySimEvent->nTracks() ; ++fsimi) {

   

    FSimTrack& myTrack = mySimEvent->track(fsimi);

    int pid = abs(myTrack.type());

    if (debug_) {
      LogDebug("FastCalorimetry") << " ===> pid = "  << pid << std::endl;      
    }
    // std::cout << " ===> pid = "  << pid << std::endl; 
    
    // Check that the particle hasn't decayed
    if(myTrack.noEndVertex() ) {
      // Simulate energy smearing for photon and electrons
      if ( pid == 11 || pid == 22  ) {
	  
	if ( myTrack.onEcal()) 
	    EMShowerSimulation(myTrack);
	else if ( myTrack.onVFcal())
	     reconstructHCAL(myTrack);
     

      } // electron or photon
      else if (pid==13)
	{
	  MuonMipSimulation(myTrack);
	}
      // Simulate energy smearing for hadrons (i.e., everything 
      // but muons... and SUSY particles that deserve a special 
      // treatment.
      else if ( pid < 1000000 ) {
	if ( myTrack.onHcal() || myTrack.onVFcal() ) { 	  
	  if(optionHDSim_ == 0 )  reconstructHCAL(myTrack);
	  else HDShowerSimulation(myTrack);
	}
      } // pid < 1000000 
    } // myTrack.noEndVertex()
  } // particle loop
  //  LogInfo("FastCalorimetry") << " Number of  hits (barrel)" << EBMapping_.size() << std::endl;
  //  LogInfo("FastCalorimetry") << " Number of  hits (Hcal)" << HMapping_.size() << std::endl;
  //  std::cout << " Nombre de hit (endcap)" << EEMapping_.size() << std::endl;

} // reconstruct

// Simulation of electromagnetic showers in PS, ECAL, HCAL 
void CalorimetryManager::EMShowerSimulation(const FSimTrack& myTrack) {
  std::vector<const RawParticle*> thePart;
  double X0depth;
  if (debug_) {
    LogDebug("FastCalorimetry") << " EMShowerSimulation "  <<myTrack << std::endl;      
  }
 

  // std::cout << " EMShowerSimulation "  << myTrack << std::endl; 
     
   

  // The Particle at ECAL entrance
  // std::cout << " Before ecalEntrance " << "Energy limit: " << myPart.e() << "PDI: "  << myTrack.type() << std::endl;
    
   myPart = myTrack.ecalEntrance(); 
   //  std::cout << " Before ecalEntrance " << myPart  << std::endl;
  // protection against infinite loop.  
  if ( myTrack.type() == 22 && myPart.e()<0.055) return; 
 
  // Barrel or Endcap ?
  int onEcal = myTrack.onEcal();
  int onHcal = myTrack.onHcal();
  int onLayer1 = myTrack.onLayer1();
  int onLayer2 = myTrack.onLayer2();

  // The entrance in ECAL
  XYZPoint ecalentrance = myPart.vertex().Vect();
  
  
  // The preshower
  PreshowerHitMaker * myPreshower = NULL ;
  if(simulatePreshower_ && (onLayer1 || onLayer2))
    {
      XYZPoint layer1entrance,layer2entrance;
      XYZVector dir1,dir2;
      if(onLayer1) 
	{
	  layer1entrance = XYZPoint(myTrack.layer1Entrance().vertex().Vect());
	  dir1 = XYZVector(myTrack.layer1Entrance().Vect().Unit());
	}
      if(onLayer2) 
	{
	  layer2entrance = XYZPoint(myTrack.layer2Entrance().vertex().Vect());
	  dir2 = XYZVector(myTrack.layer2Entrance().Vect().Unit());
	}
      //      std::cout << " Layer1entrance " << layer1entrance << std::endl;
      //      std::cout << " Layer2entrance " << layer2entrance << std::endl;
      myPreshower = new PreshowerHitMaker(myCalorimeter_,
					  layer1entrance,
					  dir1,
					  layer2entrance,
					  dir2,
					  aLandauGenerator);
      myPreshower->setMipEnergy(mipValues_[0],mipValues_[1]);
    }

  

  // The ECAL Properties
  EMECALShowerParametrization 
    showerparam(myCalorimeter_->ecalProperties(onEcal), 
		myCalorimeter_->hcalProperties(onHcal), 
		myCalorimeter_->layer1Properties(onLayer1), 
		myCalorimeter_->layer2Properties(onLayer2),
		theCoreIntervals_,
		theTailIntervals_,
		RCFactor_,
		RTFactor_);



  // Photons : create an e+e- pair
  if ( myTrack.type() == 22 ) {
    // Depth for the first e+e- pair creation (in X0)
    X0depth = -log(random->flatShoot()) * (9./7.);
    
    // Initialization
    double eMass = 0.000510998902; 
    double xe=0;
    double xm=eMass/myPart.e();
    double weight = 0.;
    
    // Generate electron energy between emass and eGamma-emass
    do {
      xe = random->flatShoot()*(1.-2.*xm) + xm;
      weight = 1. - 4./3.*xe*(1.-xe);
    

    } while ( weight < random->flatShoot() );


    // Protection agains infinite loop in Famos Shower
     if ( myPart.e()*xe < 0.055 || myPart.e()*(1.-xe) < 0.055 ) {
       if ( myPart.e() > 0.055 ) thePart.push_back(&myPart);
      
     } else {      
     
      myElec = (myPart) * xe;
      myPosi = (myPart) * (1.-xe);
      myElec.setVertex(myPart.vertex());
      myPosi.setVertex(myPart.vertex());
      thePart.push_back(&myElec);
      thePart.push_back(&myPosi);
     
     

     }
  // Electrons
  } else {  
    
    X0depth = 0.;
    
    if ( myPart.e() > 0.055 ) thePart.push_back(&myPart);
    //  std::cout << " myPart.e(): "<<  myPart.e() << std::endl;

  } 
  
  // After the different protections, this shouldn't happen. 
  if(thePart.size()==0) 
    { 
      if(myPreshower==NULL) return; 
      delete myPreshower; 
      return; 
    } 

  // find the most energetic particle
  double maxEnergy=-1.;
  for(unsigned ip=0;ip < thePart.size();++ip)
    if(thePart[ip]->e() > maxEnergy) maxEnergy = thePart[ip]->e();

  // Initialize the Grid in ECAL
  int size = gridSize_;
   if(maxEnergy>100) size=11;
   

//  if ( maxEnergy < threshold5x5 ) size = 5;
//  if ( maxEnergy < threshold3x3 ) size = 3;
   

  EMShower theShower(random,aGammaGenerator,&showerparam,&thePart);

  double maxShower = theShower.getMaximumOfShower();


  if (maxShower > 20.) maxShower = 2.; // simple pivot-searching protection 

  double depth((X0depth + maxShower) * 
	       myCalorimeter_->ecalProperties(onEcal)->radLenIncm());
  XYZPoint meanShower = ecalentrance + myPart.Vect().Unit()*depth;

 
  //std::cout<<" onEcal:" << onEcal << std::endl;
  //  if(onEcal!=1) return ; 

  // The closest crystal
  DetId pivot(myCalorimeter_->getClosestCell(meanShower, true, onEcal==1));

  if(pivot.subdetId() == 0) {   // further protection against avbsence of pivot
    edm::LogWarning("CalorimetryManager") <<  "Pivot for egamma  e = "  << myTrack.hcalEntrance().e() << " is not found at depth " << depth << " and meanShower coordinates = " << meanShower << std::endl;
    //   std::cout <<  "Pivot for egamma  e = "  << myTrack.hcalEntrance().e() << " is not found at depth " << depth << " and meanShower coordinates = " << meanShower << std::endl;
    return;
  }
  
  EcalHitMaker myGrid(myCalorimeter_,ecalentrance,pivot,onEcal,size,0,random);
  //                                             ^^^^
  //                                         for EM showers
  myGrid.setPulledPadSurvivalProbability(pulledPadSurvivalProbability_);
  myGrid.setCrackPadSurvivalProbability(crackPadSurvivalProbability_);
  myGrid.setRadiusFactor(radiusFactor_);


 //maximumdepth dependence of the radiusfactorbehindpreshower
 // adjusted by Shilpi Jain (Mar-Apr 2010)
  if(onLayer1 || onLayer2)
    {
      float b               = radiusPreshowerCorrections_[0];
      float a               = radiusFactor_*( 1.+radiusPreshowerCorrections_[1]*radiusPreshowerCorrections_[0] );
      float maxdepth        = X0depth+theShower.getMaximumOfShower();
      float newRadiusFactor = radiusFactor_;
      if(myPart.e()<=250.)
        {
	  newRadiusFactor = a/(1.+b*maxdepth); 
	}
      myGrid.setRadiusFactor(newRadiusFactor);
    }
  else // otherwise use the normal radius factor
    {
      myGrid.setRadiusFactor(radiusFactor_);
    }
  myGrid.setPreshowerPresent(simulatePreshower_);
  
  // The shower simulation
  myGrid.setTrackParameters(myPart.Vect().Unit(),X0depth,myTrack);

//  std::cout << " PS ECAL GAP HCAL X0 " << myGrid.ps1TotalX0()+myGrid.ps2TotalX0() << " " << myGrid.ecalTotalX0();
//  std::cout << " " << myGrid.ecalHcalGapTotalX0() << " " << myGrid.hcalTotalX0() << std::endl;
//  std::cout << " PS ECAL GAP HCAL L0 " << myGrid.ps1TotalL0()+myGrid.ps2TotalL0() << " " << myGrid.ecalTotalL0();
//   std::cout << " " << myGrid.ecalHcalGapTotalL0() << " " << myGrid.hcalTotalL0() << std::endl;
//   std::cout << "ECAL-HCAL " << myTrack.momentum().eta() << " " <<  myGrid.ecalHcalGapTotalL0() << std::endl;
//
//  std::cout << " Grid created " << std::endl;
  if(myPreshower) theShower.setPreshower(myPreshower);
  
  HcalHitMaker myHcalHitMaker(myGrid,(unsigned)0); 

  theShower.setGrid(&myGrid);
  theShower.setHcal(&myHcalHitMaker);
  //  std:: cout << " About to compute " << std::endl;
  theShower.compute();
  //  std::cout << " Coming back from compute" << std::endl;
  //myHistos->fill("h502", myPart->eta(),myGrid.totalX0());
  
  // Save the hits !
  std::map<uint32_t,float>::const_iterator mapitr;
  std::map<uint32_t,float>::const_iterator endmapitr=myGrid.getHits().end();
  
  for(mapitr=myGrid.getHits().begin();mapitr!=endmapitr;++mapitr)
    {
      if(onEcal==1)
	{
	  updateMap(EBDetId(mapitr->first).hashedIndex(), mapitr->second,myTrack.id(),EBMapping_,firedCellsEB_);
	  //std::cout << " Adding Barrel: " <<mapitr->first << " " << mapitr->second <<std::endl;
	  //	  std::cout << " myTrack.id(): " << myTrack.id() <<std::endl;
	}
	    
      else if(onEcal==2)
	updateMap(EEDetId(mapitr->first).hashedIndex(), mapitr->second,myTrack.id(),EEMapping_,firedCellsEE_);
      //if (mapitr->second > 0.) std::cout << " Adding EndCaps: " <<mapitr->first << " " << mapitr->second <<std::endl; 
    }

  // Now fill the HCAL hits
  endmapitr=myHcalHitMaker.getHits().end();
  for(mapitr=myHcalHitMaker.getHits().begin();mapitr!=endmapitr;++mapitr)
    {
      updateMap(HcalDetId(mapitr->first).hashed_index(),mapitr->second,myTrack.id(),HMapping_,firedCellsHCAL_);
      //if (mapitr->second > 0.) std::cout << " Adding HCAL: " <<mapitr->first << " " << mapitr->second <<std::endl; 
    }

  // delete the preshower
  if(myPreshower!=0)
    {
      endmapitr=myPreshower->getHits().end();
      for(mapitr=myPreshower->getHits().begin();mapitr!=endmapitr;++mapitr)
	{
	  updateMap(mapitr->first,mapitr->second,myTrack.id(),ESMapping_);
	  //if (mapitr->second > 0.) std::cout << " Adding PreShower: " <<mapitr->first << " " << mapitr->second <<std::endl; 
	}
      delete myPreshower;
    }
  //  std::cout << " Deleting myPreshower " << std::endl;
  
}



// Simulation of electromagnetic showers in VFCAL
void CalorimetryManager::reconstructECAL(const FSimTrack& track) {
  if(debug_) {
    XYZTLorentzVector moment = track.momentum();
    std::cout << "FASTEnergyReconstructor::reconstructECAL - " << std::endl
	 << "  eta " << moment.eta() << std::endl
         << "  phi " << moment.phi() << std::endl
         << "   et " << moment.Et()  << std::endl;
  }
  
  int hit; 
  
  bool central=track.onEcal()==1;
  
  //Reconstruct only electrons and photons. 

  //deal with different conventions
  // ParticlePropagator 1 <-> Barrel
  //                    2 <-> EC
  // whereas for Artur(this code):
  //                    0 <-> Barrel
  //                    1 <-> EC
  //                    2 <-> VF
  XYZTLorentzVector trackPosition;
  if( track.onEcal() ) {
    hit=track.onEcal()-1;
    trackPosition=track.ecalEntrance().vertex();
  } else {
    hit=2;
    trackPosition=track.vfcalEntrance().vertex();
  }
  
  double pathEta   = trackPosition.eta();
  double pathPhi   = trackPosition.phi();	
  double EGen      = track.ecalEntrance().e();
  

  double e=0.;
  double sigma=0.;
  // if full simulation and in HF, but without showering anyway...
  if(hit == 2 && optionHDSim_ == 2 ) { 
    std::pair<double,double> response =
      myHDResponse_->responseHCAL(0, EGen, pathEta, 0); // last par.= 0 = e/gamma 
    e     = response.first;
    sigma = response.second;
  }

  double emeas = 0.;
  
  if(sigma>0.)
    emeas = random->gaussShoot(e,sigma);

  if(debug_)
    std::cout << "FASTEnergyReconstructor::reconstructECAL : " 
         << "  on-calo  eta, phi  = " << pathEta << " " << pathPhi << std::endl 
	 << "  Egen  = " << EGen << std::endl 
	 << "  Eres  = " << e << std::endl 
	 << " sigma  = " << sigma << std::endl 
	 << "  Emeas = " << emeas << std::endl; 


   if(debug_)
    std::cout << "FASTEnergyReconstructor::reconstructECAL : " 
	 << " Track position - " << trackPosition.Vect() 
	 << "   bool central - " << central
         << "   hit - " << hit   << std::endl;  

  DetId detid;  
  if( hit==2 ) 
      detid = myCalorimeter_->getClosestCell(trackPosition.Vect(),false,central);
  // Check that the detid is HCAL forward
  HcalDetId hdetid(detid);
  if(!hdetid.subdetId()!=HcalForward) return;

  if(debug_)
    std::cout << "FASTEnergyReconstructor::reconstructECAL : " 
	      << " CellID - " <<  detid.rawId() << std::endl;

  if( hit != 2  || emeas > 0.) 
    if(!detid.null()) 
      {
	updateMap(hdetid.hashed_index(),emeas,track.id(),HMapping_,firedCellsHCAL_);
      }

}


void CalorimetryManager::reconstructHCAL(const FSimTrack& myTrack)
{
  int hit;
  int pid = abs(myTrack.type());
  if (debug_) {
    LogDebug("FastCalorimetry") << " reconstructHCAL "  << myTrack << std::endl;      
  }

  //  FSimTrack myTrack = mySimEvent.track(fsimi);

  //  int pid=abs(myTrack.type());
  //  std::cout << "reconstructHCAL " << std::endl;
  
  XYZTLorentzVector trackPosition;
  if (myTrack.onHcal()) {
    trackPosition=myTrack.hcalEntrance().vertex();
    hit = myTrack.onHcal()-1;
  } else {
    trackPosition=myTrack.vfcalEntrance().vertex();
    hit = 2;
  }

  double pathEta   = trackPosition.eta();
  double pathPhi   = trackPosition.phi();	
  //  double pathTheta = trackPosition.theta();

  double EGen  = myTrack.hcalEntrance().e();
  double e     = 0.;
  double sigma = 0.;
  double emeas = 0.;

  if(pid == 13) { 
    //    std::cout << " We should not be here " << std::endl;
    std::pair<double,double> response =
      myHDResponse_->responseHCAL(0, EGen, pathEta, 2); // 2=muon 
    emeas  = response.first;
    if(debug_)
      LogDebug("FastCalorimetry") << "CalorimetryManager::reconstructHCAL - MUON !!!" << std::endl;
  }
  else if( pid == 22 || pid == 11)
    {
      
      std::pair<double,double> response =
	myHDResponse_->responseHCAL(0, EGen, pathEta, 0); // last par. = 0 = e/gamma
      e     = response.first;              //
      sigma = response.second;             //
      emeas = random->gaussShoot(e,sigma); //

      //  cout <<  "CalorimetryManager::reconstructHCAL - e/gamma !!!" << std::endl;
      if(debug_)
	LogDebug("FastCalorimetry") << "CalorimetryManager::reconstructHCAL - e/gamma !!!" << std::endl;
    }
    else {
      e     = myHDResponse_->getHCALEnergyResponse(EGen,hit);
      sigma = myHDResponse_->getHCALEnergyResolution(EGen, hit);
      
      emeas = random->gaussShoot(e,sigma);  
    }
    

  if(debug_)
    LogDebug("FastCalorimetry") << "CalorimetryManager::reconstructHCAL - on-calo "   
				<< "  eta = " << pathEta 
				<< "  phi = " << pathPhi 
				<< "  Egen = " << EGen 
				<< "  Eres = " << e 
				<< "  sigma = " << sigma 
				<< "  Emeas = " << emeas << std::endl;

  if(emeas > 0.) {  
    DetId cell = myCalorimeter_->getClosestCell(trackPosition.Vect(),false,false);
    updateMap(HcalDetId(cell).hashed_index(), emeas, myTrack.id(),HMapping_,firedCellsHCAL_);
  }
}

void CalorimetryManager::HDShowerSimulation(const FSimTrack& myTrack)
{
  //  TimeMe t(" FASTEnergyReconstructor::HDShower");
  XYZTLorentzVector moment = myTrack.momentum();
  
  if(debug_)
    LogDebug("FastCalorimetry") 
      << "CalorimetryManager::HDShowerSimulation - track param."
      << std::endl
      << "  eta = " << moment.eta() << std::endl
      << "  phi = " << moment.phi() << std::endl
      << "   et = " << moment.Et()  << std::endl
      << "   e  = " << myTrack.hcalEntrance().e() << std::endl;

  if (debug_) {
      LogDebug("FastCalorimetry") << " HDShowerSimulation "  << myTrack << std::endl;      
    }


  int hit;
  //  int pid = abs(myTrack.type());

  XYZTLorentzVector trackPosition;
  if ( myTrack.onEcal() ) {
    trackPosition=myTrack.ecalEntrance().vertex();
    hit = myTrack.onEcal()-1;                               //
    myPart = myTrack.ecalEntrance();
  } else if ( myTrack.onVFcal()) {
    trackPosition=myTrack.vfcalEntrance().vertex();
    hit = 2;
    myPart = myTrack.vfcalEntrance();
  }
  else
    {
      LogDebug("FastCalorimetry") << " The particle is not in the acceptance " << std::endl;
      return;
    }

  // int onHCAL = hit + 1; - specially for myCalorimeter->hcalProperties(onHCAL)
  // (below) to get VFcal properties ...
  int onHCAL = hit + 1;
  int onECAL = myTrack.onEcal();
  
  double pathEta   = trackPosition.eta();
  double pathPhi   = trackPosition.phi();	
  //  double pathTheta = trackPosition.theta();

  double eint  = moment.e();
  double eGen  = myTrack.hcalEntrance().e();
  double e     = 0.;
  double sigma = 0.;

  double emeas = 0.;  
  
  //===========================================================================
  if(eGen > 0.) {  

    // ECAL and HCAL properties to get
    HDShowerParametrization 
      theHDShowerparam(myCalorimeter_->ecalProperties(onECAL),
		       myCalorimeter_->hcalProperties(onHCAL),
		       myHSParameters_);
    
    //Making ECAL Grid (and segments calculation)
    XYZPoint caloentrance;
    XYZVector direction;
    if(myTrack.onEcal()) 
      {	
	caloentrance = myTrack.ecalEntrance().vertex().Vect();
	direction = myTrack.ecalEntrance().Vect().Unit();
      }
    else if(myTrack.onHcal())
      {
	caloentrance = myTrack.hcalEntrance().vertex().Vect();
	direction = myTrack.hcalEntrance().Vect().Unit();
      }
    else
      {
	caloentrance = myTrack.vfcalEntrance().vertex().Vect();
	direction = myTrack.vfcalEntrance().Vect().Unit();
      }

  if(debug_)
    LogDebug("FastCalorimetry") 
      << "CalorimetryManager::HDShowerSimulation - on-calo 1 "
      << std::endl
      << "  onEcal    = " <<  myTrack.onEcal()  << std::endl
      << "  onHcal    = " <<  myTrack.onHcal()  << std::endl
      << "  onVFcal   = " <<  myTrack.onVFcal() << std::endl
      << "  position  = " << caloentrance << std::endl;


    DetId pivot;
    if(myTrack.onEcal())
      {
	pivot=myCalorimeter_->getClosestCell(caloentrance,
					     true, myTrack.onEcal()==1);
      }
    else if(myTrack.onHcal())
      {
	//	std::cout << " CalorimetryManager onHcal " <<  myTrack.onHcal() << " caloentrance" << caloentrance  << std::endl;
	pivot=myCalorimeter_->getClosestCell(caloentrance,					     
					    false, false);
      }

    EcalHitMaker myGrid(myCalorimeter_,caloentrance,pivot,
			pivot.null()? 0 : myTrack.onEcal(),hdGridSize_,1,
			random);
    // 1=HAD shower

    myGrid.setTrackParameters(direction,0,myTrack);
    // Build the FAMOS HCAL 
    HcalHitMaker myHcalHitMaker(myGrid,(unsigned)1); 
    
    // Shower simulation
    bool status = false;
    int  mip = 2;
    // Use HFShower for HF
    if ( !myTrack.onEcal() && !myTrack.onHcal() ) {
      //      std::cout << "CalorimetryManager::HDShowerSimulation(): track entrance = "
      //		<< myTrack.vfcalEntrance().vertex().X() << " "
      //		<< myTrack.vfcalEntrance().vertex().Y() << " "
      //		<< myTrack.vfcalEntrance().vertex().Z() << " "
      //		<< " , Energy (Gen/Scale) = " << eGen << " " << e << std::endl;

      // Warning : We give here the particle energy with the response
      //           but without the resolution/gaussian smearing
      //           For HF, the resolution is due to the PE statistic

      HFShower theShower(random,
			 &theHDShowerparam,
			 &myGrid,
			 &myHcalHitMaker,
			 onECAL,
			 eGen);
			 //			 eGen);
			 //			 e); // PV Warning : temporarly set the energy to the generated E

      status = theShower.compute();
    } else { 
      if(hdSimMethod_ == 0) {
	HDShower theShower(random,
			   &theHDShowerparam,
			   &myGrid,
			   &myHcalHitMaker,
			   onECAL,
			   eGen);
	status = theShower.compute();
        mip    = theShower.getmip();
      }
      else if (hdSimMethod_ == 1) {
	HDRShower theShower(random,
			    &theHDShowerparam,
			    &myGrid,
			    &myHcalHitMaker,
			    onECAL,
			    eGen);
	status = theShower.computeShower();
        mip = 2;
      }
      else if (hdSimMethod_ == 2 ) {
	//        std::cout << "Using GflashHadronShowerProfile hdSimMethod_ == 2" << std::endl;

        //dynamically loading a corresponding profile by the particle type
        int particleType = myTrack.type();
        theProfile = thePiKProfile;
        if(particleType == -2212) theProfile = theAntiProtonProfile;
        else if(particleType == 2212) theProfile = theProtonProfile;

        //input variables for GflashHadronShowerProfile
        int showerType = 99 + myTrack.onEcal();
        double globalTime = 150.0; // a temporary reference hit time in nanosecond
        float charge = (float)(myTrack.charge());
        Gflash3Vector gfpos(trackPosition.X(),trackPosition.Y(),trackPosition.Z());
        Gflash3Vector gfmom(moment.X(),moment.Y(),moment.Z());

        theProfile->initialize(showerType,eGen,globalTime,charge,gfpos,gfmom);
        theProfile->loadParameters();
        theProfile->hadronicParameterization();

        //make hits
	std::vector<GflashHit>& gflashHitList = theProfile->getGflashHitList();
	std::vector<GflashHit>::const_iterator spotIter    = gflashHitList.begin();
	std::vector<GflashHit>::const_iterator spotIterEnd = gflashHitList.end();

	Gflash::CalorimeterNumber whichCalor = Gflash::kNULL;
        bool result;

        for( ; spotIter != spotIterEnd; spotIter++){

          double pathLength = theProfile->getGflashShowino()->getPathLengthAtShower()
            + (30*100/eGen)*(spotIter->getTime() - globalTime);

          double currentDepth = std::max(0.0,pathLength - theProfile->getGflashShowino()->getPathLengthOnEcal());

          //find the the showino position at the currentDepth
          GflashTrajectoryPoint trajectoryPoint;
          theProfile->getGflashShowino()->getHelix()->getGflashTrajectoryPoint(trajectoryPoint,pathLength);
          Gflash3Vector positionAtCurrentDepth = trajectoryPoint.getPosition();
          //find radial distrance
          Gflash3Vector lateralDisplacement = positionAtCurrentDepth - spotIter->getPosition()/CLHEP::cm;
          double rShower = lateralDisplacement.r();
          double azimuthalAngle = lateralDisplacement.phi();

          whichCalor = Gflash::getCalorimeterNumber(positionAtCurrentDepth);

          if(whichCalor==Gflash::kESPM || whichCalor==Gflash::kENCA) {
            bool statusPad = myGrid.getPads(currentDepth,true);
            if(!statusPad) continue;
            myGrid.setSpotEnergy(1.2*spotIter->getEnergy()/CLHEP::GeV);
            result = myGrid.addHit(rShower/Gflash::intLength[Gflash::kESPM],azimuthalAngle,0);
          }
          else if(whichCalor==Gflash::kHB || whichCalor==Gflash::kHE) {
            bool setHDdepth = myHcalHitMaker.setDepth(currentDepth,true);
            if(!setHDdepth) continue;
            myHcalHitMaker.setSpotEnergy(1.4*spotIter->getEnergy()/CLHEP::GeV);
            result = myHcalHitMaker.addHit(rShower/Gflash::intLength[Gflash::kHB],azimuthalAngle,0);
          }
        }
        status = true;
      }
      else {
	edm::LogInfo("FastSimulationCalorimetry") << " SimMethod " << hdSimMethod_ <<" is NOT available ";
      }
    }


    if(status) {

      // Here to switch between simple formulae and parameterized response
      if(optionHDSim_ == 1) {
	e     = myHDResponse_->getHCALEnergyResponse  (eGen, hit);
	sigma = myHDResponse_->getHCALEnergyResolution(eGen, hit);
      }
      else { // optionHDsim == 2
	std::pair<double,double> response =
	  myHDResponse_->responseHCAL(mip, eGen, pathEta, 1); // 1=hadron
	e     = response.first;
	sigma = response.second;
      }
      
      emeas = random->gaussShoot(e,sigma);      
      double correction = emeas / eGen;
      
      // RespCorrP factors (ECAL and HCAL separately) calculation
      respCorr(eint);     

      if(debug_)
	LogDebug("FastCalorimetry") 
	  << "CalorimetryManager::HDShowerSimulation - on-calo 2" << std::endl
	  << "   eta  = " << pathEta << std::endl
	  << "   phi  = " << pathPhi << std::endl
	  << "  Egen  = " << eGen << std::endl
	  << "  Eres  = " << e << std::endl
	  << " sigma  = " << sigma << std::endl
	  << " Emeas  = " << emeas << std::endl
	  << "  corr  = " << correction << std::endl
	  << "   mip  = " << mip << std::endl;
      
      
      // was map<unsigned,double> but CaloHitMaker uses float
      std::map<unsigned,float>::const_iterator mapitr;
      std::map<unsigned,float>::const_iterator endmapitr;
      if(myTrack.onEcal() > 0) {
	// Save ECAL hits 
	endmapitr=myGrid.getHits().end();
	for(mapitr=myGrid.getHits().begin(); mapitr!=endmapitr; ++mapitr) {
	  double energy = mapitr->second;
          energy *= correction;              // RESCALING 
          energy *= ecorr;

	  if(energy > 0.000001) { 
	    if(onECAL==1)
		updateMap(EBDetId(mapitr->first).hashedIndex(),energy,myTrack.id(),EBMapping_,firedCellsEB_);

	    else if(onECAL==2)
	      updateMap(EEDetId(mapitr->first).hashedIndex(),energy,myTrack.id(),EEMapping_,firedCellsEE_);

	    if(debug_)
	      LogDebug("FastCalorimetry") << " ECAL cell " << mapitr->first << " added,  E = " 
		   << energy << std::endl;  
	  }
	}
      }
      
      // Save HCAL hits
      endmapitr=myHcalHitMaker.getHits().end();
      for(mapitr=myHcalHitMaker.getHits().begin(); mapitr!=endmapitr; ++mapitr) {
	double energy = mapitr->second;
	energy *= correction;               // RESCALING 
	energy *= hcorr;

	updateMap(HcalDetId(mapitr->first).hashed_index(),energy,myTrack.id(),HMapping_,firedCellsHCAL_);
	if(debug_)
	  LogDebug("FastCalorimetry") << " HCAL cell "  
	       << mapitr->first << " added    E = " 
	       << mapitr->second << std::endl;  
      }
    }      
    else {  // shower simulation failed  
//      std::cout << " Shower simulation failed " << trackPosition.Vect() << std::endl;
//      std::cout << " The FSimTrack " << myTrack << std::endl;
//      std::cout << " HF entrance on VFcal" << myTrack.onVFcal() << std::endl;
//      std::cout << " trackPosition.eta() " << trackPosition.eta() << std::endl;
      if(myTrack.onHcal() || myTrack.onVFcal())
	{
	  DetId cell = myCalorimeter_->getClosestCell(trackPosition.Vect(),false,false);
	  updateMap(HcalDetId(cell).hashed_index(),emeas,myTrack.id(),HMapping_,firedCellsHCAL_);
	  if(debug_)
	    LogDebug("FastCalorimetry") << " HCAL simple cell "   
					<< cell.rawId() << " added    E = " 
					<< emeas << std::endl;  
	}
    }

  } // e > 0. ...

  if(debug_)
    LogDebug("FastCalorimetry") << std::endl << " FASTEnergyReconstructor::HDShowerSimulation  finished "
	 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void CalorimetryManager::MuonMipSimulation(const FSimTrack& myTrack)
{
//ref:http://cmslxr.fnal.gov/lxr/source/Calibration/HcalAlCaRecoProducers/src/AlCaHOCalibProducer.cc#775
//Initialization SHP
 SteppingHelixPropagator myHelix(*&magField_,anyDirection);
 
  //  TimeMe t(" FASTEnergyReconstructor::HDShower");
  XYZTLorentzVector moment = myTrack.momentum();

  // Backward compatibility behaviour
  if(!theMuonHcalEffects) 
    {
      if(myTrack.onHcal() || myTrack.onVFcal() ) 
	reconstructHCAL(myTrack);

      return;
    }

  if(debug_)
    LogDebug("FastCalorimetry") << "CalorimetryManager::MuonMipSimulation - track param."
         << std::endl
	 << "  eta = " << moment.eta() << std::endl
         << "  phi = " << moment.phi() << std::endl
         << "   et = " << moment.Et()  << std::endl;

  int hit;
  //  int pid = abs(myTrack.type());

  // XYZTLorentzVector trackPosition;
 GlobalPoint startingPosition;
 GlobalVector startingMomentum;
 if ( myTrack.onEcal() ){
      // trackPosition=myTrack.ecalEntrance().vertex();
    hit = myTrack.onEcal()-1;                               //
    myPart = myTrack.ecalEntrance();

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first implementation of influence of Magnetic field on muon trajectory into ECAL  
// Sandro F. de Souza
     
    startingPosition = GlobalPoint(myTrack.ecalEntrance().x(),myTrack.ecalEntrance().y(),myTrack.ecalEntrance().z());



    startingMomentum = GlobalVector(myTrack.ecalEntrance().momentum().x(),myTrack.ecalEntrance().momentum().y(),myTrack.ecalEntrance().momentum().z());

     LogDebug("FastCalorimetry")<< " the Muon START Position on ECAL:" << startingPosition << std::endl;
     LogDebug("FastCalorimetry")<< " the Muon START momentum on ECAL:" << startingMomentum << std::endl;
   
  

 } else if ( myTrack.onVFcal()) {
    //    trackPosition=myTrack.vfcalEntrance().vertex();
    hit = 2;
    myPart = myTrack.vfcalEntrance();
 


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// first implementation of influence of Magnetic field on muon trajectory into  VFCAL 
// Sandro F. de Souza
     
    startingPosition = GlobalPoint(myTrack.vfcalEntrance().x(),myTrack.vfcalEntrance().y(),myTrack.vfcalEntrance().z());



    startingMomentum = GlobalVector(myTrack.vfcalEntrance().momentum().x(),myTrack.vfcalEntrance().momentum().y(),myTrack.vfcalEntrance().momentum().z());


    



    LogDebug("FastCalorimetry")<< " the Muon START Position on VFCAL:" << startingPosition << std::endl;
    LogDebug("FastCalorimetry")<< " the Muon START momentum on VFCAL:" << startingMomentum << std::endl;
   

 }
  else
    {
      LogDebug("FastCalorimetry") << " The particle is not in the acceptance " << std::endl;
      return;
    }

  //Setting propagator/////////////////////////////
//Rotation Matrix
    float r21 = 0;
    float r22 = (startingMomentum.unit()).z()/sqrt(1-pow((startingMomentum.unit()).x(),2));
    float r23 = -(startingMomentum.unit()).y()/sqrt(1-pow((startingMomentum.unit()).x(),2));
    float r31 = (startingMomentum.unit()).x();
    float r32 = (startingMomentum.unit()).y();
    float r33 = (startingMomentum.unit()).z();
    float r11 = r22*r33-r23*r32;
    float r12 = r23*r31;
    float r13 = -r22*r31;

    // std::cout << " the Muon START Position:" << startingPosition << std::endl;
    //  std::cout << " the Muon START momentum:" << startingMomentum << std::endl;
     
  int charge = myTrack.charge(); 

  Surface::RotationType rotation(r11, r12, r13,
                                  r21, r22, r23,
                                  r31, r32, r33);
  PlaneBuilder pb;
  PlaneBuilder::ReturnType target = pb.plane(startingPosition,rotation);
  GlobalTrajectoryParameters fts(startingPosition,
                                 startingMomentum,
                                 charge,
				 magField_);


  TrajectoryStateOnSurface state (fts, *target);
  TrajectoryStateOnSurface propagatedState = state;

  SteppingHelixStateInfo shsStart(*(propagatedState.freeTrajectoryState()));
  LogDebug("FastCalorimetry")<< "Propagating " << propagatedState  << std::endl;
  // std::cout << "Propagating " << propagatedState  << std::endl;
  /////////////////////////////////////////////////////////////////////////////////////


  // int onHCAL = hit + 1; - specially for myCalorimeter->hcalProperties(onHCAL)
  // (below) to get VFcal properties ...
  // not needed ? 
  //  int onHCAL = hit + 1; 
  int onECAL = myTrack.onEcal();
  
  //  double pathEta   = trackPosition.eta();
  //  double pathPhi   = trackPosition.phi();	
  //  double pathTheta = trackPosition.theta();
  
  //===========================================================================

  // ECAL and HCAL properties to get
      
  //Making ECAL Grid (and segments calculation)
  XYZPoint caloentrance;
  XYZVector direction;
  if(myTrack.onEcal()) 
    {	
      caloentrance = myTrack.ecalEntrance().vertex().Vect();
      direction = myTrack.ecalEntrance().Vect().Unit();
 




   }
  else if(myTrack.onHcal())
    {
      caloentrance = myTrack.hcalEntrance().vertex().Vect();
      direction = myTrack.hcalEntrance().Vect().Unit();
  

  }
  else
    {
      caloentrance = myTrack.vfcalEntrance().vertex().Vect();
      direction = myTrack.vfcalEntrance().Vect().Unit();
    }
  
  DetId pivot;
  if(myTrack.onEcal())
    {
      pivot=myCalorimeter_->getClosestCell(caloentrance,
					   true, myTrack.onEcal()==1);
    }
  else if(myTrack.onHcal())
    {
      std::cout << " CalorimetryManager onHcal " <<  myTrack.onHcal() << " caloentrance" << caloentrance  << std::endl;
      
       pivot=myCalorimeter_->getClosestCell(caloentrance,					     
					   false, false);
    }
  
  EcalHitMaker myGrid(myCalorimeter_,caloentrance,pivot,
		      pivot.null()? 0 : myTrack.onEcal(),hdGridSize_,0,
		      random);
  // 0 =EM shower -> Unit = X0
  
  myGrid.setTrackParameters(direction,0,myTrack);
  
  const std::vector<CaloSegment>& segments=myGrid.getSegments();
  unsigned nsegments=segments.size();

  GlobalPoint Ecal_entrance_segment;
  GlobalPoint Ecal_exit_segment;

  int ifirstHcal=-1;
  int ilastEcal=-1;
  EnergyLossSimulator* energyLossECAL = 0;
 //Muon Brem and photo production on ECAL 
//Sandro Fonseca de Souza (UERJ/Brazil) and Dilson de Jesus Damiao (CBPF/Brazil)
  //Date 12/10/2011

  MuonBremsstrahlungSimulator* muonBremECAL = 0;
 
  if (theMuonEcalEffects) energyLossECAL = theMuonEcalEffects->energyLossSimulator();
  // Muon brem in ECAL
  if (theMuonEcalEffects) muonBremECAL = theMuonEcalEffects->muonBremsstrahlungSimulator();

  for(unsigned iseg=0;iseg<nsegments&&ifirstHcal<0;++iseg)
    {
      
      // in the ECAL, there are two types of segments: PbWO4 and GAP
      float segmentSizeinX0=segments[iseg].X0length();
      
      // Martijn - insert your computations here
      float energy=0.0;
      if (segmentSizeinX0>0.001 && segments[iseg].material()==CaloSegment::PbWO4 ) {
     
//Information of particle in segments
     Ecal_entrance_segment = GlobalPoint(segments[iseg].entrance().x(),
                                   segments[iseg].entrance().y(),
                                   segments[iseg].entrance().z());
     
     Ecal_exit_segment = GlobalPoint(segments[iseg].exit().x(),
                                   segments[iseg].exit().y(),
                                   segments[iseg].exit().z());

    
      LogDebug("FastCalorimetry")<< "position_segment_entrance: "<< segments[iseg].entrance()<< std::endl;   
      LogDebug("FastCalorimetry")<< "position_segment_exit: "<< segments[iseg].exit()<< std::endl;

//    std::cout<< "length (in cm): "<< length << std::endl;
 //   std::cout<< "segmentSizeinX0: "<< segmentSizeinX0 << std::endl;
 //   std::cout<< "RadPath on SHP: "<< shsDest.radPath() << std::endl;
//    std::cout<< "Path on SHP: "<< shsDest.path() << std::endl;

//Initialization of Stepping Helix State info using the information of segments

 const SteppingHelixStateInfo& shsDest = ((const SteppingHelixPropagator*)ptr_)->propagate(shsStart,Ecal_entrance_segment,Ecal_exit_segment);

               LogDebug("FastCalorimetry")<<"Global Position calculates using SHP in ECAL: "<<  shsDest.position()<< std::endl;
               LogDebug("FastCalorimetry")<<"Global Momentum calculates using SHP in ECAL: "<<  shsDest.momentum()<< std::endl;

	     //   std::cout <<"Global Position calculates using SHP in ECAL: "<<  shsDest.position()<< std::endl;
// 	       std::cout <<"Global Momentum calculates using SHP in ECAL: "<<  shsDest.momentum()<< std::endl;

//Initialization the momentum, position charge and energy on SHS

       GlobalPoint gPos(shsDest.position().x(),shsDest.position().y(),shsDest.position().z());
       GlobalVector gMom(shsDest.momentum().x(),shsDest.momentum().y(),shsDest.momentum().z());
//Muon Mass
       double mu = 0.1056583692;
//initial muon energy
       double ini_en = std::sqrt(startingMomentum.mag2()+mu*mu);

//Muon Energy
       double en = std::sqrt(gMom.mag2()+mu*mu);

//delta energy due SHP (Contribution of the flutuation of muon energy )
       double delta_energy = en - ini_en;
       LogDebug("FastCalorimetry")<<"Flutuation of muon energy due SHP (delta_energy):  "<<  fabs(delta_energy) << std::endl;
       //  std::cout <<"Flutuation of muon energy due SHP om ECAL (delta_energy):  "<<  fabs(delta_energy) << std::endl;


       XYZTLorentzVector trackPosition_shsDest(gPos.x(),gPos.y(),gPos.z(),0.);
       XYZTLorentzVector moment_shsDest(gMom.x(),gMom.y(),gMom.z(),en);
       float charge = (float)(shsDest.charge());
 // The energy loss simulator
	ParticlePropagator theMuon(moment_shsDest,trackPosition_shsDest,charge,0);
        theMuon.setID(-(int)charge*13);
	if ( energyLossECAL ) { 
          energyLossECAL->updateState(theMuon, segmentSizeinX0);
//include of correction enrgy due effect of Mag Field´on Muon Trajectory
	  energy = (energyLossECAL->deltaMom().E() + fabs(delta_energy));
	  moment -= energyLossECAL->deltaMom();
        //std::cout<< "Energy loss on ECAL: "<< energyLossECAL->deltaMom().E()<< std::endl;

	 //Implementation of Muon Brem in ECAL (only applied if also energyLossECAL is switched on)
        if (muonBremECAL){
            muonBremECAL->updateState(theMuon, segmentSizeinX0);
            energy += muonBremECAL->deltaP_BremPhoton().E();
            moment -=  muonBremECAL->deltaP_BremMuon();
	    
          XYZTLorentzVector thePhoton4Momentum(muonBremECAL->deltaP_BremPhoton().Px(),
                                               muonBremECAL->deltaP_BremPhoton().Py(),
                                               muonBremECAL ->deltaP_BremPhoton().Pz(),
                                               muonBremECAL->deltaP_BremPhoton().E());
     
	  if(muonBremECAL->deltaP_BremPhoton().E()>0.){

         // Add a photon
          RawParticle thePhoton(22,thePhoton4Momentum);
	  thePhoton.setID(22);
         //Using new position of muon with brem photon vertex (check it) 
	
          XYZTLorentzVector vertex (myTrack.vertex().position().x(),
			      myTrack.vertex().position().y(),
			      myTrack.vertex().position().z(),
			      myTrack.vertex().position().t());



          thePhoton.setVertex(vertex);

	  //  std ::cout <<" Muon Brem Photon vertex: " << thePhoton.vertex() <<std::endl;

         _theUpdatedState.push_back(thePhoton);
        
	
         //Creating Hits of the photon 
         int ivertex = mySimEvent->addSimVertex(vertex,-1,FSimVertexType::BREM_VERTEX);
         mySimEvent->addSimTrack(&thePhoton,ivertex);

         int fsimi = ((int) mySimEvent->nTracks()-1);
         FSimTrack& myTrack = mySimEvent->track(fsimi);
	 myTrack.setEcal(thePhoton,onECAL);

    ///////////////////////////////////////////////////////
	 // Photons : create an e+e- pair
	 if( myTrack.noEndVertex() ) {

	   // Simulate energy smearing for photon and electrons
	   if (abs( myTrack.type()) == 11 ||  myTrack.type() == 22 ) {
	     LogDebug("FastCalorimetry") << "Using EMShower by muons on ECAL:" <<std::endl;  
	     if ( myTrack.onEcal()) 
	       EMShowerSimulation(myTrack);
	     else if ( myTrack.onVFcal())
	      LogDebug("FastCalorimetry") <<  "Using EMShower by muons on VFCAL" << std::endl;  
	     reconstructHCAL(myTrack);
 
	   }
	 }

  ///////////////////////////////////////////////////////////////////


	
///////////
    if(muonBremECAL->deltaP_BremPhoton().E()>0.){ LogDebug("FastCalorimetry") << "Muon Brem energy on ECAL: "<< energy << std::endl; }//mySimEvent->print(); thePhoton.print();}
           if( debug_ ) {
              float muon_mip =  energyLossECAL->deltaMom().E();
              float muon_brem_ecal = muonBremECAL->deltaP_BremPhoton().E();
              if ( muon_brem_ecal > 0.){
                LogDebug("FastCalorimetry")<< "Energy loss ONLY by Muon Brem on ECAL: "<< muon_brem_ecal<< std::endl;
                LogDebug("FastCalorimetry")<< "Energy loss ONLY by MIP ECAL: "<<  muon_mip << std::endl;
                LogDebug("FastCalorimetry")<< "Photon momentum:"<< muonBremECAL->deltaP_BremPhoton().P()<<" Photon Energy: "<< muonBremECAL->deltaP_BremPhoton().E() << std::endl;
                LogDebug("FastCalorimetry")<< "Contribution of Muon Brem by MUON Momentum on ECAL : "<< muonBremECAL->deltaP_BremMuon().P() << std::endl;
                LogDebug("FastCalorimetry")<< "Contribution of MIP  by MUON Momentum on HCAL: "<< energyLossECAL->deltaMom().P() << std::endl;
              }
            }

            }
	}
 

}
          else{
            LogDebug("FastCalorimetry")<<"Using only MIP for MUON"<<std::endl;
            LogDebug("FastCalorimetry")<< "Muon moment using ONLY the contribution of MIP "<<   moment <<"moment:"<< energyLossECAL->deltaMom() << std::endl;
          }





  }

  // that's all for ECAL, Florian
      // Save the hit only if it is a crystal
      if(segments[iseg].material()==CaloSegment::PbWO4)
	{
	  myGrid.getPads(segments[iseg].sX0Entrance()+segmentSizeinX0*0.5);
	  myGrid.setSpotEnergy(energy);
          myGrid.addHit(0.,0.);
	  ilastEcal=iseg;
	}
      // Check for end of loop:
      if(segments[iseg].material()==CaloSegment::HCAL)
	{
	  ifirstHcal=iseg;
	}
    }
  
  // Position of Ecal Exit
  math::XYZVector ecalExit;
  if(ilastEcal>=0)
    math::XYZVector ecalExit=segments[ilastEcal].exit();
  
  // Position of HCAL entrance
  math::XYZVector hcalEntrance;
  if(ifirstHcal>=0)
    hcalEntrance=segments[ifirstHcal].entrance();

  // dummy HCAL exit, just to test the FSimTrack saving mechanism
  math::XYZVector hcalExit;
  if(ifirstHcal>=0)
    hcalExit=segments[ifirstHcal].exit();


  // Build the FAMOS HCAL 
  HcalHitMaker myHcalHitMaker(myGrid,(unsigned)2);     
  // float mipenergy=0.1;
  // Create the helix with the stepping helix propagator
  // to add a hit, just do
  //  myHcalHitMaker.setSpotEnergy(mipenergy);
  // myHcalHitMaker.addHit(hcalEntrance);
  ///
  /////


  ////// TEMPORARY First attempt to include HCAL 
  int onHCAL = myTrack.onHcal();
  int ilastHcal=-1;
  float mipenergy=0.0;
  GlobalPoint HCAL_entrance_segment;
  GlobalPoint HCAL_exit_segment;

  EnergyLossSimulator* energyLossHCAL = 0;
  // Implementation of Muon Brem Effect in HCAL 
  //Sandro Fonseca de Souza (UERJ/Brazil) and Dilson de Jesus Damiao (CBPF/Brazil)
  //Date 12/10/2011
  MuonBremsstrahlungSimulator* muonBremHCAL = 0;
 
  if (theMuonHcalEffects) energyLossHCAL = theMuonHcalEffects->energyLossSimulator();
  // Muon Brem effect
  if (theMuonHcalEffects) muonBremHCAL = theMuonHcalEffects->muonBremsstrahlungSimulator(); 

  if(ifirstHcal>0 && energyLossHCAL){
    for(unsigned iseg=ifirstHcal;iseg<nsegments;++iseg)
      {
	float segmentSizeinX0=segments[iseg].X0length();
	if (segmentSizeinX0>0.001 && segments[iseg].material()==CaloSegment::HCAL ) {


//Information of particle in segments
     HCAL_entrance_segment = GlobalPoint(segments[iseg].entrance().x(),
                                   segments[iseg].entrance().y(),
                                   segments[iseg].entrance().z());
     
     HCAL_exit_segment = GlobalPoint(segments[iseg].exit().x(),
                                   segments[iseg].exit().y(),
                                   segments[iseg].exit().z());

     //  LogDebug("FastCalorimetry")<< "Propagating on HCAL " << propagatedState  << std::endl;
   LogDebug("FastCalorimetry")<< "position_segment_entrance on HCAL: "<< segments[iseg].entrance()<< std::endl;   
   LogDebug("FastCalorimetry")<< "position_segment_exit on HCAL: "<< segments[iseg].exit()<< std::endl;


//Initialization of Stepping Helix State info using the information of segments(HCAL)

 const SteppingHelixStateInfo& shsDest = ((const SteppingHelixPropagator*)ptr_)->propagate(shsStart,HCAL_entrance_segment,HCAL_exit_segment);

 LogDebug("FastCalorimetry")<<"Global Position calculates using SHP on HCAL: "<<  shsDest.position()<< std::endl;
 LogDebug("FastCalorimetry")<<"Global Momentum calculates using SHP on HCAL: "<<  shsDest.momentum()<< std::endl;


//Initialization the momentum, position charge and energy on SHS

  GlobalPoint gPos(shsDest.position().x(),shsDest.position().y(),shsDest.position().z());
  GlobalVector gMom(shsDest.momentum().x(),shsDest.momentum().y(),shsDest.momentum().z());
//Muon Mass
       double mu = 0.1056583692;
//initial muon energy
       double ini_en = std::sqrt(startingMomentum.mag2()+mu*mu);
       // std::cout <<" ini_en_hcal:  "<< ini_en_hcal << std::endl;
//Muon Energy
       double en = std::sqrt(gMom.mag2()+mu*mu);
       //   std::cout <<" en_hcal:  "<< en_hcal << std::endl;

//delta energy due SHP (Contribution of the flutuation of muon energy )
       double delta_energy = en - ini_en;
      LogDebug("FastCalorimetry")<<"Flutuation of muon energy due SHP (delta_energy):  "<<  fabs(delta_energy) << std::endl;
     
      // std::cout <<"Flutuation of muon energy due SHP om HCAL (delta_energy):  "<<  fabs(delta_energy) << std::endl;
     
      XYZTLorentzVector trackPosition_shsDest(gPos.x(),gPos.y(),gPos.z(),0.);
      XYZTLorentzVector moment_shsDest(gMom.x(),gMom.y(),gMom.z(),en);
      float charge = (float)(shsDest.charge());

      // The energy loss simulator
      // float charge = (float)(myTrack.charge());
      // ParticlePropagator theMuon(moment,trackPosition,charge,0);

      ParticlePropagator theMuon(moment_shsDest,trackPosition_shsDest,charge,0);

      theMuon.setID(-(int)charge*13);
      energyLossHCAL->updateState(theMuon, segmentSizeinX0);

      //include of the correction enrgy due effect of Mag Field´on Muon Trajectory on HCAL

      // mipenergy = energyLossHCAL->deltaMom().E();
      mipenergy = (energyLossHCAL->deltaMom().E() + fabs(delta_energy));
      moment -= energyLossHCAL->deltaMom();
      // std::cout<< "Energy loss ONLY on HCAL: "<< energyLossHCAL->deltaMom().E()<< std::endl;	
      //Implementation of Muon Brem in HCAL (only applied if also energyLossHCAL is switched on)
      if (muonBremHCAL){
	muonBremHCAL->updateState(theMuon, segmentSizeinX0);
	mipenergy += muonBremHCAL->deltaP_BremPhoton().E();
	moment -=  muonBremHCAL->deltaP_BremMuon();


	//Testing
	XYZTLorentzVector thePhoton4Momentum(muonBremHCAL->deltaP_BremPhoton().Px(),
					     muonBremHCAL->deltaP_BremPhoton().Py(),
					     muonBremHCAL ->deltaP_BremPhoton().Pz(),
					     muonBremHCAL->deltaP_BremPhoton().E());
    


	if(muonBremHCAL->deltaP_BremPhoton().E()>0.){

          // Add a photon
          RawParticle thePhoton(22,thePhoton4Momentum);
	  //  thePhoton.setT(0.1);
          thePhoton.setID(22);
	  //Using new position of muon with brem photon vertex (check it) 
	
          XYZTLorentzVector vertex (myTrack.vertex().position().x(),
				    myTrack.vertex().position().y(),
				    myTrack.vertex().position().z(),
				    myTrack.vertex().position().t());



          thePhoton.setVertex(vertex);

	  _theUpdatedState.push_back(thePhoton);
        


      //Creating Hits of the photon 
         int ivertex = mySimEvent->addSimVertex(vertex,-1,FSimVertexType::BREM_VERTEX);
         mySimEvent->addSimTrack(&thePhoton,ivertex);

         int fsimi = ((int) mySimEvent->nTracks()-1);
         FSimTrack& myTrack = mySimEvent->track(fsimi);
         myTrack.setHcal(thePhoton,onHCAL);

    ///////////////////////////////////////////////////////
         // Photons : create an e+e- pair
         if( myTrack.noEndVertex() ) {

           // Simulate energy smearing for photon and electrons
           if (abs( myTrack.type()) == 11 ||  myTrack.type() == 22 ) {
             LogDebug("FastCalorimetry") << "Using EMShower by muons on HCAL:" <<std::endl;
             if ( myTrack.onHcal())
               EMShowerSimulation(myTrack);
       //      else if ( myTrack.onVFcal())
       //       LogDebug("FastCalorimetry") <<  "Using EMShower by muons on VFCAL" << std::endl;
       //      reconstructHCAL(myTrack);

           }
         }



	    if( debug_ ) {
	      float muon_mip =  energyLossHCAL->deltaMom().E();
	      float muon_brem = muonBremHCAL->deltaP_BremPhoton().E();
	      if ( muon_brem > 0.){
		LogDebug("FastCalorimetry")<< "Energy loss ONLY by Muon Brem on HCAL: "<< muon_brem<< std::endl;
		LogDebug("FastCalorimetry")<< "Energy loss ONLY by MIP HCAL: "<<  muon_mip << std::endl;
		LogDebug("FastCalorimetry")<< "Photon momentum:"<< muonBremHCAL->deltaP_BremPhoton().P()<<" Photon Energy: "<< muonBremHCAL->deltaP_BremPhoton().E() << std::endl;
		LogDebug("FastCalorimetry")<< "Contribution of Muon Brem by MUON Momentum on HCAL : "<< muonBremHCAL->deltaP_BremMuon().P() << std::endl;
		LogDebug("FastCalorimetry")<< "Total Energy loss on HCAL by (mip + Muon brem): "<<  mipenergy << std::endl;
		LogDebug("FastCalorimetry")<< "Contribution of MIP  by MUON Momentum on HCAL: "<< energyLossHCAL->deltaMom().P() << std::endl;
		LogDebug("FastCalorimetry")<< "Muon Momentum after the contribution of ( mip + Muon brem) on HCAL: "<<   moment.P() << " Muon Energy: "<< moment.E() << std::endl;
	      }
	    }

        }
	  }
          else{
	    LogDebug("FastCalorimetry")<<"Using only MIP for MUON"<<std::endl;
	    LogDebug("FastCalorimetry")<< "Muon moment using ONLY the contribution of MIP "<<   moment <<"moment:"<< energyLossHCAL->deltaMom() << std::endl; 
	  }
      
          myHcalHitMaker.setSpotEnergy(mipenergy);
	  myHcalHitMaker.addHit(segments[iseg].entrance());
	}

 

	if(segments[iseg].material()==CaloSegment::HCAL)
	  {
	    ilastHcal=iseg;
	  }
      }
  }
  if(ilastHcal>=0) // Update hcalExit position from 'dummy' HCAL entrance to 'straight line' HCAL exit
    hcalExit=segments[ilastHcal].exit();

  //////
  /////
  ////
  ///
  //



  // Copy the muon SimTrack (for Energy loss and Muon Brem)
  FSimTrack muonTrack(myTrack);
  if(energyLossHCAL || muonBremHCAL ) {
    muonTrack.setTkPosition(hcalExit);
    muonTrack.setTkMomentum(moment);
  } else if(energyLossECAL || muonBremECAL ) {
    muonTrack.setTkPosition(ecalExit);
    muonTrack.setTkMomentum(moment);
  } // else just leave tracker surface position and momentum...  

  muonSimTracks.push_back(muonTrack);


  // no need to change below this line
  std::map<unsigned,float>::const_iterator mapitr;
  std::map<unsigned,float>::const_iterator endmapitr;
  if(myTrack.onEcal() > 0) {
    // Save ECAL hits 
    endmapitr=myGrid.getHits().end();
    for(mapitr=myGrid.getHits().begin(); mapitr!=endmapitr; ++mapitr) {
      double energy = mapitr->second;
      if(onECAL==1)
	{
	  updateMap(EBDetId(mapitr->first).hashedIndex(),energy,myTrack.id(),EBMapping_,firedCellsEB_);
	}      
      else if(onECAL==2)
	{
	  updateMap(EEDetId(mapitr->first).hashedIndex(),energy,myTrack.id(),EEMapping_,firedCellsEE_);
	}
      
      if(debug_)
	LogDebug("FastCalorimetry") << " ECAL cell " << mapitr->first << " added,  E = " 
				    << energy << std::endl;  

//if (energy>0.) std::cout<< " ECAL cell " << mapitr->first << " added,  E = "
//                                    << energy << std::endl;



    }
  }
      
  // Save HCAL hits
  endmapitr=myHcalHitMaker.getHits().end();
  for(mapitr=myHcalHitMaker.getHits().begin(); mapitr!=endmapitr; ++mapitr) {
    double energy = mapitr->second;
    {
      updateMap(HcalDetId(mapitr->first).hashed_index(),energy,myTrack.id(),HMapping_,firedCellsHCAL_);
    }
    if(debug_)
      LogDebug("FastCalorimetry") << " HCAL cell "
				  << mapitr->first << " added    E = "
				  << mapitr->second << std::endl;


 //if (energy>0.)std::cout << " HCAL cell "<< mapitr->first << " added    E = "<< mapitr->second << std::endl;


 }
  
  if(debug_)
    LogDebug("FastCalorimetry") << std::endl << " FASTEnergyReconstructor::MipShowerSimulation  finished "
	 << std::endl;
}


void CalorimetryManager::readParameters(const edm::ParameterSet& fastCalo) {

  edm::ParameterSet ECALparameters = fastCalo.getParameter<edm::ParameterSet>("ECAL");

  evtsToDebug_ = fastCalo.getUntrackedParameter<std::vector<unsigned int> >("EvtsToDebug",std::vector<unsigned>());
  debug_ = fastCalo.getUntrackedParameter<bool>("Debug");

  gridSize_ = ECALparameters.getParameter<int>("GridSize");
  spotFraction_ = ECALparameters.getParameter<double>("SpotFraction");
  pulledPadSurvivalProbability_ = ECALparameters.getParameter<double>("FrontLeakageProbability");
  crackPadSurvivalProbability_ = ECALparameters.getParameter<double>("GapLossProbability");
  theCoreIntervals_ = ECALparameters.getParameter<std::vector<double> >("CoreIntervals");
  theTailIntervals_ = ECALparameters.getParameter<std::vector<double> >("TailIntervals");
  
  RCFactor_ = ECALparameters.getParameter<double>("RCFactor");
  RTFactor_ = ECALparameters.getParameter<double>("RTFactor");
  radiusFactor_ = ECALparameters.getParameter<double>("RadiusFactor");
  radiusPreshowerCorrections_ = ECALparameters.getParameter<std::vector<double> >("RadiusPreshowerCorrections");
  mipValues_ = ECALparameters.getParameter<std::vector<double> >("MipsinGeV");
  simulatePreshower_ = ECALparameters.getParameter<bool>("SimulatePreshower");

  if(gridSize_ <1) gridSize_= 7;
  if(pulledPadSurvivalProbability_ <0. || pulledPadSurvivalProbability_>1 ) pulledPadSurvivalProbability_= 1.;
  if(crackPadSurvivalProbability_ <0. || crackPadSurvivalProbability_>1 ) crackPadSurvivalProbability_= 0.9;
  
  LogInfo("FastCalorimetry") << " Fast ECAL simulation parameters " << std::endl;
  LogInfo("FastCalorimetry") << " =============================== " << std::endl;
  if(simulatePreshower_)
    LogInfo("FastCalorimetry") << " The preshower is present " << std::endl;
  else
    LogInfo("FastCalorimetry") << " The preshower is NOT present " << std::endl;
  LogInfo("FastCalorimetry") << " Grid Size : " << gridSize_  << std::endl; 
  if(spotFraction_>0.) 
    LogInfo("FastCalorimetry") << " Spot Fraction : " << spotFraction_ << std::endl;
  else
    {
      LogInfo("FastCalorimetry") << " Core of the shower " << std::endl;
      for(unsigned ir=0; ir < theCoreIntervals_.size()/2;++ir)
	{
	  LogInfo("FastCalorimetry") << " r < " << theCoreIntervals_[ir*2] << " R_M : " << theCoreIntervals_[ir*2+1] << "        ";
	}
      LogInfo("FastCalorimetry") << std::endl;
	
      LogInfo("FastCalorimetry") << " Tail of the shower " << std::endl;
      for(unsigned ir=0; ir < theTailIntervals_.size()/2;++ir)
	{
	  LogInfo("FastCalorimetry") << " r < " << theTailIntervals_[ir*2] << " R_M : " << theTailIntervals_[ir*2+1] << "        ";
	}
      LogInfo("FastCalotimetry") << "Radius correction factor " << radiusFactor_ << std::endl;
      LogInfo("FastCalorimetry") << std::endl;
      if(mipValues_.size()>2) 	{
	LogInfo("FastCalorimetry") << "Improper number of parameters for the preshower ; using 95keV" << std::endl;
	mipValues_.clear();
	mipValues_.resize(2,0.000095);
	}
    }

  LogInfo("FastCalorimetry") << " FrontLeakageProbability : " << pulledPadSurvivalProbability_ << std::endl;
  LogInfo("FastCalorimetry") << " GapLossProbability : " << crackPadSurvivalProbability_ << std::endl;

  
  // RespCorrP: p (momentum), ECAL and HCAL corrections = f(p)
  edm::ParameterSet CalorimeterParam = fastCalo.getParameter<edm::ParameterSet>("CalorimeterProperties");

  rsp = CalorimeterParam.getParameter<std::vector<double> >("RespCorrP");
   LogInfo("FastCalorimetry") << " RespCorrP (rsp) size " << rsp.size() << std::endl;

  if( rsp.size()%3 !=0 )  {
    LogInfo("FastCalorimetry") 
      << " RespCorrP size is wrong -> no corrections applied !!!" 
      << std::endl;

      p_knots.push_back(14000.);
      k_e.push_back    (1.);
      k_h.push_back    (1.);
  }
  else {
    for(unsigned i = 0; i < rsp.size(); i += 3) { 
     LogInfo("FastCalorimetry") << "i = " << i/3 << "   p = " << rsp [i] 
				<< "   k_e(p) = " << rsp[i+1] 
				<< "   k_e(p) = " << rsp[i+2] << std::endl; 
      
      p_knots.push_back(rsp[i]);
      k_e.push_back    (rsp[i+1]);
      k_h.push_back    (rsp[i+2]); 
    }
  }  
 

  //FR
  edm::ParameterSet HCALparameters = fastCalo.getParameter<edm::ParameterSet>("HCAL");
  optionHDSim_ = HCALparameters.getParameter<int>("SimOption");
  hdGridSize_  = HCALparameters.getParameter<int>("GridSize");
  hdSimMethod_ = HCALparameters.getParameter<int>("SimMethod");
  //RF

  unfoldedMode_ = fastCalo.getUntrackedParameter<bool>("UnfoldedMode",false);
}


void CalorimetryManager::updateMap(uint32_t cellid,float energy,int id,std::map<uint32_t,std::vector<std::pair<int,float> > > & mymap)
{
  //  std::cout << " updateMap " << std::endl;
  std::map<unsigned,std::vector<std::pair<int,float> > >::iterator cellitr;
  cellitr = mymap.find(cellid);
  if(!unfoldedMode_) id=0;
  if( cellitr==mymap.end())
    {      
      std::vector<std::pair<int,float> > myElement;
      myElement.push_back(std::pair<int,float> (id,energy));
      mymap[cellid]=myElement;
    }
  else
    {
      if(!unfoldedMode_)
	{
	  cellitr->second[0].second+=energy;
	}
      else
	cellitr->second.push_back(std::pair<int,float>(id,energy));
    }
}

void CalorimetryManager::updateMap(int hi,float energy,int tid,std::vector<std::vector<std::pair<int,float> > > & mymap, std::vector<int>& firedCells)
{
  // Standard case first : one entry per cell 
  if(!unfoldedMode_)
    {
      // if new entry, update the list 
      if(mymap[hi].size()==0)
	{
	  firedCells.push_back(hi);
	  mymap[hi].push_back(std::pair<int,float>(0,energy));
	}
      else
	mymap[hi][0].second+=energy;
    }
  else
    {
      //      std::cout << "update map " << mymap[hi].size() << " " << hi << std::setw(8) << std::setprecision(6) <<  energy ;
      //      std::cout << " " << mymap[hi][0].second << std::endl;
      // the minimal size is always 1 ; in this case, no push_back 
      if(mymap[hi].size()==0)
	{
	  //	  if(tid==0) std::cout << " Shit ! " << std::endl;
	  firedCells.push_back(hi);
	}

      mymap[hi].push_back(std::pair<int,float>(tid,energy));
    }
  
}


void CalorimetryManager::respCorr(double p) {

  int sizeP = p_knots.size();

  if(sizeP <= 1) {
    ecorr = 1.;
    hcorr = 1.;
  }
  else {
    int ip = -1;    
    for (int i = 0; i < sizeP; i++) { 
      if (p < p_knots[i]) { ip = i; break;}
    }
    if (ip == 0) {
      ecorr = k_e[0];
      hcorr = k_h[0];
    }
    else {
      if(ip == -1) {
	ecorr = k_e[sizeP-1];
	hcorr = k_h[sizeP-1];
      } 
      else {
	double x1 =  p_knots[ip-1];
	double x2 =  p_knots[ip];
	double y1 =  k_e[ip-1];
	double y2 =  k_e[ip];
	
	if(x1 == x2) {
	  //        std::cout << " equal p_knots values!!! " << std::endl;
	}	
      
	ecorr = (y1 + (y2 - y1) * (p - x1)/(x2 - x1));
	
	y1 =  k_h[ip-1];
	y2 =  k_h[ip];
	hcorr = (y1 + (y2 - y1) * (p - x1)/(x2 - x1)); 
	
      }
    }
  }

  if(debug_)
    LogDebug("FastCalorimetry") << " p, ecorr, hcorr = " << p << " "  
			        << ecorr << "  " << hcorr << std::endl;
	
}


void CalorimetryManager::loadFromEcalBarrel(edm::PCaloHitContainer & c) const
{ 
  unsigned size=firedCellsEB_.size();
  //  float sum=0.;
  for(unsigned ic=0;ic<size;++ic)
    {
      int hi=firedCellsEB_[ic];
      if(!unfoldedMode_)
	{
	  c.push_back(PCaloHit(EBDetId::unhashIndex(hi),EBMapping_[hi][0].second,0.,0));
	  //	  std::cout << "Adding " << hi << " " << EBDetId::unhashIndex(hi) << " " ;
	  //	  std::cout << EBMapping_[hi][0].second << " " << EBMapping_[hi][0].first << std::endl;
	}
      else
	{
	  unsigned npart=EBMapping_[hi].size();
	  for(unsigned ip=0;ip<npart;++ip)
	    {
	      c.push_back(PCaloHit(EBDetId::unhashIndex(hi),EBMapping_[hi][ip].second,0.,
				   EBMapping_[hi][ip].first));

	    }
	}
	
      //      sum+=cellit->second;
    }
  
//  for(unsigned ic=0;ic<61200;++ic) 
//    { 
//      EBDetId myCell(EBDetId::unhashIndex(ic)); 
//      if(!myCell.null()) 
//        { 
//	  float total=0.;
//	  for(unsigned id=0;id<EBMapping_[ic].size();++id)
//	    total+=EBMapping_[ic][id].second;
//	  if(EBMapping_[ic].size()>0)
//	    std::cout << "Adding " << ic << " " << myCell << " " << std::setprecision(8) <<total << std::endl; 
//        } 
//    } 


  //  std::cout << " SUM : " << sum << std::endl;
  //  std::cout << " Added " <<c.size() << " hits " <<std::endl;
}


void CalorimetryManager::loadFromEcalEndcap(edm::PCaloHitContainer & c) const
{
  unsigned size=firedCellsEE_.size();
  //  float sum=0.;
  for(unsigned ic=0;ic<size;++ic)
    {
      int hi=firedCellsEE_[ic];
      if(!unfoldedMode_)
	c.push_back(PCaloHit(EEDetId::unhashIndex(hi),EEMapping_[hi][0].second,0.,0));
      else
	{
	  unsigned npart=EEMapping_[hi].size();
	  for(unsigned ip=0;ip<npart;++ip)
	    c.push_back(PCaloHit(EEDetId::unhashIndex(hi),EEMapping_[hi][ip].second,0.,
				 EEMapping_[hi][ip].first));
	}
	
      //      sum+=cellit->second;
    }
  //  std::cout << " SUM : " << sum << std::endl;
  //  std::cout << " Added " <<c.size() << " hits " <<std::endl;
}

void CalorimetryManager::loadFromHcal(edm::PCaloHitContainer & c) const
{
  unsigned size=firedCellsHCAL_.size();
  //  float sum=0.;
  for(unsigned ic=0;ic<size;++ic)
    {
      int hi=firedCellsHCAL_[ic];
      if(!unfoldedMode_)
	c.push_back(PCaloHit(theDetIds_[hi],HMapping_[hi][0].second,0.,0));
      else
	{
	  unsigned npart=HMapping_[hi].size();
	  for(unsigned ip=0;ip<npart;++ip)
	    c.push_back(PCaloHit(theDetIds_[hi],HMapping_[hi][ip].second,0.,
				 HMapping_[hi][ip].first));
	}
	
      //      sum+=cellit->second;
    }
  //  std::cout << " SUM : " << sum << std::endl;
  //  std::cout << " Added " <<c.size() << " hits " <<std::endl;
}

void CalorimetryManager::loadFromPreshower(edm::PCaloHitContainer & c) const
{
  std::map<uint32_t,std::vector<std::pair< int,float> > >::const_iterator cellit;
  std::map<uint32_t,std::vector<std::pair <int,float> > >::const_iterator preshEnd=ESMapping_.end();
  
  for(cellit=ESMapping_.begin();cellit!=preshEnd;++cellit)
    {
      if(!unfoldedMode_)	
	c.push_back(PCaloHit(cellit->first,cellit->second[0].second,0.,0));
      else
	{
	  unsigned npart=cellit->second.size();
	  for(unsigned ip=0;ip<npart;++ip)
	    {
	      c.push_back(PCaloHit(cellit->first,cellit->second[ip].second,0.,cellit->second[ip].first));
	    }
	}
    }
}


// The main danger in this method is to screw up to relationships between particles
// So, the muon FSimTracks created by FSimEvent.cc are simply to be updated 
void CalorimetryManager::loadMuonSimTracks(edm::SimTrackContainer &muons) const
{
  unsigned size=muons.size();
  for(unsigned i=0; i<size;++i)
    {
      int id=muons[i].trackId();
      if(abs(muons[i].type())!=13) continue;
      // identify the corresponding muon in the local collection

      std::vector<FSimTrack>::const_iterator itcheck=find_if(muonSimTracks.begin(),muonSimTracks.end(),FSimTrackEqual(id));
      if(itcheck!=muonSimTracks.end())
	{
	  muons[i].setTkPosition(itcheck->trackerSurfacePosition());
	  muons[i].setTkMomentum(itcheck->trackerSurfaceMomentum());
//	  std::cout << " Found the SimTrack " << std::endl;
//	  std::cout << *itcheck << std::endl;
//	  std::cout << "SimTrack Id "<< id << " " << muons[i] << " " << std::endl;
	}
//      else
//	{
//	  std::cout << " Calorimetery Manager : this should really not happen " << std::endl;
//	  std::cout << " Was looking for " << id << " " << muons[i] << std::endl;
//	  for(unsigned i=0;i<muonSimTracks.size();++i)
//	    std::cout << muonSimTracks[i] << std::endl;
//	}
    }

}


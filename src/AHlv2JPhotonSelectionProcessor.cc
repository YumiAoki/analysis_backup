// *****************************************************
// e+e- ------> AH ------> gamma + X
// Processor for photon selection
//                        ----Junping
// *****************************************************
#include "AHlv2JPhotonSelectionProcessor.h"
#include <iostream>
#include <sstream>
#include <iomanip>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/Cluster.h>
#include <UTIL/LCTypedVector.h>
#include <EVENT/Track.h>
#include <UTIL/LCRelationNavigator.h>
#include <EVENT/ParticleID.h>
#include <marlin/Exceptions.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "Utilities.h"

using namespace lcio ;
using namespace marlin ;
using namespace std;

using namespace mylib;

AHlv2JPhotonSelectionProcessor aAHlv2JPhotonSelectionProcessor ;


AHlv2JPhotonSelectionProcessor::AHlv2JPhotonSelectionProcessor() : Processor("AHlv2JPhotonSelectionProcessor") {
  
  // modify processor description
  _description = "AHlv2JPhotonSelectionProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticlesSkimmed") ) ;

  registerInputCollection( LCIO::LCRELATION,
			   "InputMCTruthLinkCollection" , 
			   "Name of the MCTruthLink collection"  ,
			   _colMCTL ,
			   std::string("RecoMCTruthLink") ) ;

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "InputPandoraPFOsCollection" , 
			   "Name of the PandoraPFOs collection"  ,
			   _colPFOs ,
			   std::string("PandoraPFOs") ) ;

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputNewPFOsCollection",
			    "Name of the new PFOs collection without isolated lepton ",
			    _colNewPFOs,
			    std::string("newPandoraPFOs") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "OutputPhotonsCollection",
			    "Name of collection with the selected photons",
			    _colPhotons,
			    std::string("photons") );

}

void AHlv2JPhotonSelectionProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStat = 0;
  
}

void AHlv2JPhotonSelectionProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void AHlv2JPhotonSelectionProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

#if 0
  Double_t fPhotonEnergyCut = 100.;  // energy cut for hard photon seed
  Double_t fRecoilMassCut = 600.;  // recoil mass cut
#else  
  Double_t fPhotonEnergyCut = 50.;  // energy cut for hard photon seed
  Double_t fRecoilMassCut = 0.;  // recoil mass cut
#endif
  //  Double_t fRecoilMassCut = 0.;  // recoil mass cut   
  Double_t fCosConeCut = 0.998;   // the angle of cone around the direction of photon for merging
  Double_t fCosSmallCone = 0.98;   // the angle of small cone around photon for isolation
  Double_t fCosLargeCone = 0.95;   // the angle of large cone around photon for isolation  

  const Double_t fEcm =250.;
  TLorentzVector lortzEcm = getLorentzEcm(fEcm);

  const Double_t massX = 125.;
  Double_t EPhoton = (fEcm*fEcm - massX*massX)/2./fEcm;

  // cut table
  if (!hStat) hStat = new TH1D("hStat", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  streamlog_out(DEBUG) << "Hello, Photon Selection!" << endl;


  LCCollectionVec *pNewPFOsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  LCCollectionVec *pPhotonsCollection = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
  pNewPFOsCollection->setSubset(true);
  
  // -- Read out PFO information --
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  if (!colPFO) {
    std::cerr << "No PFO Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "MCParticle and PandoraPFOs Collections found!" << ends;

  std::vector<lcio::ReconstructedParticle*> photons;
  std::vector<lcio::ReconstructedParticle*> newPFOs;
  
  Int_t nPFOs = colPFO->getNumberOfElements();

  ReconstructedParticle *thePhoton = 0;
  Double_t deltaEMin = 1000.;

  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *recPart = dynamic_cast<ReconstructedParticle*>(colPFO->getElementAt(i));
    newPFOs.push_back(recPart);
    Double_t energy = recPart->getEnergy();
    Double_t charge = recPart->getCharge();
    TVector3 momentum = TVector3(recPart->getMomentum());
    //    Double_t momentumMagnitude = momentum.Mag();
    TLorentzVector lortz = TLorentzVector(momentum,energy);
    if (TMath::Abs(charge) < 0.5 && recPart->getType() == 22) {
      photons.push_back(recPart);
      if (energy > fPhotonEnergyCut) {
	Double_t mrecoil = getRecoilMass(lortzEcm,lortz);
	if (mrecoil > fRecoilMassCut) {
	  Double_t deltaE = TMath::Abs(energy-EPhoton);
	  if (deltaE < deltaEMin) {
	    deltaEMin = deltaE;
	    thePhoton = recPart;
	  }
	}
      }
    }
  }

  if (!thePhoton)  throw marlin::SkipEventException(this);
  hStat->Fill(++selid);
  //  gCutName[(Int_t)selid] << "E_Photon > 100 && M_Recoil > 600" << ends;
  gCutName[(Int_t)selid] << "E_Photon > 50 && M_Recoil > 0" << ends;  

  // recover splited photons
  std::vector<lcio::ReconstructedParticle*> photonsMerged;
  ReconstructedParticleImpl * recoPhoton = new ReconstructedParticleImpl();
  TLorentzVector lortzPhoton = TLorentzVector(thePhoton->getMomentum(),thePhoton->getEnergy());
  TLorentzVector lortzPhotonOrig = lortzPhoton;
  recoPhoton->addParticle(thePhoton);
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator ib=photons.begin();ib<photons.end()-1;ib++) {
    if (*ib == thePhoton) continue;
    Double_t cosTheta = getCosTheta(*ib,thePhoton);
    if (cosTheta < fCosConeCut) continue;
    recoPhoton->addParticle(*ib);
    lortzPhoton += TLorentzVector((*ib)->getMomentum(),(*ib)->getEnergy());
    photonsMerged.push_back(*ib);
  }
  Double_t energy = lortzPhoton.E();
  Double_t mass   = lortzPhoton.M();
  Double_t momentum[3] = {lortzPhoton.Px(),lortzPhoton.Py(),lortzPhoton.Pz()};
  Double_t charge = thePhoton->getCharge();
  recoPhoton->setMomentum(momentum);
  recoPhoton->setEnergy(energy);
  recoPhoton->setMass(mass);
  recoPhoton->setCharge(charge);
  recoPhoton->setType(94);
  pPhotonsCollection->addElement(recoPhoton);
  
  // get photon cone energies
  Bool_t woFSR = kFALSE;
  Double_t coneEnergy0[3] = {0.,0.,0.};
  Double_t pFSR[4] = {0.,0.,0.,0.};
  Double_t pLargeCone[4]  = {0.,0.,0.,0.};
  Int_t nConePhoton = 0;
  getConeEnergy(thePhoton,colPFO,fCosSmallCone,woFSR,coneEnergy0,pFSR,fCosLargeCone,pLargeCone,nConePhoton);
  //      Double_t coneEnergy = coneEnergy0[0];
  Double_t coneEN     = coneEnergy0[1];
  Double_t coneEC     = coneEnergy0[2];
  TLorentzVector lortzLargeCone = TLorentzVector(pLargeCone[0],pLargeCone[1],pLargeCone[2],pLargeCone[3]);
  //  TVector3 momentumLargeCone = lortzLargeCone.Vect();
  Double_t cosThetaWithLargeCone = getCosTheta(lortzPhotonOrig,lortzLargeCone);
  Double_t energyRatioWithLargeCone = lortzPhotonOrig.E()/(lortzPhotonOrig.E()+lortzLargeCone.E());
  
  // save other PFOs to a new collection
  for (std::vector<lcio::ReconstructedParticle*>::const_iterator iObj=newPFOs.begin();iObj<newPFOs.end();++iObj) {
    if ((*iObj) == thePhoton) continue;
    if (isFoundInVector((*iObj),photonsMerged)) continue;
    pNewPFOsCollection->addElement(*iObj);
  }

  // save collections
  evt->addCollection(pNewPFOsCollection,_colNewPFOs.c_str());
  evt->addCollection(pPhotonsCollection,_colPhotons.c_str());


  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  last->cd();
}



void AHlv2JPhotonSelectionProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AHlv2JPhotonSelectionProcessor::end(){ 

  cerr << "AHlv2JPhotonSelectionProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  
  cerr << "  =============" << endl;
  cerr << "   Cut Summary " << endl;
  cerr << "  =============" << endl;
  cerr << "   ll+X    " << endl;
  cerr << "  =============" << endl;
  cerr << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;

}

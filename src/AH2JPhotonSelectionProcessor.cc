// *****************************************************
// e+e- ------> AH ------> gamma + X
// Processor for photon selection
//                        ----Junping
// *****************************************************
#include "AH2JPhotonSelectionProcessor.h"
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

AH2JPhotonSelectionProcessor aAH2JPhotonSelectionProcessor ;


AH2JPhotonSelectionProcessor::AH2JPhotonSelectionProcessor() : Processor("AH2JPhotonSelectionProcessor") {
  
  // modify processor description
  _description = "AH2JPhotonSelectionProcessor does whatever it does ..." ;
  

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

void AH2JPhotonSelectionProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStat = 0;
  
}

void AH2JPhotonSelectionProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void AH2JPhotonSelectionProcessor::processEvent( LCEvent * evt ) { 

    
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

  // cut table
  if (!hStat) hStat = new TH1D("hStat", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  streamlog_out(DEBUG) << "Hello, Photon Selection!" << endl;

  static TNtupleD *hGen = 0;
  if (!hGen) {
    stringstream tupstr_gen;
    tupstr_gen << "ea:mx:mr:eisr1:eisr2:"
		<< "e0:px0:py0:pz0:"
		<< "mff"
	       << ends;
    hGen = new TNtupleD("hGen","",tupstr_gen.str().data());
  }

  // -- Get the MCTruth Linker --
  LCCollection *colMCTL = evt->getCollection(_colMCTL);
  //  LCRelationNavigator *navMCTL = new LCRelationNavigator(colMCTL);


  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nMCP = colMC->getNumberOfElements();
  TLorentzVector lortzGam1MC,lortzGam2MC,lortzGam3MC,lortzXMC,lortzRecMC;
  TLorentzVector lortzf1MC, lortzf2MC;
  TLorentzVector lortzISR1MC, lortzISR2MC;
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    //    Int_t pdg = mcPart->getPDG();
    Double_t energy = mcPart->getEnergy();
    TVector3 pv = TVector3(mcPart->getMomentum());
    TLorentzVector lortz = TLorentzVector(pv,energy);
    //    Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
   #if 0
    if (i == 1){
      lortzGam1MC = lortz;
    }
   /* if (i == 1) {
      lortzGam2MC = lortz;
      lortzXMC += lortz;
    }
    if (i == 0) {
      lortzGam3MC = lortz;
      lortzXMC += lortz;
    }*/
    if (i == 2) {
      lortzISR1MC = lortz;
    }
    if (i == 3) {
      lortzISR2MC = lortz;
    }
    if (i == 0) {
      lortzXMC = lortz;
   
	 }
#else
    if (i == 0) {
      lortzISR1MC = lortz;
    }
    if (i == 1) {
      lortzISR2MC = lortz;
    }
    if (i == 2) {
      lortzf1MC = lortz;
    }
    if (i == 3) {
      lortzf2MC = lortz;
    }
#endif  
}

#if 0
  TLorentzVector lortzFinal = lortzGam1MC + lortzISR1MC + lortzISR2MC + lortzXMC;
#else 

  TLorentzVector lortzFinal = lortzISR1MC + lortzISR2MC + lortzf1MC + lortzf2MC;
#endif
 TLorentzVector lortz2f = lortzf1MC + lortzf2MC;

  const Double_t fEcm =250.;
  TLorentzVector lortzEcm = getLorentzEcm(fEcm);

  const Double_t massX = 125.;
  Double_t EPhoton = (fEcm*fEcm - massX*massX)/2./fEcm;

  Double_t data_gen[50];
  data_gen[ 0]= lortzGam1MC.E();
  //data_gen[ 1]= lortzGam2MC.E();
  //data_gen[ 2]= lortzGam3MC.E();  
  data_gen[ 1]= lortzXMC.M();
  data_gen[ 2]= getRecoilMass(lortzEcm,lortzGam1MC);
  data_gen[ 3]= lortzISR1MC.E();
  data_gen[ 4]= lortzISR2MC.E();
  data_gen[ 5]= lortzFinal.E();
  data_gen[ 6]= lortzFinal.Px();
  data_gen[ 7]= lortzFinal.Py();
  data_gen[ 8]= lortzFinal.Pz();
  data_gen[ 9]= lortz2f.M();
  hGen->Fill(data_gen);

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



void AH2JPhotonSelectionProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AH2JPhotonSelectionProcessor::end(){ 

  cerr << "AH2JPhotonSelectionProcessor::end()  " << name() 
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

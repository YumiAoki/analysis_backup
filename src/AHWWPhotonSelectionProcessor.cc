// *****************************************************
// e+e- ------> AH ------> gamma + X
// Processor for photon selection
//                        ----Junping
// *****************************************************
#include "AHWWPhotonSelectionProcessor.h"
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

AHWWPhotonSelectionProcessor aAHWWPhotonSelectionProcessor ;


AHWWPhotonSelectionProcessor::AHWWPhotonSelectionProcessor() : Processor("AHWWPhotonSelectionProcessor") {
  
  // modify processor description
  _description = "AHWWPhotonSelectionProcessor does whatever it does ..." ;
  

  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::MCPARTICLE,
			   "InputMCParticlesCollection" , 
			   "Name of the MCParticle collection"  ,
			   _colMCP ,
			   std::string("MCParticlesSkimmed") ) ;


}

void AHWWPhotonSelectionProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStatMC = 0;
  
}

void AHWWPhotonSelectionProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void AHWWPhotonSelectionProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;


  // cut table
  if (!hStatMC) hStatMC = new TH1D("hStatMC", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStatMC->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  streamlog_out(DEBUG) << "Hello, Photon Selection!" << endl;

  static TNtupleD *hGenW = 0;
  if (!hGenW) {
    stringstream tupstr_gen;
    tupstr_gen << "eisr:cosisr:mw1:mw2:mww:eisr1:eisr2"
	       << ends;
    hGenW = new TNtupleD("hGenW","",tupstr_gen.str().data());
  }


  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nMCP = colMC->getNumberOfElements();
  TLorentzVector lortzGam1MC,lortzGam2MC,lortzGam3MC,lortzXMC,lortzRecMC;
  TLorentzVector lortzW1MC, lortzW2MC;
  TLorentzVector lortzISR1MC, lortzISR2MC;
  for (Int_t i=0;i<nMCP;i++) {
    MCParticle *mcPart = dynamic_cast<MCParticle*>(colMC->getElementAt(i));
    //    Int_t pdg = mcPart->getPDG();
    Double_t energy = mcPart->getEnergy();
    TVector3 pv = TVector3(mcPart->getMomentum());
    TLorentzVector lortz = TLorentzVector(pv,energy);
    //    Int_t ioverlay = mcPart->isOverlay()? 1 : 0;
    if (i == 0) {
      lortzISR1MC = lortz;
    }
    if (i == 1) {
      lortzISR2MC = lortz;
    }
    if (i == 2 || i == 3) {
      lortzW1MC += lortz;
    }
    if (i == 4 || i == 5) {
      lortzW2MC += lortz;
    }
  }

TLorentzVector lortzISRMC = lortzISRMC.E() > lortzISRMC.E() ? lortzISR1MC: lortzISR2MC;
TLorentzVector lortzWWMC = lortzW1MC + lortzW2MC;

  Double_t data_gen[50];
  data_gen[ 0]= lortzISRMC.E();
  data_gen[ 1]= lortzISRMC.CosTheta();
  data_gen[ 2]= lortzW1MC.M();
  data_gen[ 3]= lortzW2MC.M();
  data_gen[ 4]= lortzWWMC.M();
  data_gen[ 5]= lortzISR1MC.E();
  data_gen[ 6]= lortzISR2MC.E();
 
 hGenW->Fill(data_gen);




  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  last->cd();
}



void AHWWPhotonSelectionProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AHWWPhotonSelectionProcessor::end(){ 

  cerr << "AHWWPhotonSelectionProcessor::end()  " << name() 
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
         << "  " << setw(10) << static_cast<int>(hStatMC->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;

}

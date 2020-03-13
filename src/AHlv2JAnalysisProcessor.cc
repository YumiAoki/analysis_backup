// *****************************************************
// e+e- ------> AH ------> gamma + X
// Processor for Analysis
//                        ----Junping
// *****************************************************
#include "AHlv2JAnalysisProcessor.h"
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
#include "UTIL/PIDHandler.h"

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

AHlv2JAnalysisProcessor aAHlv2JAnalysisProcessor ;


AHlv2JAnalysisProcessor::AHlv2JAnalysisProcessor() : Processor("AHlv2JAnalysisProcessor") {
  
  // modify processor description
  _description = "AHlv2JAnalysisProcessor does whatever it does ..." ;
  

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
			    "InputNewPFOsCollection",
			    "Name of the new PFOs collection after some pre-cuts",
			    _colNewPFOs,
			    std::string("newPandoraPFOs") );

  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			    "InputPhotonsCollection",
			    "Name of collection with the selected photons",
			    _colPhotons,
			    std::string("photons") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputLeptonsCollection",
			   "Name of collection with the selected leptons",
			   _colLeptons,
			   std::string("Leptons") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE, 
			   "InputJetsCollection",
			   "Name of the jets collection",
			   _colJets,
			   std::string("Refined_2Jet") );
}

void AHlv2JAnalysisProcessor::init() { 

  streamlog_out(DEBUG) << "   init called  " 
		       << std::endl ;
  
  streamlog_out(DEBUG) << "   init for h to ww called  " 
		       << std::endl ;
  
  // usually a good idea to
  printParameters() ;

  _nRun = 0 ;
  _nEvt = 0 ;

  hStatAnl = 0;

}

void AHlv2JAnalysisProcessor::processRunHeader( LCRunHeader* run) { 

  _nRun++ ;
} 

void AHlv2JAnalysisProcessor::processEvent( LCEvent * evt ) { 

    
  // this gets called for every event 
  // usually the working horse ...
  _nEvt++;

  Double_t kMassZ = 91.187;     // Z mass
  Double_t fMassZCut = 40.;     // mass cut for lepton pair from Z
  Double_t fEpsilon = 1.E-10;

  static const Double_t fEnergyCut = 0.05; // minimum energy of each pfo
  Double_t fCosConeCut = 0.998;   // the angle of cone around the direction of photon for merging
  Double_t fCosSmallCone = 0.98;   // the angle of small cone around photon for isolation
  Double_t fCosLargeCone = 0.95;   // the angle of large cone around photon for isolation  

  
  // cut table
  if (!hStatAnl) hStatAnl = new TH1D("hStatAnl", "Cut Table", 20, 0, 20);
  Double_t selid = -0.5;
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "No Cuts" << ends;

  TDirectory *last = gDirectory;
  gFile->cd("/");

  cerr << endl << "Hello, e+e- --> AH Analysis! Event No. " << _nEvt << endl;

  static TNtupleD *hAnl = 0;
  if (!hAnl) {
    stringstream tupstr_anl;
    tupstr_anl << "eamc:eisr1mc:eisr2mc:mxmc:mrmc:cosmc" << ":"
	       << "ea:mr:cos:eaorig:mrorig:cosorig"   << ":"
	       << "mcserial:ma"                       << ":"
	       << "coneen:coneec:coslarcon:energyratio" << ":"
	       << "evis"   << ":"
	       << "px1:py1:pz1:e1:px2:py2:pz2:e2:mpx:mpy:mpz:emis:mpt"                   << ":"
	       << "probj1:probj2:bmax1:bmax2"                                            << ":"
	       << "yminus:yplus:yminus4:yplus4"                                          << ":"
	       << "njets:npfos:npfos1:npfos2:npfosc:npfospt500:npfosc1:npfosc2"          << ":"
	       << "cosj1j2:cosaj:cosaj1:cosaj2:mH:cosH"   << ":"
	       << "pxa:pya:pza" <<":"
		<< "e0:px0:py0:pz0:"
		<< "mff:mw1rec:mw2rec:m2WMC:nwqq:nwlv:mw1mc:mw2mc" << ":" 
               << "pxlep:pylep:pzlep:elep:coslep"
	       << ends;
    hAnl = new TNtupleD("hAnl","",tupstr_anl.str().data());
  }

  // ------------------------------------------------
  // -- read out the MCParticles information
  // ------------------------------------------------
  LCCollection *colMC = evt->getCollection(_colMCP);
  if (!colMC) {
    std::cerr << "No MC Collection Found!" << std::endl;
    throw marlin::SkipEventException(this);
  }
  Int_t nMCP = colMC->getNumberOfElements();
  TLorentzVector lortzGam1MC,lortzGam2MC,lortzGam3MC,lortzXMC,lortzRecMC;
  TLorentzVector lortzISR1MC, lortzISR2MC, lortzf1MC, lortzf2MC;
  TLorentzVector lortzW1MC, lortzW2MC;
  Int_t nwqq = 0;
  Int_t nwlv = 0;
  double mwplusmc = 0;
  double mwminusmc = 0;
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
    if (i == 4 || i == 5) {
      MCParticleVec daughter = mcPart->getDaughters();
      
      MCParticle *daughter1 = 0;
      if (daughter.size() >0 ) daughter1 = daughter[0];
      
      int dpdg = 0;
      if ( daughter1 ) dpdg = daughter1->getPDG();
      if ( abs(dpdg)<=6){ 
          nwqq++;
      }
      if (abs(dpdg)>=11 && abs(dpdg)<=16) {
          nwlv++;
      }
        
    }
    if (i == 2) {
      lortzISR1MC = lortz;
    }
    if (i == 3) {
      lortzISR2MC = lortz;
    }

     double mw = mcPart->getMass();
      if ( i == 4){ 
          mwplusmc = mw;
          lortzW1MC = lortz;
      }
      if ( i == 5){ 
          mwminusmc = mw;
          lortzW2MC = lortz;
      }

#endif  
}

TLorentzVector lortz2WMC = lortzW1MC + lortzW2MC;
double m2WMC = lortz2WMC.M();

   /* if (i == 1) {
      lortzGam1MC = lortz;
    }
    if (i == 4) {
      lortzGam2MC = lortz;
      lortzXMC += lortz;
    }
    if (i == 5) {
      lortzGam3MC = lortz;
      lortzXMC += lortz;
    }
    if (i == 2) {
      lortzISR1MC = lortz;
    }
    if (i == 3) {
      lortzISR2MC = lortz;
    }
  }*/

#if 0
  TLorentzVector lortzFinal = lortzGam1MC + lortzISR1MC + lortzISR2MC + lortzXMC;
#else 

  TLorentzVector lortzFinal = lortzISR1MC + lortzISR2MC + lortzf1MC + lortzf2MC;
#endif
 TLorentzVector lortz2f = lortzf1MC + lortzf2MC;
  // ------------------------------------------------
  // -- read out the Thrust information
  // ------------------------------------------------
#if 0
  LCCollection *colSelRecPart = evt->getCollection("SelectedReconstructedParticle");
  Double_t principleThrust = colSelRecPart->parameters().getFloatVal("principleThrustValue");
  Double_t majorThrust = colSelRecPart->parameters().getFloatVal("majorThrustValue");
  Double_t minorThrust = colSelRecPart->parameters().getFloatVal("minorThrustValue");
  FloatVec tAxis;
  FloatVec thrustAxis = colSelRecPart->parameters().getFloatVals("principleThrustAxis",tAxis);
  TVector3 principleAxis = TVector3(thrustAxis[0],thrustAxis[1],thrustAxis[2]);
  Double_t cosThrustAxis = principleAxis.CosTheta();
#endif

  // ------------------------------------------------
  // -- read out the NewPFOs information
  // ------------------------------------------------
  LCCollection *colPFO = evt->getCollection(_colPFOs);
  LCCollection *colNewPFO = evt->getCollection(_colNewPFOs);
  if (!colNewPFO) {
    cerr << "No NewPFOs Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "NewPFOs Collection Found" << ends;
  // get the visible energy
  Int_t nPFOs = colNewPFO->getNumberOfElements();
  TLorentzVector lortzVis = TLorentzVector(0.,0.,0.,0.);
  Int_t nParticles = 0, nParticlesC = 0, nParticlesPt500 = 0;
  for (Int_t i=0;i<nPFOs;i++) {
    ReconstructedParticle *pfo = dynamic_cast<ReconstructedParticle*>(colNewPFO->getElementAt(i));
    TLorentzVector lortz = getLorentzVector(pfo);
    lortzVis += lortz;
    if (pfo->getEnergy() >= fEnergyCut) nParticles++;
    if (pfo->getEnergy() >= fEnergyCut && pfo->getCharge() != 0) nParticlesC++;
    if (lortz.Pt() > 0.5) nParticlesPt500++;
  }

  // ------------------------------------------------
  // -- read out the MCTruthLink information
  // ------------------------------------------------
  LCCollection *colMCTL = evt->getCollection(_colMCTL);

  // ------------------------------------------------
  // -- read out the Leptons information
  // ------------------------------------------------
  LCCollection *colLep = evt->getCollection(_colLeptons);
  Int_t nLeps = colLep->getNumberOfElements();
  if (!colLep || nLeps != 1) {
    cerr << "No Leptons Collection Found or nLeptons != 1 !" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "Leptons Collection Found" << ends;
  ReconstructedParticle *lepton = dynamic_cast<ReconstructedParticle*>(colLep->getElementAt(0));
  TLorentzVector lortzLepton = getLorentzVector(lepton);
  Double_t _mva_lep = colLep->getParameters().getFloatVal("ISOLepTagging");
  Double_t _lep_type = colLep->getParameters().getIntVal("ISOLepType");
  Double_t charge = lepton->getCharge();

  // ------------------------------------------------
  // -- read out the Photons information
  // ------------------------------------------------
  LCCollection *colPho = evt->getCollection(_colPhotons);
  if (!colPho) {
    cerr << "No Photons Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "Photons Collection Found" << ends;
  ReconstructedParticle *recoPhoton = dynamic_cast<ReconstructedParticle*>(colPho->getElementAt(0));
  TLorentzVector lortzPhoton = getLorentzVector(recoPhoton);
  std::vector<lcio::ReconstructedParticle*> photons = recoPhoton->getParticles();
  ReconstructedParticle *thePhoton = photons.at(0);
  TLorentzVector lortzPhotonOrig = getLorentzVector(thePhoton);
  // get photon cone energies
  Bool_t woFSR = kFALSE;
  Double_t coneEnergy0[3] = {0.,0.,0.};
  Double_t pFSR[4] = {0.,0.,0.,0.};
  Double_t pLargeCone[4]  = {0.,0.,0.,0.};
  Int_t nConePhoton = 0;
  getConeEnergy(thePhoton,colPFO,fCosSmallCone,woFSR,coneEnergy0,pFSR,fCosLargeCone,pLargeCone,nConePhoton);
  Double_t coneEN     = coneEnergy0[1];
  Double_t coneEC     = coneEnergy0[2];
  TLorentzVector lortzLargeCone = TLorentzVector(pLargeCone[0],pLargeCone[1],pLargeCone[2],pLargeCone[3]);
  Double_t cosThetaWithLargeCone = getCosTheta(lortzPhotonOrig,lortzLargeCone);
  Double_t energyRatioWithLargeCone = lortzPhotonOrig.E()/(lortzPhotonOrig.E()+lortzLargeCone.E());

  // ------------------------------------------------
  // -- read out the Jets information
  // ------------------------------------------------
  LCCollection *colJet = evt->getCollection(_colJets);
  if (!colJet) {
    cerr << "No Refined_4Jet Collection Found!" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "Refined_4Jet Collection Found" << ends;
  Int_t nJets = colJet->getNumberOfElements();
  if (nJets != 2) {
    cerr << "Number of Jets is not 2" << endl;
    throw marlin::SkipEventException(this);
  }
  hStatAnl->Fill(++selid);
  gCutName[(Int_t)selid] << "2 Jets Found" << ends;

  ReconstructedParticle *jets[2];
  // flavor tagging information
  PIDHandler pidh (colJet);
  Int_t algo = pidh.getAlgorithmID("lcfiplus");
  Double_t FLV[2][11];
  Int_t nPFOsCJ1 = 0, nPFOsCJ2 = 0, nPFOsCJ3 = 0,nPFOsCJ4 = 0;
  for (Int_t i=0;i<nJets;i++) {
    jets[i] = dynamic_cast<ReconstructedParticle*>(colJet->getElementAt(i));
    const ParticleID & jetID = pidh.getParticleID(jets[i], algo);
    FloatVec params = jetID.getParameters();
    FLV[i][0] = params[pidh.getParameterIndex(algo, "BTag")];
    FLV[i][1] = params[pidh.getParameterIndex(algo, "CTag")];
    FLV[i][2] = params[pidh.getParameterIndex(algo, "BCTag")];
    std::vector<lcio::ReconstructedParticle*> partVec = jets[i]->getParticles();
    for (std::vector<lcio::ReconstructedParticle*>::const_iterator iPart=partVec.begin();iPart!=partVec.end();++iPart) {
      TVector3 momPart = TVector3((*iPart)->getMomentum());
      Double_t pTPart = momPart.Pt();
      if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 0) nPFOsCJ1++;
      if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 1) nPFOsCJ2++;
     // if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 2) nPFOsCJ3++;
     // if ((*iPart)->getCharge() != 0 && pTPart > 0.5 && i == 3) nPFOsCJ4++;
    }
  }
  Int_t algo_y = pidh.getAlgorithmID("yth");
  const ParticleID & ythID = pidh.getParticleID(jets[0], algo_y);
  FloatVec params_y = ythID.getParameters();
  Double_t yMinus = params_y[pidh.getParameterIndex(algo_y, "y12")];
  Double_t yPlus  = params_y[pidh.getParameterIndex(algo_y, "y23")];
  Double_t yMinus4= params_y[pidh.getParameterIndex(algo_y, "y34")];
  Double_t yPlus4 = params_y[pidh.getParameterIndex(algo_y, "y45")];

  // ------------------------------------------------
  // -- get the useful physical quantities and save them to ntuple
  // ------------------------------------------------
  ReconstructedParticle *jet1 = jets[0];
  ReconstructedParticle *jet2 = jets[1];
 // ReconstructedParticle *jet3 = jets[2];
 // ReconstructedParticle *jet4 = jets[3];
  Int_t nPFOsJ1 = jet1->getParticles().size();
  Int_t nPFOsJ2 = jet2->getParticles().size();
 // Int_t nPFOsJ3 = jet3->getParticles().size();
 // Int_t nPFOsJ4 = jet4->getParticles().size();
 
  TLorentzVector lortzJ1 = getLorentzVector(jet1);
  TLorentzVector lortzJ2 = getLorentzVector(jet2);
//  TLorentzVector lortzJ3 = getLorentzVector(jet3);
//  TLorentzVector lortzJ4 = getLorentzVector(jet4);
  Double_t cosThetaJ1J2 = getCosTheta(lortzJ1,lortzJ2);
  Double_t cosThetaAJ1 = getCosTheta(lortzPhoton,lortzJ1);
  Double_t cosThetaAJ2 = getCosTheta(lortzPhoton,lortzJ2);
  Double_t cosThetaAJ = cosThetaAJ1 > cosThetaAJ2? cosThetaAJ1 : cosThetaAJ2;
  
  Double_t cosThetaLJ1 = getCosTheta(lortzLepton,lortzJ1);
  Double_t cosThetaLJ2 = getCosTheta(lortzLepton,lortzJ2);
  // get the flavor tagging information
  Double_t bProb_j1 = FLV[0][0];
  Double_t bProb_j2 = FLV[1][0];
 // Double_t bProb_j3 = FLV[2][0];
 // Double_t bProb_j4 = FLV[3][0];
  Double_t bcProb_j1 = FLV[0][2];
  Double_t bcProb_j2 = FLV[1][2];
  Double_t bcProb_j3 = FLV[2][2];
 // Double_t bcProb_j4 = FLV[3][2];
  // get the two most like b-jet
  
  
  Int_t nbmax1,nbmax2;
  Double_t bmax1=0.,bmax2=0.;
  for (Int_t i=0;i<2;i++) {
    if (FLV[i][0] >= bmax1) {
      nbmax1 = i;
      bmax1 = FLV[i][0];
    }
  }
  for (Int_t j=0;j<2;j++) {
    if (j != nbmax1 && FLV[j][0] >= bmax2) {
      nbmax2 = j;
      bmax2 = FLV[j][0];
    }
  }
  

  const Double_t fEcm =250.;
  TLorentzVector lortzEcm = getLorentzEcm(fEcm);
  TLorentzVector lortzMis = lortzEcm-lortzVis-lortzPhoton-lortzLepton;
  
  TLorentzVector lortzX = lortzJ1 + lortzJ2 + lortzLepton + lortzMis ;

  // ------------------------------------------------
  // -- reconstruct W  Yumi
  // ------------------------------------------------
  TLorentzVector lortzW1 = lortzJ1 + lortzJ2;
  TLorentzVector lortzW2 = lortzLepton + lortzMis;
  Double_t mw1rec = lortzW1.M();
  Double_t mw2rec = lortzW2.M();
//const Double_t fEcm =1000.;


  Double_t data_anl[100];
  data_anl[ 0]= lortzGam1MC.E();
  data_anl[ 1]= lortzISR1MC.E();
  data_anl[ 2]= lortzISR2MC.E();
  data_anl[ 3]= lortzXMC.M();    
  data_anl[ 4]= getRecoilMass(lortzEcm,lortzGam1MC);
  data_anl[ 5]= lortzGam1MC.CosTheta();
  data_anl[ 6]= lortzPhoton.E();
  data_anl[ 7]= getRecoilMass(lortzEcm,lortzPhoton);
  data_anl[ 8]= lortzPhoton.CosTheta();
  data_anl[ 9]= lortzPhotonOrig.E();
  data_anl[10]= getRecoilMass(lortzEcm,lortzPhotonOrig);
  data_anl[11]= lortzPhotonOrig.CosTheta();
  data_anl[12]= getMCSerial(thePhoton,colMCTL,colMC);
  data_anl[13]= lortzPhoton.M();
  data_anl[14]= coneEN;
  data_anl[15]= coneEC;
  data_anl[16]= cosThetaWithLargeCone;
  data_anl[17]= energyRatioWithLargeCone;
  data_anl[18]= lortzVis.E();
  data_anl[19] = lortzJ1.Px();
  data_anl[20] = lortzJ1.Py();
  data_anl[21] = lortzJ1.Pz();
  data_anl[22] = lortzJ1.E();
  data_anl[23] = lortzJ2.Px();
  data_anl[24] = lortzJ2.Py();
  data_anl[25] = lortzJ2.Pz();
  data_anl[26] = lortzJ2.E();
  data_anl[27] = lortzMis.Px();
  data_anl[28] = lortzMis.Py();
  data_anl[29] = lortzMis.Pz();
  data_anl[30] = lortzMis.E();
  data_anl[31] = lortzMis.Pt();  
  data_anl[32] = bProb_j1;
  data_anl[33] = bProb_j2;
  data_anl[34] = bmax1;
  data_anl[35] = bmax2;
  data_anl[36] = yMinus;
  data_anl[37] = yPlus;
  data_anl[38] = yMinus4;
  data_anl[39] = yPlus4;
  data_anl[40] = nJets;
  data_anl[41] = nParticles;
  data_anl[42] = nPFOsJ1;
  data_anl[43] = nPFOsJ2;
  data_anl[44] = nParticlesC;
  data_anl[45] = nParticlesPt500;
  data_anl[46] = nPFOsCJ1;
  data_anl[47] = nPFOsCJ2;
  data_anl[48] = cosThetaJ1J2;
  data_anl[49] = cosThetaAJ;  
  data_anl[50] = cosThetaAJ1;  
  data_anl[51] = cosThetaAJ2;  
  data_anl[52] = lortzX.M();
  data_anl[53] = lortzX.CosTheta();
  data_anl[54] = lortzPhoton.Px();
  data_anl[55] = lortzPhoton.Py();  
  data_anl[56] = lortzPhoton.Pz();
  data_anl[57]= lortzFinal.E();
  data_anl[58]= lortzFinal.Px();
  data_anl[59]= lortzFinal.Py();
  data_anl[60]= lortzFinal.Pz();
  data_anl[61]= lortz2f.M();
  data_anl[62]= mw1rec;
  data_anl[63]= mw2rec;
  data_anl[64]= m2WMC;
  data_anl[65]= nwqq;
  data_anl[66]= nwlv;
  data_anl[67]= mwplusmc;
  data_anl[68]= mwminusmc;
  data_anl[69] = lortzLepton.Px();
  data_anl[70] = lortzLepton.Py();
  data_anl[71] = lortzLepton.Pz();
  data_anl[72] = lortzLepton.E();
  data_anl[73] = lortzLepton.CosTheta();
 // data_anl[62]= lortzXNew.M();
 // data_anl[63]= cosJ1X;
 // data_anl[64]= cosJ2X;
  hAnl->Fill(data_anl);
  
  //-- note: this will not be printed if compiled w/o MARLINDEBUG=1 !

  streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
  		       << "   in run:  " << evt->getRunNumber() 
  		       << std::endl ;

  //  _nEvt ++ ;

  last->cd();
}



void AHlv2JAnalysisProcessor::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}


void AHlv2JAnalysisProcessor::end(){ 

  cerr << "AHlv2JAnalysisProcessor::end()  " << name() 
       << " processed " << _nEvt << " events in " << _nRun << " runs "
       << endl ;
  //  cerr << endl;
  cerr << "  =============" << endl;
  cerr << "   Cut Summary " << endl;
  cerr << "  =============" << endl;
  cerr << "   ll+4 Jet    " << endl;
  cerr << "  =============" << endl;
  cerr << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<20 && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStatAnl->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;
  
}

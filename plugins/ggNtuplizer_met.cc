//04.29.2020
//adding px,py,pt info for each met
//05.20.2020
//adding RAWMET info
//06.30.2020
//adding RAW PUPPI MET info

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t metFilters_;
float genMET_;
float genMETPhi_;
float pfMET_;
float pfMETPhi_;

//adding RAWMET for recoil
//float rawMET_;
float rawMET_pt_;
float rawMET_px_;
float rawMET_py_;
float rawMETPhi_;


//adding px, py, pt, sumEt info
float genMET_pt_;
float genMET_px_;
float genMET_py_;
float genMET_pz_;
float genMET_sumEt_;

float pfMET_pt_;
float pfMET_px_;
float pfMET_py_;
float pfMET_pz_;
float pfMET_sumEt_;

//metspuppi
float puppiMET_;
float puppiMETPhi_;

float puppiMET_pt_;
float puppiMET_px_;
float puppiMET_py_;
float puppiMET_pz_;
float puppiMET_sumEt_;


//RAW PUPPI MET
float rawPuppiMET_pt_;
float rawPuppiMET_px_;
float rawPuppiMET_py_;
float rawPuppiMETPhi_;

//
float pfMET_T1JERUp_;
float pfMET_T1JERDo_;
float pfMET_T1JESUp_;
float pfMET_T1JESDo_;
float pfMET_T1MESUp_;
float pfMET_T1MESDo_;
float pfMET_T1EESUp_;
float pfMET_T1EESDo_;
float pfMET_T1PESUp_;
float pfMET_T1PESDo_;
float pfMET_T1TESUp_;
float pfMET_T1TESDo_;
float pfMET_T1UESUp_;
float pfMET_T1UESDo_;
float pfMET_T1TxyPhi_;
float pfMET_T1TxyPt_;
float pfMETPhi_T1JESUp_;
float pfMETPhi_T1JESDo_;
float pfMETPhi_T1UESUp_;
float pfMETPhi_T1UESDo_;

void ggNtuplizer::branchesMET(TTree* tree) {

  if (doGenParticles_) {
    tree->Branch("genMET",      &genMET_);
    tree->Branch("genMETPhi",   &genMETPhi_);

    tree->Branch("genMET_pt",   &genMET_pt_);
    tree->Branch("genMET_px",   &genMET_px_);
    tree->Branch("genMET_py",   &genMET_py_);
    tree->Branch("genMET_pz",   &genMET_pz_);
    tree->Branch("genMET_sumEt",&genMET_sumEt_);
  }
  if (addFilterInfoMINIAOD_)
  tree->Branch("metFilters",       &metFilters_);
  tree->Branch("pfMET",            &pfMET_);
  tree->Branch("pfMETPhi",         &pfMETPhi_);

  tree->Branch("pfMET_pt",   &pfMET_pt_);
  tree->Branch("pfMET_px",   &pfMET_px_);
  tree->Branch("pfMET_py",   &pfMET_py_);
  tree->Branch("pfMET_pz",   &pfMET_pz_);
  tree->Branch("pfMET_sumEt",&pfMET_sumEt_);

  //rawmet
  //tree->Branch("rawMET",&rawMET_);
  tree->Branch("rawMET_pt",&rawMET_pt_);
  tree->Branch("rawMET_px",&rawMET_px_);
  tree->Branch("rawMET_py",&rawMET_py_);
  tree->Branch("rawMETPhi",&rawMETPhi_);


  //metspuppi//
  tree->Branch("puppiMET",        &puppiMET_);
  tree->Branch("puppiMETPhi",     &puppiMETPhi_);

  tree->Branch("puppiMET_pt",   &puppiMET_pt_);
  tree->Branch("puppiMET_px",   &puppiMET_px_);
  tree->Branch("puppiMET_py",   &puppiMET_py_);
  tree->Branch("puppiMET_pz",   &puppiMET_pz_);
  tree->Branch("puppiMET_sumEt",&puppiMET_sumEt_);

  //rawmet
  //tree->Branch("rawMET",&rawMET_);
  tree->Branch("rawPuppiMET_pt",&rawPuppiMET_pt_);
  tree->Branch("rawPuppiMET_px",&rawPuppiMET_px_);
  tree->Branch("rawPuppiMET_py",&rawPuppiMET_py_);
  tree->Branch("rawPuppiMETPhi",&rawPuppiMETPhi_);




  tree->Branch("pfMET_T1JERUp",    &pfMET_T1JERUp_);
  tree->Branch("pfMET_T1JERDo",    &pfMET_T1JERDo_);
  tree->Branch("pfMET_T1JESUp",    &pfMET_T1JESUp_);
  tree->Branch("pfMET_T1JESDo",    &pfMET_T1JESDo_);
  /*
  tree->Branch("pfMET_T1MESUp",    &pfMET_T1MESUp_);
  tree->Branch("pfMET_T1MESDo",    &pfMET_T1MESDo_);
  tree->Branch("pfMET_T1EESUp",    &pfMET_T1EESUp_);
  tree->Branch("pfMET_T1EESDo",    &pfMET_T1EESDo_);
  tree->Branch("pfMET_T1PESUp",    &pfMET_T1PESUp_);
  tree->Branch("pfMET_T1PESDo",    &pfMET_T1PESDo_);
  tree->Branch("pfMET_T1TESUp",    &pfMET_T1TESUp_);
  tree->Branch("pfMET_T1TESDo",    &pfMET_T1TESDo_);
  */
  tree->Branch("pfMET_T1UESUp",    &pfMET_T1UESUp_);
  tree->Branch("pfMET_T1UESDo",    &pfMET_T1UESDo_);
  tree->Branch("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp_);
  tree->Branch("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo_);
  tree->Branch("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp_);
  tree->Branch("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo_);

}



void ggNtuplizer::fillMET(const edm::Event& e, const edm::EventSetup& es) {

  metFilters_ = 0;

  if (addFilterInfoMINIAOD_) {
    string filterNamesToCheck[9] = {
      "Flag_HBHENoiseFilter",
      "Flag_HBHENoiseIsoFilter", 
      "Flag_globalSuperTightHalo2016Filter",
      "Flag_goodVertices",
      "Flag_eeBadScFilter",
      "Flag_EcalDeadCellTriggerPrimitiveFilter",
      "Flag_BadPFMuonFilter",
      "Flag_ecalBadCalibReducedMINIAODFilter",
      "Flag_BadChargedCandidateFilter"
    };

    edm::Handle<edm::TriggerResults> patFilterResultsHandle;
    e.getByToken(patTrgResultsLabel_, patFilterResultsHandle);
    edm::TriggerResults const& patFilterResults = *patFilterResultsHandle;
    
    auto&& filterNames = e.triggerNames(patFilterResults);

    // === the following lines allow us to find the filters stored in the event ! ===
    /*
    edm::TriggerNames const& triggerNames = e.triggerNames(patFilterResults);
    for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
	  triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
      int triggerId = triggerNames.triggerIndex(*triggerName);
      if ( triggerId >= 0 && triggerId < (int)triggerNames.size() ) {
	std::string triggerDecision = ( patFilterResultsHandle->accept(triggerId) ) ? "passed" : "failed";
	  
	std::cout << " triggerName = " << (*triggerName) << " " << triggerDecision << std::endl;
      }
    }    
    */
    
    for (unsigned iF = 0; iF < 9; ++iF) {
      unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
      if ( index == filterNames.size() ) {
	//std::cout<<filterNamesToCheck[iF] << " is missing, exiting"<<std::endl;

	edm::Handle<bool> passecalBadCalibFilterUpdate;
	e.getByToken(ecalBadCalibFilterUpdate_, passecalBadCalibFilterUpdate);
	if (passecalBadCalibFilterUpdate.isValid()) {
	  bool passecalBadCalibFilterUpdate_ = (*passecalBadCalibFilterUpdate);	
	  if (passecalBadCalibFilterUpdate_) metFilters_ += pow(2, iF+1);
	}

      } else {
	if ( !patFilterResults.accept(index) ) {
	  metFilters_ += pow(2, iF+1);
	}
      }
    }
  } 

//puppimets

  edm::Handle<edm::View<pat::MET> > metsPUPPIHandle;
  e.getByToken(metsPUPPI_, metsPUPPIHandle);

  puppiMET_  = -99;
  puppiMETPhi_ = -99;
  puppiMET_pt_ = -99;
  puppiMET_px_ = -99;
  puppiMET_py_ = -99;
  puppiMET_pz_ = -99;
  puppiMET_sumEt_ = -99;

  rawPuppiMET_pt_ = -99;
  rawPuppiMET_px_ = -99;
  rawPuppiMET_py_ = -99;
  rawPuppiMETPhi_ = -99;

  

  if (metsPUPPIHandle.isValid()){
    const pat::MET *puppiMET = 0;
    puppiMET  = &(metsPUPPIHandle->front());
    puppiMET_ = puppiMET->et();
    puppiMETPhi_ = puppiMET->phi();

    puppiMET_pt_ = puppiMET->pt();
    puppiMET_px_ = puppiMET->px();
    puppiMET_py_ = puppiMET->py();
    puppiMET_pz_ = puppiMET->pz();
    puppiMET_sumEt_ = puppiMET->sumEt();

    rawPuppiMET_pt_ = puppiMET->corPt(pat::MET::Raw);
    rawPuppiMET_px_ = puppiMET->corPx(pat::MET::Raw);
    rawPuppiMET_py_ = puppiMET->corPy(pat::MET::Raw);
    rawPuppiMETPhi_ = puppiMET->corPhi(pat::MET::Raw);
}

/*
//adding RAWMET info
edm::Handle<edm::View<pat::MET::Raw> > rawMETHandle;
e.getByToken(rawMETlabel_, rawMETHandle);

rawMET_ = -99;
rawMET_pt_ = -99;
rawMET_px_ = -99;
rawMET_py_ = -99;
rawMETPhi_ = -99;


if (rawMETHandle.isValid()) {
    const pat::MET::Raw *rawMET = 0;
    rawMET  = &(rawMETHandle->front());
    rawMET_ = rawMET->et();
    rawMET_pt_ = rawMET->pt();
    rawMET_px_ = rawMET->px();
    rawMET_py_ = rawMET->py();
    rawMETPhi_ = rawMET->phi();
}
*/

    



  
//
  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  e.getByToken(pfMETlabel_, pfMETHandle);

  genMET_    = -99;
  genMETPhi_ = -99;
  pfMET_     = -99;
  pfMETPhi_  = -99;

  pfMET_pt_ = -99;
  pfMET_px_ = -99;
  pfMET_py_ = -99;
  pfMET_pz_ = -99;
  pfMET_sumEt_ = -99;

  genMET_pt_ = -99;
  genMET_px_ = -99;
  genMET_py_ = -99;
  genMET_pz_ = -99;
  genMET_sumEt_ = -99;

    //rawMET_ = -99;
    rawMET_pt_ = -99;
    rawMET_px_ = -99;
    rawMET_py_ = -99;
    rawMETPhi_ = -99;


  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET     = &(pfMETHandle->front());
    pfMET_    = pfMET->et();
    pfMETPhi_ = pfMET->phi();

    pfMET_pt_ = pfMET->pt();
    pfMET_px_ = pfMET->px();
    pfMET_py_ = pfMET->py();
    pfMET_pz_ = pfMET->pz();
    pfMET_sumEt_ = pfMET->sumEt();
    
    // Type1MET uncertainties =======================================
    pfMET_T1JERUp_ = pfMET->shiftedPt(pat::MET::JetResUp);
    pfMET_T1JERDo_ = pfMET->shiftedPt(pat::MET::JetResDown);
    pfMET_T1JESUp_ = pfMET->shiftedPt(pat::MET::JetEnUp);
    pfMET_T1JESDo_ = pfMET->shiftedPt(pat::MET::JetEnDown);
    /*
      pfMET_T1MESUp_ = pfMET->shiftedPt(pat::MET::MuonEnUp);
      pfMET_T1MESDo_ = pfMET->shiftedPt(pat::MET::MuonEnDown);
      pfMET_T1EESUp_ = pfMET->shiftedPt(pat::MET::ElectronEnUp);
      pfMET_T1EESDo_ = pfMET->shiftedPt(pat::MET::ElectronEnDown);
      pfMET_T1PESUp_ = pfMET->shiftedPt(pat::MET::PhotonEnUp);
      pfMET_T1PESDo_ = pfMET->shiftedPt(pat::MET::PhotonEnDown);
      pfMET_T1TESUp_ = pfMET->shiftedPt(pat::MET::TauEnUp);
      pfMET_T1TESDo_ = pfMET->shiftedPt(pat::MET::TauEnDown);
    */
    pfMET_T1UESUp_ = pfMET->shiftedPt(pat::MET::UnclusteredEnUp);
    pfMET_T1UESDo_ = pfMET->shiftedPt(pat::MET::UnclusteredEnDown);
    
    pfMETPhi_T1JESUp_ = pfMET->shiftedPhi(pat::MET::JetEnUp);
    pfMETPhi_T1JESDo_ = pfMET->shiftedPhi(pat::MET::JetEnDown);
    pfMETPhi_T1UESUp_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnUp);
    pfMETPhi_T1UESDo_ = pfMET->shiftedPhi(pat::MET::UnclusteredEnDown);



    //store rawmet from uncor pf
    //rawMET_ = pfMET->corEt(pat::MET::Raw);
    rawMET_pt_ = pfMET->corPt(pat::MET::Raw);
    rawMET_px_ = pfMET->corPx(pat::MET::Raw);
    rawMET_py_ = pfMET->corPy(pat::MET::Raw);
    rawMETPhi_ = pfMET->corPhi(pat::MET::Raw);
    
    if (!e.isRealData()) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();


      genMET_pt_ = pfMET->genMET()->pt();
      genMET_px_ = pfMET->genMET()->px();
      genMET_py_ = pfMET->genMET()->py();
      genMET_pz_ = pfMET->genMET()->pz();
      genMET_sumEt_ = pfMET->genMET()->sumEt();
    }

  } 

}

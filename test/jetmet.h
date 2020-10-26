//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jul 24 15:02:01 2020 by ROOT version 6.18/04
// from TTree EventTree/Event data (tag V10_02_10_04)
// found on file: ggtree_mc2017_GJET40to100.root
//////////////////////////////////////////////////////////

#ifndef jetmet_h
#define jetmet_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "TString.h"
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class jetmet {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Long64_t        event;
   Int_t           lumis;
   Bool_t          isData;
   Int_t           nVtx;
   Int_t           nGoodVtx;
   Bool_t          isPVGood;
   Float_t         vtx;
   Float_t         vty;
   Float_t         vtz;
   Float_t         rho;
   Float_t         rhoCentral;
   Double_t        L1ECALPrefire;
   Double_t        L1ECALPrefireUp;
   Double_t        L1ECALPrefireDown;
   ULong64_t       HLTEleMuX;
   ULong64_t       HLTPho;
   ULong64_t       HLTPhoRejectedByPS;
   ULong64_t       HLTJet;
   ULong64_t       HLTEleMuXIsPrescaled;
   ULong64_t       HLTPhoIsPrescaled;
   ULong64_t       HLTJetIsPrescaled;
   vector<float>   *pdf;
   Float_t         pthat;
   Float_t         processID;
   Float_t         genWeight;
   Float_t         genHT;
   Float_t         genPho1;
   Float_t         genPho2;
   TString         *EventTag;
   Int_t           nPUInfo;
   vector<int>     *nPU;
   vector<int>     *puBX;
   vector<float>   *puTrue;
   Int_t           nLHE;
   vector<int>     *lhePID;
   vector<float>   *lhePx;
   vector<float>   *lhePy;
   vector<float>   *lhePz;
   vector<float>   *lheE;
   Int_t           nMC;
   vector<int>     *mcPID;
   vector<float>   *mcVtx;
   vector<float>   *mcVty;
   vector<float>   *mcVtz;
   vector<float>   *mcPt;
   vector<float>   *mcMass;
   vector<float>   *mcEta;
   vector<float>   *mcPhi;
   vector<float>   *mcE;
   vector<float>   *mcEt;
   vector<int>     *mcGMomPID;
   vector<int>     *mcMomPID;
   vector<float>   *mcMomPt;
   vector<float>   *mcMomMass;
   vector<float>   *mcMomEta;
   vector<float>   *mcMomPhi;
   vector<unsigned short> *mcStatusFlag;
   vector<int>     *mcParentage;
   vector<int>     *mcStatus;
   vector<float>   *mcCalIsoDR03;
   vector<float>   *mcTrkIsoDR03;
   vector<float>   *mcCalIsoDR04;
   vector<float>   *mcTrkIsoDR04;
   Float_t         genMET;
   Float_t         genMETPhi;
   Float_t         genMET_pt;
   Float_t         genMET_px;
   Float_t         genMET_py;
   Float_t         genMET_pz;
   Float_t         genMET_sumEt;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Float_t         pfMET_pt;
   Float_t         pfMET_px;
   Float_t         pfMET_py;
   Float_t         pfMET_pz;
   Float_t         pfMET_sumEt;
   Float_t         rawMET_pt;
   Float_t         rawMET_px;
   Float_t         rawMET_py;
   Float_t         rawMETPhi;
   Float_t         puppiMET;
   Float_t         puppiMETPhi;
   Float_t         puppiMET_pt;
   Float_t         puppiMET_px;
   Float_t         puppiMET_py;
   Float_t         puppiMET_pz;
   Float_t         puppiMET_sumEt;
   Float_t         rawPuppiMET_pt;
   Float_t         rawPuppiMET_px;
   Float_t         rawPuppiMET_py;
   Float_t         rawPuppiMETPhi;
   Float_t         pfMET_T1JERUp;
   Float_t         pfMET_T1JERDo;
   Float_t         pfMET_T1JESUp;
   Float_t         pfMET_T1JESDo;
   Float_t         pfMET_T1UESUp;
   Float_t         pfMET_T1UESDo;
   Float_t         pfMETPhi_T1JESUp;
   Float_t         pfMETPhi_T1JESDo;
   Float_t         pfMETPhi_T1UESUp;
   Float_t         pfMETPhi_T1UESDo;
   Int_t           nPho;
   vector<float>   *phoE;
   vector<float>   *phoSigmaE;
   vector<float>   *phoEt;
   vector<float>   *phoEta;
   vector<float>   *phoPhi;
   vector<float>   *phoCalibE;
   vector<float>   *phoCalibEt;
   vector<float>   *phoSCE;
   vector<float>   *phoSCRawE;
   vector<float>   *phoESEnP1;
   vector<float>   *phoESEnP2;
   vector<float>   *phoSCEta;
   vector<float>   *phoSCPhi;
   vector<float>   *phoSCEtaWidth;
   vector<float>   *phoSCPhiWidth;
   vector<float>   *phoSCBrem;
   vector<int>     *phohasPixelSeed;
   vector<int>     *phoEleVeto;
   vector<float>   *phoR9;
   vector<float>   *phoHoverE;
   vector<float>   *phoESEffSigmaRR;
   vector<float>   *phoSigmaIEtaIEtaFull5x5;
   vector<float>   *phoSigmaIEtaIPhiFull5x5;
   vector<float>   *phoSigmaIPhiIPhiFull5x5;
   vector<float>   *phoE2x2Full5x5;
   vector<float>   *phoE5x5Full5x5;
   vector<float>   *phoR9Full5x5;
   vector<float>   *phoPFChIso;
   vector<float>   *phoPFChPVIso;
   vector<float>   *phoPFPhoIso;
   vector<float>   *phoPFNeuIso;
   vector<float>   *phoPFChWorstIso;
   vector<float>   *phoPFChWorstVetoIso;
   vector<float>   *phoEcalPFClusterIso;
   vector<float>   *phoHcalPFClusterIso;
   vector<float>   *phoIDMVA;
   vector<ULong64_t> *phoFiredSingleTrgs;
   vector<ULong64_t> *phoFiredDoubleTrgs;
   vector<ULong64_t> *phoFiredTripleTrgs;
   vector<ULong64_t> *phoFiredL1Trgs;
   vector<float>   *phoSeedTime;
   vector<float>   *phoSeedEnergy;
   vector<float>   *phoMIPTotEnergy;
   vector<float>   *phoMIPChi2;
   vector<float>   *phoMIPSlope;
   vector<float>   *phoMIPIntercept;
   vector<float>   *phoMIPNhitCone;
   vector<float>   *phoMIPIsHalo;
   vector<unsigned short> *phoxtalBits;
   vector<unsigned short> *phoIDbit;
   vector<float>   *phoScale_stat_up;
   vector<float>   *phoScale_stat_dn;
   vector<float>   *phoScale_syst_up;
   vector<float>   *phoScale_syst_dn;
   vector<float>   *phoScale_gain_up;
   vector<float>   *phoScale_gain_dn;
   vector<float>   *phoResol_rho_up;
   vector<float>   *phoResol_rho_dn;
   vector<float>   *phoResol_phi_up;
   vector<float>   *phoResol_phi_dn;
   Int_t           nEle;
   vector<int>     *eleCharge;
   vector<int>     *eleChargeConsistent;
   vector<float>   *eleEn;
   vector<float>   *eleSCEn;
   vector<float>   *eleEcalEn;
   vector<float>   *eleESEnP1;
   vector<float>   *eleESEnP2;
   vector<float>   *eleD0;
   vector<float>   *eleDz;
   vector<float>   *eleSIP;
   vector<float>   *elePt;
   vector<float>   *elePtError;
   vector<float>   *eleEta;
   vector<float>   *elePhi;
   vector<float>   *eleR9;
   vector<float>   *eleCalibPt;
   vector<float>   *eleCalibEn;
   vector<float>   *eleSCEta;
   vector<float>   *eleSCPhi;
   vector<float>   *eleSCRawEn;
   vector<float>   *eleSCEtaWidth;
   vector<float>   *eleSCPhiWidth;
   vector<float>   *eleHoverE;
   vector<float>   *eleEoverP;
   vector<float>   *eleEoverPout;
   vector<float>   *eleEoverPInv;
   vector<float>   *eleBrem;
   vector<float>   *eledEtaAtVtx;
   vector<float>   *eledPhiAtVtx;
   vector<float>   *eleSigmaIEtaIEtaFull5x5;
   vector<float>   *eleSigmaIPhiIPhiFull5x5;
   vector<int>     *eleConvVeto;
   vector<int>     *eleMissHits;
   vector<float>   *eleESEffSigmaRR;
   vector<float>   *elePFChIso;
   vector<float>   *elePFPhoIso;
   vector<float>   *elePFNeuIso;
   vector<float>   *elePFPUIso;
   vector<float>   *elePFClusEcalIso;
   vector<float>   *elePFClusHcalIso;
   vector<float>   *eleIDMVAIso;
   vector<float>   *eleIDMVANoIso;
   vector<float>   *eleR9Full5x5;
   vector<int>     *eleEcalDrivenSeed;
   vector<float>   *eleTrkdxy;
   vector<float>   *eleKFHits;
   vector<float>   *eleKFChi2;
   vector<float>   *eleGSFChi2;
   vector<vector<float> > *eleGSFPt;
   vector<vector<float> > *eleGSFEta;
   vector<vector<float> > *eleGSFPhi;
   vector<vector<float> > *eleGSFCharge;
   vector<vector<int> > *eleGSFHits;
   vector<vector<int> > *eleGSFMissHits;
   vector<vector<int> > *eleGSFNHitsMax;
   vector<vector<float> > *eleGSFVtxProb;
   vector<vector<float> > *eleGSFlxyPV;
   vector<vector<float> > *eleGSFlxyBS;
   vector<ULong64_t> *eleFiredSingleTrgs;
   vector<ULong64_t> *eleFiredDoubleTrgs;
   vector<ULong64_t> *eleFiredL1Trgs;
   vector<unsigned short> *eleIDbit;
   vector<float>   *eleScale_stat_up;
   vector<float>   *eleScale_stat_dn;
   vector<float>   *eleScale_syst_up;
   vector<float>   *eleScale_syst_dn;
   vector<float>   *eleScale_gain_up;
   vector<float>   *eleScale_gain_dn;
   vector<float>   *eleResol_rho_up;
   vector<float>   *eleResol_rho_dn;
   vector<float>   *eleResol_phi_up;
   vector<float>   *eleResol_phi_dn;
   Int_t           nMu;
   vector<float>   *muPt;
   vector<float>   *muEn;
   vector<float>   *muEta;
   vector<float>   *muPhi;
   vector<int>     *muCharge;
   vector<int>     *muType;
   vector<int>     *muIDbit;
   vector<float>   *muD0;
   vector<float>   *muDz;
   vector<float>   *muSIP;
   vector<float>   *muChi2NDF;
   vector<float>   *muInnerD0;
   vector<float>   *muInnerDz;
   vector<int>     *muTrkLayers;
   vector<int>     *muPixelLayers;
   vector<int>     *muPixelHits;
   vector<int>     *muMuonHits;
   vector<int>     *muStations;
   vector<int>     *muMatches;
   vector<int>     *muTrkQuality;
   vector<float>   *muIsoTrk;
   vector<float>   *muPFChIso;
   vector<float>   *muPFPhoIso;
   vector<float>   *muPFNeuIso;
   vector<float>   *muPFPUIso;
   vector<float>   *muPFChIso03;
   vector<float>   *muPFPhoIso03;
   vector<float>   *muPFNeuIso03;
   vector<float>   *muPFPUIso03;
   vector<ULong64_t> *muFiredTrgs;
   vector<ULong64_t> *muFiredL1Trgs;
   vector<float>   *muInnervalidFraction;
   vector<float>   *musegmentCompatibility;
   vector<float>   *muchi2LocalPosition;
   vector<float>   *mutrkKink;
   vector<float>   *muBestTrkPtError;
   vector<float>   *muBestTrkPt;
   vector<int>     *muBestTrkType;
   Int_t           npfPho;
   vector<float>   *pfPhoEt;
   vector<float>   *pfPhoEta;
   vector<float>   *pfPhoPhi;
   vector<float>   *pfPhoEtErr;
   vector<float>   *pfPhoChIso;
   vector<float>   *pfPhoPhoIso;
   vector<float>   *pfPhoNeuIso;
   vector<int>     *pfPhoEleIndex;
   vector<int>     *pfPhoPhoIndex;
   Int_t           nJet;
   vector<float>   *jetPt;
   vector<float>   *jetEn;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetRawPt;
   vector<float>   *jetRawEn;
   vector<float>   *jetMt;
   vector<float>   *jetArea;
   vector<float>   *jetLeadTrackPt;
   vector<float>   *jetLeadTrackEta;
   vector<float>   *jetLeadTrackPhi;
   vector<int>     *jetLepTrackPID;
   vector<float>   *jetLepTrackPt;
   vector<float>   *jetLepTrackEta;
   vector<float>   *jetLepTrackPhi;
   vector<float>   *jetCSV2BJetTags;
   vector<float>   *jetDeepCSVTags_b;
   vector<float>   *jetDeepCSVTags_bb;
   vector<float>   *jetDeepCSVTags_c;
   vector<float>   *jetDeepCSVTags_udsg;
   vector<int>     *jetPartonID;
   vector<int>     *jetHadFlvr;
   vector<float>   *jetGenJetEn;
   vector<float>   *jetGenJetPt;
   vector<float>   *jetGenJetEta;
   vector<float>   *jetGenJetPhi;
   vector<int>     *jetGenPartonID;
   vector<float>   *jetGenEn;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenPhi;
   vector<int>     *jetGenPartonMomID;
   vector<float>   *jetP4Smear;
   vector<float>   *jetP4SmearUp;
   vector<float>   *jetP4SmearDo;
   vector<bool>    *jetPFLooseId;
   vector<int>     *jetID;
   vector<float>   *jetPUID;
   vector<int>     *jetPUFullID;
   vector<float>   *jetJECUnc;
   vector<ULong64_t> *jetFiredTrgs;
   vector<float>   *jetCHF;
   vector<float>   *jetNHF;
   vector<float>   *jetCEF;
   vector<float>   *jetNEF;
   vector<int>     *jetNCH;
   vector<int>     *jetNNP;
   vector<float>   *jetMUF;
   vector<float>   *jetVtxPt;
   vector<float>   *jetVtxMass;
   vector<float>   *jetVtxNtrks;
   vector<float>   *jetVtx3DVal;
   vector<float>   *jetVtx3DSig;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumis;   //!
   TBranch        *b_isData;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_nGoodVtx;   //!
   TBranch        *b_isPVGood;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vty;   //!
   TBranch        *b_vtz;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_rhoCentral;   //!
   TBranch        *b_L1ECALPrefire;   //!
   TBranch        *b_L1ECALPrefireUp;   //!
   TBranch        *b_L1ECALPrefireDown;   //!
   TBranch        *b_HLTEleMuX;   //!
   TBranch        *b_HLTPho;   //!
   TBranch        *b_HLTPhoRejectedByPS;   //!
   TBranch        *b_HLTJet;   //!
   TBranch        *b_HLTEleMuXIsPrescaled;   //!
   TBranch        *b_HLTPhoIsPrescaled;   //!
   TBranch        *b_HLTJetIsPrescaled;   //!
   TBranch        *b_pdf;   //!
   TBranch        *b_pthat;   //!
   TBranch        *b_processID;   //!
   TBranch        *b_genWeight;   //!
   TBranch        *b_genHT;   //!
   TBranch        *b_genPho1;   //!
   TBranch        *b_genPho2;   //!
   TBranch        *b_EventTag;   //!
   TBranch        *b_nPUInfo;   //!
   TBranch        *b_nPU;   //!
   TBranch        *b_puBX;   //!
   TBranch        *b_puTrue;   //!
   TBranch        *b_nLHE;   //!
   TBranch        *b_lhePID;   //!
   TBranch        *b_lhePx;   //!
   TBranch        *b_lhePy;   //!
   TBranch        *b_lhePz;   //!
   TBranch        *b_lheE;   //!
   TBranch        *b_nMC;   //!
   TBranch        *b_mcPID;   //!
   TBranch        *b_mcVtx;   //!
   TBranch        *b_mcVty;   //!
   TBranch        *b_mcVtz;   //!
   TBranch        *b_mcPt;   //!
   TBranch        *b_mcMass;   //!
   TBranch        *b_mcEta;   //!
   TBranch        *b_mcPhi;   //!
   TBranch        *b_mcE;   //!
   TBranch        *b_mcEt;   //!
   TBranch        *b_mcGMomPID;   //!
   TBranch        *b_mcMomPID;   //!
   TBranch        *b_mcMomPt;   //!
   TBranch        *b_mcMomMass;   //!
   TBranch        *b_mcMomEta;   //!
   TBranch        *b_mcMomPhi;   //!
   TBranch        *b_mcStatusFlag;   //!
   TBranch        *b_mcParentage;   //!
   TBranch        *b_mcStatus;   //!
   TBranch        *b_mcCalIsoDR03;   //!
   TBranch        *b_mcTrkIsoDR03;   //!
   TBranch        *b_mcCalIsoDR04;   //!
   TBranch        *b_mcTrkIsoDR04;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMETPhi;   //!
   TBranch        *b_genMET_pt;   //!
   TBranch        *b_genMET_px;   //!
   TBranch        *b_genMET_py;   //!
   TBranch        *b_genMET_pz;   //!
   TBranch        *b_genMET_sumEt;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_pfMET_pt;   //!
   TBranch        *b_pfMET_px;   //!
   TBranch        *b_pfMET_py;   //!
   TBranch        *b_pfMET_pz;   //!
   TBranch        *b_pfMET_sumEt;   //!
   TBranch        *b_rawMET_pt;   //!
   TBranch        *b_rawMET_px;   //!
   TBranch        *b_rawMET_py;   //!
   TBranch        *b_rawMETPhi;   //!
   TBranch        *b_puppiMET;   //!
   TBranch        *b_puppiMETPhi;   //!
   TBranch        *b_puppiMET_pt;   //!
   TBranch        *b_puppiMET_px;   //!
   TBranch        *b_puppiMET_py;   //!
   TBranch        *b_puppiMET_pz;   //!
   TBranch        *b_puppiMET_sumEt;   //!
   TBranch        *b_rawPuppiMET_pt;   //!
   TBranch        *b_rawPuppiMET_px;   //!
   TBranch        *b_rawPuppiMET_py;   //!
   TBranch        *b_rawPuppiMETPhi;   //!
   TBranch        *b_pfMET_T1JERUp;   //!
   TBranch        *b_pfMET_T1JERDo;   //!
   TBranch        *b_pfMET_T1JESUp;   //!
   TBranch        *b_pfMET_T1JESDo;   //!
   TBranch        *b_pfMET_T1UESUp;   //!
   TBranch        *b_pfMET_T1UESDo;   //!
   TBranch        *b_pfMETPhi_T1JESUp;   //!
   TBranch        *b_pfMETPhi_T1JESDo;   //!
   TBranch        *b_pfMETPhi_T1UESUp;   //!
   TBranch        *b_pfMETPhi_T1UESDo;   //!
   TBranch        *b_nPho;   //!
   TBranch        *b_phoE;   //!
   TBranch        *b_phoSigmaE;   //!
   TBranch        *b_phoEt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoCalibE;   //!
   TBranch        *b_phoCalibEt;   //!
   TBranch        *b_phoSCE;   //!
   TBranch        *b_phoSCRawE;   //!
   TBranch        *b_phoESEnP1;   //!
   TBranch        *b_phoESEnP2;   //!
   TBranch        *b_phoSCEta;   //!
   TBranch        *b_phoSCPhi;   //!
   TBranch        *b_phoSCEtaWidth;   //!
   TBranch        *b_phoSCPhiWidth;   //!
   TBranch        *b_phoSCBrem;   //!
   TBranch        *b_phohasPixelSeed;   //!
   TBranch        *b_phoEleVeto;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoHoverE;   //!
   TBranch        *b_phoESEffSigmaRR;   //!
   TBranch        *b_phoSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_phoSigmaIEtaIPhiFull5x5;   //!
   TBranch        *b_phoSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_phoE2x2Full5x5;   //!
   TBranch        *b_phoE5x5Full5x5;   //!
   TBranch        *b_phoR9Full5x5;   //!
   TBranch        *b_phoPFChIso;   //!
   TBranch        *b_phoPFChPVIso;   //!
   TBranch        *b_phoPFPhoIso;   //!
   TBranch        *b_phoPFNeuIso;   //!
   TBranch        *b_phoPFChWorstIso;   //!
   TBranch        *b_phoPFChWorstVetoIso;   //!
   TBranch        *b_phoEcalPFClusterIso;   //!
   TBranch        *b_phoHcalPFClusterIso;   //!
   TBranch        *b_phoIDMVA;   //!
   TBranch        *b_phoFiredSingleTrgs;   //!
   TBranch        *b_phoFiredDoubleTrgs;   //!
   TBranch        *b_phoFiredTripleTrgs;   //!
   TBranch        *b_phoFiredL1Trgs;   //!
   TBranch        *b_phoSeedTime;   //!
   TBranch        *b_phoSeedEnergy;   //!
   TBranch        *b_phoMIPTotEnergy;   //!
   TBranch        *b_phoMIPChi2;   //!
   TBranch        *b_phoMIPSlope;   //!
   TBranch        *b_phoMIPIntercept;   //!
   TBranch        *b_phoMIPNhitCone;   //!
   TBranch        *b_phoMIPIsHalo;   //!
   TBranch        *b_phoxtalBits;   //!
   TBranch        *b_phoIDbit;   //!
   TBranch        *b_phoScale_stat_up;   //!
   TBranch        *b_phoScale_stat_dn;   //!
   TBranch        *b_phoScale_syst_up;   //!
   TBranch        *b_phoScale_syst_dn;   //!
   TBranch        *b_phoScale_gain_up;   //!
   TBranch        *b_phoScale_gain_dn;   //!
   TBranch        *b_phoResol_rho_up;   //!
   TBranch        *b_phoResol_rho_dn;   //!
   TBranch        *b_phoResol_phi_up;   //!
   TBranch        *b_phoResol_phi_dn;   //!
   TBranch        *b_nEle;   //!
   TBranch        *b_eleCharge;   //!
   TBranch        *b_eleChargeConsistent;   //!
   TBranch        *b_eleEn;   //!
   TBranch        *b_eleSCEn;   //!
   TBranch        *b_eleEcalEn;   //!
   TBranch        *b_eleESEnP1;   //!
   TBranch        *b_eleESEnP2;   //!
   TBranch        *b_eleD0;   //!
   TBranch        *b_eleDz;   //!
   TBranch        *b_eleSIP;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_elePtError;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleR9;   //!
   TBranch        *b_eleCalibPt;   //!
   TBranch        *b_eleCalibEn;   //!
   TBranch        *b_eleSCEta;   //!
   TBranch        *b_eleSCPhi;   //!
   TBranch        *b_eleSCRawEn;   //!
   TBranch        *b_eleSCEtaWidth;   //!
   TBranch        *b_eleSCPhiWidth;   //!
   TBranch        *b_eleHoverE;   //!
   TBranch        *b_eleEoverP;   //!
   TBranch        *b_eleEoverPout;   //!
   TBranch        *b_eleEoverPInv;   //!
   TBranch        *b_eleBrem;   //!
   TBranch        *b_eledEtaAtVtx;   //!
   TBranch        *b_eledPhiAtVtx;   //!
   TBranch        *b_eleSigmaIEtaIEtaFull5x5;   //!
   TBranch        *b_eleSigmaIPhiIPhiFull5x5;   //!
   TBranch        *b_eleConvVeto;   //!
   TBranch        *b_eleMissHits;   //!
   TBranch        *b_eleESEffSigmaRR;   //!
   TBranch        *b_elePFChIso;   //!
   TBranch        *b_elePFPhoIso;   //!
   TBranch        *b_elePFNeuIso;   //!
   TBranch        *b_elePFPUIso;   //!
   TBranch        *b_elePFClusEcalIso;   //!
   TBranch        *b_elePFClusHcalIso;   //!
   TBranch        *b_eleIDMVAIso;   //!
   TBranch        *b_eleIDMVANoIso;   //!
   TBranch        *b_eleR9Full5x5;   //!
   TBranch        *b_eleEcalDrivenSeed;   //!
   TBranch        *b_eleTrkdxy;   //!
   TBranch        *b_eleKFHits;   //!
   TBranch        *b_eleKFChi2;   //!
   TBranch        *b_eleGSFChi2;   //!
   TBranch        *b_eleGSFPt;   //!
   TBranch        *b_eleGSFEta;   //!
   TBranch        *b_eleGSFPhi;   //!
   TBranch        *b_eleGSFCharge;   //!
   TBranch        *b_eleGSFHits;   //!
   TBranch        *b_eleGSFMissHits;   //!
   TBranch        *b_eleGSFNHitsMax;   //!
   TBranch        *b_eleGSFVtxProb;   //!
   TBranch        *b_eleGSFlxyPV;   //!
   TBranch        *b_eleGSFlxyBS;   //!
   TBranch        *b_eleFiredSingleTrgs;   //!
   TBranch        *b_eleFiredDoubleTrgs;   //!
   TBranch        *b_eleFiredL1Trgs;   //!
   TBranch        *b_eleIDbit;   //!
   TBranch        *b_eleScale_stat_up;   //!
   TBranch        *b_eleScale_stat_dn;   //!
   TBranch        *b_eleScale_syst_up;   //!
   TBranch        *b_eleScale_syst_dn;   //!
   TBranch        *b_eleScale_gain_up;   //!
   TBranch        *b_eleScale_gain_dn;   //!
   TBranch        *b_eleResol_rho_up;   //!
   TBranch        *b_eleResol_rho_dn;   //!
   TBranch        *b_eleResol_phi_up;   //!
   TBranch        *b_eleResol_phi_dn;   //!
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEn;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIDbit;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muSIP;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muMatches;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
   TBranch        *b_muPFChIso03;   //!
   TBranch        *b_muPFPhoIso03;   //!
   TBranch        *b_muPFNeuIso03;   //!
   TBranch        *b_muPFPUIso03;   //!
   TBranch        *b_muFiredTrgs;   //!
   TBranch        *b_muFiredL1Trgs;   //!
   TBranch        *b_muInnervalidFraction;   //!
   TBranch        *b_musegmentCompatibility;   //!
   TBranch        *b_muchi2LocalPosition;   //!
   TBranch        *b_mutrkKink;   //!
   TBranch        *b_muBestTrkPtError;   //!
   TBranch        *b_muBestTrkPt;   //!
   TBranch        *b_muBestTrkType;   //!
   TBranch        *b_npfPho;   //!
   TBranch        *b_pfPhoEt;   //!
   TBranch        *b_pfPhoEta;   //!
   TBranch        *b_pfPhoPhi;   //!
   TBranch        *b_pfPhoEtErr;   //!
   TBranch        *b_pfPhoChIso;   //!
   TBranch        *b_pfPhoPhoIso;   //!
   TBranch        *b_pfPhoNeuIso;   //!
   TBranch        *b_pfPhoEleIndex;   //!
   TBranch        *b_pfPhoPhoIndex;   //!
   TBranch        *b_nJet;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEn;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetRawPt;   //!
   TBranch        *b_jetRawEn;   //!
   TBranch        *b_jetMt;   //!
   TBranch        *b_jetArea;   //!
   TBranch        *b_jetLeadTrackPt;   //!
   TBranch        *b_jetLeadTrackEta;   //!
   TBranch        *b_jetLeadTrackPhi;   //!
   TBranch        *b_jetLepTrackPID;   //!
   TBranch        *b_jetLepTrackPt;   //!
   TBranch        *b_jetLepTrackEta;   //!
   TBranch        *b_jetLepTrackPhi;   //!
   TBranch        *b_jetCSV2BJetTags;   //!
   TBranch        *b_jetDeepCSVTags_b;   //!
   TBranch        *b_jetDeepCSVTags_bb;   //!
   TBranch        *b_jetDeepCSVTags_c;   //!
   TBranch        *b_jetDeepCSVTags_udsg;   //!
   TBranch        *b_jetPartonID;   //!
   TBranch        *b_jetHadFlvr;   //!
   TBranch        *b_jetGenJetEn;   //!
   TBranch        *b_jetGenJetPt;   //!
   TBranch        *b_jetGenJetEta;   //!
   TBranch        *b_jetGenJetPhi;   //!
   TBranch        *b_jetGenPartonID;   //!
   TBranch        *b_jetGenEn;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenPhi;   //!
   TBranch        *b_jetGenPartonMomID;   //!
   TBranch        *b_jetP4Smear;   //!
   TBranch        *b_jetP4SmearUp;   //!
   TBranch        *b_jetP4SmearDo;   //!
   TBranch        *b_jetPFLooseId;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_jetPUID;   //!
   TBranch        *b_jetPUFullID;   //!
   TBranch        *b_jetJECUnc;   //!
   TBranch        *b_jetFiredTrgs;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetCEF;   //!
   TBranch        *b_jetNEF;   //!
   TBranch        *b_jetNCH;   //!
   TBranch        *b_jetNNP;   //!
   TBranch        *b_jetMUF;   //!
   TBranch        *b_jetVtxPt;   //!
   TBranch        *b_jetVtxMass;   //!
   TBranch        *b_jetVtxNtrks;   //!
   TBranch        *b_jetVtx3DVal;   //!
   TBranch        *b_jetVtx3DSig;   //!

   jetmet(TTree *tree=0);
   virtual ~jetmet();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef jetmet_cxx
jetmet::jetmet(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ggtree_mc2017_GJET40to100.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ggtree_mc2017_GJET40to100.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ggtree_mc2017_GJET40to100.root:/ggNtuplizer");
      dir->GetObject("EventTree",tree);

   }
   Init(tree);
}

jetmet::~jetmet()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t jetmet::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t jetmet::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void jetmet::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   pdf = 0;
   EventTag = 0;
   nPU = 0;
   puBX = 0;
   puTrue = 0;
   lhePID = 0;
   lhePx = 0;
   lhePy = 0;
   lhePz = 0;
   lheE = 0;
   mcPID = 0;
   mcVtx = 0;
   mcVty = 0;
   mcVtz = 0;
   mcPt = 0;
   mcMass = 0;
   mcEta = 0;
   mcPhi = 0;
   mcE = 0;
   mcEt = 0;
   mcGMomPID = 0;
   mcMomPID = 0;
   mcMomPt = 0;
   mcMomMass = 0;
   mcMomEta = 0;
   mcMomPhi = 0;
   mcStatusFlag = 0;
   mcParentage = 0;
   mcStatus = 0;
   mcCalIsoDR03 = 0;
   mcTrkIsoDR03 = 0;
   mcCalIsoDR04 = 0;
   mcTrkIsoDR04 = 0;
   phoE = 0;
   phoSigmaE = 0;
   phoEt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoCalibE = 0;
   phoCalibEt = 0;
   phoSCE = 0;
   phoSCRawE = 0;
   phoESEnP1 = 0;
   phoESEnP2 = 0;
   phoSCEta = 0;
   phoSCPhi = 0;
   phoSCEtaWidth = 0;
   phoSCPhiWidth = 0;
   phoSCBrem = 0;
   phohasPixelSeed = 0;
   phoEleVeto = 0;
   phoR9 = 0;
   phoHoverE = 0;
   phoESEffSigmaRR = 0;
   phoSigmaIEtaIEtaFull5x5 = 0;
   phoSigmaIEtaIPhiFull5x5 = 0;
   phoSigmaIPhiIPhiFull5x5 = 0;
   phoE2x2Full5x5 = 0;
   phoE5x5Full5x5 = 0;
   phoR9Full5x5 = 0;
   phoPFChIso = 0;
   phoPFChPVIso = 0;
   phoPFPhoIso = 0;
   phoPFNeuIso = 0;
   phoPFChWorstIso = 0;
   phoPFChWorstVetoIso = 0;
   phoEcalPFClusterIso = 0;
   phoHcalPFClusterIso = 0;
   phoIDMVA = 0;
   phoFiredSingleTrgs = 0;
   phoFiredDoubleTrgs = 0;
   phoFiredTripleTrgs = 0;
   phoFiredL1Trgs = 0;
   phoSeedTime = 0;
   phoSeedEnergy = 0;
   phoMIPTotEnergy = 0;
   phoMIPChi2 = 0;
   phoMIPSlope = 0;
   phoMIPIntercept = 0;
   phoMIPNhitCone = 0;
   phoMIPIsHalo = 0;
   phoxtalBits = 0;
   phoIDbit = 0;
   phoScale_stat_up = 0;
   phoScale_stat_dn = 0;
   phoScale_syst_up = 0;
   phoScale_syst_dn = 0;
   phoScale_gain_up = 0;
   phoScale_gain_dn = 0;
   phoResol_rho_up = 0;
   phoResol_rho_dn = 0;
   phoResol_phi_up = 0;
   phoResol_phi_dn = 0;
   eleCharge = 0;
   eleChargeConsistent = 0;
   eleEn = 0;
   eleSCEn = 0;
   eleEcalEn = 0;
   eleESEnP1 = 0;
   eleESEnP2 = 0;
   eleD0 = 0;
   eleDz = 0;
   eleSIP = 0;
   elePt = 0;
   elePtError = 0;
   eleEta = 0;
   elePhi = 0;
   eleR9 = 0;
   eleCalibPt = 0;
   eleCalibEn = 0;
   eleSCEta = 0;
   eleSCPhi = 0;
   eleSCRawEn = 0;
   eleSCEtaWidth = 0;
   eleSCPhiWidth = 0;
   eleHoverE = 0;
   eleEoverP = 0;
   eleEoverPout = 0;
   eleEoverPInv = 0;
   eleBrem = 0;
   eledEtaAtVtx = 0;
   eledPhiAtVtx = 0;
   eleSigmaIEtaIEtaFull5x5 = 0;
   eleSigmaIPhiIPhiFull5x5 = 0;
   eleConvVeto = 0;
   eleMissHits = 0;
   eleESEffSigmaRR = 0;
   elePFChIso = 0;
   elePFPhoIso = 0;
   elePFNeuIso = 0;
   elePFPUIso = 0;
   elePFClusEcalIso = 0;
   elePFClusHcalIso = 0;
   eleIDMVAIso = 0;
   eleIDMVANoIso = 0;
   eleR9Full5x5 = 0;
   eleEcalDrivenSeed = 0;
   eleTrkdxy = 0;
   eleKFHits = 0;
   eleKFChi2 = 0;
   eleGSFChi2 = 0;
   eleGSFPt = 0;
   eleGSFEta = 0;
   eleGSFPhi = 0;
   eleGSFCharge = 0;
   eleGSFHits = 0;
   eleGSFMissHits = 0;
   eleGSFNHitsMax = 0;
   eleGSFVtxProb = 0;
   eleGSFlxyPV = 0;
   eleGSFlxyBS = 0;
   eleFiredSingleTrgs = 0;
   eleFiredDoubleTrgs = 0;
   eleFiredL1Trgs = 0;
   eleIDbit = 0;
   eleScale_stat_up = 0;
   eleScale_stat_dn = 0;
   eleScale_syst_up = 0;
   eleScale_syst_dn = 0;
   eleScale_gain_up = 0;
   eleScale_gain_dn = 0;
   eleResol_rho_up = 0;
   eleResol_rho_dn = 0;
   eleResol_phi_up = 0;
   eleResol_phi_dn = 0;
   muPt = 0;
   muEn = 0;
   muEta = 0;
   muPhi = 0;
   muCharge = 0;
   muType = 0;
   muIDbit = 0;
   muD0 = 0;
   muDz = 0;
   muSIP = 0;
   muChi2NDF = 0;
   muInnerD0 = 0;
   muInnerDz = 0;
   muTrkLayers = 0;
   muPixelLayers = 0;
   muPixelHits = 0;
   muMuonHits = 0;
   muStations = 0;
   muMatches = 0;
   muTrkQuality = 0;
   muIsoTrk = 0;
   muPFChIso = 0;
   muPFPhoIso = 0;
   muPFNeuIso = 0;
   muPFPUIso = 0;
   muPFChIso03 = 0;
   muPFPhoIso03 = 0;
   muPFNeuIso03 = 0;
   muPFPUIso03 = 0;
   muFiredTrgs = 0;
   muFiredL1Trgs = 0;
   muInnervalidFraction = 0;
   musegmentCompatibility = 0;
   muchi2LocalPosition = 0;
   mutrkKink = 0;
   muBestTrkPtError = 0;
   muBestTrkPt = 0;
   muBestTrkType = 0;
   pfPhoEt = 0;
   pfPhoEta = 0;
   pfPhoPhi = 0;
   pfPhoEtErr = 0;
   pfPhoChIso = 0;
   pfPhoPhoIso = 0;
   pfPhoNeuIso = 0;
   pfPhoEleIndex = 0;
   pfPhoPhoIndex = 0;
   jetPt = 0;
   jetEn = 0;
   jetEta = 0;
   jetPhi = 0;
   jetRawPt = 0;
   jetRawEn = 0;
   jetMt = 0;
   jetArea = 0;
   jetLeadTrackPt = 0;
   jetLeadTrackEta = 0;
   jetLeadTrackPhi = 0;
   jetLepTrackPID = 0;
   jetLepTrackPt = 0;
   jetLepTrackEta = 0;
   jetLepTrackPhi = 0;
   jetCSV2BJetTags = 0;
   jetDeepCSVTags_b = 0;
   jetDeepCSVTags_bb = 0;
   jetDeepCSVTags_c = 0;
   jetDeepCSVTags_udsg = 0;
   jetPartonID = 0;
   jetHadFlvr = 0;
   jetGenJetEn = 0;
   jetGenJetPt = 0;
   jetGenJetEta = 0;
   jetGenJetPhi = 0;
   jetGenPartonID = 0;
   jetGenEn = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenPhi = 0;
   jetGenPartonMomID = 0;
   jetP4Smear = 0;
   jetP4SmearUp = 0;
   jetP4SmearDo = 0;
   jetPFLooseId = 0;
   jetID = 0;
   jetPUID = 0;
   jetPUFullID = 0;
   jetJECUnc = 0;
   jetFiredTrgs = 0;
   jetCHF = 0;
   jetNHF = 0;
   jetCEF = 0;
   jetNEF = 0;
   jetNCH = 0;
   jetNNP = 0;
   jetMUF = 0;
   jetVtxPt = 0;
   jetVtxMass = 0;
   jetVtxNtrks = 0;
   jetVtx3DVal = 0;
   jetVtx3DSig = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumis", &lumis, &b_lumis);
   fChain->SetBranchAddress("isData", &isData, &b_isData);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("nGoodVtx", &nGoodVtx, &b_nGoodVtx);
   fChain->SetBranchAddress("isPVGood", &isPVGood, &b_isPVGood);
   fChain->SetBranchAddress("vtx", &vtx, &b_vtx);
   fChain->SetBranchAddress("vty", &vty, &b_vty);
   fChain->SetBranchAddress("vtz", &vtz, &b_vtz);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("rhoCentral", &rhoCentral, &b_rhoCentral);
   fChain->SetBranchAddress("L1ECALPrefire", &L1ECALPrefire, &b_L1ECALPrefire);
   fChain->SetBranchAddress("L1ECALPrefireUp", &L1ECALPrefireUp, &b_L1ECALPrefireUp);
   fChain->SetBranchAddress("L1ECALPrefireDown", &L1ECALPrefireDown, &b_L1ECALPrefireDown);
   fChain->SetBranchAddress("HLTEleMuX", &HLTEleMuX, &b_HLTEleMuX);
   fChain->SetBranchAddress("HLTPho", &HLTPho, &b_HLTPho);
   fChain->SetBranchAddress("HLTPhoRejectedByPS", &HLTPhoRejectedByPS, &b_HLTPhoRejectedByPS);
   fChain->SetBranchAddress("HLTJet", &HLTJet, &b_HLTJet);
   fChain->SetBranchAddress("HLTEleMuXIsPrescaled", &HLTEleMuXIsPrescaled, &b_HLTEleMuXIsPrescaled);
   fChain->SetBranchAddress("HLTPhoIsPrescaled", &HLTPhoIsPrescaled, &b_HLTPhoIsPrescaled);
   fChain->SetBranchAddress("HLTJetIsPrescaled", &HLTJetIsPrescaled, &b_HLTJetIsPrescaled);
   fChain->SetBranchAddress("pdf", &pdf, &b_pdf);
   fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
   fChain->SetBranchAddress("processID", &processID, &b_processID);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("genHT", &genHT, &b_genHT);
   fChain->SetBranchAddress("genPho1", &genPho1, &b_genPho1);
   fChain->SetBranchAddress("genPho2", &genPho2, &b_genPho2);
   fChain->SetBranchAddress("EventTag", &EventTag, &b_EventTag);
   fChain->SetBranchAddress("nPUInfo", &nPUInfo, &b_nPUInfo);
   fChain->SetBranchAddress("nPU", &nPU, &b_nPU);
   fChain->SetBranchAddress("puBX", &puBX, &b_puBX);
   fChain->SetBranchAddress("puTrue", &puTrue, &b_puTrue);
   fChain->SetBranchAddress("nLHE", &nLHE, &b_nLHE);
   fChain->SetBranchAddress("lhePID", &lhePID, &b_lhePID);
   fChain->SetBranchAddress("lhePx", &lhePx, &b_lhePx);
   fChain->SetBranchAddress("lhePy", &lhePy, &b_lhePy);
   fChain->SetBranchAddress("lhePz", &lhePz, &b_lhePz);
   fChain->SetBranchAddress("lheE", &lheE, &b_lheE);
   fChain->SetBranchAddress("nMC", &nMC, &b_nMC);
   fChain->SetBranchAddress("mcPID", &mcPID, &b_mcPID);
   fChain->SetBranchAddress("mcVtx", &mcVtx, &b_mcVtx);
   fChain->SetBranchAddress("mcVty", &mcVty, &b_mcVty);
   fChain->SetBranchAddress("mcVtz", &mcVtz, &b_mcVtz);
   fChain->SetBranchAddress("mcPt", &mcPt, &b_mcPt);
   fChain->SetBranchAddress("mcMass", &mcMass, &b_mcMass);
   fChain->SetBranchAddress("mcEta", &mcEta, &b_mcEta);
   fChain->SetBranchAddress("mcPhi", &mcPhi, &b_mcPhi);
   fChain->SetBranchAddress("mcE", &mcE, &b_mcE);
   fChain->SetBranchAddress("mcEt", &mcEt, &b_mcEt);
   fChain->SetBranchAddress("mcGMomPID", &mcGMomPID, &b_mcGMomPID);
   fChain->SetBranchAddress("mcMomPID", &mcMomPID, &b_mcMomPID);
   fChain->SetBranchAddress("mcMomPt", &mcMomPt, &b_mcMomPt);
   fChain->SetBranchAddress("mcMomMass", &mcMomMass, &b_mcMomMass);
   fChain->SetBranchAddress("mcMomEta", &mcMomEta, &b_mcMomEta);
   fChain->SetBranchAddress("mcMomPhi", &mcMomPhi, &b_mcMomPhi);
   fChain->SetBranchAddress("mcStatusFlag", &mcStatusFlag, &b_mcStatusFlag);
   fChain->SetBranchAddress("mcParentage", &mcParentage, &b_mcParentage);
   fChain->SetBranchAddress("mcStatus", &mcStatus, &b_mcStatus);
   fChain->SetBranchAddress("mcCalIsoDR03", &mcCalIsoDR03, &b_mcCalIsoDR03);
   fChain->SetBranchAddress("mcTrkIsoDR03", &mcTrkIsoDR03, &b_mcTrkIsoDR03);
   fChain->SetBranchAddress("mcCalIsoDR04", &mcCalIsoDR04, &b_mcCalIsoDR04);
   fChain->SetBranchAddress("mcTrkIsoDR04", &mcTrkIsoDR04, &b_mcTrkIsoDR04);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMETPhi", &genMETPhi, &b_genMETPhi);
   fChain->SetBranchAddress("genMET_pt", &genMET_pt, &b_genMET_pt);
   fChain->SetBranchAddress("genMET_px", &genMET_px, &b_genMET_px);
   fChain->SetBranchAddress("genMET_py", &genMET_py, &b_genMET_py);
   fChain->SetBranchAddress("genMET_pz", &genMET_pz, &b_genMET_pz);
   fChain->SetBranchAddress("genMET_sumEt", &genMET_sumEt, &b_genMET_sumEt);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("pfMET_pt", &pfMET_pt, &b_pfMET_pt);
   fChain->SetBranchAddress("pfMET_px", &pfMET_px, &b_pfMET_px);
   fChain->SetBranchAddress("pfMET_py", &pfMET_py, &b_pfMET_py);
   fChain->SetBranchAddress("pfMET_pz", &pfMET_pz, &b_pfMET_pz);
   fChain->SetBranchAddress("pfMET_sumEt", &pfMET_sumEt, &b_pfMET_sumEt);
   fChain->SetBranchAddress("rawMET_pt", &rawMET_pt, &b_rawMET_pt);
   fChain->SetBranchAddress("rawMET_px", &rawMET_px, &b_rawMET_px);
   fChain->SetBranchAddress("rawMET_py", &rawMET_py, &b_rawMET_py);
   fChain->SetBranchAddress("rawMETPhi", &rawMETPhi, &b_rawMETPhi);
   fChain->SetBranchAddress("puppiMET", &puppiMET, &b_puppiMET);
   fChain->SetBranchAddress("puppiMETPhi", &puppiMETPhi, &b_puppiMETPhi);
   fChain->SetBranchAddress("puppiMET_pt", &puppiMET_pt, &b_puppiMET_pt);
   fChain->SetBranchAddress("puppiMET_px", &puppiMET_px, &b_puppiMET_px);
   fChain->SetBranchAddress("puppiMET_py", &puppiMET_py, &b_puppiMET_py);
   fChain->SetBranchAddress("puppiMET_pz", &puppiMET_pz, &b_puppiMET_pz);
   fChain->SetBranchAddress("puppiMET_sumEt", &puppiMET_sumEt, &b_puppiMET_sumEt);
   fChain->SetBranchAddress("rawPuppiMET_pt", &rawPuppiMET_pt, &b_rawPuppiMET_pt);
   fChain->SetBranchAddress("rawPuppiMET_px", &rawPuppiMET_px, &b_rawPuppiMET_px);
   fChain->SetBranchAddress("rawPuppiMET_py", &rawPuppiMET_py, &b_rawPuppiMET_py);
   fChain->SetBranchAddress("rawPuppiMETPhi", &rawPuppiMETPhi, &b_rawPuppiMETPhi);
   fChain->SetBranchAddress("pfMET_T1JERUp", &pfMET_T1JERUp, &b_pfMET_T1JERUp);
   fChain->SetBranchAddress("pfMET_T1JERDo", &pfMET_T1JERDo, &b_pfMET_T1JERDo);
   fChain->SetBranchAddress("pfMET_T1JESUp", &pfMET_T1JESUp, &b_pfMET_T1JESUp);
   fChain->SetBranchAddress("pfMET_T1JESDo", &pfMET_T1JESDo, &b_pfMET_T1JESDo);
   fChain->SetBranchAddress("pfMET_T1UESUp", &pfMET_T1UESUp, &b_pfMET_T1UESUp);
   fChain->SetBranchAddress("pfMET_T1UESDo", &pfMET_T1UESDo, &b_pfMET_T1UESDo);
   fChain->SetBranchAddress("pfMETPhi_T1JESUp", &pfMETPhi_T1JESUp, &b_pfMETPhi_T1JESUp);
   fChain->SetBranchAddress("pfMETPhi_T1JESDo", &pfMETPhi_T1JESDo, &b_pfMETPhi_T1JESDo);
   fChain->SetBranchAddress("pfMETPhi_T1UESUp", &pfMETPhi_T1UESUp, &b_pfMETPhi_T1UESUp);
   fChain->SetBranchAddress("pfMETPhi_T1UESDo", &pfMETPhi_T1UESDo, &b_pfMETPhi_T1UESDo);
   fChain->SetBranchAddress("nPho", &nPho, &b_nPho);
   fChain->SetBranchAddress("phoE", &phoE, &b_phoE);
   fChain->SetBranchAddress("phoSigmaE", &phoSigmaE, &b_phoSigmaE);
   fChain->SetBranchAddress("phoEt", &phoEt, &b_phoEt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoCalibE", &phoCalibE, &b_phoCalibE);
   fChain->SetBranchAddress("phoCalibEt", &phoCalibEt, &b_phoCalibEt);
   fChain->SetBranchAddress("phoSCE", &phoSCE, &b_phoSCE);
   fChain->SetBranchAddress("phoSCRawE", &phoSCRawE, &b_phoSCRawE);
   fChain->SetBranchAddress("phoESEnP1", &phoESEnP1, &b_phoESEnP1);
   fChain->SetBranchAddress("phoESEnP2", &phoESEnP2, &b_phoESEnP2);
   fChain->SetBranchAddress("phoSCEta", &phoSCEta, &b_phoSCEta);
   fChain->SetBranchAddress("phoSCPhi", &phoSCPhi, &b_phoSCPhi);
   fChain->SetBranchAddress("phoSCEtaWidth", &phoSCEtaWidth, &b_phoSCEtaWidth);
   fChain->SetBranchAddress("phoSCPhiWidth", &phoSCPhiWidth, &b_phoSCPhiWidth);
   fChain->SetBranchAddress("phoSCBrem", &phoSCBrem, &b_phoSCBrem);
   fChain->SetBranchAddress("phohasPixelSeed", &phohasPixelSeed, &b_phohasPixelSeed);
   fChain->SetBranchAddress("phoEleVeto", &phoEleVeto, &b_phoEleVeto);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoHoverE", &phoHoverE, &b_phoHoverE);
   fChain->SetBranchAddress("phoESEffSigmaRR", &phoESEffSigmaRR, &b_phoESEffSigmaRR);
   fChain->SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5, &b_phoSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5, &b_phoSigmaIEtaIPhiFull5x5);
   fChain->SetBranchAddress("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5, &b_phoSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("phoE2x2Full5x5", &phoE2x2Full5x5, &b_phoE2x2Full5x5);
   fChain->SetBranchAddress("phoE5x5Full5x5", &phoE5x5Full5x5, &b_phoE5x5Full5x5);
   fChain->SetBranchAddress("phoR9Full5x5", &phoR9Full5x5, &b_phoR9Full5x5);
   fChain->SetBranchAddress("phoPFChIso", &phoPFChIso, &b_phoPFChIso);
   fChain->SetBranchAddress("phoPFChPVIso", &phoPFChPVIso, &b_phoPFChPVIso);
   fChain->SetBranchAddress("phoPFPhoIso", &phoPFPhoIso, &b_phoPFPhoIso);
   fChain->SetBranchAddress("phoPFNeuIso", &phoPFNeuIso, &b_phoPFNeuIso);
   fChain->SetBranchAddress("phoPFChWorstIso", &phoPFChWorstIso, &b_phoPFChWorstIso);
   fChain->SetBranchAddress("phoPFChWorstVetoIso", &phoPFChWorstVetoIso, &b_phoPFChWorstVetoIso);
   fChain->SetBranchAddress("phoEcalPFClusterIso", &phoEcalPFClusterIso, &b_phoEcalPFClusterIso);
   fChain->SetBranchAddress("phoHcalPFClusterIso", &phoHcalPFClusterIso, &b_phoHcalPFClusterIso);
   fChain->SetBranchAddress("phoIDMVA", &phoIDMVA, &b_phoIDMVA);
   fChain->SetBranchAddress("phoFiredSingleTrgs", &phoFiredSingleTrgs, &b_phoFiredSingleTrgs);
   fChain->SetBranchAddress("phoFiredDoubleTrgs", &phoFiredDoubleTrgs, &b_phoFiredDoubleTrgs);
   fChain->SetBranchAddress("phoFiredTripleTrgs", &phoFiredTripleTrgs, &b_phoFiredTripleTrgs);
   fChain->SetBranchAddress("phoFiredL1Trgs", &phoFiredL1Trgs, &b_phoFiredL1Trgs);
   fChain->SetBranchAddress("phoSeedTime", &phoSeedTime, &b_phoSeedTime);
   fChain->SetBranchAddress("phoSeedEnergy", &phoSeedEnergy, &b_phoSeedEnergy);
   fChain->SetBranchAddress("phoMIPTotEnergy", &phoMIPTotEnergy, &b_phoMIPTotEnergy);
   fChain->SetBranchAddress("phoMIPChi2", &phoMIPChi2, &b_phoMIPChi2);
   fChain->SetBranchAddress("phoMIPSlope", &phoMIPSlope, &b_phoMIPSlope);
   fChain->SetBranchAddress("phoMIPIntercept", &phoMIPIntercept, &b_phoMIPIntercept);
   fChain->SetBranchAddress("phoMIPNhitCone", &phoMIPNhitCone, &b_phoMIPNhitCone);
   fChain->SetBranchAddress("phoMIPIsHalo", &phoMIPIsHalo, &b_phoMIPIsHalo);
   fChain->SetBranchAddress("phoxtalBits", &phoxtalBits, &b_phoxtalBits);
   fChain->SetBranchAddress("phoIDbit", &phoIDbit, &b_phoIDbit);
   fChain->SetBranchAddress("phoScale_stat_up", &phoScale_stat_up, &b_phoScale_stat_up);
   fChain->SetBranchAddress("phoScale_stat_dn", &phoScale_stat_dn, &b_phoScale_stat_dn);
   fChain->SetBranchAddress("phoScale_syst_up", &phoScale_syst_up, &b_phoScale_syst_up);
   fChain->SetBranchAddress("phoScale_syst_dn", &phoScale_syst_dn, &b_phoScale_syst_dn);
   fChain->SetBranchAddress("phoScale_gain_up", &phoScale_gain_up, &b_phoScale_gain_up);
   fChain->SetBranchAddress("phoScale_gain_dn", &phoScale_gain_dn, &b_phoScale_gain_dn);
   fChain->SetBranchAddress("phoResol_rho_up", &phoResol_rho_up, &b_phoResol_rho_up);
   fChain->SetBranchAddress("phoResol_rho_dn", &phoResol_rho_dn, &b_phoResol_rho_dn);
   fChain->SetBranchAddress("phoResol_phi_up", &phoResol_phi_up, &b_phoResol_phi_up);
   fChain->SetBranchAddress("phoResol_phi_dn", &phoResol_phi_dn, &b_phoResol_phi_dn);
   fChain->SetBranchAddress("nEle", &nEle, &b_nEle);
   fChain->SetBranchAddress("eleCharge", &eleCharge, &b_eleCharge);
   fChain->SetBranchAddress("eleChargeConsistent", &eleChargeConsistent, &b_eleChargeConsistent);
   fChain->SetBranchAddress("eleEn", &eleEn, &b_eleEn);
   fChain->SetBranchAddress("eleSCEn", &eleSCEn, &b_eleSCEn);
   fChain->SetBranchAddress("eleEcalEn", &eleEcalEn, &b_eleEcalEn);
   fChain->SetBranchAddress("eleESEnP1", &eleESEnP1, &b_eleESEnP1);
   fChain->SetBranchAddress("eleESEnP2", &eleESEnP2, &b_eleESEnP2);
   fChain->SetBranchAddress("eleD0", &eleD0, &b_eleD0);
   fChain->SetBranchAddress("eleDz", &eleDz, &b_eleDz);
   fChain->SetBranchAddress("eleSIP", &eleSIP, &b_eleSIP);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("elePtError", &elePtError, &b_elePtError);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleR9", &eleR9, &b_eleR9);
   fChain->SetBranchAddress("eleCalibPt", &eleCalibPt, &b_eleCalibPt);
   fChain->SetBranchAddress("eleCalibEn", &eleCalibEn, &b_eleCalibEn);
   fChain->SetBranchAddress("eleSCEta", &eleSCEta, &b_eleSCEta);
   fChain->SetBranchAddress("eleSCPhi", &eleSCPhi, &b_eleSCPhi);
   fChain->SetBranchAddress("eleSCRawEn", &eleSCRawEn, &b_eleSCRawEn);
   fChain->SetBranchAddress("eleSCEtaWidth", &eleSCEtaWidth, &b_eleSCEtaWidth);
   fChain->SetBranchAddress("eleSCPhiWidth", &eleSCPhiWidth, &b_eleSCPhiWidth);
   fChain->SetBranchAddress("eleHoverE", &eleHoverE, &b_eleHoverE);
   fChain->SetBranchAddress("eleEoverP", &eleEoverP, &b_eleEoverP);
   fChain->SetBranchAddress("eleEoverPout", &eleEoverPout, &b_eleEoverPout);
   fChain->SetBranchAddress("eleEoverPInv", &eleEoverPInv, &b_eleEoverPInv);
   fChain->SetBranchAddress("eleBrem", &eleBrem, &b_eleBrem);
   fChain->SetBranchAddress("eledEtaAtVtx", &eledEtaAtVtx, &b_eledEtaAtVtx);
   fChain->SetBranchAddress("eledPhiAtVtx", &eledPhiAtVtx, &b_eledPhiAtVtx);
   fChain->SetBranchAddress("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5, &b_eleSigmaIEtaIEtaFull5x5);
   fChain->SetBranchAddress("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5, &b_eleSigmaIPhiIPhiFull5x5);
   fChain->SetBranchAddress("eleConvVeto", &eleConvVeto, &b_eleConvVeto);
   fChain->SetBranchAddress("eleMissHits", &eleMissHits, &b_eleMissHits);
   fChain->SetBranchAddress("eleESEffSigmaRR", &eleESEffSigmaRR, &b_eleESEffSigmaRR);
   fChain->SetBranchAddress("elePFChIso", &elePFChIso, &b_elePFChIso);
   fChain->SetBranchAddress("elePFPhoIso", &elePFPhoIso, &b_elePFPhoIso);
   fChain->SetBranchAddress("elePFNeuIso", &elePFNeuIso, &b_elePFNeuIso);
   fChain->SetBranchAddress("elePFPUIso", &elePFPUIso, &b_elePFPUIso);
   fChain->SetBranchAddress("elePFClusEcalIso", &elePFClusEcalIso, &b_elePFClusEcalIso);
   fChain->SetBranchAddress("elePFClusHcalIso", &elePFClusHcalIso, &b_elePFClusHcalIso);
   fChain->SetBranchAddress("eleIDMVAIso", &eleIDMVAIso, &b_eleIDMVAIso);
   fChain->SetBranchAddress("eleIDMVANoIso", &eleIDMVANoIso, &b_eleIDMVANoIso);
   fChain->SetBranchAddress("eleR9Full5x5", &eleR9Full5x5, &b_eleR9Full5x5);
   fChain->SetBranchAddress("eleEcalDrivenSeed", &eleEcalDrivenSeed, &b_eleEcalDrivenSeed);
   fChain->SetBranchAddress("eleTrkdxy", &eleTrkdxy, &b_eleTrkdxy);
   fChain->SetBranchAddress("eleKFHits", &eleKFHits, &b_eleKFHits);
   fChain->SetBranchAddress("eleKFChi2", &eleKFChi2, &b_eleKFChi2);
   fChain->SetBranchAddress("eleGSFChi2", &eleGSFChi2, &b_eleGSFChi2);
   fChain->SetBranchAddress("eleGSFPt", &eleGSFPt, &b_eleGSFPt);
   fChain->SetBranchAddress("eleGSFEta", &eleGSFEta, &b_eleGSFEta);
   fChain->SetBranchAddress("eleGSFPhi", &eleGSFPhi, &b_eleGSFPhi);
   fChain->SetBranchAddress("eleGSFCharge", &eleGSFCharge, &b_eleGSFCharge);
   fChain->SetBranchAddress("eleGSFHits", &eleGSFHits, &b_eleGSFHits);
   fChain->SetBranchAddress("eleGSFMissHits", &eleGSFMissHits, &b_eleGSFMissHits);
   fChain->SetBranchAddress("eleGSFNHitsMax", &eleGSFNHitsMax, &b_eleGSFNHitsMax);
   fChain->SetBranchAddress("eleGSFVtxProb", &eleGSFVtxProb, &b_eleGSFVtxProb);
   fChain->SetBranchAddress("eleGSFlxyPV", &eleGSFlxyPV, &b_eleGSFlxyPV);
   fChain->SetBranchAddress("eleGSFlxyBS", &eleGSFlxyBS, &b_eleGSFlxyBS);
   fChain->SetBranchAddress("eleFiredSingleTrgs", &eleFiredSingleTrgs, &b_eleFiredSingleTrgs);
   fChain->SetBranchAddress("eleFiredDoubleTrgs", &eleFiredDoubleTrgs, &b_eleFiredDoubleTrgs);
   fChain->SetBranchAddress("eleFiredL1Trgs", &eleFiredL1Trgs, &b_eleFiredL1Trgs);
   fChain->SetBranchAddress("eleIDbit", &eleIDbit, &b_eleIDbit);
   fChain->SetBranchAddress("eleScale_stat_up", &eleScale_stat_up, &b_eleScale_stat_up);
   fChain->SetBranchAddress("eleScale_stat_dn", &eleScale_stat_dn, &b_eleScale_stat_dn);
   fChain->SetBranchAddress("eleScale_syst_up", &eleScale_syst_up, &b_eleScale_syst_up);
   fChain->SetBranchAddress("eleScale_syst_dn", &eleScale_syst_dn, &b_eleScale_syst_dn);
   fChain->SetBranchAddress("eleScale_gain_up", &eleScale_gain_up, &b_eleScale_gain_up);
   fChain->SetBranchAddress("eleScale_gain_dn", &eleScale_gain_dn, &b_eleScale_gain_dn);
   fChain->SetBranchAddress("eleResol_rho_up", &eleResol_rho_up, &b_eleResol_rho_up);
   fChain->SetBranchAddress("eleResol_rho_dn", &eleResol_rho_dn, &b_eleResol_rho_dn);
   fChain->SetBranchAddress("eleResol_phi_up", &eleResol_phi_up, &b_eleResol_phi_up);
   fChain->SetBranchAddress("eleResol_phi_dn", &eleResol_phi_dn, &b_eleResol_phi_dn);
   fChain->SetBranchAddress("nMu", &nMu, &b_nMu);
   fChain->SetBranchAddress("muPt", &muPt, &b_muPt);
   fChain->SetBranchAddress("muEn", &muEn, &b_muEn);
   fChain->SetBranchAddress("muEta", &muEta, &b_muEta);
   fChain->SetBranchAddress("muPhi", &muPhi, &b_muPhi);
   fChain->SetBranchAddress("muCharge", &muCharge, &b_muCharge);
   fChain->SetBranchAddress("muType", &muType, &b_muType);
   fChain->SetBranchAddress("muIDbit", &muIDbit, &b_muIDbit);
   fChain->SetBranchAddress("muD0", &muD0, &b_muD0);
   fChain->SetBranchAddress("muDz", &muDz, &b_muDz);
   fChain->SetBranchAddress("muSIP", &muSIP, &b_muSIP);
   fChain->SetBranchAddress("muChi2NDF", &muChi2NDF, &b_muChi2NDF);
   fChain->SetBranchAddress("muInnerD0", &muInnerD0, &b_muInnerD0);
   fChain->SetBranchAddress("muInnerDz", &muInnerDz, &b_muInnerDz);
   fChain->SetBranchAddress("muTrkLayers", &muTrkLayers, &b_muTrkLayers);
   fChain->SetBranchAddress("muPixelLayers", &muPixelLayers, &b_muPixelLayers);
   fChain->SetBranchAddress("muPixelHits", &muPixelHits, &b_muPixelHits);
   fChain->SetBranchAddress("muMuonHits", &muMuonHits, &b_muMuonHits);
   fChain->SetBranchAddress("muStations", &muStations, &b_muStations);
   fChain->SetBranchAddress("muMatches", &muMatches, &b_muMatches);
   fChain->SetBranchAddress("muTrkQuality", &muTrkQuality, &b_muTrkQuality);
   fChain->SetBranchAddress("muIsoTrk", &muIsoTrk, &b_muIsoTrk);
   fChain->SetBranchAddress("muPFChIso", &muPFChIso, &b_muPFChIso);
   fChain->SetBranchAddress("muPFPhoIso", &muPFPhoIso, &b_muPFPhoIso);
   fChain->SetBranchAddress("muPFNeuIso", &muPFNeuIso, &b_muPFNeuIso);
   fChain->SetBranchAddress("muPFPUIso", &muPFPUIso, &b_muPFPUIso);
   fChain->SetBranchAddress("muPFChIso03", &muPFChIso03, &b_muPFChIso03);
   fChain->SetBranchAddress("muPFPhoIso03", &muPFPhoIso03, &b_muPFPhoIso03);
   fChain->SetBranchAddress("muPFNeuIso03", &muPFNeuIso03, &b_muPFNeuIso03);
   fChain->SetBranchAddress("muPFPUIso03", &muPFPUIso03, &b_muPFPUIso03);
   fChain->SetBranchAddress("muFiredTrgs", &muFiredTrgs, &b_muFiredTrgs);
   fChain->SetBranchAddress("muFiredL1Trgs", &muFiredL1Trgs, &b_muFiredL1Trgs);
   fChain->SetBranchAddress("muInnervalidFraction", &muInnervalidFraction, &b_muInnervalidFraction);
   fChain->SetBranchAddress("musegmentCompatibility", &musegmentCompatibility, &b_musegmentCompatibility);
   fChain->SetBranchAddress("muchi2LocalPosition", &muchi2LocalPosition, &b_muchi2LocalPosition);
   fChain->SetBranchAddress("mutrkKink", &mutrkKink, &b_mutrkKink);
   fChain->SetBranchAddress("muBestTrkPtError", &muBestTrkPtError, &b_muBestTrkPtError);
   fChain->SetBranchAddress("muBestTrkPt", &muBestTrkPt, &b_muBestTrkPt);
   fChain->SetBranchAddress("muBestTrkType", &muBestTrkType, &b_muBestTrkType);
   fChain->SetBranchAddress("npfPho", &npfPho, &b_npfPho);
   fChain->SetBranchAddress("pfPhoEt", &pfPhoEt, &b_pfPhoEt);
   fChain->SetBranchAddress("pfPhoEta", &pfPhoEta, &b_pfPhoEta);
   fChain->SetBranchAddress("pfPhoPhi", &pfPhoPhi, &b_pfPhoPhi);
   fChain->SetBranchAddress("pfPhoEtErr", &pfPhoEtErr, &b_pfPhoEtErr);
   fChain->SetBranchAddress("pfPhoChIso", &pfPhoChIso, &b_pfPhoChIso);
   fChain->SetBranchAddress("pfPhoPhoIso", &pfPhoPhoIso, &b_pfPhoPhoIso);
   fChain->SetBranchAddress("pfPhoNeuIso", &pfPhoNeuIso, &b_pfPhoNeuIso);
   fChain->SetBranchAddress("pfPhoEleIndex", &pfPhoEleIndex, &b_pfPhoEleIndex);
   fChain->SetBranchAddress("pfPhoPhoIndex", &pfPhoPhoIndex, &b_pfPhoPhoIndex);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEn", &jetEn, &b_jetEn);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetRawPt", &jetRawPt, &b_jetRawPt);
   fChain->SetBranchAddress("jetRawEn", &jetRawEn, &b_jetRawEn);
   fChain->SetBranchAddress("jetMt", &jetMt, &b_jetMt);
   fChain->SetBranchAddress("jetArea", &jetArea, &b_jetArea);
   fChain->SetBranchAddress("jetLeadTrackPt", &jetLeadTrackPt, &b_jetLeadTrackPt);
   fChain->SetBranchAddress("jetLeadTrackEta", &jetLeadTrackEta, &b_jetLeadTrackEta);
   fChain->SetBranchAddress("jetLeadTrackPhi", &jetLeadTrackPhi, &b_jetLeadTrackPhi);
   fChain->SetBranchAddress("jetLepTrackPID", &jetLepTrackPID, &b_jetLepTrackPID);
   fChain->SetBranchAddress("jetLepTrackPt", &jetLepTrackPt, &b_jetLepTrackPt);
   fChain->SetBranchAddress("jetLepTrackEta", &jetLepTrackEta, &b_jetLepTrackEta);
   fChain->SetBranchAddress("jetLepTrackPhi", &jetLepTrackPhi, &b_jetLepTrackPhi);
   fChain->SetBranchAddress("jetCSV2BJetTags", &jetCSV2BJetTags, &b_jetCSV2BJetTags);
   fChain->SetBranchAddress("jetDeepCSVTags_b", &jetDeepCSVTags_b, &b_jetDeepCSVTags_b);
   fChain->SetBranchAddress("jetDeepCSVTags_bb", &jetDeepCSVTags_bb, &b_jetDeepCSVTags_bb);
   fChain->SetBranchAddress("jetDeepCSVTags_c", &jetDeepCSVTags_c, &b_jetDeepCSVTags_c);
   fChain->SetBranchAddress("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg, &b_jetDeepCSVTags_udsg);
   fChain->SetBranchAddress("jetPartonID", &jetPartonID, &b_jetPartonID);
   fChain->SetBranchAddress("jetHadFlvr", &jetHadFlvr, &b_jetHadFlvr);
   fChain->SetBranchAddress("jetGenJetEn", &jetGenJetEn, &b_jetGenJetEn);
   fChain->SetBranchAddress("jetGenJetPt", &jetGenJetPt, &b_jetGenJetPt);
   fChain->SetBranchAddress("jetGenJetEta", &jetGenJetEta, &b_jetGenJetEta);
   fChain->SetBranchAddress("jetGenJetPhi", &jetGenJetPhi, &b_jetGenJetPhi);
   fChain->SetBranchAddress("jetGenPartonID", &jetGenPartonID, &b_jetGenPartonID);
   fChain->SetBranchAddress("jetGenEn", &jetGenEn, &b_jetGenEn);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenPhi", &jetGenPhi, &b_jetGenPhi);
   fChain->SetBranchAddress("jetGenPartonMomID", &jetGenPartonMomID, &b_jetGenPartonMomID);
   fChain->SetBranchAddress("jetP4Smear", &jetP4Smear, &b_jetP4Smear);
   fChain->SetBranchAddress("jetP4SmearUp", &jetP4SmearUp, &b_jetP4SmearUp);
   fChain->SetBranchAddress("jetP4SmearDo", &jetP4SmearDo, &b_jetP4SmearDo);
   fChain->SetBranchAddress("jetPFLooseId", &jetPFLooseId, &b_jetPFLooseId);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("jetPUID", &jetPUID, &b_jetPUID);
   fChain->SetBranchAddress("jetPUFullID", &jetPUFullID, &b_jetPUFullID);
   fChain->SetBranchAddress("jetJECUnc", &jetJECUnc, &b_jetJECUnc);
   fChain->SetBranchAddress("jetFiredTrgs", &jetFiredTrgs, &b_jetFiredTrgs);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetCEF", &jetCEF, &b_jetCEF);
   fChain->SetBranchAddress("jetNEF", &jetNEF, &b_jetNEF);
   fChain->SetBranchAddress("jetNCH", &jetNCH, &b_jetNCH);
   fChain->SetBranchAddress("jetNNP", &jetNNP, &b_jetNNP);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
   fChain->SetBranchAddress("jetVtxPt", &jetVtxPt, &b_jetVtxPt);
   fChain->SetBranchAddress("jetVtxMass", &jetVtxMass, &b_jetVtxMass);
   fChain->SetBranchAddress("jetVtxNtrks", &jetVtxNtrks, &b_jetVtxNtrks);
   fChain->SetBranchAddress("jetVtx3DVal", &jetVtx3DVal, &b_jetVtx3DVal);
   fChain->SetBranchAddress("jetVtx3DSig", &jetVtx3DSig, &b_jetVtx3DSig);
   Notify();
}

Bool_t jetmet::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void jetmet::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t jetmet::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef jetmet_cxx

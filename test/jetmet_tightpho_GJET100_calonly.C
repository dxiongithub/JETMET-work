//04.29.2020
//calculating the recoil of each MET (PF, PUPPI, Raw)
//for service work - JETMET group
//ref: https://github.com/Soumyatifr/nanoAOD-tools/blob/MET_NANOAOD/python/postprocessing/examples/MCeventselection.py#L84

//GJET dataset

/*
pfmet_px= pfmet.pt*m.cos(pfmet.phi)
pfmet_py= pfmet.pt*m.sin(pfmet.phi)
puppimet_px= puppimet.pt*m.cos(puppimet.phi)
puppimet_py= puppimet.pt*m.sin(puppimet.phi)
rawmet_px= rawmet.pt*m.cos(rawmet.phi)
rawmet_py= rawmet.pt*m.sin(rawmet.phi)


uparpendicular_pf= ((-pfmet_px-photon_px)*photon_py-(-pfmet_py-photon_py)*photon_px)/tp.pt
uparallal_pf= ((-pfmet_px-photon_px)*photon_px+(-pfmet_py-photon_py)*photon_py)/tp.pt
uparpendicular_puppi= ((-puppimet_px-photon_px)*photon_py-(-puppimet_py-photon_py)*photon_px)/tp.pt
uparallal_puppi= ((-puppimet_px-photon_px)*photon_px+(-puppimet_py-photon_py)*photon_py)/tp.pt
uparpendicular_raw= ((-rawmet_px-photon_px)*photon_py-(-rawmet_py-photon_py)*photon_px)/tp.pt
uparallal_raw= ((-rawmet_px-photon_px)*photon_px+(-rawmet_py-photon_py)*photon_py)/tp.pt
*/



#define jetmet_cxx
#include "jetmet.h"
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TMath.h>
#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>


using namespace std;

void jetmet::Loop()
{
//   In a ROOT session, you can do:
//      root> .L jetmet.C
//      root> jetmet t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch



        //--------------------------------------------------------------------------------------------
        //start from here, we use only tight photon ID selection
        //https://github.com/Soumyatifr/nanoAOD-tools/blob/MET_NANOAOD/python/postprocessing/examples/MCeventselection.py#L69
        /*
        selection:
        fabs eta < 1.442
        phopt > 50
        phoIDbit == 2
		*/


        /*
        (a) Loose: phoIDbit[]>>0&1 ---> gives 0 or 1. if 0--> this phoID is failed. if 1--> this phoID is passed
        (b) Medium: phoIDbit[]>>1&1
        (c) Tight: phoIDbit[]>>2&1
        */

        //--------------------------------------------------------------------------------------------
        //start from here, we use the multiple ID selection
        //select events with exactly one tight photon, at least one jet, no loose electron or muon
        //https://github.com/Soumyatifr/nanoAOD-tools/blob/MET_NANOAOD/python/postprocessing/examples/MCeventselection.py#L69
        /*
        selections:

        exactly one photon:
        fabs eta < 1.442
        phopt > 50
        phoIDbit == 2

        at least one jet:
        jetpt > 40
        jet eta < 2.4
        jetID > 1 ==> tight jet ==> == 6 <==  += pow(2,2) and for loose += pow(2,1)

        electron:
        loose electron
        ele.pt > 10

        muon:
        loose muon
        mu.pt > 10

        trigger:
        Photon50_R9Id90_HE10_IsoM
        Photon75_R9Id90_HE10_IsoM 
        Photon90_R9Id90_HE10_IsoM 
        Photon120_R9Id90_HE10_IsoM 
        Photon165_R9Id90_HE10_IsoM 

        select:
        tightjet > 0 && onetightpho > 0 && looseelectron == 0 && loosemuon == 0 && 
        flag.BadPFMuonFilter == 1 (pow(2,6+3) = 2^9 = 512)
        flag.ecalBadCalibFilter == 1 (not in ggNtuplizer)
        flag.globalTightHalo2016Filter == 1 (pow(2, 2+1) = 2^3 = 8)
        flag.HBHENoiseFilter==1 (pow(2, 0+1) = 2^1 = 2)
        flag.HBHENoiseIsoFilter ==1 (pow(2, 1+1) = 2^2 = 4)
        flag.EcalDeadCellTriggerPrimitiveFilter==1 (pow(2, 5+1) = 2^6 = 64)
        flag.eeBadScFilter==1 (pow(2, 4+1) = 2^5 = 32)

        then calculate px py and per and pen
		*/

        int onetightPho = 0;

        for (int ipho = 0; ipho < nPho; ipho++)
        {
            file2 << (*phoIDbit)[ipho] << endl;
            
            if (fabs((*phoEta)[ipho]) < 1.442 && (*phoEt)[ipho] > 50 && ((*phoIDbit)[ipho] >> 2 & 1) && nPho == 1)
            {
                onetightPho ++;
            }
        }

        
        if (onetightPho > 0)
        {
            nPho_one->Fill(onetightPho);
        }
        

       int tightJet = 0;
        
        for (int ijet = 0; ijet < nJet; ijet++)
        {
            file3 << (*jetID)[ijet] << endl;

            if (fabs((*jetEta)[ijet]) < 2.4 && (*jetPt)[ijet] > 40 && (*jetID)[ijet] == 6)
            {
                tightJet ++;
            }
        }

        nJet_one->Fill(tightJet);


        int looseelectron = 0;

        for (int iele = 0; iele < nEle; iele++)
        {
            file4 << (*eleIDbit)[iele] << endl;

            if ((*elePt)[iele] > 10 && ((*eleIDbit)[iele] >> 1 & 1))
            {
                looseelectron ++;
            }
        }

        nEle_one->Fill(looseelectron);


        int loosemuon = 0;
        
        for (int imu = 0; imu < nMu; imu++)
        {
            file5 << (*muIDbit)[imu] << endl;

            if ((*muPt)[imu] > 10 && ((*muIDbit)[imu] >> 0 & 1))
            {
                loosemuon++;
            }
        }

        nMu_one->Fill(loosemuon);



        //start recoil calculation

        Float_t puppipho_px = 0;
        Float_t puppipho_py = 0;

        Float_t rawpuppipho_px = 0;
        Float_t rawpuppipho_py = 0;
       
        Float_t pfpho_px = 0;
        Float_t pfpho_py = 0;

        Float_t rawpho_px = 0;
        Float_t rawpho_py = 0;

        Float_t puppi_pen = 0;
        Float_t puppi_par = 0;

        Float_t rawpuppi_pen = 0;
        Float_t rawpuppi_par = 0;

        Float_t pf_pen = 0;
        Float_t pf_par = 0;

        Float_t raw_pen = 0;
        Float_t raw_par = 0;




        //selections
        if (((HLTPho >> 30 & 1) == 1) || //&& HLTPho < 2147483648  //HLT_Photon50_R9Id90_HE10_IsoM_v
        ((HLTPho >> 31 & 1) == 1) ||     //&& HLTPho < 4294967296  //HLT_Photon75_R9Id90_HE10_IsoM_v
        ((HLTPho >> 32 & 1) == 1) ||     //&& HLTPho < 8589934592  //HLT_Photon90_R9Id90_HE10_IsoM_v
        ((HLTPho >> 33 & 1) == 1) ||     //&& HLTPho < 17079869084 //HLT_Photon120_R9Id90_HE10_IsoM_v
        ((HLTPho >> 34 & 1) == 1))       //&& HLTPho < 34359738368 //HLT_Photon165_R9Id90_HE10_IsoM_v
        {
            if (tightJet > 0 && onetightPho > 0 && looseelectron == 0 && loosemuon == 0)
            /*
            flag.BadPFMuonFilter == 1 (pow(2,6+3) = 2^9 = 512)
            flag.ecalBadCalibFilter == 1 (not in ggNtuplizer)
            flag.globalTightHalo2016Filter == 1 (pow(2, 2+1) = 2^3 = 8)
            flag.HBHENoiseFilter==1 (pow(2, 0+1) = 2^1 = 2)
            flag.HBHENoiseIsoFilter ==1 (pow(2, 1+1) = 2^2 = 4)
            flag.EcalDeadCellTriggerPrimitiveFilter==1 (pow(2, 5+1) = 2^6 = 64)
            flag.eeBadScFilter==1 (pow(2, 4+1) = 2^5 = 32)
            */
            //metFilters == 512 && metFilters == 8 && metFilters == 2 && metFilters == 4 && metFilters == 64 && metFilters == 32)
            //There's no metFilters in MC config in ggNtuplizer
            {
                //start recoil calculation
                //puppimet

                /*
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                */


                for (int ipuppiPho = 0; ipuppiPho < nPho; ipuppiPho++)
                {
                    puppimetphi->Fill(puppiMETPhi);

                    puppimet_tightphopt->Fill((*phoEt)[ipuppiPho]);

                    //fill the px py from ggNtuplizer
                    puppimetpt->Fill(puppiMET_pt);
                    puppimet_px->Fill(puppiMET_px);
                    puppimet_py->Fill(puppiMET_py);

                    
                    //try to calculate and also fill the one from ggNtuplizer
                    //puppimet_px = puppimet.pt*m.cos(puppimet.phi)
                    puppimet_px_cal->Fill((puppiMET_pt)*cos(puppiMETPhi));
                    puppimet_py_cal->Fill((puppiMET_pt)*sin(puppiMETPhi));
                    //==> passed test! the same as ggNtuplizer stored.
                    

                    //cal photon_px and photon_py
                    puppipho_px = puppipho_px + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                    puppipho_py = puppipho_py + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                    puppipho_pxdist->Fill(puppipho_px);
                    puppipho_pydist->Fill(puppipho_py);
                    

                    //calculate the perpe and para
                    //uparpendicular_puppi= ((-puppimet_px-photon_px)*photon_py-(-puppimet_py-photon_py)*photon_px)/tp.pt
                    //uparallal_puppi= ((-puppimet_px-photon_px)*photon_px+(-puppimet_py-photon_py)*photon_py)/tp.pt

                    puppi_pen = ((-puppiMET_px - puppipho_px) * puppipho_py - (-puppiMET_py - puppipho_py) * puppipho_px)/((*phoEt)[ipuppiPho]);
                    puppi_par = ((-puppiMET_px - puppipho_px) * puppipho_px + (-puppiMET_py - puppipho_py) * puppipho_py)/((*phoEt)[ipuppiPho]);

                    uperpen_puppi->Fill(puppi_pen);
                    uparall_puppi->Fill(puppi_par);
                
                    puppi_abs_scale->Fill(-(puppi_par)/((*phoEt)[ipuppiPho]));

                    puppi_paraqt->Fill(puppi_par + (*phoEt)[ipuppiPho]);


                    puppimet_nGoodvtx->Fill(nGoodVtx);
                
                }

                
                Float_t puppipho_px_bins = 0;
                Float_t puppipho_py_bins = 0;

                Float_t puppi_pen_bins = 0;
                Float_t puppi_par_bins = 0;

                //puppimet with different bins
                for (int ipuppiPho = 0; ipuppiPho < nPho; ipuppiPho++)
                {
                    if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 70)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale5070->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara5070->Fill(puppi_par_bins);
                        puppi_uperpen5070->Fill(puppi_pen_bins);
                        puppi_uparaqt5070->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 70 && (*phoEt)[ipuppiPho] < 90)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale7090->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara7090->Fill(puppi_par_bins);
                        puppi_uperpen7090->Fill(puppi_pen_bins);
                        puppi_uparaqt7090->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 110)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale90110->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara90110->Fill(puppi_par_bins);
                        puppi_uperpen90110->Fill(puppi_pen_bins);
                        puppi_uparaqt90110->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 110 && (*phoEt)[ipuppiPho] < 130)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale110130->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara110130->Fill(puppi_par_bins);
                        puppi_uperpen110130->Fill(puppi_pen_bins);
                        puppi_uparaqt110130->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 150)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale130150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara130150->Fill(puppi_par_bins);
                        puppi_uperpen130150->Fill(puppi_pen_bins);
                        puppi_uparaqt130150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 150 && (*phoEt)[ipuppiPho] < 170)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale150170->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara150170->Fill(puppi_par_bins);
                        puppi_uperpen150170->Fill(puppi_pen_bins);
                        puppi_uparaqt150170->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 190)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale170190->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara170190->Fill(puppi_par_bins);
                        puppi_uperpen170190->Fill(puppi_pen_bins);
                        puppi_uparaqt170190->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 190 && (*phoEt)[ipuppiPho] < 210)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale190210->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara190210->Fill(puppi_par_bins);
                        puppi_uperpen190210->Fill(puppi_pen_bins);
                        puppi_uparaqt190210->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 230)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale210230->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara210230->Fill(puppi_par_bins);
                        puppi_uperpen210230->Fill(puppi_pen_bins);
                        puppi_uparaqt210230->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 230 && (*phoEt)[ipuppiPho] < 250)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale230250->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara230250->Fill(puppi_par_bins);
                        puppi_uperpen230250->Fill(puppi_pen_bins);
                        puppi_uparaqt230250->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 270)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale250270->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara250270->Fill(puppi_par_bins);
                        puppi_uperpen250270->Fill(puppi_pen_bins);
                        puppi_uparaqt250270->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 270 && (*phoEt)[ipuppiPho] < 290)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale270290->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara270290->Fill(puppi_par_bins);
                        puppi_uperpen270290->Fill(puppi_pen_bins);
                        puppi_uparaqt270290->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 310)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale290310->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara290310->Fill(puppi_par_bins);
                        puppi_uperpen290310->Fill(puppi_pen_bins);
                        puppi_uparaqt290310->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 310 && (*phoEt)[ipuppiPho] < 330)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale310330->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara310330->Fill(puppi_par_bins);
                        puppi_uperpen310330->Fill(puppi_pen_bins);
                        puppi_uparaqt310330->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 350)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale330350->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara330350->Fill(puppi_par_bins);
                        puppi_uperpen330350->Fill(puppi_pen_bins);
                        puppi_uparaqt330350->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 350 && (*phoEt)[ipuppiPho] < 370)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale350370->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara350370->Fill(puppi_par_bins);
                        puppi_uperpen350370->Fill(puppi_pen_bins);
                        puppi_uparaqt350370->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 400)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale370400->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara370400->Fill(puppi_par_bins);
                        puppi_uperpen370400->Fill(puppi_pen_bins);
                        puppi_uparaqt370400->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }


                    //vs nvtx
                    if (nGoodVtx > 0 && nGoodVtx < 4)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx04->Fill(puppi_par_bins);
                        puppi_uperpennvtx04->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 4 && nGoodVtx < 6)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx46->Fill(puppi_par_bins);
                        puppi_uperpennvtx46->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 6 && nGoodVtx < 8)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx68->Fill(puppi_par_bins);
                        puppi_uperpennvtx68->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 8 && nGoodVtx < 12)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx812->Fill(puppi_par_bins);
                        puppi_uperpennvtx812->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 12 && nGoodVtx < 18)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx1218->Fill(puppi_par_bins);
                        puppi_uperpennvtx1218->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 18 && nGoodVtx < 24)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx1824->Fill(puppi_par_bins);
                        puppi_uperpennvtx1824->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 24 && nGoodVtx < 28)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx2428->Fill(puppi_par_bins);
                        puppi_uperpennvtx2428->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 28 && nGoodVtx < 32)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx2832->Fill(puppi_par_bins);
                        puppi_uperpennvtx2832->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 32 && nGoodVtx < 36)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx3236->Fill(puppi_par_bins);
                        puppi_uperpennvtx3236->Fill(puppi_pen_bins);
                    }

                    if (nGoodVtx > 36 && nGoodVtx < 40)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        
                        puppi_uparanvtx3640->Fill(puppi_par_bins);
                        puppi_uperpennvtx3640->Fill(puppi_pen_bins);
                    }


                    /*
                    ======--------======--======--======---------=======---------=======------===
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--------======--======--======---------=======---------=========--=====
                    ======--============--======--======--==============--================--=====
                    ======--============----------======--==============--==============------===
                    */

                    //vs puppipt < 100
                    if (puppiMET_pt > 10 && puppiMET_pt < 20)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt1020->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt1020->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt1020->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt1020->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }


                    if (puppiMET_pt > 20 && puppiMET_pt < 30)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt2030->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt2030->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt2030->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt2030->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 30 && puppiMET_pt < 40)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt3040->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt3040->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt3040->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt3040->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 40 && puppiMET_pt < 50)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt4050->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt4050->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt4050->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt4050->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 50 && puppiMET_pt < 60)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt5070->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt5070->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt5070->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt5070->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 60 && puppiMET_pt < 70)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt6070->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt6070->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt6070->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt6070->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 70 && puppiMET_pt < 80)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt7080->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt7080->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt7080->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt7080->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 80 && puppiMET_pt < 90)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt7090->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt7090->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt7090->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt7090->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 90 && puppiMET_pt < 100)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt90100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt90100->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt90100->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt90100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5070_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5070_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen5070_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt5070_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 70 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale7090_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara7090_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen7090_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt7090_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50110_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90110_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90110_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen90110_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt90110_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 110 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale110130_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110130_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110130_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110130_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130150_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130150_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen130150_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt130150_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 150 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale150170_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara150170_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen150170_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt150170_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170190_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170190_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen170190_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt170190_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 190 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale190210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara190210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen190210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt190210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210230_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210230_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210230_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210230_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 230 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale230250_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara230250_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen230250_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt230250_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250270_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250270_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen250270_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt250270_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 270 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale270290_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara270290_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen270290_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt270290_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen290310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt290310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 310 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale310330_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310330_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen310330_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt310330_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330350_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330350_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen330350_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt330350_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 350 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale350370_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara350370_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen350370_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt350370_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370400_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370400_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen370400_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt370400_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx04_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx04_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx46_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx46_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx68_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx68_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx812_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx812_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1218_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx1218_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1824_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx1824_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2428_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx2428_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2832_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx2832_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3236_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx3236_puppi100->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3640_puppi100->Fill(puppi_par_bins);
                            puppi_uperpennvtx3640_puppi100->Fill(puppi_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5070_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5070_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen5070_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt5070_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 70 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale7090_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara7090_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen7090_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt7090_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50110_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90110_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90110_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen90110_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt90110_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 110 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale110130_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110130_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110130_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110130_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130150_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130150_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen130150_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt130150_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 150 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale150170_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara150170_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen150170_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt150170_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170190_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170190_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen170190_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt170190_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 190 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale190210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara190210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen190210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt190210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210230_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210230_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210230_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210230_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 230 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale230250_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara230250_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen230250_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt230250_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250270_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250270_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen250270_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt250270_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 270 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale270290_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara270290_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen270290_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt270290_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen290310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt290310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 310 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale310330_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310330_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen310330_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt310330_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330350_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330350_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen330350_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt330350_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 350 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale350370_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara350370_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen350370_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt350370_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370400_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370400_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen370400_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt370400_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx04_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx04_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx46_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx46_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx68_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx68_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx812_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx812_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1218_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx1218_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1824_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx1824_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2428_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx2428_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2832_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx2832_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3236_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx3236_puppi150->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3640_puppi150->Fill(puppi_par_bins);
                            puppi_uperpennvtx3640_puppi150->Fill(puppi_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5070_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5070_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen5070_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt5070_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 70 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale7090_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara7090_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen7090_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt7090_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50110_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90110_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90110_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen90110_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt90110_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50110_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50110_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen50110_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt50110_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 110 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale110130_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110130_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110130_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110130_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130150_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130150_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen130150_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt130150_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 150 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale150170_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara150170_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen150170_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt150170_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170190_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170190_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen170190_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt170190_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 190 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale190210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara190210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen190210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt190210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale110210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara110210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen110210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt110210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210230_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210230_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210230_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210230_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 230 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale230250_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara230250_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen230250_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt230250_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250270_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250270_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen250270_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt250270_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 270 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale270290_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara270290_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen270290_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt270290_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen290310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt290310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210310_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210310_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210310_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210310_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 310 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale310330_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310330_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen310330_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt310330_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330350_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330350_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen330350_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt330350_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 350 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale350370_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara350370_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen350370_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt350370_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370400_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370400_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen370400_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt370400_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale310400_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara310400_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen310400_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt310400_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx04_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx04_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx46_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx46_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx68_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx68_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx812_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx812_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1218_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx1218_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx1824_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx1824_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2428_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx2428_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx2832_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx2832_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3236_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx3236_puppi200->Fill(puppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            
                            puppi_uparanvtx3640_puppi200->Fill(puppi_par_bins);
                            puppi_uperpennvtx3640_puppi200->Fill(puppi_pen_bins);
                        }
                    }



                    


                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5090_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5090_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen5090_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt5090_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50130_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90130_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90130_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen90130_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt90130_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50130_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130170_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130170_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen130170_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt130170_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen170210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt170210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210250_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210250_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210250_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210250_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale210290_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250290_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250290_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen250290_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt250290_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210290_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290330_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290330_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen290330_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt290330_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330370_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330370_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen330370_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt330370_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370410_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370450_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370450_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen370450_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt370450_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 410 && (*phoEt)[ipuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale410450_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara410450_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen410450_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt410450_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370410_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 450 && (*phoEt)[ipuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale450490_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450490_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen450490_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt450490_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 490 && (*phoEt)[ipuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale490530_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara490530_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen490530_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt490530_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi100->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi100->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi100->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi100->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5090_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5090_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen5090_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt5090_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50130_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90130_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90130_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen90130_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt90130_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50130_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130170_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130170_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen130170_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt130170_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen170210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt170210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210250_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210250_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210250_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210250_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale210290_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250290_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250290_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen250290_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt250290_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210290_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290330_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290330_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen290330_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt290330_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330370_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330370_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen330370_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt330370_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370410_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370450_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370450_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen370450_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt370450_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 410 && (*phoEt)[ipuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale410450_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara410450_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen410450_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt410450_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370410_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 450 && (*phoEt)[ipuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale450490_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450490_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen450490_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt450490_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 490 && (*phoEt)[ipuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale490530_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara490530_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen490530_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt490530_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi150->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi150->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[ipuppiPho] > 50 && (*phoEt)[ipuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale5090_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara5090_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen5090_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt5090_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale50130_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 90 && (*phoEt)[ipuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale90130_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara90130_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen90130_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt90130_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale50130_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara50130_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen50130_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt50130_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 130 && (*phoEt)[ipuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale130170_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130170_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen130170_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt130170_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 170 && (*phoEt)[ipuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale170210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara170210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen170210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt170210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale130210_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara130210_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen130210_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt130210_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 210 && (*phoEt)[ipuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale210250_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210250_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210250_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210250_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                            puppi_abs_scale210290_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);


                        }

                        if ((*phoEt)[ipuppiPho] > 250 && (*phoEt)[ipuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale250290_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara250290_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen250290_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt250290_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale210290_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara210290_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen210290_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt210290_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 290 && (*phoEt)[ipuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale290330_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290330_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen290330_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt290330_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 330 && (*phoEt)[ipuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale330370_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara330370_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen330370_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt330370_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale290370_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara290370_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen290370_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt290370_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 370 && (*phoEt)[ipuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale370410_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370450_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370450_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen370450_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt370450_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 410 && (*phoEt)[ipuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale410450_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara410450_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen410450_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt410450_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale370410_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara370410_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen370410_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt370410_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 450 && (*phoEt)[ipuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale450490_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450490_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen450490_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt450490_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }

                        if ((*phoEt)[ipuppiPho] > 490 && (*phoEt)[ipuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                            puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                            puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                            puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                            
                            puppi_abs_scale490530_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara490530_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen490530_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt490530_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);

                            puppi_abs_scale450530_puppi200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                            puppi_upara450530_puppi200->Fill(puppi_par_bins);
                            puppi_uperpen450530_puppi200->Fill(puppi_pen_bins);
                            puppi_uparaqt450530_puppi200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                        }
                    }





                }

                /*
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                ===========================================================
                */


                //rawpuppimet
                for (int irawpuppiPho = 0; irawpuppiPho < nPho; irawpuppiPho++)
                {
                    rawpuppimetphi->Fill(rawPuppiMETPhi);

                    rawpuppimet_tightphopt->Fill((*phoEt)[irawpuppiPho]);

                    //fill the px py from ggNtuplizer
                    rawpuppimetpt->Fill(rawPuppiMET_pt);
                    rawPuppimet_px->Fill(rawPuppiMET_px);
                    rawPuppimet_py->Fill(rawPuppiMET_py);

                    
                    //try to calculate and also fill the one from ggNtuplizer
                    //rawPuppimet_px = rawpuppimet.pt*m.cos(rawpuppimet.phi)
                    rawPuppimet_px_cal->Fill((rawPuppiMET_pt)*cos(rawPuppiMETPhi));
                    rawPuppimet_py_cal->Fill((rawPuppiMET_pt)*sin(rawPuppiMETPhi));
                    //==> passed test! the same as ggNtuplizer stored.
                    

                    //cal photon_px and photon_py
                    rawpuppipho_px = rawpuppipho_px + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                    rawpuppipho_py = rawpuppipho_py + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                    rawpuppipho_pxdist->Fill(rawpuppipho_px);
                    rawpuppipho_pydist->Fill(rawpuppipho_py);
                    

                    //calculate the perpe and para
                    //uparpendicular_rawpuppi= ((-rawPuppimet_px-photon_px)*photon_py-(-rawPuppimet_py-photon_py)*photon_px)/tp.pt
                    //uparallal_rawpuppi= ((-rawPuppimet_px-photon_px)*photon_px+(-rawPuppimet_py-photon_py)*photon_py)/tp.pt

                    rawpuppi_pen = ((-rawPuppiMET_px - rawpuppipho_px) * rawpuppipho_py - (-rawPuppiMET_py - rawpuppipho_py) * rawpuppipho_px)/((*phoEt)[irawpuppiPho]);
                    rawpuppi_par = ((-rawPuppiMET_px - rawpuppipho_px) * rawpuppipho_px + (-rawPuppiMET_py - rawpuppipho_py) * rawpuppipho_py)/((*phoEt)[irawpuppiPho]);

                    uperpen_rawpuppi->Fill(rawpuppi_pen);
                    uparall_rawpuppi->Fill(rawpuppi_par);
                
                    rawpuppi_abs_scale->Fill(-(rawpuppi_par)/((*phoEt)[irawpuppiPho]));

                    rawpuppi_paraqt->Fill(rawpuppi_par + (*phoEt)[irawpuppiPho]);


                    rawpuppimet_nGoodvtx->Fill(nGoodVtx);
                
                }


                Float_t rawpuppipho_px_bins = 0;
                Float_t rawpuppipho_py_bins = 0;

                Float_t rawpuppi_pen_bins = 0;
                Float_t rawpuppi_par_bins = 0;

                //rawpuppimet with different bins
                for (int irawpuppiPho = 0; irawpuppiPho < nPho; irawpuppiPho++)
                {
                    if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 70)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale5070->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara5070->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen5070->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt5070->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 70 && (*phoEt)[irawpuppiPho] < 90)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale7090->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara7090->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen7090->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt7090->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 110)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale90110->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara90110->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen90110->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt90110->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 110 && (*phoEt)[irawpuppiPho] < 130)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale110130->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara110130->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen110130->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt110130->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 150)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale130150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara130150->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen130150->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt130150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 150 && (*phoEt)[irawpuppiPho] < 170)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale150170->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara150170->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen150170->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt150170->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 190)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale170190->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara170190->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen170190->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt170190->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 190 && (*phoEt)[irawpuppiPho] < 210)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale190210->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara190210->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen190210->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt190210->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 230)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale210230->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara210230->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen210230->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt210230->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 230 && (*phoEt)[irawpuppiPho] < 250)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale230250->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara230250->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen230250->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt230250->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 270)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale250270->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara250270->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen250270->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt250270->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 270 && (*phoEt)[irawpuppiPho] < 290)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale270290->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara270290->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen270290->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt270290->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 310)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale290310->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara290310->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen290310->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt290310->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 310 && (*phoEt)[irawpuppiPho] < 330)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale310330->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara310330->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen310330->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt310330->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 350)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale330350->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara330350->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen330350->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt330350->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 350 && (*phoEt)[irawpuppiPho] < 370)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale350370->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara350370->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen350370->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt350370->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 400)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale370400->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara370400->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen370400->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt370400->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }


                    //vs nvtx
                    if (nGoodVtx > 0 && nGoodVtx < 4)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx04->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx04->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 4 && nGoodVtx < 6)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx46->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx46->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 6 && nGoodVtx < 8)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx68->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx68->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 8 && nGoodVtx < 12)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx812->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx812->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 12 && nGoodVtx < 18)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx1218->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx1218->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 18 && nGoodVtx < 24)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx1824->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx1824->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 24 && nGoodVtx < 28)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx2428->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx2428->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 28 && nGoodVtx < 32)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx2832->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx2832->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 32 && nGoodVtx < 36)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx3236->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx3236->Fill(rawpuppi_pen_bins);
                    }

                    if (nGoodVtx > 36 && nGoodVtx < 40)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        
                        rawpuppi_uparanvtx3640->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpennvtx3640->Fill(rawpuppi_pen_bins);
                    }



                    /*
                    ======--------======--======--======---------=======---------=======------===
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--------======--======--======---------=======---------=========--=====
                    ======--============--======--======--==============--================--=====
                    ======--============----------======--==============--==============------===
                    */

                    //vs puppipt < 100
                    if (puppiMET_pt > 10 && puppiMET_pt < 20)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt1020->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt1020->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt1020->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt1020->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }


                    if (puppiMET_pt > 20 && puppiMET_pt < 30)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt2030->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt2030->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt2030->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt2030->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 30 && puppiMET_pt < 40)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt3040->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt3040->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt3040->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt3040->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 40 && puppiMET_pt < 50)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt4050->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt4050->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt4050->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt4050->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 50 && puppiMET_pt < 60)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt5070->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt5070->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt5070->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt5070->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 60 && puppiMET_pt < 70)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt6070->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt6070->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt6070->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt6070->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 70 && puppiMET_pt < 80)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt7080->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt7080->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt7080->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt7080->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 80 && puppiMET_pt < 90)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt7090->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt7090->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt7090->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt7090->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (puppiMET_pt > 90 && puppiMET_pt < 100)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt90100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt90100->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt90100->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt90100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }


                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5070_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5070_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5070_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5070_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 70 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale7090_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara7090_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen7090_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt7090_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50110_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90110_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90110_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90110_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90110_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 110 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale110130_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110130_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110130_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110130_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130150_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130150_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130150_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130150_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 150 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale150170_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara150170_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen150170_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt150170_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170190_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170190_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170190_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170190_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 190 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale190210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara190210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen190210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt190210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210230_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210230_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210230_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210230_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 230 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale230250_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara230250_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen230250_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt230250_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250270_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250270_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250270_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250270_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 270 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale270290_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara270290_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen270290_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt270290_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 310 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale310330_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310330_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310330_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310330_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330350_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330350_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330350_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330350_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 350 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale350370_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara350370_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen350370_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt350370_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370400_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370400_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370400_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370400_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx04_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx04_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx46_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx46_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx68_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx68_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx812_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx812_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1218_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1218_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1824_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1824_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2428_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2428_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2832_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2832_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3236_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3236_puppi100->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3640_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3640_puppi100->Fill(rawpuppi_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5070_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5070_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5070_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5070_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 70 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale7090_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara7090_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen7090_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt7090_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50110_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90110_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90110_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90110_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90110_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 110 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale110130_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110130_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110130_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110130_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130150_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130150_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130150_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130150_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 150 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale150170_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara150170_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen150170_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt150170_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170190_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170190_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170190_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170190_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 190 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale190210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara190210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen190210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt190210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210230_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210230_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210230_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210230_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 230 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale230250_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara230250_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen230250_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt230250_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250270_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250270_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250270_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250270_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 270 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale270290_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara270290_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen270290_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt270290_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 310 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale310330_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310330_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310330_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310330_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330350_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330350_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330350_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330350_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 350 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale350370_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara350370_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen350370_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt350370_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370400_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370400_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370400_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370400_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx04_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx04_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx46_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx46_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx68_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx68_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx812_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx812_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1218_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1218_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1824_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1824_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2428_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2428_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2832_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2832_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3236_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3236_puppi150->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3640_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3640_puppi150->Fill(rawpuppi_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5070_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5070_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5070_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5070_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 70 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale7090_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara7090_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen7090_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt7090_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50110_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90110_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90110_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90110_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90110_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50110_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50110_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50110_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50110_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 110 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale110130_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110130_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110130_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110130_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130150_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130150_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130150_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130150_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 150 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale150170_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara150170_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen150170_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt150170_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170190_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170190_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170190_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170190_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 190 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale190210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara190210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen190210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt190210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale110210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara110210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen110210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt110210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210230_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210230_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210230_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210230_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 230 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale230250_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara230250_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen230250_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt230250_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250270_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250270_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250270_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250270_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 270 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale270290_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara270290_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen270290_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt270290_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210310_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210310_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210310_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210310_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 310 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale310330_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310330_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310330_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310330_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330350_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330350_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330350_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330350_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 350 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale350370_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara350370_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen350370_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt350370_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370400_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370400_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370400_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370400_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale310400_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara310400_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen310400_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt310400_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx04_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx04_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx46_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx46_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx68_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx68_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx812_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx812_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1218_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1218_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx1824_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx1824_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2428_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2428_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx2832_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx2832_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3236_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3236_puppi200->Fill(rawpuppi_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            
                            rawpuppi_uparanvtx3640_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpennvtx3640_puppi200->Fill(rawpuppi_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5090_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5090_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5090_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5090_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50130_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90130_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90130_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90130_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90130_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50130_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130170_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130170_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130170_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130170_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210250_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210250_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210250_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210250_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale210290_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250290_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250290_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250290_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250290_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210290_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290330_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290330_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290330_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290330_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330370_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330370_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330370_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330370_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370410_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370450_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370450_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370450_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370450_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 410 && (*phoEt)[irawpuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale410450_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara410450_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen410450_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt410450_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370410_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 450 && (*phoEt)[irawpuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale450490_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450490_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450490_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450490_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 490 && (*phoEt)[irawpuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale490530_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara490530_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen490530_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt490530_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi100->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi100->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi100->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi100->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5090_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5090_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5090_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5090_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50130_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90130_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90130_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90130_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90130_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50130_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130170_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130170_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130170_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130170_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210250_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210250_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210250_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210250_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale210290_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250290_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250290_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250290_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250290_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210290_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290330_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290330_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290330_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290330_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330370_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330370_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330370_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330370_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370410_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370450_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370450_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370450_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370450_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 410 && (*phoEt)[irawpuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale410450_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara410450_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen410450_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt410450_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370410_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 450 && (*phoEt)[irawpuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale450490_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450490_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450490_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450490_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 490 && (*phoEt)[irawpuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale490530_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara490530_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen490530_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt490530_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi150->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi150->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[irawpuppiPho] > 50 && (*phoEt)[irawpuppiPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale5090_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara5090_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen5090_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt5090_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale50130_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 90 && (*phoEt)[irawpuppiPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale90130_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara90130_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen90130_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt90130_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale50130_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara50130_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen50130_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt50130_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 130 && (*phoEt)[irawpuppiPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale130170_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130170_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130170_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130170_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 170 && (*phoEt)[irawpuppiPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale170210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara170210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen170210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt170210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale130210_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara130210_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen130210_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt130210_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 210 && (*phoEt)[irawpuppiPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale210250_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210250_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210250_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210250_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                            rawpuppi_abs_scale210290_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);


                        }

                        if ((*phoEt)[irawpuppiPho] > 250 && (*phoEt)[irawpuppiPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale250290_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara250290_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen250290_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt250290_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale210290_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara210290_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen210290_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt210290_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 290 && (*phoEt)[irawpuppiPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale290330_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290330_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290330_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290330_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 330 && (*phoEt)[irawpuppiPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale330370_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara330370_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen330370_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt330370_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale290370_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara290370_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen290370_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt290370_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 370 && (*phoEt)[irawpuppiPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale370410_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370450_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370450_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370450_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370450_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 410 && (*phoEt)[irawpuppiPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale410450_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara410450_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen410450_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt410450_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale370410_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara370410_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen370410_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt370410_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 450 && (*phoEt)[irawpuppiPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale450490_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450490_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450490_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450490_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }

                        if ((*phoEt)[irawpuppiPho] > 490 && (*phoEt)[irawpuppiPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                            rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                            rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                            rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                            
                            rawpuppi_abs_scale490530_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara490530_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen490530_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt490530_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);

                            rawpuppi_abs_scale450530_puppi200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                            rawpuppi_upara450530_puppi200->Fill(rawpuppi_par_bins);
                            rawpuppi_uperpen450530_puppi200->Fill(rawpuppi_pen_bins);
                            rawpuppi_uparaqt450530_puppi200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                        }
                    }

                }
		
               
                //pfmet
                for (int ipfPho = 0; ipfPho < nPho; ipfPho++)
                {
                    pfmetphi->Fill(pfMETPhi);

                    pfmet_tightphopt->Fill((*phoEt)[ipfPho]);

                    //fill the px py from ggNtuplizer
                    pfmetpt->Fill(pfMET_pt);
                    pfmet_px->Fill(pfMET_px);
                    pfmet_py->Fill(pfMET_py);


                    
                    //try to calculate and also fill the one from ggNtuplizer
                    //pfmet_px = pfmet.pt*m.cos(pfmet.phi)
                    pfmet_px_cal->Fill((pfMET_pt)*cos(pfMETPhi));
                    pfmet_py_cal->Fill((pfMET_pt)*sin(pfMETPhi));
                    //==> passed test! the same as ggNtuplizer stored.
                    

                    //cal photon_px and photon_py
                    pfpho_px = pfpho_px + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                    pfpho_py = pfpho_py + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                    pfpho_pxdist->Fill(pfpho_px);
                    pfpho_pydist->Fill(pfpho_py);

                    //calculate the perpe and para
                    //uparpendicular_pf= ((-pfmet_px-photon_px)*photon_py-(-pfmet_py-photon_py)*photon_px)/tp.pt
                    //uparallal_pf= ((-pfmet_px-photon_px)*photon_px+(-pfmet_py-photon_py)*photon_py)/tp.pt

                    pf_pen = ((-pfMET_px - pfpho_px) * pfpho_py - (-pfMET_py - pfpho_py) * pfpho_px)/((*phoEt)[ipfPho]);
                    pf_par = ((-pfMET_px - pfpho_px) * pfpho_px + (-pfMET_py - pfpho_py) * pfpho_py)/((*phoEt)[ipfPho]);
                    
                    uperpen_pf->Fill(pf_pen);
                    uparall_pf->Fill(pf_par);

                    pf_abs_scale->Fill(-(pf_par)/((*phoEt)[ipfPho]));

                    pf_paraqt->Fill(pf_par + (*phoEt)[ipfPho]);


                    pf_nGoodvtx->Fill(nGoodVtx);
                }


                Float_t pfpho_px_bins = 0;
                Float_t pfpho_py_bins = 0;

                Float_t pf_pen_bins = 0;
                Float_t pf_par_bins = 0;

                //pfmet with different bins
                for (int ipfPho = 0; ipfPho < nPho; ipfPho++)
                {
                    if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 70)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale5070->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara5070->Fill(pf_par_bins);
                        pf_uperpen5070->Fill(pf_pen_bins);
                        pf_uparaqt5070->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 70 && (*phoEt)[ipfPho] < 90)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale7090->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara7090->Fill(pf_par_bins);
                        pf_uperpen7090->Fill(pf_pen_bins);
                        pf_uparaqt7090->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 110)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale90110->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara90110->Fill(pf_par_bins);
                        pf_uperpen90110->Fill(pf_pen_bins);
                        pf_uparaqt90110->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 110 && (*phoEt)[ipfPho] < 130)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale110130->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara110130->Fill(pf_par_bins);
                        pf_uperpen110130->Fill(pf_pen_bins);
                        pf_uparaqt110130->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 150)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale130150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara130150->Fill(pf_par_bins);
                        pf_uperpen130150->Fill(pf_pen_bins);
                        pf_uparaqt130150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 150 && (*phoEt)[ipfPho] < 170)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale150170->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara150170->Fill(pf_par_bins);
                        pf_uperpen150170->Fill(pf_pen_bins);
                        pf_uparaqt150170->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 190)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale170190->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara170190->Fill(pf_par_bins);
                        pf_uperpen170190->Fill(pf_pen_bins);
                        pf_uparaqt170190->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 190 && (*phoEt)[ipfPho] < 210)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale190210->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara190210->Fill(pf_par_bins);
                        pf_uperpen190210->Fill(pf_pen_bins);
                        pf_uparaqt190210->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 230)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale210230->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara210230->Fill(pf_par_bins);
                        pf_uperpen210230->Fill(pf_pen_bins);
                        pf_uparaqt210230->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 230 && (*phoEt)[ipfPho] < 250)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale230250->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara230250->Fill(pf_par_bins);
                        pf_uperpen230250->Fill(pf_pen_bins);
                        pf_uparaqt230250->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 270)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale250270->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara250270->Fill(pf_par_bins);
                        pf_uperpen250270->Fill(pf_pen_bins);
                        pf_uparaqt250270->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 270 && (*phoEt)[ipfPho] < 290)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale270290->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara270290->Fill(pf_par_bins);
                        pf_uperpen270290->Fill(pf_pen_bins);
                        pf_uparaqt270290->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 310)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale290310->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara290310->Fill(pf_par_bins);
                        pf_uperpen290310->Fill(pf_pen_bins);
                        pf_uparaqt290310->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 310 && (*phoEt)[ipfPho] < 330)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale310330->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara310330->Fill(pf_par_bins);
                        pf_uperpen310330->Fill(pf_pen_bins);
                        pf_uparaqt310330->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 350)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale330350->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara330350->Fill(pf_par_bins);
                        pf_uperpen330350->Fill(pf_pen_bins);
                        pf_uparaqt330350->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 350 && (*phoEt)[ipfPho] < 370)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale350370->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara350370->Fill(pf_par_bins);
                        pf_uperpen350370->Fill(pf_pen_bins);
                        pf_uparaqt350370->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 400)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale370400->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara370400->Fill(pf_par_bins);
                        pf_uperpen370400->Fill(pf_pen_bins);
                        pf_uparaqt370400->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }



                    //vs nvtx
                    if (nGoodVtx > 0 && nGoodVtx < 4)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx04->Fill(pf_par_bins);
                        pf_uperpennvtx04->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 4 && nGoodVtx < 6)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx46->Fill(pf_par_bins);
                        pf_uperpennvtx46->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 6 && nGoodVtx < 8)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx68->Fill(pf_par_bins);
                        pf_uperpennvtx68->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 8 && nGoodVtx < 12)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx812->Fill(pf_par_bins);
                        pf_uperpennvtx812->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 12 && nGoodVtx < 18)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx1218->Fill(pf_par_bins);
                        pf_uperpennvtx1218->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 18 && nGoodVtx < 24)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx1824->Fill(pf_par_bins);
                        pf_uperpennvtx1824->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 24 && nGoodVtx < 28)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx2428->Fill(pf_par_bins);
                        pf_uperpennvtx2428->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 28 && nGoodVtx < 32)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx2832->Fill(pf_par_bins);
                        pf_uperpennvtx2832->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 32 && nGoodVtx < 36)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx3236->Fill(pf_par_bins);
                        pf_uperpennvtx3236->Fill(pf_pen_bins);
                    }

                    if (nGoodVtx > 36 && nGoodVtx < 40)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        
                        pf_uparanvtx3640->Fill(pf_par_bins);
                        pf_uperpennvtx3640->Fill(pf_pen_bins);
                    }


                    /*
                    ======--------======--======--======---------=======---------=======------===
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--------======--======--======---------=======---------=========--=====
                    ======--============--======--======--==============--================--=====
                    ======--============----------======--==============--==============------===
                    */

                    //vs puppipt < 100
                    if (puppiMET_pt > 10 && puppiMET_pt < 20)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt1020->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt1020->Fill(pf_par_bins);
                        pf_uperpen_puppipt1020->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt1020->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }


                    if (puppiMET_pt > 20 && puppiMET_pt < 30)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt2030->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt2030->Fill(pf_par_bins);
                        pf_uperpen_puppipt2030->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt2030->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 30 && puppiMET_pt < 40)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt3040->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt3040->Fill(pf_par_bins);
                        pf_uperpen_puppipt3040->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt3040->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 40 && puppiMET_pt < 50)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt4050->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt4050->Fill(pf_par_bins);
                        pf_uperpen_puppipt4050->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt4050->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 50 && puppiMET_pt < 60)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt5070->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt5070->Fill(pf_par_bins);
                        pf_uperpen_puppipt5070->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt5070->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 60 && puppiMET_pt < 70)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt6070->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt6070->Fill(pf_par_bins);
                        pf_uperpen_puppipt6070->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt6070->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 70 && puppiMET_pt < 80)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt7080->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt7080->Fill(pf_par_bins);
                        pf_uperpen_puppipt7080->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt7080->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 80 && puppiMET_pt < 90)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt7090->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt7090->Fill(pf_par_bins);
                        pf_uperpen_puppipt7090->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt7090->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 90 && puppiMET_pt < 100)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt90100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt90100->Fill(pf_par_bins);
                        pf_uperpen_puppipt90100->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt90100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 70)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5070_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5070_puppi100->Fill(pf_par_bins);
                            pf_uperpen5070_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt5070_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi100->Fill(pf_par_bins);
                            pf_uperpen50110_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 70 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale7090_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara7090_puppi100->Fill(pf_par_bins);
                            pf_uperpen7090_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt7090_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50110_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi100->Fill(pf_par_bins);
                            pf_uperpen50110_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 110)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90110_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90110_puppi100->Fill(pf_par_bins);
                            pf_uperpen90110_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt90110_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi100->Fill(pf_par_bins);
                            pf_uperpen50110_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 110 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale110130_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110130_puppi100->Fill(pf_par_bins);
                            pf_uperpen110130_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110130_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi100->Fill(pf_par_bins);
                            pf_uperpen110210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 150)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130150_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130150_puppi100->Fill(pf_par_bins);
                            pf_uperpen130150_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt130150_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi100->Fill(pf_par_bins);
                            pf_uperpen110210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 150 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale150170_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara150170_puppi100->Fill(pf_par_bins);
                            pf_uperpen150170_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt150170_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi100->Fill(pf_par_bins);
                            pf_uperpen110210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 190)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170190_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170190_puppi100->Fill(pf_par_bins);
                            pf_uperpen170190_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt170190_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi100->Fill(pf_par_bins);
                            pf_uperpen110210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 190 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale190210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara190210_puppi100->Fill(pf_par_bins);
                            pf_uperpen190210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt190210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi100->Fill(pf_par_bins);
                            pf_uperpen110210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 230)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210230_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210230_puppi100->Fill(pf_par_bins);
                            pf_uperpen210230_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210230_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi100->Fill(pf_par_bins);
                            pf_uperpen210310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 230 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale230250_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara230250_puppi100->Fill(pf_par_bins);
                            pf_uperpen230250_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt230250_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi100->Fill(pf_par_bins);
                            pf_uperpen210310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 270)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250270_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250270_puppi100->Fill(pf_par_bins);
                            pf_uperpen250270_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt250270_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi100->Fill(pf_par_bins);
                            pf_uperpen210310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 270 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale270290_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara270290_puppi100->Fill(pf_par_bins);
                            pf_uperpen270290_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt270290_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi100->Fill(pf_par_bins);
                            pf_uperpen210310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 310)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290310_puppi100->Fill(pf_par_bins);
                            pf_uperpen290310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt290310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi100->Fill(pf_par_bins);
                            pf_uperpen210310_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 310 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale310330_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310330_puppi100->Fill(pf_par_bins);
                            pf_uperpen310330_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt310330_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi100->Fill(pf_par_bins);
                            pf_uperpen310400_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 350)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330350_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330350_puppi100->Fill(pf_par_bins);
                            pf_uperpen330350_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt330350_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi100->Fill(pf_par_bins);
                            pf_uperpen310400_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 350 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale350370_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara350370_puppi100->Fill(pf_par_bins);
                            pf_uperpen350370_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt350370_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi100->Fill(pf_par_bins);
                            pf_uperpen310400_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 400)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370400_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370400_puppi100->Fill(pf_par_bins);
                            pf_uperpen370400_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt370400_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi100->Fill(pf_par_bins);
                            pf_uperpen310400_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx04_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx04_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx46_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx46_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx68_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx68_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx812_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx812_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1218_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx1218_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1824_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx1824_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2428_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx2428_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2832_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx2832_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3236_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx3236_puppi100->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3640_puppi100->Fill(pf_par_bins);
                            pf_uperpennvtx3640_puppi100->Fill(pf_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 70)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5070_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5070_puppi150->Fill(pf_par_bins);
                            pf_uperpen5070_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt5070_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi150->Fill(pf_par_bins);
                            pf_uperpen50110_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 70 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale7090_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara7090_puppi150->Fill(pf_par_bins);
                            pf_uperpen7090_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt7090_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50110_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi150->Fill(pf_par_bins);
                            pf_uperpen50110_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 110)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90110_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90110_puppi150->Fill(pf_par_bins);
                            pf_uperpen90110_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt90110_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi150->Fill(pf_par_bins);
                            pf_uperpen50110_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 110 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale110130_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110130_puppi150->Fill(pf_par_bins);
                            pf_uperpen110130_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110130_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi150->Fill(pf_par_bins);
                            pf_uperpen110210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 150)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130150_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130150_puppi150->Fill(pf_par_bins);
                            pf_uperpen130150_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt130150_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi150->Fill(pf_par_bins);
                            pf_uperpen110210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 150 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale150170_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara150170_puppi150->Fill(pf_par_bins);
                            pf_uperpen150170_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt150170_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi150->Fill(pf_par_bins);
                            pf_uperpen110210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 190)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170190_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170190_puppi150->Fill(pf_par_bins);
                            pf_uperpen170190_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt170190_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi150->Fill(pf_par_bins);
                            pf_uperpen110210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 190 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale190210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara190210_puppi150->Fill(pf_par_bins);
                            pf_uperpen190210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt190210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi150->Fill(pf_par_bins);
                            pf_uperpen110210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 230)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210230_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210230_puppi150->Fill(pf_par_bins);
                            pf_uperpen210230_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210230_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi150->Fill(pf_par_bins);
                            pf_uperpen210310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 230 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale230250_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara230250_puppi150->Fill(pf_par_bins);
                            pf_uperpen230250_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt230250_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi150->Fill(pf_par_bins);
                            pf_uperpen210310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 270)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250270_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250270_puppi150->Fill(pf_par_bins);
                            pf_uperpen250270_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt250270_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi150->Fill(pf_par_bins);
                            pf_uperpen210310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 270 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale270290_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara270290_puppi150->Fill(pf_par_bins);
                            pf_uperpen270290_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt270290_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi150->Fill(pf_par_bins);
                            pf_uperpen210310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 310)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290310_puppi150->Fill(pf_par_bins);
                            pf_uperpen290310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt290310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi150->Fill(pf_par_bins);
                            pf_uperpen210310_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 310 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale310330_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310330_puppi150->Fill(pf_par_bins);
                            pf_uperpen310330_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt310330_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi150->Fill(pf_par_bins);
                            pf_uperpen310400_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 350)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330350_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330350_puppi150->Fill(pf_par_bins);
                            pf_uperpen330350_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt330350_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi150->Fill(pf_par_bins);
                            pf_uperpen310400_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 350 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale350370_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara350370_puppi150->Fill(pf_par_bins);
                            pf_uperpen350370_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt350370_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi150->Fill(pf_par_bins);
                            pf_uperpen310400_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 400)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370400_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370400_puppi150->Fill(pf_par_bins);
                            pf_uperpen370400_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt370400_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi150->Fill(pf_par_bins);
                            pf_uperpen310400_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx04_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx04_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx46_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx46_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx68_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx68_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx812_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx812_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1218_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx1218_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1824_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx1824_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2428_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx2428_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2832_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx2832_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3236_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx3236_puppi150->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3640_puppi150->Fill(pf_par_bins);
                            pf_uperpennvtx3640_puppi150->Fill(pf_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 70)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5070_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5070_puppi200->Fill(pf_par_bins);
                            pf_uperpen5070_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt5070_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi200->Fill(pf_par_bins);
                            pf_uperpen50110_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 70 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale7090_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara7090_puppi200->Fill(pf_par_bins);
                            pf_uperpen7090_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt7090_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50110_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi200->Fill(pf_par_bins);
                            pf_uperpen50110_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 110)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90110_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90110_puppi200->Fill(pf_par_bins);
                            pf_uperpen90110_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt90110_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50110_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50110_puppi200->Fill(pf_par_bins);
                            pf_uperpen50110_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt50110_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 110 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale110130_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110130_puppi200->Fill(pf_par_bins);
                            pf_uperpen110130_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110130_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi200->Fill(pf_par_bins);
                            pf_uperpen110210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 150)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130150_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130150_puppi200->Fill(pf_par_bins);
                            pf_uperpen130150_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt130150_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi200->Fill(pf_par_bins);
                            pf_uperpen110210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 150 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale150170_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara150170_puppi200->Fill(pf_par_bins);
                            pf_uperpen150170_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt150170_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi200->Fill(pf_par_bins);
                            pf_uperpen110210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 190)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170190_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170190_puppi200->Fill(pf_par_bins);
                            pf_uperpen170190_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt170190_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi200->Fill(pf_par_bins);
                            pf_uperpen110210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 190 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale190210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara190210_puppi200->Fill(pf_par_bins);
                            pf_uperpen190210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt190210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale110210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara110210_puppi200->Fill(pf_par_bins);
                            pf_uperpen110210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt110210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 230)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210230_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210230_puppi200->Fill(pf_par_bins);
                            pf_uperpen210230_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210230_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi200->Fill(pf_par_bins);
                            pf_uperpen210310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 230 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale230250_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara230250_puppi200->Fill(pf_par_bins);
                            pf_uperpen230250_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt230250_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi200->Fill(pf_par_bins);
                            pf_uperpen210310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 270)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250270_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250270_puppi200->Fill(pf_par_bins);
                            pf_uperpen250270_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt250270_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi200->Fill(pf_par_bins);
                            pf_uperpen210310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 270 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale270290_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara270290_puppi200->Fill(pf_par_bins);
                            pf_uperpen270290_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt270290_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi200->Fill(pf_par_bins);
                            pf_uperpen210310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 310)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290310_puppi200->Fill(pf_par_bins);
                            pf_uperpen290310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt290310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210310_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210310_puppi200->Fill(pf_par_bins);
                            pf_uperpen210310_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210310_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 310 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale310330_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310330_puppi200->Fill(pf_par_bins);
                            pf_uperpen310330_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt310330_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi200->Fill(pf_par_bins);
                            pf_uperpen310400_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 350)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330350_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330350_puppi200->Fill(pf_par_bins);
                            pf_uperpen330350_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt330350_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi200->Fill(pf_par_bins);
                            pf_uperpen310400_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 350 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale350370_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara350370_puppi200->Fill(pf_par_bins);
                            pf_uperpen350370_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt350370_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi200->Fill(pf_par_bins);
                            pf_uperpen310400_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 400)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370400_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370400_puppi200->Fill(pf_par_bins);
                            pf_uperpen370400_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt370400_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale310400_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara310400_puppi200->Fill(pf_par_bins);
                            pf_uperpen310400_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt310400_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx04_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx04_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx46_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx46_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx68_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx68_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx812_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx812_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1218_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx1218_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx1824_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx1824_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2428_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx2428_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx2832_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx2832_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3236_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx3236_puppi200->Fill(pf_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            
                            pf_uparanvtx3640_puppi200->Fill(pf_par_bins);
                            pf_uperpennvtx3640_puppi200->Fill(pf_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5090_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5090_puppi100->Fill(pf_par_bins);
                            pf_uperpen5090_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt5090_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50130_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi100->Fill(pf_par_bins);
                            pf_uperpen50130_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90130_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90130_puppi100->Fill(pf_par_bins);
                            pf_uperpen90130_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt90130_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50130_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi100->Fill(pf_par_bins);
                            pf_uperpen50130_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130170_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130170_puppi100->Fill(pf_par_bins);
                            pf_uperpen130170_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt130170_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi100->Fill(pf_par_bins);
                            pf_uperpen130210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170210_puppi100->Fill(pf_par_bins);
                            pf_uperpen170210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt170210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi100->Fill(pf_par_bins);
                            pf_uperpen130210_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210250_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210250_puppi100->Fill(pf_par_bins);
                            pf_uperpen210250_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210250_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale210290_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi100->Fill(pf_par_bins);
                            pf_uperpen210290_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250290_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250290_puppi100->Fill(pf_par_bins);
                            pf_uperpen250290_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt250290_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210290_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi100->Fill(pf_par_bins);
                            pf_uperpen210290_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290330_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290330_puppi100->Fill(pf_par_bins);
                            pf_uperpen290330_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt290330_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi100->Fill(pf_par_bins);
                            pf_uperpen290370_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330370_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330370_puppi100->Fill(pf_par_bins);
                            pf_uperpen330370_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt330370_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi100->Fill(pf_par_bins);
                            pf_uperpen290370_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 410)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370410_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi100->Fill(pf_par_bins);
                            pf_uperpen370410_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370450_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370450_puppi100->Fill(pf_par_bins);
                            pf_uperpen370450_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt370450_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 410 && (*phoEt)[ipfPho] < 450)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale410450_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara410450_puppi100->Fill(pf_par_bins);
                            pf_uperpen410450_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt410450_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370410_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi100->Fill(pf_par_bins);
                            pf_uperpen370410_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 450 && (*phoEt)[ipfPho] < 490)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale450490_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450490_puppi100->Fill(pf_par_bins);
                            pf_uperpen450490_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt450490_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi100->Fill(pf_par_bins);
                            pf_uperpen450530_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 490 && (*phoEt)[ipfPho] < 530)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale490530_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara490530_puppi100->Fill(pf_par_bins);
                            pf_uperpen490530_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt490530_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi100->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi100->Fill(pf_par_bins);
                            pf_uperpen450530_puppi100->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi100->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5090_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5090_puppi150->Fill(pf_par_bins);
                            pf_uperpen5090_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt5090_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50130_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi150->Fill(pf_par_bins);
                            pf_uperpen50130_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90130_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90130_puppi150->Fill(pf_par_bins);
                            pf_uperpen90130_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt90130_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50130_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi150->Fill(pf_par_bins);
                            pf_uperpen50130_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130170_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130170_puppi150->Fill(pf_par_bins);
                            pf_uperpen130170_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt130170_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi150->Fill(pf_par_bins);
                            pf_uperpen130210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170210_puppi150->Fill(pf_par_bins);
                            pf_uperpen170210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt170210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi150->Fill(pf_par_bins);
                            pf_uperpen130210_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210250_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210250_puppi150->Fill(pf_par_bins);
                            pf_uperpen210250_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210250_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale210290_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi150->Fill(pf_par_bins);
                            pf_uperpen210290_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250290_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250290_puppi150->Fill(pf_par_bins);
                            pf_uperpen250290_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt250290_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210290_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi150->Fill(pf_par_bins);
                            pf_uperpen210290_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290330_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290330_puppi150->Fill(pf_par_bins);
                            pf_uperpen290330_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt290330_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi150->Fill(pf_par_bins);
                            pf_uperpen290370_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330370_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330370_puppi150->Fill(pf_par_bins);
                            pf_uperpen330370_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt330370_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi150->Fill(pf_par_bins);
                            pf_uperpen290370_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 410)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370410_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi150->Fill(pf_par_bins);
                            pf_uperpen370410_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370450_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370450_puppi150->Fill(pf_par_bins);
                            pf_uperpen370450_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt370450_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 410 && (*phoEt)[ipfPho] < 450)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale410450_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara410450_puppi150->Fill(pf_par_bins);
                            pf_uperpen410450_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt410450_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370410_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi150->Fill(pf_par_bins);
                            pf_uperpen370410_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 450 && (*phoEt)[ipfPho] < 490)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale450490_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450490_puppi150->Fill(pf_par_bins);
                            pf_uperpen450490_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt450490_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi150->Fill(pf_par_bins);
                            pf_uperpen450530_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 490 && (*phoEt)[ipfPho] < 530)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale490530_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara490530_puppi150->Fill(pf_par_bins);
                            pf_uperpen490530_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt490530_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi150->Fill(pf_par_bins);
                            pf_uperpen450530_puppi150->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[ipfPho] > 50 && (*phoEt)[ipfPho] < 90)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale5090_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara5090_puppi200->Fill(pf_par_bins);
                            pf_uperpen5090_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt5090_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale50130_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi200->Fill(pf_par_bins);
                            pf_uperpen50130_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 90 && (*phoEt)[ipfPho] < 130)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale90130_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara90130_puppi200->Fill(pf_par_bins);
                            pf_uperpen90130_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt90130_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale50130_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara50130_puppi200->Fill(pf_par_bins);
                            pf_uperpen50130_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt50130_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 130 && (*phoEt)[ipfPho] < 170)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale130170_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130170_puppi200->Fill(pf_par_bins);
                            pf_uperpen130170_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt130170_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi200->Fill(pf_par_bins);
                            pf_uperpen130210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 170 && (*phoEt)[ipfPho] < 210)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale170210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara170210_puppi200->Fill(pf_par_bins);
                            pf_uperpen170210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt170210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale130210_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara130210_puppi200->Fill(pf_par_bins);
                            pf_uperpen130210_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt130210_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 210 && (*phoEt)[ipfPho] < 250)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale210250_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210250_puppi200->Fill(pf_par_bins);
                            pf_uperpen210250_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210250_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                            pf_abs_scale210290_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi200->Fill(pf_par_bins);
                            pf_uperpen210290_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);


                        }

                        if ((*phoEt)[ipfPho] > 250 && (*phoEt)[ipfPho] < 290)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale250290_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara250290_puppi200->Fill(pf_par_bins);
                            pf_uperpen250290_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt250290_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale210290_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara210290_puppi200->Fill(pf_par_bins);
                            pf_uperpen210290_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt210290_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 290 && (*phoEt)[ipfPho] < 330)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale290330_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290330_puppi200->Fill(pf_par_bins);
                            pf_uperpen290330_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt290330_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi200->Fill(pf_par_bins);
                            pf_uperpen290370_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 330 && (*phoEt)[ipfPho] < 370)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale330370_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara330370_puppi200->Fill(pf_par_bins);
                            pf_uperpen330370_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt330370_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale290370_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara290370_puppi200->Fill(pf_par_bins);
                            pf_uperpen290370_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt290370_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 370 && (*phoEt)[ipfPho] < 410)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale370410_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi200->Fill(pf_par_bins);
                            pf_uperpen370410_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370450_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370450_puppi200->Fill(pf_par_bins);
                            pf_uperpen370450_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt370450_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 410 && (*phoEt)[ipfPho] < 450)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale410450_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara410450_puppi200->Fill(pf_par_bins);
                            pf_uperpen410450_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt410450_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale370410_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara370410_puppi200->Fill(pf_par_bins);
                            pf_uperpen370410_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt370410_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 450 && (*phoEt)[ipfPho] < 490)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale450490_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450490_puppi200->Fill(pf_par_bins);
                            pf_uperpen450490_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt450490_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi200->Fill(pf_par_bins);
                            pf_uperpen450530_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }

                        if ((*phoEt)[ipfPho] > 490 && (*phoEt)[ipfPho] < 530)
                        {
                            //cal photon_px and photon_py
                            pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                            pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                            pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                            pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                            
                            pf_abs_scale490530_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara490530_puppi200->Fill(pf_par_bins);
                            pf_uperpen490530_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt490530_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);

                            pf_abs_scale450530_puppi200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                            pf_upara450530_puppi200->Fill(pf_par_bins);
                            pf_uperpen450530_puppi200->Fill(pf_pen_bins);
                            pf_uparaqt450530_puppi200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                        }
                    }
                }


                


                //rawmet with all bins
                for (int irawPho = 0; irawPho < nPho; irawPho++)
                {
                    rawmetphi->Fill(rawMETPhi);

                    rawmet_tightphopt->Fill((*phoEt)[irawPho]);

                    //fill the px py from ggNtuplizer
                    rawmetpt->Fill(rawMET_pt);
                    rawmet_px->Fill(rawMET_px);
                    rawmet_py->Fill(rawMET_py);


                    
                    //try to calculate and also fill the one from ggNtuplizer
                    //rawmet_px = rawmet.pt*m.cos(rawmet.phi)
                    rawmet_px_cal->Fill((rawMET_pt)*cos(rawMETPhi));
                    rawmet_py_cal->Fill((rawMET_pt)*sin(rawMETPhi));
                    //==> passed test! the same as ggNtuplizer stored.
                    

                    //cal photon_px and photon_py
                    rawpho_px = rawpho_px + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                    rawpho_py = rawpho_py + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                    rawpho_pxdist->Fill(rawpho_px);
                    rawpho_pydist->Fill(rawpho_py);

                    //calculate the perpe and para
                    //uparpendicular_raw= ((-rawmet_px-photon_px)*photon_py-(-rawmet_py-photon_py)*photon_px)/tp.pt
                    //uparallal_raw= ((-rawmet_px-photon_px)*photon_px+(-rawmet_py-photon_py)*photon_py)/tp.pt

                    raw_pen = ((-rawMET_px - rawpho_px) * rawpho_py - (-rawMET_py - rawpho_py) * rawpho_px)/((*phoEt)[irawPho]);
                    raw_par = ((-rawMET_px - rawpho_px) * rawpho_px + (-rawMET_py - rawpho_py) * rawpho_py)/((*phoEt)[irawPho]);
                    
                    uperpen_raw->Fill(raw_pen);
                    uparall_raw->Fill(raw_par);

                    raw_abs_scale->Fill(-(raw_par)/((*phoEt)[irawPho]));

                    raw_paraqt->Fill(raw_par + (*phoEt)[irawPho]);


                    raw_nGoodvtx->Fill(nGoodVtx);
                }

                Float_t rawpho_px_bins = 0;
                Float_t rawpho_py_bins = 0;

                Float_t raw_pen_bins = 0;
                Float_t raw_par_bins = 0;

                //rawmet with different bins
                for (int irawPho = 0; irawPho < nPho; irawPho++)
                {
                    if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 70)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale5070->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara5070->Fill(raw_par_bins);
                        raw_uperpen5070->Fill(raw_pen_bins);
                        raw_uparaqt5070->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 70 && (*phoEt)[irawPho] < 90)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale7090->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara7090->Fill(raw_par_bins);
                        raw_uperpen7090->Fill(raw_pen_bins);
                        raw_uparaqt7090->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 110)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale90110->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara90110->Fill(raw_par_bins);
                        raw_uperpen90110->Fill(raw_pen_bins);
                        raw_uparaqt90110->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 110 && (*phoEt)[irawPho] < 130)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale110130->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara110130->Fill(raw_par_bins);
                        raw_uperpen110130->Fill(raw_pen_bins);
                        raw_uparaqt110130->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 150)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale130150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara130150->Fill(raw_par_bins);
                        raw_uperpen130150->Fill(raw_pen_bins);
                        raw_uparaqt130150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 150 && (*phoEt)[irawPho] < 170)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale150170->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara150170->Fill(raw_par_bins);
                        raw_uperpen150170->Fill(raw_pen_bins);
                        raw_uparaqt150170->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 190)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale170190->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara170190->Fill(raw_par_bins);
                        raw_uperpen170190->Fill(raw_pen_bins);
                        raw_uparaqt170190->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 190 && (*phoEt)[irawPho] < 210)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale190210->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara190210->Fill(raw_par_bins);
                        raw_uperpen190210->Fill(raw_pen_bins);
                        raw_uparaqt190210->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 230)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale210230->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara210230->Fill(raw_par_bins);
                        raw_uperpen210230->Fill(raw_pen_bins);
                        raw_uparaqt210230->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 230 && (*phoEt)[irawPho] < 250)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale230250->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara230250->Fill(raw_par_bins);
                        raw_uperpen230250->Fill(raw_pen_bins);
                        raw_uparaqt230250->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 270)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale250270->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara250270->Fill(raw_par_bins);
                        raw_uperpen250270->Fill(raw_pen_bins);
                        raw_uparaqt250270->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 270 && (*phoEt)[irawPho] < 290)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale270290->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara270290->Fill(raw_par_bins);
                        raw_uperpen270290->Fill(raw_pen_bins);
                        raw_uparaqt270290->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 310)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale290310->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara290310->Fill(raw_par_bins);
                        raw_uperpen290310->Fill(raw_pen_bins);
                        raw_uparaqt290310->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 310 && (*phoEt)[irawPho] < 330)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale310330->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara310330->Fill(raw_par_bins);
                        raw_uperpen310330->Fill(raw_pen_bins);
                        raw_uparaqt310330->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 350)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale330350->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara330350->Fill(raw_par_bins);
                        raw_uperpen330350->Fill(raw_pen_bins);
                        raw_uparaqt330350->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 350 && (*phoEt)[irawPho] < 370)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale350370->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara350370->Fill(raw_par_bins);
                        raw_uperpen350370->Fill(raw_pen_bins);
                        raw_uparaqt350370->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 400)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale370400->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara370400->Fill(raw_par_bins);
                        raw_uperpen370400->Fill(raw_pen_bins);
                        raw_uparaqt370400->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }


                    //vs nvtx
                    if (nGoodVtx > 0 && nGoodVtx < 4)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx04->Fill(raw_par_bins);
                        raw_uperpennvtx04->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 4 && nGoodVtx < 6)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx46->Fill(raw_par_bins);
                        raw_uperpennvtx46->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 6 && nGoodVtx < 8)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx68->Fill(raw_par_bins);
                        raw_uperpennvtx68->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 8 && nGoodVtx < 12)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx812->Fill(raw_par_bins);
                        raw_uperpennvtx812->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 12 && nGoodVtx < 18)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx1218->Fill(raw_par_bins);
                        raw_uperpennvtx1218->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 18 && nGoodVtx < 24)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx1824->Fill(raw_par_bins);
                        raw_uperpennvtx1824->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 24 && nGoodVtx < 28)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx2428->Fill(raw_par_bins);
                        raw_uperpennvtx2428->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 28 && nGoodVtx < 32)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx2832->Fill(raw_par_bins);
                        raw_uperpennvtx2832->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 32 && nGoodVtx < 36)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx3236->Fill(raw_par_bins);
                        raw_uperpennvtx3236->Fill(raw_pen_bins);
                    }

                    if (nGoodVtx > 36 && nGoodVtx < 40)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        
                        raw_uparanvtx3640->Fill(raw_par_bins);
                        raw_uperpennvtx3640->Fill(raw_pen_bins);
                    }

                    /*
                    ======--------======--======--======---------=======---------=======------===
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--====--======--======--======--=====--=======--=====--=========--=====
                    ======--------======--======--======---------=======---------=========--=====
                    ======--============--======--======--==============--================--=====
                    ======--============----------======--==============--==============------===
                    */

                    //vs puppipt < 100
                    if (puppiMET_pt > 10 && puppiMET_pt < 20)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt1020->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt1020->Fill(raw_par_bins);
                        raw_uperpen_puppipt1020->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt1020->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }


                    if (puppiMET_pt > 20 && puppiMET_pt < 30)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt2030->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt2030->Fill(raw_par_bins);
                        raw_uperpen_puppipt2030->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt2030->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 30 && puppiMET_pt < 40)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt3040->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt3040->Fill(raw_par_bins);
                        raw_uperpen_puppipt3040->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt3040->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 40 && puppiMET_pt < 50)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt4050->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt4050->Fill(raw_par_bins);
                        raw_uperpen_puppipt4050->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt4050->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 50 && puppiMET_pt < 60)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt5070->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt5070->Fill(raw_par_bins);
                        raw_uperpen_puppipt5070->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt5070->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 60 && puppiMET_pt < 70)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt6070->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt6070->Fill(raw_par_bins);
                        raw_uperpen_puppipt6070->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt6070->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 70 && puppiMET_pt < 80)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt7080->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt7080->Fill(raw_par_bins);
                        raw_uperpen_puppipt7080->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt7080->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 80 && puppiMET_pt < 90)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt7090->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt7090->Fill(raw_par_bins);
                        raw_uperpen_puppipt7090->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt7090->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 90 && puppiMET_pt < 100)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt90100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt90100->Fill(raw_par_bins);
                        raw_uperpen_puppipt90100->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt90100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5070_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5070_puppi100->Fill(raw_par_bins);
                            raw_uperpen5070_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt5070_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi100->Fill(raw_par_bins);
                            raw_uperpen50110_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 70 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale7090_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara7090_puppi100->Fill(raw_par_bins);
                            raw_uperpen7090_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt7090_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50110_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi100->Fill(raw_par_bins);
                            raw_uperpen50110_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90110_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90110_puppi100->Fill(raw_par_bins);
                            raw_uperpen90110_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt90110_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi100->Fill(raw_par_bins);
                            raw_uperpen50110_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 110 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale110130_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110130_puppi100->Fill(raw_par_bins);
                            raw_uperpen110130_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110130_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi100->Fill(raw_par_bins);
                            raw_uperpen110210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130150_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130150_puppi100->Fill(raw_par_bins);
                            raw_uperpen130150_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt130150_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi100->Fill(raw_par_bins);
                            raw_uperpen110210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 150 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale150170_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara150170_puppi100->Fill(raw_par_bins);
                            raw_uperpen150170_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt150170_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi100->Fill(raw_par_bins);
                            raw_uperpen110210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170190_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170190_puppi100->Fill(raw_par_bins);
                            raw_uperpen170190_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt170190_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi100->Fill(raw_par_bins);
                            raw_uperpen110210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 190 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale190210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara190210_puppi100->Fill(raw_par_bins);
                            raw_uperpen190210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt190210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi100->Fill(raw_par_bins);
                            raw_uperpen110210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210230_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210230_puppi100->Fill(raw_par_bins);
                            raw_uperpen210230_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210230_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi100->Fill(raw_par_bins);
                            raw_uperpen210310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 230 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale230250_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara230250_puppi100->Fill(raw_par_bins);
                            raw_uperpen230250_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt230250_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi100->Fill(raw_par_bins);
                            raw_uperpen210310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250270_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250270_puppi100->Fill(raw_par_bins);
                            raw_uperpen250270_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt250270_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi100->Fill(raw_par_bins);
                            raw_uperpen210310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 270 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale270290_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara270290_puppi100->Fill(raw_par_bins);
                            raw_uperpen270290_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt270290_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi100->Fill(raw_par_bins);
                            raw_uperpen210310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290310_puppi100->Fill(raw_par_bins);
                            raw_uperpen290310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt290310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi100->Fill(raw_par_bins);
                            raw_uperpen210310_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 310 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale310330_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310330_puppi100->Fill(raw_par_bins);
                            raw_uperpen310330_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt310330_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi100->Fill(raw_par_bins);
                            raw_uperpen310400_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330350_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330350_puppi100->Fill(raw_par_bins);
                            raw_uperpen330350_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt330350_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi100->Fill(raw_par_bins);
                            raw_uperpen310400_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 350 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale350370_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara350370_puppi100->Fill(raw_par_bins);
                            raw_uperpen350370_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt350370_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi100->Fill(raw_par_bins);
                            raw_uperpen310400_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370400_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370400_puppi100->Fill(raw_par_bins);
                            raw_uperpen370400_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt370400_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi100->Fill(raw_par_bins);
                            raw_uperpen310400_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx04_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx04_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx46_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx46_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx68_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx68_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx812_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx812_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1218_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx1218_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1824_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx1824_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2428_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx2428_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2832_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx2832_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3236_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx3236_puppi100->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3640_puppi100->Fill(raw_par_bins);
                            raw_uperpennvtx3640_puppi100->Fill(raw_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5070_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5070_puppi150->Fill(raw_par_bins);
                            raw_uperpen5070_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt5070_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi150->Fill(raw_par_bins);
                            raw_uperpen50110_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 70 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale7090_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara7090_puppi150->Fill(raw_par_bins);
                            raw_uperpen7090_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt7090_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50110_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi150->Fill(raw_par_bins);
                            raw_uperpen50110_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90110_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90110_puppi150->Fill(raw_par_bins);
                            raw_uperpen90110_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt90110_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi150->Fill(raw_par_bins);
                            raw_uperpen50110_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 110 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale110130_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110130_puppi150->Fill(raw_par_bins);
                            raw_uperpen110130_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110130_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi150->Fill(raw_par_bins);
                            raw_uperpen110210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130150_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130150_puppi150->Fill(raw_par_bins);
                            raw_uperpen130150_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt130150_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi150->Fill(raw_par_bins);
                            raw_uperpen110210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 150 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale150170_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara150170_puppi150->Fill(raw_par_bins);
                            raw_uperpen150170_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt150170_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi150->Fill(raw_par_bins);
                            raw_uperpen110210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170190_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170190_puppi150->Fill(raw_par_bins);
                            raw_uperpen170190_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt170190_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi150->Fill(raw_par_bins);
                            raw_uperpen110210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 190 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale190210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara190210_puppi150->Fill(raw_par_bins);
                            raw_uperpen190210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt190210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi150->Fill(raw_par_bins);
                            raw_uperpen110210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210230_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210230_puppi150->Fill(raw_par_bins);
                            raw_uperpen210230_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210230_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi150->Fill(raw_par_bins);
                            raw_uperpen210310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 230 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale230250_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara230250_puppi150->Fill(raw_par_bins);
                            raw_uperpen230250_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt230250_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi150->Fill(raw_par_bins);
                            raw_uperpen210310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250270_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250270_puppi150->Fill(raw_par_bins);
                            raw_uperpen250270_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt250270_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi150->Fill(raw_par_bins);
                            raw_uperpen210310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 270 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale270290_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara270290_puppi150->Fill(raw_par_bins);
                            raw_uperpen270290_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt270290_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi150->Fill(raw_par_bins);
                            raw_uperpen210310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290310_puppi150->Fill(raw_par_bins);
                            raw_uperpen290310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt290310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi150->Fill(raw_par_bins);
                            raw_uperpen210310_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 310 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale310330_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310330_puppi150->Fill(raw_par_bins);
                            raw_uperpen310330_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt310330_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi150->Fill(raw_par_bins);
                            raw_uperpen310400_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330350_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330350_puppi150->Fill(raw_par_bins);
                            raw_uperpen330350_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt330350_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi150->Fill(raw_par_bins);
                            raw_uperpen310400_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 350 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale350370_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara350370_puppi150->Fill(raw_par_bins);
                            raw_uperpen350370_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt350370_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi150->Fill(raw_par_bins);
                            raw_uperpen310400_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370400_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370400_puppi150->Fill(raw_par_bins);
                            raw_uperpen370400_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt370400_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi150->Fill(raw_par_bins);
                            raw_uperpen310400_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx04_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx04_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx46_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx46_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx68_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx68_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx812_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx812_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1218_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx1218_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1824_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx1824_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2428_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx2428_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2832_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx2832_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3236_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx3236_puppi150->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3640_puppi150->Fill(raw_par_bins);
                            raw_uperpennvtx3640_puppi150->Fill(raw_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 70)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5070_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5070_puppi200->Fill(raw_par_bins);
                            raw_uperpen5070_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt5070_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi200->Fill(raw_par_bins);
                            raw_uperpen50110_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 70 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale7090_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara7090_puppi200->Fill(raw_par_bins);
                            raw_uperpen7090_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt7090_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50110_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi200->Fill(raw_par_bins);
                            raw_uperpen50110_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 110)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90110_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90110_puppi200->Fill(raw_par_bins);
                            raw_uperpen90110_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt90110_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50110_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50110_puppi200->Fill(raw_par_bins);
                            raw_uperpen50110_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt50110_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 110 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale110130_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110130_puppi200->Fill(raw_par_bins);
                            raw_uperpen110130_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110130_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi200->Fill(raw_par_bins);
                            raw_uperpen110210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 150)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130150_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130150_puppi200->Fill(raw_par_bins);
                            raw_uperpen130150_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt130150_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi200->Fill(raw_par_bins);
                            raw_uperpen110210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 150 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale150170_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara150170_puppi200->Fill(raw_par_bins);
                            raw_uperpen150170_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt150170_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi200->Fill(raw_par_bins);
                            raw_uperpen110210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 190)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170190_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170190_puppi200->Fill(raw_par_bins);
                            raw_uperpen170190_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt170190_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi200->Fill(raw_par_bins);
                            raw_uperpen110210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 190 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale190210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara190210_puppi200->Fill(raw_par_bins);
                            raw_uperpen190210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt190210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale110210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara110210_puppi200->Fill(raw_par_bins);
                            raw_uperpen110210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt110210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 230)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210230_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210230_puppi200->Fill(raw_par_bins);
                            raw_uperpen210230_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210230_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi200->Fill(raw_par_bins);
                            raw_uperpen210310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 230 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale230250_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara230250_puppi200->Fill(raw_par_bins);
                            raw_uperpen230250_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt230250_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi200->Fill(raw_par_bins);
                            raw_uperpen210310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 270)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250270_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250270_puppi200->Fill(raw_par_bins);
                            raw_uperpen250270_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt250270_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi200->Fill(raw_par_bins);
                            raw_uperpen210310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 270 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale270290_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara270290_puppi200->Fill(raw_par_bins);
                            raw_uperpen270290_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt270290_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi200->Fill(raw_par_bins);
                            raw_uperpen210310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 310)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290310_puppi200->Fill(raw_par_bins);
                            raw_uperpen290310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt290310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210310_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210310_puppi200->Fill(raw_par_bins);
                            raw_uperpen210310_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210310_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 310 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale310330_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310330_puppi200->Fill(raw_par_bins);
                            raw_uperpen310330_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt310330_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi200->Fill(raw_par_bins);
                            raw_uperpen310400_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 350)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330350_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330350_puppi200->Fill(raw_par_bins);
                            raw_uperpen330350_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt330350_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi200->Fill(raw_par_bins);
                            raw_uperpen310400_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 350 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale350370_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara350370_puppi200->Fill(raw_par_bins);
                            raw_uperpen350370_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt350370_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi200->Fill(raw_par_bins);
                            raw_uperpen310400_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 400)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370400_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370400_puppi200->Fill(raw_par_bins);
                            raw_uperpen370400_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt370400_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale310400_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara310400_puppi200->Fill(raw_par_bins);
                            raw_uperpen310400_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt310400_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }


                        //vs nvtx
                        if (nGoodVtx > 0 && nGoodVtx < 4)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx04_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx04_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 4 && nGoodVtx < 6)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx46_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx46_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 6 && nGoodVtx < 8)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx68_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx68_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 8 && nGoodVtx < 12)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx812_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx812_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 12 && nGoodVtx < 18)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1218_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx1218_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 18 && nGoodVtx < 24)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx1824_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx1824_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 24 && nGoodVtx < 28)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2428_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx2428_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 28 && nGoodVtx < 32)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx2832_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx2832_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 32 && nGoodVtx < 36)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3236_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx3236_puppi200->Fill(raw_pen_bins);
                        }

                        if (nGoodVtx > 36 && nGoodVtx < 40)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            
                            raw_uparanvtx3640_puppi200->Fill(raw_par_bins);
                            raw_uperpennvtx3640_puppi200->Fill(raw_pen_bins);
                        }
                    }

                    if (puppiMET_pt > 100)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5090_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5090_puppi100->Fill(raw_par_bins);
                            raw_uperpen5090_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt5090_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50130_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi100->Fill(raw_par_bins);
                            raw_uperpen50130_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90130_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90130_puppi100->Fill(raw_par_bins);
                            raw_uperpen90130_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt90130_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50130_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi100->Fill(raw_par_bins);
                            raw_uperpen50130_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130170_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130170_puppi100->Fill(raw_par_bins);
                            raw_uperpen130170_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt130170_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi100->Fill(raw_par_bins);
                            raw_uperpen130210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170210_puppi100->Fill(raw_par_bins);
                            raw_uperpen170210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt170210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi100->Fill(raw_par_bins);
                            raw_uperpen130210_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210250_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210250_puppi100->Fill(raw_par_bins);
                            raw_uperpen210250_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210250_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale210290_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi100->Fill(raw_par_bins);
                            raw_uperpen210290_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);


                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250290_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250290_puppi100->Fill(raw_par_bins);
                            raw_uperpen250290_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt250290_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210290_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi100->Fill(raw_par_bins);
                            raw_uperpen210290_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290330_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290330_puppi100->Fill(raw_par_bins);
                            raw_uperpen290330_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt290330_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi100->Fill(raw_par_bins);
                            raw_uperpen290370_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330370_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330370_puppi100->Fill(raw_par_bins);
                            raw_uperpen330370_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt330370_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi100->Fill(raw_par_bins);
                            raw_uperpen290370_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370410_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi100->Fill(raw_par_bins);
                            raw_uperpen370410_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370450_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370450_puppi100->Fill(raw_par_bins);
                            raw_uperpen370450_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt370450_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 410 && (*phoEt)[irawPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale410450_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara410450_puppi100->Fill(raw_par_bins);
                            raw_uperpen410450_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt410450_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370410_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi100->Fill(raw_par_bins);
                            raw_uperpen370410_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 450 && (*phoEt)[irawPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale450490_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450490_puppi100->Fill(raw_par_bins);
                            raw_uperpen450490_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt450490_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi100->Fill(raw_par_bins);
                            raw_uperpen450530_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 490 && (*phoEt)[irawPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale490530_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara490530_puppi100->Fill(raw_par_bins);
                            raw_uperpen490530_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt490530_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi100->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi100->Fill(raw_par_bins);
                            raw_uperpen450530_puppi100->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi100->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }
                    }

                    if (puppiMET_pt > 150)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5090_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5090_puppi150->Fill(raw_par_bins);
                            raw_uperpen5090_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt5090_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50130_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi150->Fill(raw_par_bins);
                            raw_uperpen50130_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90130_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90130_puppi150->Fill(raw_par_bins);
                            raw_uperpen90130_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt90130_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50130_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi150->Fill(raw_par_bins);
                            raw_uperpen50130_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130170_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130170_puppi150->Fill(raw_par_bins);
                            raw_uperpen130170_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt130170_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi150->Fill(raw_par_bins);
                            raw_uperpen130210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170210_puppi150->Fill(raw_par_bins);
                            raw_uperpen170210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt170210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi150->Fill(raw_par_bins);
                            raw_uperpen130210_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210250_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210250_puppi150->Fill(raw_par_bins);
                            raw_uperpen210250_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210250_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale210290_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi150->Fill(raw_par_bins);
                            raw_uperpen210290_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);


                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250290_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250290_puppi150->Fill(raw_par_bins);
                            raw_uperpen250290_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt250290_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210290_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi150->Fill(raw_par_bins);
                            raw_uperpen210290_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290330_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290330_puppi150->Fill(raw_par_bins);
                            raw_uperpen290330_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt290330_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi150->Fill(raw_par_bins);
                            raw_uperpen290370_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330370_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330370_puppi150->Fill(raw_par_bins);
                            raw_uperpen330370_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt330370_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi150->Fill(raw_par_bins);
                            raw_uperpen290370_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370410_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi150->Fill(raw_par_bins);
                            raw_uperpen370410_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370450_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370450_puppi150->Fill(raw_par_bins);
                            raw_uperpen370450_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt370450_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 410 && (*phoEt)[irawPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale410450_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara410450_puppi150->Fill(raw_par_bins);
                            raw_uperpen410450_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt410450_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370410_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi150->Fill(raw_par_bins);
                            raw_uperpen370410_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 450 && (*phoEt)[irawPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale450490_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450490_puppi150->Fill(raw_par_bins);
                            raw_uperpen450490_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt450490_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi150->Fill(raw_par_bins);
                            raw_uperpen450530_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 490 && (*phoEt)[irawPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale490530_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara490530_puppi150->Fill(raw_par_bins);
                            raw_uperpen490530_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt490530_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi150->Fill(raw_par_bins);
                            raw_uperpen450530_puppi150->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }
                    }

                    if (puppiMET_pt > 200)
                    {
                        if ((*phoEt)[irawPho] > 50 && (*phoEt)[irawPho] < 90)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale5090_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara5090_puppi200->Fill(raw_par_bins);
                            raw_uperpen5090_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt5090_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale50130_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi200->Fill(raw_par_bins);
                            raw_uperpen50130_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 90 && (*phoEt)[irawPho] < 130)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale90130_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara90130_puppi200->Fill(raw_par_bins);
                            raw_uperpen90130_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt90130_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale50130_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara50130_puppi200->Fill(raw_par_bins);
                            raw_uperpen50130_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt50130_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 130 && (*phoEt)[irawPho] < 170)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale130170_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130170_puppi200->Fill(raw_par_bins);
                            raw_uperpen130170_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt130170_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi200->Fill(raw_par_bins);
                            raw_uperpen130210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 170 && (*phoEt)[irawPho] < 210)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale170210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara170210_puppi200->Fill(raw_par_bins);
                            raw_uperpen170210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt170210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale130210_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara130210_puppi200->Fill(raw_par_bins);
                            raw_uperpen130210_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt130210_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 210 && (*phoEt)[irawPho] < 250)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale210250_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210250_puppi200->Fill(raw_par_bins);
                            raw_uperpen210250_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210250_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);


                            raw_abs_scale210290_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi200->Fill(raw_par_bins);
                            raw_uperpen210290_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);


                        }

                        if ((*phoEt)[irawPho] > 250 && (*phoEt)[irawPho] < 290)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale250290_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara250290_puppi200->Fill(raw_par_bins);
                            raw_uperpen250290_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt250290_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale210290_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara210290_puppi200->Fill(raw_par_bins);
                            raw_uperpen210290_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt210290_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 290 && (*phoEt)[irawPho] < 330)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale290330_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290330_puppi200->Fill(raw_par_bins);
                            raw_uperpen290330_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt290330_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi200->Fill(raw_par_bins);
                            raw_uperpen290370_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 330 && (*phoEt)[irawPho] < 370)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale330370_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara330370_puppi200->Fill(raw_par_bins);
                            raw_uperpen330370_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt330370_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale290370_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara290370_puppi200->Fill(raw_par_bins);
                            raw_uperpen290370_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt290370_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 370 && (*phoEt)[irawPho] < 410)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale370410_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi200->Fill(raw_par_bins);
                            raw_uperpen370410_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370450_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370450_puppi200->Fill(raw_par_bins);
                            raw_uperpen370450_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt370450_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 410 && (*phoEt)[irawPho] < 450)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale410450_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara410450_puppi200->Fill(raw_par_bins);
                            raw_uperpen410450_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt410450_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale370410_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara370410_puppi200->Fill(raw_par_bins);
                            raw_uperpen370410_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt370410_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 450 && (*phoEt)[irawPho] < 490)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale450490_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450490_puppi200->Fill(raw_par_bins);
                            raw_uperpen450490_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt450490_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi200->Fill(raw_par_bins);
                            raw_uperpen450530_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }

                        if ((*phoEt)[irawPho] > 490 && (*phoEt)[irawPho] < 530)
                        {
                            //cal photon_px and photon_py
                            rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                            rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                            raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                            raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                            
                            raw_abs_scale490530_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara490530_puppi200->Fill(raw_par_bins);
                            raw_uperpen490530_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt490530_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);

                            raw_abs_scale450530_puppi200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                            raw_upara450530_puppi200->Fill(raw_par_bins);
                            raw_uperpen450530_puppi200->Fill(raw_pen_bins);
                            raw_uparaqt450530_puppi200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                        }
                    }


                }


                for (int ipho = 0; ipho < nPho; ipho++)
                {
                    //METs 2D compairson
                    if (puppiMET_pt > 100)
                    {
                        pf_vs_puppi100->Fill(pfMET_pt, puppiMET_pt);

                        pf_vs_rawpuppi100->Fill(pfMET_pt, rawPuppiMET_pt);
                        
                        raw_vs_puppi100->Fill(rawMET_pt, puppiMET_pt);

                        raw_vs_rawpuppi100->Fill(rawMET_pt, rawPuppiMET_pt);

                        rawpuppi_vs_puppi100->Fill(rawPuppiMET_pt, puppiMET_pt);

                        pf_vs_raw100->Fill(pfMET_pt, rawPuppiMET_pt);
                    }

                    if (puppiMET_pt > 200)
                    {
                        pf_vs_puppi200->Fill(pfMET_pt, puppiMET_pt);

                        pf_vs_rawpuppi200->Fill(pfMET_pt, rawPuppiMET_pt);
                        
                        raw_vs_puppi200->Fill(rawMET_pt, puppiMET_pt);

                        raw_vs_rawpuppi200->Fill(rawMET_pt, rawPuppiMET_pt);

                        rawpuppi_vs_puppi200->Fill(rawPuppiMET_pt, puppiMET_pt);

                        pf_vs_raw200->Fill(pfMET_pt, rawPuppiMET_pt);
                    }
                }
            }
        }
    }



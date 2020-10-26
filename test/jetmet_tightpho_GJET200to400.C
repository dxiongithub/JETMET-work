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

    TFile *JETMETPuppi = new TFile("jetmet_tightpho_GJET200to400.root", "RECREATE");

    //puppimet with eaxctly one tight photon, at least one jet, no loose electron or muon
    TH1F *puppimetphi = new TH1F("puppimetphi", "puppimetphi", 200, 25, 25);
	puppimetphi->GetXaxis()->SetTitle("phi");
	puppimetphi->GetYaxis()->SetTitle("Entries");

    TH1F *puppimetpt = new TH1F("puppimetpt", "puppimetpt", 200, 0, 400);
	puppimetpt->GetXaxis()->SetTitle("p_{T} (GeV)");
	puppimetpt->GetYaxis()->SetTitle("Entries");

    TH1F *puppimet_px_cal = new TH1F("puppimet_px_cal", "puppimet_px_cal", 200, 25, 25);
	puppimet_px_cal->GetXaxis()->SetTitle("px");
	puppimet_px_cal->GetYaxis()->SetTitle("Entries");

    TH1F *puppimet_py_cal = new TH1F("puppimet_py_cal", "puppimet_py_cal", 200, 25, 25);
	puppimet_py_cal->GetXaxis()->SetTitle("py");
	puppimet_py_cal->GetYaxis()->SetTitle("Entries");

    TH1F *puppimet_px = new TH1F("puppimet_px", "puppimet_px", 200, 25, 25);
	puppimet_px->GetXaxis()->SetTitle("px");
	puppimet_px->GetYaxis()->SetTitle("Entries");

    TH1F *puppimet_py = new TH1F("puppimet_py", "puppimet_py", 200, 25, 25);
	puppimet_py->GetXaxis()->SetTitle("py");
	puppimet_py->GetYaxis()->SetTitle("Entries");

    TH1F *uperpen_puppi = new TH1F("uperpen_puppi", "uperpen_puppi", 200, 25, 25);
	uperpen_puppi->GetXaxis()->SetTitle("u#perp");
	uperpen_puppi->GetYaxis()->SetTitle("Entries");

    TH1F *uparall_puppi = new TH1F("uparall_puppi", "uparall_puppi", 200, 25, 25);
	uparall_puppi->GetXaxis()->SetTitle("u#parallel");
	uparall_puppi->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale = new TH1F("puppi_abs_scale", "puppi_abs_scale", 200, 25, 25);
	puppi_abs_scale->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_paraqt = new TH1F("puppi_paraqt", "puppi_paraqt", 200, 25, 25);
	puppi_paraqt->GetXaxis()->SetTitle("u#parallel + q_{T}");
	puppi_paraqt->GetYaxis()->SetTitle("Entries");

    TH1F *puppimet_tightphopt = new TH1F("puppimet_tightphopt", "puppimet_tightphopt", 200, 0, 400);
	puppimet_tightphopt->GetXaxis()->SetTitle("p_{T} (GeV)");
	puppimet_tightphopt->GetYaxis()->SetTitle("Entries");

    TH1F *puppipho_pxdist = new TH1F("puppipho_pxdist", "puppipho_pxdist", 200, 25, 25);
	puppipho_pxdist->GetXaxis()->SetTitle("p_{T} (GeV)");
	puppipho_pxdist->GetYaxis()->SetTitle("Entries");

    TH1F *puppipho_pydist = new TH1F("puppipho_pydist", "puppipho_pydist", 200, 25, 25);
	puppipho_pydist->GetXaxis()->SetTitle("p_{T} (GeV)");
	puppipho_pydist->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_response = new TH1F("puppi_response", "puppi_response", 200, 25, 25);
	puppi_response->GetXaxis()->SetTitle("q_{T} (GeV)");
	puppi_response->GetYaxis()->SetTitle("- <u#parallel> / <q_{T}>");


    TH1F *puppimet_nGoodvtx = new TH1F("puppimet_nvtx", "puppimet_nvtx", 200, 25, 25);
	puppimet_nGoodvtx->GetXaxis()->SetTitle("nvtx");
	puppimet_nGoodvtx->GetYaxis()->SetTitle("Entries");




    //puppimet absolute scale with different bins
    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */
    TH1F *puppi_abs_scale5060 = new TH1F("puppi_abs_scale5060", "puppi_abs_scale5060", 200, 25, 25);
	puppi_abs_scale5060->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale5060->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale8090 = new TH1F("puppi_abs_scale8090", "puppi_abs_scale8090", 200, 25, 25);
	puppi_abs_scale8090->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale8090->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale100110 = new TH1F("puppi_abs_scale100110", "puppi_abs_scale100110", 200, 25, 25);
	puppi_abs_scale100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale120130 = new TH1F("puppi_abs_scale120130", "puppi_abs_scale120130", 200, 25, 25);
	puppi_abs_scale120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale140150 = new TH1F("puppi_abs_scale140150", "puppi_abs_scale140150", 200, 25, 25);
	puppi_abs_scale140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale160170 = new TH1F("puppi_abs_scale160170", "puppi_abs_scale160170", 200, 25, 25);
	puppi_abs_scale160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale190200 = new TH1F("puppi_abs_scale190200", "puppi_abs_scale190200", 200, 25, 25);
	puppi_abs_scale190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale190200->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale220230 = new TH1F("puppi_abs_scale220230", "puppi_abs_scale220230", 200, 25, 25);
	puppi_abs_scale220230->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale220230->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale250260 = new TH1F("puppi_abs_scale250260", "puppi_abs_scale250260", 200, 25, 25);
	puppi_abs_scale250260->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale250260->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scaled260270 = new TH1F("puppi_abs_scaled260270", "puppi_abs_scaled260270", 200, 25, 25);
	puppi_abs_scaled260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scaled260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale280290 = new TH1F("puppi_abs_scale280290", "puppi_abs_scale280290", 200, 25, 25);
	puppi_abs_scale280290->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale280290->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale300310 = new TH1F("puppi_abs_scale300310", "puppi_abs_scale300310", 200, 25, 25);
	puppi_abs_scale300310->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale300310->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale310320 = new TH1F("puppi_abs_scale310320", "puppi_abs_scale310320", 200, 25, 25);
	puppi_abs_scale310320->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale310320->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale320330 = new TH1F("puppi_abs_scale320330", "puppi_abs_scale320330", 200, 25, 25);
	puppi_abs_scale320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale350360 = new TH1F("puppi_abs_scale350360", "puppi_abs_scale350360", 200, 25, 25);
	puppi_abs_scale350360->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale350360->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale380390 = new TH1F("puppi_abs_scale380390", "puppi_abs_scale380390", 200, 25, 25);
	puppi_abs_scale380390->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale380390->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale390400 = new TH1F("puppi_abs_scale390400", "puppi_abs_scale390400", 200, 25, 25);
	puppi_abs_scale390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale390400->GetYaxis()->SetTitle("Entries");


    //puppimet u_para with different bins
    TH1F *puppi_upara5060 = new TH1F("puppi_upara5060", "puppi_upara5060", 200, 25, 25);
	puppi_upara5060->GetXaxis()->SetTitle("u#parallel");
	puppi_upara5060->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara8090 = new TH1F("puppi_upara8090", "puppi_upara8090", 200, 25, 25);
	puppi_upara8090->GetXaxis()->SetTitle("u#parallel");
	puppi_upara8090->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara100110 = new TH1F("puppi_upara100110", "puppi_upara100110", 200, 25, 25);
	puppi_upara100110->GetXaxis()->SetTitle("u#parallel");
	puppi_upara100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara120130 = new TH1F("puppi_upara120130", "puppi_upara120130", 200, 25, 25);
	puppi_upara120130->GetXaxis()->SetTitle("u#parallel");
	puppi_upara120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara140150 = new TH1F("puppi_upara140150", "puppi_upara140150", 200, 25, 25);
	puppi_upara140150->GetXaxis()->SetTitle("u#parallel");
	puppi_upara140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara160170 = new TH1F("puppi_upara160170", "puppi_upara160170", 200, 25, 25);
	puppi_upara160170->GetXaxis()->SetTitle("u#parallel");
	puppi_upara160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara190200 = new TH1F("puppi_upara190200", "puppi_upara190200", 200, 25, 25);
	puppi_upara190200->GetXaxis()->SetTitle("u#parallel");
	puppi_upara190200->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara220230 = new TH1F("puppi_upara220230", "puppi_upara220230", 200, 25, 25);
	puppi_upara220230->GetXaxis()->SetTitle("u#parallel");
	puppi_upara220230->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara250260 = new TH1F("puppi_upara250260", "puppi_upara250260", 200, 25, 25);
	puppi_upara250260->GetXaxis()->SetTitle("u#parallel");
	puppi_upara250260->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara260270 = new TH1F("puppi_upara260270", "puppi_upara260270", 200, 25, 25);
	puppi_upara260270->GetXaxis()->SetTitle("u#parallel");
	puppi_upara260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara280290 = new TH1F("puppi_upara280290", "puppi_upara280290", 200, 25, 25);
	puppi_upara280290->GetXaxis()->SetTitle("u#parallel");
	puppi_upara280290->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara300310 = new TH1F("puppi_upara300310", "puppi_upara300310", 200, 25, 25);
	puppi_upara300310->GetXaxis()->SetTitle("u#parallel");
	puppi_upara300310->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara310320 = new TH1F("puppi_upara310320", "puppi_upara310320", 200, 25, 25);
	puppi_upara310320->GetXaxis()->SetTitle("u#parallel");
	puppi_upara310320->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara320330 = new TH1F("puppi_upara320330", "puppi_upara320330", 200, 25, 25);
	puppi_upara320330->GetXaxis()->SetTitle("u#parallel");
	puppi_upara320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara350360 = new TH1F("puppi_upara350360", "puppi_upara350360", 200, 25, 25);
	puppi_upara350360->GetXaxis()->SetTitle("u#parallel");
	puppi_upara350360->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara380390 = new TH1F("puppi_upara380390", "puppi_upara380390", 200, 25, 25);
	puppi_upara380390->GetXaxis()->SetTitle("u#parallel");
	puppi_upara380390->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara390400 = new TH1F("puppi_upara390400", "puppi_upara390400", 200, 25, 25);
	puppi_upara390400->GetXaxis()->SetTitle("u#parallel");
	puppi_upara390400->GetYaxis()->SetTitle("Entries");


    //puppimet u_perpen with different bins
    TH1F *puppi_uperpen5060 = new TH1F("puppi_uperpen5060", "puppi_uperpen5060", 200, 25, 25);
	puppi_uperpen5060->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen5060->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen8090 = new TH1F("puppi_uperpen8090", "puppi_uperpen8090", 200, 25, 25);
	puppi_uperpen8090->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen8090->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen100110 = new TH1F("puppi_uperpen100110", "puppi_uperpen100110", 200, 25, 25);
	puppi_uperpen100110->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen120130 = new TH1F("puppi_uperpen120130", "puppi_uperpen120130", 200, 25, 25);
	puppi_uperpen120130->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen140150 = new TH1F("puppi_uperpen140150", "puppi_uperpen140150", 200, 25, 25);
	puppi_uperpen140150->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen160170 = new TH1F("puppi_uperpen160170", "puppi_uperpen160170", 200, 25, 25);
	puppi_uperpen160170->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen190200 = new TH1F("puppi_uperpen190200", "puppi_uperpen190200", 200, 25, 25);
	puppi_uperpen190200->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen190200->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen220230 = new TH1F("puppi_uperpen220230", "puppi_uperpen220230", 200, 25, 25);
	puppi_uperpen220230->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen220230->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen250260 = new TH1F("puppi_uperpen250260", "puppi_uperpen250260", 200, 25, 25);
	puppi_uperpen250260->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen250260->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen260270 = new TH1F("puppi_uperpen260270", "puppi_uperpen260270", 200, 25, 25);
	puppi_uperpen260270->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen280290 = new TH1F("puppi_uperpen280290", "puppi_uperpen280290", 200, 25, 25);
	puppi_uperpen280290->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen280290->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen300310 = new TH1F("puppi_uperpen300310", "puppi_uperpen300310", 200, 25, 25);
	puppi_uperpen300310->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen300310->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen310320 = new TH1F("puppi_uperpen310320", "puppi_uperpen310320", 200, 25, 25);
	puppi_uperpen310320->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen310320->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen320330 = new TH1F("puppi_uperpen320330", "puppi_uperpen320330", 200, 25, 25);
	puppi_uperpen320330->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen350360 = new TH1F("puppi_uperpen350360", "puppi_uperpen350360", 200, 25, 25);
	puppi_uperpen350360->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen350360->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen380390 = new TH1F("puppi_uperpen380390", "puppi_uperpen380390", 200, 25, 25);
	puppi_uperpen380390->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen380390->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen390400 = new TH1F("puppi_uperpen390400", "puppi_uperpen390400", 200, 25, 25);
	puppi_uperpen390400->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen390400->GetYaxis()->SetTitle("Entries");


    //puppimet u_para with different nGoodvtx bins
    TH1F *puppi_uparanvtx04 = new TH1F("puppi_uparanvtx04", "puppi_uparanvtx04", 200, 25, 25);
	puppi_uparanvtx04->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx46 = new TH1F("puppi_uparanvtx46", "puppi_uparanvtx46", 200, 25, 25);
	puppi_uparanvtx46->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx68 = new TH1F("puppi_uparanvtx68", "puppi_uparanvtx68", 200, 25, 25);
	puppi_uparanvtx68->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx812 = new TH1F("puppi_uparanvtx812", "puppi_uparanvtx812", 200, 25, 25);
	puppi_uparanvtx812->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx1218 = new TH1F("puppi_uparanvtx1218", "puppi_uparanvtx1218", 200, 25, 25);
	puppi_uparanvtx1218->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx1824 = new TH1F("puppi_uparanvtx1824", "puppi_uparanvtx1824", 200, 25, 25);
	puppi_uparanvtx1824->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx2428 = new TH1F("puppi_uparanvtx2428", "puppi_uparanvtx2428", 200, 25, 25);
	puppi_uparanvtx2428->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx2832 = new TH1F("puppi_uparanvtx2832", "puppi_uparanvtx2832", 200, 25, 25);
	puppi_uparanvtx2832->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx3236 = new TH1F("puppi_uparanvtx3236", "puppi_uparanvtx3236", 200, 25, 25);
	puppi_uparanvtx3236->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparanvtx3640 = new TH1F("puppi_uparanvtx3640", "puppi_uparanvtx3640", 200, 25, 25);
	puppi_uparanvtx3640->GetXaxis()->SetTitle("u#parallel");
	puppi_uparanvtx3640->GetYaxis()->SetTitle("Entries");

    

    //puppimet u_perpen with different nGoodvtx bins
    TH1F *puppi_uperpennvtx04 = new TH1F("puppi_uperpennvtx04", "puppi_uperpennvtx04", 200, 25, 25);
	puppi_uperpennvtx04->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx46 = new TH1F("puppi_uperpennvtx46", "puppi_uperpennvtx46", 200, 25, 25);
	puppi_uperpennvtx46->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx68 = new TH1F("puppi_uperpennvtx68", "puppi_uperpennvtx68", 200, 25, 25);
	puppi_uperpennvtx68->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx812 = new TH1F("puppi_uperpennvtx812", "puppi_uperpennvtx812", 200, 25, 25);
	puppi_uperpennvtx812->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx1218 = new TH1F("puppi_uperpennvtx1218", "puppi_uperpennvtx1218", 200, 25, 25);
	puppi_uperpennvtx1218->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx1824 = new TH1F("puppi_uperpennvtx1824", "puppi_uperpennvtx1824", 200, 25, 25);
	puppi_uperpennvtx1824->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx2428 = new TH1F("puppi_uperpennvtx2428", "puppi_uperpennvtx2428", 200, 25, 25);
	puppi_uperpennvtx2428->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx2832 = new TH1F("puppi_uperpennvtx2832", "puppi_uperpennvtx2832", 200, 25, 25);
	puppi_uperpennvtx2832->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx3236 = new TH1F("puppi_uperpennvtx3236", "puppi_uperpennvtx3236", 200, 25, 25);
	puppi_uperpennvtx3236->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpennvtx3640 = new TH1F("puppi_uperpennvtx3640", "puppi_uperpennvtx3640", 200, 25, 25);
	puppi_uperpennvtx3640->GetXaxis()->SetTitle("u#perp");
	puppi_uperpennvtx3640->GetYaxis()->SetTitle("Entries");


    //puppimet u_paraqt with different bins
    TH1F *puppi_uparaqt5060 = new TH1F("puppi_uparaqt5060", "puppi_uparaqt5060", 200, 25, 25);
	puppi_uparaqt5060->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt5060->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt8090 = new TH1F("puppi_uparaqt8090", "puppi_uparaqt8090", 200, 25, 25);
	puppi_uparaqt8090->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt8090->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt100110 = new TH1F("puppi_uparaqt100110", "puppi_uparaqt100110", 200, 25, 25);
	puppi_uparaqt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt120130 = new TH1F("puppi_uparaqt120130", "puppi_uparaqt120130", 200, 25, 25);
	puppi_uparaqt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt140150 = new TH1F("puppi_uparaqt140150", "puppi_uparaqt140150", 200, 25, 25);
	puppi_uparaqt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt160170 = new TH1F("puppi_uparaqt160170", "puppi_uparaqt160170", 200, 25, 25);
	puppi_uparaqt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt190200 = new TH1F("puppi_uparaqt190200", "puppi_uparaqt190200", 200, 25, 25);
	puppi_uparaqt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt190200->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt220230 = new TH1F("puppi_uparaqt220230", "puppi_uparaqt220230", 200, 25, 25);
	puppi_uparaqt220230->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt220230->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt250260 = new TH1F("puppi_uparaqt250260", "puppi_uparaqt250260", 200, 25, 25);
	puppi_uparaqt250260->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt250260->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt260270 = new TH1F("puppi_uparaqt260270", "puppi_uparaqt260270", 200, 25, 25);
	puppi_uparaqt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt280290 = new TH1F("puppi_uparaqt280290", "puppi_uparaqt280290", 200, 25, 25);
	puppi_uparaqt280290->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt280290->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt300310 = new TH1F("puppi_uparaqt300310", "puppi_uparaqt300310", 200, 25, 25);
	puppi_uparaqt300310->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt300310->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt310320 = new TH1F("puppi_uparaqt310320", "puppi_uparaqt310320", 200, 25, 25);
	puppi_uparaqt310320->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt310320->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt320330 = new TH1F("puppi_uparaqt320330", "puppi_uparaqt320330", 200, 25, 25);
	puppi_uparaqt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt350360 = new TH1F("puppi_uparaqt350360", "puppi_uparaqt350360", 200, 25, 25);
	puppi_uparaqt350360->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt350360->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt380390 = new TH1F("puppi_uparaqt380390", "puppi_uparaqt380390", 200, 25, 25);
	puppi_uparaqt380390->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt380390->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt390400 = new TH1F("puppi_uparaqt390400", "puppi_uparaqt390400", 200, 25, 25);
	puppi_uparaqt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt390400->GetYaxis()->SetTitle("Entries");

    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */


    /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */

    //puppi vs puppi pt > 100
    TH1F *puppi_abs_scale_puppipt100110 = new TH1F("puppi_abs_scale_puppipt100110", "puppi_abs_scale_puppipt100110", 200, 25, 25);
	puppi_abs_scale_puppipt100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt110120 = new TH1F("puppi_abs_scale_puppipt110120", "puppi_abs_scale_puppipt110120", 200, 25, 25);
	puppi_abs_scale_puppipt110120->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt120130 = new TH1F("puppi_abs_scale_puppipt120130", "puppi_abs_scale_puppipt120130", 200, 25, 25);
	puppi_abs_scale_puppipt120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt130140 = new TH1F("puppi_abs_scale_puppipt130140", "puppi_abs_scale_puppipt130140", 200, 25, 25);
	puppi_abs_scale_puppipt130140->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt140150 = new TH1F("puppi_abs_scale_puppipt140150", "puppi_abs_scale_puppipt140150", 200, 25, 25);
	puppi_abs_scale_puppipt140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt150160 = new TH1F("puppi_abs_scale_puppipt150160", "puppi_abs_scale_puppipt150160", 200, 25, 25);
	puppi_abs_scale_puppipt150160->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt160170 = new TH1F("puppi_abs_scale_puppipt160170", "puppi_abs_scale_puppipt160170", 200, 25, 25);
	puppi_abs_scale_puppipt160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt170180 = new TH1F("puppi_abs_scale_puppipt170180", "puppi_abs_scale_puppipt170180", 200, 25, 25);
	puppi_abs_scale_puppipt170180->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt180190 = new TH1F("puppi_abs_scale_puppipt180190", "puppi_abs_scale_puppipt180190", 200, 25, 25);
	puppi_abs_scale_puppipt180190->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt190200 = new TH1F("puppi_abs_scale_puppipt190200", "puppi_abs_scale_puppipt190200", 200, 25, 25);
	puppi_abs_scale_puppipt190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt190200->GetYaxis()->SetTitle("Entries");


    //puppi vs puppi pt > 200
    TH1F *puppi_abs_scale_puppipt200210 = new TH1F("puppi_abs_scale_puppipt200210", "puppi_abs_scale_puppipt200210", 200, 25, 25);
	puppi_abs_scale_puppipt200210->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt230240 = new TH1F("puppi_abs_scale_puppipt230240", "puppi_abs_scale_puppipt230240", 200, 25, 25);
	puppi_abs_scale_puppipt230240->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt260270 = new TH1F("puppi_abs_scale_puppipt260270", "puppi_abs_scale_puppipt260270", 200, 25, 25);
	puppi_abs_scale_puppipt260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt290300 = new TH1F("puppi_abs_scale_puppipt290300", "puppi_abs_scale_puppipt290300", 200, 25, 25);
	puppi_abs_scale_puppipt290300->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt320330 = new TH1F("puppi_abs_scale_puppipt320330", "puppi_abs_scale_puppipt320330", 200, 25, 25);
	puppi_abs_scale_puppipt320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt330340 = new TH1F("puppi_abs_scale_puppipt330340", "puppi_abs_scale_puppipt330340", 200, 25, 25);
	puppi_abs_scale_puppipt330340->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt360370 = new TH1F("puppi_abs_scale_puppipt360370", "puppi_abs_scale_puppipt360370", 200, 25, 25);
	puppi_abs_scale_puppipt360370->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt390400 = new TH1F("puppi_abs_scale_puppipt390400", "puppi_abs_scale_puppipt390400", 200, 25, 25);
	puppi_abs_scale_puppipt390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_abs_scale_puppipt420430 = new TH1F("puppi_abs_scale_puppipt420430", "puppi_abs_scale_puppipt420430", 200, 25, 25);
	puppi_abs_scale_puppipt420430->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	puppi_abs_scale_puppipt420430->GetYaxis()->SetTitle("Entries");



    //puppiupara vs puppi pt > 100
    TH1F *puppi_upara_puppipt100110 = new TH1F("puppi_upara_puppipt100110", "puppi_upara_puppipt100110", 200, 25, 25);
	puppi_upara_puppipt100110->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt110120 = new TH1F("puppi_upara_puppipt110120", "puppi_upara_puppipt110120", 200, 25, 25);
	puppi_upara_puppipt110120->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt120130 = new TH1F("puppi_upara_puppipt120130", "puppi_upara_puppipt120130", 200, 25, 25);
	puppi_upara_puppipt120130->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt130140 = new TH1F("puppi_upara_puppipt130140", "puppi_upara_puppipt130140", 200, 25, 25);
	puppi_upara_puppipt130140->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt140150 = new TH1F("puppi_upara_puppipt140150", "puppi_upara_puppipt140150", 200, 25, 25);
	puppi_upara_puppipt140150->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt150160 = new TH1F("puppi_upara_puppipt150160", "puppi_upara_puppipt150160", 200, 25, 25);
	puppi_upara_puppipt150160->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt160170 = new TH1F("puppi_upara_puppipt160170", "puppi_upara_puppipt160170", 200, 25, 25);
	puppi_upara_puppipt160170->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt170180 = new TH1F("puppi_upara_puppipt170180", "puppi_upara_puppipt170180", 200, 25, 25);
	puppi_upara_puppipt170180->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt180190 = new TH1F("puppi_upara_puppipt180190", "puppi_upara_puppipt180190", 200, 25, 25);
	puppi_upara_puppipt180190->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt190200 = new TH1F("puppi_upara_puppipt190200", "puppi_upara_puppipt190200", 200, 25, 25);
	puppi_upara_puppipt190200->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt190200->GetYaxis()->SetTitle("Entries");


    //puppiupara vs puppi pt > 200
    TH1F *puppi_upara_puppipt200210 = new TH1F("puppi_upara_puppipt200210", "puppi_upara_puppipt200210", 200, 25, 25);
	puppi_upara_puppipt200210->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt230240 = new TH1F("puppi_upara_puppipt230240", "puppi_upara_puppipt230240", 200, 25, 25);
	puppi_upara_puppipt230240->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt260270 = new TH1F("puppi_upara_puppipt260270", "puppi_upara_puppipt260270", 200, 25, 25);
	puppi_upara_puppipt260270->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt290300 = new TH1F("puppi_upara_puppipt290300", "puppi_upara_puppipt290300", 200, 25, 25);
	puppi_upara_puppipt290300->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt320330 = new TH1F("puppi_upara_puppipt320330", "puppi_upara_puppipt320330", 200, 25, 25);
	puppi_upara_puppipt320330->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt330340 = new TH1F("puppi_upara_puppipt330340", "puppi_upara_puppipt330340", 200, 25, 25);
	puppi_upara_puppipt330340->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt360370 = new TH1F("puppi_upara_puppipt360370", "puppi_upara_puppipt360370", 200, 25, 25);
	puppi_upara_puppipt360370->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt390400 = new TH1F("puppi_upara_puppipt390400", "puppi_upara_puppipt390400", 200, 25, 25);
	puppi_upara_puppipt390400->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_upara_puppipt420430 = new TH1F("puppi_upara_puppipt420430", "puppi_upara_puppipt420430", 200, 25, 25);
	puppi_upara_puppipt420430->GetXaxis()->SetTitle("u#parallel");
	puppi_upara_puppipt420430->GetYaxis()->SetTitle("Entries");


    //puppiuperpen vs puppi pt > 100
    TH1F *puppi_uperpen_puppipt100110 = new TH1F("puppi_uperpen_puppipt100110", "puppi_uperpen_puppipt100110", 200, 25, 25);
	puppi_uperpen_puppipt100110->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt110120 = new TH1F("puppi_uperpen_puppipt110120", "puppi_uperpen_puppipt110120", 200, 25, 25);
	puppi_uperpen_puppipt110120->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt120130 = new TH1F("puppi_uperpen_puppipt120130", "puppi_uperpen_puppipt120130", 200, 25, 25);
	puppi_uperpen_puppipt120130->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt130140 = new TH1F("puppi_uperpen_puppipt130140", "puppi_uperpen_puppipt130140", 200, 25, 25);
	puppi_uperpen_puppipt130140->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt140150 = new TH1F("puppi_uperpen_puppipt140150", "puppi_uperpen_puppipt140150", 200, 25, 25);
	puppi_uperpen_puppipt140150->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt150160 = new TH1F("puppi_uperpen_puppipt150160", "puppi_uperpen_puppipt150160", 200, 25, 25);
	puppi_uperpen_puppipt150160->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt160170 = new TH1F("puppi_uperpen_puppipt160170", "puppi_uperpen_puppipt160170", 200, 25, 25);
	puppi_uperpen_puppipt160170->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt170180 = new TH1F("puppi_uperpen_puppipt170180", "puppi_uperpen_puppipt170180", 200, 25, 25);
	puppi_uperpen_puppipt170180->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt180190 = new TH1F("puppi_uperpen_puppipt180190", "puppi_uperpen_puppipt180190", 200, 25, 25);
	puppi_uperpen_puppipt180190->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt190200 = new TH1F("puppi_uperpen_puppipt190200", "puppi_uperpen_puppipt190200", 200, 25, 25);
	puppi_uperpen_puppipt190200->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt190200->GetYaxis()->SetTitle("Entries");


    //puppiuperpen vs puppi pt > 200
    TH1F *puppi_uperpen_puppipt200210 = new TH1F("puppi_uperpen_puppipt200210", "puppi_uperpen_puppipt200210", 200, 25, 25);
	puppi_uperpen_puppipt200210->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt230240 = new TH1F("puppi_uperpen_puppipt230240", "puppi_uperpen_puppipt230240", 200, 25, 25);
	puppi_uperpen_puppipt230240->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt260270 = new TH1F("puppi_uperpen_puppipt260270", "puppi_uperpen_puppipt260270", 200, 25, 25);
	puppi_uperpen_puppipt260270->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt290300 = new TH1F("puppi_uperpen_puppipt290300", "puppi_uperpen_puppipt290300", 200, 25, 25);
	puppi_uperpen_puppipt290300->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt320330 = new TH1F("puppi_uperpen_puppipt320330", "puppi_uperpen_puppipt320330", 200, 25, 25);
	puppi_uperpen_puppipt320330->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt330340 = new TH1F("puppi_uperpen_puppipt330340", "puppi_uperpen_puppipt330340", 200, 25, 25);
	puppi_uperpen_puppipt330340->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt360370 = new TH1F("puppi_uperpen_puppipt360370", "puppi_uperpen_puppipt360370", 200, 25, 25);
	puppi_uperpen_puppipt360370->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt390400 = new TH1F("puppi_uperpen_puppipt390400", "puppi_uperpen_puppipt390400", 200, 25, 25);
	puppi_uperpen_puppipt390400->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uperpen_puppipt420430 = new TH1F("puppi_uperpen_puppipt420430", "puppi_uperpen_puppipt420430", 200, 25, 25);
	puppi_uperpen_puppipt420430->GetXaxis()->SetTitle("u#perp");
	puppi_uperpen_puppipt420430->GetYaxis()->SetTitle("Entries");



    //puppiuparaqt vs puppi pt > 100
    TH1F *puppi_uparaqt_puppipt100110 = new TH1F("puppi_uparaqt_puppipt100110", "puppi_uparaqt_puppipt100110", 200, 25, 25);
	puppi_uparaqt_puppipt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt110120 = new TH1F("puppi_uparaqt_puppipt110120", "puppi_uparaqt_puppipt110120", 200, 25, 25);
	puppi_uparaqt_puppipt110120->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt120130 = new TH1F("puppi_uparaqt_puppipt120130", "puppi_uparaqt_puppipt120130", 200, 25, 25);
	puppi_uparaqt_puppipt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt130140 = new TH1F("puppi_uparaqt_puppipt130140", "puppi_uparaqt_puppipt130140", 200, 25, 25);
	puppi_uparaqt_puppipt130140->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt140150 = new TH1F("puppi_uparaqt_puppipt140150", "puppi_uparaqt_puppipt140150", 200, 25, 25);
	puppi_uparaqt_puppipt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt150160 = new TH1F("puppi_uparaqt_puppipt150160", "puppi_uparaqt_puppipt150160", 200, 25, 25);
	puppi_uparaqt_puppipt150160->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt160170 = new TH1F("puppi_uparaqt_puppipt160170", "puppi_uparaqt_puppipt160170", 200, 25, 25);
	puppi_uparaqt_puppipt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt170180 = new TH1F("puppi_uparaqt_puppipt170180", "puppi_uparaqt_puppipt170180", 200, 25, 25);
	puppi_uparaqt_puppipt170180->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt180190 = new TH1F("puppi_uparaqt_puppipt180190", "puppi_uparaqt_puppipt180190", 200, 25, 25);
	puppi_uparaqt_puppipt180190->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt190200 = new TH1F("puppi_uparaqt_puppipt190200", "puppi_uparaqt_puppipt190200", 200, 25, 25);
	puppi_uparaqt_puppipt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt190200->GetYaxis()->SetTitle("Entries");


    //puppiuparaqt vs puppi pt > 200
    TH1F *puppi_uparaqt_puppipt200210 = new TH1F("puppi_uparaqt_puppipt200210", "puppi_uparaqt_puppipt200210", 200, 25, 25);
	puppi_uparaqt_puppipt200210->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt230240 = new TH1F("puppi_uparaqt_puppipt230240", "puppi_uparaqt_puppipt230240", 200, 25, 25);
	puppi_uparaqt_puppipt230240->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt260270 = new TH1F("puppi_uparaqt_puppipt260270", "puppi_uparaqt_puppipt260270", 200, 25, 25);
	puppi_uparaqt_puppipt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt290300 = new TH1F("puppi_uparaqt_puppipt290300", "puppi_uparaqt_puppipt290300", 200, 25, 25);
	puppi_uparaqt_puppipt290300->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt320330 = new TH1F("puppi_uparaqt_puppipt320330", "puppi_uparaqt_puppipt320330", 200, 25, 25);
	puppi_uparaqt_puppipt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt330340 = new TH1F("puppi_uparaqt_puppipt330340", "puppi_uparaqt_puppipt330340", 200, 25, 25);
	puppi_uparaqt_puppipt330340->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt360370 = new TH1F("puppi_uparaqt_puppipt360370", "puppi_uparaqt_puppipt360370", 200, 25, 25);
	puppi_uparaqt_puppipt360370->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt390400 = new TH1F("puppi_uparaqt_puppipt390400", "puppi_uparaqt_puppipt390400", 200, 25, 25);
	puppi_uparaqt_puppipt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *puppi_uparaqt_puppipt420430 = new TH1F("puppi_uparaqt_puppipt420430", "puppi_uparaqt_puppipt420430", 200, 25, 25);
	puppi_uparaqt_puppipt420430->GetXaxis()->SetTitle("u#parallel+q_{T}");
	puppi_uparaqt_puppipt420430->GetYaxis()->SetTitle("Entries");

    
     /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */



    //RAWpuppimet with eaxctly one tight photon, at least one jet, no loose electron or muon
    TH1F *rawpuppimetphi = new TH1F("rawpuppimetphi", "rawpuppimetphi", 200, 25, 25);
	rawpuppimetphi->GetXaxis()->SetTitle("phi");
	rawpuppimetphi->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimetpt = new TH1F("rawpuppimetpt", "rawpuppimetpt", 200, 0, 400);
	rawpuppimetpt->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpuppimetpt->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimet_px_cal = new TH1F("rawpuppimet_px_cal", "rawpuppimet_px_cal", 200, 25, 25);
	rawpuppimet_px_cal->GetXaxis()->SetTitle("px");
	rawpuppimet_px_cal->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimet_py_cal = new TH1F("rawpuppimet_py_cal", "rawpuppimet_py_cal", 200, 25, 25);
	rawpuppimet_py_cal->GetXaxis()->SetTitle("py");
	rawpuppimet_py_cal->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimet_px = new TH1F("rawpuppimet_px", "rawpuppimet_px", 200, 25, 25);
	rawpuppimet_px->GetXaxis()->SetTitle("px");
	rawpuppimet_px->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimet_py = new TH1F("rawpuppimet_py", "rawpuppimet_py", 200, 25, 25);
	rawpuppimet_py->GetXaxis()->SetTitle("py");
	rawpuppimet_py->GetYaxis()->SetTitle("Entries");

    TH1F *uperpen_rawpuppi = new TH1F("uperpen_rawpuppi", "uperpen_rawpuppi", 200, 25, 25);
	uperpen_rawpuppi->GetXaxis()->SetTitle("u#perp");
	uperpen_rawpuppi->GetYaxis()->SetTitle("Entries");

    TH1F *uparall_rawpuppi = new TH1F("uparall_rawpuppi", "uparall_rawpuppi", 200, 25, 25);
	uparall_rawpuppi->GetXaxis()->SetTitle("u#parallel");
	uparall_rawpuppi->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale = new TH1F("rawpuppi_abs_scale", "rawpuppi_abs_scale", 200, 25, 25);
	rawpuppi_abs_scale->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_paraqt = new TH1F("rawpuppi_paraqt", "rawpuppi_paraqt", 200, 25, 25);
	rawpuppi_paraqt->GetXaxis()->SetTitle("u#parallel + q_{T}");
	rawpuppi_paraqt->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppimet_tightphopt = new TH1F("rawpuppimet_tightphopt", "rawpuppimet_tightphopt", 200, 0, 400);
	rawpuppimet_tightphopt->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpuppimet_tightphopt->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppipho_pxdist = new TH1F("rawpuppipho_pxdist", "rawpuppipho_pxdist", 200, 25, 25);
	rawpuppipho_pxdist->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpuppipho_pxdist->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppipho_pydist = new TH1F("rawpuppipho_pydist", "rawpuppipho_pydist", 200, 25, 25);
	rawpuppipho_pydist->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpuppipho_pydist->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_response = new TH1F("rawpuppi_response", "rawpuppi_response", 200, 25, 25);
	rawpuppi_response->GetXaxis()->SetTitle("q_{T} (GeV)");
	rawpuppi_response->GetYaxis()->SetTitle("- <u#parallel> / <q_{T}>");


    TH1F *rawpuppimet_nGoodvtx = new TH1F("rawpuppimet_nvtx", "rawpuppimet_nvtx", 200, 25, 25);
	rawpuppimet_nGoodvtx->GetXaxis()->SetTitle("nvtx");
	rawpuppimet_nGoodvtx->GetYaxis()->SetTitle("Entries");




    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */

    //rawpuppimet absolute scale with different bins
    TH1F *rawpuppi_abs_scale5060 = new TH1F("rawpuppi_abs_scale5060", "rawpuppi_abs_scale5060", 200, 25, 25);
	rawpuppi_abs_scale5060->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale5060->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale8090 = new TH1F("rawpuppi_abs_scale8090", "rawpuppi_abs_scale8090", 200, 25, 25);
	rawpuppi_abs_scale8090->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale8090->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale100110 = new TH1F("rawpuppi_abs_scale100110", "rawpuppi_abs_scale100110", 200, 25, 25);
	rawpuppi_abs_scale100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale120130 = new TH1F("rawpuppi_abs_scale120130", "rawpuppi_abs_scale120130", 200, 25, 25);
	rawpuppi_abs_scale120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale140150 = new TH1F("rawpuppi_abs_scale140150", "rawpuppi_abs_scale140150", 200, 25, 25);
	rawpuppi_abs_scale140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale160170 = new TH1F("rawpuppi_abs_scale160170", "rawpuppi_abs_scale160170", 200, 25, 25);
	rawpuppi_abs_scale160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale190200 = new TH1F("rawpuppi_abs_scale190200", "rawpuppi_abs_scale190200", 200, 25, 25);
	rawpuppi_abs_scale190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale190200->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale220230 = new TH1F("rawpuppi_abs_scale220230", "rawpuppi_abs_scale220230", 200, 25, 25);
	rawpuppi_abs_scale220230->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale220230->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale250260 = new TH1F("rawpuppi_abs_scale250260", "rawpuppi_abs_scale250260", 200, 25, 25);
	rawpuppi_abs_scale250260->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale250260->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scaled260270 = new TH1F("rawpuppi_abs_scaled260270", "rawpuppi_abs_scaled260270", 200, 25, 25);
	rawpuppi_abs_scaled260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scaled260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale280290 = new TH1F("rawpuppi_abs_scale280290", "rawpuppi_abs_scale280290", 200, 25, 25);
	rawpuppi_abs_scale280290->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale280290->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale300310 = new TH1F("rawpuppi_abs_scale300310", "rawpuppi_abs_scale300310", 200, 25, 25);
	rawpuppi_abs_scale300310->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale300310->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale310320 = new TH1F("rawpuppi_abs_scale310320", "rawpuppi_abs_scale310320", 200, 25, 25);
	rawpuppi_abs_scale310320->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale310320->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale320330 = new TH1F("rawpuppi_abs_scale320330", "rawpuppi_abs_scale320330", 200, 25, 25);
	rawpuppi_abs_scale320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale350360 = new TH1F("rawpuppi_abs_scale350360", "rawpuppi_abs_scale350360", 200, 25, 25);
	rawpuppi_abs_scale350360->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale350360->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale380390 = new TH1F("rawpuppi_abs_scale380390", "rawpuppi_abs_scale380390", 200, 25, 25);
	rawpuppi_abs_scale380390->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale380390->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale390400 = new TH1F("rawpuppi_abs_scale390400", "rawpuppi_abs_scale390400", 200, 25, 25);
	rawpuppi_abs_scale390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale390400->GetYaxis()->SetTitle("Entries");


    //rawpuppimet u_para with different bins
    TH1F *rawpuppi_upara5060 = new TH1F("rawpuppi_upara5060", "rawpuppi_upara5060", 200, 25, 25);
	rawpuppi_upara5060->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara5060->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara8090 = new TH1F("rawpuppi_upara8090", "rawpuppi_upara8090", 200, 25, 25);
	rawpuppi_upara8090->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara8090->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara100110 = new TH1F("rawpuppi_upara100110", "rawpuppi_upara100110", 200, 25, 25);
	rawpuppi_upara100110->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara120130 = new TH1F("rawpuppi_upara120130", "rawpuppi_upara120130", 200, 25, 25);
	rawpuppi_upara120130->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara140150 = new TH1F("rawpuppi_upara140150", "rawpuppi_upara140150", 200, 25, 25);
	rawpuppi_upara140150->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara160170 = new TH1F("rawpuppi_upara160170", "rawpuppi_upara160170", 200, 25, 25);
	rawpuppi_upara160170->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara190200 = new TH1F("rawpuppi_upara190200", "rawpuppi_upara190200", 200, 25, 25);
	rawpuppi_upara190200->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara190200->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara220230 = new TH1F("rawpuppi_upara220230", "rawpuppi_upara220230", 200, 25, 25);
	rawpuppi_upara220230->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara220230->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara250260 = new TH1F("rawpuppi_upara250260", "rawpuppi_upara250260", 200, 25, 25);
	rawpuppi_upara250260->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara250260->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara260270 = new TH1F("rawpuppi_upara260270", "rawpuppi_upara260270", 200, 25, 25);
	rawpuppi_upara260270->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara280290 = new TH1F("rawpuppi_upara280290", "rawpuppi_upara280290", 200, 25, 25);
	rawpuppi_upara280290->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara280290->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara300310 = new TH1F("rawpuppi_upara300310", "rawpuppi_upara300310", 200, 25, 25);
	rawpuppi_upara300310->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara300310->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara310320 = new TH1F("rawpuppi_upara310320", "rawpuppi_upara310320", 200, 25, 25);
	rawpuppi_upara310320->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara310320->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara320330 = new TH1F("rawpuppi_upara320330", "rawpuppi_upara320330", 200, 25, 25);
	rawpuppi_upara320330->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara350360 = new TH1F("rawpuppi_upara350360", "rawpuppi_upara350360", 200, 25, 25);
	rawpuppi_upara350360->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara350360->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara380390 = new TH1F("rawpuppi_upara380390", "rawpuppi_upara380390", 200, 25, 25);
	rawpuppi_upara380390->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara380390->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara390400 = new TH1F("rawpuppi_upara390400", "rawpuppi_upara390400", 200, 25, 25);
	rawpuppi_upara390400->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara390400->GetYaxis()->SetTitle("Entries");


    //rawpuppimet u_perpen with different bins
    TH1F *rawpuppi_uperpen5060 = new TH1F("rawpuppi_uperpen5060", "rawpuppi_uperpen5060", 200, 25, 25);
	rawpuppi_uperpen5060->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen5060->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen8090 = new TH1F("rawpuppi_uperpen8090", "rawpuppi_uperpen8090", 200, 25, 25);
	rawpuppi_uperpen8090->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen8090->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen100110 = new TH1F("rawpuppi_uperpen100110", "rawpuppi_uperpen100110", 200, 25, 25);
	rawpuppi_uperpen100110->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen120130 = new TH1F("rawpuppi_uperpen120130", "rawpuppi_uperpen120130", 200, 25, 25);
	rawpuppi_uperpen120130->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen140150 = new TH1F("rawpuppi_uperpen140150", "rawpuppi_uperpen140150", 200, 25, 25);
	rawpuppi_uperpen140150->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen160170 = new TH1F("rawpuppi_uperpen160170", "rawpuppi_uperpen160170", 200, 25, 25);
	rawpuppi_uperpen160170->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen190200 = new TH1F("rawpuppi_uperpen190200", "rawpuppi_uperpen190200", 200, 25, 25);
	rawpuppi_uperpen190200->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen190200->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen220230 = new TH1F("rawpuppi_uperpen220230", "rawpuppi_uperpen220230", 200, 25, 25);
	rawpuppi_uperpen220230->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen220230->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen250260 = new TH1F("rawpuppi_uperpen250260", "rawpuppi_uperpen250260", 200, 25, 25);
	rawpuppi_uperpen250260->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen250260->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen260270 = new TH1F("rawpuppi_uperpen260270", "rawpuppi_uperpen260270", 200, 25, 25);
	rawpuppi_uperpen260270->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen280290 = new TH1F("rawpuppi_uperpen280290", "rawpuppi_uperpen280290", 200, 25, 25);
	rawpuppi_uperpen280290->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen280290->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen300310 = new TH1F("rawpuppi_uperpen300310", "rawpuppi_uperpen300310", 200, 25, 25);
	rawpuppi_uperpen300310->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen300310->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen310320 = new TH1F("rawpuppi_uperpen310320", "rawpuppi_uperpen310320", 200, 25, 25);
	rawpuppi_uperpen310320->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen310320->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen320330 = new TH1F("rawpuppi_uperpen320330", "rawpuppi_uperpen320330", 200, 25, 25);
	rawpuppi_uperpen320330->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen350360 = new TH1F("rawpuppi_uperpen350360", "rawpuppi_uperpen350360", 200, 25, 25);
	rawpuppi_uperpen350360->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen350360->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen380390 = new TH1F("rawpuppi_uperpen380390", "rawpuppi_uperpen380390", 200, 25, 25);
	rawpuppi_uperpen380390->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen380390->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen390400 = new TH1F("rawpuppi_uperpen390400", "rawpuppi_uperpen390400", 200, 25, 25);
	rawpuppi_uperpen390400->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen390400->GetYaxis()->SetTitle("Entries");


    //rawpuppimet u_para with different nGoodvtx bins
    TH1F *rawpuppi_uparanvtx04 = new TH1F("rawpuppi_uparanvtx04", "rawpuppi_uparanvtx04", 200, 25, 25);
	rawpuppi_uparanvtx04->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx46 = new TH1F("rawpuppi_uparanvtx46", "rawpuppi_uparanvtx46", 200, 25, 25);
	rawpuppi_uparanvtx46->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx68 = new TH1F("rawpuppi_uparanvtx68", "rawpuppi_uparanvtx68", 200, 25, 25);
	rawpuppi_uparanvtx68->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx812 = new TH1F("rawpuppi_uparanvtx812", "rawpuppi_uparanvtx812", 200, 25, 25);
	rawpuppi_uparanvtx812->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx1218 = new TH1F("rawpuppi_uparanvtx1218", "rawpuppi_uparanvtx1218", 200, 25, 25);
	rawpuppi_uparanvtx1218->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx1824 = new TH1F("rawpuppi_uparanvtx1824", "rawpuppi_uparanvtx1824", 200, 25, 25);
	rawpuppi_uparanvtx1824->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx2428 = new TH1F("rawpuppi_uparanvtx2428", "rawpuppi_uparanvtx2428", 200, 25, 25);
	rawpuppi_uparanvtx2428->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx2832 = new TH1F("rawpuppi_uparanvtx2832", "rawpuppi_uparanvtx2832", 200, 25, 25);
	rawpuppi_uparanvtx2832->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx3236 = new TH1F("rawpuppi_uparanvtx3236", "rawpuppi_uparanvtx3236", 200, 25, 25);
	rawpuppi_uparanvtx3236->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparanvtx3640 = new TH1F("rawpuppi_uparanvtx3640", "rawpuppi_uparanvtx3640", 200, 25, 25);
	rawpuppi_uparanvtx3640->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_uparanvtx3640->GetYaxis()->SetTitle("Entries");

    

    //rawpuppimet u_perpen with different nGoodvtx bins
    TH1F *rawpuppi_uperpennvtx04 = new TH1F("rawpuppi_uperpennvtx04", "rawpuppi_uperpennvtx04", 200, 25, 25);
	rawpuppi_uperpennvtx04->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx46 = new TH1F("rawpuppi_uperpennvtx46", "rawpuppi_uperpennvtx46", 200, 25, 25);
	rawpuppi_uperpennvtx46->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx68 = new TH1F("rawpuppi_uperpennvtx68", "rawpuppi_uperpennvtx68", 200, 25, 25);
	rawpuppi_uperpennvtx68->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx812 = new TH1F("rawpuppi_uperpennvtx812", "rawpuppi_uperpennvtx812", 200, 25, 25);
	rawpuppi_uperpennvtx812->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx1218 = new TH1F("rawpuppi_uperpennvtx1218", "rawpuppi_uperpennvtx1218", 200, 25, 25);
	rawpuppi_uperpennvtx1218->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx1824 = new TH1F("rawpuppi_uperpennvtx1824", "rawpuppi_uperpennvtx1824", 200, 25, 25);
	rawpuppi_uperpennvtx1824->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx2428 = new TH1F("rawpuppi_uperpennvtx2428", "rawpuppi_uperpennvtx2428", 200, 25, 25);
	rawpuppi_uperpennvtx2428->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx2832 = new TH1F("rawpuppi_uperpennvtx2832", "rawpuppi_uperpennvtx2832", 200, 25, 25);
	rawpuppi_uperpennvtx2832->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx3236 = new TH1F("rawpuppi_uperpennvtx3236", "rawpuppi_uperpennvtx3236", 200, 25, 25);
	rawpuppi_uperpennvtx3236->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpennvtx3640 = new TH1F("rawpuppi_uperpennvtx3640", "rawpuppi_uperpennvtx3640", 200, 25, 25);
	rawpuppi_uperpennvtx3640->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpennvtx3640->GetYaxis()->SetTitle("Entries");


    //rawpuppimet u_paraqt with different bins
    TH1F *rawpuppi_uparaqt5060 = new TH1F("rawpuppi_uparaqt5060", "rawpuppi_uparaqt5060", 200, 25, 25);
	rawpuppi_uparaqt5060->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt5060->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt8090 = new TH1F("rawpuppi_uparaqt8090", "rawpuppi_uparaqt8090", 200, 25, 25);
	rawpuppi_uparaqt8090->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt8090->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt100110 = new TH1F("rawpuppi_uparaqt100110", "rawpuppi_uparaqt100110", 200, 25, 25);
	rawpuppi_uparaqt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt120130 = new TH1F("rawpuppi_uparaqt120130", "rawpuppi_uparaqt120130", 200, 25, 25);
	rawpuppi_uparaqt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt140150 = new TH1F("rawpuppi_uparaqt140150", "rawpuppi_uparaqt140150", 200, 25, 25);
	rawpuppi_uparaqt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt160170 = new TH1F("rawpuppi_uparaqt160170", "rawpuppi_uparaqt160170", 200, 25, 25);
	rawpuppi_uparaqt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt190200 = new TH1F("rawpuppi_uparaqt190200", "rawpuppi_uparaqt190200", 200, 25, 25);
	rawpuppi_uparaqt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt190200->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt220230 = new TH1F("rawpuppi_uparaqt220230", "rawpuppi_uparaqt220230", 200, 25, 25);
	rawpuppi_uparaqt220230->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt220230->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt250260 = new TH1F("rawpuppi_uparaqt250260", "rawpuppi_uparaqt250260", 200, 25, 25);
	rawpuppi_uparaqt250260->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt250260->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt260270 = new TH1F("rawpuppi_uparaqt260270", "rawpuppi_uparaqt260270", 200, 25, 25);
	rawpuppi_uparaqt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt280290 = new TH1F("rawpuppi_uparaqt280290", "rawpuppi_uparaqt280290", 200, 25, 25);
	rawpuppi_uparaqt280290->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt280290->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt300310 = new TH1F("rawpuppi_uparaqt300310", "rawpuppi_uparaqt300310", 200, 25, 25);
	rawpuppi_uparaqt300310->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt300310->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt310320 = new TH1F("rawpuppi_uparaqt310320", "rawpuppi_uparaqt310320", 200, 25, 25);
	rawpuppi_uparaqt310320->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt310320->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt320330 = new TH1F("rawpuppi_uparaqt320330", "rawpuppi_uparaqt320330", 200, 25, 25);
	rawpuppi_uparaqt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt350360 = new TH1F("rawpuppi_uparaqt350360", "rawpuppi_uparaqt350360", 200, 25, 25);
	rawpuppi_uparaqt350360->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt350360->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt380390 = new TH1F("rawpuppi_uparaqt380390", "rawpuppi_uparaqt380390", 200, 25, 25);
	rawpuppi_uparaqt380390->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt380390->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt390400 = new TH1F("rawpuppi_uparaqt390400", "rawpuppi_uparaqt390400", 200, 25, 25);
	rawpuppi_uparaqt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt390400->GetYaxis()->SetTitle("Entries");

    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */


    /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */

    //rawpuppi vs puppi pt > 100
    TH1F *rawpuppi_abs_scale_puppipt100110 = new TH1F("rawpuppi_abs_scale_puppipt100110", "rawpuppi_abs_scale_puppipt100110", 200, 25, 25);
	rawpuppi_abs_scale_puppipt100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt110120 = new TH1F("rawpuppi_abs_scale_puppipt110120", "rawpuppi_abs_scale_puppipt110120", 200, 25, 25);
	rawpuppi_abs_scale_puppipt110120->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt120130 = new TH1F("rawpuppi_abs_scale_puppipt120130", "rawpuppi_abs_scale_puppipt120130", 200, 25, 25);
	rawpuppi_abs_scale_puppipt120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt130140 = new TH1F("rawpuppi_abs_scale_puppipt130140", "rawpuppi_abs_scale_puppipt130140", 200, 25, 25);
	rawpuppi_abs_scale_puppipt130140->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt140150 = new TH1F("rawpuppi_abs_scale_puppipt140150", "rawpuppi_abs_scale_puppipt140150", 200, 25, 25);
	rawpuppi_abs_scale_puppipt140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt150160 = new TH1F("rawpuppi_abs_scale_puppipt150160", "rawpuppi_abs_scale_puppipt150160", 200, 25, 25);
	rawpuppi_abs_scale_puppipt150160->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt160170 = new TH1F("rawpuppi_abs_scale_puppipt160170", "rawpuppi_abs_scale_puppipt160170", 200, 25, 25);
	rawpuppi_abs_scale_puppipt160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt170180 = new TH1F("rawpuppi_abs_scale_puppipt170180", "rawpuppi_abs_scale_puppipt170180", 200, 25, 25);
	rawpuppi_abs_scale_puppipt170180->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt180190 = new TH1F("rawpuppi_abs_scale_puppipt180190", "rawpuppi_abs_scale_puppipt180190", 200, 25, 25);
	rawpuppi_abs_scale_puppipt180190->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt190200 = new TH1F("rawpuppi_abs_scale_puppipt190200", "rawpuppi_abs_scale_puppipt190200", 200, 25, 25);
	rawpuppi_abs_scale_puppipt190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawpuppi vs puppi pt > 200
    TH1F *rawpuppi_abs_scale_puppipt200210 = new TH1F("rawpuppi_abs_scale_puppipt200210", "rawpuppi_abs_scale_puppipt200210", 200, 25, 25);
	rawpuppi_abs_scale_puppipt200210->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt230240 = new TH1F("rawpuppi_abs_scale_puppipt230240", "rawpuppi_abs_scale_puppipt230240", 200, 25, 25);
	rawpuppi_abs_scale_puppipt230240->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt260270 = new TH1F("rawpuppi_abs_scale_puppipt260270", "rawpuppi_abs_scale_puppipt260270", 200, 25, 25);
	rawpuppi_abs_scale_puppipt260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt290300 = new TH1F("rawpuppi_abs_scale_puppipt290300", "rawpuppi_abs_scale_puppipt290300", 200, 25, 25);
	rawpuppi_abs_scale_puppipt290300->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt320330 = new TH1F("rawpuppi_abs_scale_puppipt320330", "rawpuppi_abs_scale_puppipt320330", 200, 25, 25);
	rawpuppi_abs_scale_puppipt320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt330340 = new TH1F("rawpuppi_abs_scale_puppipt330340", "rawpuppi_abs_scale_puppipt330340", 200, 25, 25);
	rawpuppi_abs_scale_puppipt330340->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt360370 = new TH1F("rawpuppi_abs_scale_puppipt360370", "rawpuppi_abs_scale_puppipt360370", 200, 25, 25);
	rawpuppi_abs_scale_puppipt360370->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt390400 = new TH1F("rawpuppi_abs_scale_puppipt390400", "rawpuppi_abs_scale_puppipt390400", 200, 25, 25);
	rawpuppi_abs_scale_puppipt390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_abs_scale_puppipt420430 = new TH1F("rawpuppi_abs_scale_puppipt420430", "rawpuppi_abs_scale_puppipt420430", 200, 25, 25);
	rawpuppi_abs_scale_puppipt420430->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	rawpuppi_abs_scale_puppipt420430->GetYaxis()->SetTitle("Entries");



    //rawpuppiupara vs puppi pt > 100
    TH1F *rawpuppi_upara_puppipt100110 = new TH1F("rawpuppi_upara_puppipt100110", "rawpuppi_upara_puppipt100110", 200, 25, 25);
	rawpuppi_upara_puppipt100110->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt110120 = new TH1F("rawpuppi_upara_puppipt110120", "rawpuppi_upara_puppipt110120", 200, 25, 25);
	rawpuppi_upara_puppipt110120->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt120130 = new TH1F("rawpuppi_upara_puppipt120130", "rawpuppi_upara_puppipt120130", 200, 25, 25);
	rawpuppi_upara_puppipt120130->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt130140 = new TH1F("rawpuppi_upara_puppipt130140", "rawpuppi_upara_puppipt130140", 200, 25, 25);
	rawpuppi_upara_puppipt130140->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt140150 = new TH1F("rawpuppi_upara_puppipt140150", "rawpuppi_upara_puppipt140150", 200, 25, 25);
	rawpuppi_upara_puppipt140150->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt150160 = new TH1F("rawpuppi_upara_puppipt150160", "rawpuppi_upara_puppipt150160", 200, 25, 25);
	rawpuppi_upara_puppipt150160->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt160170 = new TH1F("rawpuppi_upara_puppipt160170", "rawpuppi_upara_puppipt160170", 200, 25, 25);
	rawpuppi_upara_puppipt160170->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt170180 = new TH1F("rawpuppi_upara_puppipt170180", "rawpuppi_upara_puppipt170180", 200, 25, 25);
	rawpuppi_upara_puppipt170180->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt180190 = new TH1F("rawpuppi_upara_puppipt180190", "rawpuppi_upara_puppipt180190", 200, 25, 25);
	rawpuppi_upara_puppipt180190->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt190200 = new TH1F("rawpuppi_upara_puppipt190200", "rawpuppi_upara_puppipt190200", 200, 25, 25);
	rawpuppi_upara_puppipt190200->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawpuppiupara vs puppi pt > 200
    TH1F *rawpuppi_upara_puppipt200210 = new TH1F("rawpuppi_upara_puppipt200210", "rawpuppi_upara_puppipt200210", 200, 25, 25);
	rawpuppi_upara_puppipt200210->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt230240 = new TH1F("rawpuppi_upara_puppipt230240", "rawpuppi_upara_puppipt230240", 200, 25, 25);
	rawpuppi_upara_puppipt230240->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt260270 = new TH1F("rawpuppi_upara_puppipt260270", "rawpuppi_upara_puppipt260270", 200, 25, 25);
	rawpuppi_upara_puppipt260270->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt290300 = new TH1F("rawpuppi_upara_puppipt290300", "rawpuppi_upara_puppipt290300", 200, 25, 25);
	rawpuppi_upara_puppipt290300->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt320330 = new TH1F("rawpuppi_upara_puppipt320330", "rawpuppi_upara_puppipt320330", 200, 25, 25);
	rawpuppi_upara_puppipt320330->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt330340 = new TH1F("rawpuppi_upara_puppipt330340", "rawpuppi_upara_puppipt330340", 200, 25, 25);
	rawpuppi_upara_puppipt330340->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt360370 = new TH1F("rawpuppi_upara_puppipt360370", "rawpuppi_upara_puppipt360370", 200, 25, 25);
	rawpuppi_upara_puppipt360370->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt390400 = new TH1F("rawpuppi_upara_puppipt390400", "rawpuppi_upara_puppipt390400", 200, 25, 25);
	rawpuppi_upara_puppipt390400->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_upara_puppipt420430 = new TH1F("rawpuppi_upara_puppipt420430", "rawpuppi_upara_puppipt420430", 200, 25, 25);
	rawpuppi_upara_puppipt420430->GetXaxis()->SetTitle("u#parallel");
	rawpuppi_upara_puppipt420430->GetYaxis()->SetTitle("Entries");


    //rawpuppiuperpen vs puppi pt > 100
    TH1F *rawpuppi_uperpen_puppipt100110 = new TH1F("rawpuppi_uperpen_puppipt100110", "rawpuppi_uperpen_puppipt100110", 200, 25, 25);
	rawpuppi_uperpen_puppipt100110->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt110120 = new TH1F("rawpuppi_uperpen_puppipt110120", "rawpuppi_uperpen_puppipt110120", 200, 25, 25);
	rawpuppi_uperpen_puppipt110120->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt120130 = new TH1F("rawpuppi_uperpen_puppipt120130", "rawpuppi_uperpen_puppipt120130", 200, 25, 25);
	rawpuppi_uperpen_puppipt120130->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt130140 = new TH1F("rawpuppi_uperpen_puppipt130140", "rawpuppi_uperpen_puppipt130140", 200, 25, 25);
	rawpuppi_uperpen_puppipt130140->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt140150 = new TH1F("rawpuppi_uperpen_puppipt140150", "rawpuppi_uperpen_puppipt140150", 200, 25, 25);
	rawpuppi_uperpen_puppipt140150->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt150160 = new TH1F("rawpuppi_uperpen_puppipt150160", "rawpuppi_uperpen_puppipt150160", 200, 25, 25);
	rawpuppi_uperpen_puppipt150160->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt160170 = new TH1F("rawpuppi_uperpen_puppipt160170", "rawpuppi_uperpen_puppipt160170", 200, 25, 25);
	rawpuppi_uperpen_puppipt160170->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt170180 = new TH1F("rawpuppi_uperpen_puppipt170180", "rawpuppi_uperpen_puppipt170180", 200, 25, 25);
	rawpuppi_uperpen_puppipt170180->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt180190 = new TH1F("rawpuppi_uperpen_puppipt180190", "rawpuppi_uperpen_puppipt180190", 200, 25, 25);
	rawpuppi_uperpen_puppipt180190->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt190200 = new TH1F("rawpuppi_uperpen_puppipt190200", "rawpuppi_uperpen_puppipt190200", 200, 25, 25);
	rawpuppi_uperpen_puppipt190200->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawpuppiuperpen vs puppi pt > 200
    TH1F *rawpuppi_uperpen_puppipt200210 = new TH1F("rawpuppi_uperpen_puppipt200210", "rawpuppi_uperpen_puppipt200210", 200, 25, 25);
	rawpuppi_uperpen_puppipt200210->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt230240 = new TH1F("rawpuppi_uperpen_puppipt230240", "rawpuppi_uperpen_puppipt230240", 200, 25, 25);
	rawpuppi_uperpen_puppipt230240->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt260270 = new TH1F("rawpuppi_uperpen_puppipt260270", "rawpuppi_uperpen_puppipt260270", 200, 25, 25);
	rawpuppi_uperpen_puppipt260270->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt290300 = new TH1F("rawpuppi_uperpen_puppipt290300", "rawpuppi_uperpen_puppipt290300", 200, 25, 25);
	rawpuppi_uperpen_puppipt290300->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt320330 = new TH1F("rawpuppi_uperpen_puppipt320330", "rawpuppi_uperpen_puppipt320330", 200, 25, 25);
	rawpuppi_uperpen_puppipt320330->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt330340 = new TH1F("rawpuppi_uperpen_puppipt330340", "rawpuppi_uperpen_puppipt330340", 200, 25, 25);
	rawpuppi_uperpen_puppipt330340->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt360370 = new TH1F("rawpuppi_uperpen_puppipt360370", "rawpuppi_uperpen_puppipt360370", 200, 25, 25);
	rawpuppi_uperpen_puppipt360370->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt390400 = new TH1F("rawpuppi_uperpen_puppipt390400", "rawpuppi_uperpen_puppipt390400", 200, 25, 25);
	rawpuppi_uperpen_puppipt390400->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uperpen_puppipt420430 = new TH1F("rawpuppi_uperpen_puppipt420430", "rawpuppi_uperpen_puppipt420430", 200, 25, 25);
	rawpuppi_uperpen_puppipt420430->GetXaxis()->SetTitle("u#perp");
	rawpuppi_uperpen_puppipt420430->GetYaxis()->SetTitle("Entries");



    //rawpuppiuparaqt vs puppi pt > 100
    TH1F *rawpuppi_uparaqt_puppipt100110 = new TH1F("rawpuppi_uparaqt_puppipt100110", "rawpuppi_uparaqt_puppipt100110", 200, 25, 25);
	rawpuppi_uparaqt_puppipt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt110120 = new TH1F("rawpuppi_uparaqt_puppipt110120", "rawpuppi_uparaqt_puppipt110120", 200, 25, 25);
	rawpuppi_uparaqt_puppipt110120->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt120130 = new TH1F("rawpuppi_uparaqt_puppipt120130", "rawpuppi_uparaqt_puppipt120130", 200, 25, 25);
	rawpuppi_uparaqt_puppipt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt130140 = new TH1F("rawpuppi_uparaqt_puppipt130140", "rawpuppi_uparaqt_puppipt130140", 200, 25, 25);
	rawpuppi_uparaqt_puppipt130140->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt140150 = new TH1F("rawpuppi_uparaqt_puppipt140150", "rawpuppi_uparaqt_puppipt140150", 200, 25, 25);
	rawpuppi_uparaqt_puppipt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt150160 = new TH1F("rawpuppi_uparaqt_puppipt150160", "rawpuppi_uparaqt_puppipt150160", 200, 25, 25);
	rawpuppi_uparaqt_puppipt150160->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt160170 = new TH1F("rawpuppi_uparaqt_puppipt160170", "rawpuppi_uparaqt_puppipt160170", 200, 25, 25);
	rawpuppi_uparaqt_puppipt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt170180 = new TH1F("rawpuppi_uparaqt_puppipt170180", "rawpuppi_uparaqt_puppipt170180", 200, 25, 25);
	rawpuppi_uparaqt_puppipt170180->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt180190 = new TH1F("rawpuppi_uparaqt_puppipt180190", "rawpuppi_uparaqt_puppipt180190", 200, 25, 25);
	rawpuppi_uparaqt_puppipt180190->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt190200 = new TH1F("rawpuppi_uparaqt_puppipt190200", "rawpuppi_uparaqt_puppipt190200", 200, 25, 25);
	rawpuppi_uparaqt_puppipt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawpuppiuparaqt vs puppi pt > 200
    TH1F *rawpuppi_uparaqt_puppipt200210 = new TH1F("rawpuppi_uparaqt_puppipt200210", "rawpuppi_uparaqt_puppipt200210", 200, 25, 25);
	rawpuppi_uparaqt_puppipt200210->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt230240 = new TH1F("rawpuppi_uparaqt_puppipt230240", "rawpuppi_uparaqt_puppipt230240", 200, 25, 25);
	rawpuppi_uparaqt_puppipt230240->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt260270 = new TH1F("rawpuppi_uparaqt_puppipt260270", "rawpuppi_uparaqt_puppipt260270", 200, 25, 25);
	rawpuppi_uparaqt_puppipt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt290300 = new TH1F("rawpuppi_uparaqt_puppipt290300", "rawpuppi_uparaqt_puppipt290300", 200, 25, 25);
	rawpuppi_uparaqt_puppipt290300->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt320330 = new TH1F("rawpuppi_uparaqt_puppipt320330", "rawpuppi_uparaqt_puppipt320330", 200, 25, 25);
	rawpuppi_uparaqt_puppipt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt330340 = new TH1F("rawpuppi_uparaqt_puppipt330340", "rawpuppi_uparaqt_puppipt330340", 200, 25, 25);
	rawpuppi_uparaqt_puppipt330340->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt360370 = new TH1F("rawpuppi_uparaqt_puppipt360370", "rawpuppi_uparaqt_puppipt360370", 200, 25, 25);
	rawpuppi_uparaqt_puppipt360370->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt390400 = new TH1F("rawpuppi_uparaqt_puppipt390400", "rawpuppi_uparaqt_puppipt390400", 200, 25, 25);
	rawpuppi_uparaqt_puppipt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *rawpuppi_uparaqt_puppipt420430 = new TH1F("rawpuppi_uparaqt_puppipt420430", "rawpuppi_uparaqt_puppipt420430", 200, 25, 25);
	rawpuppi_uparaqt_puppipt420430->GetXaxis()->SetTitle("u#parallel+q_{T}");
	rawpuppi_uparaqt_puppipt420430->GetYaxis()->SetTitle("Entries");



     /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */


    //pfmet with eaxctly one tight photon, at least one jet, no loose electron or muon
    TH1F *pfmetphi = new TH1F("pfmetphi", "pfmetphi", 200, 25, 25);
	pfmetphi->GetXaxis()->SetTitle("phi");
	pfmetphi->GetYaxis()->SetTitle("Entries");

    TH1F *pfmetpt = new TH1F("pfmetpt", "pfmetpt", 200, 0, 400);
	pfmetpt->GetXaxis()->SetTitle("p_{T} (GeV)");
	pfmetpt->GetYaxis()->SetTitle("Entries");

    TH1F *pfmet_px_cal = new TH1F("pfmet_px_cal", "pfmet_px_cal", 200, 25, 25);
	pfmet_px_cal->GetXaxis()->SetTitle("px");
	pfmet_px_cal->GetYaxis()->SetTitle("Entries");

    TH1F *pfmet_py_cal = new TH1F("pfmet_py_cal", "pfmet_py_cal", 200, 25, 25);
	pfmet_py_cal->GetXaxis()->SetTitle("py");
	pfmet_py_cal->GetYaxis()->SetTitle("Entries");

    TH1F *pfmet_px = new TH1F("pfmet_px", "pfmet_px", 200, 25, 25);
	pfmet_px->GetXaxis()->SetTitle("px");
	pfmet_px->GetYaxis()->SetTitle("Entries");

    TH1F *pfmet_py = new TH1F("pfmet_py", "pfmet_py", 200, 25, 25);
	pfmet_py->GetXaxis()->SetTitle("py");
	pfmet_py->GetYaxis()->SetTitle("Entries");

    TH1F *uperpen_pf = new TH1F("uperpen_pf", "uperpen_pf", 200, 25, 25);
	uperpen_pf->GetXaxis()->SetTitle("u#perp");
	uperpen_pf->GetYaxis()->SetTitle("Entries");

    TH1F *uparall_pf = new TH1F("uparall_pf", "uparall_pf", 200, 25, 25);
	uparall_pf->GetXaxis()->SetTitle("u#parallel");
	uparall_pf->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale = new TH1F("pf_abs_scale", "pf_abs_scale", 200, 25, 25);
	pf_abs_scale->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale->GetYaxis()->SetTitle("Entries");

    TH1F *pf_paraqt = new TH1F("pf_paraqt", "pf_paraqt", 200, 25, 25);
	pf_paraqt->GetXaxis()->SetTitle("u#parallel + q_{T}");
	pf_paraqt->GetYaxis()->SetTitle("Entries");

    TH1F *pfmet_tightphopt = new TH1F("pfmet_tightphopt", "pfmet_tightphopt", 200, 0, 400);
	pfmet_tightphopt->GetXaxis()->SetTitle("p_{T} (GeV)");
	pfmet_tightphopt->GetYaxis()->SetTitle("Entries");

    TH1F *pfpho_pxdist = new TH1F("pfpho_pxdist", "pfpho_pxdist", 200, 25, 25);
	pfpho_pxdist->GetXaxis()->SetTitle("p_{T} (GeV)");
	pfpho_pxdist->GetYaxis()->SetTitle("Entries");

    TH1F *pfpho_pydist = new TH1F("pfpho_pydist", "pfpho_pydist", 200, 25, 25);
	pfpho_pydist->GetXaxis()->SetTitle("p_{T} (GeV)");
	pfpho_pydist->GetYaxis()->SetTitle("Entries");

    TH1F *pf_nGoodvtx = new TH1F("pf_nvtx", "pf_nvtx", 200, 25, 25);
	pf_nGoodvtx->GetXaxis()->SetTitle("nvtx");
	pf_nGoodvtx->GetYaxis()->SetTitle("Entries");

    TH1F *pf_response = new TH1F("pf_response", "pf_response", 200, 25, 25);
	pf_response->GetXaxis()->SetTitle("q_{T} (GeV)");
	pf_response->GetYaxis()->SetTitle("- <u#parallel> / <q_{T}>");


    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */


    //pfmet absolute scale with different bins
    TH1F *pf_abs_scale5060 = new TH1F("pf_abs_scale5060", "pf_abs_scale5060", 200, 25, 25);
	pf_abs_scale5060->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale5060->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale8090 = new TH1F("pf_abs_scale8090", "pf_abs_scale8090", 200, 25, 25);
	pf_abs_scale8090->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale8090->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale100110 = new TH1F("pf_abs_scale100110", "pf_abs_scale100110", 200, 25, 25);
	pf_abs_scale100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale120130 = new TH1F("pf_abs_scale120130", "pf_abs_scale120130", 200, 25, 25);
	pf_abs_scale120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale140150 = new TH1F("pf_abs_scale140150", "pf_abs_scale140150", 200, 25, 25);
	pf_abs_scale140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale160170 = new TH1F("pf_abs_scale160170", "pf_abs_scale160170", 200, 25, 25);
	pf_abs_scale160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale190200 = new TH1F("pf_abs_scale190200", "pf_abs_scale190200", 200, 25, 25);
	pf_abs_scale190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale190200->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale220230 = new TH1F("pf_abs_scale220230", "pf_abs_scale220230", 200, 25, 25);
	pf_abs_scale220230->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale220230->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale250260 = new TH1F("pf_abs_scale250260", "pf_abs_scale250260", 200, 25, 25);
	pf_abs_scale250260->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale250260->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scaled260270 = new TH1F("pf_abs_scaled260270", "pf_abs_scaled260270", 200, 25, 25);
	pf_abs_scaled260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scaled260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale280290 = new TH1F("pf_abs_scale280290", "pf_abs_scale280290", 200, 25, 25);
	pf_abs_scale280290->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale280290->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale300310 = new TH1F("pf_abs_scale300310", "pf_abs_scale300310", 200, 25, 25);
	pf_abs_scale300310->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale300310->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale310320 = new TH1F("pf_abs_scale310320", "pf_abs_scale310320", 200, 25, 25);
	pf_abs_scale310320->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale310320->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale320330 = new TH1F("pf_abs_scale320330", "pf_abs_scale320330", 200, 25, 25);
	pf_abs_scale320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale350360 = new TH1F("pf_abs_scale350360", "pf_abs_scale350360", 200, 25, 25);
	pf_abs_scale350360->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale350360->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale380390 = new TH1F("pf_abs_scale380390", "pf_abs_scale380390", 200, 25, 25);
	pf_abs_scale380390->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale380390->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale390400 = new TH1F("pf_abs_scale390400", "pf_abs_scale390400", 200, 25, 25);
	pf_abs_scale390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale390400->GetYaxis()->SetTitle("Entries");


    //pfmet u_para with different bins
    TH1F *pf_upara5060 = new TH1F("pf_upara5060", "pf_upara5060", 200, 25, 25);
	pf_upara5060->GetXaxis()->SetTitle("u#parallel");
	pf_upara5060->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara8090 = new TH1F("pf_upara8090", "pf_upara8090", 200, 25, 25);
	pf_upara8090->GetXaxis()->SetTitle("u#parallel");
	pf_upara8090->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara100110 = new TH1F("pf_upara100110", "pf_upara100110", 200, 25, 25);
	pf_upara100110->GetXaxis()->SetTitle("u#parallel");
	pf_upara100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara120130 = new TH1F("pf_upara120130", "pf_upara120130", 200, 25, 25);
	pf_upara120130->GetXaxis()->SetTitle("u#parallel");
	pf_upara120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara140150 = new TH1F("pf_upara140150", "pf_upara140150", 200, 25, 25);
	pf_upara140150->GetXaxis()->SetTitle("u#parallel");
	pf_upara140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara160170 = new TH1F("pf_upara160170", "pf_upara160170", 200, 25, 25);
	pf_upara160170->GetXaxis()->SetTitle("u#parallel");
	pf_upara160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara190200 = new TH1F("pf_upara190200", "pf_upara190200", 200, 25, 25);
	pf_upara190200->GetXaxis()->SetTitle("u#parallel");
	pf_upara190200->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara220230 = new TH1F("pf_upara220230", "pf_upara220230", 200, 25, 25);
	pf_upara220230->GetXaxis()->SetTitle("u#parallel");
	pf_upara220230->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara250260 = new TH1F("pf_upara250260", "pf_upara250260", 200, 25, 25);
	pf_upara250260->GetXaxis()->SetTitle("u#parallel");
	pf_upara250260->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara260270 = new TH1F("pf_upara260270", "pf_upara260270", 200, 25, 25);
	pf_upara260270->GetXaxis()->SetTitle("u#parallel");
	pf_upara260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara280290 = new TH1F("pf_upara280290", "pf_upara280290", 200, 25, 25);
	pf_upara280290->GetXaxis()->SetTitle("u#parallel");
	pf_upara280290->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara300310 = new TH1F("pf_upara300310", "pf_upara300310", 200, 25, 25);
	pf_upara300310->GetXaxis()->SetTitle("u#parallel");
	pf_upara300310->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara310320 = new TH1F("pf_upara310320", "pf_upara310320", 200, 25, 25);
	pf_upara310320->GetXaxis()->SetTitle("u#parallel");
	pf_upara310320->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara320330 = new TH1F("pf_upara320330", "pf_upara320330", 200, 25, 25);
	pf_upara320330->GetXaxis()->SetTitle("u#parallel");
	pf_upara320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara350360 = new TH1F("pf_upara350360", "pf_upara350360", 200, 25, 25);
	pf_upara350360->GetXaxis()->SetTitle("u#parallel");
	pf_upara350360->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara380390 = new TH1F("pf_upara380390", "pf_upara380390", 200, 25, 25);
	pf_upara380390->GetXaxis()->SetTitle("u#parallel");
	pf_upara380390->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara390400 = new TH1F("pf_upara390400", "pf_upara390400", 200, 25, 25);
	pf_upara390400->GetXaxis()->SetTitle("u#parallel");
	pf_upara390400->GetYaxis()->SetTitle("Entries");


    //pfmet u_perpen with different bins
    TH1F *pf_uperpen5060 = new TH1F("pf_uperpen5060", "pf_uperpen5060", 200, 25, 25);
	pf_uperpen5060->GetXaxis()->SetTitle("u#perp");
	pf_uperpen5060->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen8090 = new TH1F("pf_uperpen8090", "pf_uperpen8090", 200, 25, 25);
	pf_uperpen8090->GetXaxis()->SetTitle("u#perp");
	pf_uperpen8090->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen100110 = new TH1F("pf_uperpen100110", "pf_uperpen100110", 200, 25, 25);
	pf_uperpen100110->GetXaxis()->SetTitle("u#perp");
	pf_uperpen100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen120130 = new TH1F("pf_uperpen120130", "pf_uperpen120130", 200, 25, 25);
	pf_uperpen120130->GetXaxis()->SetTitle("u#perp");
	pf_uperpen120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen140150 = new TH1F("pf_uperpen140150", "pf_uperpen140150", 200, 25, 25);
	pf_uperpen140150->GetXaxis()->SetTitle("u#perp");
	pf_uperpen140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen160170 = new TH1F("pf_uperpen160170", "pf_uperpen160170", 200, 25, 25);
	pf_uperpen160170->GetXaxis()->SetTitle("u#perp");
	pf_uperpen160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen190200 = new TH1F("pf_uperpen190200", "pf_uperpen190200", 200, 25, 25);
	pf_uperpen190200->GetXaxis()->SetTitle("u#perp");
	pf_uperpen190200->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen220230 = new TH1F("pf_uperpen220230", "pf_uperpen220230", 200, 25, 25);
	pf_uperpen220230->GetXaxis()->SetTitle("u#perp");
	pf_uperpen220230->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen250260 = new TH1F("pf_uperpen250260", "pf_uperpen250260", 200, 25, 25);
	pf_uperpen250260->GetXaxis()->SetTitle("u#perp");
	pf_uperpen250260->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen260270 = new TH1F("pf_uperpen260270", "pf_uperpen260270", 200, 25, 25);
	pf_uperpen260270->GetXaxis()->SetTitle("u#perp");
	pf_uperpen260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen280290 = new TH1F("pf_uperpen280290", "pf_uperpen280290", 200, 25, 25);
	pf_uperpen280290->GetXaxis()->SetTitle("u#perp");
	pf_uperpen280290->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen300310 = new TH1F("pf_uperpen300310", "pf_uperpen300310", 200, 25, 25);
	pf_uperpen300310->GetXaxis()->SetTitle("u#perp");
	pf_uperpen300310->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen310320 = new TH1F("pf_uperpen310320", "pf_uperpen310320", 200, 25, 25);
	pf_uperpen310320->GetXaxis()->SetTitle("u#perp");
	pf_uperpen310320->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen320330 = new TH1F("pf_uperpen320330", "pf_uperpen320330", 200, 25, 25);
	pf_uperpen320330->GetXaxis()->SetTitle("u#perp");
	pf_uperpen320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen350360 = new TH1F("pf_uperpen350360", "pf_uperpen350360", 200, 25, 25);
	pf_uperpen350360->GetXaxis()->SetTitle("u#perp");
	pf_uperpen350360->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen380390 = new TH1F("pf_uperpen380390", "pf_uperpen380390", 200, 25, 25);
	pf_uperpen380390->GetXaxis()->SetTitle("u#perp");
	pf_uperpen380390->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen390400 = new TH1F("pf_uperpen390400", "pf_uperpen390400", 200, 25, 25);
	pf_uperpen390400->GetXaxis()->SetTitle("u#perp");
	pf_uperpen390400->GetYaxis()->SetTitle("Entries");


    //pfmet u_para with different nGoodvtx bins
    TH1F *pf_uparanvtx04 = new TH1F("pf_uparanvtx04", "pf_uparanvtx04", 200, 25, 25);
	pf_uparanvtx04->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx46 = new TH1F("pf_uparanvtx46", "pf_uparanvtx46", 200, 25, 25);
	pf_uparanvtx46->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx68 = new TH1F("pf_uparanvtx68", "pf_uparanvtx68", 200, 25, 25);
	pf_uparanvtx68->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx812 = new TH1F("pf_uparanvtx812", "pf_uparanvtx812", 200, 25, 25);
	pf_uparanvtx812->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx1218 = new TH1F("pf_uparanvtx1218", "pf_uparanvtx1218", 200, 25, 25);
	pf_uparanvtx1218->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx1824 = new TH1F("pf_uparanvtx1824", "pf_uparanvtx1824", 200, 25, 25);
	pf_uparanvtx1824->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx2428 = new TH1F("pf_uparanvtx2428", "pf_uparanvtx2428", 200, 25, 25);
	pf_uparanvtx2428->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx2832 = new TH1F("pf_uparanvtx2832", "pf_uparanvtx2832", 200, 25, 25);
	pf_uparanvtx2832->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx3236 = new TH1F("pf_uparanvtx3236", "pf_uparanvtx3236", 200, 25, 25);
	pf_uparanvtx3236->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparanvtx3640 = new TH1F("pf_uparanvtx3640", "pf_uparanvtx3640", 200, 25, 25);
	pf_uparanvtx3640->GetXaxis()->SetTitle("u#parallel");
	pf_uparanvtx3640->GetYaxis()->SetTitle("Entries");

    

    //pfmet u_perpen with different nGoodvtx bins
    TH1F *pf_uperpennvtx04 = new TH1F("pf_uperpennvtx04", "pf_uperpennvtx04", 200, 25, 25);
	pf_uperpennvtx04->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx46 = new TH1F("pf_uperpennvtx46", "pf_uperpennvtx46", 200, 25, 25);
	pf_uperpennvtx46->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx68 = new TH1F("pf_uperpennvtx68", "pf_uperpennvtx68", 200, 25, 25);
	pf_uperpennvtx68->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx812 = new TH1F("pf_uperpennvtx812", "pf_uperpennvtx812", 200, 25, 25);
	pf_uperpennvtx812->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx1218 = new TH1F("pf_uperpennvtx1218", "pf_uperpennvtx1218", 200, 25, 25);
	pf_uperpennvtx1218->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx1824 = new TH1F("pf_uperpennvtx1824", "pf_uperpennvtx1824", 200, 25, 25);
	pf_uperpennvtx1824->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx2428 = new TH1F("pf_uperpennvtx2428", "pf_uperpennvtx2428", 200, 25, 25);
	pf_uperpennvtx2428->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx2832 = new TH1F("pf_uperpennvtx2832", "pf_uperpennvtx2832", 200, 25, 25);
	pf_uperpennvtx2832->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx3236 = new TH1F("pf_uperpennvtx3236", "pf_uperpennvtx3236", 200, 25, 25);
	pf_uperpennvtx3236->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpennvtx3640 = new TH1F("pf_uperpennvtx3640", "pf_uperpennvtx3640", 200, 25, 25);
	pf_uperpennvtx3640->GetXaxis()->SetTitle("u#perp");
	pf_uperpennvtx3640->GetYaxis()->SetTitle("Entries");


    //pfmet u_paraqt with different bins
    TH1F *pf_uparaqt5060 = new TH1F("pf_uparaqt5060", "pf_uparaqt5060", 200, 25, 25);
	pf_uparaqt5060->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt5060->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt8090 = new TH1F("pf_uparaqt8090", "pf_uparaqt8090", 200, 25, 25);
	pf_uparaqt8090->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt8090->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt100110 = new TH1F("pf_uparaqt100110", "pf_uparaqt100110", 200, 25, 25);
	pf_uparaqt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt120130 = new TH1F("pf_uparaqt120130", "pf_uparaqt120130", 200, 25, 25);
	pf_uparaqt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt140150 = new TH1F("pf_uparaqt140150", "pf_uparaqt140150", 200, 25, 25);
	pf_uparaqt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt160170 = new TH1F("pf_uparaqt160170", "pf_uparaqt160170", 200, 25, 25);
	pf_uparaqt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt190200 = new TH1F("pf_uparaqt190200", "pf_uparaqt190200", 200, 25, 25);
	pf_uparaqt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt190200->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt220230 = new TH1F("pf_uparaqt220230", "pf_uparaqt220230", 200, 25, 25);
	pf_uparaqt220230->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt220230->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt250260 = new TH1F("pf_uparaqt250260", "pf_uparaqt250260", 200, 25, 25);
	pf_uparaqt250260->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt250260->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt260270 = new TH1F("pf_uparaqt260270", "pf_uparaqt260270", 200, 25, 25);
	pf_uparaqt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt280290 = new TH1F("pf_uparaqt280290", "pf_uparaqt280290", 200, 25, 25);
	pf_uparaqt280290->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt280290->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt300310 = new TH1F("pf_uparaqt300310", "pf_uparaqt300310", 200, 25, 25);
	pf_uparaqt300310->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt300310->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt310320 = new TH1F("pf_uparaqt310320", "pf_uparaqt310320", 200, 25, 25);
	pf_uparaqt310320->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt310320->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt320330 = new TH1F("pf_uparaqt320330", "pf_uparaqt320330", 200, 25, 25);
	pf_uparaqt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt350360 = new TH1F("pf_uparaqt350360", "pf_uparaqt350360", 200, 25, 25);
	pf_uparaqt350360->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt350360->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt380390 = new TH1F("pf_uparaqt380390", "pf_uparaqt380390", 200, 25, 25);
	pf_uparaqt380390->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt380390->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt390400 = new TH1F("pf_uparaqt390400", "pf_uparaqt390400", 200, 25, 25);
	pf_uparaqt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt390400->GetYaxis()->SetTitle("Entries");

    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */

   /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */

    //pf vs puppi pt > 100
    TH1F *pf_abs_scale_puppipt100110 = new TH1F("pf_abs_scale_puppipt100110", "pf_abs_scale_puppipt100110", 200, 25, 25);
	pf_abs_scale_puppipt100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt110120 = new TH1F("pf_abs_scale_puppipt110120", "pf_abs_scale_puppipt110120", 200, 25, 25);
	pf_abs_scale_puppipt110120->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt120130 = new TH1F("pf_abs_scale_puppipt120130", "pf_abs_scale_puppipt120130", 200, 25, 25);
	pf_abs_scale_puppipt120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt130140 = new TH1F("pf_abs_scale_puppipt130140", "pf_abs_scale_puppipt130140", 200, 25, 25);
	pf_abs_scale_puppipt130140->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt140150 = new TH1F("pf_abs_scale_puppipt140150", "pf_abs_scale_puppipt140150", 200, 25, 25);
	pf_abs_scale_puppipt140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt150160 = new TH1F("pf_abs_scale_puppipt150160", "pf_abs_scale_puppipt150160", 200, 25, 25);
	pf_abs_scale_puppipt150160->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt160170 = new TH1F("pf_abs_scale_puppipt160170", "pf_abs_scale_puppipt160170", 200, 25, 25);
	pf_abs_scale_puppipt160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt170180 = new TH1F("pf_abs_scale_puppipt170180", "pf_abs_scale_puppipt170180", 200, 25, 25);
	pf_abs_scale_puppipt170180->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt180190 = new TH1F("pf_abs_scale_puppipt180190", "pf_abs_scale_puppipt180190", 200, 25, 25);
	pf_abs_scale_puppipt180190->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt190200 = new TH1F("pf_abs_scale_puppipt190200", "pf_abs_scale_puppipt190200", 200, 25, 25);
	pf_abs_scale_puppipt190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt190200->GetYaxis()->SetTitle("Entries");


    //pf vs puppi pt > 200
    TH1F *pf_abs_scale_puppipt200210 = new TH1F("pf_abs_scale_puppipt200210", "pf_abs_scale_puppipt200210", 200, 25, 25);
	pf_abs_scale_puppipt200210->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt230240 = new TH1F("pf_abs_scale_puppipt230240", "pf_abs_scale_puppipt230240", 200, 25, 25);
	pf_abs_scale_puppipt230240->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt260270 = new TH1F("pf_abs_scale_puppipt260270", "pf_abs_scale_puppipt260270", 200, 25, 25);
	pf_abs_scale_puppipt260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt290300 = new TH1F("pf_abs_scale_puppipt290300", "pf_abs_scale_puppipt290300", 200, 25, 25);
	pf_abs_scale_puppipt290300->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt320330 = new TH1F("pf_abs_scale_puppipt320330", "pf_abs_scale_puppipt320330", 200, 25, 25);
	pf_abs_scale_puppipt320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt330340 = new TH1F("pf_abs_scale_puppipt330340", "pf_abs_scale_puppipt330340", 200, 25, 25);
	pf_abs_scale_puppipt330340->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt360370 = new TH1F("pf_abs_scale_puppipt360370", "pf_abs_scale_puppipt360370", 200, 25, 25);
	pf_abs_scale_puppipt360370->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt390400 = new TH1F("pf_abs_scale_puppipt390400", "pf_abs_scale_puppipt390400", 200, 25, 25);
	pf_abs_scale_puppipt390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *pf_abs_scale_puppipt420430 = new TH1F("pf_abs_scale_puppipt420430", "pf_abs_scale_puppipt420430", 200, 25, 25);
	pf_abs_scale_puppipt420430->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	pf_abs_scale_puppipt420430->GetYaxis()->SetTitle("Entries");


    //pfupara vs puppi pt > 100
    TH1F *pf_upara_puppipt100110 = new TH1F("pf_upara_puppipt100110", "pf_upara_puppipt100110", 200, 25, 25);
	pf_upara_puppipt100110->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt110120 = new TH1F("pf_upara_puppipt110120", "pf_upara_puppipt110120", 200, 25, 25);
	pf_upara_puppipt110120->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt120130 = new TH1F("pf_upara_puppipt120130", "pf_upara_puppipt120130", 200, 25, 25);
	pf_upara_puppipt120130->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt130140 = new TH1F("pf_upara_puppipt130140", "pf_upara_puppipt130140", 200, 25, 25);
	pf_upara_puppipt130140->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt140150 = new TH1F("pf_upara_puppipt140150", "pf_upara_puppipt140150", 200, 25, 25);
	pf_upara_puppipt140150->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt150160 = new TH1F("pf_upara_puppipt150160", "pf_upara_puppipt150160", 200, 25, 25);
	pf_upara_puppipt150160->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt160170 = new TH1F("pf_upara_puppipt160170", "pf_upara_puppipt160170", 200, 25, 25);
	pf_upara_puppipt160170->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt170180 = new TH1F("pf_upara_puppipt170180", "pf_upara_puppipt170180", 200, 25, 25);
	pf_upara_puppipt170180->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt180190 = new TH1F("pf_upara_puppipt180190", "pf_upara_puppipt180190", 200, 25, 25);
	pf_upara_puppipt180190->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt190200 = new TH1F("pf_upara_puppipt190200", "pf_upara_puppipt190200", 200, 25, 25);
	pf_upara_puppipt190200->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt190200->GetYaxis()->SetTitle("Entries");


    //pfupara vs puppi pt > 200
    TH1F *pf_upara_puppipt200210 = new TH1F("pf_upara_puppipt200210", "pf_upara_puppipt200210", 200, 25, 25);
	pf_upara_puppipt200210->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt230240 = new TH1F("pf_upara_puppipt230240", "pf_upara_puppipt230240", 200, 25, 25);
	pf_upara_puppipt230240->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt260270 = new TH1F("pf_upara_puppipt260270", "pf_upara_puppipt260270", 200, 25, 25);
	pf_upara_puppipt260270->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt290300 = new TH1F("pf_upara_puppipt290300", "pf_upara_puppipt290300", 200, 25, 25);
	pf_upara_puppipt290300->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt320330 = new TH1F("pf_upara_puppipt320330", "pf_upara_puppipt320330", 200, 25, 25);
	pf_upara_puppipt320330->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt330340 = new TH1F("pf_upara_puppipt330340", "pf_upara_puppipt330340", 200, 25, 25);
	pf_upara_puppipt330340->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt360370 = new TH1F("pf_upara_puppipt360370", "pf_upara_puppipt360370", 200, 25, 25);
	pf_upara_puppipt360370->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt390400 = new TH1F("pf_upara_puppipt390400", "pf_upara_puppipt390400", 200, 25, 25);
	pf_upara_puppipt390400->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *pf_upara_puppipt420430 = new TH1F("pf_upara_puppipt420430", "pf_upara_puppipt420430", 200, 25, 25);
	pf_upara_puppipt420430->GetXaxis()->SetTitle("u#parallel");
	pf_upara_puppipt420430->GetYaxis()->SetTitle("Entries");


    //pfuperpen vs puppi pt > 100
    TH1F *pf_uperpen_puppipt100110 = new TH1F("pf_uperpen_puppipt100110", "pf_uperpen_puppipt100110", 200, 25, 25);
	pf_uperpen_puppipt100110->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt110120 = new TH1F("pf_uperpen_puppipt110120", "pf_uperpen_puppipt110120", 200, 25, 25);
	pf_uperpen_puppipt110120->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt120130 = new TH1F("pf_uperpen_puppipt120130", "pf_uperpen_puppipt120130", 200, 25, 25);
	pf_uperpen_puppipt120130->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt130140 = new TH1F("pf_uperpen_puppipt130140", "pf_uperpen_puppipt130140", 200, 25, 25);
	pf_uperpen_puppipt130140->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt140150 = new TH1F("pf_uperpen_puppipt140150", "pf_uperpen_puppipt140150", 200, 25, 25);
	pf_uperpen_puppipt140150->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt150160 = new TH1F("pf_uperpen_puppipt150160", "pf_uperpen_puppipt150160", 200, 25, 25);
	pf_uperpen_puppipt150160->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt160170 = new TH1F("pf_uperpen_puppipt160170", "pf_uperpen_puppipt160170", 200, 25, 25);
	pf_uperpen_puppipt160170->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt170180 = new TH1F("pf_uperpen_puppipt170180", "pf_uperpen_puppipt170180", 200, 25, 25);
	pf_uperpen_puppipt170180->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt180190 = new TH1F("pf_uperpen_puppipt180190", "pf_uperpen_puppipt180190", 200, 25, 25);
	pf_uperpen_puppipt180190->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt190200 = new TH1F("pf_uperpen_puppipt190200", "pf_uperpen_puppipt190200", 200, 25, 25);
	pf_uperpen_puppipt190200->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt190200->GetYaxis()->SetTitle("Entries");


    //pfuperpen vs puppi pt > 200
    TH1F *pf_uperpen_puppipt200210 = new TH1F("pf_uperpen_puppipt200210", "pf_uperpen_puppipt200210", 200, 25, 25);
	pf_uperpen_puppipt200210->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt230240 = new TH1F("pf_uperpen_puppipt230240", "pf_uperpen_puppipt230240", 200, 25, 25);
	pf_uperpen_puppipt230240->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt260270 = new TH1F("pf_uperpen_puppipt260270", "pf_uperpen_puppipt260270", 200, 25, 25);
	pf_uperpen_puppipt260270->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt290300 = new TH1F("pf_uperpen_puppipt290300", "pf_uperpen_puppipt290300", 200, 25, 25);
	pf_uperpen_puppipt290300->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt320330 = new TH1F("pf_uperpen_puppipt320330", "pf_uperpen_puppipt320330", 200, 25, 25);
	pf_uperpen_puppipt320330->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt330340 = new TH1F("pf_uperpen_puppipt330340", "pf_uperpen_puppipt330340", 200, 25, 25);
	pf_uperpen_puppipt330340->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt360370 = new TH1F("pf_uperpen_puppipt360370", "pf_uperpen_puppipt360370", 200, 25, 25);
	pf_uperpen_puppipt360370->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt390400 = new TH1F("pf_uperpen_puppipt390400", "pf_uperpen_puppipt390400", 200, 25, 25);
	pf_uperpen_puppipt390400->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uperpen_puppipt420430 = new TH1F("pf_uperpen_puppipt420430", "pf_uperpen_puppipt420430", 200, 25, 25);
	pf_uperpen_puppipt420430->GetXaxis()->SetTitle("u#perp");
	pf_uperpen_puppipt420430->GetYaxis()->SetTitle("Entries");



    //pfuparaqt vs puppi pt > 100
    TH1F *pf_uparaqt_puppipt100110 = new TH1F("pf_uparaqt_puppipt100110", "pf_uparaqt_puppipt100110", 200, 25, 25);
	pf_uparaqt_puppipt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt110120 = new TH1F("pf_uparaqt_puppipt110120", "pf_uparaqt_puppipt110120", 200, 25, 25);
	pf_uparaqt_puppipt110120->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt120130 = new TH1F("pf_uparaqt_puppipt120130", "pf_uparaqt_puppipt120130", 200, 25, 25);
	pf_uparaqt_puppipt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt130140 = new TH1F("pf_uparaqt_puppipt130140", "pf_uparaqt_puppipt130140", 200, 25, 25);
	pf_uparaqt_puppipt130140->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt140150 = new TH1F("pf_uparaqt_puppipt140150", "pf_uparaqt_puppipt140150", 200, 25, 25);
	pf_uparaqt_puppipt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt150160 = new TH1F("pf_uparaqt_puppipt150160", "pf_uparaqt_puppipt150160", 200, 25, 25);
	pf_uparaqt_puppipt150160->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt160170 = new TH1F("pf_uparaqt_puppipt160170", "pf_uparaqt_puppipt160170", 200, 25, 25);
	pf_uparaqt_puppipt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt170180 = new TH1F("pf_uparaqt_puppipt170180", "pf_uparaqt_puppipt170180", 200, 25, 25);
	pf_uparaqt_puppipt170180->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt180190 = new TH1F("pf_uparaqt_puppipt180190", "pf_uparaqt_puppipt180190", 200, 25, 25);
	pf_uparaqt_puppipt180190->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt190200 = new TH1F("pf_uparaqt_puppipt190200", "pf_uparaqt_puppipt190200", 200, 25, 25);
	pf_uparaqt_puppipt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt190200->GetYaxis()->SetTitle("Entries");


    //pfuparaqt vs puppi pt > 200
    TH1F *pf_uparaqt_puppipt200210 = new TH1F("pf_uparaqt_puppipt200210", "pf_uparaqt_puppipt200210", 200, 25, 25);
	pf_uparaqt_puppipt200210->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt230240 = new TH1F("pf_uparaqt_puppipt230240", "pf_uparaqt_puppipt230240", 200, 25, 25);
	pf_uparaqt_puppipt230240->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt260270 = new TH1F("pf_uparaqt_puppipt260270", "pf_uparaqt_puppipt260270", 200, 25, 25);
	pf_uparaqt_puppipt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt290300 = new TH1F("pf_uparaqt_puppipt290300", "pf_uparaqt_puppipt290300", 200, 25, 25);
	pf_uparaqt_puppipt290300->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt320330 = new TH1F("pf_uparaqt_puppipt320330", "pf_uparaqt_puppipt320330", 200, 25, 25);
	pf_uparaqt_puppipt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt330340 = new TH1F("pf_uparaqt_puppipt330340", "pf_uparaqt_puppipt330340", 200, 25, 25);
	pf_uparaqt_puppipt330340->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt360370 = new TH1F("pf_uparaqt_puppipt360370", "pf_uparaqt_puppipt360370", 200, 25, 25);
	pf_uparaqt_puppipt360370->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt390400 = new TH1F("pf_uparaqt_puppipt390400", "pf_uparaqt_puppipt390400", 200, 25, 25);
	pf_uparaqt_puppipt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *pf_uparaqt_puppipt420430 = new TH1F("pf_uparaqt_puppipt420430", "pf_uparaqt_puppipt420430", 200, 25, 25);
	pf_uparaqt_puppipt420430->GetXaxis()->SetTitle("u#parallel+q_{T}");
	pf_uparaqt_puppipt420430->GetYaxis()->SetTitle("Entries");


     /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */
    
 
    //rawmet with eaxctly one tight photon, at least one jet, no loose electron or muon
    TH1F *rawmetphi = new TH1F("rawmetphi", "rawmetphi", 200, 25, 25);
	rawmetphi->GetXaxis()->SetTitle("phi");
	rawmetphi->GetYaxis()->SetTitle("Entries");

    TH1F *rawmetpt = new TH1F("rawmetpt", "rawmetpt", 200, 0, 400);
	rawmetpt->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawmetpt->GetYaxis()->SetTitle("Entries");

    TH1F *rawmet_px_cal = new TH1F("rawmet_px_cal", "rawmet_px_cal", 200, 25, 25);
	rawmet_px_cal->GetXaxis()->SetTitle("px");
	rawmet_px_cal->GetYaxis()->SetTitle("Entries");

    TH1F *rawmet_py_cal = new TH1F("rawmet_py_cal", "rawmet_py_cal", 200, 25, 25);
	rawmet_py_cal->GetXaxis()->SetTitle("py");
	rawmet_py_cal->GetYaxis()->SetTitle("Entries");

    TH1F *rawmet_px = new TH1F("rawmet_px", "rawmet_px", 200, 25, 25);
	rawmet_px->GetXaxis()->SetTitle("px");
	rawmet_px->GetYaxis()->SetTitle("Entries");

    TH1F *rawmet_py = new TH1F("rawmet_py", "rawmet_py", 200, 25, 25);
	rawmet_py->GetXaxis()->SetTitle("py");
	rawmet_py->GetYaxis()->SetTitle("Entries");

    TH1F *uperpen_raw = new TH1F("uperpen_raw", "uperpen_raw", 200, 25, 25);
	uperpen_raw->GetXaxis()->SetTitle("u#perp");
	uperpen_raw->GetYaxis()->SetTitle("Entries");

    TH1F *uparall_raw = new TH1F("uparall_raw", "uparall_raw", 200, 25, 25);
	uparall_raw->GetXaxis()->SetTitle("u#parallel");
	uparall_raw->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale = new TH1F("raw_abs_scale", "raw_abs_scale", 200, 25, 25);
	raw_abs_scale->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale->GetYaxis()->SetTitle("Entries");


    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */

    //rawmet absolute scale with different bins
    TH1F *raw_abs_scale5060 = new TH1F("raw_abs_scale5060", "raw_abs_scale5060", 200, 25, 25);
	raw_abs_scale5060->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale5060->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale8090 = new TH1F("raw_abs_scale8090", "raw_abs_scale8090", 200, 25, 25);
	raw_abs_scale8090->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale8090->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale100110 = new TH1F("raw_abs_scale100110", "raw_abs_scale100110", 200, 25, 25);
	raw_abs_scale100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale120130 = new TH1F("raw_abs_scale120130", "raw_abs_scale120130", 200, 25, 25);
	raw_abs_scale120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale140150 = new TH1F("raw_abs_scale140150", "raw_abs_scale140150", 200, 25, 25);
	raw_abs_scale140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale160170 = new TH1F("raw_abs_scale160170", "raw_abs_scale160170", 200, 25, 25);
	raw_abs_scale160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale190200 = new TH1F("raw_abs_scale190200", "raw_abs_scale190200", 200, 25, 25);
	raw_abs_scale190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale190200->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale220230 = new TH1F("raw_abs_scale220230", "raw_abs_scale220230", 200, 25, 25);
	raw_abs_scale220230->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale220230->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale250260 = new TH1F("raw_abs_scale250260", "raw_abs_scale250260", 200, 25, 25);
	raw_abs_scale250260->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale250260->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scaled260270 = new TH1F("raw_abs_scaled260270", "raw_abs_scaled260270", 200, 25, 25);
	raw_abs_scaled260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scaled260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale280290 = new TH1F("raw_abs_scale280290", "raw_abs_scale280290", 200, 25, 25);
	raw_abs_scale280290->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale280290->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale300310 = new TH1F("raw_abs_scale300310", "raw_abs_scale300310", 200, 25, 25);
	raw_abs_scale300310->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale300310->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale310320 = new TH1F("raw_abs_scale310320", "raw_abs_scale310320", 200, 25, 25);
	raw_abs_scale310320->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale310320->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale320330 = new TH1F("raw_abs_scale320330", "raw_abs_scale320330", 200, 25, 25);
	raw_abs_scale320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale350360 = new TH1F("raw_abs_scale350360", "raw_abs_scale350360", 200, 25, 25);
	raw_abs_scale350360->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale350360->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale380390 = new TH1F("raw_abs_scale380390", "raw_abs_scale380390", 200, 25, 25);
	raw_abs_scale380390->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale380390->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale390400 = new TH1F("raw_abs_scale390400", "raw_abs_scale390400", 200, 25, 25);
	raw_abs_scale390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale390400->GetYaxis()->SetTitle("Entries");


    TH1F *raw_nGoodvtx = new TH1F("raw_nvtx", "raw_nvtx", 200, 25, 25);
	raw_nGoodvtx->GetXaxis()->SetTitle("nvtx");
	raw_nGoodvtx->GetYaxis()->SetTitle("Entries");

    //rawmet u_para with different bins
    TH1F *raw_upara5060 = new TH1F("raw_upara5060", "raw_upara5060", 200, 25, 25);
	raw_upara5060->GetXaxis()->SetTitle("u#parallel");
	raw_upara5060->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara8090 = new TH1F("raw_upara8090", "raw_upara8090", 200, 25, 25);
	raw_upara8090->GetXaxis()->SetTitle("u#parallel");
	raw_upara8090->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara100110 = new TH1F("raw_upara100110", "raw_upara100110", 200, 25, 25);
	raw_upara100110->GetXaxis()->SetTitle("u#parallel");
	raw_upara100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara120130 = new TH1F("raw_upara120130", "raw_upara120130", 200, 25, 25);
	raw_upara120130->GetXaxis()->SetTitle("u#parallel");
	raw_upara120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara140150 = new TH1F("raw_upara140150", "raw_upara140150", 200, 25, 25);
	raw_upara140150->GetXaxis()->SetTitle("u#parallel");
	raw_upara140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara160170 = new TH1F("raw_upara160170", "raw_upara160170", 200, 25, 25);
	raw_upara160170->GetXaxis()->SetTitle("u#parallel");
	raw_upara160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara190200 = new TH1F("raw_upara190200", "raw_upara190200", 200, 25, 25);
	raw_upara190200->GetXaxis()->SetTitle("u#parallel");
	raw_upara190200->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara220230 = new TH1F("raw_upara220230", "raw_upara220230", 200, 25, 25);
	raw_upara220230->GetXaxis()->SetTitle("u#parallel");
	raw_upara220230->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara250260 = new TH1F("raw_upara250260", "raw_upara250260", 200, 25, 25);
	raw_upara250260->GetXaxis()->SetTitle("u#parallel");
	raw_upara250260->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara260270 = new TH1F("raw_upara260270", "raw_upara260270", 200, 25, 25);
	raw_upara260270->GetXaxis()->SetTitle("u#parallel");
	raw_upara260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara280290 = new TH1F("raw_upara280290", "raw_upara280290", 200, 25, 25);
	raw_upara280290->GetXaxis()->SetTitle("u#parallel");
	raw_upara280290->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara300310 = new TH1F("raw_upara300310", "raw_upara300310", 200, 25, 25);
	raw_upara300310->GetXaxis()->SetTitle("u#parallel");
	raw_upara300310->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara310320 = new TH1F("raw_upara310320", "raw_upara310320", 200, 25, 25);
	raw_upara310320->GetXaxis()->SetTitle("u#parallel");
	raw_upara310320->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara320330 = new TH1F("raw_upara320330", "raw_upara320330", 200, 25, 25);
	raw_upara320330->GetXaxis()->SetTitle("u#parallel");
	raw_upara320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara350360 = new TH1F("raw_upara350360", "raw_upara350360", 200, 25, 25);
	raw_upara350360->GetXaxis()->SetTitle("u#parallel");
	raw_upara350360->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara380390 = new TH1F("raw_upara380390", "raw_upara380390", 200, 25, 25);
	raw_upara380390->GetXaxis()->SetTitle("u#parallel");
	raw_upara380390->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara390400 = new TH1F("raw_upara390400", "raw_upara390400", 200, 25, 25);
	raw_upara390400->GetXaxis()->SetTitle("u#parallel");
	raw_upara390400->GetYaxis()->SetTitle("Entries");


    //rawmet u_perpen with different bins
    TH1F *raw_uperpen5060 = new TH1F("raw_uperpen5060", "raw_uperpen5060", 200, 25, 25);
	raw_uperpen5060->GetXaxis()->SetTitle("u#perp");
	raw_uperpen5060->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen8090 = new TH1F("raw_uperpen8090", "raw_uperpen8090", 200, 25, 25);
	raw_uperpen8090->GetXaxis()->SetTitle("u#perp");
	raw_uperpen8090->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen100110 = new TH1F("raw_uperpen100110", "raw_uperpen100110", 200, 25, 25);
	raw_uperpen100110->GetXaxis()->SetTitle("u#perp");
	raw_uperpen100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen120130 = new TH1F("raw_uperpen120130", "raw_uperpen120130", 200, 25, 25);
	raw_uperpen120130->GetXaxis()->SetTitle("u#perp");
	raw_uperpen120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen140150 = new TH1F("raw_uperpen140150", "raw_uperpen140150", 200, 25, 25);
	raw_uperpen140150->GetXaxis()->SetTitle("u#perp");
	raw_uperpen140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen160170 = new TH1F("raw_uperpen160170", "raw_uperpen160170", 200, 25, 25);
	raw_uperpen160170->GetXaxis()->SetTitle("u#perp");
	raw_uperpen160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen190200 = new TH1F("raw_uperpen190200", "raw_uperpen190200", 200, 25, 25);
	raw_uperpen190200->GetXaxis()->SetTitle("u#perp");
	raw_uperpen190200->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen220230 = new TH1F("raw_uperpen220230", "raw_uperpen220230", 200, 25, 25);
	raw_uperpen220230->GetXaxis()->SetTitle("u#perp");
	raw_uperpen220230->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen250260 = new TH1F("raw_uperpen250260", "raw_uperpen250260", 200, 25, 25);
	raw_uperpen250260->GetXaxis()->SetTitle("u#perp");
	raw_uperpen250260->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen260270 = new TH1F("raw_uperpen260270", "raw_uperpen260270", 200, 25, 25);
	raw_uperpen260270->GetXaxis()->SetTitle("u#perp");
	raw_uperpen260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen280290 = new TH1F("raw_uperpen280290", "raw_uperpen280290", 200, 25, 25);
	raw_uperpen280290->GetXaxis()->SetTitle("u#perp");
	raw_uperpen280290->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen300310 = new TH1F("raw_uperpen300310", "raw_uperpen300310", 200, 25, 25);
	raw_uperpen300310->GetXaxis()->SetTitle("u#perp");
	raw_uperpen300310->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen310320 = new TH1F("raw_uperpen310320", "raw_uperpen310320", 200, 25, 25);
	raw_uperpen310320->GetXaxis()->SetTitle("u#perp");
	raw_uperpen310320->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen320330 = new TH1F("raw_uperpen320330", "raw_uperpen320330", 200, 25, 25);
	raw_uperpen320330->GetXaxis()->SetTitle("u#perp");
	raw_uperpen320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen350360 = new TH1F("raw_uperpen350360", "raw_uperpen350360", 200, 25, 25);
	raw_uperpen350360->GetXaxis()->SetTitle("u#perp");
	raw_uperpen350360->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen380390 = new TH1F("raw_uperpen380390", "raw_uperpen380390", 200, 25, 25);
	raw_uperpen380390->GetXaxis()->SetTitle("u#perp");
	raw_uperpen380390->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen390400 = new TH1F("raw_uperpen390400", "raw_uperpen390400", 200, 25, 25);
	raw_uperpen390400->GetXaxis()->SetTitle("u#perp");
	raw_uperpen390400->GetYaxis()->SetTitle("Entries");


    //rawmet u_para with different nGoodvtx bins
    TH1F *raw_uparanvtx04 = new TH1F("raw_uparanvtx04", "raw_uparanvtx04", 200, 25, 25);
	raw_uparanvtx04->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx46 = new TH1F("raw_uparanvtx46", "raw_uparanvtx46", 200, 25, 25);
	raw_uparanvtx46->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx68 = new TH1F("raw_uparanvtx68", "raw_uparanvtx68", 200, 25, 25);
	raw_uparanvtx68->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx812 = new TH1F("raw_uparanvtx812", "raw_uparanvtx812", 200, 25, 25);
	raw_uparanvtx812->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx1218 = new TH1F("raw_uparanvtx1218", "raw_uparanvtx1218", 200, 25, 25);
	raw_uparanvtx1218->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx1824 = new TH1F("raw_uparanvtx1824", "raw_uparanvtx1824", 200, 25, 25);
	raw_uparanvtx1824->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx2428 = new TH1F("raw_uparanvtx2428", "raw_uparanvtx2428", 200, 25, 25);
	raw_uparanvtx2428->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx2832 = new TH1F("raw_uparanvtx2832", "raw_uparanvtx2832", 200, 25, 25);
	raw_uparanvtx2832->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx3236 = new TH1F("raw_uparanvtx3236", "raw_uparanvtx3236", 200, 25, 25);
	raw_uparanvtx3236->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparanvtx3640 = new TH1F("raw_uparanvtx3640", "raw_uparanvtx3640", 200, 25, 25);
	raw_uparanvtx3640->GetXaxis()->SetTitle("u#parallel");
	raw_uparanvtx3640->GetYaxis()->SetTitle("Entries");

    

    //rawmet u_perpen with different nGoodvtx bins
    TH1F *raw_uperpennvtx04 = new TH1F("raw_uperpennvtx04", "raw_uperpennvtx04", 200, 25, 25);
	raw_uperpennvtx04->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx04->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx46 = new TH1F("raw_uperpennvtx46", "raw_uperpennvtx46", 200, 25, 25);
	raw_uperpennvtx46->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx46->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx68 = new TH1F("raw_uperpennvtx68", "raw_uperpennvtx68", 200, 25, 25);
	raw_uperpennvtx68->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx68->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx812 = new TH1F("raw_uperpennvtx812", "raw_uperpennvtx812", 200, 25, 25);
	raw_uperpennvtx812->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx812->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx1218 = new TH1F("raw_uperpennvtx1218", "raw_uperpennvtx1218", 200, 25, 25);
	raw_uperpennvtx1218->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx1218->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx1824 = new TH1F("raw_uperpennvtx1824", "raw_uperpennvtx1824", 200, 25, 25);
	raw_uperpennvtx1824->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx1824->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx2428 = new TH1F("raw_uperpennvtx2428", "raw_uperpennvtx2428", 200, 25, 25);
	raw_uperpennvtx2428->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx2428->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx2832 = new TH1F("raw_uperpennvtx2832", "raw_uperpennvtx2832", 200, 25, 25);
	raw_uperpennvtx2832->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx2832->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx3236 = new TH1F("raw_uperpennvtx3236", "raw_uperpennvtx3236", 200, 25, 25);
	raw_uperpennvtx3236->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx3236->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpennvtx3640 = new TH1F("raw_uperpennvtx3640", "raw_uperpennvtx3640", 200, 25, 25);
	raw_uperpennvtx3640->GetXaxis()->SetTitle("u#perp");
	raw_uperpennvtx3640->GetYaxis()->SetTitle("Entries");



    //rawmet u_paraqt with different bins
    TH1F *raw_uparaqt5060 = new TH1F("raw_uparaqt5060", "raw_uparaqt5060", 200, 25, 25);
	raw_uparaqt5060->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt5060->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt8090 = new TH1F("raw_uparaqt8090", "raw_uparaqt8090", 200, 25, 25);
	raw_uparaqt8090->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt8090->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt100110 = new TH1F("raw_uparaqt100110", "raw_uparaqt100110", 200, 25, 25);
	raw_uparaqt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt120130 = new TH1F("raw_uparaqt120130", "raw_uparaqt120130", 200, 25, 25);
	raw_uparaqt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt140150 = new TH1F("raw_uparaqt140150", "raw_uparaqt140150", 200, 25, 25);
	raw_uparaqt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt160170 = new TH1F("raw_uparaqt160170", "raw_uparaqt160170", 200, 25, 25);
	raw_uparaqt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt190200 = new TH1F("raw_uparaqt190200", "raw_uparaqt190200", 200, 25, 25);
	raw_uparaqt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt190200->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt220230 = new TH1F("raw_uparaqt220230", "raw_uparaqt220230", 200, 25, 25);
	raw_uparaqt220230->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt220230->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt250260 = new TH1F("raw_uparaqt250260", "raw_uparaqt250260", 200, 25, 25);
	raw_uparaqt250260->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt250260->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt260270 = new TH1F("raw_uparaqt260270", "raw_uparaqt260270", 200, 25, 25);
	raw_uparaqt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt280290 = new TH1F("raw_uparaqt280290", "raw_uparaqt280290", 200, 25, 25);
	raw_uparaqt280290->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt280290->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt300310 = new TH1F("raw_uparaqt300310", "raw_uparaqt300310", 200, 25, 25);
	raw_uparaqt300310->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt300310->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt310320 = new TH1F("raw_uparaqt310320", "raw_uparaqt310320", 200, 25, 25);
	raw_uparaqt310320->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt310320->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt320330 = new TH1F("raw_uparaqt320330", "raw_uparaqt320330", 200, 25, 25);
	raw_uparaqt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt350360 = new TH1F("raw_uparaqt350360", "raw_uparaqt350360", 200, 25, 25);
	raw_uparaqt350360->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt350360->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt380390 = new TH1F("raw_uparaqt380390", "raw_uparaqt380390", 200, 25, 25);
	raw_uparaqt380390->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt380390->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt390400 = new TH1F("raw_uparaqt390400", "raw_uparaqt390400", 200, 25, 25);
	raw_uparaqt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt390400->GetYaxis()->SetTitle("Entries");

    /*
    =================--------======--======--======---------=================
    =================--====--======--======--======--=====--=================
    =================--====--======----------======--=====--=================
    =================--------======----------======--=====--=================
    =================--============--======--======--=====--=================
    =================--============--======--======---------=================
    */

   /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */

    //raw vs puppi pt > 100
    TH1F *raw_abs_scale_puppipt100110 = new TH1F("raw_abs_scale_puppipt100110", "raw_abs_scale_puppipt100110", 200, 25, 25);
	raw_abs_scale_puppipt100110->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt110120 = new TH1F("raw_abs_scale_puppipt110120", "raw_abs_scale_puppipt110120", 200, 25, 25);
	raw_abs_scale_puppipt110120->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt120130 = new TH1F("raw_abs_scale_puppipt120130", "raw_abs_scale_puppipt120130", 200, 25, 25);
	raw_abs_scale_puppipt120130->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt130140 = new TH1F("raw_abs_scale_puppipt130140", "raw_abs_scale_puppipt130140", 200, 25, 25);
	raw_abs_scale_puppipt130140->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt140150 = new TH1F("raw_abs_scale_puppipt140150", "raw_abs_scale_puppipt140150", 200, 25, 25);
	raw_abs_scale_puppipt140150->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt150160 = new TH1F("raw_abs_scale_puppipt150160", "raw_abs_scale_puppipt150160", 200, 25, 25);
	raw_abs_scale_puppipt150160->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt160170 = new TH1F("raw_abs_scale_puppipt160170", "raw_abs_scale_puppipt160170", 200, 25, 25);
	raw_abs_scale_puppipt160170->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt170180 = new TH1F("raw_abs_scale_puppipt170180", "raw_abs_scale_puppipt170180", 200, 25, 25);
	raw_abs_scale_puppipt170180->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt180190 = new TH1F("raw_abs_scale_puppipt180190", "raw_abs_scale_puppipt180190", 200, 25, 25);
	raw_abs_scale_puppipt180190->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt190200 = new TH1F("raw_abs_scale_puppipt190200", "raw_abs_scale_puppipt190200", 200, 25, 25);
	raw_abs_scale_puppipt190200->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt190200->GetYaxis()->SetTitle("Entries");


    //raw vs puppi pt > 200
    TH1F *raw_abs_scale_puppipt200210 = new TH1F("raw_abs_scale_puppipt200210", "raw_abs_scale_puppipt200210", 200, 25, 25);
	raw_abs_scale_puppipt200210->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt230240 = new TH1F("raw_abs_scale_puppipt230240", "raw_abs_scale_puppipt230240", 200, 25, 25);
	raw_abs_scale_puppipt230240->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt260270 = new TH1F("raw_abs_scale_puppipt260270", "raw_abs_scale_puppipt260270", 200, 25, 25);
	raw_abs_scale_puppipt260270->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt290300 = new TH1F("raw_abs_scale_puppipt290300", "raw_abs_scale_puppipt290300", 200, 25, 25);
	raw_abs_scale_puppipt290300->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt320330 = new TH1F("raw_abs_scale_puppipt320330", "raw_abs_scale_puppipt320330", 200, 25, 25);
	raw_abs_scale_puppipt320330->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt330340 = new TH1F("raw_abs_scale_puppipt330340", "raw_abs_scale_puppipt330340", 200, 25, 25);
	raw_abs_scale_puppipt330340->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt360370 = new TH1F("raw_abs_scale_puppipt360370", "raw_abs_scale_puppipt360370", 200, 25, 25);
	raw_abs_scale_puppipt360370->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt390400 = new TH1F("raw_abs_scale_puppipt390400", "raw_abs_scale_puppipt390400", 200, 25, 25);
	raw_abs_scale_puppipt390400->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *raw_abs_scale_puppipt420430 = new TH1F("raw_abs_scale_puppipt420430", "raw_abs_scale_puppipt420430", 200, 25, 25);
	raw_abs_scale_puppipt420430->GetXaxis()->SetTitle("- <u#parallel> / <q_{T}>");
	raw_abs_scale_puppipt420430->GetYaxis()->SetTitle("Entries");


    //rawupara vs puppi pt > 100
    TH1F *raw_upara_puppipt100110 = new TH1F("raw_upara_puppipt100110", "raw_upara_puppipt100110", 200, 25, 25);
	raw_upara_puppipt100110->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt110120 = new TH1F("raw_upara_puppipt110120", "raw_upara_puppipt110120", 200, 25, 25);
	raw_upara_puppipt110120->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt120130 = new TH1F("raw_upara_puppipt120130", "raw_upara_puppipt120130", 200, 25, 25);
	raw_upara_puppipt120130->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt130140 = new TH1F("raw_upara_puppipt130140", "raw_upara_puppipt130140", 200, 25, 25);
	raw_upara_puppipt130140->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt140150 = new TH1F("raw_upara_puppipt140150", "raw_upara_puppipt140150", 200, 25, 25);
	raw_upara_puppipt140150->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt150160 = new TH1F("raw_upara_puppipt150160", "raw_upara_puppipt150160", 200, 25, 25);
	raw_upara_puppipt150160->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt160170 = new TH1F("raw_upara_puppipt160170", "raw_upara_puppipt160170", 200, 25, 25);
	raw_upara_puppipt160170->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt170180 = new TH1F("raw_upara_puppipt170180", "raw_upara_puppipt170180", 200, 25, 25);
	raw_upara_puppipt170180->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt180190 = new TH1F("raw_upara_puppipt180190", "raw_upara_puppipt180190", 200, 25, 25);
	raw_upara_puppipt180190->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt190200 = new TH1F("raw_upara_puppipt190200", "raw_upara_puppipt190200", 200, 25, 25);
	raw_upara_puppipt190200->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawupara vs puppi pt > 200
    TH1F *raw_upara_puppipt200210 = new TH1F("raw_upara_puppipt200210", "raw_upara_puppipt200210", 200, 25, 25);
	raw_upara_puppipt200210->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt230240 = new TH1F("raw_upara_puppipt230240", "raw_upara_puppipt230240", 200, 25, 25);
	raw_upara_puppipt230240->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt260270 = new TH1F("raw_upara_puppipt260270", "raw_upara_puppipt260270", 200, 25, 25);
	raw_upara_puppipt260270->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt290300 = new TH1F("raw_upara_puppipt290300", "raw_upara_puppipt290300", 200, 25, 25);
	raw_upara_puppipt290300->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt320330 = new TH1F("raw_upara_puppipt320330", "raw_upara_puppipt320330", 200, 25, 25);
	raw_upara_puppipt320330->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt330340 = new TH1F("raw_upara_puppipt330340", "raw_upara_puppipt330340", 200, 25, 25);
	raw_upara_puppipt330340->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt360370 = new TH1F("raw_upara_puppipt360370", "raw_upara_puppipt360370", 200, 25, 25);
	raw_upara_puppipt360370->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt390400 = new TH1F("raw_upara_puppipt390400", "raw_upara_puppipt390400", 200, 25, 25);
	raw_upara_puppipt390400->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *raw_upara_puppipt420430 = new TH1F("raw_upara_puppipt420430", "raw_upara_puppipt420430", 200, 25, 25);
	raw_upara_puppipt420430->GetXaxis()->SetTitle("u#parallel");
	raw_upara_puppipt420430->GetYaxis()->SetTitle("Entries");


    //rawuperpen vs puppi pt > 100
    TH1F *raw_uperpen_puppipt100110 = new TH1F("raw_uperpen_puppipt100110", "raw_uperpen_puppipt100110", 200, 25, 25);
	raw_uperpen_puppipt100110->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt110120 = new TH1F("raw_uperpen_puppipt110120", "raw_uperpen_puppipt110120", 200, 25, 25);
	raw_uperpen_puppipt110120->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt120130 = new TH1F("raw_uperpen_puppipt120130", "raw_uperpen_puppipt120130", 200, 25, 25);
	raw_uperpen_puppipt120130->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt130140 = new TH1F("raw_uperpen_puppipt130140", "raw_uperpen_puppipt130140", 200, 25, 25);
	raw_uperpen_puppipt130140->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt140150 = new TH1F("raw_uperpen_puppipt140150", "raw_uperpen_puppipt140150", 200, 25, 25);
	raw_uperpen_puppipt140150->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt150160 = new TH1F("raw_uperpen_puppipt150160", "raw_uperpen_puppipt150160", 200, 25, 25);
	raw_uperpen_puppipt150160->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt160170 = new TH1F("raw_uperpen_puppipt160170", "raw_uperpen_puppipt160170", 200, 25, 25);
	raw_uperpen_puppipt160170->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt170180 = new TH1F("raw_uperpen_puppipt170180", "raw_uperpen_puppipt170180", 200, 25, 25);
	raw_uperpen_puppipt170180->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt180190 = new TH1F("raw_uperpen_puppipt180190", "raw_uperpen_puppipt180190", 200, 25, 25);
	raw_uperpen_puppipt180190->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt190200 = new TH1F("raw_uperpen_puppipt190200", "raw_uperpen_puppipt190200", 200, 25, 25);
	raw_uperpen_puppipt190200->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawuperpen vs puppi pt > 200
    TH1F *raw_uperpen_puppipt200210 = new TH1F("raw_uperpen_puppipt200210", "raw_uperpen_puppipt200210", 200, 25, 25);
	raw_uperpen_puppipt200210->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt230240 = new TH1F("raw_uperpen_puppipt230240", "raw_uperpen_puppipt230240", 200, 25, 25);
	raw_uperpen_puppipt230240->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt260270 = new TH1F("raw_uperpen_puppipt260270", "raw_uperpen_puppipt260270", 200, 25, 25);
	raw_uperpen_puppipt260270->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt290300 = new TH1F("raw_uperpen_puppipt290300", "raw_uperpen_puppipt290300", 200, 25, 25);
	raw_uperpen_puppipt290300->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt320330 = new TH1F("raw_uperpen_puppipt320330", "raw_uperpen_puppipt320330", 200, 25, 25);
	raw_uperpen_puppipt320330->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt330340 = new TH1F("raw_uperpen_puppipt330340", "raw_uperpen_puppipt330340", 200, 25, 25);
	raw_uperpen_puppipt330340->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt360370 = new TH1F("raw_uperpen_puppipt360370", "raw_uperpen_puppipt360370", 200, 25, 25);
	raw_uperpen_puppipt360370->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt390400 = new TH1F("raw_uperpen_puppipt390400", "raw_uperpen_puppipt390400", 200, 25, 25);
	raw_uperpen_puppipt390400->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uperpen_puppipt420430 = new TH1F("raw_uperpen_puppipt420430", "raw_uperpen_puppipt420430", 200, 25, 25);
	raw_uperpen_puppipt420430->GetXaxis()->SetTitle("u#perp");
	raw_uperpen_puppipt420430->GetYaxis()->SetTitle("Entries");



    //rawuparaqt vs puppi pt > 100
    TH1F *raw_uparaqt_puppipt100110 = new TH1F("raw_uparaqt_puppipt100110", "raw_uparaqt_puppipt100110", 200, 25, 25);
	raw_uparaqt_puppipt100110->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt100110->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt110120 = new TH1F("raw_uparaqt_puppipt110120", "raw_uparaqt_puppipt110120", 200, 25, 25);
	raw_uparaqt_puppipt110120->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt110120->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt120130 = new TH1F("raw_uparaqt_puppipt120130", "raw_uparaqt_puppipt120130", 200, 25, 25);
	raw_uparaqt_puppipt120130->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt120130->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt130140 = new TH1F("raw_uparaqt_puppipt130140", "raw_uparaqt_puppipt130140", 200, 25, 25);
	raw_uparaqt_puppipt130140->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt130140->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt140150 = new TH1F("raw_uparaqt_puppipt140150", "raw_uparaqt_puppipt140150", 200, 25, 25);
	raw_uparaqt_puppipt140150->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt140150->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt150160 = new TH1F("raw_uparaqt_puppipt150160", "raw_uparaqt_puppipt150160", 200, 25, 25);
	raw_uparaqt_puppipt150160->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt150160->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt160170 = new TH1F("raw_uparaqt_puppipt160170", "raw_uparaqt_puppipt160170", 200, 25, 25);
	raw_uparaqt_puppipt160170->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt160170->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt170180 = new TH1F("raw_uparaqt_puppipt170180", "raw_uparaqt_puppipt170180", 200, 25, 25);
	raw_uparaqt_puppipt170180->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt170180->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt180190 = new TH1F("raw_uparaqt_puppipt180190", "raw_uparaqt_puppipt180190", 200, 25, 25);
	raw_uparaqt_puppipt180190->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt180190->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt190200 = new TH1F("raw_uparaqt_puppipt190200", "raw_uparaqt_puppipt190200", 200, 25, 25);
	raw_uparaqt_puppipt190200->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt190200->GetYaxis()->SetTitle("Entries");


    //rawuparaqt vs puppi pt > 200
    TH1F *raw_uparaqt_puppipt200210 = new TH1F("raw_uparaqt_puppipt200210", "raw_uparaqt_puppipt200210", 200, 25, 25);
	raw_uparaqt_puppipt200210->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt200210->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt230240 = new TH1F("raw_uparaqt_puppipt230240", "raw_uparaqt_puppipt230240", 200, 25, 25);
	raw_uparaqt_puppipt230240->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt230240->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt260270 = new TH1F("raw_uparaqt_puppipt260270", "raw_uparaqt_puppipt260270", 200, 25, 25);
	raw_uparaqt_puppipt260270->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt260270->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt290300 = new TH1F("raw_uparaqt_puppipt290300", "raw_uparaqt_puppipt290300", 200, 25, 25);
	raw_uparaqt_puppipt290300->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt290300->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt320330 = new TH1F("raw_uparaqt_puppipt320330", "raw_uparaqt_puppipt320330", 200, 25, 25);
	raw_uparaqt_puppipt320330->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt320330->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt330340 = new TH1F("raw_uparaqt_puppipt330340", "raw_uparaqt_puppipt330340", 200, 25, 25);
	raw_uparaqt_puppipt330340->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt330340->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt360370 = new TH1F("raw_uparaqt_puppipt360370", "raw_uparaqt_puppipt360370", 200, 25, 25);
	raw_uparaqt_puppipt360370->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt360370->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt390400 = new TH1F("raw_uparaqt_puppipt390400", "raw_uparaqt_puppipt390400", 200, 25, 25);
	raw_uparaqt_puppipt390400->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt390400->GetYaxis()->SetTitle("Entries");

    TH1F *raw_uparaqt_puppipt420430 = new TH1F("raw_uparaqt_puppipt420430", "raw_uparaqt_puppipt420430", 200, 25, 25);
	raw_uparaqt_puppipt420430->GetXaxis()->SetTitle("u#parallel+q_{T}");
	raw_uparaqt_puppipt420430->GetYaxis()->SetTitle("Entries");



     /*
    ======--------======--======--======---------=======---------=======------===
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--====--======--======--======--=====--=======--=====--=========--=====
    ======--------======--======--======---------=======---------=========--=====
    ======--============--======--======--==============--================--=====
    ======--============----------======--==============--==============------===
    */




    TH1F *raw_paraqt = new TH1F("raw_paraqt", "raw_paraqt", 200, 25, 25);
	raw_paraqt->GetXaxis()->SetTitle("u#parallel + q_{T}");
	raw_paraqt->GetYaxis()->SetTitle("Entries");

    TH1F *rawmet_tightphopt = new TH1F("rawmet_tightphopt", "rawmet_tightphopt", 200, 0, 400);
	rawmet_tightphopt->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawmet_tightphopt->GetYaxis()->SetTitle("Entries");

    TH1F *rawpho_pxdist = new TH1F("rawpho_pxdist", "rawpho_pxdist", 200, 25, 25);
	rawpho_pxdist->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpho_pxdist->GetYaxis()->SetTitle("Entries");

    TH1F *rawpho_pydist = new TH1F("rawpho_pydist", "rawpho_pydist", 200, 25, 25);
	rawpho_pydist->GetXaxis()->SetTitle("p_{T} (GeV)");
	rawpho_pydist->GetYaxis()->SetTitle("Entries");

    TH1F *raw_response = new TH1F("raw_response", "raw_response", 200, 25, 25);
	raw_response->GetXaxis()->SetTitle("q_{T} (GeV)");
	raw_response->GetYaxis()->SetTitle("- <u#parallel> / <q_{T}>");

    
    //--------------------------2D METs compairson-----------------------
    TH2F *pf_vs_puppi100 = new TH2F("pf_vs_puppi100", "pf_vs_puppi100", 200, 25, 25, 200, 25, 25);
	pf_vs_puppi100->GetXaxis()->SetTitle("PFMET");
	pf_vs_puppi100->GetYaxis()->SetTitle("PuppiMET");

    TH2F *pf_vs_rawpuppi100 = new TH2F("pf_vs_rawpuppi100", "pf_vs_rawpuppi100", 200, 25, 25, 200, 25, 25);
	pf_vs_rawpuppi100->GetXaxis()->SetTitle("PFMET");
	pf_vs_rawpuppi100->GetYaxis()->SetTitle("rawPuppiMET");

    TH2F *raw_vs_puppi100 = new TH2F("raw_vs_puppi100", "raw_vs_puppi100", 200, 25, 25, 200, 25, 25);
	raw_vs_puppi100->GetXaxis()->SetTitle("rawMET");
	raw_vs_puppi100->GetYaxis()->SetTitle("PuppiMET");

    TH2F *raw_vs_rawpuppi100 = new TH2F("raw_vs_rawpuppi100", "raw_vs_rawpuppi100", 200, 25, 25, 200, 25, 25);
	raw_vs_rawpuppi100->GetXaxis()->SetTitle("rawMET");
	raw_vs_rawpuppi100->GetYaxis()->SetTitle("rawPuppiMET");

    TH2F *rawpuppi_vs_puppi100 = new TH2F("rawpuppi_vs_puppi100", "rawpuppi_vs_puppi100", 200, 25, 25, 200, 25, 25);
	rawpuppi_vs_puppi100->GetXaxis()->SetTitle("rawPuppiMET");
	rawpuppi_vs_puppi100->GetYaxis()->SetTitle("PuppiMET");

    TH2F *pf_vs_raw100 = new TH2F("pf_vs_raw100", "pf_vs_raw100", 200, 25, 25, 200, 25, 25);
	pf_vs_raw100->GetXaxis()->SetTitle("pfMET");
	pf_vs_raw100->GetYaxis()->SetTitle("rawMET");


    TH2F *pf_vs_puppi200 = new TH2F("pf_vs_puppi200", "pf_vs_puppi200", 200, 25, 25, 200, 25, 25);
	pf_vs_puppi200->GetXaxis()->SetTitle("PFMET");
	pf_vs_puppi200->GetYaxis()->SetTitle("PuppiMET");

    TH2F *pf_vs_rawpuppi200 = new TH2F("pf_vs_rawpuppi200", "pf_vs_rawpuppi200", 200, 25, 25, 200, 25, 25);
	pf_vs_rawpuppi200->GetXaxis()->SetTitle("PFMET");
	pf_vs_rawpuppi200->GetYaxis()->SetTitle("rawPuppiMET");

    TH2F *raw_vs_puppi200 = new TH2F("raw_vs_puppi200", "raw_vs_puppi200", 200, 25, 25, 200, 25, 25);
	raw_vs_puppi200->GetXaxis()->SetTitle("rawMET");
	raw_vs_puppi200->GetYaxis()->SetTitle("PuppiMET");

    TH2F *raw_vs_rawpuppi200 = new TH2F("raw_vs_rawpuppi200", "raw_vs_rawpuppi200", 200, 25, 25, 200, 25, 25);
	raw_vs_rawpuppi200->GetXaxis()->SetTitle("rawMET");
	raw_vs_rawpuppi200->GetYaxis()->SetTitle("rawPuppiMET");

    TH2F *rawpuppi_vs_puppi200 = new TH2F("rawpuppi_vs_puppi200", "rawpuppi_vs_puppi200", 200, 25, 25, 200, 25, 25);
	rawpuppi_vs_puppi200->GetXaxis()->SetTitle("rawPuppiMET");
	rawpuppi_vs_puppi200->GetYaxis()->SetTitle("PuppiMET");

    TH2F *pf_vs_raw200 = new TH2F("pf_vs_raw200", "pf_vs_raw200", 200, 25, 25, 200, 25, 25);
	pf_vs_raw200->GetXaxis()->SetTitle("pfMET");
	pf_vs_raw200->GetYaxis()->SetTitle("rawMET");

    


    //--------------------------test area-------------------------------
    TH1F *nPho_one = new TH1F("nPho_one", "nPho_one", 200, 25, 25);
	nPho_one->GetXaxis()->SetTitle("nPho");
	nPho_one->GetYaxis()->SetTitle("Entries");

    TH1F *nJet_one = new TH1F("nJet_one", "nJet_one", 200, 25, 25);
	nJet_one->GetXaxis()->SetTitle("nJet");
	nJet_one->GetYaxis()->SetTitle("Entries");

    TH1F *nEle_one = new TH1F("nEle_one", "nEle_one", 200, 25, 25);
	nEle_one->GetXaxis()->SetTitle("nEle");
	nEle_one->GetYaxis()->SetTitle("Entries");

    TH1F *nMu_one = new TH1F("nMu_one", "nMu_one", 200, 25, 25);
	nMu_one->GetXaxis()->SetTitle("nMu");
	nMu_one->GetYaxis()->SetTitle("Entries");


    ofstream file1, file2, file3, file4, file5;

    file2.open("photonID.txt", ios::out);

    file3.open("jetID.txt", ios::out);

    file4.open("electronID.txt", ios::out);

    file5.open("muonID.txt", ios::out);



    cout << "Start analysis on jetmet_tightpho_GJET200to400.C " << endl;

	auto timenow_start = chrono::system_clock::to_time_t(chrono::system_clock::now());

	cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;


    //========================================================//
    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) 
    {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        if ((jentry % 10000) == 0)
			// to print the number of processed entries
			std::cout << "Processed: " << jentry << std::endl;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        //========================================================//


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
        ((HLTPho >> 33 & 1) == 1) ||     //&& HLTPho < 17179869184 //HLT_Photon120_R9Id90_HE10_IsoM_v
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

                    puppi_response->Fill(((*phoEt)[ipuppiPho]), (-(puppi_par)/((*phoEt)[ipuppiPho])));

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
                        
                        puppi_abs_scale5060->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara5060->Fill(puppi_par_bins);
                        puppi_uperpen5060->Fill(puppi_pen_bins);
                        puppi_uparaqt5060->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 71 && (*phoEt)[ipuppiPho] < 90)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale8090->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara8090->Fill(puppi_par_bins);
                        puppi_uperpen8090->Fill(puppi_pen_bins);
                        puppi_uparaqt8090->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 91 && (*phoEt)[ipuppiPho] < 110)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale100110->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara100110->Fill(puppi_par_bins);
                        puppi_uperpen100110->Fill(puppi_pen_bins);
                        puppi_uparaqt100110->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 111 && (*phoEt)[ipuppiPho] < 130)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale120130->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara120130->Fill(puppi_par_bins);
                        puppi_uperpen120130->Fill(puppi_pen_bins);
                        puppi_uparaqt120130->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 131 && (*phoEt)[ipuppiPho] < 150)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale140150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara140150->Fill(puppi_par_bins);
                        puppi_uperpen140150->Fill(puppi_pen_bins);
                        puppi_uparaqt140150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 151 && (*phoEt)[ipuppiPho] < 170)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale160170->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara160170->Fill(puppi_par_bins);
                        puppi_uperpen160170->Fill(puppi_pen_bins);
                        puppi_uparaqt160170->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 171 && (*phoEt)[ipuppiPho] < 190)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale190200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara190200->Fill(puppi_par_bins);
                        puppi_uperpen190200->Fill(puppi_pen_bins);
                        puppi_uparaqt190200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 191 && (*phoEt)[ipuppiPho] < 210)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale220230->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara220230->Fill(puppi_par_bins);
                        puppi_uperpen220230->Fill(puppi_pen_bins);
                        puppi_uparaqt220230->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 211 && (*phoEt)[ipuppiPho] < 230)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale250260->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara250260->Fill(puppi_par_bins);
                        puppi_uperpen250260->Fill(puppi_pen_bins);
                        puppi_uparaqt250260->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 231 && (*phoEt)[ipuppiPho] < 250)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scaled260270->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara260270->Fill(puppi_par_bins);
                        puppi_uperpen260270->Fill(puppi_pen_bins);
                        puppi_uparaqt260270->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 251 && (*phoEt)[ipuppiPho] < 270)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale280290->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara280290->Fill(puppi_par_bins);
                        puppi_uperpen280290->Fill(puppi_pen_bins);
                        puppi_uparaqt280290->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 271 && (*phoEt)[ipuppiPho] < 290)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale300310->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara300310->Fill(puppi_par_bins);
                        puppi_uperpen300310->Fill(puppi_pen_bins);
                        puppi_uparaqt300310->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 291 && (*phoEt)[ipuppiPho] < 310)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale310320->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara310320->Fill(puppi_par_bins);
                        puppi_uperpen310320->Fill(puppi_pen_bins);
                        puppi_uparaqt310320->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 311 && (*phoEt)[ipuppiPho] < 330)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale320330->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara320330->Fill(puppi_par_bins);
                        puppi_uperpen320330->Fill(puppi_pen_bins);
                        puppi_uparaqt320330->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 331 && (*phoEt)[ipuppiPho] < 350)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale350360->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara350360->Fill(puppi_par_bins);
                        puppi_uperpen350360->Fill(puppi_pen_bins);
                        puppi_uparaqt350360->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 351 && (*phoEt)[ipuppiPho] < 370)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale380390->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara380390->Fill(puppi_par_bins);
                        puppi_uperpen380390->Fill(puppi_pen_bins);
                        puppi_uparaqt380390->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if ((*phoEt)[ipuppiPho] > 371 && (*phoEt)[ipuppiPho] < 400)
                    {
                        //cal photon_px and photon_py
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale390400->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara390400->Fill(puppi_par_bins);
                        puppi_uperpen390400->Fill(puppi_pen_bins);
                        puppi_uparaqt390400->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
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

                    //vs puppipt > 100
                    if (puppiMET_pt > 100 && puppiMET_pt < 110)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt100110->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt100110->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt100110->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt100110->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }


                    if (puppiMET_pt > 110 && puppiMET_pt < 120)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt110120->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt110120->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt110120->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt110120->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 120 && puppiMET_pt < 130)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt120130->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt120130->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt120130->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt120130->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 130 && puppiMET_pt < 140)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt130140->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt130140->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt130140->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt130140->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 140 && puppiMET_pt < 150)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt140150->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt140150->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt140150->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt140150->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 150 && puppiMET_pt < 160)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt150160->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt150160->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt150160->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt150160->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 160 && puppiMET_pt < 170)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt160170->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt160170->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt160170->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt160170->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 170 && puppiMET_pt < 180)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt170180->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt170180->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt170180->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt170180->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 180 && puppiMET_pt < 190)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt180190->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt180190->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt180190->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt180190->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 190 && puppiMET_pt < 200)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt190200->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt190200->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt190200->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt190200->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }


                    if (puppiMET_pt > 200 && puppiMET_pt < 210)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt200210->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt200210->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt200210->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt200210->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 230 && puppiMET_pt < 240)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt230240->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt230240->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt230240->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt230240->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 260 && puppiMET_pt < 270)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt260270->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt260270->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt260270->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt260270->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 290 && puppiMET_pt < 300)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt290300->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt290300->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt290300->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt290300->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 320 && puppiMET_pt < 330)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt320330->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt320330->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt320330->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt320330->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 330 && puppiMET_pt < 340)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt330340->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt330340->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt330340->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt330340->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 360 && puppiMET_pt < 370)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt360370->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt360370->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt360370->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt360370->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 390 && puppiMET_pt < 400)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt390400->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt390400->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt390400->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt390400->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }

                    if (puppiMET_pt > 420 && puppiMET_pt < 430)
                    {
                        puppipho_px_bins = puppipho_px_bins + ((*phoEt)[ipuppiPho]*cos((*phoPhi)[ipuppiPho]));
                        puppipho_py_bins = puppipho_py_bins + ((*phoEt)[ipuppiPho]*sin((*phoPhi)[ipuppiPho]));

                        puppi_pen_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_py_bins - (-puppiMET_py - puppipho_py_bins) * puppipho_px_bins)/((*phoEt)[ipuppiPho]);
                        puppi_par_bins = ((-puppiMET_px - puppipho_px_bins) * puppipho_px_bins + (-puppiMET_py - puppipho_py_bins) * puppipho_py_bins)/((*phoEt)[ipuppiPho]);
                        
                        puppi_abs_scale_puppipt420430->Fill(-(puppi_par_bins)/((*phoEt)[ipuppiPho]));
                        puppi_upara_puppipt420430->Fill(puppi_par_bins);
                        puppi_uperpen_puppipt420430->Fill(puppi_pen_bins);
                        puppi_uparaqt_puppipt420430->Fill(puppi_par_bins + (*phoEt)[ipuppiPho]);
                    }
                }


                //rawpuppimet
                for (int irawpuppiPho = 0; irawpuppiPho < nPho; irawpuppiPho++)
                {
                    rawpuppimetphi->Fill(rawPuppiMETPhi);

                    rawpuppimet_tightphopt->Fill((*phoEt)[irawpuppiPho]);

                    //fill the px py from ggNtuplizer
                    rawpuppimetpt->Fill(rawPuppiMET_pt);
                    rawpuppimet_px->Fill(rawPuppiMET_px);
                    rawpuppimet_py->Fill(rawPuppiMET_py);

                    
                    //try to calculate and also fill the one from ggNtuplizer
                    //rawPuppimet_px = rawpuppimet.pt*m.cos(rawpuppimet.phi)
                    rawpuppimet_px_cal->Fill((rawPuppiMET_pt)*cos(rawPuppiMETPhi));
                    rawpuppimet_py_cal->Fill((rawPuppiMET_pt)*sin(rawPuppiMETPhi));
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

                    rawpuppi_response->Fill(((*phoEt)[irawpuppiPho]), (-(rawpuppi_par)/((*phoEt)[irawpuppiPho])));

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
                        
                        rawpuppi_abs_scale5060->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara5060->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen5060->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt5060->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 71 && (*phoEt)[irawpuppiPho] < 90)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale8090->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara8090->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen8090->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt8090->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 91 && (*phoEt)[irawpuppiPho] < 110)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale100110->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara100110->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen100110->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt100110->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 111 && (*phoEt)[irawpuppiPho] < 130)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale120130->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara120130->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen120130->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt120130->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 131 && (*phoEt)[irawpuppiPho] < 150)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale140150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara140150->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen140150->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt140150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 151 && (*phoEt)[irawpuppiPho] < 170)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale160170->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara160170->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen160170->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt160170->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 171 && (*phoEt)[irawpuppiPho] < 190)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale190200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara190200->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen190200->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt190200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 191 && (*phoEt)[irawpuppiPho] < 210)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale220230->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara220230->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen220230->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt220230->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 211 && (*phoEt)[irawpuppiPho] < 230)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale250260->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara250260->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen250260->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt250260->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 231 && (*phoEt)[irawpuppiPho] < 250)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scaled260270->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara260270->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen260270->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt260270->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 251 && (*phoEt)[irawpuppiPho] < 270)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale280290->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara280290->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen280290->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt280290->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 271 && (*phoEt)[irawpuppiPho] < 290)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale300310->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara300310->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen300310->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt300310->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 291 && (*phoEt)[irawpuppiPho] < 310)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale310320->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara310320->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen310320->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt310320->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 311 && (*phoEt)[irawpuppiPho] < 330)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale320330->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara320330->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen320330->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt320330->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 331 && (*phoEt)[irawpuppiPho] < 350)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale350360->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara350360->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen350360->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt350360->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 351 && (*phoEt)[irawpuppiPho] < 370)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale380390->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara380390->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen380390->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt380390->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if ((*phoEt)[irawpuppiPho] > 371 && (*phoEt)[irawpuppiPho] < 400)
                    {
                        //cal photon_px and photon_py
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale390400->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara390400->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen390400->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt390400->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
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

                    //vs puppipt > 100
                    if (rawPuppiMET_pt > 100 && rawPuppiMET_pt < 110)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt100110->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt100110->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt100110->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt100110->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }


                    if (rawPuppiMET_pt > 110 && rawPuppiMET_pt < 120)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt110120->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt110120->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt110120->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt110120->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 120 && rawPuppiMET_pt < 130)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt120130->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt120130->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt120130->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt120130->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 130 && rawPuppiMET_pt < 140)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt130140->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt130140->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt130140->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt130140->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 140 && rawPuppiMET_pt < 150)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt140150->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt140150->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt140150->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt140150->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 150 && rawPuppiMET_pt < 160)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt150160->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt150160->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt150160->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt150160->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 160 && rawPuppiMET_pt < 170)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt160170->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt160170->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt160170->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt160170->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 170 && rawPuppiMET_pt < 180)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt170180->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt170180->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt170180->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt170180->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 180 && rawPuppiMET_pt < 190)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt180190->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt180190->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt180190->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt180190->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 190 && rawPuppiMET_pt < 200)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt190200->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt190200->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt190200->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt190200->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }


                    if (rawPuppiMET_pt > 200 && rawPuppiMET_pt < 210)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt200210->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt200210->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt200210->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt200210->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 230 && rawPuppiMET_pt < 240)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt230240->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt230240->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt230240->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt230240->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 260 && rawPuppiMET_pt < 270)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt260270->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt260270->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt260270->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt260270->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 290 && rawPuppiMET_pt < 300)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt290300->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt290300->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt290300->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt290300->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 320 && rawPuppiMET_pt < 330)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt320330->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt320330->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt320330->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt320330->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 330 && rawPuppiMET_pt < 340)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt330340->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt330340->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt330340->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt330340->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 360 && rawPuppiMET_pt < 370)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt360370->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt360370->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt360370->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt360370->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 390 && rawPuppiMET_pt < 400)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt390400->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt390400->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt390400->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt390400->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
                    }

                    if (rawPuppiMET_pt > 420 && rawPuppiMET_pt < 430)
                    {
                        rawpuppipho_px_bins = rawpuppipho_px_bins + ((*phoEt)[irawpuppiPho]*cos((*phoPhi)[irawpuppiPho]));
                        rawpuppipho_py_bins = rawpuppipho_py_bins + ((*phoEt)[irawpuppiPho]*sin((*phoPhi)[irawpuppiPho]));

                        rawpuppi_pen_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_py_bins - (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_px_bins)/((*phoEt)[irawpuppiPho]);
                        rawpuppi_par_bins = ((-rawPuppiMET_px - rawpuppipho_px_bins) * rawpuppipho_px_bins + (-rawPuppiMET_py - rawpuppipho_py_bins) * rawpuppipho_py_bins)/((*phoEt)[irawpuppiPho]);
                        
                        rawpuppi_abs_scale_puppipt420430->Fill(-(rawpuppi_par_bins)/((*phoEt)[irawpuppiPho]));
                        rawpuppi_upara_puppipt420430->Fill(rawpuppi_par_bins);
                        rawpuppi_uperpen_puppipt420430->Fill(rawpuppi_pen_bins);
                        rawpuppi_uparaqt_puppipt420430->Fill(rawpuppi_par_bins + (*phoEt)[irawpuppiPho]);
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

                    pf_response->Fill(((*phoEt)[ipfPho]), (-(pf_par)/((*phoEt)[ipfPho])));

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
                        
                        pf_abs_scale5060->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara5060->Fill(pf_par_bins);
                        pf_uperpen5060->Fill(pf_pen_bins);
                        pf_uparaqt5060->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 71 && (*phoEt)[ipfPho] < 90)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale8090->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara8090->Fill(pf_par_bins);
                        pf_uperpen8090->Fill(pf_pen_bins);
                        pf_uparaqt8090->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 91 && (*phoEt)[ipfPho] < 110)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale100110->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara100110->Fill(pf_par_bins);
                        pf_uperpen100110->Fill(pf_pen_bins);
                        pf_uparaqt100110->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 111 && (*phoEt)[ipfPho] < 130)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale120130->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara120130->Fill(pf_par_bins);
                        pf_uperpen120130->Fill(pf_pen_bins);
                        pf_uparaqt120130->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 131 && (*phoEt)[ipfPho] < 150)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale140150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara140150->Fill(pf_par_bins);
                        pf_uperpen140150->Fill(pf_pen_bins);
                        pf_uparaqt140150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 151 && (*phoEt)[ipfPho] < 170)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale160170->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara160170->Fill(pf_par_bins);
                        pf_uperpen160170->Fill(pf_pen_bins);
                        pf_uparaqt160170->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 171 && (*phoEt)[ipfPho] < 190)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale190200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara190200->Fill(pf_par_bins);
                        pf_uperpen190200->Fill(pf_pen_bins);
                        pf_uparaqt190200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 191 && (*phoEt)[ipfPho] < 210)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale220230->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara220230->Fill(pf_par_bins);
                        pf_uperpen220230->Fill(pf_pen_bins);
                        pf_uparaqt220230->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 211 && (*phoEt)[ipfPho] < 230)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale250260->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara250260->Fill(pf_par_bins);
                        pf_uperpen250260->Fill(pf_pen_bins);
                        pf_uparaqt250260->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 231 && (*phoEt)[ipfPho] < 250)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scaled260270->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara260270->Fill(pf_par_bins);
                        pf_uperpen260270->Fill(pf_pen_bins);
                        pf_uparaqt260270->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 251 && (*phoEt)[ipfPho] < 270)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale280290->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara280290->Fill(pf_par_bins);
                        pf_uperpen280290->Fill(pf_pen_bins);
                        pf_uparaqt280290->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 271 && (*phoEt)[ipfPho] < 290)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale300310->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara300310->Fill(pf_par_bins);
                        pf_uperpen300310->Fill(pf_pen_bins);
                        pf_uparaqt300310->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 291 && (*phoEt)[ipfPho] < 310)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale310320->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara310320->Fill(pf_par_bins);
                        pf_uperpen310320->Fill(pf_pen_bins);
                        pf_uparaqt310320->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 311 && (*phoEt)[ipfPho] < 330)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale320330->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara320330->Fill(pf_par_bins);
                        pf_uperpen320330->Fill(pf_pen_bins);
                        pf_uparaqt320330->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 331 && (*phoEt)[ipfPho] < 350)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale350360->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara350360->Fill(pf_par_bins);
                        pf_uperpen350360->Fill(pf_pen_bins);
                        pf_uparaqt350360->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 351 && (*phoEt)[ipfPho] < 370)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale380390->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara380390->Fill(pf_par_bins);
                        pf_uperpen380390->Fill(pf_pen_bins);
                        pf_uparaqt380390->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if ((*phoEt)[ipfPho] > 371 && (*phoEt)[ipfPho] < 400)
                    {
                        //cal photon_px and photon_py
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale390400->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara390400->Fill(pf_par_bins);
                        pf_uperpen390400->Fill(pf_pen_bins);
                        pf_uparaqt390400->Fill(pf_par_bins + (*phoEt)[ipfPho]);
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

                    //vs puppipt > 100
                    if (pfMET_pt > 100 && pfMET_pt < 110)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt100110->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt100110->Fill(pf_par_bins);
                        pf_uperpen_puppipt100110->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt100110->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }


                    if (pfMET_pt > 110 && pfMET_pt < 120)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt110120->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt110120->Fill(pf_par_bins);
                        pf_uperpen_puppipt110120->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt110120->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 120 && pfMET_pt < 130)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt120130->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt120130->Fill(pf_par_bins);
                        pf_uperpen_puppipt120130->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt120130->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 130 && pfMET_pt < 140)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt130140->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt130140->Fill(pf_par_bins);
                        pf_uperpen_puppipt130140->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt130140->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 140 && pfMET_pt < 150)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt140150->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt140150->Fill(pf_par_bins);
                        pf_uperpen_puppipt140150->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt140150->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 150 && pfMET_pt < 160)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt150160->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt150160->Fill(pf_par_bins);
                        pf_uperpen_puppipt150160->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt150160->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 160 && pfMET_pt < 170)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt160170->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt160170->Fill(pf_par_bins);
                        pf_uperpen_puppipt160170->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt160170->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 170 && pfMET_pt < 180)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt170180->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt170180->Fill(pf_par_bins);
                        pf_uperpen_puppipt170180->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt170180->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 180 && pfMET_pt < 190)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt180190->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt180190->Fill(pf_par_bins);
                        pf_uperpen_puppipt180190->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt180190->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 190 && pfMET_pt < 200)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt190200->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt190200->Fill(pf_par_bins);
                        pf_uperpen_puppipt190200->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt190200->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }


                    if (pfMET_pt > 200 && pfMET_pt < 210)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt200210->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt200210->Fill(pf_par_bins);
                        pf_uperpen_puppipt200210->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt200210->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 230 && pfMET_pt < 240)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt230240->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt230240->Fill(pf_par_bins);
                        pf_uperpen_puppipt230240->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt230240->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 260 && pfMET_pt < 270)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt260270->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt260270->Fill(pf_par_bins);
                        pf_uperpen_puppipt260270->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt260270->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 290 && pfMET_pt < 300)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt290300->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt290300->Fill(pf_par_bins);
                        pf_uperpen_puppipt290300->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt290300->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 320 && pfMET_pt < 330)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt320330->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt320330->Fill(pf_par_bins);
                        pf_uperpen_puppipt320330->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt320330->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 330 && pfMET_pt < 340)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt330340->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt330340->Fill(pf_par_bins);
                        pf_uperpen_puppipt330340->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt330340->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 360 && pfMET_pt < 370)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt360370->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt360370->Fill(pf_par_bins);
                        pf_uperpen_puppipt360370->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt360370->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 390 && pfMET_pt < 400)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt390400->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt390400->Fill(pf_par_bins);
                        pf_uperpen_puppipt390400->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt390400->Fill(pf_par_bins + (*phoEt)[ipfPho]);
                    }

                    if (pfMET_pt > 420 && pfMET_pt < 430)
                    {
                        pfpho_px_bins = pfpho_px_bins + ((*phoEt)[ipfPho]*cos((*phoPhi)[ipfPho]));
                        pfpho_py_bins = pfpho_py_bins + ((*phoEt)[ipfPho]*sin((*phoPhi)[ipfPho]));

                        pf_pen_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_py_bins - (-pfMET_py - pfpho_py_bins) * pfpho_px_bins)/((*phoEt)[ipfPho]);
                        pf_par_bins = ((-pfMET_px - pfpho_px_bins) * pfpho_px_bins + (-pfMET_py - pfpho_py_bins) * pfpho_py_bins)/((*phoEt)[ipfPho]);
                        
                        pf_abs_scale_puppipt420430->Fill(-(pf_par_bins)/((*phoEt)[ipfPho]));
                        pf_upara_puppipt420430->Fill(pf_par_bins);
                        pf_uperpen_puppipt420430->Fill(pf_pen_bins);
                        pf_uparaqt_puppipt420430->Fill(pf_par_bins + (*phoEt)[ipfPho]);
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

                    raw_response->Fill(((*phoEt)[irawPho]), (-(raw_par)/((*phoEt)[irawPho])));

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
                        
                        raw_abs_scale5060->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara5060->Fill(raw_par_bins);
                        raw_uperpen5060->Fill(raw_pen_bins);
                        raw_uparaqt5060->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 71 && (*phoEt)[irawPho] < 90)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale8090->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara8090->Fill(raw_par_bins);
                        raw_uperpen8090->Fill(raw_pen_bins);
                        raw_uparaqt8090->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 91 && (*phoEt)[irawPho] < 110)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale100110->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara100110->Fill(raw_par_bins);
                        raw_uperpen100110->Fill(raw_pen_bins);
                        raw_uparaqt100110->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 111 && (*phoEt)[irawPho] < 130)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale120130->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara120130->Fill(raw_par_bins);
                        raw_uperpen120130->Fill(raw_pen_bins);
                        raw_uparaqt120130->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 131 && (*phoEt)[irawPho] < 150)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale140150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara140150->Fill(raw_par_bins);
                        raw_uperpen140150->Fill(raw_pen_bins);
                        raw_uparaqt140150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 151 && (*phoEt)[irawPho] < 170)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale160170->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara160170->Fill(raw_par_bins);
                        raw_uperpen160170->Fill(raw_pen_bins);
                        raw_uparaqt160170->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 171 && (*phoEt)[irawPho] < 190)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale190200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara190200->Fill(raw_par_bins);
                        raw_uperpen190200->Fill(raw_pen_bins);
                        raw_uparaqt190200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 191 && (*phoEt)[irawPho] < 210)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale220230->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara220230->Fill(raw_par_bins);
                        raw_uperpen220230->Fill(raw_pen_bins);
                        raw_uparaqt220230->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 211 && (*phoEt)[irawPho] < 230)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale250260->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara250260->Fill(raw_par_bins);
                        raw_uperpen250260->Fill(raw_pen_bins);
                        raw_uparaqt250260->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 231 && (*phoEt)[irawPho] < 250)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scaled260270->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara260270->Fill(raw_par_bins);
                        raw_uperpen260270->Fill(raw_pen_bins);
                        raw_uparaqt260270->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 251 && (*phoEt)[irawPho] < 270)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale280290->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara280290->Fill(raw_par_bins);
                        raw_uperpen280290->Fill(raw_pen_bins);
                        raw_uparaqt280290->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 271 && (*phoEt)[irawPho] < 290)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale300310->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara300310->Fill(raw_par_bins);
                        raw_uperpen300310->Fill(raw_pen_bins);
                        raw_uparaqt300310->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 291 && (*phoEt)[irawPho] < 310)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale310320->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara310320->Fill(raw_par_bins);
                        raw_uperpen310320->Fill(raw_pen_bins);
                        raw_uparaqt310320->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 311 && (*phoEt)[irawPho] < 330)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale320330->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara320330->Fill(raw_par_bins);
                        raw_uperpen320330->Fill(raw_pen_bins);
                        raw_uparaqt320330->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 331 && (*phoEt)[irawPho] < 350)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale350360->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara350360->Fill(raw_par_bins);
                        raw_uperpen350360->Fill(raw_pen_bins);
                        raw_uparaqt350360->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 351 && (*phoEt)[irawPho] < 370)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale380390->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara380390->Fill(raw_par_bins);
                        raw_uperpen380390->Fill(raw_pen_bins);
                        raw_uparaqt380390->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if ((*phoEt)[irawPho] > 371 && (*phoEt)[irawPho] < 400)
                    {
                        //cal photon_px and photon_py
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale390400->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara390400->Fill(raw_par_bins);
                        raw_uperpen390400->Fill(raw_pen_bins);
                        raw_uparaqt390400->Fill(raw_par_bins + (*phoEt)[irawPho]);
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

                    //vs puppipt > 100
                    if (rawMET_pt > 100 && rawMET_pt < 110)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt100110->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt100110->Fill(raw_par_bins);
                        raw_uperpen_puppipt100110->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt100110->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }


                    if (rawMET_pt > 110 && rawMET_pt < 120)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt110120->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt110120->Fill(raw_par_bins);
                        raw_uperpen_puppipt110120->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt110120->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 120 && rawMET_pt < 130)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt120130->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt120130->Fill(raw_par_bins);
                        raw_uperpen_puppipt120130->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt120130->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 130 && rawMET_pt < 140)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt130140->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt130140->Fill(raw_par_bins);
                        raw_uperpen_puppipt130140->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt130140->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 140 && rawMET_pt < 150)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt140150->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt140150->Fill(raw_par_bins);
                        raw_uperpen_puppipt140150->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt140150->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 150 && rawMET_pt < 160)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt150160->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt150160->Fill(raw_par_bins);
                        raw_uperpen_puppipt150160->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt150160->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 160 && rawMET_pt < 170)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt160170->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt160170->Fill(raw_par_bins);
                        raw_uperpen_puppipt160170->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt160170->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 170 && rawMET_pt < 180)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt170180->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt170180->Fill(raw_par_bins);
                        raw_uperpen_puppipt170180->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt170180->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 180 && rawMET_pt < 190)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt180190->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt180190->Fill(raw_par_bins);
                        raw_uperpen_puppipt180190->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt180190->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 190 && rawMET_pt < 200)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt190200->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt190200->Fill(raw_par_bins);
                        raw_uperpen_puppipt190200->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt190200->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }


                    if (rawMET_pt > 200 && rawMET_pt < 210)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt200210->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt200210->Fill(raw_par_bins);
                        raw_uperpen_puppipt200210->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt200210->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 230 && rawMET_pt < 240)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt230240->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt230240->Fill(raw_par_bins);
                        raw_uperpen_puppipt230240->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt230240->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 260 && rawMET_pt < 270)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt260270->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt260270->Fill(raw_par_bins);
                        raw_uperpen_puppipt260270->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt260270->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 290 && rawMET_pt < 300)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt290300->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt290300->Fill(raw_par_bins);
                        raw_uperpen_puppipt290300->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt290300->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 320 && rawMET_pt < 330)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt320330->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt320330->Fill(raw_par_bins);
                        raw_uperpen_puppipt320330->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt320330->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 330 && rawMET_pt < 340)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt330340->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt330340->Fill(raw_par_bins);
                        raw_uperpen_puppipt330340->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt330340->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 360 && rawMET_pt < 370)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt360370->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt360370->Fill(raw_par_bins);
                        raw_uperpen_puppipt360370->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt360370->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 390 && rawMET_pt < 400)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt390400->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt390400->Fill(raw_par_bins);
                        raw_uperpen_puppipt390400->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt390400->Fill(raw_par_bins + (*phoEt)[irawPho]);
                    }

                    if (rawMET_pt > 420 && rawMET_pt < 430)
                    {
                        rawpho_px_bins = rawpho_px_bins + ((*phoEt)[irawPho]*cos((*phoPhi)[irawPho]));
                        rawpho_py_bins = rawpho_py_bins + ((*phoEt)[irawPho]*sin((*phoPhi)[irawPho]));

                        raw_pen_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_py_bins - (-rawMET_py - rawpho_py_bins) * rawpho_px_bins)/((*phoEt)[irawPho]);
                        raw_par_bins = ((-rawMET_px - rawpho_px_bins) * rawpho_px_bins + (-rawMET_py - rawpho_py_bins) * rawpho_py_bins)/((*phoEt)[irawPho]);
                        
                        raw_abs_scale_puppipt420430->Fill(-(raw_par_bins)/((*phoEt)[irawPho]));
                        raw_upara_puppipt420430->Fill(raw_par_bins);
                        raw_uperpen_puppipt420430->Fill(raw_pen_bins);
                        raw_uparaqt_puppipt420430->Fill(raw_par_bins + (*phoEt)[irawPho]);
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



    //puppimet with one tight pho id, one jet, no loose ele and muon
    puppimetphi->SetLineColor(1);
    puppimetphi->Write();

    puppimetpt->SetLineColor(1);
    puppimetpt->Write();

    puppimet_px_cal->SetLineColor(1);
    puppimet_px_cal->Write();

    puppimet_py_cal->SetLineColor(1);
    puppimet_py_cal->Write();

    puppimet_px->SetLineColor(1);
    puppimet_px->Write();

    puppimet_py->SetLineColor(1);
    puppimet_py->Write();

    uperpen_puppi->SetLineColor(1);
    uperpen_puppi->Write();

    uparall_puppi->SetLineColor(1);
    uparall_puppi->Write();

    puppi_abs_scale->SetLineColor(1);
    puppi_abs_scale->Write();

    puppi_paraqt->SetLineColor(1);
    puppi_paraqt->Write();

    puppimet_tightphopt->SetLineColor(1);
    puppimet_tightphopt->Write();

    puppipho_pxdist->SetLineColor(1);
    puppipho_pxdist->Write();

    puppipho_pydist->SetLineColor(1);
    puppipho_pydist->Write();

    puppi_response->SetLineColor(1);
    puppi_response->Write();

    puppimet_nGoodvtx->SetLineColor(1);
    puppimet_nGoodvtx->Write();


    puppi_abs_scale5060->SetLineColor(1);
    puppi_abs_scale5060->Write();

    puppi_abs_scale8090->SetLineColor(1);
    puppi_abs_scale8090->Write();

    puppi_abs_scale100110->SetLineColor(1);
    puppi_abs_scale100110->Write();

    puppi_abs_scale120130->SetLineColor(1);
    puppi_abs_scale120130->Write();

    puppi_abs_scale140150->SetLineColor(1);
    puppi_abs_scale140150->Write();

    puppi_abs_scale160170->SetLineColor(1);
    puppi_abs_scale160170->Write();

    puppi_abs_scale190200->SetLineColor(1);
    puppi_abs_scale190200->Write();

    puppi_abs_scale220230->SetLineColor(1);
    puppi_abs_scale220230->Write();

    puppi_abs_scale250260->SetLineColor(1);
    puppi_abs_scale250260->Write();

    puppi_abs_scaled260270->SetLineColor(1);
    puppi_abs_scaled260270->Write();

    puppi_abs_scale280290->SetLineColor(1);
    puppi_abs_scale280290->Write();

    puppi_abs_scale300310->SetLineColor(1);
    puppi_abs_scale300310->Write();

    puppi_abs_scale310320->SetLineColor(1);
    puppi_abs_scale310320->Write();

    puppi_abs_scale320330->SetLineColor(1);
    puppi_abs_scale320330->Write();

    puppi_abs_scale350360->SetLineColor(1);
    puppi_abs_scale350360->Write();

    puppi_abs_scale380390->SetLineColor(1);
    puppi_abs_scale380390->Write();

    puppi_abs_scale390400->SetLineColor(1);
    puppi_abs_scale390400->Write();


    
    puppi_upara5060->SetLineColor(1);
    puppi_upara5060->Write();

    puppi_upara8090->SetLineColor(1);
    puppi_upara8090->Write();

    puppi_upara100110->SetLineColor(1);
    puppi_upara100110->Write();

    puppi_upara120130->SetLineColor(1);
    puppi_upara120130->Write();

    puppi_upara140150->SetLineColor(1);
    puppi_upara140150->Write();

    puppi_upara160170->SetLineColor(1);
    puppi_upara160170->Write();

    puppi_upara190200->SetLineColor(1);
    puppi_upara190200->Write();

    puppi_upara220230->SetLineColor(1);
    puppi_upara220230->Write();

    puppi_upara250260->SetLineColor(1);
    puppi_upara250260->Write();

    puppi_upara260270->SetLineColor(1);
    puppi_upara260270->Write();

    puppi_upara280290->SetLineColor(1);
    puppi_upara280290->Write();

    puppi_upara300310->SetLineColor(1);
    puppi_upara300310->Write();

    puppi_upara310320->SetLineColor(1);
    puppi_upara310320->Write();

    puppi_upara320330->SetLineColor(1);
    puppi_upara320330->Write();

    puppi_upara350360->SetLineColor(1);
    puppi_upara350360->Write();

    puppi_upara380390->SetLineColor(1);
    puppi_upara380390->Write();

    puppi_upara390400->SetLineColor(1);
    puppi_upara390400->Write();


    
    puppi_uperpen5060->SetLineColor(1);
    puppi_uperpen5060->Write();

    puppi_uperpen8090->SetLineColor(1);
    puppi_uperpen8090->Write();

    puppi_uperpen100110->SetLineColor(1);
    puppi_uperpen100110->Write();

    puppi_uperpen120130->SetLineColor(1);
    puppi_uperpen120130->Write();

    puppi_uperpen140150->SetLineColor(1);
    puppi_uperpen140150->Write();

    puppi_uperpen160170->SetLineColor(1);
    puppi_uperpen160170->Write();

    puppi_uperpen190200->SetLineColor(1);
    puppi_uperpen190200->Write();

    puppi_uperpen220230->SetLineColor(1);
    puppi_uperpen220230->Write();

    puppi_uperpen250260->SetLineColor(1);
    puppi_uperpen250260->Write();

    puppi_uperpen260270->SetLineColor(1);
    puppi_uperpen260270->Write();

    puppi_uperpen280290->SetLineColor(1);
    puppi_uperpen280290->Write();

    puppi_uperpen300310->SetLineColor(1);
    puppi_uperpen300310->Write();

    puppi_uperpen310320->SetLineColor(1);
    puppi_uperpen310320->Write();

    puppi_uperpen320330->SetLineColor(1);
    puppi_uperpen320330->Write();

    puppi_uperpen350360->SetLineColor(1);
    puppi_uperpen350360->Write();

    puppi_uperpen380390->SetLineColor(1);
    puppi_uperpen380390->Write();

    puppi_uperpen390400->SetLineColor(1);
    puppi_uperpen390400->Write();


    puppi_uparaqt5060->SetLineColor(1);
    puppi_uparaqt5060->Write();

    puppi_uparaqt8090->SetLineColor(1);
    puppi_uparaqt8090->Write();

    puppi_uparaqt100110->SetLineColor(1);
    puppi_uparaqt100110->Write();

    puppi_uparaqt120130->SetLineColor(1);
    puppi_uparaqt120130->Write();

    puppi_uparaqt140150->SetLineColor(1);
    puppi_uparaqt140150->Write();

    puppi_uparaqt160170->SetLineColor(1);
    puppi_uparaqt160170->Write();

    puppi_uparaqt190200->SetLineColor(1);
    puppi_uparaqt190200->Write();

    puppi_uparaqt220230->SetLineColor(1);
    puppi_uparaqt220230->Write();

    puppi_uparaqt250260->SetLineColor(1);
    puppi_uparaqt250260->Write();

    puppi_uparaqt260270->SetLineColor(1);
    puppi_uparaqt260270->Write();

    puppi_uparaqt280290->SetLineColor(1);
    puppi_uparaqt280290->Write();

    puppi_uparaqt300310->SetLineColor(1);
    puppi_uparaqt300310->Write();

    puppi_uparaqt310320->SetLineColor(1);
    puppi_uparaqt310320->Write();

    puppi_uparaqt320330->SetLineColor(1);
    puppi_uparaqt320330->Write();

    puppi_uparaqt350360->SetLineColor(1);
    puppi_uparaqt350360->Write();

    puppi_uparaqt380390->SetLineColor(1);
    puppi_uparaqt380390->Write();

    puppi_uparaqt390400->SetLineColor(1);
    puppi_uparaqt390400->Write();



    puppi_uparanvtx04->SetLineColor(1);
    puppi_uparanvtx04->Write();

    puppi_uparanvtx46->SetLineColor(1);
    puppi_uparanvtx46->Write();

    puppi_uparanvtx68->SetLineColor(1);
    puppi_uparanvtx68->Write();

    puppi_uparanvtx812->SetLineColor(1);
    puppi_uparanvtx812->Write();
    
    puppi_uparanvtx1218->SetLineColor(1);
    puppi_uparanvtx1218->Write();

    puppi_uparanvtx1824->SetLineColor(1);
    puppi_uparanvtx1824->Write();

    puppi_uparanvtx2428->SetLineColor(1);
    puppi_uparanvtx2428->Write();

    puppi_uparanvtx2832->SetLineColor(1);
    puppi_uparanvtx2832->Write();

    puppi_uparanvtx3236->SetLineColor(1);
    puppi_uparanvtx3236->Write();

    puppi_uparanvtx3640->SetLineColor(1);
    puppi_uparanvtx3640->Write();


    puppi_uperpennvtx04->SetLineColor(1);
    puppi_uperpennvtx04->Write();

    puppi_uperpennvtx46->SetLineColor(1);
    puppi_uperpennvtx46->Write();

    puppi_uperpennvtx68->SetLineColor(1);
    puppi_uperpennvtx68->Write();

    puppi_uperpennvtx812->SetLineColor(1);
    puppi_uperpennvtx812->Write();
    
    puppi_uperpennvtx1218->SetLineColor(1);
    puppi_uperpennvtx1218->Write();

    puppi_uperpennvtx1824->SetLineColor(1);
    puppi_uperpennvtx1824->Write();

    puppi_uperpennvtx2428->SetLineColor(1);
    puppi_uperpennvtx2428->Write();

    puppi_uperpennvtx2832->SetLineColor(1);
    puppi_uperpennvtx2832->Write();

    puppi_uperpennvtx3236->SetLineColor(1);
    puppi_uperpennvtx3236->Write();

    puppi_uperpennvtx3640->SetLineColor(1);
    puppi_uperpennvtx3640->Write();


    

    puppi_abs_scale_puppipt100110->SetLineColor(1);
    puppi_abs_scale_puppipt100110->Write();

    puppi_abs_scale_puppipt110120->SetLineColor(1);
    puppi_abs_scale_puppipt110120->Write();

    puppi_abs_scale_puppipt120130->SetLineColor(1);
    puppi_abs_scale_puppipt120130->Write();

    puppi_abs_scale_puppipt130140->SetLineColor(1);
    puppi_abs_scale_puppipt130140->Write();

    puppi_abs_scale_puppipt140150->SetLineColor(1);
    puppi_abs_scale_puppipt140150->Write();

    puppi_abs_scale_puppipt150160->SetLineColor(1);
    puppi_abs_scale_puppipt150160->Write();

    puppi_abs_scale_puppipt160170->SetLineColor(1);
    puppi_abs_scale_puppipt160170->Write();

    puppi_abs_scale_puppipt170180->SetLineColor(1);
    puppi_abs_scale_puppipt170180->Write();

    puppi_abs_scale_puppipt180190->SetLineColor(1);
    puppi_abs_scale_puppipt180190->Write();

    puppi_abs_scale_puppipt190200->SetLineColor(1);
    puppi_abs_scale_puppipt190200->Write();

    puppi_abs_scale_puppipt200210->SetLineColor(1);
    puppi_abs_scale_puppipt200210->Write();

    puppi_abs_scale_puppipt230240->SetLineColor(1);
    puppi_abs_scale_puppipt230240->Write();

    puppi_abs_scale_puppipt260270->SetLineColor(1);
    puppi_abs_scale_puppipt260270->Write();

    puppi_abs_scale_puppipt290300->SetLineColor(1);
    puppi_abs_scale_puppipt290300->Write();

    puppi_abs_scale_puppipt320330->SetLineColor(1);
    puppi_abs_scale_puppipt320330->Write();

    puppi_abs_scale_puppipt330340->SetLineColor(1);
    puppi_abs_scale_puppipt330340->Write();

    puppi_abs_scale_puppipt360370->SetLineColor(1);
    puppi_abs_scale_puppipt360370->Write();

    puppi_abs_scale_puppipt390400->SetLineColor(1);
    puppi_abs_scale_puppipt390400->Write();

    puppi_abs_scale_puppipt420430->SetLineColor(1);
    puppi_abs_scale_puppipt420430->Write();


    puppi_upara_puppipt100110->SetLineColor(1);
    puppi_upara_puppipt100110->Write();

    puppi_upara_puppipt110120->SetLineColor(1);
    puppi_upara_puppipt110120->Write();

    puppi_upara_puppipt120130->SetLineColor(1);
    puppi_upara_puppipt120130->Write();

    puppi_upara_puppipt130140->SetLineColor(1);
    puppi_upara_puppipt130140->Write();

    puppi_upara_puppipt140150->SetLineColor(1);
    puppi_upara_puppipt140150->Write();

    puppi_upara_puppipt150160->SetLineColor(1);
    puppi_upara_puppipt150160->Write();

    puppi_upara_puppipt160170->SetLineColor(1);
    puppi_upara_puppipt160170->Write();

    puppi_upara_puppipt170180->SetLineColor(1);
    puppi_upara_puppipt170180->Write();

    puppi_upara_puppipt180190->SetLineColor(1);
    puppi_upara_puppipt180190->Write();

    puppi_upara_puppipt190200->SetLineColor(1);
    puppi_upara_puppipt190200->Write();

    puppi_upara_puppipt200210->SetLineColor(1);
    puppi_upara_puppipt200210->Write();

    puppi_upara_puppipt230240->SetLineColor(1);
    puppi_upara_puppipt230240->Write();

    puppi_upara_puppipt260270->SetLineColor(1);
    puppi_upara_puppipt260270->Write();

    puppi_upara_puppipt290300->SetLineColor(1);
    puppi_upara_puppipt290300->Write();

    puppi_upara_puppipt320330->SetLineColor(1);
    puppi_upara_puppipt320330->Write();

    puppi_upara_puppipt330340->SetLineColor(1);
    puppi_upara_puppipt330340->Write();

    puppi_upara_puppipt360370->SetLineColor(1);
    puppi_upara_puppipt360370->Write();

    puppi_upara_puppipt390400->SetLineColor(1);
    puppi_upara_puppipt390400->Write();

    puppi_upara_puppipt420430->SetLineColor(1);
    puppi_upara_puppipt420430->Write();


    puppi_uperpen_puppipt100110->SetLineColor(1);
    puppi_uperpen_puppipt100110->Write();

    puppi_uperpen_puppipt110120->SetLineColor(1);
    puppi_uperpen_puppipt110120->Write();

    puppi_uperpen_puppipt120130->SetLineColor(1);
    puppi_uperpen_puppipt120130->Write();

    puppi_uperpen_puppipt130140->SetLineColor(1);
    puppi_uperpen_puppipt130140->Write();

    puppi_uperpen_puppipt140150->SetLineColor(1);
    puppi_uperpen_puppipt140150->Write();

    puppi_uperpen_puppipt150160->SetLineColor(1);
    puppi_uperpen_puppipt150160->Write();

    puppi_uperpen_puppipt160170->SetLineColor(1);
    puppi_uperpen_puppipt160170->Write();

    puppi_uperpen_puppipt170180->SetLineColor(1);
    puppi_uperpen_puppipt170180->Write();

    puppi_uperpen_puppipt180190->SetLineColor(1);
    puppi_uperpen_puppipt180190->Write();

    puppi_uperpen_puppipt190200->SetLineColor(1);
    puppi_uperpen_puppipt190200->Write();

    puppi_uperpen_puppipt200210->SetLineColor(1);
    puppi_uperpen_puppipt200210->Write();

    puppi_uperpen_puppipt230240->SetLineColor(1);
    puppi_uperpen_puppipt230240->Write();

    puppi_uperpen_puppipt260270->SetLineColor(1);
    puppi_uperpen_puppipt260270->Write();

    puppi_uperpen_puppipt290300->SetLineColor(1);
    puppi_uperpen_puppipt290300->Write();

    puppi_uperpen_puppipt320330->SetLineColor(1);
    puppi_uperpen_puppipt320330->Write();

    puppi_uperpen_puppipt330340->SetLineColor(1);
    puppi_uperpen_puppipt330340->Write();

    puppi_uperpen_puppipt360370->SetLineColor(1);
    puppi_uperpen_puppipt360370->Write();

    puppi_uperpen_puppipt390400->SetLineColor(1);
    puppi_uperpen_puppipt390400->Write();

    puppi_uperpen_puppipt420430->SetLineColor(1);
    puppi_uperpen_puppipt420430->Write();


    puppi_uparaqt_puppipt100110->SetLineColor(1);
    puppi_uparaqt_puppipt100110->Write();

    puppi_uparaqt_puppipt110120->SetLineColor(1);
    puppi_uparaqt_puppipt110120->Write();

    puppi_uparaqt_puppipt120130->SetLineColor(1);
    puppi_uparaqt_puppipt120130->Write();

    puppi_uparaqt_puppipt130140->SetLineColor(1);
    puppi_uparaqt_puppipt130140->Write();

    puppi_uparaqt_puppipt140150->SetLineColor(1);
    puppi_uparaqt_puppipt140150->Write();

    puppi_uparaqt_puppipt150160->SetLineColor(1);
    puppi_uparaqt_puppipt150160->Write();

    puppi_uparaqt_puppipt160170->SetLineColor(1);
    puppi_uparaqt_puppipt160170->Write();

    puppi_uparaqt_puppipt170180->SetLineColor(1);
    puppi_uparaqt_puppipt170180->Write();

    puppi_uparaqt_puppipt180190->SetLineColor(1);
    puppi_uparaqt_puppipt180190->Write();

    puppi_uparaqt_puppipt190200->SetLineColor(1);
    puppi_uparaqt_puppipt190200->Write();

    puppi_uparaqt_puppipt200210->SetLineColor(1);
    puppi_uparaqt_puppipt200210->Write();

    puppi_uparaqt_puppipt230240->SetLineColor(1);
    puppi_uparaqt_puppipt230240->Write();

    puppi_uparaqt_puppipt260270->SetLineColor(1);
    puppi_uparaqt_puppipt260270->Write();

    puppi_uparaqt_puppipt290300->SetLineColor(1);
    puppi_uparaqt_puppipt290300->Write();

    puppi_uparaqt_puppipt320330->SetLineColor(1);
    puppi_uparaqt_puppipt320330->Write();

    puppi_uparaqt_puppipt330340->SetLineColor(1);
    puppi_uparaqt_puppipt330340->Write();

    puppi_uparaqt_puppipt360370->SetLineColor(1);
    puppi_uparaqt_puppipt360370->Write();

    puppi_uparaqt_puppipt390400->SetLineColor(1);
    puppi_uparaqt_puppipt390400->Write();

    puppi_uparaqt_puppipt420430->SetLineColor(1);
    puppi_uparaqt_puppipt420430->Write();






    //rawpuppi
    rawpuppimetphi->SetLineColor(1);
    rawpuppimetphi->Write();

    rawpuppimetpt->SetLineColor(1);
    rawpuppimetpt->Write();

    rawpuppimet_px_cal->SetLineColor(1);
    rawpuppimet_px_cal->Write();

    rawpuppimet_py_cal->SetLineColor(1);
    rawpuppimet_py_cal->Write();

    rawpuppimet_px->SetLineColor(1);
    rawpuppimet_px->Write();

    rawpuppimet_py->SetLineColor(1);
    rawpuppimet_py->Write();

    uperpen_rawpuppi->SetLineColor(1);
    uperpen_rawpuppi->Write();

    uparall_rawpuppi->SetLineColor(1);
    uparall_rawpuppi->Write();

    rawpuppi_abs_scale->SetLineColor(1);
    rawpuppi_abs_scale->Write();

    rawpuppi_paraqt->SetLineColor(1);
    rawpuppi_paraqt->Write();

    rawpuppimet_tightphopt->SetLineColor(1);
    rawpuppimet_tightphopt->Write();

    rawpuppipho_pxdist->SetLineColor(1);
    rawpuppipho_pxdist->Write();

    rawpuppipho_pydist->SetLineColor(1);
    rawpuppipho_pydist->Write();

    rawpuppi_response->SetLineColor(1);
    rawpuppi_response->Write();

    rawpuppimet_nGoodvtx->SetLineColor(1);
    rawpuppimet_nGoodvtx->Write();


    rawpuppi_abs_scale5060->SetLineColor(1);
    rawpuppi_abs_scale5060->Write();

    rawpuppi_abs_scale8090->SetLineColor(1);
    rawpuppi_abs_scale8090->Write();

    rawpuppi_abs_scale100110->SetLineColor(1);
    rawpuppi_abs_scale100110->Write();

    rawpuppi_abs_scale120130->SetLineColor(1);
    rawpuppi_abs_scale120130->Write();

    rawpuppi_abs_scale140150->SetLineColor(1);
    rawpuppi_abs_scale140150->Write();

    rawpuppi_abs_scale160170->SetLineColor(1);
    rawpuppi_abs_scale160170->Write();

    rawpuppi_abs_scale190200->SetLineColor(1);
    rawpuppi_abs_scale190200->Write();

    rawpuppi_abs_scale220230->SetLineColor(1);
    rawpuppi_abs_scale220230->Write();

    rawpuppi_abs_scale250260->SetLineColor(1);
    rawpuppi_abs_scale250260->Write();

    rawpuppi_abs_scaled260270->SetLineColor(1);
    rawpuppi_abs_scaled260270->Write();

    rawpuppi_abs_scale280290->SetLineColor(1);
    rawpuppi_abs_scale280290->Write();

    rawpuppi_abs_scale300310->SetLineColor(1);
    rawpuppi_abs_scale300310->Write();

    rawpuppi_abs_scale310320->SetLineColor(1);
    rawpuppi_abs_scale310320->Write();

    rawpuppi_abs_scale320330->SetLineColor(1);
    rawpuppi_abs_scale320330->Write();

    rawpuppi_abs_scale350360->SetLineColor(1);
    rawpuppi_abs_scale350360->Write();

    rawpuppi_abs_scale380390->SetLineColor(1);
    rawpuppi_abs_scale380390->Write();

    rawpuppi_abs_scale390400->SetLineColor(1);
    rawpuppi_abs_scale390400->Write();


    
    rawpuppi_upara5060->SetLineColor(1);
    rawpuppi_upara5060->Write();

    rawpuppi_upara8090->SetLineColor(1);
    rawpuppi_upara8090->Write();

    rawpuppi_upara100110->SetLineColor(1);
    rawpuppi_upara100110->Write();

    rawpuppi_upara120130->SetLineColor(1);
    rawpuppi_upara120130->Write();

    rawpuppi_upara140150->SetLineColor(1);
    rawpuppi_upara140150->Write();

    rawpuppi_upara160170->SetLineColor(1);
    rawpuppi_upara160170->Write();

    rawpuppi_upara190200->SetLineColor(1);
    rawpuppi_upara190200->Write();

    rawpuppi_upara220230->SetLineColor(1);
    rawpuppi_upara220230->Write();

    rawpuppi_upara250260->SetLineColor(1);
    rawpuppi_upara250260->Write();

    rawpuppi_upara260270->SetLineColor(1);
    rawpuppi_upara260270->Write();

    rawpuppi_upara280290->SetLineColor(1);
    rawpuppi_upara280290->Write();

    rawpuppi_upara300310->SetLineColor(1);
    rawpuppi_upara300310->Write();

    rawpuppi_upara310320->SetLineColor(1);
    rawpuppi_upara310320->Write();

    rawpuppi_upara320330->SetLineColor(1);
    rawpuppi_upara320330->Write();

    rawpuppi_upara350360->SetLineColor(1);
    rawpuppi_upara350360->Write();

    rawpuppi_upara380390->SetLineColor(1);
    rawpuppi_upara380390->Write();

    rawpuppi_upara390400->SetLineColor(1);
    rawpuppi_upara390400->Write();


    
    rawpuppi_uperpen5060->SetLineColor(1);
    rawpuppi_uperpen5060->Write();

    rawpuppi_uperpen8090->SetLineColor(1);
    rawpuppi_uperpen8090->Write();

    rawpuppi_uperpen100110->SetLineColor(1);
    rawpuppi_uperpen100110->Write();

    rawpuppi_uperpen120130->SetLineColor(1);
    rawpuppi_uperpen120130->Write();

    rawpuppi_uperpen140150->SetLineColor(1);
    rawpuppi_uperpen140150->Write();

    rawpuppi_uperpen160170->SetLineColor(1);
    rawpuppi_uperpen160170->Write();

    rawpuppi_uperpen190200->SetLineColor(1);
    rawpuppi_uperpen190200->Write();

    rawpuppi_uperpen220230->SetLineColor(1);
    rawpuppi_uperpen220230->Write();

    rawpuppi_uperpen250260->SetLineColor(1);
    rawpuppi_uperpen250260->Write();

    rawpuppi_uperpen260270->SetLineColor(1);
    rawpuppi_uperpen260270->Write();

    rawpuppi_uperpen280290->SetLineColor(1);
    rawpuppi_uperpen280290->Write();

    rawpuppi_uperpen300310->SetLineColor(1);
    rawpuppi_uperpen300310->Write();

    rawpuppi_uperpen310320->SetLineColor(1);
    rawpuppi_uperpen310320->Write();

    rawpuppi_uperpen320330->SetLineColor(1);
    rawpuppi_uperpen320330->Write();

    rawpuppi_uperpen350360->SetLineColor(1);
    rawpuppi_uperpen350360->Write();

    rawpuppi_uperpen380390->SetLineColor(1);
    rawpuppi_uperpen380390->Write();

    rawpuppi_uperpen390400->SetLineColor(1);
    rawpuppi_uperpen390400->Write();


    rawpuppi_uparaqt5060->SetLineColor(1);
    rawpuppi_uparaqt5060->Write();

    rawpuppi_uparaqt8090->SetLineColor(1);
    rawpuppi_uparaqt8090->Write();

    rawpuppi_uparaqt100110->SetLineColor(1);
    rawpuppi_uparaqt100110->Write();

    rawpuppi_uparaqt120130->SetLineColor(1);
    rawpuppi_uparaqt120130->Write();

    rawpuppi_uparaqt140150->SetLineColor(1);
    rawpuppi_uparaqt140150->Write();

    rawpuppi_uparaqt160170->SetLineColor(1);
    rawpuppi_uparaqt160170->Write();

    rawpuppi_uparaqt190200->SetLineColor(1);
    rawpuppi_uparaqt190200->Write();

    rawpuppi_uparaqt220230->SetLineColor(1);
    rawpuppi_uparaqt220230->Write();

    rawpuppi_uparaqt250260->SetLineColor(1);
    rawpuppi_uparaqt250260->Write();

    rawpuppi_uparaqt260270->SetLineColor(1);
    rawpuppi_uparaqt260270->Write();

    rawpuppi_uparaqt280290->SetLineColor(1);
    rawpuppi_uparaqt280290->Write();

    rawpuppi_uparaqt300310->SetLineColor(1);
    rawpuppi_uparaqt300310->Write();

    rawpuppi_uparaqt310320->SetLineColor(1);
    rawpuppi_uparaqt310320->Write();

    rawpuppi_uparaqt320330->SetLineColor(1);
    rawpuppi_uparaqt320330->Write();

    rawpuppi_uparaqt350360->SetLineColor(1);
    rawpuppi_uparaqt350360->Write();

    rawpuppi_uparaqt380390->SetLineColor(1);
    rawpuppi_uparaqt380390->Write();

    rawpuppi_uparaqt390400->SetLineColor(1);
    rawpuppi_uparaqt390400->Write();



    rawpuppi_uparanvtx04->SetLineColor(1);
    rawpuppi_uparanvtx04->Write();

    rawpuppi_uparanvtx46->SetLineColor(1);
    rawpuppi_uparanvtx46->Write();

    rawpuppi_uparanvtx68->SetLineColor(1);
    rawpuppi_uparanvtx68->Write();

    rawpuppi_uparanvtx812->SetLineColor(1);
    rawpuppi_uparanvtx812->Write();
    
    rawpuppi_uparanvtx1218->SetLineColor(1);
    rawpuppi_uparanvtx1218->Write();

    rawpuppi_uparanvtx1824->SetLineColor(1);
    rawpuppi_uparanvtx1824->Write();

    rawpuppi_uparanvtx2428->SetLineColor(1);
    rawpuppi_uparanvtx2428->Write();

    rawpuppi_uparanvtx2832->SetLineColor(1);
    rawpuppi_uparanvtx2832->Write();

    rawpuppi_uparanvtx3236->SetLineColor(1);
    rawpuppi_uparanvtx3236->Write();

    rawpuppi_uparanvtx3640->SetLineColor(1);
    rawpuppi_uparanvtx3640->Write();


    rawpuppi_uperpennvtx04->SetLineColor(1);
    rawpuppi_uperpennvtx04->Write();

    rawpuppi_uperpennvtx46->SetLineColor(1);
    rawpuppi_uperpennvtx46->Write();

    rawpuppi_uperpennvtx68->SetLineColor(1);
    rawpuppi_uperpennvtx68->Write();

    rawpuppi_uperpennvtx812->SetLineColor(1);
    rawpuppi_uperpennvtx812->Write();
    
    rawpuppi_uperpennvtx1218->SetLineColor(1);
    rawpuppi_uperpennvtx1218->Write();

    rawpuppi_uperpennvtx1824->SetLineColor(1);
    rawpuppi_uperpennvtx1824->Write();

    rawpuppi_uperpennvtx2428->SetLineColor(1);
    rawpuppi_uperpennvtx2428->Write();

    rawpuppi_uperpennvtx2832->SetLineColor(1);
    rawpuppi_uperpennvtx2832->Write();

    rawpuppi_uperpennvtx3236->SetLineColor(1);
    rawpuppi_uperpennvtx3236->Write();

    rawpuppi_uperpennvtx3640->SetLineColor(1);
    rawpuppi_uperpennvtx3640->Write();



    

    rawpuppi_abs_scale_puppipt100110->SetLineColor(1);
    rawpuppi_abs_scale_puppipt100110->Write();

    rawpuppi_abs_scale_puppipt110120->SetLineColor(1);
    rawpuppi_abs_scale_puppipt110120->Write();

    rawpuppi_abs_scale_puppipt120130->SetLineColor(1);
    rawpuppi_abs_scale_puppipt120130->Write();

    rawpuppi_abs_scale_puppipt130140->SetLineColor(1);
    rawpuppi_abs_scale_puppipt130140->Write();

    rawpuppi_abs_scale_puppipt140150->SetLineColor(1);
    rawpuppi_abs_scale_puppipt140150->Write();

    rawpuppi_abs_scale_puppipt150160->SetLineColor(1);
    rawpuppi_abs_scale_puppipt150160->Write();

    rawpuppi_abs_scale_puppipt160170->SetLineColor(1);
    rawpuppi_abs_scale_puppipt160170->Write();

    rawpuppi_abs_scale_puppipt170180->SetLineColor(1);
    rawpuppi_abs_scale_puppipt170180->Write();

    rawpuppi_abs_scale_puppipt180190->SetLineColor(1);
    rawpuppi_abs_scale_puppipt180190->Write();

    rawpuppi_abs_scale_puppipt190200->SetLineColor(1);
    rawpuppi_abs_scale_puppipt190200->Write();

    rawpuppi_abs_scale_puppipt200210->SetLineColor(1);
    rawpuppi_abs_scale_puppipt200210->Write();

    rawpuppi_abs_scale_puppipt230240->SetLineColor(1);
    rawpuppi_abs_scale_puppipt230240->Write();

    rawpuppi_abs_scale_puppipt260270->SetLineColor(1);
    rawpuppi_abs_scale_puppipt260270->Write();

    rawpuppi_abs_scale_puppipt290300->SetLineColor(1);
    rawpuppi_abs_scale_puppipt290300->Write();

    rawpuppi_abs_scale_puppipt320330->SetLineColor(1);
    rawpuppi_abs_scale_puppipt320330->Write();

    rawpuppi_abs_scale_puppipt330340->SetLineColor(1);
    rawpuppi_abs_scale_puppipt330340->Write();

    rawpuppi_abs_scale_puppipt360370->SetLineColor(1);
    rawpuppi_abs_scale_puppipt360370->Write();

    rawpuppi_abs_scale_puppipt390400->SetLineColor(1);
    rawpuppi_abs_scale_puppipt390400->Write();

    rawpuppi_abs_scale_puppipt420430->SetLineColor(1);
    rawpuppi_abs_scale_puppipt420430->Write();

    rawpuppi_upara_puppipt100110->SetLineColor(1);
    rawpuppi_upara_puppipt100110->Write();

    rawpuppi_upara_puppipt110120->SetLineColor(1);
    rawpuppi_upara_puppipt110120->Write();

    rawpuppi_upara_puppipt120130->SetLineColor(1);
    rawpuppi_upara_puppipt120130->Write();

    rawpuppi_upara_puppipt130140->SetLineColor(1);
    rawpuppi_upara_puppipt130140->Write();

    rawpuppi_upara_puppipt140150->SetLineColor(1);
    rawpuppi_upara_puppipt140150->Write();

    rawpuppi_upara_puppipt150160->SetLineColor(1);
    rawpuppi_upara_puppipt150160->Write();

    rawpuppi_upara_puppipt160170->SetLineColor(1);
    rawpuppi_upara_puppipt160170->Write();

    rawpuppi_upara_puppipt170180->SetLineColor(1);
    rawpuppi_upara_puppipt170180->Write();

    rawpuppi_upara_puppipt180190->SetLineColor(1);
    rawpuppi_upara_puppipt180190->Write();

    rawpuppi_upara_puppipt190200->SetLineColor(1);
    rawpuppi_upara_puppipt190200->Write();

    rawpuppi_upara_puppipt200210->SetLineColor(1);
    rawpuppi_upara_puppipt200210->Write();

    rawpuppi_upara_puppipt230240->SetLineColor(1);
    rawpuppi_upara_puppipt230240->Write();

    rawpuppi_upara_puppipt260270->SetLineColor(1);
    rawpuppi_upara_puppipt260270->Write();

    rawpuppi_upara_puppipt290300->SetLineColor(1);
    rawpuppi_upara_puppipt290300->Write();

    rawpuppi_upara_puppipt320330->SetLineColor(1);
    rawpuppi_upara_puppipt320330->Write();

    rawpuppi_upara_puppipt330340->SetLineColor(1);
    rawpuppi_upara_puppipt330340->Write();

    rawpuppi_upara_puppipt360370->SetLineColor(1);
    rawpuppi_upara_puppipt360370->Write();

    rawpuppi_upara_puppipt390400->SetLineColor(1);
    rawpuppi_upara_puppipt390400->Write();

    rawpuppi_upara_puppipt420430->SetLineColor(1);
    rawpuppi_upara_puppipt420430->Write();


    rawpuppi_uperpen_puppipt100110->SetLineColor(1);
    rawpuppi_uperpen_puppipt100110->Write();

    rawpuppi_uperpen_puppipt110120->SetLineColor(1);
    rawpuppi_uperpen_puppipt110120->Write();

    rawpuppi_uperpen_puppipt120130->SetLineColor(1);
    rawpuppi_uperpen_puppipt120130->Write();

    rawpuppi_uperpen_puppipt130140->SetLineColor(1);
    rawpuppi_uperpen_puppipt130140->Write();

    rawpuppi_uperpen_puppipt140150->SetLineColor(1);
    rawpuppi_uperpen_puppipt140150->Write();

    rawpuppi_uperpen_puppipt150160->SetLineColor(1);
    rawpuppi_uperpen_puppipt150160->Write();

    rawpuppi_uperpen_puppipt160170->SetLineColor(1);
    rawpuppi_uperpen_puppipt160170->Write();

    rawpuppi_uperpen_puppipt170180->SetLineColor(1);
    rawpuppi_uperpen_puppipt170180->Write();

    rawpuppi_uperpen_puppipt180190->SetLineColor(1);
    rawpuppi_uperpen_puppipt180190->Write();

    rawpuppi_uperpen_puppipt190200->SetLineColor(1);
    rawpuppi_uperpen_puppipt190200->Write();

    rawpuppi_uperpen_puppipt200210->SetLineColor(1);
    rawpuppi_uperpen_puppipt200210->Write();

    rawpuppi_uperpen_puppipt230240->SetLineColor(1);
    rawpuppi_uperpen_puppipt230240->Write();

    rawpuppi_uperpen_puppipt260270->SetLineColor(1);
    rawpuppi_uperpen_puppipt260270->Write();

    rawpuppi_uperpen_puppipt290300->SetLineColor(1);
    rawpuppi_uperpen_puppipt290300->Write();

    rawpuppi_uperpen_puppipt320330->SetLineColor(1);
    rawpuppi_uperpen_puppipt320330->Write();

    rawpuppi_uperpen_puppipt330340->SetLineColor(1);
    rawpuppi_uperpen_puppipt330340->Write();

    rawpuppi_uperpen_puppipt360370->SetLineColor(1);
    rawpuppi_uperpen_puppipt360370->Write();

    rawpuppi_uperpen_puppipt390400->SetLineColor(1);
    rawpuppi_uperpen_puppipt390400->Write();

    rawpuppi_uperpen_puppipt420430->SetLineColor(1);
    rawpuppi_uperpen_puppipt420430->Write();


    rawpuppi_uparaqt_puppipt100110->SetLineColor(1);
    rawpuppi_uparaqt_puppipt100110->Write();

    rawpuppi_uparaqt_puppipt110120->SetLineColor(1);
    rawpuppi_uparaqt_puppipt110120->Write();

    rawpuppi_uparaqt_puppipt120130->SetLineColor(1);
    rawpuppi_uparaqt_puppipt120130->Write();

    rawpuppi_uparaqt_puppipt130140->SetLineColor(1);
    rawpuppi_uparaqt_puppipt130140->Write();

    rawpuppi_uparaqt_puppipt140150->SetLineColor(1);
    rawpuppi_uparaqt_puppipt140150->Write();

    rawpuppi_uparaqt_puppipt150160->SetLineColor(1);
    rawpuppi_uparaqt_puppipt150160->Write();

    rawpuppi_uparaqt_puppipt160170->SetLineColor(1);
    rawpuppi_uparaqt_puppipt160170->Write();

    rawpuppi_uparaqt_puppipt170180->SetLineColor(1);
    rawpuppi_uparaqt_puppipt170180->Write();

    rawpuppi_uparaqt_puppipt180190->SetLineColor(1);
    rawpuppi_uparaqt_puppipt180190->Write();

    rawpuppi_uparaqt_puppipt190200->SetLineColor(1);
    rawpuppi_uparaqt_puppipt190200->Write();

    rawpuppi_uparaqt_puppipt200210->SetLineColor(1);
    rawpuppi_uparaqt_puppipt200210->Write();

    rawpuppi_uparaqt_puppipt230240->SetLineColor(1);
    rawpuppi_uparaqt_puppipt230240->Write();

    rawpuppi_uparaqt_puppipt260270->SetLineColor(1);
    rawpuppi_uparaqt_puppipt260270->Write();

    rawpuppi_uparaqt_puppipt290300->SetLineColor(1);
    rawpuppi_uparaqt_puppipt290300->Write();

    rawpuppi_uparaqt_puppipt320330->SetLineColor(1);
    rawpuppi_uparaqt_puppipt320330->Write();

    rawpuppi_uparaqt_puppipt330340->SetLineColor(1);
    rawpuppi_uparaqt_puppipt330340->Write();

    rawpuppi_uparaqt_puppipt360370->SetLineColor(1);
    rawpuppi_uparaqt_puppipt360370->Write();

    rawpuppi_uparaqt_puppipt390400->SetLineColor(1);
    rawpuppi_uparaqt_puppipt390400->Write();

    rawpuppi_uparaqt_puppipt420430->SetLineColor(1);
    rawpuppi_uparaqt_puppipt420430->Write();



    

    //pfmet with only tight pho id
    pfmetphi->SetLineColor(1);
    pfmetphi->Write();

    pfmetpt->SetLineColor(1);
    pfmetpt->Write();

    pfmet_px->SetLineColor(1);
    pfmet_px->Write();

    pfmet_py->SetLineColor(1);
    pfmet_py->Write();

    uperpen_pf->SetLineColor(1);
    uperpen_pf->Write();

    uparall_pf->SetLineColor(1);
    uparall_pf->Write();

    pf_abs_scale->SetLineColor(1);
    pf_abs_scale->Write();

    pf_paraqt->SetLineColor(1);
    pf_paraqt->Write();

    




    //pfmet with one tight pho id, one jet, no loose ele and muon
    pfmetphi->SetLineColor(1);
    pfmetphi->Write();

    pfmetpt->SetLineColor(1);
    pfmetpt->Write();

    pfmet_px_cal->SetLineColor(1);
    pfmet_px_cal->Write();

    pfmet_py_cal->SetLineColor(1);
    pfmet_py_cal->Write();

    pfmet_px->SetLineColor(1);
    pfmet_px->Write();

    pfmet_py->SetLineColor(1);
    pfmet_py->Write();

    uperpen_pf->SetLineColor(1);
    uperpen_pf->Write();

    uparall_pf->SetLineColor(1);
    uparall_pf->Write();

    pf_abs_scale->SetLineColor(1);
    pf_abs_scale->Write();

    pf_paraqt->SetLineColor(1);
    pf_paraqt->Write();

    pfmet_tightphopt->SetLineColor(1);
    pfmet_tightphopt->Write();

    pfpho_pxdist->SetLineColor(1);
    pfpho_pxdist->Write();

    pfpho_pydist->SetLineColor(1);
    pfpho_pydist->Write();

    pf_response->SetLineColor(1);
    pf_response->Write();

    pf_nGoodvtx->SetLineColor(1);
    pf_nGoodvtx->Write();




    pf_abs_scale5060->SetLineColor(1);
    pf_abs_scale5060->Write();

    pf_abs_scale8090->SetLineColor(1);
    pf_abs_scale8090->Write();

    pf_abs_scale100110->SetLineColor(1);
    pf_abs_scale100110->Write();

    pf_abs_scale120130->SetLineColor(1);
    pf_abs_scale120130->Write();

    pf_abs_scale140150->SetLineColor(1);
    pf_abs_scale140150->Write();

    pf_abs_scale160170->SetLineColor(1);
    pf_abs_scale160170->Write();

    pf_abs_scale190200->SetLineColor(1);
    pf_abs_scale190200->Write();

    pf_abs_scale220230->SetLineColor(1);
    pf_abs_scale220230->Write();

    pf_abs_scale250260->SetLineColor(1);
    pf_abs_scale250260->Write();

    pf_abs_scaled260270->SetLineColor(1);
    pf_abs_scaled260270->Write();

    pf_abs_scale280290->SetLineColor(1);
    pf_abs_scale280290->Write();

    pf_abs_scale300310->SetLineColor(1);
    pf_abs_scale300310->Write();

    pf_abs_scale310320->SetLineColor(1);
    pf_abs_scale310320->Write();

    pf_abs_scale320330->SetLineColor(1);
    pf_abs_scale320330->Write();

    pf_abs_scale350360->SetLineColor(1);
    pf_abs_scale350360->Write();

    pf_abs_scale380390->SetLineColor(1);
    pf_abs_scale380390->Write();

    pf_abs_scale390400->SetLineColor(1);
    pf_abs_scale390400->Write();


    
    pf_upara5060->SetLineColor(1);
    pf_upara5060->Write();

    pf_upara8090->SetLineColor(1);
    pf_upara8090->Write();

    pf_upara100110->SetLineColor(1);
    pf_upara100110->Write();

    pf_upara120130->SetLineColor(1);
    pf_upara120130->Write();

    pf_upara140150->SetLineColor(1);
    pf_upara140150->Write();

    pf_upara160170->SetLineColor(1);
    pf_upara160170->Write();

    pf_upara190200->SetLineColor(1);
    pf_upara190200->Write();

    pf_upara220230->SetLineColor(1);
    pf_upara220230->Write();

    pf_upara250260->SetLineColor(1);
    pf_upara250260->Write();

    pf_upara260270->SetLineColor(1);
    pf_upara260270->Write();

    pf_upara280290->SetLineColor(1);
    pf_upara280290->Write();

    pf_upara300310->SetLineColor(1);
    pf_upara300310->Write();

    pf_upara310320->SetLineColor(1);
    pf_upara310320->Write();

    pf_upara320330->SetLineColor(1);
    pf_upara320330->Write();

    pf_upara350360->SetLineColor(1);
    pf_upara350360->Write();

    pf_upara380390->SetLineColor(1);
    pf_upara380390->Write();

    pf_upara390400->SetLineColor(1);
    pf_upara390400->Write();


    
    pf_uperpen5060->SetLineColor(1);
    pf_uperpen5060->Write();

    pf_uperpen8090->SetLineColor(1);
    pf_uperpen8090->Write();

    pf_uperpen100110->SetLineColor(1);
    pf_uperpen100110->Write();

    pf_uperpen120130->SetLineColor(1);
    pf_uperpen120130->Write();

    pf_uperpen140150->SetLineColor(1);
    pf_uperpen140150->Write();

    pf_uperpen160170->SetLineColor(1);
    pf_uperpen160170->Write();

    pf_uperpen190200->SetLineColor(1);
    pf_uperpen190200->Write();

    pf_uperpen220230->SetLineColor(1);
    pf_uperpen220230->Write();

    pf_uperpen250260->SetLineColor(1);
    pf_uperpen250260->Write();

    pf_uperpen260270->SetLineColor(1);
    pf_uperpen260270->Write();

    pf_uperpen280290->SetLineColor(1);
    pf_uperpen280290->Write();

    pf_uperpen300310->SetLineColor(1);
    pf_uperpen300310->Write();

    pf_uperpen310320->SetLineColor(1);
    pf_uperpen310320->Write();

    pf_uperpen320330->SetLineColor(1);
    pf_uperpen320330->Write();

    pf_uperpen350360->SetLineColor(1);
    pf_uperpen350360->Write();

    pf_uperpen380390->SetLineColor(1);
    pf_uperpen380390->Write();

    pf_uperpen390400->SetLineColor(1);
    pf_uperpen390400->Write();


    pf_uparaqt5060->SetLineColor(1);
    pf_uparaqt5060->Write();

    pf_uparaqt8090->SetLineColor(1);
    pf_uparaqt8090->Write();

    pf_uparaqt100110->SetLineColor(1);
    pf_uparaqt100110->Write();

    pf_uparaqt120130->SetLineColor(1);
    pf_uparaqt120130->Write();

    pf_uparaqt140150->SetLineColor(1);
    pf_uparaqt140150->Write();

    pf_uparaqt160170->SetLineColor(1);
    pf_uparaqt160170->Write();

    pf_uparaqt190200->SetLineColor(1);
    pf_uparaqt190200->Write();

    pf_uparaqt220230->SetLineColor(1);
    pf_uparaqt220230->Write();

    pf_uparaqt250260->SetLineColor(1);
    pf_uparaqt250260->Write();

    pf_uparaqt260270->SetLineColor(1);
    pf_uparaqt260270->Write();

    pf_uparaqt280290->SetLineColor(1);
    pf_uparaqt280290->Write();

    pf_uparaqt300310->SetLineColor(1);
    pf_uparaqt300310->Write();

    pf_uparaqt310320->SetLineColor(1);
    pf_uparaqt310320->Write();

    pf_uparaqt320330->SetLineColor(1);
    pf_uparaqt320330->Write();

    pf_uparaqt350360->SetLineColor(1);
    pf_uparaqt350360->Write();

    pf_uparaqt380390->SetLineColor(1);
    pf_uparaqt380390->Write();

    pf_uparaqt390400->SetLineColor(1);
    pf_uparaqt390400->Write();



    pf_uparanvtx04->SetLineColor(1);
    pf_uparanvtx04->Write();

    pf_uparanvtx46->SetLineColor(1);
    pf_uparanvtx46->Write();

    pf_uparanvtx68->SetLineColor(1);
    pf_uparanvtx68->Write();

    pf_uparanvtx812->SetLineColor(1);
    pf_uparanvtx812->Write();
    
    pf_uparanvtx1218->SetLineColor(1);
    pf_uparanvtx1218->Write();

    pf_uparanvtx1824->SetLineColor(1);
    pf_uparanvtx1824->Write();

    pf_uparanvtx2428->SetLineColor(1);
    pf_uparanvtx2428->Write();

    pf_uparanvtx2832->SetLineColor(1);
    pf_uparanvtx2832->Write();

    pf_uparanvtx3236->SetLineColor(1);
    pf_uparanvtx3236->Write();

    pf_uparanvtx3640->SetLineColor(1);
    pf_uparanvtx3640->Write();


    pf_uperpennvtx04->SetLineColor(1);
    pf_uperpennvtx04->Write();

    pf_uperpennvtx46->SetLineColor(1);
    pf_uperpennvtx46->Write();

    pf_uperpennvtx68->SetLineColor(1);
    pf_uperpennvtx68->Write();

    pf_uperpennvtx812->SetLineColor(1);
    pf_uperpennvtx812->Write();
    
    pf_uperpennvtx1218->SetLineColor(1);
    pf_uperpennvtx1218->Write();

    pf_uperpennvtx1824->SetLineColor(1);
    pf_uperpennvtx1824->Write();

    pf_uperpennvtx2428->SetLineColor(1);
    pf_uperpennvtx2428->Write();

    pf_uperpennvtx2832->SetLineColor(1);
    pf_uperpennvtx2832->Write();

    pf_uperpennvtx3236->SetLineColor(1);
    pf_uperpennvtx3236->Write();

    pf_uperpennvtx3640->SetLineColor(1);
    pf_uperpennvtx3640->Write();



    

    pf_abs_scale_puppipt100110->SetLineColor(1);
    pf_abs_scale_puppipt100110->Write();

    pf_abs_scale_puppipt110120->SetLineColor(1);
    pf_abs_scale_puppipt110120->Write();

    pf_abs_scale_puppipt120130->SetLineColor(1);
    pf_abs_scale_puppipt120130->Write();

    pf_abs_scale_puppipt130140->SetLineColor(1);
    pf_abs_scale_puppipt130140->Write();

    pf_abs_scale_puppipt140150->SetLineColor(1);
    pf_abs_scale_puppipt140150->Write();

    pf_abs_scale_puppipt150160->SetLineColor(1);
    pf_abs_scale_puppipt150160->Write();

    pf_abs_scale_puppipt160170->SetLineColor(1);
    pf_abs_scale_puppipt160170->Write();

    pf_abs_scale_puppipt170180->SetLineColor(1);
    pf_abs_scale_puppipt170180->Write();

    pf_abs_scale_puppipt180190->SetLineColor(1);
    pf_abs_scale_puppipt180190->Write();

    pf_abs_scale_puppipt190200->SetLineColor(1);
    pf_abs_scale_puppipt190200->Write();

    pf_abs_scale_puppipt200210->SetLineColor(1);
    pf_abs_scale_puppipt200210->Write();

    pf_abs_scale_puppipt230240->SetLineColor(1);
    pf_abs_scale_puppipt230240->Write();

    pf_abs_scale_puppipt260270->SetLineColor(1);
    pf_abs_scale_puppipt260270->Write();

    pf_abs_scale_puppipt290300->SetLineColor(1);
    pf_abs_scale_puppipt290300->Write();

    pf_abs_scale_puppipt320330->SetLineColor(1);
    pf_abs_scale_puppipt320330->Write();

    pf_abs_scale_puppipt330340->SetLineColor(1);
    pf_abs_scale_puppipt330340->Write();

    pf_abs_scale_puppipt360370->SetLineColor(1);
    pf_abs_scale_puppipt360370->Write();

    pf_abs_scale_puppipt390400->SetLineColor(1);
    pf_abs_scale_puppipt390400->Write();

    pf_abs_scale_puppipt420430->SetLineColor(1);
    pf_abs_scale_puppipt420430->Write();


    pf_upara_puppipt100110->SetLineColor(1);
    pf_upara_puppipt100110->Write();

    pf_upara_puppipt110120->SetLineColor(1);
    pf_upara_puppipt110120->Write();

    pf_upara_puppipt120130->SetLineColor(1);
    pf_upara_puppipt120130->Write();

    pf_upara_puppipt130140->SetLineColor(1);
    pf_upara_puppipt130140->Write();

    pf_upara_puppipt140150->SetLineColor(1);
    pf_upara_puppipt140150->Write();

    pf_upara_puppipt150160->SetLineColor(1);
    pf_upara_puppipt150160->Write();

    pf_upara_puppipt160170->SetLineColor(1);
    pf_upara_puppipt160170->Write();

    pf_upara_puppipt170180->SetLineColor(1);
    pf_upara_puppipt170180->Write();

    pf_upara_puppipt180190->SetLineColor(1);
    pf_upara_puppipt180190->Write();

    pf_upara_puppipt190200->SetLineColor(1);
    pf_upara_puppipt190200->Write();

    pf_upara_puppipt200210->SetLineColor(1);
    pf_upara_puppipt200210->Write();

    pf_upara_puppipt230240->SetLineColor(1);
    pf_upara_puppipt230240->Write();

    pf_upara_puppipt260270->SetLineColor(1);
    pf_upara_puppipt260270->Write();

    pf_upara_puppipt290300->SetLineColor(1);
    pf_upara_puppipt290300->Write();

    pf_upara_puppipt320330->SetLineColor(1);
    pf_upara_puppipt320330->Write();

    pf_upara_puppipt330340->SetLineColor(1);
    pf_upara_puppipt330340->Write();

    pf_upara_puppipt360370->SetLineColor(1);
    pf_upara_puppipt360370->Write();

    pf_upara_puppipt390400->SetLineColor(1);
    pf_upara_puppipt390400->Write();

    pf_upara_puppipt420430->SetLineColor(1);
    pf_upara_puppipt420430->Write();


    pf_uperpen_puppipt100110->SetLineColor(1);
    pf_uperpen_puppipt100110->Write();

    pf_uperpen_puppipt110120->SetLineColor(1);
    pf_uperpen_puppipt110120->Write();

    pf_uperpen_puppipt120130->SetLineColor(1);
    pf_uperpen_puppipt120130->Write();

    pf_uperpen_puppipt130140->SetLineColor(1);
    pf_uperpen_puppipt130140->Write();

    pf_uperpen_puppipt140150->SetLineColor(1);
    pf_uperpen_puppipt140150->Write();

    pf_uperpen_puppipt150160->SetLineColor(1);
    pf_uperpen_puppipt150160->Write();

    pf_uperpen_puppipt160170->SetLineColor(1);
    pf_uperpen_puppipt160170->Write();

    pf_uperpen_puppipt170180->SetLineColor(1);
    pf_uperpen_puppipt170180->Write();

    pf_uperpen_puppipt180190->SetLineColor(1);
    pf_uperpen_puppipt180190->Write();

    pf_uperpen_puppipt190200->SetLineColor(1);
    pf_uperpen_puppipt190200->Write();

    pf_uperpen_puppipt200210->SetLineColor(1);
    pf_uperpen_puppipt200210->Write();

    pf_uperpen_puppipt230240->SetLineColor(1);
    pf_uperpen_puppipt230240->Write();

    pf_uperpen_puppipt260270->SetLineColor(1);
    pf_uperpen_puppipt260270->Write();

    pf_uperpen_puppipt290300->SetLineColor(1);
    pf_uperpen_puppipt290300->Write();

    pf_uperpen_puppipt320330->SetLineColor(1);
    pf_uperpen_puppipt320330->Write();

    pf_uperpen_puppipt330340->SetLineColor(1);
    pf_uperpen_puppipt330340->Write();

    pf_uperpen_puppipt360370->SetLineColor(1);
    pf_uperpen_puppipt360370->Write();

    pf_uperpen_puppipt390400->SetLineColor(1);
    pf_uperpen_puppipt390400->Write();

    pf_uperpen_puppipt420430->SetLineColor(1);
    pf_uperpen_puppipt420430->Write();


    pf_uparaqt_puppipt100110->SetLineColor(1);
    pf_uparaqt_puppipt100110->Write();

    pf_uparaqt_puppipt110120->SetLineColor(1);
    pf_uparaqt_puppipt110120->Write();

    pf_uparaqt_puppipt120130->SetLineColor(1);
    pf_uparaqt_puppipt120130->Write();

    pf_uparaqt_puppipt130140->SetLineColor(1);
    pf_uparaqt_puppipt130140->Write();

    pf_uparaqt_puppipt140150->SetLineColor(1);
    pf_uparaqt_puppipt140150->Write();

    pf_uparaqt_puppipt150160->SetLineColor(1);
    pf_uparaqt_puppipt150160->Write();

    pf_uparaqt_puppipt160170->SetLineColor(1);
    pf_uparaqt_puppipt160170->Write();

    pf_uparaqt_puppipt170180->SetLineColor(1);
    pf_uparaqt_puppipt170180->Write();

    pf_uparaqt_puppipt180190->SetLineColor(1);
    pf_uparaqt_puppipt180190->Write();

    pf_uparaqt_puppipt190200->SetLineColor(1);
    pf_uparaqt_puppipt190200->Write();

    pf_uparaqt_puppipt200210->SetLineColor(1);
    pf_uparaqt_puppipt200210->Write();

    pf_uparaqt_puppipt230240->SetLineColor(1);
    pf_uparaqt_puppipt230240->Write();

    pf_uparaqt_puppipt260270->SetLineColor(1);
    pf_uparaqt_puppipt260270->Write();

    pf_uparaqt_puppipt290300->SetLineColor(1);
    pf_uparaqt_puppipt290300->Write();

    pf_uparaqt_puppipt320330->SetLineColor(1);
    pf_uparaqt_puppipt320330->Write();

    pf_uparaqt_puppipt330340->SetLineColor(1);
    pf_uparaqt_puppipt330340->Write();

    pf_uparaqt_puppipt360370->SetLineColor(1);
    pf_uparaqt_puppipt360370->Write();

    pf_uparaqt_puppipt390400->SetLineColor(1);
    pf_uparaqt_puppipt390400->Write();

    pf_uparaqt_puppipt420430->SetLineColor(1);
    pf_uparaqt_puppipt420430->Write();

 




    //rawmet with only tight pho id
    rawmetphi->SetLineColor(1);
    rawmetphi->Write();

    rawmetpt->SetLineColor(1);
    rawmetpt->Write();

    rawmet_px->SetLineColor(1);
    rawmet_px->Write();

    rawmet_py->SetLineColor(1);
    rawmet_py->Write();

    uperpen_raw->SetLineColor(1);
    uperpen_raw->Write();

    uparall_raw->SetLineColor(1);
    uparall_raw->Write();

    raw_abs_scale->SetLineColor(1);
    raw_abs_scale->Write();

    raw_paraqt->SetLineColor(1);
    raw_paraqt->Write();

    




    //rawmet with one tight pho id, one jet, no loose ele and muon
    rawmetphi->SetLineColor(1);
    rawmetphi->Write();

    rawmetpt->SetLineColor(1);
    rawmetpt->Write();

    rawmet_px_cal->SetLineColor(1);
    rawmet_px_cal->Write();

    rawmet_py_cal->SetLineColor(1);
    rawmet_py_cal->Write();

    rawmet_px->SetLineColor(1);
    rawmet_px->Write();

    rawmet_py->SetLineColor(1);
    rawmet_py->Write();

    uperpen_raw->SetLineColor(1);
    uperpen_raw->Write();

    uparall_raw->SetLineColor(1);
    uparall_raw->Write();

    raw_abs_scale->SetLineColor(1);
    raw_abs_scale->Write();


    raw_nGoodvtx->SetLineColor(1);
    raw_nGoodvtx->Write();


    raw_abs_scale5060->SetLineColor(1);
    raw_abs_scale5060->Write();

    raw_abs_scale8090->SetLineColor(1);
    raw_abs_scale8090->Write();

    raw_abs_scale100110->SetLineColor(1);
    raw_abs_scale100110->Write();

    raw_abs_scale120130->SetLineColor(1);
    raw_abs_scale120130->Write();

    raw_abs_scale140150->SetLineColor(1);
    raw_abs_scale140150->Write();

    raw_abs_scale160170->SetLineColor(1);
    raw_abs_scale160170->Write();

    raw_abs_scale190200->SetLineColor(1);
    raw_abs_scale190200->Write();

    raw_abs_scale220230->SetLineColor(1);
    raw_abs_scale220230->Write();

    raw_abs_scale250260->SetLineColor(1);
    raw_abs_scale250260->Write();

    raw_abs_scaled260270->SetLineColor(1);
    raw_abs_scaled260270->Write();

    raw_abs_scale280290->SetLineColor(1);
    raw_abs_scale280290->Write();

    raw_abs_scale300310->SetLineColor(1);
    raw_abs_scale300310->Write();

    raw_abs_scale310320->SetLineColor(1);
    raw_abs_scale310320->Write();

    raw_abs_scale320330->SetLineColor(1);
    raw_abs_scale320330->Write();

    raw_abs_scale350360->SetLineColor(1);
    raw_abs_scale350360->Write();

    raw_abs_scale380390->SetLineColor(1);
    raw_abs_scale380390->Write();

    raw_abs_scale390400->SetLineColor(1);
    raw_abs_scale390400->Write();


    
    raw_upara5060->SetLineColor(1);
    raw_upara5060->Write();

    raw_upara8090->SetLineColor(1);
    raw_upara8090->Write();

    raw_upara100110->SetLineColor(1);
    raw_upara100110->Write();

    raw_upara120130->SetLineColor(1);
    raw_upara120130->Write();

    raw_upara140150->SetLineColor(1);
    raw_upara140150->Write();

    raw_upara160170->SetLineColor(1);
    raw_upara160170->Write();

    raw_upara190200->SetLineColor(1);
    raw_upara190200->Write();

    raw_upara220230->SetLineColor(1);
    raw_upara220230->Write();

    raw_upara250260->SetLineColor(1);
    raw_upara250260->Write();

    raw_upara260270->SetLineColor(1);
    raw_upara260270->Write();

    raw_upara280290->SetLineColor(1);
    raw_upara280290->Write();

    raw_upara300310->SetLineColor(1);
    raw_upara300310->Write();

    raw_upara310320->SetLineColor(1);
    raw_upara310320->Write();

    raw_upara320330->SetLineColor(1);
    raw_upara320330->Write();

    raw_upara350360->SetLineColor(1);
    raw_upara350360->Write();

    raw_upara380390->SetLineColor(1);
    raw_upara380390->Write();

    raw_upara390400->SetLineColor(1);
    raw_upara390400->Write();


    
    raw_uperpen5060->SetLineColor(1);
    raw_uperpen5060->Write();

    raw_uperpen8090->SetLineColor(1);
    raw_uperpen8090->Write();

    raw_uperpen100110->SetLineColor(1);
    raw_uperpen100110->Write();

    raw_uperpen120130->SetLineColor(1);
    raw_uperpen120130->Write();

    raw_uperpen140150->SetLineColor(1);
    raw_uperpen140150->Write();

    raw_uperpen160170->SetLineColor(1);
    raw_uperpen160170->Write();

    raw_uperpen190200->SetLineColor(1);
    raw_uperpen190200->Write();

    raw_uperpen220230->SetLineColor(1);
    raw_uperpen220230->Write();

    raw_uperpen250260->SetLineColor(1);
    raw_uperpen250260->Write();

    raw_uperpen260270->SetLineColor(1);
    raw_uperpen260270->Write();

    raw_uperpen280290->SetLineColor(1);
    raw_uperpen280290->Write();

    raw_uperpen300310->SetLineColor(1);
    raw_uperpen300310->Write();

    raw_uperpen310320->SetLineColor(1);
    raw_uperpen310320->Write();

    raw_uperpen320330->SetLineColor(1);
    raw_uperpen320330->Write();

    raw_uperpen350360->SetLineColor(1);
    raw_uperpen350360->Write();

    raw_uperpen380390->SetLineColor(1);
    raw_uperpen380390->Write();

    raw_uperpen390400->SetLineColor(1);
    raw_uperpen390400->Write();


    raw_uparaqt5060->SetLineColor(1);
    raw_uparaqt5060->Write();

    raw_uparaqt8090->SetLineColor(1);
    raw_uparaqt8090->Write();

    raw_uparaqt100110->SetLineColor(1);
    raw_uparaqt100110->Write();

    raw_uparaqt120130->SetLineColor(1);
    raw_uparaqt120130->Write();

    raw_uparaqt140150->SetLineColor(1);
    raw_uparaqt140150->Write();

    raw_uparaqt160170->SetLineColor(1);
    raw_uparaqt160170->Write();

    raw_uparaqt190200->SetLineColor(1);
    raw_uparaqt190200->Write();

    raw_uparaqt220230->SetLineColor(1);
    raw_uparaqt220230->Write();

    raw_uparaqt250260->SetLineColor(1);
    raw_uparaqt250260->Write();

    raw_uparaqt260270->SetLineColor(1);
    raw_uparaqt260270->Write();

    raw_uparaqt280290->SetLineColor(1);
    raw_uparaqt280290->Write();

    raw_uparaqt300310->SetLineColor(1);
    raw_uparaqt300310->Write();

    raw_uparaqt310320->SetLineColor(1);
    raw_uparaqt310320->Write();

    raw_uparaqt320330->SetLineColor(1);
    raw_uparaqt320330->Write();

    raw_uparaqt350360->SetLineColor(1);
    raw_uparaqt350360->Write();

    raw_uparaqt380390->SetLineColor(1);
    raw_uparaqt380390->Write();

    raw_uparaqt390400->SetLineColor(1);
    raw_uparaqt390400->Write();



    raw_uparanvtx04->SetLineColor(1);
    raw_uparanvtx04->Write();

    raw_uparanvtx46->SetLineColor(1);
    raw_uparanvtx46->Write();

    raw_uparanvtx68->SetLineColor(1);
    raw_uparanvtx68->Write();

    raw_uparanvtx812->SetLineColor(1);
    raw_uparanvtx812->Write();
    
    raw_uparanvtx1218->SetLineColor(1);
    raw_uparanvtx1218->Write();

    raw_uparanvtx1824->SetLineColor(1);
    raw_uparanvtx1824->Write();

    raw_uparanvtx2428->SetLineColor(1);
    raw_uparanvtx2428->Write();

    raw_uparanvtx2832->SetLineColor(1);
    raw_uparanvtx2832->Write();

    raw_uparanvtx3236->SetLineColor(1);
    raw_uparanvtx3236->Write();

    raw_uparanvtx3640->SetLineColor(1);
    raw_uparanvtx3640->Write();


    raw_uperpennvtx04->SetLineColor(1);
    raw_uperpennvtx04->Write();

    raw_uperpennvtx46->SetLineColor(1);
    raw_uperpennvtx46->Write();

    raw_uperpennvtx68->SetLineColor(1);
    raw_uperpennvtx68->Write();

    raw_uperpennvtx812->SetLineColor(1);
    raw_uperpennvtx812->Write();
    
    raw_uperpennvtx1218->SetLineColor(1);
    raw_uperpennvtx1218->Write();

    raw_uperpennvtx1824->SetLineColor(1);
    raw_uperpennvtx1824->Write();

    raw_uperpennvtx2428->SetLineColor(1);
    raw_uperpennvtx2428->Write();

    raw_uperpennvtx2832->SetLineColor(1);
    raw_uperpennvtx2832->Write();

    raw_uperpennvtx3236->SetLineColor(1);
    raw_uperpennvtx3236->Write();

    raw_uperpennvtx3640->SetLineColor(1);
    raw_uperpennvtx3640->Write();



    raw_abs_scale_puppipt100110->SetLineColor(1);
    raw_abs_scale_puppipt100110->Write();

    raw_abs_scale_puppipt110120->SetLineColor(1);
    raw_abs_scale_puppipt110120->Write();

    raw_abs_scale_puppipt120130->SetLineColor(1);
    raw_abs_scale_puppipt120130->Write();

    raw_abs_scale_puppipt130140->SetLineColor(1);
    raw_abs_scale_puppipt130140->Write();

    raw_abs_scale_puppipt140150->SetLineColor(1);
    raw_abs_scale_puppipt140150->Write();

    raw_abs_scale_puppipt150160->SetLineColor(1);
    raw_abs_scale_puppipt150160->Write();

    raw_abs_scale_puppipt160170->SetLineColor(1);
    raw_abs_scale_puppipt160170->Write();

    raw_abs_scale_puppipt170180->SetLineColor(1);
    raw_abs_scale_puppipt170180->Write();

    raw_abs_scale_puppipt180190->SetLineColor(1);
    raw_abs_scale_puppipt180190->Write();

    raw_abs_scale_puppipt190200->SetLineColor(1);
    raw_abs_scale_puppipt190200->Write();

    raw_abs_scale_puppipt200210->SetLineColor(1);
    raw_abs_scale_puppipt200210->Write();

    raw_abs_scale_puppipt230240->SetLineColor(1);
    raw_abs_scale_puppipt230240->Write();

    raw_abs_scale_puppipt260270->SetLineColor(1);
    raw_abs_scale_puppipt260270->Write();

    raw_abs_scale_puppipt290300->SetLineColor(1);
    raw_abs_scale_puppipt290300->Write();

    raw_abs_scale_puppipt320330->SetLineColor(1);
    raw_abs_scale_puppipt320330->Write();

    raw_abs_scale_puppipt330340->SetLineColor(1);
    raw_abs_scale_puppipt330340->Write();

    raw_abs_scale_puppipt360370->SetLineColor(1);
    raw_abs_scale_puppipt360370->Write();

    raw_abs_scale_puppipt390400->SetLineColor(1);
    raw_abs_scale_puppipt390400->Write();

    raw_abs_scale_puppipt420430->SetLineColor(1);
    raw_abs_scale_puppipt420430->Write();



    raw_upara_puppipt100110->SetLineColor(1);
    raw_upara_puppipt100110->Write();

    raw_upara_puppipt110120->SetLineColor(1);
    raw_upara_puppipt110120->Write();

    raw_upara_puppipt120130->SetLineColor(1);
    raw_upara_puppipt120130->Write();

    raw_upara_puppipt130140->SetLineColor(1);
    raw_upara_puppipt130140->Write();

    raw_upara_puppipt140150->SetLineColor(1);
    raw_upara_puppipt140150->Write();

    raw_upara_puppipt150160->SetLineColor(1);
    raw_upara_puppipt150160->Write();

    raw_upara_puppipt160170->SetLineColor(1);
    raw_upara_puppipt160170->Write();

    raw_upara_puppipt170180->SetLineColor(1);
    raw_upara_puppipt170180->Write();

    raw_upara_puppipt180190->SetLineColor(1);
    raw_upara_puppipt180190->Write();

    raw_upara_puppipt190200->SetLineColor(1);
    raw_upara_puppipt190200->Write();

    raw_upara_puppipt200210->SetLineColor(1);
    raw_upara_puppipt200210->Write();

    raw_upara_puppipt230240->SetLineColor(1);
    raw_upara_puppipt230240->Write();

    raw_upara_puppipt260270->SetLineColor(1);
    raw_upara_puppipt260270->Write();

    raw_upara_puppipt290300->SetLineColor(1);
    raw_upara_puppipt290300->Write();

    raw_upara_puppipt320330->SetLineColor(1);
    raw_upara_puppipt320330->Write();

    raw_upara_puppipt330340->SetLineColor(1);
    raw_upara_puppipt330340->Write();

    raw_upara_puppipt360370->SetLineColor(1);
    raw_upara_puppipt360370->Write();

    raw_upara_puppipt390400->SetLineColor(1);
    raw_upara_puppipt390400->Write();

    raw_upara_puppipt420430->SetLineColor(1);
    raw_upara_puppipt420430->Write();


    raw_uperpen_puppipt100110->SetLineColor(1);
    raw_uperpen_puppipt100110->Write();

    raw_uperpen_puppipt110120->SetLineColor(1);
    raw_uperpen_puppipt110120->Write();

    raw_uperpen_puppipt120130->SetLineColor(1);
    raw_uperpen_puppipt120130->Write();

    raw_uperpen_puppipt130140->SetLineColor(1);
    raw_uperpen_puppipt130140->Write();

    raw_uperpen_puppipt140150->SetLineColor(1);
    raw_uperpen_puppipt140150->Write();

    raw_uperpen_puppipt150160->SetLineColor(1);
    raw_uperpen_puppipt150160->Write();

    raw_uperpen_puppipt160170->SetLineColor(1);
    raw_uperpen_puppipt160170->Write();

    raw_uperpen_puppipt170180->SetLineColor(1);
    raw_uperpen_puppipt170180->Write();

    raw_uperpen_puppipt180190->SetLineColor(1);
    raw_uperpen_puppipt180190->Write();

    raw_uperpen_puppipt190200->SetLineColor(1);
    raw_uperpen_puppipt190200->Write();

    raw_uperpen_puppipt200210->SetLineColor(1);
    raw_uperpen_puppipt200210->Write();

    raw_uperpen_puppipt230240->SetLineColor(1);
    raw_uperpen_puppipt230240->Write();

    raw_uperpen_puppipt260270->SetLineColor(1);
    raw_uperpen_puppipt260270->Write();

    raw_uperpen_puppipt290300->SetLineColor(1);
    raw_uperpen_puppipt290300->Write();

    raw_uperpen_puppipt320330->SetLineColor(1);
    raw_uperpen_puppipt320330->Write();

    raw_uperpen_puppipt330340->SetLineColor(1);
    raw_uperpen_puppipt330340->Write();

    raw_uperpen_puppipt360370->SetLineColor(1);
    raw_uperpen_puppipt360370->Write();

    raw_uperpen_puppipt390400->SetLineColor(1);
    raw_uperpen_puppipt390400->Write();

    raw_uperpen_puppipt420430->SetLineColor(1);
    raw_uperpen_puppipt420430->Write();


    raw_uparaqt_puppipt100110->SetLineColor(1);
    raw_uparaqt_puppipt100110->Write();

    raw_uparaqt_puppipt110120->SetLineColor(1);
    raw_uparaqt_puppipt110120->Write();

    raw_uparaqt_puppipt120130->SetLineColor(1);
    raw_uparaqt_puppipt120130->Write();

    raw_uparaqt_puppipt130140->SetLineColor(1);
    raw_uparaqt_puppipt130140->Write();

    raw_uparaqt_puppipt140150->SetLineColor(1);
    raw_uparaqt_puppipt140150->Write();

    raw_uparaqt_puppipt150160->SetLineColor(1);
    raw_uparaqt_puppipt150160->Write();

    raw_uparaqt_puppipt160170->SetLineColor(1);
    raw_uparaqt_puppipt160170->Write();

    raw_uparaqt_puppipt170180->SetLineColor(1);
    raw_uparaqt_puppipt170180->Write();

    raw_uparaqt_puppipt180190->SetLineColor(1);
    raw_uparaqt_puppipt180190->Write();

    raw_uparaqt_puppipt190200->SetLineColor(1);
    raw_uparaqt_puppipt190200->Write();

    raw_uparaqt_puppipt200210->SetLineColor(1);
    raw_uparaqt_puppipt200210->Write();

    raw_uparaqt_puppipt230240->SetLineColor(1);
    raw_uparaqt_puppipt230240->Write();

    raw_uparaqt_puppipt260270->SetLineColor(1);
    raw_uparaqt_puppipt260270->Write();

    raw_uparaqt_puppipt290300->SetLineColor(1);
    raw_uparaqt_puppipt290300->Write();

    raw_uparaqt_puppipt320330->SetLineColor(1);
    raw_uparaqt_puppipt320330->Write();

    raw_uparaqt_puppipt330340->SetLineColor(1);
    raw_uparaqt_puppipt330340->Write();

    raw_uparaqt_puppipt360370->SetLineColor(1);
    raw_uparaqt_puppipt360370->Write();

    raw_uparaqt_puppipt390400->SetLineColor(1);
    raw_uparaqt_puppipt390400->Write();

    raw_uparaqt_puppipt420430->SetLineColor(1);
    raw_uparaqt_puppipt420430->Write();

    


    raw_paraqt->SetLineColor(1);
    raw_paraqt->Write();

    rawmet_tightphopt->SetLineColor(1);
    rawmet_tightphopt->Write();

    rawpho_pxdist->SetLineColor(1);
    rawpho_pxdist->Write();

    rawpho_pydist->SetLineColor(1);
    rawpho_pydist->Write();

    raw_response->SetLineColor(1);
    raw_response->Write();


    //--------------METs comparison
    pf_vs_puppi100->SetOption("COLZ");
	pf_vs_puppi100->Write();

    pf_vs_rawpuppi100->SetOption("COLZ");
	pf_vs_rawpuppi100->Write();
    
    raw_vs_puppi100->SetOption("COLZ");
	raw_vs_puppi100->Write();

    raw_vs_rawpuppi100->SetOption("COLZ");
	raw_vs_rawpuppi100->Write();

    rawpuppi_vs_puppi100->SetOption("COLZ");
	rawpuppi_vs_puppi100->Write();

    pf_vs_raw100->SetOption("COLZ");
	pf_vs_raw100->Write();


    pf_vs_puppi200->SetOption("COLZ");
	pf_vs_puppi200->Write();

    pf_vs_rawpuppi200->SetOption("COLZ");
	pf_vs_rawpuppi200->Write();
    
    raw_vs_puppi200->SetOption("COLZ");
	raw_vs_puppi200->Write();

    raw_vs_rawpuppi200->SetOption("COLZ");
	raw_vs_rawpuppi200->Write();

    rawpuppi_vs_puppi200->SetOption("COLZ");
	rawpuppi_vs_puppi200->Write();

    pf_vs_raw200->SetOption("COLZ");
	pf_vs_raw200->Write();




    //-------------------------test area---------------------------
    nPho_one->SetLineColor(1);
    nPho_one->Write();

    nJet_one->SetLineColor(1);
    nJet_one->Write();

    nEle_one->SetLineColor(1);
    nEle_one->Write();

    nMu_one->SetLineColor(1);
    nMu_one->Write();


    cout << "done filling histogram!" << endl;
    cout << "analysis completed!" << endl;
    cout << "output root file: " << "jetmet_tightpho_GJET200to400.root" << endl;

    auto timenow_end = chrono::system_clock::to_time_t(chrono::system_clock::now());

    cout << "start time: " << ctime(&timenow_start) << " (CT)" << endl;

	cout << "End time: " << ctime(&timenow_end) << " (CT)" << endl;

    file1.close();
    file2.close();
    file3.close();
    file4.close();
    file5.close();


}


//see _june2020.C for previous works

//07.24.2020 "talk on June 25 2020"

//remove _all and the test code (with only tight ID)

//adding RAWPUPPIMET for everything

//adding events with PUPPI > 100 GeV and PUPPI > 200GeV

/*response (scale) for PUPPI MET, PF MET, Raw PUPPI MET and 
Raw PF MET as a function of the PUPPI MET for 
events with PUPPI MET > 100 GeV and similarly for the resolution.*/

//resolution vs nGoodvtx => Resolution of uparallel and uperpendicular vs nVtx

//Resolution of uparallel and uperpendicular vs SumEt (I do have sumEt for each MET stored) -> next time

//make new C and h.

/*2D plots that you had shown before for PF MET vs
PUPPI MET or GenMET vs PUPPI MET  etc for Gamma+Jets events with PUPPI
MET > 100 GeV and also > 200 GeV.*/


//TO-DO (move up following if finished): -> all done 07.29.2020









#ifndef H2Zy_H2ZyAnalysis_H
#define H2Zy_H2ZyAnalysis_H

#include <EventLoop/Worker.h>
#include "HGamAnalysisFramework/HGamCommon.h"
#include "HGamAnalysisFramework/HgammaAnalysis.h"

#include "H2Zy/TruthSelect.h"
#include "H2Zy/EventFill.h"
#include "H2Zy/Variables.h"
#include "H2Zy/Object_llg.h"
#include "H2Zy/Container_llg.h"
#include "H2Zy/HZgammaHelper.h"
#include "H2Zy/MCLumi.h"
#include "H2Zy/OverlapRemovalTree.h"

#include "xAODEgamma/PhotonFwd.h"
#include "xAODEgamma/ElectronFwd.h"

#include "LHAPDF/PDF.h"

#ifndef __CINT__
#include "xAODMuon/Muon.h"
#endif 

#include <vector>
#include <set>
#include <utility>

//track handler header
#include "H2Zy/TrackHandler.h"
#include "H2Zy/MergedElectronID.h"
#include "H2Zy/MergedElectronID_v2.h"
#include "H2Zy/MergedElectronID_v2F.h"
#include "H2Zy/MergedElectronID_v3.h"

#include "IsolationSelection/IsolationCloseByCorrectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"


class H2ZyAnalysis : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)

private:
#ifndef __CINT__
  HG::TrackHandler *m_trackHandler; //! 
  HG::MergedElectronID * m_mergedElectronID; //!
  HG::MergedElectronID_v2 * m_mergedElectronID_v2; //!
  HG::MergedElectronID_v2F * m_mergedElectronID_v2F; //!
  HG::MergedElectronID_v3 * m_mergedElectronID_v3; //!
#endif // __CINT__ 
    
protected:
  inline virtual HG::TrackHandler *trackHandler() { return m_trackHandler; }


private:
  
  enum CutEnum { xAOD=0, DxAOD=1,ALLEVTS=2,
		 ALLEVTS_NOPU=3, GRL=4, PV=5, EVENT_QUALITY=6, Triggers=7,  Initial_sel=8,  
		 twolepton=9, mll_threshold=10, twolepton_onephoton=11, pre_sel=12, llgcut=13, TRIGGER_MATCH=14,
		 mllcut=15, Ph_ptcut1=16, Ph_ID=17, Ph_ISO=18, MASSCUT=19, Ph_ptcut2=20, PASSALL=21 };
  
  

  /*
 enum CutEnum { xAOD=0, DxAOD=1,ALLEVTS=2,
                 ALLEVTS_NOPU=3, GRL=4, PV=5, EVENT_QUALITY=6, Triggers=7,  Initial_sel=8, z_boson_assigment=9,
                 twolepton=10, mll_threshold=11, twolepton_onephoton=12, pre_sel=13, llgcut=14, TRIGGER_MATCH=15,
                 mllcut=16, Ph_ptcut1=17, Ph_ID=18, Ph_ISO=19, MASSCUT=20, Ph_ptcut2=21, PASSALL=22 };
*/

  CutEnum m_cutFlow; //!
  int istep_event ; //!
  int istep_weight ; //!
  xAOD::PhotonContainer photons; //!
  xAOD::ElectronContainer electrons; //!
  xAOD::ElectronContainer electrons_isoCorrTool;//fabio
  xAOD::MuonContainer muons0 ; //!
  xAOD::MuonContainer muons_isoCorrTool;//fabio 
  xAOD::MuonContainer muons ; //!
  xAOD::MuonContainer muons_corr; //fla
  xAOD::JetContainer jets   ; //!
  xAOD::JetContainer loose_jets   ; //!
  xAOD::MissingETContainer met;//!
  xAOD::MissingETContainer met2;//!
  xAOD::PhotonContainer all_correctedphotons; //!
  xAOD::ElectronContainer all_correctedelectrons; //!
  xAOD::ElectronContainer all_recoEl_test; //!


  //fla
  xAOD::TrackParticleContainer m_allTracks; //!
  xAOD::TrackParticleContainer m_preSelTracks; //!
  xAOD::TrackParticleContainer m_selTracks; //!
  xAOD::ElectronContainer m_allElectrons; //! 
  


  // output trees
  TTree *m_outputTree;      //!
  TTree *m_summaryTree;     //!
  TTree *m_truthTree;       //!
  
  // output histo with sum of weight info
  TH1D *m_sumW;            //!

  std::map<TString, TTree*> map_outputTree;          //!

  HZG::OverlapRemovalTree *m_ParticlesTreeBeforeOR; //!
  HZG::OverlapRemovalTree *m_ParticlesTreeAfterOR;  //!

#ifndef __CINT__
  std::vector< xAOD::Muon* > m_selMuons;  //!
#endif // __CINT__    
  std::vector< xAOD::Electron* > m_selElectrons;  //!
  std::vector< xAOD::Photon* > m_selPhotons;      //!

  HZG::TruthSelect m_truthselector; //!
  HZG::EventFill m_eventfill; //!

  // bool to keep track of correct initialization order
  bool m_initialized;          //!

  HG::StrV m_list_of_requiredTriggers; //!

  // information on initial and final events and sum of weights for DAODs
  double m_finalSumOfWeights;   //!
  double m_initialSumOfWeights; //!
  double m_finalSumOfWeightsSquared;   //!
  double m_initialSumOfWeightsSquared; //!
  double m_finalEvents;         //!
  double m_initialEvents;       //!
  double m_initialWeight;       //!
  double m_mcweight;            //!

  bool m_newFile; //!
  unsigned int nEventsProcessed ;  //!
  double sumOfWeights           ;  //!
  double sumOfWeightsSquared    ;  //!

  unsigned int nEventsDxAOD       ;  //!
  double sumOfWeightsDxAOD        ;  //!
  double sumOfWeightsSquaredDxAOD ;  //!

  TH1F* m_histo_cutflow;                        //!  
  TH1F* m_histo_cutflow_wt;                     //!  
  
  TH1F* m_histo_lep1_caloIso = 0;
  TH1F* m_histo_lep2_caloIso = 0;

  TH1F* m_histo_lep1_neflowisol20 = 0;
  TH1F* m_histo_lep2_neflowisol20 = 0;

  TH1F* m_histo_lep1_caloIsoEnergy = 0;
  TH1F* m_histo_lep2_caloIsoEnergy = 0;  

  TH1F* m_histo_lep2_caloIsoEnergy_lowTail = 0;

  TH1F* m_histo_l1_OnLeadLeptORnot = 0;
  TH1F* m_histo_l2_OnLeadLeptORnot = 0;
  TH1F* m_histo_nLeptons = 0; 

  //dR muons default
  TH1F* m_histo_dR_muons_default = 0;
  //dR muons corr
  TH1F* m_histo_dR_muons_corr = 0;


  TH1F* m_histo_lep1_trackIso = 0;
  TH1F* m_histo_lep2_trackIso = 0;

  TH1F* m_histo_lep1_caloIso_pass2lep = 0;
  TH1F* m_histo_lep2_caloIso_pass2lep = 0;

  TH1F* m_histo_lep1_trackIso_pass2lep = 0;
  TH1F* m_histo_lep2_trackIso_pass2lep = 0;

  TH1F* m_histo_lep1_neflowisol20_pass2lep = 0;
  TH1F* m_histo_lep2_neflowisol20_pass2lep = 0;

  TH1F* m_dilepton_mass_passTriggers = 0;
  TH1F* m_dilepton_mass_2leptonSelection = 0; 

  TH1F* m_histo_Nelectrons = 0;
  TH1F* m_histo_Nelectrons_corr = 0;

  TH1F* m_histo_Nelectrons_after = 0;
  TH1F* m_histo_Nelectrons_corr_after = 0;

//EG variables check
TH1F* m_clusterE_lep1_passTriggers = 0;
TH1F* m_clusterE_lep2_passTriggers = 0;
TH1F* m_clusterE_dR_passTriggers = 0;
TH1F* m_MassclusterE_passTriggers = 0;
TH1F* m_lep1Trk_passTriggers = 0;
TH1F* m_lep2Trk_passTriggers = 0;
TH1F* m_lep1_ccl_E0_passTriggers = 0;
TH1F* m_lep1_ccl_E1_passTriggers = 0;
TH1F* m_lep1_ccl_E2_passTriggers = 0;
TH1F* m_lep1_ccl_E3_passTriggers = 0;
TH1F* m_lep2_ccl_E0_passTriggers = 0;
TH1F* m_lep2_ccl_E1_passTriggers = 0;
TH1F* m_lep2_ccl_E2_passTriggers = 0;
TH1F* m_lep2_ccl_E3_passTriggers = 0;

TH1F* m_clusterE_lep1_2Lepton = 0;
TH1F* m_clusterE_lep2_2Lepton = 0;
TH1F* m_clusterE_dR_2Lepton = 0;
TH1F* m_MassclusterE_2Lepton = 0;
TH1F* m_lep1Trk_2Lepton = 0;
TH1F* m_lep2Trk_2Lepton = 0;
TH1F* m_lep1_ccl_E0_2Lepton = 0;
TH1F* m_lep1_ccl_E1_2Lepton = 0;
TH1F* m_lep1_ccl_E2_2Lepton = 0;
TH1F* m_lep1_ccl_E3_2Lepton = 0;
TH1F* m_lep2_ccl_E0_2Lepton = 0;
TH1F* m_lep2_ccl_E1_2Lepton = 0;
TH1F* m_lep2_ccl_E2_2Lepton = 0;
TH1F* m_lep2_ccl_E3_2Lepton = 0;

TH1F* m_iso_el_allElectrons = 0;
TH1F* m_iso_el_Electrons = 0;

 TH1F* m_l1_reco = 0;
 TH1F* m_l2_reco = 0;  

 TH1F* m_l1_truth = 0; 
 TH1F* m_l2_truth = 0; 

 TH1F* m_sum_l1_l2_reco = 0;

 TH1F* m_Nelectrons_reco = 0;
 TH1F* m_Ntracks_reco = 0;  

 TH1I* m_channels = 0;
 TH2F* merged_channel = 0;
 TH1F* m_eOverP = 0;
 TH1F* m_Rhad = 0;
 TH1F* m_Eratio = 0;
 TH1F* m_Rphi = 0;
 TH1F* m_Reta = 0;
 TH1F* m_f3 = 0;
 TH1F* m_wEta2 = 0;
 TH1F* m_wToTS1 = 0;
 TH1F* m_trk_TRT_PID1 = 0;
 TH1F* m_trk_TRT_PID2 = 0;
 TH1F* m_trk_dEta1 = 0;

 TH1F* m_merged_channel_el_pt = 0;
 TH1F* m_merged_channel_el_eta = 0;
 TH1F* m_merged_channel_el_phi = 0;

 TH1F* m_merged_channel_trk1_pt = 0;
 TH1F* m_merged_channel_trk2_pt = 0;

 TH1F* m_merged_channel_trk1_eta = 0;
 TH1F* m_merged_channel_trk2_eta = 0;
 TH1F* m_merged_channel_mtrktrk = 0;

 TH1F* m_merged_channel_el_pt_v2 = 0;


  TH1F* m_elSize = 0;
  TH1F* m_histo_dR = 0;
  TH1F* m_elSize_After2ndORL = 0;
  TH1F* m_histo_dR_After2ndORL = 0; 
  TH2F* NpreEl_vs_nTracks = 0;

  TH1F* m_elSize_ZeroElectronsEvents = 0;
  TH1F* m_elSize_ZeroElectronsEvents_leadingpT = 0;
  TH1F* m_elSize_ZeroElectronsEvents_subleadingpT = 0;
  TH1F* m_elSize_ZeroElectronsEvents_leadingEta = 0;
  TH1F* m_elSize_ZeroElectronsEvents_subleadingEta = 0;
  TH1F* m_elSize_ZeroElectronsEvents_passOQCut = 0;
  TH1F* m_elSize_ZeroElectronsEvents_passPtEtaCuts = 0;
  TH1F* m_elSize_ZeroElectronsEvents_passIPCuts = 0;
  TH1F* m_elSize_ZeroElectronsEvents_passHVCut = 0;
  TH1F* m_elSize_ZeroElectronsEvents_passPIDCut = 0;


  TH1F* m_HumElectronEvents_TruthleptonpT = 0;
  TH1F* m_HumElectronEvents_TruthleptonEta = 0;
  TH1F* m_HumElectronEvents_TruthleptonPhi = 0;


  TH1F* m_HumElectronEvents_RecoleptonpT = 0;
  TH1F* m_HumElectronEvents_RecoleptonEta = 0;
  TH1F* m_HumElectronEvents_RecoleptonPhi = 0;
 
  TH1F* m_HumElectronEvents_TruthStatus = 0;
  TH1F* m_HumElectronEvents_TruthPDGID = 0;
  TH1F* m_HumElectronEvents_TruthChilds = 0;
  TH1F* m_HumElectronEvents_pTResolution = 0;

  TH1F* m_ElectronEvents_TruthLeptonLeadpT = 0;
  TH1F* m_ElectronEvents_TruthLeptonLeadEta = 0;
  TH1F* m_ElectronEvents_RecoLeptonLeadpT  = 0;
  TH1F* m_ElectronEvents_RecoLeptonLeadEta = 0;
  TH1F* m_ElectronEvents_TruthLeptonSubLeadpT = 0;
  TH1F* m_ElectronEvents_TruthLeptonSubLeadEta = 0;
  TH1F* m_ElectronEvents_RecoLeptonSubLeadpT = 0;
  TH1F* m_ElectronEvents_RecoLeptonSubLeadEta = 0;
  TH1F* m_ElectronEvents_pTResolution_leading = 0;
  TH1F* m_ElectronEvents_pTResolution_subleading = 0;

  TH1F* m_ElectronEvents_tracking_pT = 0;
  TH1F* m_ElectronEvents_trackingEta = 0;

  TH1F* m_elSize_rawCont = 0;
  TH1F* m_elSize_preSelCont = 0;  
  TH1F* m_allLeptons_RecoLeptonpT = 0;
  TH1F* m_allLeptons_RecoLeptonEta = 0;
  TH1F* m_allLeptons_TruthLeptonpT = 0;
  TH1F* m_allLeptons_TruthLeptonpT_ZIsFather = 0;
  TH1F* m_allLeptons_TruthLeptonEta = 0;
  TH1F* m_allLeptons_TruthLeptonEta_ZIsFather = 0;
  TH1F* m_allLeptons_TruthLeptonPDGID = 0;
  TH1F* m_allLeptons_NotTruthElectronPDGID = 0;
  TH1F* m_allLeptons_NotTruthElectronBrokenLink = 0;
  TH1F* m_allLeptons_TruthLeptonpT_ZIsFather_NonPrompt = 0;
  TH1F* m_allLeptons_TruthLeptonEta_ZIsFather_NonPrompt = 0;
  TH1F* m_allLeptons_TruthLeptonpT_NonPrompt_ZIsNotFather = 0;
  TH1F* m_allLeptons_ZtrueElectrons = 0;
  TH1F* m_allLeptons_truthZleadingpT = 0;
  TH1F* m_allLeptons_truthZSubleadingpT = 0;
  TH1F* m_Zleptons = 0;
  TH1F* m_truthZleptons_mass = 0;

  TH1F* m_noPrecut_truthZleadingpT = 0;
  TH1F* m_noPrecut_truthZSubleadingpT = 0;
  TH1F* m_withPrecut_truthZleadingpT = 0;
  TH1F* m_withPrecut_truthZSubleadingpT = 0;
  TH1F* m_recoleadingpT = 0;
  TH1F* m_recoSubleadingpT = 0;

  TH1F* m_subLead_JetsOnlyEvts_LinkedTruthParticle = 0;
  TH1F* m_subLead_JetsOnlyEvts_pTresolution = 0;

  TH1F* m_subLead_JetsOnlyEvts_truthSubleadpT = 0;
  TH1F* m_subLead_JetsOnlyEvts_matchedJetpT = 0;
  TH1F* m_subLead_JetsOnlyEvts_Mljet = 0;

  TH1F* dR_el_lead_matchedHisto = 0;
  TH1F* dR_el_sublead_matchedHisto = 0;
  TH1F* leadingPt_Matchedreco  = 0;
  TH1F* subleadingPt_Matchedreco = 0;
  TH1F* histo_pl1_res  = 0;
  TH1F* histo_pl2_res = 0;

  TH1F* leadingPt_truth = 0;
  TH1F* subleadingPt_truth = 0;

  TH1F* leadingPt_MatchedrecoPresel  = 0;
  TH1F* subleadingPt_MatchedrecoPresel = 0;

  TH1F* histo_pl1_res_Presel  = 0;
  TH1F* histo_pl2_res_Presel = 0;
 
  TH1F* recoMatched_mll = 0;
  TH1F* recoMatched_m_true_inv_ll = 0;

  TH1F* histo_deltaR_leadingAndSubReco_rawCont = 0;
  TH1F* histo_deltaR_leadingAndSubReco_preSelCont = 0;
  TH1F* histo_isolation_leadingReco_caloOnly_preSelCont = 0;
  TH1F* histo_isolation_leadingReco_trackOnly_preSelCont = 0;
  TH1F* histo_isolation_subleadingReco_caloOnly_preSelCont = 0;
  TH1F* histo_isolation_subleadingReco_trackOnly_preSelCont = 0;

  TH1F* histo_isolation_leadingReco_caloOnly_preSelCont_IsoCorr = 0;
  TH1F* histo_isolation_leadingReco_trackOnly_preSelCont_IsoCorr = 0;
  TH1F* histo_isolation_subleadingReco_caloOnly_preSelCont_IsoCorr = 0;
  TH1F* histo_isolation_subleadingReco_trackOnly_preSelCont_IsoCorr = 0;

  TH1F* histo_isolation_leadingReco_caloOnly_preSelCont_passIso = 0;
  TH1F* histo_isolation_leadingReco_trackOnly_preSelCont_passIso = 0;
  TH1F* histo_isolation_subleadingReco_caloOnly_preSelCont_passIso = 0;
  TH1F* histo_isolation_subleadingReco_trackOnly_preSelCont_passIso = 0;

  //check isolation order
  TH1F* histo_isolationOrder_caloOnly = 0;
  TH1F* histo_isolationOrder_trackOnly = 0;
  TH1F* histo_isolationOrder_caloOnly_corr = 0;
  TH1F* histo_isolationOrder_trackOnly_corr = 0;

  TH1F* truthLepPt = 0; TH1F* truthLepEta = 0; TH1F* truthLepPt_preCut = 0; TH1F* truthLepEta_preCut = 0; TH1F* recoLepPt = 0;  TH1F* recoLepEta = 0; TH1F* recoLead_Ptres = 0;
   TH1F* recoSubLead_Ptres = 0; TH1F* recoMatched_dR = 0; TH1F*  fcloose_calo = 0; TH1F*  fcloose_trk = 0; TH1F* trkOnly_fixed = 0; TH1F* trkOnly_var = 0; TH1F* passFCLoose_calo = 0; TH1F* passFCLoose_trk = 0; TH1F*  passTightTrkFixed  = 0; TH1F*  passTightTrkVar = 0;  TH1F* fcloose_calo_corr = 0;  TH1F* fcloose_trk_corr = 0; TH1F*  trkOnly_fixed_corr = 0;
  TH1F* trkOnly_var_corr = 0;  TH1F* passFCLoose_calo_corr = 0; TH1F*  passFCLoose_trk_corr = 0; TH1F* passTightTrkFixed_corr = 0; TH1F*  passTightTrkVar_corr = 0;

  TH1F* categoryLeptons_leadingTruth = 0;
  TH1F* categoryJets_leadingTruth  = 0;
  TH1F* categoryNoMatch_leadingTruth = 0;
  TH1F* categoryLeptonsAndJets_leadingTruth = 0;
  TH1F* categoryLeptons_SubleadingTruth = 0;
  TH1F* categoryJets_SubleadingTruth = 0;
  TH1F* categoryNoMatch_SubleadingTruth = 0;
  TH1F* categoryLeptonsAndJets_SubleadingTruth = 0;
  TH1F* jetMass = 0;
  TH1F* jetPt = 0;

  TH1F* dRTwoTrueLeptons_catGoesIntoJet = 0;
  TH1F* Ptleading_catGoesIntoJet = 0;
  TH1F* PtSubleading_catGoesIntoJet = 0;
  TH1F* dRl1Jet_catGoesIntoJet = 0;
  TH1F* dRl2Jet_catGoesIntoJet = 0;

  TH1F* leptonsGoesIntoLeptons_leading_resolution = 0; 
  TH1F* leptonsGoesIntoLeptons_subleading_resolution = 0;   
  TH1F* dR_twoTrueLeptons = 0;
  TH2F* truthLeading_2Dresolution = 0;
  TH2F* truthSubLeading_2Dresolution = 0;
  TH1F* leadGoesLeptons_cluster = 0;
  TH1F* subleadGoesLeptons_cluster = 0;
  TH1F* leadGoesLeptons_RecoCluster = 0;
  TH1F* bothTrueLeptons_cluster = 0;

  TH1F* trkFoundLeading_NoMatch = 0; 
  TH1F* trkFoundSubLeading_NoMatch = 0;
  TH1F* dR_diTrks = 0;

  TH1F* jetTrackLeading_GoesIntoJet = 0;
  TH1F* jetTrackSubLeading_GoesIntoJet = 0;
  TH1F* dR_jetTracks_bothTracks = 0;

  TH1F* leadingLinkedTruthParticle_GoesIntoLeptons_Pt = 0;
  TH1F* leadingLinkedTruthParticle_GoesIntoLeptons_Eta = 0;
  TH1F* leadingTruthParticle_GoesIntoLeptons_Pt = 0;
  TH1F* leadingTruthParticle_GoesIntoLeptons_Eta = 0;
  TH1F* mtrk1trk2 = 0;
  TH1F* h_1dResolution_NoMatch = 0;
  TH2F* h_2dResolution_NoMatch  = 0; 
  TH1F* h_sumOfTrue_NoMatch = 0;

  TH1F* m_passingPrecut_truthZleadingpT  = 0;
  TH1F* m_passingPrecut_truthZSubleadingpT = 0;
  TH1F* dR_twoTrueLeptons_PassPrecuts = 0;


  //for truth leading 
  TH1F* hLead_matchViaFunction_recoMatched_Pt = 0;
  TH1F* hLead_matchLeptons_trueLead_Pt = 0;
  TH1F* hLead_matchLeptons_recoMatched_Pt = 0;
  TH1F* hLead_matchLeptons_1d_resolution = 0;
  TH2F* hLead_matchLeptons_2d_resolution = 0; 

  //for subleading
  TH1F* hSubLead_matchViaFunction_recoMatched_Pt = 0;
  TH1F* hSubLead_matchLeptons_trueSubLead_Pt = 0;
  TH1F* hSubLead_matchLeptons_recoMatched_Pt = 0;
  TH1F* hSubLead_matchLeptons_1d_resolution = 0;
  TH2F* hSubLead_matchLeptons_2d_resolution = 0;

  TH1F* hLead_matchLeptons_1d_InsideresolutionWindow = 0;
  TH1F* hSubLead_matchLeptons_1d_InsideresolutionWindow = 0;
  TH1F* hLeptons_matchLeptons_1d_BothInsideresolutionWindow = 0;
  TH1F* hLeptons_matchLeptons_Mll_BothInsideresolutionWindow = 0;

  //outside the resolution window
  TH1F* hLead_matchLeptons_outsideResWindow_Pt = 0;
  TH1F* hLead_matchLeptons_RecoEl_outsideResWindow_Cluster = 0;
  TH1F* hLead_matchLeptons_1d_resolution_outsideResWindow  = 0;
  TH2F* hLead_matchLeptons_2d_resolution_outsideResWindow = 0;
  TH1F* hSubLead_matchLeptons_outsideResWindow_Pt  = 0;
  TH1F* hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster = 0;
  TH1F* hSubLead_matchLeptons_1d_resolution_outsideResWindow = 0;
  TH2F* hSubLead_matchLeptons_2d_resolution_outsideResWindow = 0;
  TH1F* dRCluster_twoRecoElMatched_NoResContraint = 0;
  TH1F* dRCluster_twoRecoElMatched_BothInsideresolutionWindow = 0;

  TH1F* hLead_matchLeptons_1d_resolution_Above1p01 = 0;
  TH2F* hLead_matchLeptons_2d_resolution_Above1p01 = 0;
  TH1F* hSubLead_matchLeptons_1d_resolution_Above1p01 = 0;
  TH2F* hSubLead_matchLeptons_2d_resolution_Above1p01 = 0;

  TH1F* hLead_matchLeptons_1d_resolution_Below0p99 = 0;
  TH2F* hLead_matchLeptons_2d_resolution_Below0p99 = 0;
  TH1F* hSubLead_matchLeptons_1d_resolution_Below0p99 = 0;
  TH2F* hSubLead_matchLeptons_2d_resolution_Below0p99  = 0;

  //for double counting
  TH1F* hLead_matchLeptons_1d_resolution_doubleCount = 0;
  TH2F* hLead_matchLeptons_2d_resolution_doubleCount = 0;
  TH1F* hSubLead_matchLeptons_1d_resolution_doubleCount = 0;
  TH2F* hSubLead_matchLeptons_2d_resolution_doubleCount = 0;

  TH1F* hLead_matchLeptons_1d_resolution_Strip = 0;
  TH2F* hLead_matchLeptons_2d_resolution_Strip = 0;

  TH1F* m_dR_TrueLeptons = 0;
  TH1F* m_TrueLeptons_LeadingPt = 0;
  TH1F* m_TrueLeptons_SubLeadingPt = 0;
  TH1F* m_TrueLeptons_LeadingEta = 0;
  TH1F* m_TrueLeptons_SubLeadingEta = 0;
  TH1F* m_dR_TrueLeptons_PreCuts = 0;
  TH1F* m_TrueLeptons_LeadingPt_PreCuts = 0;
  TH1F* m_TrueLeptons_SubLeadingPt_PreCuts = 0;
  TH1F* m_TrueLeptons_LeadingEta_PreCuts = 0;
  TH1F* m_TrueLeptons_SubLeadingEta_PreCuts = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_TrueLeadingPt = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_TrueLeadingEta = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_RecoElMatchedPt = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_RecoElMatchedEta = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_1d_resolution = 0;
  TH2F* m_TrueLeadGoesIntoLeptons_2d_resolution = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_dR_l1AndReco = 0;
  
  TH1F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d  = 0;
 
  TH2F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d = 0;

  TH1F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR = 0;

  //std::map<TString,TH1F*> m_histoTH1F;          //!
  std::map<TString,TH1F*> m_cutflowhistoTH1F;          //!
  std::map<int, TString> m_event_cutflow_name;  //!
  std::map<int, TString> m_event_weight_name;   //!
  std::map<int, int> m_mcchannel_state;  //!
  TString histweight_name; //!
  TString hist_cutflow_name; //!

  // Trigger pass
  TH1F* m_histo_trigger_pass;                     //!
  TH1F* m_histo_trigger_pass_wt;                  //!
  std::map<int, TString> m_event_trigger_name;   //!
  TH1F* m_histo_trigger_denom;                   //!
  TH1F* m_histo_trigger_denom_wt;                //!

  // Trigger match
  TH1F* m_histo_trigger_pass_match;                     //!
  TH1F* m_histo_trigger_pass_match_wt;                  //!
  TH1F* m_histo_trigSF;                                 //!

  int m_set_outputtree; //!
  Dmap m_deventmap; //!
  Fmap m_feventmap; //!
  Imap m_ieventmap; //!
  UImap m_uieventmap; //!
  ULLmap m_ulleventmap; //!
  Bmap m_beventmap; //!
  Bmap m_btreemap; //!

  IVmap m_iveventmap;
  FVmap m_fveventmap;
  TLVVmap m_tlvveventmap;

  vector<TString> m_vbn;

  HZG::Container_llg V_llg; //!
  HZG::Object_llg *cand_llg; //!

  CP::IsolationSelectionTool           *m_isoTool_loose;  //!
  CP::IsolationSelectionTool           *m_isoTool_medium; //!
  CP::IsolationSelectionTool           *m_isoTool_tight;  //!
  CP::IsolationHelper* m_isoHelper; //!
  vector<xAOD::Iso::IsolationType> IsoORtypes; //!

  vector<LHAPDF::PDF*> m_pdfs; //!
  vector<LHAPDF::PDF*> m_pdfs_mstw;
  vector<LHAPDF::PDF*> m_pdfs_mmht;
  vector<LHAPDF::PDF*> m_pdfs_nn;

  //--- event info
  bool m_isDerivation;
  int m_channel;                         //! 1: eegamma; 2: mumugamma
  int m_hasDphoton;
  int m_hasttyphoton;
  unsigned int m_runNumber;              //! run number
  unsigned int m_eventNumber;            //! event number
  unsigned int m_lbn;                    //! lumi block number
  unsigned int m_nGoodCands;             //! the number of good candidates (llg triplets, photons, ...)

  //--- event selection flag --
  bool m_pass_grl;
  bool m_pass_pv;
  bool m_pass_quality;
  int m_nMuons;
  int m_nElectrons;
  int m_nPhotons;

  // - mc info
  unsigned int m_mc_channel_number;      //! mc channel number
  int m_mc_Year;                         //! mc Year
  float m_mc_weight_xs;                  //! cross-section*Br
  float m_mc_totevent;                   //! totevent
  float m_mc_weight_gen;                 //! generator weight
  float m_mc_weight_ph;                  //! photon scale factor
  float m_mc_weight_l1;                  //! lepton scale factors
  float m_mc_weight_l2;                  //! lepton scale factors

  // - data 
  unsigned int m_data_channel_number; //! data channel number

  // Apply GRL
  bool m_checkGRL;

  // trigger info
  bool m_trigger_passed;			//!
  unsigned int m_trigger_passed_items;		//!
  unsigned int m_trigger_matched_items;		//!
  bool m_checkTrig;
  bool m_checkTrigMatch;

  // Overlap Removal
  bool m_doOverlapRemoval;

  bool m_MxAODinput;//!
  TString m_photonContainerName;//!
  TString m_jetContainerName;//! 
  TString m_elecContainerName; //!
  TString m_muonContainerName; //!
  TString m_evtInfoName; //!
  bool first; //! 
  
  // events, filter efficiency and luminosity of MC samples
  MCLumi* m_mclumi; //!

  bool m_theorySyst; //!
  HG::StrV m_expSyst; //!

public:
  // standard constructors and destructor
  H2ZyAnalysis() { }
  H2ZyAnalysis(const char *name);
  virtual ~H2ZyAnalysis();

  // functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode fileExecute();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode histFinalize ();
  virtual EL::StatusCode finalize();


  CutEnum cutflow(TString sysname);
  CutEnum initialcutflow();
  CutEnum precutflow();
  CutEnum finalcutflow();
  void fillCutFlow(CutEnum cut, TString sysname, double w);


  void SetupOutputTree( HZG::Container_llg &_V_llg, TString sysname="");
  void SetupSummaryTree();
  void FillVars(TString varname, HZG::Object_llg &_llg);
  void FillVars(vector<TString> varname, HZG::Object_llg &_llg);
  void ResetVars();

  void SetEventInfo();
  void RecSave(HZG::Object_llg *&_llg);
  void setcutflowname();

  //Trigger related variables & functions
  bool m_applyTrig;
  bool m_applyTrigMatch;

  void addBookKeeping();
  void createTriggerHists();
  void setupTrigger();
  bool fillTriggerInfo(int channel, float weight);
  void finalizeTrigger();
  void triggerEfficiencyPerItem();
  // VBF truth match

  // Fabio's
  int nEvents = 0;;
  double norm_yield = 0;
  int n_WeidEvts = 0;
  double sum_w = 0;
  int nTracksIgualHum = 0;
  int nTracksAboveHum = 0;
  int N1lepevts_recoHasNoTruthLinked = 0;
  int N2lepevtsLeading_recoHasNoTruthLinked = 0;
  int N2lepevtsSubLeading_recoHasNoTruthLinked = 0;
  int nTracksLower_pt = 0;
  int nTracksLower_Eta = 0;
  int n_totRecoPreSelElectrons = 0;
  int n_totTruePreSelElectrons = 0;
  int n_totTruePreSelElectrons_ZIsFather = 0;
  int evtsWithTrueLeadingLeptons = 0;
  int evtsWithTrueSubLeadingLeptons = 0;
  int nevts = 0;
  int n_trueElIsStable = 0;
  int n_truthElFatherIs11 = 0;
  int n_trueElNoParents = 0; 
  int nEvts_matchLeptons = 0;
  int nEvts_matchJets = 0;
  int nEvts_NoMatch = 0;
  int nEvts_matchLeptonsAndJets = 0;
 
  int nEvtsSub_matchLeptons = 0;
  int nEvtsSub_matchJets = 0;
  int nEvtsSub_NoMatch = 0;
  int nEvtsSub_matchLeptonsAndJets = 0;
  int nEvtsSub_IsmatchedLeptonOnly = 0;

  int nEvts_empty = 0;
  int nEvts_1 = 0;
  int nEvts_2 = 0;
  int nEvts_IsmatchedLeptonOnly = 0;
  int nEvts_Less10GeV = 0;
  int nEvts_HigherEta2p47 = 0;
  int evts_noPrecut = 0;
  int nevts_before = 0;
  int nevts_after_lead = 0;
  int nevts_after_sublead = 0;
  int nevts_atleast2Truth = 0;
  int nevts_atleast1Truth = 0;
  int nevts_nTrkJet = 0; 
  int nevts_LeadIsFoundsubLeadIntoJet = 0;
  int nevts_LeadIsFoundsubLeadNoMatch = 0;
  int n_BothTruthLeptonFound = 0;
  int n_1TruthLeptonFound = 0;
  int n_1TruthLeptonFoundAnd1JetFound = 0;
  int n_1TruthLeptonNoSubLead = 0;
  int nevtsBoth_before = 0;
  int nevtsAfter_before = 0;
  int n_BothTruthLeptonNotFound = 0;
  int nEvts_bothPreSelRecoFound = 0;
  int n_PreSelTracksIsOne = 0;
  int n_PreSelTracksIsLargerOne = 0;
  int n_EvtsDoubleCounting = 0;

  int n_LeadIntoLeptons_l2RecoMatch = 0;
  int n_SubLeadIntoLeptons_l1RecoMatch = 0;

  int n_evts_strip = 0;
  int n_LeadIntoLeptons_Strip_l2RecoMatch = 0;
  int n_evts_TrueSub_Out = 0;
  int n_evts_TrueSub_In = 0;
  int nEvts_strip_trueSubPassPreCut = 0;
  


  //fla - merged electron ID
  
  asg::AnaToolHandle<CP::IIsolationSelectionTool> m_iso_tool;
  asg::AnaToolHandle<CP::IIsolationCloseByCorrectionTool> m_isoCloseByFCLoose;
  void AddElectronDecorations(xAOD::ElectronContainer& electrons);
  std::map<HG::Iso::IsolationType, CP::IsolationCloseByCorrectionTool*> m_isoCloseByTools_Ele; //!
  std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_eleIsoAccCorr; //!
  HG::Iso::IsolationType m_eleResolvedIsoWP_FCLoose;
  HG::Iso::IsolationType m_eleResolvedIsoWP_TightTrkFixed;
  HG::Iso::IsolationType m_eleResolvedIsoWP_TightTrkVar;
  void  FindZboson_ElectronChannelAware(xAOD::TrackParticleContainer* inTracks,
                                                  xAOD::TrackParticle*& sel_trk1,
                                                  xAOD::TrackParticle*& sel_trk2,
                                                  double& return_mll,
                                                  const HG::TrackElectronMap& trkEleMap,
                                                  xAOD::ElectronContainer* inEleCont,
                                                  xAOD::ElectronContainer* outEleCont);

  
  vector<xAOD::Electron*> match_DR_electrons(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  vector<xAOD::Jet*> match_DR_jet(const xAOD::TruthParticle* true_lepton, xAOD::JetContainer recoJetCont);
  vector<xAOD::Electron*> match_DR_electrons_InFunction(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  //void match_DR_electrons_InFunction(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  xAOD::ElectronContainer skipRecoMatchedToTrueLeading(xAOD::Electron* el_matchedLeading, xAOD::ElectronContainer recoElCont);
  //check double counting
  bool  checkDoubleCountingEvents(xAOD::Electron* el_matchedLeading, const xAOD::TruthParticle* true_lepton);
  
  bool truePassPreCut(const xAOD::TruthParticle* true_lepton);

  TString m_eleIDPreselection;
  const xAOD::TruthParticle* recoTruthMatching(const xAOD::Electron *el);
  bool isTruthLepton(const xAOD::Electron *el);
  void decorateCorrectedIsoCut(xAOD::ElectronContainer & electrons);
  void CheckIsoToolsOrder(TH1F *h_isoCalo, TH1F *h_isoTrack, TH1F *h_isoCalo_IsoCorr, TH1F *h_isoTrack_IsoCorr);
  void ResolvedChannelEvts(TH1F *h_truthLepPt,TH1F *h_truthLepEta, TH1F *h_truthLepPt_preCut, TH1F *h_truthLepEta_preCut, TH1F *h_recoLepPt, TH1F *h_recoLepEta, TH1F *h_recoLead_Ptres, TH1F *h_recoSubLead_Ptres, TH1F *h_recoMatched_dR, TH1F *h_fcloose_calo, TH1F *h_fcloose_trk, TH1F *h_trkOnly_fixed, TH1F *h_trkOnly_var, TH1F *h_passFCLoose_calo, TH1F *h_passFCLoose_trk, TH1F *h_passTightTrkFixed, TH1F *h_passTightTrkVar, TH1F *h_fcloose_calo_corr, TH1F *h_fcloose_trk_corr, TH1F *h_trkOnly_fixed_corr, TH1F *h_trkOnly_var_corr, TH1F *h_passFCLoose_calo_corr, TH1F *h_passFCLoose_trk_corr, TH1F *h_passTightTrkFixed_corr, TH1F *h_passTightTrkVar_corr);
  // this is needed to distribute the algorithm to the workers
  bool m_applySystematicLoop;

  ClassDef(H2ZyAnalysis, 1);
};

#endif // Zy_H2ZyAnalysis_H

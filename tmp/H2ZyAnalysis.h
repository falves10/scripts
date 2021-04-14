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
  enum CutEnum { xAOD=0, DxAOD=1,ALLEVTS=2,
		 ALLEVTS_NOPU=3, GRL=4, PV=5, EVENT_QUALITY=6, Triggers=7,  Initial_sel=8, 
		 twolepton=9, mll_threshold=10, twolepton_onephoton=11, pre_sel=12, llgcut=13, TRIGGER_MATCH=14,
		 mllcut=15, Ph_ptcut1=16, Ph_ID=17, Ph_ISO=18, MASSCUT=19, Ph_ptcut2=20, PASSALL=21 };
  
  CutEnum m_cutFlow; //!
  int istep_event ; //!
  int istep_weight ; //!
  xAOD::PhotonContainer photons; //!
  xAOD::ElectronContainer electrons; //!
  xAOD::MuonContainer muons0 ; //!
  xAOD::MuonContainer muons ; //!
  xAOD::JetContainer jets   ; //!
  xAOD::JetContainer loose_jets   ; //!
  xAOD::MissingETContainer met;//!
  xAOD::MissingETContainer met2;//!
  xAOD::PhotonContainer all_correctedphotons; //!
  xAOD::ElectronContainer all_correctedelectrons; //!

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
  
  // Reco Studies
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
  
  TH1F* m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt = 0;
  TH1F* m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta = 0;
  TH1F* m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt = 0;
  TH1F* m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta = 0;
  TH1F* m_TrueSubLeadGoesIntoLeptons_1d_resolution = 0;
  TH2F* m_TrueSubLeadGoesIntoLeptons_2d_resolution = 0;
  
  TH1F* m_TrueLeadGoesIntoLeptons_dR_l1AndReco = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d  = 0;
  TH2F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d = 0;
  TH1F* m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR = 0;

  TH1F* m_el_truthMatched = 0;
  TH1F* m_GSFtrk_pT = 0;
  TH1F* m_GSFtrk_eta = 0;
  TH1F* m_IDtrk_pT = 0;
  TH1F* m_IDtrk_eta = 0;
  
  TH1F* m_dEta_trueLeptons = 0;
  TH1F* m_dPhi_trueLeptons = 0;
  TH1F* m_trueSubLead_pT_noMatched = 0;
  TH1F* m_dR_realLeptons = 0;
  

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
  virtual EL::StatusCode finalize();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode histFinalize ();

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
  void createRecoStudiesHists();//fla
  void setupTrigger();
  bool fillTriggerInfo(int channel, float weight);
  void finalizeTrigger();
  void triggerEfficiencyPerItem();
  // VBF truth match

  //Reco Studies variables
  int nEvents = 0;
  int nEvts_matchLeptons = 0;
  int nEvts_matchJets = 0;
  int nEvts_NoMatch = 0;
  int nEvts_matchLeptonsAndJets = 0;
 
  int nEvtsSub_matchLeptons = 0;
  int nEvtsSub_matchJets = 0;
  int nEvtsSub_NoMatch = 0;
  int nEvtsSub_matchLeptonsAndJets = 0;
  
  int nevts_leadAfterPrecut = 0;
  int nevts_SubleadAfterPrecut = 0;
  int n_evts_TrueSub_Out = 0;
  int n_evts_TrueSub_In = 0;
  int n_LeadIntoLeptons_l2RecoMatch = 0;
  int n_SubLeadIntoLeptons_l1RecoMatch = 0;
  int n_evts_strip = 0;
  int nEvts_strip_trueSubPassPreCut = 0;
  int nEvts_IsmatchedLeptonOnly = 0;
  int n_LeadIntoLeptons_Strip_l2RecoMatch = 0; 
  int nEvts_bothTrueMatchSameRecoEl_Strip = 0;
  int nEvts_trueLeadMatchARecoViaFunc = 0;
  int nEvts_trueSubLeadMatchARecoViaFunc = 0;
  int nEvts_trueSubLead_BrokenLink = 0;
  int nEvts_trueSubLead_FoundTrack = 0;
  int nEvts_trueSubLead_FoundTrack_IDTrk = 0;
  int nEvts_trksMatchSameEl = 0;
  int nEvts_zero = 0;
  int nEvts_hum = 0;
  int nEvts_dois = 0;
  
  
  //Functions used in the Reconstruction type studies
  vector<xAOD::Electron*> match_DR_electrons(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  vector<xAOD::Jet*> match_DR_jet(const xAOD::TruthParticle* true_lepton, xAOD::JetContainer recoJetCont);
  xAOD::ElectronContainer skipRecoMatchedToTrueLeading(xAOD::Electron* el_matchedLeading, xAOD::ElectronContainer recoElCont);
  bool checkDoubleCountingEvents(xAOD::Electron* el_matchedLeading, const xAOD::TruthParticle* true_lepton);
  bool isTruthLepton(const xAOD::Electron *el);
  bool truePassPreCut(const xAOD::TruthParticle* true_lepton);
  bool trueMatch_toolFunction(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  xAOD::ElectronContainer getRecoMatched( vector<const xAOD::TruthParticle*> lepton_true, xAOD::ElectronContainer recoElCont );
  bool trueMatch_BrokenLink(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  vector<const xAOD::TruthParticle*> trueMatch_noMatched(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont);
  bool findTracksMatchingSameElectron(const xAOD::TruthParticle* truelepton_l1, const xAOD::TruthParticle* truelepton_l2, xAOD::ElectronContainer recoElCont);
  const xAOD::IParticle* getTheOriginalPointer(const xAOD::IParticle& part);
  float getTruthMatchProbability(const xAOD::TrackParticle* trackParticle);
  xAOD::TrackParticle* foundIDTrack(const xAOD::TruthParticle* tlep);


  // this is needed to distribute the algorithm to the workers
  bool m_applySystematicLoop;

  ClassDef(H2ZyAnalysis, 1);
};

#endif // Zy_H2ZyAnalysis_H

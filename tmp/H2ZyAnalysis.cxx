//#include <sys/time.h>  
#include "H2Zy/H2ZyAnalysis.h"
#include "HGamAnalysisFramework/HGamCommon.h"
#include "HGamAnalysisFramework/HGamVariables.h"
#include "TTreeFormula.h"
#include "H2Zy/OverlapRemovalTree.h"
#include "LHAPDF/LHAPDF.h"
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/Reweighting.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "xAODEgamma/ElectronxAODHelpers.h"

// this is needed to distribute the algorithm to the workers
ClassImp(H2ZyAnalysis)



H2ZyAnalysis::H2ZyAnalysis(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  m_set_outputtree                     = 0;
  m_histo_cutflow                      = 0;
  m_histo_cutflow_wt                   = 0;
  m_histo_trigger_pass                 = 0;
  m_histo_trigger_pass_wt              = 0;
  m_histo_trigger_pass_match           = 0;
  m_histo_trigger_pass_match_wt        = 0;
  m_histo_trigger_denom                = 0;
  m_histo_trigger_denom_wt             = 0;
  m_histo_trigSF                       = 0;
  m_sumW                               = 0;
  m_checkGRL                           = false;
  m_checkTrig                          = false;
  m_checkTrigMatch                     = false;
  m_mclumi                             = 0;
  m_doOverlapRemoval                   = false;
  m_theorySyst                         = false;
}



H2ZyAnalysis::~H2ZyAnalysis()
{
  // Here you delete any memory you allocated during your analysis.
}


EL::StatusCode H2ZyAnalysis::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  //histoStore()->createTH1F("m_yy", 60, 110, 140);
  //m_eventHandler = new HG::EventHandler(event(), store());
  //m_eventHandler->initialize(config());
  //m_truthselector.initialize(*(config()),  event(), &m_beventmap, &m_ieventmap, &m_feventmap, &m_fveventmap, &m_tlvveventmap);


  // create tree containing truth-level information
  m_truthTree = m_truthselector.settruthtree(*(config()));
  wk()->addOutput (m_truthTree);
  
  //config output MxAOD collections
  m_photonContainerName = "HGam"+config()->getStr("PhotonHandler.ContainerName");
  m_jetContainerName    = "HGam"+config()->getStr("JetHandler.ContainerName");
  m_elecContainerName   = "HGam"+config()->getStr("ElectronHandler.ContainerName");
  m_muonContainerName   = "HGam"+config()->getStr("MuonHandler.ContainerName");
  m_evtInfoName         = "EventInfo"; // fix this

  if (config()->isDefined("MxAOD.Variables.Photon"))
	  event()->setAuxItemList( (m_photonContainerName+"Aux.").Data(), config()->getStr("MxAOD.Variables.Photon").Data() );
  if (config()->isDefined("MxAOD.Variables.Jet"))
	  event()->setAuxItemList( (m_jetContainerName+"Aux.").Data(), config()->getStr("MxAOD.Variables.Jet").Data() );
  if (config()->isDefined("MxAOD.Variables.Electron"))
	  event()->setAuxItemList( (m_elecContainerName+"Aux.").Data(), config()->getStr("MxAOD.Variables.Electron").Data() );
  if (config()->isDefined("MxAOD.Variables.Muon"))
	  event()->setAuxItemList( (m_muonContainerName+"Aux.").Data(), config()->getStr("MxAOD.Variables.Muon").Data() );
  if (config()->isDefined("MxAOD.Variables.EventInfo"))
	  event()->setAuxItemList( (m_evtInfoName+"Aux.").Data(), config()->getStr("MxAOD.Variables.EventInfo").Data() );
  //--------- 


  m_MxAODinput = false;
  first = true;
  m_set_outputtree=0;
  setcutflowname();

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode H2ZyAnalysis::histInitialize()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  HgammaAnalysis::histInitialize();

  m_finalSumOfWeights = 0;
  m_initialSumOfWeights = 0;
  m_finalSumOfWeightsSquared = 0;
  m_initialSumOfWeightsSquared = 0;
  m_finalEvents = 0;
  m_initialEvents = 0;

  // create and add to output small tree with info on MC xsection, filter efficieny, ...
//hyp  SetupSummaryTree();

  // create and add to output histogram for cutflow (unweighted)
  m_histo_cutflow = new TH1F("histo_cutflow", "histo_cutflow", 20, 0, 20.);
  wk()->addOutput(m_histo_cutflow);

  // create and add to output histogram for cutflow (weighted)
  m_histo_cutflow_wt = new TH1F("histo_cutflow_wt", "histo_cutflow_wt", 20, 0, 20.);
  wk()->addOutput(m_histo_cutflow_wt);
    
  

  // create and add to output trigger related histograms
  createTriggerHists();
  // create and add to output reco studies related histograms
  createRecoStudiesHists();

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode H2ZyAnalysis::fileExecute()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  // We use this function to extract initial and final sum of weights
  // and initial and final events from the metadata (for MC and data derivations)
  // For data AOD, we just sum the number of events in the file
  // For MC AOD, the sum of weights is computed in FillMCInfo

  // this function is called before initialize, unfortunately

  // file name
  Info("fileExecute", "Extracting bookkeeping information from file \"%s\"",wk()->inputFile()->GetName());
  Info("fileExecute", "This file has %lli entries", wk()->xaodEvent()->getEntries());

  HgammaAnalysis::fileExecute();  

  m_newFile=true;

  // read CutBookKeeper information from MetaData
  // recipe from https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AnalysisMetadata#Examples_MC15
  // get the MetaData tree once a new file is opened, with
  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!MetaData) {
    Error("fileExecute()", "MetaData not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  
  bool isAOD = MetaData->GetBranch("StreamAOD");
  bool isMAOD = !MetaData->GetBranch("TriggerMenu");
  HG::setAndLock_InputType(isAOD,isMAOD);

  MetaData->LoadTree(0);

  nEventsProcessed         = 0;
  sumOfWeights             = 0;
  sumOfWeightsSquared      = 0;
  nEventsDxAOD             = 0;
  sumOfWeightsDxAOD        = 0;
  sumOfWeightsSquaredDxAOD = 0;


  if (!HG::isAOD()){
	  bool m_bookkeeper = false;
	  if (config()->getBool("EventHandler.bookkeep", false)) m_bookkeeper = true;    

	  nEventsProcessed = 0;
	  sumOfWeights           = 0;
	  sumOfWeightsSquared    = 0;
	  nEventsDxAOD       = 0;
	  sumOfWeightsDxAOD        = 0;
	  sumOfWeightsSquaredDxAOD = 0;

	  if(m_bookkeeper){
		  // check for corruption
		  const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
		  if(!wk()->xaodEvent()->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
			  Error("fileExecute()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
			  return EL::StatusCode::FAILURE;
		  }
		  if ( incompleteCBC->size() != 0 ) {
			  //Error("fileExecute()","Found incomplete Bookkeepers! Check file for corruption.");
			  //return EL::StatusCode::FAILURE;
		  }

		  // Now, let's find the actual information
		  const xAOD::CutBookkeeperContainer* completeCBC = 0;
		  if(!wk()->xaodEvent()->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
			  Error("fileExecute()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
			  return EL::StatusCode::FAILURE;
		  }

		  
		  // First, let's find the smallest cycle number,
		  // i.e., the original first processing step/cycle
		  int minCycle = 10000;
		  for ( auto cbk : *completeCBC ) {
			  if ( ! cbk->name().empty()  && minCycle > cbk->cycle() ){ minCycle = cbk->cycle(); }
		  }
		  // Now, let's actually find the right one that contains all the needed info...
		  const xAOD::CutBookkeeper* allEventsCBK=0;
		  const xAOD::CutBookkeeper* DxAODEventsCBK=0;
		  std::string derivationName = "HIGG1D2Kernel";
		  int maxCycle = -1;
		  for (const auto& cbk: *completeCBC) {
		    if (cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD") {
		      allEventsCBK = cbk;
		      maxCycle = cbk->cycle();
		    }
		    if ( cbk->name() == derivationName){
		      DxAODEventsCBK = cbk;
		    }
		  }
		  if(!allEventsCBK) {
			  Error("fileExecute()","Failed to find AllExecutedEvents CutBookkeeper in MetaData! Exiting.");
			  return EL::StatusCode::FAILURE;
		  }
		  if(!DxAODEventsCBK) {
			  Error("fileExecute()",Form("Failed to find %s CutBookkeeper in MetaData! Exiting.",derivationName.c_str()));
			  return EL::StatusCode::FAILURE;
		  }

		  nEventsProcessed = allEventsCBK->nAcceptedEvents();
		  sumOfWeights           = allEventsCBK->sumOfEventWeights();
		  sumOfWeightsSquared    = allEventsCBK->sumOfEventWeightsSquared();

		  nEventsDxAOD       = DxAODEventsCBK->nAcceptedEvents();
		  sumOfWeightsDxAOD        = DxAODEventsCBK->sumOfEventWeights();
		  sumOfWeightsSquaredDxAOD = DxAODEventsCBK->sumOfEventWeightsSquared();
		  cout<<nEventsProcessed<<" "<<sumOfWeights<<" "<<sumOfWeightsSquared<<endl;
		  cout<<nEventsDxAOD<<" "<<sumOfWeightsDxAOD<<" "<<sumOfWeightsSquaredDxAOD<<endl;
	  }

  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode H2ZyAnalysis::changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  //Warning in <H2ZyAnalysis::changeInput()>: Could not find branch EventBookkeepers.m_nWeightedAcceptedEvents in MetaData
  Info("changeInput", "Zy inside changeinput");
  HgammaAnalysis::changeInput(true);

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode H2ZyAnalysis::finalize()
{
	
    
    cout<<" nEvts_matchLeptons " << nEvts_matchLeptons << endl;
    cout<<" nEvts_matchJets " << nEvts_matchJets << endl;
    cout<<" nEvts_NoMatch " << nEvts_NoMatch << endl;
    cout<<" nEvts_matchLeptonsAndJets " << nEvts_matchLeptonsAndJets << endl;
    
    cout<<" nEvtsSub_matchLeptons " << nEvtsSub_matchLeptons << endl;
    cout<<" nEvtsSub_matchJets " << nEvtsSub_matchJets << endl;
    cout<<" nEvtsSub_NoMatch " << nEvtsSub_NoMatch << endl;
    cout<<" nEvtsSub_matchLeptonsAndJets " << nEvtsSub_matchLeptonsAndJets << endl;
    cout<<" nevts_leadAfterPrecut " << nevts_leadAfterPrecut << endl;
    cout<<" nevts_SubleadAfterPrecut " << nevts_SubleadAfterPrecut << endl;
    
    cout<<" n_LeadIntoLeptons_l2RecoMatch " << n_LeadIntoLeptons_l2RecoMatch << endl;
    cout<<" n_SubLeadIntoLeptons_l1RecoMatch " << n_SubLeadIntoLeptons_l1RecoMatch << endl;
    cout<<" Events in the vertical strip "<< n_evts_strip << endl;
    cout<<" Events in the vertical strip (TSL matches same reco) "<< n_LeadIntoLeptons_Strip_l2RecoMatch << endl;
    cout<< " n_evts_TrueSub_Out " << n_evts_TrueSub_Out << endl;
    cout<< " n_evts_TrueSub_In " << n_evts_TrueSub_In << endl;
    cout<<" nEvts_strip_trueSubPassPreCut " << nEvts_strip_trueSubPassPreCut << endl;
    cout<<" nEvts_bothTrueMatchSameRecoEl_Strip " << nEvts_bothTrueMatchSameRecoEl_Strip << endl;
    cout<<" nEvts_trueLeadMatchARecoViaFunc " << nEvts_trueLeadMatchARecoViaFunc << endl;
    cout<<" nEvts_trueSubLeadMatchARecoViaFunc " << nEvts_trueSubLeadMatchARecoViaFunc << endl;
    cout<<" nEvts_trueSubLead_BrokenLink " << nEvts_trueSubLead_BrokenLink << endl;
    cout<<" nEvts_trueSubLead_FoundTrack " << nEvts_trueSubLead_FoundTrack << endl;
    cout<<" nEvts_trueSubLead_FoundTrack_IDTrk " << nEvts_trueSubLead_FoundTrack_IDTrk << endl;
    cout<<" nEvts_trksMatchSameEl " << nEvts_trksMatchSameEl << endl;
    cout<<" howMany evts found zero " << nEvts_zero << endl;
    cout<<" howMany evts found Hum " << nEvts_hum << endl;
    cout<<" howMany evts found Dois ou mais " << nEvts_dois << endl;
    
    
    ofstream myfile;
    myfile.open ("recoTypeCategories.txt");
    myfile << "leading\t" << nEvts_matchLeptons << "\t" << nEvts_matchJets << "\t" << nEvts_NoMatch << "\t" << nEvts_matchLeptonsAndJets <<  "\tnevts_leadAfterPrecut\t" << nevts_leadAfterPrecut << endl;
    myfile << "subleading\t" << nEvtsSub_matchLeptons << "\t" << nEvtsSub_matchJets << "\t" << nEvtsSub_NoMatch << "\t" << nEvtsSub_matchLeptonsAndJets<< "\tnevts_after_sublead\t" << nevts_SubleadAfterPrecut << endl;
    myfile.close();
    
	return EL::StatusCode::SUCCESS;
}

EL::StatusCode H2ZyAnalysis::initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  Info("initialize", "starting HgammaAnalysis::initialize");
  HgammaAnalysis::initialize();


  m_applySystematicLoop = config()->getBool("EventHandler.SystematicLoop");

  if(config()->isDefined("H2ZyAnalysis.ExpSyst"))   {
    m_expSyst=config()->getStrV("H2ZyAnalysis.ExpSyst");
  }

  // check if it is DATA or MC
  Info("initialize()", "Sample is %s", HG::isMC() ? "MC" : "DATA");
  V_llg.setConfig(*(config()), HG::isMC());

  // check whether we should apply the GRL or not
  if (config()->getBool("EventHandler.CheckGRL", false)) m_checkGRL = true;

  // setup trigger-related tools
  setupTrigger();

  // check whether we should do overlap removal (default is yes)
  if (config()->getBool("EventHandler.DoOverlapRemoval", true)) m_doOverlapRemoval = true;

  Info("initialize()", "Sample is %s", HG::isAOD() ? "AOD" : "DAOD");

  std::string upmenu = "/../";

  // initialize tool for MC equivalent luminosities
  if (HG::isMC()) {
    // initialize tool for Higgs xsection
    m_mclumi = new MCLumi((std::string)PathResolverFindCalibFile("H2Zy/MC16_lumi.txt"));
    m_mc_weight_xs = 0.0; // initialize xs_weight to 0 so that execute will compute it only the 1st time needed

    // saving sample type
    asg::AsgMetadataTool amdt("MetaDataTool");

    if (!amdt.inputMetaStore()->contains<xAOD::FileMetaData>("FileMetaData")) {
      m_mc_Year = -999;
    } else {
    const xAOD::FileMetaData *fmd = nullptr;
    amdt.inputMetaStore()->retrieve(fmd, "FileMetaData");
    std::string productionRelease, amiTag;
    fmd->value(xAOD::FileMetaData::productionRelease, productionRelease); // Release that was used to make the file
    fmd->value(xAOD::FileMetaData::amiTag, amiTag); // AMI tag used to process the file the last time
    TString sample_Type = HG::get_mcTypeUsingAmiTag(amiTag);

    m_mc_Year = -999;
    if( sample_Type.Contains("MC15") || sample_Type == "MC16a") m_mc_Year = 2016;
    else if (sample_Type == "MC16c" || sample_Type == "MC16d") m_mc_Year = 2017;
    else if (sample_Type == "MC16e") m_mc_Year = 2018;
    }
  }

  // create trees for info about overlap removal and add them to output
  m_ParticlesTreeBeforeOR = new HZG::OverlapRemovalTree("NOR");
  m_ParticlesTreeAfterOR  = new HZG::OverlapRemovalTree("OR");
  wk()->addOutput((TTree*)m_ParticlesTreeBeforeOR->GetTree());
  wk()->addOutput((TTree*)m_ParticlesTreeAfterOR->GetTree());

  // initialise tool for PDF reweighting
  m_pdfs.clear();
  if (HG::isMC()) {
    if (!(config()->getStr("H2ZyAnalysis.AddPDFWeights","NO")=="NO")) {
      const LHAPDF::PDFSet pdfset(config()->getStr("H2ZyAnalysis.AddPDFWeights").Data());
      const LHAPDF::PDFSet pdfset_mstw("MSTW2008nlo68cl");
      const LHAPDF::PDFSet pdfset_mmht("MMHT2014nlo68cl");
      //const LHAPDF::PDFSet pdfset_nn("NNPDF30_nlo_nf_5_pdfas");
      m_pdfs = pdfset.mkPDFs(); // pointers to PDF set members
      m_pdfs_mstw = pdfset_mstw.mkPDFs();
      m_pdfs_mmht = pdfset_mmht.mkPDFs();
      //m_pdfs_nn = pdfset_nn.mkPDFs();
    }
    m_theorySyst=config()->getBool("H2ZyAnalysis.AddTheorySystematicWeights");
  }

  m_truthselector.setMVA();

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode H2ZyAnalysis::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  //Info("execute", "starting HgammaAnalysis::execute");
  HgammaAnalysis::execute();
  //Info("execute", "end HgammaAnalysis::execute");

  // Retrieve event-level information
  SetEventInfo();
  //Info("execute", "end SetEventInfo");
  if (m_newFile){
	  addBookKeeping();
	  m_newFile=false;
  }

  // Define containers for storing good particles
  //xAOD::PhotonContainer S_photons(SG::VIEW_ELEMENTS);
  //xAOD::ElectronContainer S_electrons(SG::VIEW_ELEMENTS);
  //xAOD::MuonContainer S_muons(SG::VIEW_ELEMENTS);
  //xAOD::JetContainer S_jets(SG::VIEW_ELEMENTS);

  //-------- Getting different object containers -------------------
  // Loop over systematics for LLG selection
  bool is_goodevent = false ;
  for (auto sys: getSystematics()) {

	  TString sysname = sys.name();
	  if(!HG::isMC() && sysname!="") continue; 
	  if(!m_applySystematicLoop && sysname!="") continue;
	  //if(sysname.BeginsWith("MET_") || sysname.BeginsWith("JET_") || sysname.BeginsWith("FT_")) continue;
	  bool runSys=false;
	  for(int iSys=0;iSys<m_expSyst.size();++iSys)  {
		  if(sysname.BeginsWith(m_expSyst[iSys]))   {
			  runSys=true;
			  break;
		  }
	  }

	  if(m_applySystematicLoop&&sysname!=""&&!runSys)  continue;

	  //if(sysname.BeginsWith("MET_") || sysname.BeginsWith("FT_")) continue;
	  //if(!sysname.BeginsWith("JET_")) continue;
	  if(sysname!="") applySystematicVariation(sys);

	  m_cutFlow = cutflow(sysname);
	  //Info("execute", "end cutflow");
	  var::cutFlow.setValue(m_cutFlow);

	  if(HG::isMC()&&m_theorySyst) m_eventfill.SetTheoryUn(eventHandler()); 
	  //!{
	  //!        //xAOD::HiggsWeights hw=eventHandler()->higgsWeights();
	  //!        TruthWeightTools::HiggsWeights hw=eventHandler()->higgsWeights();

	  //!        double n = hw.nominal;

	  //!        double ratio=m_feventmap["mc_weight_final"]/n;

	  //!        m_feventmap["mc_weight_nominal"]=n*ratio;

	  //!        for(int iSyst=0;iSyst<hw.pdf4lhc_unc.size();++iSyst)  {
	  //!      	  m_feventmap[Form("mc_weight_pdf4lhc_unc_%d",iSyst)]=hw.pdf4lhc_unc[iSyst];
	  //!        }

	  //!        m_feventmap["mc_weight_alphaS_up"]=hw.alphaS_up*ratio;
	  //!        m_feventmap["mc_weight_alphaS_dn"]=hw.alphaS_dn*ratio;

	  //!        for(int iSyst=0;iSyst<hw.qcd.size();++iSyst)  {
	  //!      	  m_feventmap[Form("mc_weight_qcd_%d",iSyst)]=hw.qcd[iSyst]*ratio;
	  //!        }

	  //!        /*
	  //!           for(int iSyst=0;iSyst<hw.ggF_qcd_wg1().size();++iSyst)  {
	  //!           m_feventmap[Form("mc_weight_ggF_qcd_wg1_%d",iSyst)]=hw.ggF_qcd_wg1()[iSyst]*ratio;
	  //!           }
	  //!           */

	  //!        for(int iSyst=0;iSyst<hw.ggF_qcd_stxs.size();++iSyst)  {
	  //!      	  m_feventmap[Form("mc_weight_ggF_qcd_stxs_%d",iSyst)]=hw.ggF_qcd_stxs[iSyst]*ratio;
	  //!        }

	  //!        for(int iSyst=0;iSyst<hw.ggF_qcd_2017.size();++iSyst)  {
	  //!      	  m_feventmap[Form("mc_weight_ggF_qcd_2017_%d",iSyst)]=hw.ggF_qcd_2017[iSyst]*ratio;
	  //!        }
	  //!}

	  {
		  double w = 1.0;
		  for (int cut=ALLEVTS; cut<m_cutFlow; ++cut) {
			  if(HG::isMC()) w= eventHandler()->mcWeight()*eventHandler()->vertexWeight()*eventHandler()->pileupWeight();
			  //w = var::weightInitial();
			  if(cut==ALLEVTS_NOPU && HG::isMC() ) w= eventHandler()->mcWeight()*eventHandler()->vertexWeight();
			  if (cut > pre_sel && cut<=TRIGGER_MATCH) {
				  w *= m_feventmap["mc_weight_leptonscalefactor"];
				  w *= m_feventmap["mc_weight_triggerscalefactor"];
			  }
			  if (cut >TRIGGER_MATCH) w = m_feventmap["mc_weight_final"];
			  fillCutFlow(CutEnum(cut),sysname, w);
		  }
	  }

	  // vas::isPassed is used to decide whether to store the event in the output or not..
	  // we keep only events that pass all cuts except the very final ones (ph ID and isolation; 
	  // needed for purity studies; and mllg, used for sidebands in final fit)
	  var::isPassed.setValue(true);
	  if(m_cutFlow<=TRIGGER_MATCH) { var::isPassed.setValue(false); continue; }
	  //if(m_cutFlow<Ph_ptcut2+1) { var::isPassed.setValue(false); continue; }
	  // if we are here, there's at least one llg candidate

	  m_eventfill.retrieveVtxSumPt(electrons, muons);

	  //------------ llg systematic selection ----------------
	  if(!m_MxAODinput){
		  // Apply Truth matching ----
		  if (HG::isMC()) 	  m_truthselector.Reco_Truth_Matching( *cand_llg );
		  //----- save llg output information (tree) ------
		  //--------------
		  m_ieventmap["EventInfo.cutflow"]=m_cutFlow;

		  if(!m_btreemap.count(sysname)){
			  SetupOutputTree(V_llg, sysname);
			  m_btreemap[sysname]=1;
		  }
		  if(m_btreemap[sysname]){
			  FillVars(m_vbn, *cand_llg);
			  map_outputTree[sysname]->Fill();
		  }
	  }

	  if(sys.name()=="")
	  {
		  //for (auto ptcl : photons) S_photons.push_back(ptcl);
		  //for (auto ptcl : electrons) S_electrons.push_back(ptcl);
		  //for (auto ptcl : muons) S_muons.push_back(ptcl);
		  //for (auto ptcl : jets) S_jets.push_back(ptcl);
	  }

	  is_goodevent = true;	  
	  // example on how to store quantities in HGamEventInfo of MxAOD
	  // eventHandler()->storeVar("pT_yy1", 1.0);
	  // eventHandler()->storeVar("pT_yy2", 1.0);

  }//end of the systematics loop


  if(!is_goodevent) { store()->clear();  return EL::StatusCode::SUCCESS; }


  //---- write MxAOD
  for (auto sys: getSystematics()) {
    TString sysname = sys.name();
    if(!HG::isMC() && sys.name()!="") continue; 
    if(!m_applySystematicLoop && sys.name()!="") continue;
    //if(sysname.BeginsWith("MET_") || sysname.BeginsWith("JET_") || sysname.BeginsWith("FT_")) continue;
      bool runSys=false;
      for(int iSys=0;iSys<m_expSyst.size();++iSys)  {
        if(sysname.BeginsWith(m_expSyst[iSys]))   {
            runSys=true;
            break;
        }
      }

      if(m_applySystematicLoop&&sysname!=""&&!runSys)  continue;
    //if(sysname.BeginsWith("MET_") || sysname.BeginsWith("FT_")) continue;
    if(sysname!="") applySystematicVariation(sys);

    bool ispass = var::isPassed();
    if(!ispass )
     {
      // example on how to store quantities in HGamEventInfo of MxAOD
      //eventHandler()->storeVar("pT_yy1", -1.0);
     }
    //HG::VarHandler::getInstance()->write(); // do not call without calling fill()
  }

  // write object containers - do not call without calling fill() ! - this is why this is commented.
  //CP_CHECK("execute()", photonHandler  ()->writeContainer(S_photons , m_photonContainerName ));
  //CP_CHECK("execute()", electronHandler()->writeContainer(S_electrons, m_elecContainerName));
  //CP_CHECK("execute()", jetHandler     ()->writeContainer(S_jets, m_jetContainerName     ));
  //CP_CHECK("execute()", muonHandler    ()->writeContainer(S_muons, m_muonContainerName    ));

  // write event info container - do not call without calling fill()
  //eventHandler()->writeEventInfo();

  // do the real write
  //event()->fill();

  // cleanup temporary objects
  delete cand_llg->m_electron_viewcontainer;
  delete cand_llg->m_muon_viewcontainer;
  store()->clear();

  return EL::StatusCode::SUCCESS;

}

EL::StatusCode H2ZyAnalysis :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.

  m_histo_cutflow->SetFillStyle(3013);
  m_histo_cutflow->SetFillColor(13);
  ResetCounterName(m_event_cutflow_name,m_histo_cutflow);
  m_histo_cutflow_wt->SetFillStyle(3013);
  m_histo_cutflow_wt->SetFillColor(13);
  ResetCounterName(m_event_cutflow_name,m_histo_cutflow_wt);

  // Final process with trigger hist
  finalizeTrigger();
  triggerEfficiencyPerItem();


  for (std::map<Str,TH1F*>::iterator it_map = m_cutflowhistoTH1F.begin();it_map!=m_cutflowhistoTH1F.end();++it_map)
  {
	  (it_map->second)->SetFillStyle(3013);
	  (it_map->second)->SetFillColor(13);
	  ResetCounterName(m_event_cutflow_name,it_map->second);
  }

//hyp  if (m_summaryTree) {
//hyp    m_initialEvents = nEventsProcessed;
//hyp    m_finalEvents = nEventsDxAOD ;
//hyp    m_initialSumOfWeights = sumOfWeights;
//hyp    m_finalSumOfWeights = sumOfWeightsDxAOD;
//hyp    m_initialSumOfWeightsSquared = sumOfWeightsSquared;
//hyp    m_finalSumOfWeightsSquared = sumOfWeightsSquaredDxAOD;
//hyp    m_summaryTree->Fill();
//hyp  }

  return EL::StatusCode::SUCCESS;
}

void H2ZyAnalysis :: ResetVars ()
{   
  m_pass_grl= false;
  m_pass_pv = false;
  m_pass_quality = false;
  m_nMuons = 0;
  m_nElectrons = 0;
  m_nPhotons = 0;
}

void H2ZyAnalysis :: FillVars(vector<TString> varname, HZG::Object_llg &_llg)
{ 
  for(int i=0; i<(int)varname.size(); i++)
  {
    FillVars(varname[i], _llg); 
  }

}

void H2ZyAnalysis :: FillVars(TString varname, HZG::Object_llg &_llg)
{ 
  if((_llg.Map_float).find(varname) != (_llg.Map_float).end()){
    append_map(varname, _llg.Map_float, m_feventmap );
  }
  else if((_llg.Map_int).find(varname) != (_llg.Map_int).end()){
    append_map(varname, _llg.Map_int, m_ieventmap );
  }
  else if((_llg.Map_uint).find(varname) != (_llg.Map_uint).end()){
    append_map(varname, _llg.Map_uint, m_uieventmap );
  }
  else if((_llg.Map_bool).find(varname) != (_llg.Map_bool).end()){
    append_map(varname, _llg.Map_bool, m_beventmap );
  }
  else if((_llg.Map_ull).find(varname) != (_llg.Map_ull).end()){
    append_map(varname, _llg.Map_ull, m_ulleventmap );
  }
  else if((_llg.Map_IV).find(varname) != (_llg.Map_IV).end()){
    append_map(varname, _llg.Map_IV, m_iveventmap );
  }
  else if((_llg.Map_FV).find(varname) != (_llg.Map_FV).end()){
    append_map(varname, _llg.Map_FV, m_fveventmap );
  }
  else if((_llg.Map_TLVV).find(varname) != (_llg.Map_TLVV).end()){
    append_map(varname, _llg.Map_TLVV, m_tlvveventmap );
  }
  else {
    //std::cout<<"something not in llg object class"<<std::endl;
    // something not in llg object class
  }


}   


void H2ZyAnalysis :: SetupOutputTree (HZG::Container_llg &_V_llg, TString sysname)
{   
        // create HZG_Tree for systematic variation
	TString treename="HZG_Tree";
	if(sysname!="") treename = treename+"_"+sysname;
	m_outputTree = new TTree(treename,treename);
        map_outputTree[sysname] = m_outputTree;

	//
	std::set <TString> Set_variablename;  
	Set_variablename.insert(cand_llg->Set_variablename.begin(), cand_llg->Set_variablename.end());

	for (std::set <TString>::iterator si=Set_variablename.begin(); si!=Set_variablename.end(); si++)
		FillVars(*si, *cand_llg); 

	//m_feventmap
	for (std::map<Str,float>::iterator it_map = m_feventmap.begin();it_map!=m_feventmap.end();++it_map)
		MakeSingleBranch(it_map->first, m_feventmap, m_outputTree,m_vbn);

	//m_ieventmap
	for (std::map<Str,int>::iterator it_map = m_ieventmap.begin();it_map!=m_ieventmap.end();++it_map)
		MakeSingleBranch(it_map->first, m_ieventmap, m_outputTree,m_vbn);

	//m_uieventmap
	for (std::map<Str,unsigned int>::iterator it_map = m_uieventmap.begin();it_map!=m_uieventmap.end();++it_map)
		MakeSingleBranch(it_map->first, m_uieventmap, m_outputTree,m_vbn);

	//m_beventmap
	for (std::map<Str,bool>::iterator it_map = m_beventmap.begin();it_map!=m_beventmap.end();++it_map)
		MakeSingleBranch(it_map->first, m_beventmap, m_outputTree,m_vbn);

	//m_ulleventmap
	for (std::map<Str,unsigned long long>::iterator it_map = m_ulleventmap.begin();it_map!=m_ulleventmap.end();++it_map)
		MakeSingleBranch(it_map->first, m_ulleventmap, m_outputTree,m_vbn);

	//m_iveventmap
	for (std::map<Str,vector<int>>::iterator it_map = m_iveventmap.begin();it_map!=m_iveventmap.end();++it_map)
		MakeVectorBranch(it_map->first, m_iveventmap, m_outputTree,m_vbn);

	//m_fveventmap
	for (std::map<Str,vector<float>>::iterator it_map = m_fveventmap.begin();it_map!=m_fveventmap.end();++it_map)
		MakeVectorBranch(it_map->first, m_fveventmap, m_outputTree,m_vbn);

	//m_tlveventmap
	for (std::map<Str,vector<TLorentzVector>>::iterator it_map = m_tlvveventmap.begin();it_map!=m_tlvveventmap.end();++it_map)
		MakeVectorBranch(it_map->first, m_tlvveventmap, m_outputTree,m_vbn);


	wk()->addOutput (m_outputTree);
}

void H2ZyAnalysis :: SetupSummaryTree ()
{
  m_summaryTree = new TTree("summary","summary");
  //m_summaryTree->Branch("isMC",&m_isMC);
  //m_summaryTree->Branch("isDerivation",&m_isDerivation);
  m_summaryTree->Branch("initialSumW",&m_initialSumOfWeights);
  m_summaryTree->Branch("finalSumW",&m_finalSumOfWeights);
  m_summaryTree->Branch("initialSumW2",&m_initialSumOfWeightsSquared);
  m_summaryTree->Branch("finalSumW2",&m_finalSumOfWeightsSquared);
  m_summaryTree->Branch("initialEvents",&m_initialEvents);
  m_summaryTree->Branch("finalEvents",&m_finalEvents);
  // should add these 2 only for MC
  m_summaryTree->Branch("mc_channel_number",&m_mc_channel_number);
  m_summaryTree->Branch("mc_weight_xs",&m_mc_weight_xs);
  // should add RunNumber for data

  // Add summary info to the output file
  wk()->addOutput (m_summaryTree);


}

//-------
// save the event level information and MC truth information ------
//-- define the information for different weights, including MCweight, vertexweight, pileupweight
//-- define the generator information: Xsec, filter
//-- create the histograms to save cutflow information--
//------------------------------------------------------
void H2ZyAnalysis :: SetEventInfo()
{
	int mctype=-1;
	m_initialWeight=1.;
	float HiggsResMass = -100;
	m_hasDphoton=0;
        //------------------------- initialize mc truth --------------
	// search if there is a prompt photon in the event (for Z+gamma/Z+jet overlap removal)
    if(HG::isMC())  {
	    m_truthselector.initialize(*(config()),  event(), &m_beventmap, &m_ieventmap, &m_feventmap, &m_iveventmap, &m_fveventmap, &m_tlvveventmap);
	    m_truthselector.searchDphoton();
	    m_hasDphoton = m_truthselector.hasDphoton;
	    m_hasttyphoton = m_truthselector.hasttyphoton;
	    m_ieventmap["mc_hasPromptPhoton"]=m_hasDphoton;
	    m_ieventmap["mc_hasttgammaPhoton"]=m_hasttyphoton;
    }

	// saving event weights
	m_eventfill.initialize(*(config()), eventInfo(),  event(), &m_beventmap, &m_ieventmap, &m_uieventmap, &m_ulleventmap, &m_feventmap, &m_iveventmap, &m_fveventmap, &m_tlvveventmap);
	m_eventfill.seteventInfo(eventHandler());
	if(HG::isMC()) m_eventfill.SaveMCWeights(m_pdfs, m_pdfs_mstw, m_pdfs_mmht);

	if (HG::isMC()) {
		// retrieve mc_channel_number
		m_mc_channel_number = eventInfo()->mcChannelNumber();		

		// retrieve initial weight (generator*pile-up*vertex)
		//m_initialWeight=var::weightInitial();			
		//! m_initialWeight= eventHandler()->mcWeight()*eventHandler()->vertexWeight()*eventHandler()->pileupWeight();

		// retrieve MC generator weight(s)
		const std::vector< float > weights = eventInfo()->mcEventWeights();
		if( weights.size() > 0 )
			m_mc_weight_gen = weights[0];
		else
			m_mc_weight_gen = 1.0;

		// compute cross-section weight (only once)
		if(m_mcchannel_state.count((int)m_mc_channel_number) ==0)
		{
			m_mcchannel_state[(int)m_mc_channel_number]=1;
			float mc_genxsec  = m_mclumi->getXsec(m_mc_channel_number)*1.e6; // in fb
			float mc_genfeff  = m_mclumi->getFEff(m_mc_channel_number);
			m_mc_weight_xs = mc_genxsec*mc_genfeff; 
			m_mc_totevent = m_mclumi->getEvts(m_mc_channel_number);

			//long  mc_gennevts = m_mclumi->getEvts(D3PDdataset);
			//float mc_genlumi  = m_mclumi->getLumi(D3PDdataset)*1e-6; // in fb-1
		}


		// save truth info in ntuple
		//if(eventInfo()->eventNumber()==100) {m_truthselector.dumptruth(); getchar();}
		m_eventNumber = eventInfo()->eventNumber();
		if(!m_MxAODinput) m_truthselector.truthsave(m_eventNumber, m_mc_channel_number);

	}

	m_feventmap["mc_weight_xs"]=m_mc_weight_xs;
	m_feventmap["mc_weight_gen"]=m_mc_weight_gen;

	// about lumi nomorlization
	m_ieventmap["EventInfo.Year"] = m_mc_Year;
	if(!HG::isMC()) {
		if (eventInfo()->runNumber()<320000) m_ieventmap["EventInfo.Year"] = 2016;
		else if (eventInfo()->runNumber()<341000) m_ieventmap["EventInfo.Year"] = 2017;
		else m_ieventmap["EventInfo.Year"] = 2018;
	}
	m_feventmap["mc_weight_lumi"] = 1.0;
	if(m_ieventmap["EventInfo.Year"] == 2016) m_feventmap["mc_weight_lumi"]= 36.20766;
	else if(m_ieventmap["EventInfo.Year"] == 2017) m_feventmap["mc_weight_lumi"]= 44.3074;
	else if(m_ieventmap["EventInfo.Year"] == 2018) m_feventmap["mc_weight_lumi"]= 58.4501;

	//----- define the histogram to fill the weight information for each mc_channel_number ---
	if(HG::isMC()) {
		hist_cutflow_name = Form("cutflow_%d",(int) m_mc_channel_number);
	}
	else {
		hist_cutflow_name = Form("cutflow_%d",(int) eventInfo()->runNumber());
	}

	if(m_cutflowhistoTH1F.count(hist_cutflow_name)==0)
	{
		// the histogram doesn't exist - create it and add it to the output file!
		for (auto sys: getSystematics()) {

			TString sysname = sys.name();
			if(!HG::isMC() && sysname!="") continue;
			if(!m_applySystematicLoop && sysname!="") continue;
			//if(sysname.BeginsWith("MET_") || sysname.BeginsWith("JET_")|| sysname.BeginsWith("FT_")) continue;
			bool runSys=false;
			for(int iSys=0;iSys<m_expSyst.size();++iSys)  {
				if(sysname.BeginsWith(m_expSyst[iSys]))   {
					runSys=true;
					break;
				}
			}

			if(m_applySystematicLoop&&sysname!=""&&!runSys)  continue;
			//if(sysname.BeginsWith("MET_") || sysname.BeginsWith("FT_")) continue;
			//if(!sysname.BeginsWith("JET_")) continue;
			if (sysname!="") sysname = "_"+sysname;

			m_cutflowhistoTH1F[Form("%s%s", hist_cutflow_name.Data(), sysname.Data())] = new TH1F(Form("%s%s", hist_cutflow_name.Data(), sysname.Data()), Form("%s%s", hist_cutflow_name.Data(),sysname.Data()), 30, 0, 30.);
			wk()->addOutput(m_cutflowhistoTH1F[Form("%s%s", hist_cutflow_name.Data(),sysname.Data())]);
			m_cutflowhistoTH1F[Form("%s%s_w", hist_cutflow_name.Data(),sysname.Data())] = new TH1F(Form("%s%s_w", hist_cutflow_name.Data(),sysname.Data()), Form("%s%s_w", hist_cutflow_name.Data(),sysname.Data()), 30, 0, 30.);
			wk()->addOutput(m_cutflowhistoTH1F[Form("%s%s_w", hist_cutflow_name.Data(),sysname.Data())]);
			m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(),sysname.Data())] = new TH1F(Form("%s%s_w2", hist_cutflow_name.Data(),sysname.Data()), Form("%s%s_w2", hist_cutflow_name.Data(),sysname.Data()), 30, 0, 30.);
			wk()->addOutput(m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(),sysname.Data())]);
		}

	}

}

//-------------------------------------------------------------------
//  TRIGGER RELATED FUNCTIONS
//-------------------------------------------------------------------

void H2ZyAnalysis::createTriggerHists() {

  // Get list of triggers
  m_list_of_requiredTriggers = config()->getStrV("EventHandler.RequiredTriggers");
  const int nbins_trigger_items = m_list_of_requiredTriggers.size();
  const float nbins_low = -0.5;
  const float nbins_high = nbins_trigger_items-0.5;

  // create and add to output histogram of trigger pass (unweighted)
  m_histo_trigger_pass = new TH1F("histo_trigger_pass", "histo_trigger_pass", nbins_trigger_items, nbins_low, nbins_high);
  wk()->addOutput(m_histo_trigger_pass);

  // create and add to output histogram of trigger pass (weighted)
  m_histo_trigger_pass_wt = new TH1F("histo_trigger_pass_wt", "histo_trigger_pass_wt", nbins_trigger_items, nbins_low, nbins_high);
   wk()->addOutput(m_histo_trigger_pass_wt);

   // create and add to output histogram of trigger pass+match (unweighted)
   m_histo_trigger_pass_match = new TH1F("histo_trigger_pass_match", "histo_trigger_pass_match", nbins_trigger_items, nbins_low, nbins_high);
   wk()->addOutput(m_histo_trigger_pass_match);

   // create and add to output histogram of trigger pass+match (weighted)
   m_histo_trigger_pass_match_wt = new TH1F("histo_trigger_pass_match_wt", "histo_trigger_pass_match_wt", nbins_trigger_items, nbins_low, nbins_high);
   wk()->addOutput(m_histo_trigger_pass_match_wt);

   // histogram for trigger efficiency (denominator) (unweighted)
   m_histo_trigger_denom = new TH1F("histo_trigger_denom", "histo_trigger_denom", nbins_trigger_items, nbins_low, nbins_high);
   wk()->addOutput(m_histo_trigger_denom);

   // histogram for trigger efficiency (denominator) (weighted)
   m_histo_trigger_denom_wt = new TH1F("histo_trigger_denom_wt", "histo_trigger_denom_wt", nbins_trigger_items, nbins_low, nbins_high);
   wk()->addOutput(m_histo_trigger_denom_wt);

   m_histo_trigSF = new TH1F("histo_trigSF", "histo_trigSF", 50, -.05, 1.55);
   wk()->addOutput(m_histo_trigSF);

}

void H2ZyAnalysis::setupTrigger() {


  // Get list of triggers
  //m_list_of_requiredTriggers = config()->getStrV("EventHandler.RequiredTriggers");
  
  // Decide on applying trigger
  if (config()->getBool("EventHandler.CheckTriggers", false)) m_checkTrig = true;
  m_applyTrig = ( m_checkTrig && (m_list_of_requiredTriggers.size() > 0) );
  
  // Decide on applying trigger matching
  if (config()->getBool("EventHandler.CheckTriggerMatching", false)) m_checkTrigMatch = true;
  m_applyTrigMatch = ( m_checkTrigMatch && (m_list_of_requiredTriggers.size() > 0) );

  return;

}

bool H2ZyAnalysis::fillTriggerInfo(int channel, float weight) {

  int istep_trigger = 0;
  m_trigger_passed_items = 0;
  m_trigger_matched_items = 0;
  bool is_passed_matched = false;
  bool is_passed_trigger = false;
  if ( !m_MxAODinput) {
      //m_histo_trigger_denom->Fill(m_list_of_requiredTriggers.size());
    for (auto trig: m_list_of_requiredTriggers) {
      m_event_trigger_name[istep_trigger] = trig.Data();
      float bin_number =  istep_trigger;
      if ( (channel==2 && trig.Contains("mu")) || (channel==1 && !trig.Contains("mu")) || trig.Contains("HLT_g") ) {
         m_histo_trigger_denom->Fill(bin_number);
         m_histo_trigger_denom_wt->Fill(bin_number, weight);

         if ( eventHandler()->passTrigger(trig) ) {
             is_passed_trigger=true;
           m_trigger_passed_items |=  (1 << istep_trigger);
           m_histo_trigger_pass->Fill(bin_number);
           m_histo_trigger_pass_wt->Fill(bin_number, weight);

           if ( m_applyTrigMatch && passTriggerMatch(trig, cand_llg->m_photon_viewcontainer , cand_llg->m_electron_viewcontainer, cand_llg->m_muon_viewcontainer, NULL) ) {
              is_passed_matched = true;
              m_trigger_matched_items |= (1 << istep_trigger);
              m_histo_trigger_pass_match->Fill(bin_number);
              m_histo_trigger_pass_match_wt->Fill(bin_number, weight);
           }
         }
      }
      istep_trigger++ ;
     }
    //if(is_passed_trigger)   {
    //    m_histo_trigger_pass->Fill(m_list_of_requiredTriggers.size());
    //}
    //if(is_passed_matched)   {
    //    m_histo_trigger_pass_match->Fill(m_list_of_requiredTriggers.size());
    //}
  }

  return is_passed_matched;

}


void H2ZyAnalysis::triggerEfficiencyPerItem() {

  if ( ! m_applyTrig ) return;

  //Trigger efficiency
  TH1F*  trigeff = (TH1F*) m_histo_trigger_pass->Clone();
  trigeff->Divide( m_histo_trigger_pass, m_histo_trigger_denom, 1., 1., "cl=0.683 b(1,1) mode");
  trigeff->Scale(100.);
  trigeff->SetName("TriggerEfficienciesPassed");
  trigeff->GetYaxis()->SetTitle( "Trigger efficiency (%)" );
  wk()->addOutput( trigeff );

  //Trigger efficiency (passed + matched)
  TH1F*  trigeffmatched = (TH1F*) m_histo_trigger_pass_match->Clone();
  trigeffmatched->Divide( m_histo_trigger_pass_match, m_histo_trigger_denom, 1., 1., "cl=0.683 b(1,1) mode");
  trigeffmatched->Scale(100.);
  trigeffmatched->SetName("TriggerEfficienciesMatched");
  trigeffmatched->GetYaxis()->SetTitle( "Trigger efficiency (passed+matched) (%)" );
  wk()->addOutput( trigeffmatched );
 
  return;

}


void H2ZyAnalysis::finalizeTrigger() {

  ResetCounterName(m_event_trigger_name, m_histo_trigger_pass);
  ResetCounterName(m_event_trigger_name, m_histo_trigger_pass_wt);

  ResetCounterName(m_event_trigger_name, m_histo_trigger_pass_match);
  ResetCounterName(m_event_trigger_name, m_histo_trigger_pass_match_wt);

  ResetCounterName(m_event_trigger_name, m_histo_trigger_denom);
  ResetCounterName(m_event_trigger_name, m_histo_trigger_denom_wt);  

  return;

}

void H2ZyAnalysis::RecSave(HZG::Object_llg *&_llg){


	m_eventfill.RecSave(eventHandler(), _llg, photonHandler(), electronHandler(), muonHandler(), photons, electrons, muons);

        //------ trigger scale factor --------
	m_histo_trigSF->Fill(m_feventmap["mc_weight_triggerscalefactor"]);

}

void H2ZyAnalysis::fillCutFlow(CutEnum cut,  TString sysname, double w) 
{
	if(sysname==""){
		m_histo_cutflow->Fill(cut);
		m_histo_cutflow_wt->Fill(cut,w);
	}
	if (sysname=="") {
	  m_cutflowhistoTH1F[Form("%s", hist_cutflow_name.Data())]->Fill(cut);
	  m_cutflowhistoTH1F[Form("%s_w", hist_cutflow_name.Data())]->Fill(cut,w);
	  m_cutflowhistoTH1F[Form("%s_w2", hist_cutflow_name.Data())]->Fill(cut,w*w);
	}
	else {
	  m_cutflowhistoTH1F[Form("%s_%s", hist_cutflow_name.Data(), sysname.Data())]->Fill(cut);
	  m_cutflowhistoTH1F[Form("%s_%s_w", hist_cutflow_name.Data(), sysname.Data())]->Fill(cut,w);
	  m_cutflowhistoTH1F[Form("%s_%s_w2", hist_cutflow_name.Data(),sysname.Data())]->Fill(cut,w*w);
	}
}

H2ZyAnalysis::CutEnum H2ZyAnalysis::cutflow(TString sysname)
{

	//timeval t1;
	//gettimeofday(&t1, NULL);
	//std::cout<<"what's the time ? before cutflow = " << t1.tv_usec << std::endl; 

	m_cutFlow = initialcutflow();
	if(m_cutFlow!=Initial_sel) {return m_cutFlow;}
	//return m_cutFlow;

	all_correctedphotons = m_eventfill.all_photon_container(photonHandler());
	photons = m_eventfill.photon_container(photonHandler());

	all_correctedelectrons = m_eventfill.all_electron_container(electronHandler());
	electrons = m_eventfill.electron_container(electronHandler());

	muons0 = m_eventfill.muon_container(muonHandler());

	loose_jets = m_eventfill.loose_jet_container(jetHandler());

    //Reconstruction studies
    
    //Basic reconstructed electron pre-selection 
    xAOD::ElectronContainer m_preSelElectrons(SG::VIEW_ELEMENTS);
  	for (auto electron : all_correctedelectrons) 
    {
    		if (!electronHandler()->passOQCut(electron)) { continue; }
    		if (!electronHandler()->passPtEtaCuts(electron)) { continue; }
            //if (!electronHandler()->passIPCuts(electron)) { continue; } //will be used later, after y*y classification
    		if (!electronHandler()->passHVCut(electron)) { continue; }
    		m_preSelElectrons.push_back(electron);
  	}
    
    vector<const xAOD::TruthParticle*> m_new_truthLeptonsFromZ = m_truthselector.TruthLeptonsFromZ();
    
    
    //Reco type with leading and sub-leading lepton
    if(m_new_truthLeptonsFromZ.size() > 0)
    {

        //dR for two true leptons (no precut)
        double dEta = m_new_truthLeptonsFromZ.at(0)->eta() - m_new_truthLeptonsFromZ.at(1)->eta();
        double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
        dPhi = (dPhi <=M_PI)? dPhi : 2*M_PI-dPhi;
        double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
        m_dR_TrueLeptons->Fill(dR);
        
        //two true leptons kinematics (no precut)
        m_TrueLeptons_LeadingPt->Fill(m_new_truthLeptonsFromZ.at(0)->pt()/1000);
        m_TrueLeptons_SubLeadingPt->Fill(m_new_truthLeptonsFromZ.at(1)->pt()/1000);
        
        m_TrueLeptons_LeadingEta->Fill(m_new_truthLeptonsFromZ.at(0)->eta());
        m_TrueLeptons_SubLeadingEta->Fill(m_new_truthLeptonsFromZ.at(1)->eta());
        
        //Precut true leptons
        bool truthLeadPassPrecut = truePassPreCut(m_new_truthLeptonsFromZ.at(0));
        bool truthSubLeadPassPrecut = truePassPreCut(m_new_truthLeptonsFromZ.at(1));
        
        if(truthLeadPassPrecut && truthSubLeadPassPrecut)
        {
            m_dR_TrueLeptons_PreCuts->Fill(dR);
            
            //two true leptons kinematics (no precut)
            m_TrueLeptons_LeadingPt_PreCuts->Fill(m_new_truthLeptonsFromZ.at(0)->pt()/1000);
            m_TrueLeptons_SubLeadingPt_PreCuts->Fill(m_new_truthLeptonsFromZ.at(1)->pt()/1000);
        
            m_TrueLeptons_LeadingEta_PreCuts->Fill(m_new_truthLeptonsFromZ.at(0)->eta());
            m_TrueLeptons_SubLeadingEta_PreCuts->Fill(m_new_truthLeptonsFromZ.at(1)->eta());
        }
        
        
        double matchLep_lead = false;
        double matchLep_sublead = false;
        double matchLep_sublead_jets = false;
        double matchLep_sublead_nomatch = false;
        
        xAOD::ElectronContainer removedMatchToLeading (SG::VIEW_ELEMENTS);

        
        xAOD::Electron* recoHasMatchedLeading;
        bool hasRecoMatchedLead = false;
        bool isCat_LeadingGoesIntoLeptons =  false;
        
        double res_l1 = 0;
        double res_2d_l1 = 0;
        
        //leading lepton
        if( truthLeadPassPrecut )
        {
            cout<<"LEADING LEPTON pT " << m_new_truthLeptonsFromZ.at(0)->pt()/1000 << " pdgID " << m_new_truthLeptonsFromZ.at(0)->pdgId() << " barcode " << m_new_truthLeptonsFromZ.at(0)->barcode() << " status " << m_new_truthLeptonsFromZ.at(0)->status() << endl;
            nevts_leadAfterPrecut++;
            
            bool matchLeptons_dR = false;
            bool Nomatch = false;
            bool MatchedJet = false;
            
            
            cout<<" CHECK FOR TRUE LEADING " << endl;
            
            //==========================
            // CHECK MATCHING LEPTONS
            //==========================
            vector<xAOD::Electron*> elTruthMatched = match_DR_electrons(m_new_truthLeptonsFromZ.at(0), all_correctedelectrons);    
            
            if(elTruthMatched.size())
                matchLeptons_dR = true;
            
            //==========================
            // CHECK MATCHING JETS
            //==========================
            vector<xAOD::Jet*> jetTruthMatched = match_DR_jet(m_new_truthLeptonsFromZ.at(0), loose_jets);
            if(jetTruthMatched.size())
                MatchedJet = true;
        
            //if( elTruthMatched.size() ) removedMatchToLeading = skipRecoMatchedToTrueLeading(elTruthMatched.at(0), all_correctedelectrons);
            //else cout<<" truth el matching with true leading electron not found! " << endl;
           
            
            //Check with getTruthParticle function
            bool matchingFunction = trueMatch_toolFunction(m_new_truthLeptonsFromZ.at(0), all_correctedelectrons);
            if(matchingFunction) nEvts_trueLeadMatchARecoViaFunc++;
            
            //=======================================================
            // COUNTING NUMBER OF EVENTS IN EACH INDIVIDUAL RECO TYPE
            //=======================================================
              
            //True leading goes into leptons
            if(matchLeptons_dR)
            {
                
                matchLep_lead = true;
                isCat_LeadingGoesIntoLeptons = true;
                
                
                nEvts_matchLeptons++;
                
                //category fraction vs dR histogram
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryLeptons_leadingTruth->Fill(dR);
                
                
                //truth leading electron pT
                double trueLep_Pt = m_new_truthLeptonsFromZ.at(0)->pt()/1000;
                double trueLep_Eta = m_new_truthLeptonsFromZ.at(0)->eta();
                m_TrueLeadGoesIntoLeptons_TrueLeadingPt->Fill(trueLep_Pt);
                m_TrueLeadGoesIntoLeptons_TrueLeadingEta->Fill(trueLep_Eta);
                
                //reco Matched electron pT
                double recoMatched_Pt = elTruthMatched.at(0)->pt()/1000;
                double recoMatched_Eta = elTruthMatched.at(0)->eta();
                m_TrueLeadGoesIntoLeptons_RecoElMatchedPt->Fill(recoMatched_Pt);
                m_TrueLeadGoesIntoLeptons_RecoElMatchedEta->Fill(recoMatched_Eta);
                
                
                //leading electron resolution
                double resolution = recoMatched_Pt/trueLep_Pt;                
                
                
                //2d fraction
                TLorentzVector lep1, lep2;
                
                
                lep1.SetPtEtaPhiM(m_new_truthLeptonsFromZ.at(0)->pt(), m_new_truthLeptonsFromZ.at(0)->eta(), m_new_truthLeptonsFromZ.at(0)->phi(), m_new_truthLeptonsFromZ.at(0)->m());
                
                lep2.SetPtEtaPhiM(m_new_truthLeptonsFromZ.at(1)->pt(), m_new_truthLeptonsFromZ.at(1)->eta(), m_new_truthLeptonsFromZ.at(1)->phi(), m_new_truthLeptonsFromZ.at(1)->m());
                
                TLorentzVector sumTrueLeptons = lep1 + lep2;
                double res2d = recoMatched_Pt/(sumTrueLeptons.Pt()/1000);
               
                //cout<<" leading resolution " << resolution << " res2d  " << res2d  << endl;
            
                m_TrueLeadGoesIntoLeptons_1d_resolution->Fill(resolution);
                m_TrueLeadGoesIntoLeptons_2d_resolution->Fill(res2d, resolution);
                
                //cout<<" Sum 4-mom true leptons " << sumTrueLeptons.Pt()/1000 << " 1d res " <<  resolution << " 2d res " << res2d << endl;
                
                double dEta_l1 = m_new_truthLeptonsFromZ.at(0)->eta() -  elTruthMatched.at(0)->eta();
                double dPhi_l1 = m_new_truthLeptonsFromZ.at(0)->phi() - elTruthMatched.at(0)->phi();
                dPhi_l1 = (dPhi_l1<=M_PI)? dPhi_l1 : 2*M_PI-dPhi_l1;
                double dR_l1 = sqrt( pow(dEta_l1,2) + pow(dPhi_l1,2) );
                
                m_TrueLeadGoesIntoLeptons_dR_l1AndReco->Fill(dR_l1);
                
                //How often the truth matched electron (with truth leading) ALSO truth match the truth sub-leading
                n_evts_TrueSub_Out++;
                bool isPassPreCut = truePassPreCut(m_new_truthLeptonsFromZ.at(1));
                if(isPassPreCut)
                {
                    n_evts_TrueSub_In++;
                    double dEta_l2RecoMatch = m_new_truthLeptonsFromZ.at(1)->eta() -  elTruthMatched.at(0)->eta();
                    double dPhi_l2RecoMatch = m_new_truthLeptonsFromZ.at(1)->phi() -  elTruthMatched.at(0)->phi();
                    dPhi_l2RecoMatch = (dPhi_l2RecoMatch <= M_PI)? dPhi_l2RecoMatch : 2*M_PI-dPhi_l2RecoMatch;
                    double dR_l2RecoMatch = sqrt( pow(dEta_l2RecoMatch,2) + pow(dPhi_l2RecoMatch,2) );
                    
                    if( dR_l2RecoMatch < 0.1)
                    {
                        n_LeadIntoLeptons_l2RecoMatch++;
                        
                        m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d->Fill(resolution);
                        m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d->Fill(res2d, resolution);
                        m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR->Fill(dR_l2RecoMatch);
                        
                        if( res2d > 0.95 && res2d < 1.05) nEvts_bothTrueMatchSameRecoEl_Strip++;
                        
                    } 
                }
                
                //How much events fall into the strip between 0.95-1.05 (2d x-axis ratio)
                if( res2d > 0.95 && res2d < 1.05)
                {
                    n_evts_strip++;
                    if(isPassPreCut) nEvts_strip_trueSubPassPreCut++;
                    
                    /*
                    double dEta_strip = m_new_truthLeptonsFromZ.at(1)->eta() -  elTruthMatched.at(0)->eta();
                    double dPhi_l2RecoMatch = m_new_truthLeptonsFromZ.at(1)->phi() -  elTruthMatched.at(0)->phi();
                    dPhi_l2RecoMatch = (dPhi_l2RecoMatch <= M_PI)? dPhi_l2RecoMatch : 2*M_PI-dPhi_l2RecoMatch;
                    double dR_l2RecoMatch = sqrt( pow(dEta_l2RecoMatch,2) + pow(dPhi_l2RecoMatch,2) );
                    
                    if( dR_l2RecoMatch < 0.1)
                    {
                        n_LeadIntoLeptons_Strip_l2RecoMatch++;
                        
                        hLead_matchLeptons_1d_resolution_Strip->Fill(resolution);
                        hLead_matchLeptons_2d_resolution_Strip->Fill(res2d, resolution);
                        
                    } 
                    */
                    
                }
                
                
                if(resolution > 0.99 && resolution < 1.01)  hLead_matchLeptons_1d_InsideresolutionWindow->Fill(resolution);
                else
                {
                    hLead_matchLeptons_outsideResWindow_Pt->Fill(trueLep_Pt);  
                    hLead_matchLeptons_1d_resolution_outsideResWindow->Fill(resolution);
                    hLead_matchLeptons_2d_resolution_outsideResWindow->Fill(res2d, resolution);
                    hLead_matchLeptons_RecoEl_outsideResWindow_Cluster->Fill(elTruthMatched.at(0)->caloCluster()->e()/1000);
                } 
                
                if(resolution > 1.01)
                {
                    hLead_matchLeptons_1d_resolution_Above1p01->Fill(resolution);
                    hLead_matchLeptons_2d_resolution_Above1p01->Fill(res2d, resolution);
                }
                if(resolution < 0.99 )
                {
                    hLead_matchLeptons_1d_resolution_Below0p99->Fill(resolution);
                    hLead_matchLeptons_2d_resolution_Below0p99->Fill(res2d, resolution);
                }
                
                
                
                
                
            } 
            
            //True leading goes into leptons only
            if(matchLeptons_dR && !MatchedJet )
                nEvts_IsmatchedLeptonOnly++;
            
            //True leading goes into jets
            if(!matchLeptons_dR && MatchedJet)
            {
                nEvts_matchJets++;
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryJets_leadingTruth->Fill(dR);
            } 
            
            //True leading goes nothing (noMatch)
            if(!matchLeptons_dR && !MatchedJet)
            {
                nEvts_NoMatch++;
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryNoMatch_leadingTruth->Fill(dR);
                
            } 
            
            //True leading goes leptons + jets
            if(matchLeptons_dR && MatchedJet)
            {
                nEvts_matchLeptonsAndJets++;
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryLeptonsAndJets_leadingTruth->Fill(dR);
                
            } 
            
            
            
        }
        
        //subleading lepton
        if( truthSubLeadPassPrecut )
        {
            cout<<"SUB-LEADING LEPTON pT " << m_new_truthLeptonsFromZ.at(1)->pt()/1000 << " pdgID " << m_new_truthLeptonsFromZ.at(1)->pdgId() << " barcode " << m_new_truthLeptonsFromZ.at(1)->barcode() << " status " << m_new_truthLeptonsFromZ.at(1)->status() << endl;
            
            nevts_SubleadAfterPrecut++;
            
            bool matchLeptons_dR = false;
            bool Nomatch = false;
            bool MatchedJet = false;
            
            cout<<" CHECK FOR TRUE SUB-LEADING " << endl;
            
            //==========================
            // CHECK MATCHING LEPTONS
            //==========================
            //vector<xAOD::Electron*> elTruthMatched = match_DR_electrons(m_new_truthLeptonsFromZ.at(1), removedMatchToLeading);
            
            vector<xAOD::Electron*> elTruthMatched = match_DR_electrons(m_new_truthLeptonsFromZ.at(1), all_correctedelectrons);
            
            if(elTruthMatched.size())
            {
                //cout<<" truth matched el pT " << elTruthMatched.at(0)->pt()/1000 << endl;
                matchLeptons_dR = true;
            } 
            
            //==========================
            // CHECK MATCHING JETS
            //==========================
            vector<xAOD::Jet*> jetTruthMatched = match_DR_jet(m_new_truthLeptonsFromZ.at(1), loose_jets);
            if(jetTruthMatched.size())
            {
                //cout<<" truth matched jet pT " << jetTruthMatched.at(0)->pt()/1000 << endl;   
                MatchedJet = true;
            } 
        
            //Check with getTruthParticle function
            bool matchingFunction = trueMatch_toolFunction(m_new_truthLeptonsFromZ.at(1), all_correctedelectrons);
            if(matchingFunction) nEvts_trueSubLeadMatchARecoViaFunc++;
            
            //check broken links
            bool brokenLink = trueMatch_BrokenLink(m_new_truthLeptonsFromZ.at(1), all_correctedelectrons);
            if(!brokenLink) nEvts_trueSubLead_BrokenLink++;
            
            //check who are those non-matched true particles
            vector<const xAOD::TruthParticle*> nonRealLeptons = trueMatch_noMatched(m_new_truthLeptonsFromZ.at(1), all_correctedelectrons);
            
            cout<<" non-real leptons linked size " << nonRealLeptons.size() << endl;
            
            //Tracking truth match
            xAOD::ElectronContainer recoTrk(SG::VIEW_ELEMENTS);
            if(!matchingFunction)
            {
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() - m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi <=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                    
                m_dEta_trueLeptons->Fill(dEta);
                m_dPhi_trueLeptons->Fill(dPhi);
                m_dR_realLeptons->Fill(dR);
                m_trueSubLead_pT_noMatched->Fill(m_new_truthLeptonsFromZ.at(1)->pt()/1000);
                 
                cout<<"=================="<<endl;
                cout<<" GSF CONTAINER " << endl;
                cout<<"=================="<<endl;
                
                const xAOD::TrackParticleContainer* tracks = 0;
                if(!event()->retrieve(tracks, "GSFTrackParticles" ).isSuccess())
                {
                    Error("execute()", "Failed to retrieve GSFTrackParticles container." );
                }
                
                cout<<" GSFTrackParticles size " << tracks->size() << endl;
                for(const auto track : *tracks)
                {
                    
                    float truth_prob = track->auxdata<float>("truthMatchProbability");
                    cout<<" get prob " << truth_prob << endl;
                    /*
                    xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*trk);
                    const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*tp_p);
                    if(!truthlepton) continue;
                    
                    if(truthlepton)
                    {
                        if(abs(truthlepton->pdgId()) == abs(m_new_truthLeptonsFromZ.at(1)->pdgId()) && truthlepton->barcode() == m_new_truthLeptonsFromZ.at(1)->barcode() && truthlepton->status() == m_new_truthLeptonsFromZ.at(1)->status())
                        {
                            cout<<" Found track, missing electron " << endl; 
                            m_GSFtrk_pT->Fill(trk->pt()/1000);
                            m_GSFtrk_eta->Fill(trk->eta());
                            nEvts_trueSubLead_FoundTrack++;
                            
                        }
                          
                    }
                    */
                }
                
                
                //xAOD::TrackParticle* matchedTracksToSubLeading = foundIDTrack(m_new_truthLeptonsFromZ.at(1));
                
                /*if( matchedTracksToSubLeading.size() == 0)
                {
                    
                    nEvts_zero++;
                    for(auto trk : *IDtracks)
                    {
                        
                        xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*trk);
                        const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*tp_p);
                        
                        if(!truthlepton) continue;
                        
                        const xAOD::TruthParticle* truthPart_p = (xAOD::TruthParticle*)getTheOriginalPointer(*truthlepton);
                        
                        cout<<"DOES IT EXIST NOW?"<< abs(truthPart_p->pdgId()) << " " << truthPart_p->barcode() << endl;
                        
                        if(truthPart_p->barcode() == m_new_truthLeptonsFromZ.at(1)->barcode())
                            cout<<" LET'S HAVE A LOOK! " << endl;
                        
                    
                        if(abs(truthlepton->pdgId()) == abs(m_new_truthLeptonsFromZ.at(1)->pdgId()) && truthlepton->barcode() == m_new_truthLeptonsFromZ.at(1)->barcode() && truthlepton->status() == m_new_truthLeptonsFromZ.at(1)->status())
                        {
                                cout<<" Found ID track, missing electron " << endl;
                                nEvts_trueSubLead_FoundTrack_IDTrk++;
                                m_IDtrk_pT->Fill(trk->pt()/1000);
                                m_IDtrk_eta->Fill(trk->eta());
                                moreThanOne++;
                        }
                        
                    }
                }*/
                //if( matchedTracksToSubLeading ) nEvts_hum++;
                
                
                //bool tracksMatchingSameEl = findTracksMatchingSameElectron(m_new_truthLeptonsFromZ.at(0), m_new_truthLeptonsFromZ.at(1), all_correctedelectrons);
                
                //if(tracksMatchingSameEl) nEvts_trksMatchSameEl++;
                
            }
            else cout<<" ALREADY MATCHED! " << endl;
            
            
            

            //=======================================================
            // COUNTING NUMBER OF EVENTS IN EACH INDIVIDUAL RECO TYPE
            //=======================================================
            
            
            if(matchLeptons_dR)
            {
                matchLep_sublead = true;
                nEvtsSub_matchLeptons++;
                
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                
                categoryLeptons_SubleadingTruth->Fill(dR);
                
                
                
                
                //truth subleading electron pT
                double trueLep_Pt = m_new_truthLeptonsFromZ.at(1)->pt()/1000;
                double trueLep_Eta = m_new_truthLeptonsFromZ.at(1)->eta();
                m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt->Fill(trueLep_Pt);
                m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta->Fill(trueLep_Eta);
                
                //reco Matched electron pT
                double recoMatched_Pt = elTruthMatched.at(0)->pt()/1000;
                double recoMatched_Eta = elTruthMatched.at(0)->eta();
                m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt->Fill(recoMatched_Pt);
                m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta->Fill(recoMatched_Eta);
            
                
                //subleading electron resolution
                double resolution = recoMatched_Pt/trueLep_Pt;                
                
                //2d fraction
                TLorentzVector lep1, lep2;
                
                lep1.SetPtEtaPhiM(m_new_truthLeptonsFromZ.at(0)->pt(), m_new_truthLeptonsFromZ.at(0)->eta(), m_new_truthLeptonsFromZ.at(0)->phi(), m_new_truthLeptonsFromZ.at(0)->m());
                
                lep2.SetPtEtaPhiM(m_new_truthLeptonsFromZ.at(1)->pt(), m_new_truthLeptonsFromZ.at(1)->eta(), m_new_truthLeptonsFromZ.at(1)->phi(), m_new_truthLeptonsFromZ.at(1)->m());
                
                TLorentzVector sumTrueLeptons = lep1 + lep2;
                
                double res2d = recoMatched_Pt/(sumTrueLeptons.Pt()/1000);
                
            
                m_TrueSubLeadGoesIntoLeptons_1d_resolution->Fill(resolution);
                m_TrueSubLeadGoesIntoLeptons_2d_resolution->Fill(res2d, resolution);
            
                /*
                if(resolution > 0.99 && resolution < 1.01) hSubLead_matchLeptons_1d_InsideresolutionWindow->Fill(resolution);
                else
                {
                    hSubLead_matchLeptons_outsideResWindow_Pt->Fill(trueLep_Pt);
                    hSubLead_matchLeptons_1d_resolution_outsideResWindow->Fill(resolution);
                    hSubLead_matchLeptons_2d_resolution_outsideResWindow->Fill(res2d, resolution);
                    hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster->Fill(elTruthMatched.at(0)->caloCluster()->e()/1000);
                }
                
                if(resolution > 1.01)
                {
                    hSubLead_matchLeptons_1d_resolution_Above1p01->Fill(resolution);
                    hSubLead_matchLeptons_2d_resolution_Above1p01->Fill(res2d, resolution);
                }
                if(resolution < 0.99 )
                {
                    hSubLead_matchLeptons_1d_resolution_Below0p99->Fill(resolution);
                    hSubLead_matchLeptons_2d_resolution_Below0p99->Fill(res2d, resolution);
                }
                */
                
                
            } 
            if(!matchLeptons_dR && MatchedJet)
            {
                matchLep_sublead_jets = true;
                nEvtsSub_matchJets++;
                
                
            } 
            if(!matchLeptons_dR && !MatchedJet)
            {
                
                cout<<" NOT MATCH SUB-LEAD EVTS " << endl;
                matchLep_sublead_nomatch = true;
                nEvtsSub_NoMatch++;
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryNoMatch_SubleadingTruth->Fill(dR);
                
                /*
                m_allTracks = trackHandler()->getCorrectedContainer();
                //m_allTracks.sort(HG::TrackHandler::comparePt);
                
                
                 const xAOD::IParticle* trk_l1;
                 const xAOD::IParticle* trk_l2;
                  
                 bool foundTrkLead = false;
                 bool foundTrkSubLead = false;
                
                
                for(auto trk : m_allTracks)
                {
                    const xAOD::TruthParticle* true_ptr = xAOD::TruthHelpers::getTruthParticle(*trk);
                    if(true_ptr)
                    {
                         if( (abs(true_ptr->pdgId()) == abs(m_new_truthLeptonsFromZ.at(0)->pdgId())) && (true_ptr->barcode() == m_new_truthLeptonsFromZ.at(0)->barcode()) && (true_ptr->status() == m_new_truthLeptonsFromZ.at(0)->status())  )
                        {
                        
                            trkFoundLeading_NoMatch->Fill(trk->pt()/1000);
                            trk_l1 = trk;
                            cout<<" found trk for leading " << trk->pt()/1000 << endl;
                            foundTrkLead = true;
                        }
                    
                        if( (abs(true_ptr->pdgId()) == abs(m_new_truthLeptonsFromZ.at(1)->pdgId())) && (true_ptr->barcode() == m_new_truthLeptonsFromZ.at(1)->barcode()) && (true_ptr->status() == m_new_truthLeptonsFromZ.at(1)->status())  )
                        {
                        
                            trkFoundSubLeading_NoMatch->Fill(trk->pt()/1000);
                            trk_l2 = trk;
                            cout<<" found trk for subleading " << trk->pt()/1000 << endl;
                            foundTrkSubLead = true;
                        }
                    }
                    
                }
                
                double trkRes_leading = (trk_l1->pt()/1000)/(m_new_truthLeptonsFromZ.at(0)->pt()/1000);
                double trkRes_leadingAndSubleading =  (trk_l1->pt()/1000)/((m_new_truthLeptonsFromZ.at(0)->pt()/1000) + (m_new_truthLeptonsFromZ.at(1)->pt()/1000));
                
                h_1dResolution_NoMatch->Fill(trkRes_leading);
                h_2dResolution_NoMatch->Fill(trkRes_leadingAndSubleading, trkRes_leading);
                
                h_sumOfTrue_NoMatch->Fill((m_new_truthLeptonsFromZ.at(0)->pt()/1000) + (m_new_truthLeptonsFromZ.at(1)->pt()/1000) );
                
                
                
                
                cout<<" leading Pt" << m_new_truthLeptonsFromZ.at(0)->pt()/1000 << " sub leading Pt " << m_new_truthLeptonsFromZ.at(1)->pt()/1000 << " sum truth leptons Pt " << (m_new_truthLeptonsFromZ.at(0)->pt()/1000) + (m_new_truthLeptonsFromZ.at(1)->pt()/1000) << endl;
                
                //dR between two found tracks inside the jet
                if(foundTrkLead && foundTrkSubLead)
                {
                    double dEta_trk = trk_l1->eta() - trk_l2->eta();
                    double dPhi_trk = trk_l1->phi() - trk_l2->phi();
                    dPhi_trk = (dPhi_trk<=M_PI)? dPhi_trk : 2*M_PI-dPhi_trk;
                    double dR_trk = sqrt( pow(dEta_trk,2) + pow(dPhi_trk,2) );
                    dR_diTrks->Fill(dR_trk);   
                    
                }
                else if(foundTrkLead && !foundTrkSubLead)
                {
                    cout<<" trk associated to leading only found! " << endl;
                    
                } 
                else if (!foundTrkLead && foundTrkSubLead)
                {
                    cout<<" trk associated to sub-leading only found! " << endl;
                } 
                else
                {
                    cout<<" no trks found! " << endl;   
                } 
                */
                
                
            } 
            if(matchLeptons_dR && MatchedJet)
            {
                nEvtsSub_matchLeptonsAndJets++;
                
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() -  m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                categoryLeptonsAndJets_SubleadingTruth->Fill(dR);
                
            }
            
            /*
            if(!matchLeptons_dR && MatchedJet)
            {
                
                cout<<" SUBLEP GOES INTO JET " << endl;
                
                double pT_res = ( ((jetMatchedMinDR->pt()/1000) - (m_new_truthLeptonsFromZ.at(1)->pt()/1000))/(m_new_truthLeptonsFromZ.at(1)->pt()/1000) );
                
                //cout<< " matched jet pT " << matched_jet->pt()/1000 << " truth sub-leading pT " << m_new_truthLeptonsFromZ.at(1)->pt()/1000 << " pT res " <<  pT_res << endl;
                
                jetMass->Fill(jetMatchedMinDR->m()/1000);
                jetPt->Fill(jetMatchedMinDR->pt()/1000);

                m_subLead_JetsOnlyEvts_pTresolution->Fill(pT_res);
                m_subLead_JetsOnlyEvts_truthSubleadpT->Fill(m_new_truthLeptonsFromZ.at(1)->pt()/1000);
                m_subLead_JetsOnlyEvts_matchedJetpT->Fill(jetMatchedMinDR->pt()/1000);

                TLorentzVector le1, le2;
                
                le1.SetPtEtaPhiM(m_new_truthLeptonsFromZ.at(1)->pt(), m_new_truthLeptonsFromZ.at(1)->eta(), m_new_truthLeptonsFromZ.at(1)->phi(), m_new_truthLeptonsFromZ.at(1)->m());
                le2.SetPtEtaPhiM(jetMatchedMinDR->pt(), jetMatchedMinDR->eta(), jetMatchedMinDR->phi(), jetMatchedMinDR->m());
                
                double n_ml1l2 = (le1 + le2).M();
                m_subLead_JetsOnlyEvts_Mljet->Fill(n_ml1l2/1000);
                
                std::vector<const xAOD::IParticle*> jettracks;
                jetMatchedMinDR->getAssociatedObjects<xAOD::IParticle>(xAOD::JetAttribute::GhostTrack,jettracks);
                //cout<<" jettracks " << jettracks.size() << endl;
                if( jettracks.size() < 5 ) nevts_nTrkJet++;
                
                 const xAOD::IParticle*  jetTrk_l1;
                 const xAOD::IParticle*  jetTrk_l2;
                
                 bool jettrkLead = false;
                 bool jettrkSubLead = false;
                
                
                double trk_l1_max = 0;
                double trk_l2_max = 0;
                
                for(auto trk : jettracks)
                {
                    const xAOD::TruthParticle* true_trk = xAOD::TruthHelpers::getTruthParticle(*trk);
                    
                    
                    if(true_trk)
                    {
                        if( (abs(true_trk->pdgId()) == abs(m_new_truthLeptonsFromZ.at(0)->pdgId())) && (true_trk->barcode() == m_new_truthLeptonsFromZ.at(0)->barcode()) && (true_trk->status() == m_new_truthLeptonsFromZ.at(0)->status())  )
                        {
                            if( trk->pt()/1000 > trk_l1_max )
                            {
                                cout<<" track Inside Jet Match true leading electron " << trk->pt()/1000 << endl;
                                trk_l1_max = trk->pt()/1000;
                                jetTrk_l1 = trk;
                                jettrkLead = true;   
                            }
                            
                        }
                            
                        if( (abs(true_trk->pdgId()) == abs(m_new_truthLeptonsFromZ.at(1)->pdgId())) && (true_trk->barcode() == m_new_truthLeptonsFromZ.at(1)->barcode()) && (true_trk->status() == m_new_truthLeptonsFromZ.at(1)->status())  )
                        {
                            if(trk->pt()/1000 > trk_l2_max)
                            {
                                cout<<" track Inside Jet Match true sub-leading electron " << trk->pt()/1000 << endl;
                                jetTrk_l2 = trk;
                                jettrkSubLead = true;
                            }
                            
                            
                        }
                            
                    }
                    
                }
                
                
                
                
                if(jettrkLead && jettrkSubLead)
                {
                    //dR between two found tracks inside the jet
                    double dEta_trk = jetTrk_l1->eta() - jetTrk_l2->eta();
                    double dPhi_trk = jetTrk_l1->phi() - jetTrk_l2->phi();
                    dPhi_trk = (dPhi_trk<=M_PI)? dPhi_trk : 2*M_PI-dPhi_trk;
                    double dR_trk = sqrt( pow(dEta_trk,2) + pow(dPhi_trk,2) );
                    dR_jetTracks_bothTracks->Fill(dR_trk);  
                    
                    //Mtrktrk
                    TLorentzVector trk1, trk2;

                    trk1.SetPtEtaPhiM(jetTrk_l1->pt(), jetTrk_l1->eta(), jetTrk_l1->phi(), jetTrk_l1->m());

                    trk2.SetPtEtaPhiM(jetTrk_l2->pt(), jetTrk_l2->eta(), jetTrk_l2->phi(), jetTrk_l2->m());

                    double m_mtrk1trk2 = (trk1 + trk2).M();

                    mtrk1trk2->Fill(m_mtrk1trk2/1000);
                    
                    cout<<" both associated " << jetTrk_l1->pt()/1000 << "\t" << jetTrk_l2->pt()/1000 << endl;
                    
                    jetTrackLeading_GoesIntoJet->Fill(jetTrk_l1->pt()/1000);
                    jetTrackSubLeading_GoesIntoJet->Fill(jetTrk_l2->pt()/1000);
                    
                }
                else if(jettrkLead && !jettrkSubLead)
                {
                    cout<<" only trk associated to leading found with mass: " << jetTrk_l1->m()/1000 << endl;
                    jetTrackLeading_GoesIntoJet->Fill(jetTrk_l1->pt()/1000);
                }
                else if(!jettrkLead && jettrkSubLead)
                {
                    cout<<" only trk associated to subleading found with mass: " << jetTrk_l2->m()/1000 << endl;
                    jetTrackSubLeading_GoesIntoJet->Fill(jetTrk_l2->pt()/1000);
                }
                else
                {
                    cout<<" no trk inside Jet match truth leptons! " << endl;
                }
                
                
                
                //dR between true leading and sub-leading 
                double dEta = m_new_truthLeptonsFromZ.at(0)->eta() - m_new_truthLeptonsFromZ.at(1)->eta();
                double dPhi = m_new_truthLeptonsFromZ.at(0)->phi() - m_new_truthLeptonsFromZ.at(1)->phi();
                dPhi = (dPhi<=M_PI)? dPhi : 2*M_PI-dPhi;
                double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                
                dRTwoTrueLeptons_catGoesIntoJet->Fill(dR);
                //pT leading and subleading
                Ptleading_catGoesIntoJet->Fill(m_new_truthLeptonsFromZ.at(0)->pt()/1000);
                PtSubleading_catGoesIntoJet->Fill(m_new_truthLeptonsFromZ.at(1)->pt()/1000);
                
                //dR leading with matched jet
                
                double dEta_l1Jet = m_new_truthLeptonsFromZ.at(0)->eta() - jetMatchedMinDR->eta();
                double dPhi_l1Jet = m_new_truthLeptonsFromZ.at(0)->phi() - jetMatchedMinDR->phi();
                dPhi_l1Jet = (dPhi_l1Jet <= M_PI)? dPhi_l1Jet : 2*M_PI-dPhi_l1Jet;
                double dR_l1Jet = sqrt( pow(dEta_l1Jet,2) + pow(dPhi_l1Jet,2) );
                
                dRl1Jet_catGoesIntoJet->Fill(dR_l1Jet);
                
                //dR subleading with matched jet
                
                double dEta_l2Jet = m_new_truthLeptonsFromZ.at(1)->eta() - jetMatchedMinDR->eta();
                double dPhi_l2Jet = m_new_truthLeptonsFromZ.at(1)->phi() - jetMatchedMinDR->phi();
                dPhi_l2Jet = (dPhi_l2Jet <= M_PI)? dPhi_l2Jet : 2*M_PI-dPhi_l2Jet;
                double dR_l2Jet = sqrt( pow(dEta_l2Jet,2) + pow(dPhi_l2Jet,2) );
                
                dRl2Jet_catGoesIntoJet->Fill(dR_l2Jet);
                
                
                
            }
            */
        }
        
      
       
        
    }
    
    //ALL TRUE LEPTONS - GetTruthParticle Function
    xAOD::ElectronContainer el_truthMatched(SG::VIEW_ELEMENTS);
    el_truthMatched = getRecoMatched( m_new_truthLeptonsFromZ, all_correctedelectrons );
    m_el_truthMatched->Fill(el_truthMatched.size());
    
    
	if (m_doOverlapRemoval)
	{
		if(config()->getBool("EventHandler.FillORTree",false)&&sysname=="") m_ParticlesTreeBeforeOR->Fill(electrons, muons0, photons, loose_jets);
		ZGam_Overlapremoval(overlapHandler(), *(config()), electrons, muons0, photons, loose_jets);
		if(config()->getBool("EventHandler.FillORTree",false)&&sysname=="") m_ParticlesTreeAfterOR->Fill(electrons, muons0, photons, loose_jets);
	}

	jets = m_eventfill.jet_container(loose_jets);

	muons = m_eventfill.muon_clean_container(muonHandler(),muons0);

    xAOD::JetContainer allPFlowJets = jetHandlerPFlow()->getCorrectedContainer();

    //xAOD::TauJetContainer m_allTaus = tauHandler()->getCorrectedContainer();
    //xAOD::TauJetContainer m_selTaus = tauHandler()->applySelection(m_allTaus);
    xAOD::TauJetContainer m_selTaus(SG::VIEW_ELEMENTS);

	met = m_eventfill.met_container(etmissHandler(), jetHandler(), photons, electrons, muons,&allPFlowJets,&m_selTaus);
	//met = m_eventfill.met_container(etmissHandler(), jetHandler(), photons, electrons, muons,&allPFlowJets);
	//met = m_eventfill.met_container(etmissHandler(), jetHandler(), photons, electrons, muons);

	//m_nPhotons = photons.size();
	m_nElectrons = electrons.size();
	m_nMuons = muons.size();

	// jet event clean	
	//m_eventfill.JetCleaningSave(jetHandler());

	// Truth container
	//m_eventfill.TruthMETSave(truthHandler(), sysname);

	//----------------------------------------------------------------------------------------------------
	//---- cut flow ------
	if(!m_MxAODinput){	
		double pre_mll=-999;
		V_llg.initialize(photons, electrons, muons, jets, loose_jets, met, all_correctedphotons, all_correctedelectrons);
		pre_mll=V_llg.get_maxmll();

		// overlap removal with the target di-lepton
		if(V_llg.m_container_llg.size() >= 1) {
			cand_llg = &((V_llg.m_container_llg)[0]);
			//std::cout<<"index = " << cand_llg->m_channel << " " << cand_llg->index_lepton1<< " " <<cand_llg->index_lepton2<<std::endl;
			//std::cout<<"size = " << (cand_llg->m_electron_viewcontainer)->size() << " " << (cand_llg->m_muon_viewcontainer)->size() << std::endl;

			xAOD::ElectronContainer S_electrons(SG::VIEW_ELEMENTS);
			xAOD::MuonContainer S_muons(SG::VIEW_ELEMENTS);

			if(cand_llg->m_channel==1) {
				S_electrons.push_back(electrons[cand_llg->index_lepton1]);
				S_electrons.push_back(electrons[cand_llg->index_lepton2]);
			}
			else if(cand_llg->m_channel==2){
				S_muons.push_back(muons[cand_llg->index_lepton1]);
				S_muons.push_back(muons[cand_llg->index_lepton2]);
			}

			ZGam_Overlapremoval2(overlapHandler(), *(config()), S_electrons, S_muons, photons, loose_jets);

			if(config()->getBool("H2ZyAnalysis.findMoreLepton",false))    {
				ZGam_Overlapremoval2(overlapHandler(), *(config()), electrons, muons, photons, loose_jets);
			}

			m_nPhotons = photons.size();

			jets = m_eventfill.jet_container(loose_jets);

            muonHandler()->decorateDeltaRJet(S_muons, jets);
			//met2 = m_eventfill.met_container(etmissHandler(), jetHandler(), photons, electrons, muons);

			m_eventfill.JetCleaningSave(jetHandler());

			m_eventfill.TruthMETSave(truthHandler(), sysname);

			//ZGam_Overlapremoval2(overlapHandler(), *(config()), *(cand_llg->m_electron_viewcontainer) , *(cand_llg->m_muon_viewcontainer), photons, loose_jets);

			// have to save the kinematics again after overlap removal
			cand_llg->initialize_g0();
			cand_llg->Store_llg();
		}

		//=======
		m_cutFlow = precutflow();
		if(m_cutFlow!=pre_sel) return m_cutFlow;
		m_cutFlow = finalcutflow();
	}
	return m_cutFlow;
}

//----------------------------------------
//--basic event selection criteria   -----
//----------------------------------------
H2ZyAnalysis::CutEnum H2ZyAnalysis::initialcutflow()
{

	first=true;
	if (first) {
		m_MxAODinput = false;
		//m_MxAODinput = event()->contains<xAOD::PhotonContainer>(m_photonContainerName.Data());
		first = false;
	}

	m_cutFlow = ALLEVTS;
	m_cutFlow = ALLEVTS_NOPU;
	istep_event=0;

        //------------------------------------ apply the separation criteria at the truth level ---------------
	if (config()->getBool("H2ZyAnalysis.TruthSample",false)) return ALLEVTS;
	else if (config()->getBool("H2ZyAnalysis.electrononly",false)) {if(m_ieventmap["mc_Z_decay_topo"]!=1111)  return ALLEVTS;}
	else if (config()->getBool("H2ZyAnalysis.muononly",false)) {if(m_ieventmap["mc_Z_decay_topo"]!=1313)  return ALLEVTS;}

	if(config()->getBool("H2ZyAnalysis.truthMllCut",false)) {
		double minTruMll=config()->getNum("H2ZyAnalysis.truthMinMll",999.);
		double maxTruMll=config()->getNum("H2ZyAnalysis.truthMaxMll",-999.);
		if(minTruMll>998.)  {
			if(!(m_feventmap["Z_truth_mass"]<=maxTruMll))    return ALLEVTS;
		}
		else if(maxTruMll<-998.)    {
			if(!(m_feventmap["Z_truth_mass"]>minTruMll))    	return ALLEVTS;

		}
		else    {
			if(!(m_feventmap["Z_truth_mass"]>minTruMll&&m_feventmap["Z_truth_mass"]<=maxTruMll))        return ALLEVTS;
		}
	}
	//----- ALL ----

	//------ GRL -----
	if (m_checkGRL)
	{
		if (HG::isMC())
			m_pass_grl = true;
		else
			m_pass_grl = eventHandler()->passGRL(eventInfo());

		if(!m_pass_grl) {  return GRL;} 
	}

	//----  PV checking ----
	if(!m_MxAODinput) m_pass_pv = m_eventfill.retrievePV();
	else m_pass_pv = true;
	if(!m_pass_pv)  { return PV; }

	//----  EVENT_QUALITY checking ----
	if(HG::isMC() || m_MxAODinput) m_pass_quality=true;
	else m_pass_quality=eventHandler()->passTile(eventInfo()) && eventHandler()->passLAr(eventInfo()) && eventHandler()->passCore(eventInfo());
	if(!m_pass_quality)  { return EVENT_QUALITY; }
	//---------------------------------

	//----  TRIGGER checking ----
	// apply trigger selection : passTriggers checks if at least one required triggers is satisfied
	m_trigger_passed = true;
	if ( m_applyTrig && !m_MxAODinput) {
		m_trigger_passed = (eventHandler()->passTriggers());
		if(!m_trigger_passed) { return Triggers;} 
	}
	m_beventmap["trigger_passed"]  = m_trigger_passed;

	return Initial_sel;
}

H2ZyAnalysis::CutEnum H2ZyAnalysis::precutflow()
{
	if(!(m_nMuons>=2 || m_nElectrons>=2)) return twolepton;
	if(V_llg.get_maxmll()<45) return mll_threshold;
	if (config()->getBool("H2ZyAnalysis.TwoLeptonSelection",false)) {if (!(m_nPhotons>=0 && (m_nMuons>=2 || m_nElectrons>=2))) return twolepton_onephoton;}
	else if (!(m_nPhotons>=1 && (m_nMuons>=2 || m_nElectrons>=2))) return twolepton_onephoton;
	return pre_sel;
}

H2ZyAnalysis::CutEnum H2ZyAnalysis::finalcutflow()
{
	// From here on we need at least one llg candidate
	if (V_llg.m_container_llg.size() == 0 ) {
		return llgcut;
	}

	//cand_llg = &((V_llg.m_container_llg)[0]);

	//std::cout<<"higss cand." << eventInfo()->eventNumber()<<" "<<m_feventmap["mc_weight_leptonscalefactor"]<< " " <<                                                 m_feventmap["mc_weight_triggerscalefactor"]<<" " <<  cand_llg->m_channel  <<std::endl;
	m_ieventmap["N_ll"] = (V_llg.m_container_llg).size();

	// Fill the trigger info
	bool trig_passed_matched = fillTriggerInfo(cand_llg->m_channel, m_initialWeight);
	m_uieventmap["trigger_passed_items"]  = m_trigger_passed_items;
	m_ieventmap["trigger_matched_items"] = m_trigger_matched_items;

	// Apply trigger matching
	if ( m_applyTrigMatch && ! trig_passed_matched ) return TRIGGER_MATCH;

	// Compute various quantities of the final state; they will eventually
	// be recorded if the event passes the desired cutflow step
	m_ieventmap["n_electron"]=m_nElectrons;
	m_ieventmap["n_muon"]=m_nMuons;

	RecSave(cand_llg);

	m_beventmap["llg_passallcuts"] = false;

	// m(ll) selection
	float minMll = config()->getNum("H2ZyAnalysis.FinalMinMllCut",-1.);
	if (minMll>0) {
		if (cand_llg->Map_float["ll_m_Zmassconstraint"]<minMll) return mllcut;
	}
	float maxMll = config()->getNum("H2ZyAnalysis.FinalMaxMllCut",-1.);
	if (maxMll>0) {
		if (cand_llg->Map_float["ll_m_Zmassconstraint"]>maxMll) return mllcut;
	}

	// loose photon pT cut
	double ph_minpt = config()->getNum("H2ZyAnalysis.LoosePhotonPtCut",-1.);
	if (ph_minpt>0.) {
		if (ph_minpt<1.) {
			if (cand_llg->Map_float["ph_pt"]/cand_llg->Map_float["llg_m_Zmassconstraint"]<ph_minpt) return Ph_ptcut1;
		}
		else {
			if (cand_llg->Map_float["ph_pt"]<ph_minpt) return Ph_ptcut1;
		}
	}

	// tight photon ID
	if( !m_beventmap["ph_istight"]) return Ph_ID;

	// photon isolation
	if(!m_beventmap[Form("ph_passiso_%s",config()->getStr("H2ZyAnalysis.FinalPhotonIsolationCut").Data())]) return Ph_ISO;

	// m(llg) cut
	float minMllg = config()->getNum("H2ZyAnalysis.FinalMinMllgCut",-1.);
	if (minMllg>0) {
		if (cand_llg->Map_float["llg_m_Zmassconstraint"]<minMllg) return MASSCUT;
	}
	float maxMllg = config()->getNum("H2ZyAnalysis.FinalMaxMllgCut",-1.);
	if (maxMllg>0) {
		if (cand_llg->Map_float["llg_m_Zmassconstraint"]>maxMllg) return MASSCUT;
	}

	// final photon pT cut
	ph_minpt = config()->getNum("H2ZyAnalysis.FinalPhotonPtCut");
	if (ph_minpt>0.) {
		if (ph_minpt<1.) {
			if (cand_llg->Map_float["ph_pt"]/cand_llg->Map_float["llg_m_Zmassconstraint"]<ph_minpt) return Ph_ptcut2;
		}
		else {
			if (cand_llg->Map_float["ph_pt"]<ph_minpt) return Ph_ptcut2;
		}
	}

	m_beventmap["llg_passallcuts"] = true;
	return PASSALL;
}

void H2ZyAnalysis::addBookKeeping()
{
	for (auto sys: getSystematics()) {

		TString sysname = sys.name();
		if(!HG::isMC() && sysname!="") continue;
		//if(sysname.BeginsWith("MET_") || sysname.BeginsWith("JET_")|| sysname.BeginsWith("FT_")) continue;
		bool runSys=false;
		for(int iSys=0;iSys<m_expSyst.size();++iSys)  {
			if(sysname.BeginsWith(m_expSyst[iSys]))   {
				runSys=true;
				break;
			}
		}

		if(m_applySystematicLoop&&sysname!=""&&!runSys)  continue;
		//if(sysname.BeginsWith("MET_") || sysname.BeginsWith("FT_")) continue;
		//if(!sysname.BeginsWith("JET_")) continue;
		if(!m_applySystematicLoop && sysname!="") continue;

		if (sysname!="") sysname = "_"+sysname;

		m_cutflowhistoTH1F[Form("%s%s", hist_cutflow_name.Data(), sysname.Data())]->AddBinContent(1, nEventsProcessed);
		m_cutflowhistoTH1F[Form("%s%s", hist_cutflow_name.Data(), sysname.Data())]->AddBinContent(2, nEventsDxAOD);

		m_cutflowhistoTH1F[Form("%s%s_w", hist_cutflow_name.Data(), sysname.Data())]->AddBinContent(1, sumOfWeights);
		m_cutflowhistoTH1F[Form("%s%s_w", hist_cutflow_name.Data(), sysname.Data())]->AddBinContent(2, sumOfWeightsDxAOD);

		m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(), sysname.Data())]->SetBinContent(1, sqrt( pow(m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(), sysname.Data())]->GetBinContent(1),2) + pow(sumOfWeightsSquared,2) ));
		m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(), sysname.Data())]->SetBinContent(2, sqrt( pow(m_cutflowhistoTH1F[Form("%s%s_w2", hist_cutflow_name.Data(), sysname.Data())]->GetBinContent(2),2) + pow(sumOfWeightsSquaredDxAOD,2) ));
	}
}

void H2ZyAnalysis::setcutflowname()
{
	m_event_cutflow_name[xAOD]="xAOD";
	m_event_cutflow_name[DxAOD]="DxAOD";
	m_event_cutflow_name[ALLEVTS]="ALL";
	m_event_cutflow_name[ALLEVTS_NOPU]="ALL_NOPU";
	m_event_cutflow_name[GRL]="GRL";
	m_event_cutflow_name[PV]="PV";
	m_event_cutflow_name[EVENT_QUALITY] = "EVENT_QUALITY";
	m_event_cutflow_name[Triggers]="TRIGGER";
	m_event_cutflow_name[Initial_sel]="initial_cut";
	m_event_cutflow_name[twolepton]="2l"; 
	m_event_cutflow_name[mll_threshold]="mll_loose_cut";
	m_event_cutflow_name[twolepton_onephoton]="2l+#gamma";
	m_event_cutflow_name[pre_sel]="pre_cut";
	m_event_cutflow_name[llgcut]="llg_size>0";
	m_event_cutflow_name[TRIGGER_MATCH]="TRIGGER_MATCH";
	m_event_cutflow_name[mllcut]="mll_final_cut";
	m_event_cutflow_name[Ph_ptcut1]="p_{T}(#gamma) first cut";
	m_event_cutflow_name[Ph_ID]="ID(#gamma)";
	m_event_cutflow_name[Ph_ISO]="Isolation(#gamma)";
	//m_event_cutflow_name[mllcut]=Form("mll>%lf",config()->getNum("H2ZyAnalysis.FinalMllCut",81.18));
	m_event_cutflow_name[MASSCUT]="m(ll#gamma)";
	m_event_cutflow_name[Ph_ptcut2]="p_{T}(#gamma) second cut";
}

//-------------------------------------------------------------------
//  Initialize histos for reconstruction studies
//-------------------------------------------------------------------
void H2ZyAnalysis::createRecoStudiesHists()
{
    //True leptons (no precut)
  m_dR_TrueLeptons = new TH1F("m_dR_TrueLeptons", "m_dR_TrueLeptons", 100, 0, 1.);
  wk()->addOutput(m_dR_TrueLeptons);
    
  m_TrueLeptons_LeadingPt = new TH1F("m_TrueLeptons_LeadingPt", "m_TrueLeptons_LeadingPt", 50, 0., 5000.);
  wk()->addOutput(m_TrueLeptons_LeadingPt);  
    
  m_TrueLeptons_SubLeadingPt = new TH1F("m_TrueLeptons_SubLeadingPt", "m_TrueLeptons_SubLeadingPt", 50, 0., 3000.);
  wk()->addOutput(m_TrueLeptons_SubLeadingPt);  

  m_TrueLeptons_LeadingEta = new TH1F("m_TrueLeptons_LeadingEta", "m_TrueLeptons_LeadingEta", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueLeptons_LeadingEta); 
    
  m_TrueLeptons_SubLeadingEta = new TH1F("m_TrueLeptons_SubLeadingEta", "m_TrueLeptons_SubLeadingEta", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueLeptons_SubLeadingEta); 
    
  //True leptons (with precut)    
  m_dR_TrueLeptons_PreCuts = new TH1F("m_dR_TrueLeptons_PreCuts", "m_dR_TrueLeptons_PreCuts", 100, 0, 1.);
  wk()->addOutput(m_dR_TrueLeptons_PreCuts);
    
  m_TrueLeptons_LeadingPt_PreCuts = new TH1F("m_TrueLeptons_LeadingPt_PreCuts", "m_TrueLeptons_LeadingPt_PreCuts", 50, 0., 5000.);
  wk()->addOutput(m_TrueLeptons_LeadingPt_PreCuts);  
    
  m_TrueLeptons_SubLeadingPt_PreCuts = new TH1F("m_TrueLeptons_SubLeadingPt_PreCuts", "m_TrueLeptons_SubLeadingPt_PreCuts", 50, 0., 3000.);
  wk()->addOutput(m_TrueLeptons_SubLeadingPt_PreCuts);  

  m_TrueLeptons_LeadingEta_PreCuts = new TH1F("m_TrueLeptons_LeadingEta_PreCuts", "m_TrueLeptons_LeadingEta_PreCuts", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueLeptons_LeadingEta_PreCuts); 
    
  m_TrueLeptons_SubLeadingEta_PreCuts = new TH1F("m_TrueLeptons_SubLeadingEta_PreCuts", "m_TrueLeptons_SubLeadingEta_PreCuts", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueLeptons_SubLeadingEta_PreCuts); 

    
  //Fractions vs true leptons dR
categoryLeptons_leadingTruth = new TH1F("categoryLeptons_leadingTruth", "categoryLeptons_leadingTruth", 50, 0., 1.);
wk()->addOutput(categoryLeptons_leadingTruth);
    
categoryJets_leadingTruth = new TH1F("categoryJets_leadingTruth", "categoryJets_leadingTruth", 50, 0., 1.);
wk()->addOutput(categoryJets_leadingTruth);
    
categoryNoMatch_leadingTruth = new TH1F("categoryNoMatch_leadingTruth", "categoryNoMatch_leadingTruth", 50, 0., 1.);
wk()->addOutput(categoryNoMatch_leadingTruth);
    
categoryLeptonsAndJets_leadingTruth = new TH1F("categoryLeptonsAndJets_leadingTruth", "categoryLeptonsAndJets_leadingTruth", 50, 0., 1.);
wk()->addOutput(categoryLeptonsAndJets_leadingTruth);
    
    
categoryLeptons_SubleadingTruth = new TH1F("categoryLeptons_SubleadingTruth", "categoryLeptons_SubleadingTruth", 50, 0., 1.);
wk()->addOutput(categoryLeptons_SubleadingTruth);
    
categoryJets_SubleadingTruth = new TH1F("categoryJets_SubleadingTruth", "categoryJets_SubleadingTruth", 50, 0., 1.);
wk()->addOutput(categoryJets_SubleadingTruth);
    
    
categoryNoMatch_SubleadingTruth = new TH1F("categoryNoMatch_SubleadingTruth", "categoryNoMatch_SubleadingTruth", 50, 0., 1.);
wk()->addOutput(categoryNoMatch_SubleadingTruth);
    
categoryLeptonsAndJets_SubleadingTruth = new TH1F("categoryLeptonsAndJets_SubleadingTruth", "categoryLeptonsAndJets_SubleadingTruth", 50, 0., 1.);
wk()->addOutput(categoryLeptonsAndJets_SubleadingTruth);
    
//Kinematics of true leptons and reconstructed matched electron (truth leading goes into leptons)
m_TrueLeadGoesIntoLeptons_TrueLeadingPt = new TH1F("m_TrueLeadGoesIntoLeptons_TrueLeadingPt", "m_TrueLeadGoesIntoLeptons_TrueLeadingPt", 50, 0., 5000.);
  wk()->addOutput(m_TrueLeadGoesIntoLeptons_TrueLeadingPt);  

m_TrueLeadGoesIntoLeptons_TrueLeadingEta = new TH1F("m_TrueLeadGoesIntoLeptons_TrueLeadingEta", "m_TrueLeadGoesIntoLeptons_TrueLeadingEta", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueLeadGoesIntoLeptons_TrueLeadingEta); 
    
m_TrueLeadGoesIntoLeptons_RecoElMatchedPt = new TH1F("m_TrueLeadGoesIntoLeptons_RecoElMatchedPt", "m_TrueLeadGoesIntoLeptons_RecoElMatchedPt", 50, 0., 5000.);
  wk()->addOutput(m_TrueLeadGoesIntoLeptons_RecoElMatchedPt); 
    
m_TrueLeadGoesIntoLeptons_RecoElMatchedEta = new TH1F("m_TrueLeadGoesIntoLeptons_RecoElMatchedEta", "m_TrueLeadGoesIntoLeptons_RecoElMatchedEta", 25, -2.5, 2.5);
wk()->addOutput(m_TrueLeadGoesIntoLeptons_RecoElMatchedEta); 
    
//Kinematics of true leptons and reconstructed matched electron (truth subleading goes into leptons)
m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt = new TH1F("m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt", "m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt", 50, 0., 3000.);
  wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingPt);  

m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta = new TH1F("m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta", "m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta", 25, -2.5, 2.5);
  wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_TrueSubLeadingEta); 
    
m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt = new TH1F("m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt", "m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt", 50, 0., 5000.);
  wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_RecoElMatchedPt); 
    
m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta = new TH1F("m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta", "m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta", 25, -2.5, 2.5);
wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_RecoElMatchedEta); 

//resolutions (truth leading goes into leptons)
m_TrueLeadGoesIntoLeptons_1d_resolution = new TH1F("m_TrueLeadGoesIntoLeptons_1d_resolution", "m_TrueLeadGoesIntoLeptons_1d_resolution", 500, 0., 2.);
wk()->addOutput(m_TrueLeadGoesIntoLeptons_1d_resolution);
    
m_TrueLeadGoesIntoLeptons_2d_resolution = new TH2F("m_TrueLeadGoesIntoLeptons_2d_resolution", "m_TrueLeadGoesIntoLeptons_2d_resolution", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(m_TrueLeadGoesIntoLeptons_2d_resolution); 
    
//resolutions (truth subleading goes into leptons)
m_TrueSubLeadGoesIntoLeptons_1d_resolution = new TH1F("m_TrueSubLeadGoesIntoLeptons_1d_resolution", "m_TrueSubLeadGoesIntoLeptons_1d_resolution", 500, 0., 2.);
wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_1d_resolution);
    
m_TrueSubLeadGoesIntoLeptons_2d_resolution = new TH2F("m_TrueSubLeadGoesIntoLeptons_2d_resolution", "m_TrueSubLeadGoesIntoLeptons_2d_resolution", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(m_TrueSubLeadGoesIntoLeptons_2d_resolution); 
    

m_TrueLeadGoesIntoLeptons_dR_l1AndReco = new TH1F("m_TrueLeadGoesIntoLeptons_dR_l1AndReco", "m_TrueLeadGoesIntoLeptons_dR_l1AndReco", 100, 0, 1.);
  wk()->addOutput(m_TrueLeadGoesIntoLeptons_dR_l1AndReco);
    
m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d = new TH1F("m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d", "m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d", 500, 0., 2.);
wk()->addOutput(m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_1d);
    
m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d = new TH2F("m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d", "m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_2d);
 
m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR =  new TH1F("m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR", "m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR", 100, 0, 1.);
  wk()->addOutput(m_TrueLeadGoesIntoLeptons_recoTruthMatchedByTrueSub_dR);
    

    
//new ones    
hSubLead_matchViaFunction_recoMatched_Pt = new TH1F("hSubLead_matchViaFunction_recoMatched_Pt", "hSubLead_matchViaFunction_recoMatched_Pt", 100, 0., 5000.);
wk()->addOutput(hSubLead_matchViaFunction_recoMatched_Pt);

hSubLead_matchLeptons_trueSubLead_Pt = new TH1F("hSubLead_matchLeptons_trueSubLead_Pt", "hSubLead_matchLeptons_trueSubLead_Pt", 100, 0., 5000.);
wk()->addOutput(hSubLead_matchLeptons_trueSubLead_Pt);
    
hSubLead_matchLeptons_recoMatched_Pt = new TH1F("hSubLead_matchLeptons_recoMatched_Pt", "hSubLead_matchLeptons_recoMatched_Pt", 100, 0., 5000.);
wk()->addOutput(hSubLead_matchLeptons_recoMatched_Pt);
    
hSubLead_matchLeptons_1d_resolution = new TH1F("hSubLead_matchLeptons_1d_resolution", "hSubLead_matchLeptons_1d_resolution", 500, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_resolution);
    
hSubLead_matchLeptons_2d_resolution = new TH2F("hSubLead_matchLeptons_2d_resolution", "hSubLead_matchLeptons_2d_resolution", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_2d_resolution); 
    
//more
hLead_matchLeptons_1d_InsideresolutionWindow = new TH1F("hLead_matchLeptons_1d_InsideresolutionWindow", "hLead_matchLeptons_1d_InsideresolutionWindow", 500, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_InsideresolutionWindow);
    
hSubLead_matchLeptons_1d_InsideresolutionWindow = new TH1F("hSubLead_matchLeptons_1d_InsideresolutionWindow", "hSubLead_matchLeptons_1d_InsideresolutionWindow", 500, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_InsideresolutionWindow);
    
hLeptons_matchLeptons_1d_BothInsideresolutionWindow = new TH1F("hLeptons_matchLeptons_1d_BothInsideresolutionWindow", "hLeptons_matchLeptons_1d_BothInsideresolutionWindow", 500, 0., 2.);
wk()->addOutput(hLeptons_matchLeptons_1d_BothInsideresolutionWindow);
    
hLeptons_matchLeptons_Mll_BothInsideresolutionWindow = new TH1F("hLeptons_matchLeptons_Mll_BothInsideresolutionWindow", "hLeptons_matchLeptons_Mll_BothInsideresolutionWindow", 100, 0., 200.);
wk()->addOutput(hLeptons_matchLeptons_Mll_BothInsideresolutionWindow);
  
//events outside the resolution window 0.99-1.01 (LEAD)
hLead_matchLeptons_outsideResWindow_Pt = new TH1F("hLead_matchLeptons_outsideResWindow_Pt", "hLead_matchLeptons_outsideResWindow_Pt", 100, 0., 5000.);
wk()->addOutput(hLead_matchLeptons_outsideResWindow_Pt); //Pt
    
hLead_matchLeptons_RecoEl_outsideResWindow_Cluster = new TH1F("hLead_matchLeptons_RecoEl_outsideResWindow_Cluster", "hLead_matchLeptons_RecoEl_outsideResWindow_Cluster", 100, 0., 5000.);
wk()->addOutput(hLead_matchLeptons_RecoEl_outsideResWindow_Cluster); //Pt
    
hLead_matchLeptons_1d_resolution_outsideResWindow = new TH1F("hLead_matchLeptons_1d_resolution_outsideResWindow", "hLead_matchLeptons_1d_resolution_outsideResWindow", 500, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_resolution_outsideResWindow); //1d resolution
    
hLead_matchLeptons_2d_resolution_outsideResWindow = new TH2F("hLead_matchLeptons_2d_resolution_outsideResWindow", "hLead_matchLeptons_2d_resolution_outsideResWindow", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hLead_matchLeptons_2d_resolution_outsideResWindow);  //2d resolution
    

    
//events outside the resolution window 0.99-1.01 (SUBLEAD)
hSubLead_matchLeptons_outsideResWindow_Pt = new TH1F("hSubLead_matchLeptons_outsideResWindow_Pt", "hSubLead_matchLeptons_outsideResWindow_Pt", 100, 0., 5000.);
wk()->addOutput(hSubLead_matchLeptons_outsideResWindow_Pt); //Pt
    
hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster = new TH1F("hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster", "hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster", 100, 0., 5000.);
wk()->addOutput(hSubLead_matchLeptons_RecoEl_outsideResWindow_Cluster); //Pt
    
hSubLead_matchLeptons_1d_resolution_outsideResWindow = new TH1F("hSubLead_matchLeptons_1d_resolution_outsideResWindow", "hSubLead_matchLeptons_1d_resolution_outsideResWindow", 500, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_resolution_outsideResWindow); //1d resolution
    
hSubLead_matchLeptons_2d_resolution_outsideResWindow = new TH2F("hSubLead_matchLeptons_2d_resolution_outsideResWindow", "hSubLead_matchLeptons_2d_resolution_outsideResWindow", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_2d_resolution_outsideResWindow);  //2d resolution

//cluster energy dR (leading and sub-leading matched)
dRCluster_twoRecoElMatched_NoResContraint = new TH1F("dRCluster_twoRecoElMatched_NoResContraint", "dRCluster_twoRecoElMatched_NoResContraint", 100, 0., 1.);
wk()->addOutput(dRCluster_twoRecoElMatched_NoResContraint);
    
dRCluster_twoRecoElMatched_BothInsideresolutionWindow = new TH1F("dRCluster_twoRecoElMatched_BothInsideresolutionWindow", "dRCluster_twoRecoElMatched_BothInsideresolutionWindow", 100, 0., 1.);
wk()->addOutput(dRCluster_twoRecoElMatched_BothInsideresolutionWindow);
    
//resolution leading and sub-leading > 1.01
hLead_matchLeptons_1d_resolution_Above1p01 = new TH1F("hLead_matchLeptons_1d_resolution_Above1p01", "hLead_matchLeptons_1d_resolution_Above1p01", 500, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_resolution_Above1p01); //1d resolution
    
hLead_matchLeptons_2d_resolution_Above1p01 = new TH2F("hLead_matchLeptons_2d_resolution_Above1p01", "hLead_matchLeptons_2d_resolution_Above1p01", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hLead_matchLeptons_2d_resolution_Above1p01);  //2d resolution
    
hSubLead_matchLeptons_1d_resolution_Above1p01 = new TH1F("hSubLead_matchLeptons_1d_resolution_Above1p01", "hSubLead_matchLeptons_1d_resolution_Above1p01", 500, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_resolution_Above1p01); //1d resolution
    
hSubLead_matchLeptons_2d_resolution_Above1p01 = new TH2F("hSubLead_matchLeptons_2d_resolution_Above1p01", "hSubLead_matchLeptons_2d_resolution_Above1p01", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_2d_resolution_Above1p01);  //2d resolution
    
//resolution leading and sub-leading < 0.99
hLead_matchLeptons_1d_resolution_Below0p99 = new TH1F("hLead_matchLeptons_1d_resolution_Below0p99", "hLead_matchLeptons_1d_resolution_Below0p99", 500, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_resolution_Below0p99); //1d resolution
    
hLead_matchLeptons_2d_resolution_Below0p99 = new TH2F("hLead_matchLeptons_2d_resolution_Below0p99", "hLead_matchLeptons_2d_resolution_Below0p99", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hLead_matchLeptons_2d_resolution_Below0p99);  //2d resolution
    
hSubLead_matchLeptons_1d_resolution_Below0p99 = new TH1F("hSubLead_matchLeptons_1d_resolution_Below0p99", "hSubLead_matchLeptons_1d_resolution_Below0p99", 500, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_resolution_Below0p99); //1d resolution
    
hSubLead_matchLeptons_2d_resolution_Below0p99 = new TH2F("hSubLead_matchLeptons_2d_resolution_Below0p99", "hSubLead_matchLeptons_2d_resolution_Below0p99", 50, 0., 2., 50, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_2d_resolution_Below0p99);
  
//DOUBLE COUNTING HISTOGRAMS
hLead_matchLeptons_1d_resolution_doubleCount = new TH1F("hLead_matchLeptons_1d_resolution_doubleCount", "hLead_matchLeptons_1d_resolution_doubleCount", 100, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_resolution_doubleCount); //1d resolution
    
hLead_matchLeptons_2d_resolution_doubleCount = new TH2F("hLead_matchLeptons_2d_resolution_doubleCount", "hLead_matchLeptons_2d_resolution_doubleCount", 100, 0., 2., 10, 0., 2.);
wk()->addOutput(hLead_matchLeptons_2d_resolution_doubleCount);  //2d resolution
        
    
hSubLead_matchLeptons_1d_resolution_doubleCount = new TH1F("hSubLead_matchLeptons_1d_resolution_doubleCount", "hSubLead_matchLeptons_1d_resolution_doubleCount", 100, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_1d_resolution_doubleCount); //1d resolution
    
hSubLead_matchLeptons_2d_resolution_doubleCount = new TH2F("hSubLead_matchLeptons_2d_resolution_doubleCount", "hSubLead_matchLeptons_2d_resolution_doubleCount", 100, 0., 2., 100, 0., 2.);
wk()->addOutput(hSubLead_matchLeptons_2d_resolution_doubleCount);  //2d resolution
    
hLead_matchLeptons_1d_resolution_Strip = new TH1F("hLead_matchLeptons_1d_resolution_Strip", "hLead_matchLeptons_1d_resolution_Strip", 100, 0., 2.);
wk()->addOutput(hLead_matchLeptons_1d_resolution_Strip); //1d resolution
    
hLead_matchLeptons_2d_resolution_Strip = new TH2F("hLead_matchLeptons_2d_resolution_Strip", "hLead_matchLeptons_2d_resolution_Strip", 100, 0., 2., 100, 0., 2.);
wk()->addOutput(hLead_matchLeptons_2d_resolution_Strip);  //2d resolution

m_el_truthMatched = new TH1F("m_el_truthMatched", "m_el_truthMatched", 4, -0.5, 3.5);
wk()->addOutput(m_el_truthMatched);   
    
//GSF track collection
m_GSFtrk_pT = new TH1F("m_GSFtrk_pT", "m_GSFtrk_pT", 50, 0., 3000);
wk()->addOutput(m_GSFtrk_pT); 
    
m_GSFtrk_eta = new TH1F("m_GSFtrk_eta", "m_GSFtrk_eta", 25, -2.5, 2.5);
wk()->addOutput(m_GSFtrk_eta); 
    
//ID track collection
m_IDtrk_pT = new TH1F("m_IDtrk_pT", "m_IDtrk_pT", 50, 0., 3000);
wk()->addOutput(m_IDtrk_pT); 
    
m_IDtrk_eta = new TH1F("m_IDtrk_eta", "m_IDtrk_eta", 25, -2.5, 2.5);
wk()->addOutput(m_IDtrk_eta); 
    
m_dEta_trueLeptons = new TH1F("m_dEta_trueLeptons", "m_dEta_trueLeptons", 40, -1., 1.);
wk()->addOutput(m_dEta_trueLeptons); 
    
m_dPhi_trueLeptons = new TH1F("m_dPhi_trueLeptons", "m_dPhi_trueLeptons", 40, -1., 1.);
wk()->addOutput(m_dPhi_trueLeptons); 
    
m_dR_realLeptons = new TH1F("m_dR_realLeptons", "m_dR_realLeptons", 100, 0.0, 1.0);
wk()->addOutput(m_dR_realLeptons);     
    
    
m_trueSubLead_pT_noMatched = new TH1F("m_trueSubLead_pT_noMatched", "m_trueSubLead_pT_noMatched", 50, 0., 3000);
wk()->addOutput(m_trueSubLead_pT_noMatched); 

    
    
    
    
}

//=======================
//RECONSTRUCTION STUDIES
//======================
vector<xAOD::Electron*> H2ZyAnalysis::match_DR_electrons(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont)
{

    double drMin = 10000;
    xAOD::Electron* minDRElectron;
                
    for ( auto el : recoElCont ) //matching to reco electrons
    {
        double dEta = true_lepton->eta() - el->eta();
        double dPhi = true_lepton->phi() - el->phi();
                    
        dPhi = (dPhi <= M_PI) ? dPhi : 2*M_PI-dPhi;
        double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );

        if( dR < drMin )
        {
            drMin = dR;
            minDRElectron = el;
        } 
    }
                
    vector<xAOD::Electron*> tmp;
    if(drMin < 0.1) tmp.push_back(minDRElectron);
                
    return tmp;
}

vector<xAOD::Jet*> H2ZyAnalysis::match_DR_jet(const xAOD::TruthParticle* true_lepton, xAOD::JetContainer recoJetCont)
{

    double drMin = 10000;
    xAOD::Jet* minDRJet;
                
    for ( auto jet : recoJetCont ) //matching to reco electrons
    {
        double dEta = true_lepton->eta() - jet->eta();
        double dPhi = true_lepton->phi() - jet->phi();
                    
        dPhi = (dPhi <= M_PI) ? dPhi : 2*M_PI-dPhi;
        double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
                            
        if( dR < drMin )
        {
            drMin = dR;
            minDRJet = jet;
        } 
    }
             
    vector<xAOD::Jet*> tmp;
    if(drMin < 0.3) tmp.push_back(minDRJet);
                
    return tmp;
}

xAOD::ElectronContainer H2ZyAnalysis::skipRecoMatchedToTrueLeading(xAOD::Electron* el_matchedLeading, xAOD::ElectronContainer recoElCont)
{
    //cout<<" el container BEFORE remove truth matched el: " << recoElCont.size() << endl;
    xAOD::ElectronContainer el_skip(SG::VIEW_ELEMENTS);
    for ( auto el : recoElCont )
    {
        if( el == el_matchedLeading ) continue;
        else el_skip.push_back(el);
    }

    //cout<<" el container AFTER remove truth matched el: " << el_skip.size() << endl;
    return el_skip;
}

bool H2ZyAnalysis::checkDoubleCountingEvents(xAOD::Electron* el_matchedLeading, const xAOD::TruthParticle* true_lepton)
{
    
    
    double dEta = true_lepton->eta() - el_matchedLeading->eta();
    double dPhi = true_lepton->phi() - el_matchedLeading->phi(); 
    dPhi = (dPhi <= M_PI) ? dPhi : 2*M_PI - dPhi;
    double dR = sqrt( pow(dEta,2) + pow(dPhi,2) );
    
    
    
    if(dR < 0.1)
    {
        return true;
    } 
    else return false;
    
    
    
}
bool H2ZyAnalysis::isTruthLepton(const xAOD::Electron *el)
{
  const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*el);
  if(!truthlepton ||abs(truthlepton->pdgId())!=11)
  {
    //cout<<" broken link or no truth electron "<<endl;
    return false;    
  } 
  else return true;
}
bool H2ZyAnalysis::truePassPreCut(const xAOD::TruthParticle* true_lepton)
{
    if( true_lepton->pt()/1000 > 10. && fabs(true_lepton->eta()) < 2.47 )
        return true;
    else return false;
}

bool H2ZyAnalysis::trueMatch_toolFunction(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont)
{
 
    bool tmp_el = false;
    
    for ( auto el : recoElCont ) //matching to reco electrons
    {
        const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*el);
        if(truthlepton)
        {
            if(abs(truthlepton->pdgId()) == abs(true_lepton->pdgId()) && truthlepton->barcode() == true_lepton->barcode() && truthlepton->status() == true_lepton->status())
            {
                tmp_el = true;
                break;
            }
        }
    }
    
    if(tmp_el) return true;
    else return false;
}

bool H2ZyAnalysis::trueMatch_BrokenLink(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont)
{
 
    bool broken_el = false;
    
    for ( auto el : recoElCont ) //matching to reco electrons
    {
        const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*el);
        if(truthlepton) broken_el = true;
    }
    
    if(broken_el) return true;
    else return false;
}

vector<const xAOD::TruthParticle*> H2ZyAnalysis::trueMatch_noMatched(const xAOD::TruthParticle* true_lepton, xAOD::ElectronContainer recoElCont)
{
 
    bool tmp_match_el = false;
    vector<const xAOD::TruthParticle*> tmp_vec;
    
    for ( auto el : recoElCont ) //matching to reco electrons
    {
        const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*el);
        if(truthlepton)
        {
            if(abs(truthlepton->pdgId()) == abs(true_lepton->pdgId()) && truthlepton->barcode() == true_lepton->barcode() && truthlepton->status() == true_lepton->status())
            {
                tmp_match_el = true;
                break;
            }
            else
            {
                cout<<" link exists, no matched " << " true info " << truthlepton->pdgId() << " " << truthlepton->barcode() << " " << truthlepton->status() << endl;
                tmp_vec.push_back(truthlepton);
            }
        }
    }
    
    return tmp_vec;
}

xAOD::ElectronContainer H2ZyAnalysis::getRecoMatched( vector<const xAOD::TruthParticle*> lepton_true, xAOD::ElectronContainer recoElCont )
{
 
    xAOD::ElectronContainer recoMatched(SG::VIEW_ELEMENTS);
    
    for(auto true_el : lepton_true)
    {
        for ( auto el : recoElCont )
        {
            const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*el);
            
            if(truthlepton)
            {
                if(abs(truthlepton->pdgId()) == abs(true_el->pdgId()) && truthlepton->barcode() == true_el->barcode() && truthlepton->status() == true_el->status())
                {
                    recoMatched.push_back(el);
                } 
            }
        }
    }

    return recoMatched;
}

bool H2ZyAnalysis::findTracksMatchingSameElectron(const xAOD::TruthParticle* truelepton_l1, const xAOD::TruthParticle* truelepton_l2, xAOD::ElectronContainer recoElCont)
{
    bool tmp_l1 = false;
    bool tmp_l2 = false;
    bool matchBothTrueEl = false;
    
    for( auto el : recoElCont ) //matching to reco electrons
    {
        for(unsigned int i = 0; i < el->nTrackParticles(); ++i)
        {
            const xAOD::TrackParticle* el_tp = el->trackParticle(i);   
            if (!el_tp) continue;
            
            const xAOD::TruthParticle* link_truep = xAOD::TruthHelpers::getTruthParticle(*el_tp);
            if(link_truep)
            {
                if(abs(truelepton_l1->pdgId()) == abs(link_truep->pdgId()) && truelepton_l1->barcode() == link_truep->barcode() && truelepton_l1->status() == link_truep->status())
                {
                    cout<<" Matching true leading " << endl;
                    tmp_l1 = true;
                }

                if(abs(truelepton_l2->pdgId()) == abs(link_truep->pdgId()) && truelepton_l2->barcode() == link_truep->barcode() && truelepton_l2->status() == link_truep->status())
                {
                    cout<<" Matching true sub-leading " << endl;
                    tmp_l2 = true;
                }

                if(tmp_l1 && tmp_l2)
                {
                    cout<<" both true leptons match the same track! " << endl;
                    matchBothTrueEl = true;
                } 
            }
            
        } 
    }
    
    if(matchBothTrueEl) return true;
    else return false;
}

const xAOD::IParticle* H2ZyAnalysis::getTheOriginalPointer(const xAOD::IParticle& part)
{
  // Find the original pointer (if any) -- in order to work with deep copies
  const xAOD::IParticle* compPointer = &part;
  if (xAOD::getOriginalObject(part)) {
    compPointer = xAOD::getOriginalObject(part);
  }
  return compPointer;
}

float H2ZyAnalysis::getTruthMatchProbability(const xAOD::TrackParticle* trackParticle)
{
  float truthProb = 0.;
  if (trackParticle->isAvailable<float>("truthMatchProbability")) {
    truthProb = trackParticle->auxdata<float>("truthMatchProbability");
  }
  return truthProb;

}



xAOD::TrackParticle* H2ZyAnalysis::foundIDTrack(const xAOD::TruthParticle* tlep)
{
    const xAOD::TrackParticleContainer* IDtracks = 0;
    if(!event()->retrieve(IDtracks, "InDetTrackParticles" ).isSuccess())
    {
        Error("execute()", "Failed to retrieve InDetTrackParticles container." );
    }
                
    cout<<" InDetTrackParticles size " << IDtracks->size() << endl;
    vector<xAOD::TrackParticle*> tmp_trk;
    int nMatch = 0;
    for(auto trk : *IDtracks)
    {
        //float truth_prob = trk->auxdata<float>("truthMatchProbability");
        //cout<<" get prob " << truth_prob << endl;
        xAOD::TrackParticle* tp_p = (xAOD::TrackParticle*)getTheOriginalPointer(*trk);
        const xAOD::TruthParticle* truthlepton = xAOD::TruthHelpers::getTruthParticle(*tp_p);
        if(!truthlepton) continue;
                    
        const xAOD::TruthParticle* truthPart_p = (xAOD::TruthParticle*)getTheOriginalPointer(*truthlepton);
        const xAOD::TruthParticle* truthLep_p = (xAOD::TruthParticle*)getTheOriginalPointer(*tlep);
        
        
        cout<<" truthPart_p " << truthPart_p << " tlep " << tlep << " truthLep_p " << truthLep_p << endl;
                    
        if(abs(truthPart_p->pdgId()) == abs(tlep->pdgId()) && truthPart_p->barcode() == tlep->barcode() && truthPart_p->status() == tlep->status())
        {
            tmp_trk.push_back(tp_p);
            nMatch++;
        }
    }
    
    cout<<" nMatch " << nMatch << endl;
    
    xAOD::TrackParticle* trk;
    double large_pT = 0;
    
    if(tmp_trk.size() == 1)
    {
        trk = tmp_trk.at(0);
    }
    else if(tmp_trk.size() >= 2)
    {
        for(auto tk : tmp_trk)
        {
            cout<<"tk pt " << tk->pt()/1000 << endl;
            if( tk->pt()/1000 > large_pT )
            {
                large_pT = tk->pt()/1000;
                trk = tk;
            }
        }
    }

    cout<<" track pT " << trk->pt()/1000 << endl;
    return trk;
    
}


#include "H2Zy/TruthSelect.h"
#include "MCUtils/PIDUtils.h"

namespace HZG {

	TruthSelect::TruthSelect() {
	  useGeV(true);
	}

	void TruthSelect::setMVA(){
	  reader = new TMVA::Reader( "!Color:!Silent" );
	  TString infile_mva = (TString)PathResolverFindCalibFile("H2Zy/TMVAClassification_BDTG.weights.xml");
	  reader->AddVariable("VBF_Dphi_Zy_jj", &Dphi_Zy_jj);
	  reader->AddVariable("VBF_Zepp", &eta_Zepp);
	  reader->AddVariable("VBF_DRmin_y_j", &DRmin_Z_j);
	  reader->AddVariable("VBF_m_jj", &m_jj);	  
	  reader->AddVariable("llg_pTt_Zmassconstraint", &pTt);
	  reader->AddVariable("VBF_Dy_j_j", &Dy_j_j);
	  reader->AddVariable("llg_dphi_Zy_Zmassconstraint", &dphi_Zy);

	  reader->BookMVA("BDTG",infile_mva);

	}
	//---------------- initialize ----------------------------------
	void TruthSelect::initialize(const HG::Config &config,  xAOD::TEvent *_event,  Bmap *_beventmap, Imap *_ieventmap, Fmap *_feventmap, IVmap *_iveventmap, FVmap *_fveventmap,  TLVVmap *_tlvveventmap){
		m_config = config;
		m_event = _event;
		m_beventmap = _beventmap;
		m_ieventmap = _ieventmap;
		m_feventmap = _feventmap;
		m_iveventmap = _iveventmap;
		m_fveventmap = _fveventmap;
		m_tlvveventmap = _tlvveventmap;

		TString m_truthParticlesname = m_config.getStr("TruthHandler.truthParticles", "TruthParticles");
		if(!m_event->retrieve(m_truthParticles, m_truthParticlesname.Data() ).isSuccess()){
			Error("execute()", "Failed to retrieve TruthParticles container." );
		}

		TString m_TruthVertexname = m_config.getStr("TruthHandler.Vertices", "TruthVertices");
		if(!m_event->retrieve(m_Vertices, m_TruthVertexname.Data() ).isSuccess()){
			Error("execute()", "Failed to retrieve TruthVertices container." );
		}

		TString m_TruthJetname = m_config.getStr("TruthHandler.JetContainerName", "AntiKt4TruthJets");
		if(!m_event->retrieve(m_antiKt4TruthJets, m_TruthJetname.Data() ).isSuccess()){
		  Error("execute()", "Failed to retrieve AntiKt4TruthJets container." );
		}

		//    if(!m_event->retrieve(m_antiKt4TruthWZJets, "AntiKt4TruthWZJets" ).isSuccess()){
		//      Error("execute()", "Failed to retrieve AntiKt4TruthWZJets container. Setting use WZ jets = false." );
		//    }

	}
	
	bool TruthSelect::notFromHadronic(const xAOD::TruthParticle *ptcl) {
		int ID = ptcl->pdgId();

		// if the particle is a hadron, return false
		if (MCUtils::PID::isHadron(ID)) return false;
		//if (MC::PID::isHadron(ID)) return false;
		if (abs(ID)>100) return false;

		// if there are no parents, not from hadron
		if (ptcl->nParents()==0) return true;

		const xAOD::TruthParticle *parent = ptcl->parent(0);
		int parentID = parent->pdgId();
		if (MCUtils::PID::isHadron(parentID)) return false; // from hadron!
		//if (MC::PID::isHadron(parentID)) return false; // from hadron!
		if (abs(parentID)>100) return false; // from hadron!
		if (abs(parentID)==15||parentID==ID) return notFromHadronic(parent);

		// if we get here, all is good
		return true;
	}
	//---------------------- dumptruth ----------------------------------


    void TruthSelect::dumptruth(){

                int iparticle=0;
                for ( const xAOD::TruthParticle *ptcl : *m_truthParticles)
                {       
                        //    if (ptcl->pdgId()==25&&ptcl->status()>60) HG::printTruthPtcl(ptcl,"The higgs");
                        //if(ptcl->barcode() > 20000) continue;
                        //if(ptcl->status() != 1) continue;
                        //if(fabs(ptcl->pdgId())>100) continue;
                        //if(! notFromHadronic(ptcl)) continue; 
                        //if (ptcl->nParents()!=0) {if(ptcl->parent(0)->pdgId()>100) continue;}
                        
                        cout<<iparticle<<endl;
                        cout<<"pid="<<ptcl->pdgId()<<" pt="<<ptcl->pt()<<" eta="<<ptcl->eta()<<" phi="<<ptcl->phi()<<" barcode="<<ptcl->barcode()<<" status="<<ptcl->status()<<" m="<<ptcl->m()<<" nChild="<<ptcl->nChildren()<<" e="<<ptcl->e()<<" p="<<ptcl->p4().P()<<" et="<<ptcl->p4().Et()<<" e/p="<<ptcl->e()/ptcl->p4().P()<<" et/pt="<<ptcl->p4().Et()/ptcl->pt();
                        
                        /*
                        if(ptcl->barcode()==806&&ptcl->nChildren()==2) {
                        const xAOD::TruthParticle *child0 = ptcl->child(0);
                        const xAOD::TruthParticle *child1 = ptcl->child(1);
                        if(child0!=nullptr) std::cout<<" child0 barcode=" << child0->barcode() <<" pid="<<child0->pdgId()<<" status="<<child0->status() ;
                        if(child1!=nullptr) std::cout<<" child1 barcode=" << child1->barcode() <<" pid="<<child1->pdgId()<<" status="<<child1->status() ;                       
                        }
                        */

                        
                        if (ptcl->nParents()==0) cout<<" parent=-1"<<endl;
                        else {
                          const xAOD::TruthParticle *parent = ptcl->parent(0);
                          if(parent!=nullptr){
                            float decayR = -1;
                            if(parent->decayVtx()){
                              decayR = parent->decayVtx()->perp() ;
                            }
                            std::cout<<" parent barcode="<<parent->barcode()<<" pid="<<parent->pdgId()<<" status="<<parent->status()<< " decayR=" << decayR <<std::endl;
                          }
                        }

                        
                        iparticle++;
                }

        }
    
    vector<const xAOD::TruthParticle*> TruthSelect::TruthLeptonsFromZ()
    {
    
        int iparticle=0;
        vector<const xAOD::TruthParticle*> tmp_vector;
        vector<const xAOD::TruthParticle*> truthElsZ;
        for ( const xAOD::TruthParticle *ptcl : *m_truthParticles)
        {      
                
            double isStable = false;
            if( ptcl->status() == 1 && ptcl->barcode() < 200000 ) isStable = true;
            
            if(isStable)
            {
                if (abs(ptcl->pdgId())== 11 && ptcl->nParents()!=0) 
                {
                    if(fabs(ptcl->parent(0)->pdgId()) > 100) continue;
                            
                    if(ptcl->parent(0)->pdgId() == 23)
                    {
                        tmp_vector.push_back(ptcl);
                    }
                    else
                    {
                        int ancestorID = AncestorID(ptcl);
                        
                        if( ancestorID == 23 )
                        {
                            tmp_vector.push_back(ptcl);
                                        
                        }
                        else cout<<" no anscestor ID 23 " << endl;
                    }
                } 
            }
                        iparticle++;
        }
        
        vector<const xAOD::TruthParticle*> temp = GetLargestpTLeptons(tmp_vector);
        truthElsZ = temp;
        
        return truthElsZ; 
    }
    
    vector<const xAOD::TruthParticle*> TruthSelect::GetLargestpTLeptons(vector<const xAOD::TruthParticle*> vec)
    {
        
        vector<const xAOD::TruthParticle*> out_vec;
        
        double l1 = 0; double l2 = 0;
        int pos_l1 = 0; int pos_l2 = 0;
        
        for(int i = 0; i < vec.size(); i++)
        {   
            if (vec.at(i)->pt()/1000 > l1)
              {     
                    l1 = vec.at(i)->pt()/1000;
                    pos_l1 = i;
              }
        }
        for(int i = 0; i < vec.size(); i++)
        {   
            if (vec.at(i)->pt()/1000 > l2)
              {
                    
                    if (vec.at(i)->pt()/1000 == l1) continue;
                    l2 = vec.at(i)->pt()/1000;
                    pos_l2 = i;
              }
        }
        
        out_vec.push_back(vec.at(pos_l1));
        out_vec.push_back(vec.at(pos_l2));
        
        return out_vec;
     
    }
    
    
	void TruthSelect::searchDphoton()
	{
		hasDphoton=0;
		hasttyphoton=0;

		for ( const xAOD::TruthParticle *ptcl : *m_truthParticles)
		{
			if(ptcl->barcode() > 200000) continue;
			if(ptcl->pdgId() != 22 ) continue; 
			if(ptcl->parent(0)!=nullptr && ptcl->pt() > 10e3) {
			  if (abs(ptcl->parent(0)->pdgId()) == 24 || abs(ptcl->parent(0)->pdgId()) == 6)
			    hasttyphoton = 1;
			}
			if(!notFromHadronic(ptcl)) continue;
			if(HG::isGoodTruthPhoton(ptcl)) hasDphoton=1;
		}
	}

	void TruthSelect::truthsave(unsigned int _eventNumber, unsigned int _mc_channel_number){

		m_eventNumber = _eventNumber;
	        
	        //----- save truth hard scattering vertex information
		truthsave_vertex();

		//----- save truth Higgs information (if present) -----
		(*m_feventmap)["mc_Higgs_m"] = 100;
		(*m_feventmap)["mc_Higgs_pt"] = -99;
		(*m_feventmap)["mc_Higgs_eta"] = -99;
		(*m_feventmap)["mc_Higgs_phi"] = -99;
		HG::TruthPtcls ys  = findHiggs();
		if(ys.size()>0) 
		{
			(*m_feventmap)["mc_Higgs_m"]=ys[0]->p4().M()/mass_scale;
			(*m_feventmap)["mc_Higgs_pt"]=ys[0]->p4().Pt()/mass_scale;
			(*m_feventmap)["mc_Higgs_eta"]=ys[0]->p4().Eta();
			(*m_feventmap)["mc_Higgs_phi"]=ys[0]->p4().Phi();
		}
		higgs_truth_mass = (*m_feventmap)["mc_Higgs_m"];
		higgs_truth_pt = (*m_feventmap)["mc_Higgs_pt"];
		higgs_truth_eta = (*m_feventmap)["mc_Higgs_eta"];
		higgs_truth_phi = (*m_feventmap)["mc_Higgs_phi"];

		// ---- save the stable particle information and the true Z mass
		(*m_tlvveventmap)["mc_SP_p4"].clear();
		(*m_iveventmap)["mc_SP_barcode"].clear();
		(*m_iveventmap)["mc_SP_pdgId"].clear();
		(*m_feventmap)["mc_Z_m"] = -999;
		TLorentzVector p_leadinglepton, p_sublepton;

		//dumptruth();
		bool isSherpa3 = true; // whether use status=3 parton in sherpa or not
		isSherpa3 = false;

		bool photon = false; bool lm = false; bool lp = false; bool higgs = false;
		p_leadinglepton.SetPtEtaPhiM(0,0,0,0);
		p_sublepton.SetPtEtaPhiM(0,0,0,0);
		p_ph.SetPtEtaPhiM(0,0,0,0);
		j1.SetPtEtaPhiM(0,0,0,0);
		j2.SetPtEtaPhiM(0,0,0,0);

		int Njet_status3 = 0;
		for ( const xAOD::TruthParticle *ptcl : *m_truthParticles)
		 {
		  if(ptcl->barcode() > 200000) continue;
		  if((ptcl->status() == 11 || ptcl->status() == 22) && (ptcl->pdgId()==23) ) 
		   {
		   bool isfromhiggs = ZIsFromHiggs(ptcl);
		   if(!isfromhiggs) continue;
		   (*m_feventmap)["mc_Z_m"] =ptcl->p4().M()/mass_scale;
		   Z_truth_pt = ptcl->p4().Pt()/mass_scale;
		   Z_truth_eta = ptcl->p4().Eta();
		   Z_truth_phi = ptcl->p4().Phi();
		   Z_truth_mass = ptcl->p4().M()/mass_scale;
		   retrieve_Zdecay_topo(ptcl);
		   if(p1.Pt()>p2.Pt()) {p_leadinglepton=p1, p_sublepton = p2;}
		   else {p_leadinglepton=p2, p_sublepton = p1; int pid = l1_truth_pdgId; l1_truth_pdgId = l2_truth_pdgId; l2_truth_pdgId = pid;}
		   l1_truth_pt = p_leadinglepton.Pt()/mass_scale, l1_truth_eta = p_leadinglepton.Eta(), l1_truth_phi = p_leadinglepton.Phi(), l1_truth_mass = p_leadinglepton.M()/mass_scale;
		   l2_truth_pt = p_sublepton.Pt()/mass_scale, l2_truth_eta = p_sublepton.Eta(), l2_truth_phi = p_sublepton.Phi(), l2_truth_mass = p_sublepton.M()/mass_scale;
		   lm = true; lp = true; 
		  }

		 if(ptcl->pdgId() == 25 && (ptcl->status() == 11||ptcl->status()==62))
		  {
		   retrieve_Higgsdecay_topo(ptcl);
		   ph_truth_pt = p_ph.Pt()/mass_scale, ph_truth_eta = p_ph.Eta(), ph_truth_phi = p_ph.Phi();
		   photon = true;
		   higgs = true;
		  }
		 
		 if(ptcl->status() ==3 && fabs(ptcl->pdgId())!=11 && fabs(ptcl->pdgId())!=13 && fabs(ptcl->pdgId())!=22&& ptcl->pt()!=0) {
		   //std::cout<<"a special parton:" << ptcl->pdgId() << " " << ptcl->pt()<< " " << ptcl->eta()<<" " <<ptcl->phi() << std::endl;
		   Njet_status3++;

		   if(Njet_status3==1&&ptcl->pt()>j1.Pt()) j1 = ptcl->p4();
		   else if (Njet_status3>=2&&ptcl->pt()>j1.Pt()) {j2=j1; j1 = ptcl->p4();}
		   else if (Njet_status3>=2&&ptcl->pt()>j2.Pt()) j2 = ptcl->p4();		   
		 }


		 //if(isSherpa3&&ptcl->status() != 3) continue;
		 //else 
		   if(ptcl->status() != 1) continue;

		 if(fabs(ptcl->pdgId())>100) continue;

		 //find and save objects for background MC
		 //std::cout << "test = " << p_ph.Pt() <<  " " << p_leadinglepton.Pt() << " " << p_sublepton.Pt() <<std::endl;

		 if(p_ph.Pt()==0&&!photon && ptcl->pdgId() == 22&& HG::isGoodTruthPhoton(ptcl)) {
		   p_ph.SetPtEtaPhiM(ptcl->pt(), ptcl->eta(), ptcl->phi(), ptcl->p4().M());
		   photon = true;
		 }
		 if(p_leadinglepton.Pt()==0&&!lp && (ptcl->pdgId() == 11 || ptcl->pdgId() == 13)) {
		   l1_truth_pdgId = ptcl->pdgId();
		   p_leadinglepton.SetPtEtaPhiM(ptcl->pt(), ptcl->eta(), ptcl->phi(), ptcl->p4().M());
		   lp = true;
		 }
		 if(p_sublepton.Pt()==0&&!lm && (ptcl->pdgId() == -11 || ptcl->pdgId() == -13)) {
		   l2_truth_pdgId = ptcl->pdgId();
		   p_sublepton.SetPtEtaPhiM(ptcl->pt(), ptcl->eta(), ptcl->phi(), ptcl->p4().M());
		   lm = true;
		 }

		 if(p1.Pt()<p2.Pt()) {
		   TLorentzVector p_lepton; p_lepton = p_sublepton;
		   p_sublepton = p_leadinglepton; p_leadinglepton = p_lepton;
		   int l_pdgId; l_pdgId = l2_truth_pdgId;
		   l2_truth_pdgId = l1_truth_pdgId; l1_truth_pdgId = l_pdgId;
		 }
		 		 
		 // giving values if it's really background MC
		 if(!higgs){
		   l1_truth_pt = p_leadinglepton.Pt()/mass_scale, l1_truth_eta = p_leadinglepton.Eta(), l1_truth_phi = p_leadinglepton.Phi(), l1_truth_mass = p_leadinglepton.M()/mass_scale;
		   l2_truth_pt = p_sublepton.Pt()/mass_scale, l2_truth_eta = p_sublepton.Eta(), l2_truth_phi = p_sublepton.Phi(), l2_truth_mass = p_sublepton.M()/mass_scale;
		   ph_truth_pt = p_ph.Pt()/mass_scale, ph_truth_eta = p_ph.Eta(), ph_truth_phi = p_ph.Phi(), ph_truth_mass = p_ph.M()/mass_scale;

		   TLorentzVector p_Z = p_leadinglepton + p_sublepton;		 
		   Z_truth_pt = p_Z.Pt()/mass_scale, Z_truth_eta = p_Z.Eta(), Z_truth_phi = p_Z.Phi(), Z_truth_mass = p_Z.M()/mass_scale;

		   TLorentzVector p_H = p_Z + p_ph;
		   higgs_truth_pt = p_H.Pt()/mass_scale, higgs_truth_eta =p_H.Eta() ,higgs_truth_phi = p_H.Phi(),higgs_truth_mass = p_H.M()/mass_scale;
		 }

		 (*m_feventmap)["Z_truth_mass"] = -999.;
		 (*m_feventmap)["Z_truth_mass"] = Z_truth_mass;
		 //??
		 if(! HG::notFromHadron(ptcl)) continue;
		 TLorentzVector TLV_SP_Moment;
		 TLV_SP_Moment.SetPtEtaPhiM(ptcl->pt(), ptcl->eta(), ptcl->phi(), ptcl->p4().M());
		 (*m_tlvveventmap)["mc_SP_p4"].push_back(TLV_SP_Moment);
		 (*m_iveventmap)["mc_SP_barcode"].push_back(ptcl->barcode());
		 (*m_iveventmap)["mc_SP_pdgId"].push_back(ptcl->pdgId());
		}

		N_jets = Njet_status3;
		j1_truth_pt = j1.Pt()/mass_scale;
		j2_truth_pt = j2.Pt()/mass_scale;

		Dphi_Zy_jj = -999;
		pTt = -999;
		eta_Zepp = -999;
		m_jj = -999;
		DRmin_Z_j = -999;
		Dy_j_j = -999;
		mvaout = -999;
		dphi_Zy = -999;

		if(!isSherpa3) retrieve_truthjets(p_leadinglepton, p_sublepton, p_ph);

		(*m_ieventmap)["Truth_N_j"] = -999;
		(*m_ieventmap)["Truth_N_j"] = N_jets;
		(*m_feventmap)["j1_pt_Truth"] = -999;
		(*m_feventmap)["j2_pt_Truth"] = -999;
		(*m_feventmap)["j1_eta_Truth"] = -999;
		(*m_feventmap)["j2_eta_Truth"] = -999;
		(*m_feventmap)["j1_phi_Truth"] = -999;
		(*m_feventmap)["j2_phi_Truth"] = -999;
		(*m_feventmap)["j1_m_Truth"] = -999;
		(*m_feventmap)["j2_m_Truth"] = -999;

		if(N_jets>=1){
		  (*m_feventmap)["j1_pt_Truth"] = j1.Pt()/mass_scale;
		  (*m_feventmap)["j1_eta_Truth"] = j1.Eta();
		  (*m_feventmap)["j1_phi_Truth"] = j1.Phi();
		  (*m_feventmap)["j1_m_Truth"] = j1.M()/mass_scale;
		}

		if(doDphi) {
		  (*m_feventmap)["j2_pt_Truth"] = j2.Pt()/mass_scale;
		  (*m_feventmap)["j2_eta_Truth"] = j2.Eta();
		  (*m_feventmap)["j2_phi_Truth"] = j2.Phi();
		  (*m_feventmap)["j2_m_Truth"] = j2.M()/mass_scale;

		  Dphi_Zy_jj = fabs((p_leadinglepton + p_sublepton + p_ph).DeltaPhi(j1 + j2));
		  if(Dphi_Zy_jj>2.94) Dphi_Zy_jj=2.94;
		  m_jj = (j1 + j2).M()/mass_scale;
		  Dy_j_j = fabs(j1.Eta()-j2.Eta());
		  pTt = fabs((p_leadinglepton + p_sublepton).Px()*p_ph.Py()-p_ph.Px()*(p_leadinglepton + p_sublepton).Py())/(p_leadinglepton + p_sublepton - p_ph).Pt()*2.0/mass_scale;
		  eta_Zepp = fabs((p_leadinglepton + p_sublepton + p_ph).Eta() - (j1.Eta() + j2.Eta())/2.0);
		  DRmin_Z_j = TMath::Min((p_leadinglepton + p_sublepton).DeltaR(j1), TMath::Min((p_leadinglepton + p_sublepton).DeltaR(j2), TMath::Min(p_ph.DeltaR(j1),p_ph.DeltaR(j2))));
		  dphi_Zy = fabs((p_leadinglepton + p_sublepton).DeltaPhi(p_ph));
		}

		mvaout = reader->EvaluateMVA("BDTG");

                HZG::Object_llg::CalculateCommonKinematic(p_ph, p_leadinglepton, p_sublepton, (*m_feventmap), "_Truth");
		HZG::Object_llg::CalculateCommonKinematic_LL(p_leadinglepton, p_sublepton, (*m_feventmap), "_Truth");
                HZG::Object_llg::CalculateCMSKinematic(p_ph, p_leadinglepton, p_sublepton, (*m_feventmap), "_Truth");
		m_truthTree->Fill();
	}

        void TruthSelect::useGeV(bool _useGeV) {
	  if (_useGeV)
	    mass_scale = GEV;
	  else
	    mass_scale = MEV;
        }

	void TruthSelect::print_child(const xAOD::TruthParticle *ptcl){
		cout<<"n child : "<<ptcl->nChildren()<<endl; 
		for(int ichild=0; ichild<(int)ptcl->nChildren(); ichild++){
			cout<<ptcl->child(ichild)->pdgId()<<" "<<ptcl->child(ichild)->barcode()<<endl;
			//if(ptcl->child(ichild)->barcode()>10000) continue;
		}
	}

	bool TruthSelect::ZIsFromHiggs(const xAOD::TruthParticle *ptcl){

	  int ID = ptcl->pdgId();
	  const xAOD::TruthParticle *father;
	  const xAOD::TruthParticle *grandfather;
	  if (ptcl->nParents()!=0) { father = ptcl->parent(0); ID = father->pdgId();} else return ID;
	  if (father->nParents()!=0) { grandfather = father->parent(0); ID = grandfather->pdgId();} else return ID;
	  while (grandfather->nParents()!=0){
	    father = grandfather->parent(0);
	    grandfather = father ;
	    if (ID==25) break;
	    ID = grandfather->pdgId();
	  }
	  bool ZIsFromHiggs = false;
	  if (ID == 25) ZIsFromHiggs = true;

	  return ZIsFromHiggs;
	}

	void TruthSelect::retrieve_Higgsdecay_topo(const xAOD::TruthParticle *ptcl){
		(*m_ieventmap)["mc_isDalitz"] = 0;
		int parent_pid = ptcl->pdgId();
		for(int ichild=0; ichild<(int)ptcl->nChildren(); ichild++){
			//if(ptcl->child(ichild)->barcode()>10000) continue;
			if(ptcl->child(ichild)->barcode()>20000) continue;
			int _child_pid = ptcl->child(ichild)->pdgId();
			if(_child_pid==parent_pid)
			{
				retrieve_Higgsdecay_topo(ptcl->child(ichild));
				break;
			}
			//std::cout<<"higgs child status "<<ptcl->child(ichild)->status()<<" "<<ptcl->child(ichild)->pdgId()<<" "<<ptcl->child(ichild)->barcode()<<std::endl;
			if(_child_pid==22 && (ptcl->child(ichild)->status()==23 || ptcl->child(ichild)->status()==11)) (*m_ieventmap)["mc_isDalitz"] = 1;
			//if(_child_pid==22 && (ptcl->child(ichild)->status()==23)) (*m_ieventmap)["mc_isDalitz"] = 1;
			if(_child_pid==22) p_ph.SetPtEtaPhiM(ptcl->child(ichild)->pt(), ptcl->child(ichild)->eta(), ptcl->child(ichild)->phi(), ptcl->child(ichild)->p4().M());
		}

	}

	void TruthSelect::retrieve_Zdecay_topo(const xAOD::TruthParticle *ptcl){
		(*m_ieventmap)["mc_Z_decay_topo"] = 0;
		int parent_pid = ptcl->pdgId();
		for(int ichild=0; ichild<(int)ptcl->nChildren(); ichild++){
			//if(ptcl->child(ichild)->barcode()>10000) continue; 
			if(ptcl->child(ichild)->barcode()>20000) continue;
			int _child_pid = ptcl->child(ichild)->pdgId();
			if(_child_pid==parent_pid) 
			{
				(*m_ieventmap)["mc_Z_decay_topo"] = 0;
				retrieve_Zdecay_topo(ptcl->child(ichild));
				break;
			}
			//std::cout<<"Z child status "<<ptcl->child(ichild)->status()<<" "<<ptcl->child(ichild)->pdgId()<<" "<<ptcl->child(ichild)->barcode()<<std::endl;
			if(_child_pid>0) {
			  p1.SetPtEtaPhiM(ptcl->child(ichild)->pt(), ptcl->child(ichild)->eta(), ptcl->child(ichild)->phi(), ptcl->child(ichild)->p4().M());
			  l1_truth_pdgId = _child_pid;
			}
			else {
			  p2.SetPtEtaPhiM(ptcl->child(ichild)->pt(), ptcl->child(ichild)->eta(), ptcl->child(ichild)->phi(), ptcl->child(ichild)->p4().M());
			  l2_truth_pdgId = _child_pid;
			}
			if(_child_pid>0) _child_pid = _child_pid*100;
			(*m_ieventmap)["mc_Z_decay_topo"] += fabs(_child_pid);

		}
	}

      void TruthSelect::retrieve_truthjets(TLorentzVector p1, TLorentzVector p2, TLorentzVector ph){
          doDphi = false;
          HG::TruthJets selected_truthJets(SG::VIEW_ELEMENTS);
          HG::TruthJets central_truthJets(SG::VIEW_ELEMENTS);
          for(const xAOD::Jet *jet : *m_antiKt4TruthJets){
            if(jet->pt()/1000 < 10||fabs(jet->eta())>4.4) continue;
            if(fabs(jet->p4().DeltaR(p1))<0.2||fabs(jet->p4().DeltaR(p2))<0.2||fabs(jet->p4().DeltaR(ph))<0.2) continue;
            selected_truthJets.push_back(jet);
            if(fabs(jet->eta())<2.5) central_truthJets.push_back(jet);
          }

          N_jets = selected_truthJets.size();
          Central_N_jets = central_truthJets.size();

          selected_truthJets.sort(HG::JetHandler::comparePt);

          if(selected_truthJets.size()>0) {
            j1_truth_pt = selected_truthJets[0]->pt()/1000;
	    j1 = selected_truthJets[0]->p4();
          }

            if(selected_truthJets.size()>1) {
            j2_truth_pt = selected_truthJets[1]->pt()/1000;
            j2 = selected_truthJets[1]->p4(); doDphi = true;
          }
        }

	void TruthSelect::truthsave_vertex(){

		const xAOD::TruthVertex* pVertex = 0;
		double PVz = 0;
		double PVx = 0;
		double PVy = 0;
		double PVt = 0;

		if( m_Vertices->size() > 0 ){
			pVertex = (*m_Vertices)[0];
			PVz = pVertex->z();
			PVx = pVertex->x();
			PVy = pVertex->y();
			PVt = pVertex->t();
		}

		(*m_feventmap)["mc_PVz"] = PVz;
		(*m_feventmap)["mc_PVx"] = PVx;
		(*m_feventmap)["mc_PVy"] = PVy;
		(*m_feventmap)["mc_PVt"] = PVt;

	}

	HG::TruthPtcls TruthSelect::findHiggs() {
		HG::TruthPtcls ys(SG::VIEW_ELEMENTS);
		for (auto ptcl : *m_truthParticles)
			if(ptcl->pdgId() == 25 && (ptcl->status()==62||ptcl->status()==11)) {ys.push_back(ptcl); }
		//        if(ptcl->pdgId() == 25 && (ptcl->status()==62)) {ys.push_back(ptcl); }
		//if ( isGoodTruthPhoton(ptcl) && isFromHiggs(ptcl) ) ys.push_back(ptcl);
		return ys;
	}

	void TruthSelect::settree_branch(){
	  m_truthTree->Branch("m_eventNumber", &m_eventNumber);
	  m_truthTree->Branch("m_mc_channel_number", &m_mc_channel_number);
	  m_truthTree->Branch("higgs_truth_pt",&higgs_truth_pt);
	  m_truthTree->Branch("higgs_truth_eta",&higgs_truth_eta);
	  m_truthTree->Branch("higgs_truth_phi",&higgs_truth_phi);
	  m_truthTree->Branch("higgs_truth_mass",&higgs_truth_mass); 
	  m_truthTree->Branch("l1_truth_pt",&l1_truth_pt); 
	  m_truthTree->Branch("l1_truth_eta",&l1_truth_eta); 
	  m_truthTree->Branch("l1_truth_phi",&l1_truth_phi); 
	  m_truthTree->Branch("l1_truth_mass",&l1_truth_mass); 
	  m_truthTree->Branch("l1_truth_pdgId",&l1_truth_pdgId); 
	  m_truthTree->Branch("l2_truth_pt",&l2_truth_pt); 
	  m_truthTree->Branch("l2_truth_eta",&l2_truth_eta); 
	  m_truthTree->Branch("l2_truth_phi",&l2_truth_phi); 
	  m_truthTree->Branch("l2_truth_mass",&l2_truth_mass); 
	  m_truthTree->Branch("l2_truth_pdgId",&l2_truth_pdgId); 
	  m_truthTree->Branch("ph_truth_pt",&ph_truth_pt); 
	  m_truthTree->Branch("ph_truth_eta",&ph_truth_eta); 
	  m_truthTree->Branch("ph_truth_phi",&ph_truth_phi); 
	  m_truthTree->Branch("ph_truth_mass",&ph_truth_mass);
	  m_truthTree->Branch("Z_truth_pt",&Z_truth_pt); 
	  m_truthTree->Branch("Z_truth_eta",&Z_truth_eta);
	  m_truthTree->Branch("Z_truth_phi",&Z_truth_phi);
	  m_truthTree->Branch("Z_truth_mass",&Z_truth_mass);
	  m_truthTree->Branch("Dphi_Zy_jj", &Dphi_Zy_jj);
	  m_truthTree->Branch("pTt", &pTt);
	  m_truthTree->Branch("dphi_Zy", &dphi_Zy);
	  m_truthTree->Branch("eta_Zepp", &eta_Zepp);
	  m_truthTree->Branch("m_jj", &m_jj);
	  m_truthTree->Branch("DRmin_Z_j", &DRmin_Z_j);
	  m_truthTree->Branch("Dy_j_j", &Dy_j_j);
	  m_truthTree->Branch("BDTG", &mvaout);
	  m_truthTree->Branch("N_jets", &N_jets);
	  m_truthTree->Branch("Central_N_jets", &Central_N_jets);
	  m_truthTree->Branch("j1_truth_pt",&j1_truth_pt);
	  m_truthTree->Branch("j2_truth_pt",&j2_truth_pt);
	}

	void TruthSelect::Reco_Truth_Matching(Object_llg &_llg)
	{
		if(_llg.index_photon==-1) {(*m_beventmap)["ph_istruth"]=0;  (*m_feventmap)["DR_ph_istruth"]=999;    (*m_beventmap)["ph_good_truth"]=false;}
		else { (*m_beventmap)["ph_istruth"]= isTruthPhoton((*_llg.m_selected_photons)[_llg.index_photon]);
			(*m_feventmap)["DR_ph_istruth"]= p_ph.DeltaR((_llg.TLV_photon)) ;

            const xAOD::Photon *ph=(*_llg.m_selected_photons)[_llg.index_photon];
            const xAOD::TruthParticle* truthph = xAOD::TruthHelpers::getTruthParticle(*ph);
            if(!truthph)    {
                (*m_ieventmap)["ph_truth_pdgID"]=0;
                (*m_beventmap)["ph_good_truth"]=false;
            }
            else    {
                (*m_ieventmap)["ph_truth_pdgID"]=truthph->pdgId();
                //cout<<"barcode "<<truthph->barcode();
                //getchar();
                (*m_beventmap)["ph_good_truth"]=HG::isGoodTruthPhoton(truthph);
            }
		}
		if(_llg.m_channel==1){
			(*m_beventmap)["l1_istruth"] = isTruthElectron((*_llg.m_selected_electrons)[_llg.index_lepton1]); 
			(*m_beventmap)["l2_istruth"] = isTruthElectron((*_llg.m_selected_electrons)[_llg.index_lepton2]); 
		} else{
			(*m_beventmap)["l1_istruth"] = isTruthMuon((*_llg.m_selected_muons)[_llg.index_lepton1]);
			(*m_beventmap)["l2_istruth"] = isTruthMuon((*_llg.m_selected_muons)[_llg.index_lepton2]);
		}
	}

	const xAOD::TruthVertex* TruthSelect::getMotherVert(const xAOD::TruthParticle* p)
	 {
	  const xAOD::TruthVertex* v = p->prodVtx();
	  while (v && v->nIncomingParticles() == 1 &&
	      v->incomingParticle(0)->pdgId() == p->pdgId())
	   {
	    p = v->incomingParticle(0);
	    v = p->prodVtx();
	   }
	  return v;
	 }
	const xAOD::TruthParticle* TruthSelect::getMother(const xAOD::TruthParticle* p)
	 {
	  const xAOD::TruthVertex* v = getMotherVert (p);
	  if (!v || v->nIncomingParticles() == 0)
	    return 0;
	  return v->incomingParticle(0);
	 }

	int TruthSelect::AncestorID(const xAOD::TruthParticle *ptcl){
	  int ID = ptcl->pdgId();
	  const xAOD::TruthParticle *father;
	  const xAOD::TruthParticle *grandfather;
	  if (ptcl->nParents()!=0) { father = ptcl->parent(0); ID = father->pdgId();} else return ID;
	  if (father->nParents()!=0) { grandfather = father->parent(0); ID = grandfather->pdgId();} else return ID;
	  while (grandfather->nParents()!=0){
	    father = grandfather->parent(0);
	    grandfather = father ;
	    if (ID==25||ID==23) break;
	    ID = grandfather->pdgId();
	  }
	  return ID;
	}


}


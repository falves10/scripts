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

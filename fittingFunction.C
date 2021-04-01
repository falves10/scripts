///#include "atlas/CommonHead.h"
#include <TTreeFormula.h>

#include "atlas/AtlasStyle.C"
#include "atlas/AtlasStyle.h"
#include "atlas/AtlasUtils.C"
#include "atlas/AtlasUtils.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TLorentzVector.h"

#include "TMath.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "TChain.h"
#include "TROOT.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TSystem.h"

#include <math.h> // for “fabs”

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <vector>
//RooFit include
#include "RooDataHist.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include <stdio.h>
#include <sys/stat.h>

#include "TMinuit.h"
using namespace RooFit ;
using namespace std;


//int m_nEtaBin(13);
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.37, 1.52,
//     1.8, 1.9, 2.1, 2.3, 2.47});

int m_nEtaBin(16);
//int m_nEtaBin(14);
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});
std::vector<double> m_etaBinning({0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30,1.35,1.40});
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.775, 0.825, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});
//std::vector<double> m_etaBinning({0.0, 0.2, 0.4, 0.6, 0.7, 0.825, 1., 1.2, 1.37, 1.52,1.8, 1.9, 2.1, 2.3, 2.47});

//EL:
//int m_nEtaBin(2);
//std::vector<double> m_etaBinning({0.0, 1.7});

//N:
//int m_nEtaBin(2);
//std::vector<double> m_etaBinning({0.0, 1.5});

//IBL:
//int m_nEtaBin(3);
//std::vector<double> m_etaBinning({0.0, 1.0, 2.0});

//FMX:
//int m_nEtaBin(1);
//std::vector<double> m_etaBinning({0.0, 1.4});

//PP0
//int m_nEtaBin(3);
//std::vector<double> m_etaBinning({0.0, 1.0, 2.2});

//A:
//int m_nEtaBin(5);
//std::vector<double> m_etaBinning({0.0, 0.6, 0.8, 1.3, 1.7});

TString type = "nominal";
TString name = "data";


//int m_nEnergyBin(10);
int m_nEnergyBin(10);
//std::vector<double> m_energyBinning({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0});
std::vector<double> m_energyBinning({0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0});


bool calibrated(true);//Adele true
float lowcut = 87.;
float highcut = 94.5;

void PrintProgressBar( int index, int total )
{
  if( index%200000 == 0 )
    {
      TString print_bar = " [";
      for( int bar = 0; bar < 20; bar++ )
	{
	  double current_fraction = double(bar) / 20.0;
	  if( double(index)/double(total) > current_fraction )
	    print_bar.Append("/");
	  else
	    print_bar.Append(".");
	}
      print_bar.Append("] ");
      std::cout << print_bar << 100.*(double(index)/double(total)) << "%\r" << std::flush;
    }
}





void infile2chain(TString _infilelist, TChain *&_fchain, TString chainname)
{
  _fchain = new TChain(chainname);
  ifstream infile(_infilelist, ios::in);
  string line;
  while (getline(infile, line)){
    _fchain->Add(line.c_str());
  }
  infile.close();
}

//convolution function for the mass
std::vector<float> convolution(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, double input_width,double input_sigma, TString path, int m_layer);

std::vector<float> convolution_eOverp(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, TString path, int m_layer, TString campaign, TString type, TString name);

//methods to save plots
void mySave(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, int iEnergyBin,TString vartagname, TString varname,TString path, int m_layer);
void mySave2D(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, TString vartagname, TString varname, TString path,TString parameter, int m_layer);

/// Get the bin number of a given eta value
unsigned int binNbOfEta(double eta) 
{
  double localEta = eta;
  localEta = fabs(eta);
  unsigned int bin = 0;
  while (bin < m_nEtaBin && localEta > m_etaBinning[bin+1])
    bin++;
  return bin;
} // end binNbOfEta


/// Get the bin number of a given E1/E2 value
unsigned int binNbOfEnergy(double E1E2) 
{
  
  double localEnergy = E1E2;
  localEnergy = fabs(E1E2);
  unsigned int bin = 0;
  while (bin < m_nEnergyBin && localEnergy > m_energyBinning[bin+1])
    bin++;
  return bin;
} // end binNbOfE1E2

vector<double>params(int iBias, int iEtaBin, int iEnergyBin, TString type, TString name);



//MAIN
//int invariantMassFitCal2(TString input) {
int fittingFunction(TString campaign,TString layer, TString variable) {
 
  SetAtlasStyle();  
   
  vector<TString> varname;
  varname.push_back(variable);
  
  vector<TString> vartagname;
  if(variable.Contains("mass"))
    vartagname.push_back("m_{ee} [GeV]");
  else if (variable.Contains("eOverp"))
    vartagname.push_back("E/p");

  int m_layer = atoi(layer);
  std::cout<<"m_layer = "<<m_layer<<std::endl;
    
  for(int ivar = 0; ivar < varname.size(); ivar++){

    std::cout<<"variable = "<<varname[ivar]<<std::endl;

    //TString path = Form("/publicfs/atlas/atlasnew/higgs/hgg/donofrioadele/EG_condor/EGAMMAStep2Code/output_"+varname[ivar]+"_"+campaign+"_layer%i_%f_%f/",m_layer,lowcut,highcut);
      
    TString path = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/flaOut_dataNominal_twoLeptons/";
    gSystem->mkdir(path, kTRUE);

    std::cout<<"path_out = "<<path<<std::endl;

    TH1F* map_histog[17];
    TH1F* vhisto2D_data[17];
    TH1F* vhisto2DMean_data[17];
    TH1F* vhisto2DSigma_data[17];
    
    std::map<pair<int, int>, map<int,TH1F*>> m_histo; 
    
    int biasSize = 16;
    //int biasSize = 2;
    
    int binning; float lowedge; float highedge; float lowcut; float highcut;// mass
    if(varname[ivar].Contains("mass")) {binning = 50; lowedge = 80; highedge = 100; lowcut = 80; highcut = 100;}// mass: 80-100
    if(varname[ivar].Contains("m_eOverp")) {binning = 50; lowedge = 0.6; highedge = 2; lowcut = 0.6; highcut = 2;}// E/p
    
    unsigned int maxEta = 0;
    
    maxEta = m_nEtaBin;
    //get input files
    //TString path_in = Form("/publicfs/atlas/atlasnew/higgs/hgg/donofrioadele/EG_condor/EGAMMAStep1Code/adele_ntuple/calibrated_"+varname[ivar]+"/Histograms/layer%i/",m_layer);
      
    TString path_in = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";

    
     
      
      
    float array[11] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.,2.0};
    
    TH1F* histo2D_mc = new TH1F ("mass_etabin_mc", "", m_nEnergyBin,array);
    TH1F* histo2D_data = new TH1F ("mass_etabin_data", "", m_nEnergyBin,array);
    //Mean mass method
    TH1F* histo2DMean_mc = new TH1F ("mass_etabin_mcMean", "", m_nEnergyBin,array);
    TH1F* histo2DMean_data = new TH1F ("mass_etabin_dataMean", "", m_nEnergyBin,array);
    //CB Parameters values
    TH1F* histo2DSigma_data = new TH1F ("mass_etabin_dataSigma", "", m_nEnergyBin,array);
    TH1F* histo2DSigma_mc = new TH1F ("mass_etabin_mcSigma", "", m_nEnergyBin,array);    
    
    histo2D_mc->Reset();
    histo2D_data->Reset();
    histo2DMean_mc->Reset();
    histo2DMean_data->Reset();
    histo2DSigma_mc->Reset();
    histo2DSigma_data->Reset();
    
    std::cout<<"biasSize = "<<biasSize<<std::endl;

    //weight_histo_eta7_energy9_bias3_m_eOverp_mergedOutput.root   
      
    cout<<"type"<<type<<"name"<<name<<endl;  
    TString forMaterial = Form("_%s_%s", type.Data(), name.Data());
      
    TFile *f_in;
      
    if( type.Contains("nominal") && name.Contains("data") )
    {
        f_in = new TFile(path_in+"weight_histo_"+campaign+"_"+varname[ivar]+".root");
        if ( !(f_in->IsOpen()) )continue;  
    }
    else
    {
        f_in = new TFile(path_in+"weight_histo_"+campaign+"_"+varname[ivar]+forMaterial+".root");
        if ( !(f_in->IsOpen()) )continue;
        std::cout<<"path_in = "<<path_in+"weight_histo_"+campaign+"_"+varname[ivar]+forMaterial+".root"<<std::endl;   
    }
    
      
    for (unsigned int iEtaBin = 0; iEtaBin < maxEta; iEtaBin++) {
        
      //if( iEtaBin != 2 ) continue;
      cout<<" Fabio_lep_eta " << iEtaBin << endl;
      for(int i = 2; i < biasSize+2 ; i++){//Adele
          
	   vhisto2D_data[i-2] = new TH1F (Form("mass_etabin%d_data_bias%i",iEtaBin,i), "", m_nEnergyBin,array);
	   vhisto2DMean_data[i-2] = new TH1F (Form("massMean_etabin%d_data_bias%i",iEtaBin,i), "", m_nEnergyBin,array);
	   vhisto2DSigma_data[i-2] = new TH1F (Form("mass_etabin%d_data_biasSigma%i",iEtaBin,i), "", m_nEnergyBin,array);
          
          
      }
      
      for (unsigned int iEnergyBin = 0; iEnergyBin < m_nEnergyBin; iEnergyBin++) {
          
        //if(iEnergyBin != 3) continue;
          
        map_histog[1] = new TH1F(Form("data_%d_%d",iEtaBin,iEnergyBin), "",    binning,  lowedge, highedge);
        map_histog[0] = new TH1F(Form("mc_%d_%d",iEtaBin,iEnergyBin), "", binning,  lowedge, highedge);

        map_histog[0] = (TH1F*) f_in->Get(Form("mc_eta_%d_energy_%d",iEtaBin,iEnergyBin));
        map_histog[1] = (TH1F*) f_in->Get(Form("data_eta_%d_energy_%d",iEtaBin,iEnergyBin));

        if(map_histog[0]->Integral()==0) continue;
        if(map_histog[1]->Integral()==0) continue;
	
	   for(int i = 0; i < biasSize ; i++)
        {//Adele
	       map_histog[i+2] = (TH1F*) f_in->Get(Form("data_eta_%d_energy_%d_bias%d",iEtaBin,iEnergyBin,i));
	       map_histog[i+2]->Sumw2();
	   }
	  
	  histo2D_mc->Sumw2();
	  histo2D_data->Sumw2();
	  
	  ///////////////////////////////////////////////////////////////////////////////
	  float mean_mass_mc = -1.;
	  float mean_mass_data = -1.;
	  float emean_mass_mc = -1.;
	  float emean_mass_data = -1.;
	  float width_mass_mc = 2.49;
	  float width_mass_data = -1.;
	  float ewidth_mass_mc = -1.;
	  float ewidth_mass_data = -1.;
	  float sigma_mass_mc = 2.6;
	  float sigma_mass_data = -1.;
	  float esigma_mass_mc = -1.;
	  float esigma_mass_data = -1.;
	  
	  std::vector<float> result_data;
	  std::vector<float> result_mc;

	  if(varname[ivar].Contains("mass"))
	    result_mc = convolution(map_histog[0],-1,iEtaBin,iEnergyBin,1,map_histog[0]->GetMean(),width_mass_mc,sigma_mass_mc,path, m_layer);
	  else if(varname[ivar].Contains("eOverp"))
      {
          if(map_histog[0]->GetEntries() <= 1000) break; //for new binning, data/MC
          //if(map_histog[0]->GetEntries() <= 00) break; 
          else result_mc = convolution_eOverp(map_histog[0],-1,iEtaBin,iEnergyBin,1,map_histog[0]->GetMean(),path, m_layer, campaign, type, name);
      }
           
	
	  mean_mass_mc = result_mc.at(0);
	  emean_mass_mc = result_mc.at(1);
          
	  if(varname[ivar].Contains("mass"))
      {
	    sigma_mass_mc = result_mc.at(3);
	    esigma_mass_mc = result_mc.at(4);
	  }
          
	  if(varname[ivar].Contains("mass"))
	    result_data = convolution(map_histog[1],-1,iEtaBin,iEnergyBin,0,mean_mass_mc,width_mass_mc,sigma_mass_mc,path, m_layer);
	  else if(varname[ivar].Contains("eOverp"))
      {
          if(map_histog[1]->GetEntries() <= 1000 ) break;////for new binning, data/MC
          //if(map_histog[1]->GetEntries() <= 100) break;
          else result_data = convolution_eOverp(map_histog[1],-1,iEtaBin,iEnergyBin,0,mean_mass_mc,path, m_layer, campaign, type, name);
      }
	      
	  
          
      mean_mass_data = result_data.at(0);
	  emean_mass_data = result_data.at(1);
	  
      if(varname[ivar].Contains("mass"))
      {
	    sigma_mass_data = result_data.at(3);
	    esigma_mass_data = result_data.at(4);
	  }
          
	  histo2D_mc->SetBinContent(iEnergyBin+1,mean_mass_mc);
	  histo2D_data->SetBinContent(iEnergyBin+1,mean_mass_data);
	  histo2D_mc->SetBinError(iEnergyBin+1,emean_mass_mc);
	  histo2D_data->SetBinError(iEnergyBin+1,emean_mass_data);

	  if(varname[ivar].Contains("mass")){
	    //Mean mass method
	    histo2DMean_mc->SetBinContent(iEnergyBin+1,map_histog[0]->GetMean());
	    histo2DMean_mc->SetBinError(iEnergyBin+1,map_histog[0]->GetMeanError());
	    histo2DMean_data->SetBinContent(iEnergyBin+1,map_histog[1]->GetMean());
	    histo2DMean_data->SetBinError(iEnergyBin+1,map_histog[1]->GetMeanError());
	    //Plot all the CB+BW function parameters
	    histo2DSigma_data->SetBinContent(iEnergyBin+1,result_data.at(3));
	    histo2DSigma_data->SetBinError(iEnergyBin+1,result_data.at(4));
	    histo2DSigma_mc->SetBinContent(iEnergyBin+1,result_mc.at(3));
	    histo2DSigma_mc->SetBinError(iEnergyBin+1,result_mc.at(4));
	  }

	  for(int i = 2; i < biasSize+2 ; i++)
      {//Adele
          
        
          
	    float vmean_mass_data = -1;
	    float vemean_mass_data = -1;
	    
	    std::vector<float> vresult;
	    vresult.clear();
	    int binmax = map_histog[i]->GetMaximumBin(); 
	    double x = map_histog[i]->GetXaxis()->GetBinCenter(binmax);
        //map_histog[i]->Print("all");
	    if(varname[ivar].Contains("mass"))
	      vresult = convolution(map_histog[i],i,iEtaBin,iEnergyBin,0,mean_mass_mc,width_mass_mc,sigma_mass_mc,path, m_layer);
          
	    if(varname[ivar].Contains("eOverp"))
            //if(iEtaBin == 0 && i == 2)
            {
                if(map_histog[i]->GetEntries() <= 1000) break;//for new binning, data/MC
                //if(map_histog[i]->GetEntries() < 100) break;
                else vresult = convolution_eOverp(map_histog[i],i,iEtaBin,iEnergyBin,0,mean_mass_mc,path, m_layer, campaign, type, name);  
                
                vmean_mass_data = vresult.at(0);
	             vemean_mass_data = vresult.at(1);
            }
	        
	    //vmean_mass_data = vresult.at(0);
	    //vemean_mass_data = vresult.at(1);
        cout<< vmean_mass_data << " " << vemean_mass_data << endl;
        
	    vhisto2D_data[i-2]->SetBinContent(iEnergyBin+1,vresult.at(0));
	    vhisto2D_data[i-2]->SetBinError(iEnergyBin+1,vresult.at(1));
        vhisto2D_data[i-2]->Print("all");
	    if(varname[ivar].Contains("mass")){
	      //Mean mass method
	      vhisto2DMean_data[i-1]->SetBinContent(iEnergyBin+1,map_histog[i]->GetMean());
	      vhisto2DMean_data[i-1]->SetBinError(iEnergyBin+1,map_histog[i]->GetMeanError());
	      //Plot all the CB+BW function parameters
	      vhisto2DSigma_data[i-1]->SetBinContent(iEnergyBin+1,vresult.at(3));
	      vhisto2DSigma_data[i-1]->SetBinError(iEnergyBin+1,vresult.at(4));
	    }
	  }
          
	  mySave(map_histog[0],map_histog[1],-1,iEtaBin,iEnergyBin,vartagname[ivar],varname[ivar],path, m_layer);
	  for(int i = 0; i < biasSize ; i++)
      {//Adele
	    mySave(map_histog[0],map_histog[i+2],i+2,iEtaBin,iEnergyBin,vartagname[ivar],varname[ivar],path, m_layer);
	  }
          
	}//end loop energy bins
      
    
	mySave2D(histo2D_mc,histo2D_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"", m_layer);
	if(varname[ivar].Contains("mass")){
	  mySave2D(histo2DMean_mc,histo2DMean_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"Mean", m_layer);
	  mySave2D(histo2DSigma_mc,histo2DSigma_data,-1,iEtaBin,vartagname[ivar],varname[ivar],path,"Sigma", m_layer);
	}
    
	for(int i = 0; i < biasSize ; i++)
    {//Adele
      
	  mySave2D(histo2D_mc,vhisto2D_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"", m_layer);
	    if(varname[ivar].Contains("mass")){
	      mySave2D(histo2DMean_mc,vhisto2DMean_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"Mean", m_layer);
	      mySave2D(histo2DSigma_mc,vhisto2DSigma_data[i],i+2,iEtaBin,vartagname[ivar],varname[ivar],path,"Sigma", m_layer);
	    }
	}
	//saving 2D histos
    }//end loop eta bins
      
    delete histo2D_mc;
    delete histo2D_data;
    delete histo2DSigma_mc;
    delete histo2DSigma_data;
    delete histo2DMean_mc;
    delete histo2DSigma_mc;
  }//end loop variables
  return 0;
}



std::vector<float> convolution(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, double input_width, double input_sigma, TString path, int m_layer){	

  std::vector<float> parameters;
  parameters.clear();

  //change the range of the plots to get the mean!!!
  RooRealVar x( "x", "x", 87., 94.5 );//84,98//80-100
  x.setBins(10000,"cache") ;
  x.setMin("cache",64.) ;
  x.setMax("cache",118.) ;
  
  //std::cout<<"mean_bw = "<<mean_bw<<std::endl;
  // Breit-Wigner                        
  RooRealVar m0( "m0", "m0", mean_bw, 87., 94.5);//80-100
  RooRealVar width( "width", "width", 2.49);//2.49,1.,4.
  if(iperiod == 0){
      width.setConstant(kTRUE);
  }
  RooBreitWigner bw( "bw", "bw", x, m0, width);

  // Crystal-Ball                                                                                                                         
  RooRealVar mean( "mean", "mean", 0. );
  RooRealVar sigma( "sigma", "sigma", input_sigma, 1., 5.);//2.6,1.,5.
  if(iperiod == 0){
    sigma.setConstant(kTRUE);
  }
  RooRealVar alpha( "alpha", "alpha", 1.3 );
  RooRealVar n( "n", "n", 5.1 );
  RooCBShape cb( "cb", "cb", x, mean, sigma, alpha, n );
  // convolution                                                                                                                         
  RooFFTConvPdf pdf_sig( "pdf_sig", "pdf_sig", x, bw, cb );
  //add polynomial bkg
  // Build Chebychev polynomial p.d.f.
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",0.,0.,1.) ;
  RooChebychev bkg1("bkg1","bkg1",x,RooArgSet(a0,a1));
  // Sum the composite signal and background
  RooRealVar bkgfrac("bkgfrac","fraction of background",0.05,0.,1.);
  RooAddPdf  pdf("pdf","pdf",RooArgList(bkg1,pdf_sig),bkgfrac) ;

  RooDataHist histo("histo","histo",x,Import(*map_histog));

  if(iperiod == 0)
    pdf.fitTo(histo,SumW2Error(kTRUE)) ;
  if(iperiod == 1)
    pdf_sig.fitTo(histo,SumW2Error(kTRUE)) ;

  TCanvas canv( "canv", "canv", 800., 600. );
  RooPlot* frame1 = x.frame(Bins(100),Title("Convolution of a Breit-Wigner and a Crystal-Ball, Chebychev pol. bkg")) ;
  histo.plotOn(frame1,Name("Data")) ;
  if(iperiod == 0){
    pdf.plotOn(frame1,Name("pdf"),LineColor(kRed)) ;
    pdf.paramOn(frame1,Layout(0.60));
    pdf.plotOn(frame1,Components("bkg1"),LineStyle(kDotted),LineColor(kBlue));
  }
  if(iperiod == 1){
    pdf_sig.plotOn(frame1,Name("pdf"),LineColor(kRed)) ;
    pdf_sig.paramOn(frame1,Layout(0.60));
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
  canvas->cd() ; 
  TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
  TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
  pad1->Draw();pad2->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.16);
  pad2->SetBottomMargin(0.24);
  frame1->GetYaxis()->SetTitleOffset(1.4) ;
  frame1->GetXaxis()->SetTitle((Form("m_{ee} [GeV], #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i ", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));
  frame1-> Draw();
  pad1->Modified();
  pad1->RedrawAxis();
  pad1->Update();
  pad2->cd();

  TF1* tf1_model_data = pdf.asTF(x);
  TF1* tf1_model = pdf_sig.asTF(x);

  TH1* clone_data = (TH1*)histo.createHistogram("clone_data",x,Binning(50));
  RooDataHist *pdfHisto_data = pdf.generateBinned(x,1000000);
  TH1* clone_fit_data = (TH1*)pdfHisto_data->createHistogram("clone_fit_data",x,Binning(50));
  clone_fit_data->Scale(clone_data->Integral()/clone_fit_data->Integral());

  RooDataHist *pdfHisto = pdf_sig.generateBinned(x,1000000);
  TH1* clone_fit = (TH1*)pdfHisto->createHistogram("clone_fit",x,Binning(50));
  clone_fit->Scale(clone_data->Integral()/clone_fit->Integral());

  if(iperiod == 0)
    clone_data->Divide(clone_fit_data);

  else if(iperiod == 1)
    clone_data->Divide(clone_fit);

  clone_data->GetXaxis()->SetTitle((Form("m_{ee} [GeV], #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));
  clone_data->GetYaxis()->SetTitle("DATA / FIT");
  clone_data->GetXaxis()->SetRangeUser(lowcut, highcut);//80-100
  clone_data->GetYaxis()->SetRangeUser(0.7,1.3);
  clone_data->GetXaxis()->SetLabelSize(0.1);
  clone_data->GetYaxis()->SetLabelSize(0.08);
  clone_data->GetXaxis()->SetTitleSize(0.08);
  clone_data->GetYaxis()->SetTitleSize(0.09);
  clone_data->GetYaxis()->SetTitleOffset(0.6);
  clone_data->GetXaxis()->SetTitleOffset(1.2);
  clone_data->Draw("E1");
  pad2->Modified();
  pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();

  TString output_folder = path+Form("FitPlots/layer%i/",m_layer);
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } 
  else {
    gSystem->mkdir(path+Form("FitPlots/layer%i/",m_layer), kTRUE);
  }

  canvas->SaveAs(path+Form("FitPlots/layer%i/weight_period%i_invariantMass_eta%i_energy%i_bias%i.pdf",m_layer,iperiod,iEtaBin,iEnergyBin,iBias));
  canvas->SaveAs(path+Form("FitPlots/layer%i/weight_period%i_invariantMass_eta%i_energy%i_bias%i.root",m_layer,iperiod,iEtaBin,iEnergyBin,iBias));

   
  double modelMean_data = tf1_model_data->GetMaximumX();
  double modelMean = tf1_model->GetMaximumX();

  if(iperiod == 0)
    parameters.push_back(modelMean_data);//0
  if(iperiod == 1)
    parameters.push_back(modelMean);//0
  parameters.push_back(m0.getError());//1
  parameters.push_back(m0.getVal());//2
  parameters.push_back(sigma.getVal());//3
  parameters.push_back(sigma.getError());//4
  parameters.push_back(mean.getVal());//5
  parameters.push_back(mean.getError());//6
  parameters.push_back(alpha.getVal());//7
  parameters.push_back(alpha.getError());//8
  parameters.push_back(n.getVal());//9
  parameters.push_back(n.getError());//10
  parameters.push_back(width.getVal());//11
  parameters.push_back(width.getError());//12

  cout << "width = "<<width.getVal() << endl;

  delete frame1;
  delete clone_data;
  delete clone_fit;
  delete pad1;
  delete pad2;
  delete canvas;

  return parameters;

}


//////////////////////////////////////////eOverP////////////////////////////////////////
 std::vector<float> convolution_eOverp(TH1F* map_histog, int iBias, int iEtaBin, int iEnergyBin, int iperiod, double mean_bw, TString path, int m_layer, TString campaign, TString type, TString name)
 {
cout<<" ETA " << iEtaBin << " ENERGY " << iEnergyBin << " BIAS "<< iBias << endl;

std::vector<TString> parname_ws    = {"Constant", "Mean","Sigma","Alpha", "N"};
std::vector<float> parameters;
parameters.clear();
    
//TF1 *f2 = new TF1("f2","crystalball",0.9,1.3); //Nominal fit range 
//TF1 *f2 = new TF1("f2","crystalball",0.95,1.5); //Alternative fit range 
     
TF1 *f2 = new TF1("f2","crystalball",0.9,1.3);


f2->SetLineColor(2);
bool teste = false;
bool noParms = false; //for other fitrange
bool isNomFitRange = true;
bool isPP0 = false;
bool isNomElectronID = true; //used false for etabin4 data 2017 agains pp0 config
 
  
     
if( iperiod == 0 )
{
    cout<<" BEGIN DATA "<<endl;
    
    
    if ( iBias == -1 )
    {
        /*Material 2 leptons new binning
        if(type.Contains("data17") && name.Contains("A"))
		  f2->SetParameters(1e3, 0.99 ,0.07, -0.8, 3);		
        else if(iEtaBin == 2 || iEtaBin == 3) // all except data17A
	    {
            vector<double> par = params( 0,  iEtaBin,  0, type, name);
            f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
        }
        else
        {
            
            vector<double> par = params( 0,  iEtaBin,  iEnergyBin, type, name);
            f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
            
        }
        */
        vector<double> par = params( 0,  iEtaBin,  iEnergyBin, type, name);
        f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
        
    } 
    else
    {
        if(type.Contains("data17") && name.Contains("N"))
        {
            vector<double> par = params( 0,  iEtaBin,  0, type, name);
            f2->SetParameters( par[0], par[1], par[2], par[3], par[4]); 
        }
        else if (type.Contains("data17") && name.Contains("A"))
        {
		if( iEtaBin == 2 || iEtaBin == 3)
		{
			//if(iEnergyBin <= 7)
			f2->SetParameters(1e3, 0.99 ,0.07, -0.8, 3);
			///else f2->SetParameters(1e3, 1.0 ,0.99, -2.4046, 3);
		}
		else
		{
			vector<double> par = params( 0,  iEtaBin,  0, type, name);
                    	f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
		}
	}
        else
        {
            //Material 2 leptons new binning
            /*
            if(iEtaBin == 3)
            {
                if( iEnergyBin == 8)
                {
                    cout<<"enterhere"<<endl;   
                    f2->SetParameters( 5.03505e+04 , 1.01248e+00, 1.13589e-01 , -5.46450e-01, 6.10820e+00);
                }
                else
                {
                    vector<double> par = params( 0,  iEtaBin,  0, type, name);
                    f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);   
                }
            
            }
            else
            {
                vector<double> par = params( iBias-2,  iEtaBin,  iEnergyBin, type, name);
                f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
            }
            */
            
            vector<double> par = params( iBias-2,  iEtaBin,  iEnergyBin, type, name);
            f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
        }  
    } 
    
    
    //data nominal and material
    f2->SetParLimits(0,0.,1e10);
    f2->SetParLimits(1,0.8,1.2);
    f2->SetParLimits(2,0.01,0.5);
    f2->SetParLimits(3,-10.,10.);
    f2->SetParLimits(4,1.,10);

    cout<<" END DATA "<<endl;   
}
else if(iperiod == 1)
{
    cout<<" BEGIN MC "<<endl;
    /*if( (iEtaBin == 0 && iEnergyBin == 4) || (iEtaBin == 0 && iEnergyBin == 7)  )  
        f2->SetParameters(1, map_histog->GetMean()  , 0.01 , -0.5 , 2 );
    else if ( (iEtaBin == 1 && iEnergyBin == 9) )
        f2->SetParameters(1, map_histog->GetMean()  , 0.01 , -0.2 , 2 );
    else
        f2->SetParameters(1, map_histog->GetMean()  , 0.01 , -1 , 2 );
    */
    
    /*
    if( (iEtaBin == 0 && iEnergyBin == 4) || (iEtaBin == 0 && iEnergyBin == 7) )  
        f2->SetParameters(1, 1.  , 0.1 , -0.5 , 2 );
    else if ( (iEtaBin == 1 && iEnergyBin == 9) )
        f2->SetParameters(1, 1.  , 0.1 , -0.2 , 2 );
    else
        f2->SetParameters(1, 1.  , 0.08 , -1 , 2 );
    */
    
    if( iEnergyBin == 0 )
        f2->SetParameters( 24.491, 0.996921, 0.041169, -1.02395, 2.79162 );
    else if( iEnergyBin == 1 )
        f2->SetParameters( 71.8089, 1.00044, 0.0411136, -0.893548, 2.61719 );
    else if( iEnergyBin == 2 )
        f2->SetParameters( 116.533, 1.00307, 0.041934, -0.800426, 2.39107 );
    else if( iEnergyBin == 3 )
        f2->SetParameters( 119.281, 1.00593, 0.0434437, -0.727336, 2.23322 );
    else if( iEnergyBin== 4 )
        f2->SetParameters( 88.0382, 1.0081, 0.0452797, -0.670262, 2.13681);
    else if( iEnergyBin == 5 )
        f2->SetParameters( 50.8651, 1.01055, 0.0474722, -0.641028, 1.95129);
    else if( iEnergyBin == 6 )
    {
        if( iEtaBin == 13 ) f2->SetParameters( 1, 1.01263, 0.049631, -0.618114, 1.80426);
        else f2->SetParameters( 24.2464, 1.01263, 0.049631, -0.618114, 1.80426);
    }
    else if(iEnergyBin == 7 )
    {
        if( iEtaBin == 8 ) f2->SetParameters( 1, 1.01458, 0.0519279, -0.610171, 1.64167);
        else f2->SetParameters( 9.89067, 1.01458, 0.0519279, -0.610171, 1.64167);

    }
    else if( iEnergyBin == 8 )
        f2->SetParameters( 3.54385, 1.01532, 0.053865, -0.581497, 1.64641);
    else if( iEnergyBin == 9 )
        f2->SetParameters( 1.62573, 1.0182, 0.0581871, -0.592639, 1.50921);

    //MC nominal
    f2->SetParLimits(0, 0.,1000);
    f2->SetParLimits(1,0.8,1.2);
    f2->SetParLimits(2,0.01,0.5);
    f2->SetParLimits(3,-10.,10.);
    f2->SetParLimits(4,1.,10);

    




    //MC material
    /*
    f2->SetParLimits(0,0.,1e6);
    f2->SetParLimits(1,0.8,1.2);
    f2->SetParLimits(2,0.01, 0.5);
    f2->SetParLimits(3,-10.,10.);
    f2->SetParLimits(4,1.,10);
    */
    cout<<" END MC "<<endl;

    
}
     

 

TH1F *hRatio = new TH1F(Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i",iperiod,iEtaBin,iEnergyBin,iBias), "", 50, 0.6, 2.); 
hRatio->Sumw2();

gStyle->SetOptFit(1);
gROOT->ForceStyle();

TCanvas* canvas = new TCanvas("canvas","canvas",800,600) ;
canvas->cd() ; 
TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
pad1->Draw();pad2->Draw();
pad1->cd();
pad1->SetBottomMargin(0.16);
pad2->SetBottomMargin(0.24);
map_histog->SetMarkerStyle(20);
map_histog->SetMarkerSize(1.);
map_histog->GetYaxis()->SetTitle("Events");
map_histog->GetXaxis()->SetTitle("EoverP");
map_histog->Draw("E");
map_histog->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.root",iperiod,iEtaBin,iEnergyBin,iBias));
map_histog->Fit("f2", "r");

for (int ibin = 11; ibin <= 22; ++ibin) {// nominal fit range
    
    int low, high;
    if( map_histog->GetBinContent(ibin) > 0 )
    {
        double res =  (map_histog->GetBinContent(ibin)- f2->Eval( map_histog->GetBinCenter(ibin) ) )/map_histog->GetBinError(ibin);
        hRatio->SetBinContent(ibin, res  );
        hRatio->SetBinError(ibin, map_histog->GetBinError(ibin)  );
        
    }
}       

  pad2->cd();
  hRatio->GetXaxis()->SetTitle((Form("EoverP, #eta bin %i: abs(#eta) in range [%g, %g], E1/E2 bin %i", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1],iEnergyBin )));

    
  hRatio->GetYaxis()->SetTitle("DATA / FIT");
  hRatio->GetYaxis()->SetRangeUser(-50.,50);
  hRatio->GetXaxis()->SetLabelSize(0.1);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetXaxis()->SetTitleSize(0.08);
  hRatio->GetYaxis()->SetTitleSize(0.09);
  hRatio->GetYaxis()->SetTitleOffset(0.6);
  hRatio->GetXaxis()->SetTitleOffset(1.2);
  hRatio->Draw("E1");
  //pad2->Modified();
  //pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();   
 
  canvas->SaveAs(path+Form("fabio_EoverP_period%i_eta%i_energy%i_bias%i.pdf",iperiod,iEtaBin,iEnergyBin,iBias));

 
  delete pad1;
  delete pad2;
  delete canvas;

 double meanValue, errMeanValue;
   
if(iperiod == 0)
{
    if(iEtaBin == 13 && iEnergyBin > 6 )
    {
        meanValue = 0;
        errMeanValue  = 0;
    }
    else
    {
        meanValue = f2->GetParameter(1);
        errMeanValue = f2->GetParError(1);
    }
}
else if(iperiod == 1)
{
    meanValue = f2->GetParameter(1);
    errMeanValue = f2->GetParError(1);
}
     

 parameters.push_back(meanValue);//0
 parameters.push_back(errMeanValue);//1

return parameters;

 
     
 }
//////////////////////////////////////////SAVING HISTOS////////////////////////////////////////

void mySave(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, int iEnergyBin,TString vartagname, TString varname,TString path, int m_layer){	

  gStyle->SetOptFit(0);
  TCanvas *canvas = new TCanvas("scanning","scanning", 800,800);
  canvas->Clear();
  canvas->cd();
  TPad *pad1 =  new TPad("pad1","pad1name",0.01,0.31,0.99,0.99);
  TPad *pad2 =  new TPad("pad2","pad2name",0.01,0.01,0.99,0.41);
  pad1->Draw();pad2->Draw();
  pad1->cd();
  pad1->SetBottomMargin(0.16);
  pad2->SetBottomMargin(0.24);
  //gPad->SetLogy();

  
  map_histog_data->GetYaxis()->SetTitle("Normalised Events");
  map_histog_data->GetXaxis()->SetRangeUser(lowcut, highcut);
  map_histog_data->GetYaxis()->SetRangeUser(-0.1, 2*map_histog_data->GetMaximum());
  map_histog_data->GetYaxis()->SetTitleSize(0.07);
  map_histog_data->GetYaxis()->SetTitleOffset(0.9);
  map_histog_data->DrawNormalized("E1");
  map_histog_mc->SetLineColor(2);
  map_histog_mc->SetMarkerColor(2);
  map_histog_mc->DrawNormalized("histsame");
  
  TLegend *leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->SetLineColor(0);
  leg->SetShadowColor(0);
  leg->SetTextSize(0.035);
  leg->SetHeader(Form("#eta bin %i: abs(#eta) in range [%g, %g]", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1]));
  

  leg->AddEntry(map_histog_mc,"MC","pl");
  if(iBias < 2)
    leg->AddEntry(map_histog_data, "DATA","l");
  else
    leg->AddEntry(map_histog_data, Form("DATA_bias%i",iBias-2),"l");

  pad1->Modified();
  pad1->RedrawAxis();
  pad1->Update();
  pad2->cd();
  
  TH1F* clone_mc = (TH1F*)map_histog_mc;
  TH1F* clone_data = (TH1F*)map_histog_data;
  clone_mc->Scale(clone_data->Integral()/clone_mc->Integral());
  clone_data->Divide(clone_mc);
  
  clone_data->GetXaxis()->SetTitle(vartagname);
  clone_data->GetYaxis()->SetTitle("DATA / MC");
  clone_data->GetXaxis()->SetRangeUser(lowcut, highcut);
  clone_data->GetYaxis()->SetRangeUser(0.7,1.3);
  clone_data->GetXaxis()->SetLabelSize(0.1);
  clone_data->GetYaxis()->SetLabelSize(0.08);
  clone_data->GetXaxis()->SetTitleSize(0.08);
  clone_data->GetYaxis()->SetTitleSize(0.09);
  clone_data->GetYaxis()->SetTitleOffset(0.6);
  clone_data->GetXaxis()->SetTitleOffset(1.2);
  clone_data->Draw("E1");
  pad2->Modified();
  pad2->SetGridy();
  pad2->RedrawAxis();
  pad2->Update();

  TString output_folder = path+Form("calibrated_%s",varname.Data());
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } 
  else {
    gSystem->mkdir(path+Form("calibrated_%s/",varname.Data()), kTRUE);
  }

  if (!(stat(output_folder+Form("/Histograms/layer%i/",m_layer), &sb) == 0 && S_ISDIR(sb.st_mode))) 
    gSystem->mkdir(output_folder+Form("/Histograms/layer%i/",m_layer), kTRUE);

  if(iBias<2){
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s.pdf",m_layer,iEtaBin, iEnergyBin, varname.Data()));
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s.root",m_layer,iEtaBin, iEnergyBin, varname.Data()));
  }
  else {
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s_bias%i.pdf",m_layer,iEtaBin, iEnergyBin, varname.Data(),iBias-2));
    canvas->SaveAs(output_folder+Form("/Histograms/layer%i/weight_graph_eta%i_energy%i_%s_bias%i.root",m_layer,iEtaBin, iEnergyBin, varname.Data(),iBias-2));
  }

  delete canvas;      
  
}



void mySave2D(TH1F* map_histog_mc, TH1F* map_histog_data, int iBias, int iEtaBin, TString vartagname, TString varname, TString path,TString parameter, int m_layer){	


  
  
    
  gStyle->SetOptFit(0);
  TCanvas *canvas2D = new TCanvas("scanning","scanning", 800,800);
  canvas2D->Clear();
  canvas2D->cd();
  if(varname=="m_mass2el") map_histog_mc->GetYaxis()->SetTitle("m_{ee} peak [GeV]");
  if(varname=="m_eOverp") map_histog_mc->GetYaxis()->SetTitle("E/p peak");
  map_histog_mc->GetXaxis()->SetTitle("E1/E2");
  map_histog_mc->SetLineColor(2);
  map_histog_mc->SetMarkerColor(2);
  map_histog_data->SetMarkerColor(1);
  map_histog_data->SetLineColor(1);
  if(varname=="m_mass2el") map_histog_mc->GetYaxis()->SetRangeUser(lowcut, highcut);
  if(varname=="m_eOverp")  map_histog_mc->GetYaxis()->SetRangeUser(0., 5.);

  map_histog_mc->Draw("E1");
  map_histog_data->Draw("E1 samese");
  TLegend *leg2D = new TLegend(0.5,0.6,0.9,0.9);
  if(varname=="m_mass2el") leg2D->SetHeader(Form("#eta bin %i: abs(#eta) in range [%g, %g]", iEtaBin, m_etaBinning[iEtaBin], m_etaBinning[iEtaBin+1]));
  leg2D->SetLineColor(0);
  leg2D->SetShadowColor(0);
  leg2D->SetTextSize(0.035);
  leg2D->AddEntry(map_histog_data,"DATA");
  leg2D->AddEntry(map_histog_mc, "MC");
  leg2D->Draw("");


  TString output_folder = path+Form("calibrated_%s",varname.Data());
  struct stat sb;
  if (stat(output_folder, &sb) == 0 && S_ISDIR(sb.st_mode)) {
    gSystem->cd(output_folder);    
  } else {
    gSystem->mkdir(path+Form("calibrated_%s/",varname.Data()), kTRUE);
  }

  if (!(stat(output_folder+Form("/Histograms2D/layer%i/",m_layer), &sb) == 0 && S_ISDIR(sb.st_mode))) 
    gSystem->mkdir(output_folder+Form("/Histograms2D/layer%i/",m_layer), kTRUE);


  if(iBias<2){  
    map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo%s_mc_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo%s_data_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_%s.pdf",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_%s.root",m_layer,parameter.Data(),iEtaBin,varname.Data()));
    
    //save the parameters for MC
    if(varname=="m_mass2el" || varname=="m_eOverp"){      
      map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_mc"+parameter+"_%i_%s.root",m_layer,iEtaBin,varname.Data()));
      map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_data"+parameter+"_%i_%s.root",m_layer,iEtaBin,varname.Data()));
    }
  }
	    
  else{  
      
    cout<<"FABIOLIFE"<<endl; 
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_bias%i_%s.pdf",m_layer,parameter.Data(),iEtaBin,iBias,varname.Data()));
    canvas2D->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_graph%s_%i_bias%i_%s.root",m_layer,parameter.Data(),iEtaBin,iBias,varname.Data()));
    
    //save the parameters for MC
    if(varname=="m_mass2el" || varname=="m_eOverp"){      
      map_histog_mc->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_mc"+parameter+"_bias%i_%i_%s.root",m_layer,iEtaBin,iBias,varname.Data()));
      map_histog_data->SaveAs(output_folder+Form("/Histograms2D/layer%i/weight_histo_data"+parameter+"_bias%i_%i_%s.root",m_layer,iEtaBin,iBias,varname.Data()));
    }
  }  
}

vector<double>params(int iBias, int iEtaBin, int iEnergyBin, TString type, TString name)
{

TString eta = Form("eta%d",iEtaBin);
TString energy = Form("energy%d",iEnergyBin);
TString bias = Form("bias%d",iBias);

TString path;    
 
if(type == "data17AndMaterial")
{
    cout<<"testefile"<<endl;
    path = Form("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/material_reBinning_txt/input_data17_%s.txt", name.Data());
}
else
{
    path = Form("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/material_reBinning_txt/input_%s_%s.txt", type.Data(), name.Data());
    
    cout<<path<<endl;
}
    
    
cout<<"teste1"<<endl;
/*if(type == "data17AndMaterial")
    path = Form("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/material_reBinning_txt/out_data17_%s.txt",name.Data());
else
    path = Form("/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/material_reBinning_txt/input_%s_%s.txt", type.Data(), name.Data());
*/

cout<<eta<<energy<<bias<<endl;
    
ifstream file( path  );

        vector<double> parameters;

        std::string line;
        //std::vector<std::string> myLines;
        string str;
        //cout<<line<<endl;
        while (std::getline(file, line))
        {
                TString found = line;
                //cout<<line<<endl;
                if(found.Contains(eta) && found.Contains(energy) && found.Contains(bias) )
                {
                        str = line;
                        break;
                }
        }

        //cout<<line<<endl;

    std::istringstream ss(str);
    std::vector<string> tmp;
    std::vector<double> params;

    string v;
    while (ss)
    {
        ss >> v;
        //cout<<" v " << v << endl;
        tmp.push_back(v);
    }

    int it = 6;
    while(it < 11)
    {
        double convertedValue = std::stod(tmp.at(it));
        params.push_back(convertedValue);
        it++;
    }

    for(int i  = 0; i < params.size(); i++) cout<< "pars " << params.at(i) << endl;

    return params;
}
                         

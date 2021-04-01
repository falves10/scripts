#include "TMinuit.h"
void runMC(TFile *f, int i, int j, ofstream & file);
vector<double>  params(int i, int j, int k);


void testFit()
{

TString path = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/";


TFile *f_in;  
      
cout<<" TESTE " << endl;
f_in = new TFile(path+"weight_histo_data_m_eOverp.root");

double c;
double mean;
double sigma;
double alpha;
double n;
    
TH1F *h;

ofstream myfile;
myfile.open ("new_example.txt");

bool useParamsFromFile = true;
TString type = "data";

Int iEta = 4;
Int iEnergy = 10;
Int iBias =  16;
    
//c->Divide(2,2);
for(int i = 3;i < iEta; i++)
{
    for(int j = 0; j < iEnergy; j++)
    {
        for(int k = 0; k < iBias ;k++)
        {
            TCanvas *c = new TCanvas();
            gStyle->SetOptFit(1111);
            
            if(type = "data")
                h = (TH1F *)f_in->Get(Form("%s_eta_%d_energy_%d_bias%d",type.Data(), i, j, k)); // data
            else h = (TH1F *)f_in->Get(Form("mc_eta_%d_energy_%d_bias0",i, j));
                

            cout<<"eta, energy, bias " << i << " " << j << " " << k << endl; 
            cout<<h->GetMean()<<" " << h->GetRMS() << endl;
	       
            h->Draw();

	       TF1 *f2 = new TF1("f2","crystalball",0.9,1.3);

	       if(h->GetEntries() < 1000) continue;
	       else
	       {
               if(useParamsFromFile)
               {
                    vector<double> par = params(i, j, k);
	               f2->SetParameters( par[0], par[1], par[2], par[3], par[4]);
               }
               else
               {
                   f2->SetParameters(2751.33, 1.00918, 0.0653444, -0.700741, 1.66549);
               }
               
               if(type = "data")
               {
                   f2->SetParLimits(0, 0.,1e10);      
                   f2->SetParLimits(1,0.8,1.2);
                   f2->SetParLimits(2,0.01,0.5);
                   f2->SetParLimits(3,-10.,10.);
                   f2->SetParLimits(4,1.,10);
               }
               else
               {
                   f2->SetParLimits(0, 0.,1000);
                   f2->SetParLimits(1,0.8,1.2);
                   f2->SetParLimits(2,0.01,0.5);
                   f2->SetParLimits(3,-10.,10.);
                   f2->SetParLimits(4,1.,10);

               }
	           
           }
	
	       h->Fit("f2", "r");

	       c->SaveAs(Form("./out_nominal/mc_eta%d_energy%d_bias%d.pdf", i, j, k));

 	       TString status_fit = gMinuit->fCstatu;           
 	       myfile << " fit status " << status_fit << " eta" << i << " energy" << j << " bias" << k << " " << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " << f2->GetParameter(3) << " " << f2->GetParameter(4) << endl;             
        }
    }
      
}

myfile.close();

}

void runMC(TFile *f, int i, int j, ofstream &file)
{
    TH1F *h;
    TCanvas *c1 = new TCanvas();
    h = (TH1F *)f->Get(Form("mc_eta_%d_energy_%d_bias0",i, j));
    cout<<"MC: eta, energy, bias " << i << " " << j << " "  << endl; 
    h->Draw();
   
 
    TF1 *f2 = new TF1("f2","crystalball",0.9,1.3);
  
    /*if( (i == 0 && j == 2))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );
    else*/
   

    /* MATERIAL EL 
    if( (i == 0 && j == 4) || (i == 0 && j == 7)  )  
    	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
    else if ( (i == 1 && j == 9) )
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.2 , 2 );
    else
	f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
    */

    //MATERIAL FMX
    /*
    if((i == 0 && j == 8) || (i == 0 && j == 9) )
    	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    else if((i == 1 && j == 6) )
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   else if ((i == 2 && j == 1)) 
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 ); 
   else if ((i == 2 && j == 9))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -7 , 2 );
   else 
	f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   */	

   //MATERIAL N
   /*
 if( (i == 0 && j == 4) ||(i == 0 && j == 8) || (i == 0 && j == 9) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    else if((i == 1 && j == 6) || (i == 1 && j == 4) || (i == 1 && j == 9)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   else if ((i == 2 && j == 1))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 9))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -7 , 2 );
   else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   */

  //MATERIAL IBL
  /*
   if( (i == 0 && j == 0) || (i == 0 && j == 4) ||(i == 0 && j == 8) || (i == 0 && j == 9) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   else if((i == 1 && j == 5)  )
   	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.4 , 2 );
   else if ((i == 2 && j == 1))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 7))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else if ((i == 2 && j == 8))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else if ((i == 2 && j == 9))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else
   	f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   */


   //MATERIAL PP0
   /*if( (i == 0 && j == 0) || (i == 0 && j == 4) ||(i == 0 && j == 8) || (i == 0 && j == 9) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   else if((i == 1 && j == 5)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.4 , 2 );
   else if ((i == 2 && j == 1))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 7))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else if ((i == 2 && j == 8))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.45 , 2 );
   else if ((i == 2 && j == 9))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );

    */

    // MATERIAL A	
    /*if( i == 0 ){
    
    if( (i == 0 && j == 3) || (i == 0 && j == 4))
   	f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );    
    else if ((i == 0 && j == 5))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );

    else if ((i == 0 && j == 7))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );    
    else
	f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );
	
    }
    else if ( i == 1 )
    {
	if (j == 2 || j == 3 || j == 4)
        	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
	else if (j == 6 || j == 7 || j == 9)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    	else
        	f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );	
    }
    else if ( i == 2 )
    {
	//if (j == 0 || j == 1)
	//	f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );
	//else if (j == 2 || j == 3 || j == 4 ) 
	//	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 ); 
	//else if (j == 5 || j == 6 || j == 7 || j == 8 || j == 9)
		//f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 ); 
    	if (j == 5 || j == 7)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.07 , 2 );
	else if ( j == 9 )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.2 , 2 );
	else 
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    }	
    else if ( i == 3 )
    {
	if(j == 4 || j == 7)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.05 , 2 );
	else
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    }
   else if ( i == 4 )
   {
        //if(j == 4 || j == 7)
        //        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.05 , 2 );
        //else
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
    } 
   else if ( i == 5 )
   {
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   }
   */

   //MATERIAL A MC16D   
   /*if( i == 0 ){

    if( (i == 0 && j == 3) || (i == 0 && j == 4) || (i == 0 && j == 5))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
    else if ((i == 0 && j == 9))
    	f2->SetParameters(1, h->GetMean()  , 0.01 , -2.2 , 2 );

    else if ((i == 0 && j == 7))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
    else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );

    }
    else if ( i == 1 )
    {
        if (j == 2 || j == 3 || j == 4)
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
        else if (j == 6 || j == 7)
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
	else if (j == 8 || j == 9)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -2.1 , 2 );
        else
                f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 );
    }
    else if ( i == 2 )
    {
        if (j == 5 || j == 7 )
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.2 , 2 );
	else if (j == 6 || j == 8)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
	else if (j == 9)
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.45 , 2 ); 
        else
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    }
    else if ( i == 3 )
    {
        if(j == 1 || j == 2 || j == 4 || j == 7)
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.05 , 2 );
        else
                f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    }
   else if ( i == 4 )
   {
      	if(j == 3 || j == 6 )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
	else if (j == 4 ||  j == 7 )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -1. , 2 );
	else if (  j == 8 )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.8 , 2 );
	else 
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    }
    */

    //MATERIAL EL (MC16D)
    /*if( (i == 0 && j == 1) || (i == 0 && j == 8) || (i == 0 && j == 9)  )  
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
    else if ( (i == 1 && j == 0) || (i == 1 && j == 8) )
    	f2->SetParameters(1, h->GetMean()  , 0.01 , -1.5 , 2 );
    else if ( (i == 1 && j == 4) )
	f2->SetParameters(1, h->GetMean()  , 0.01 , -1.5 , 2 );
    else if ( (i == 1 && j == 9) )
	f2->SetParameters(1, h->GetMean()  , 0.01 , -2.1 , 2 );
    else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );  

    */

    //MATERIAL FMX (MC16D) 
    /* if((i == 0 && j == 4) || (i == 0 && j == 7) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    else if((i == 1 && j == 0) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -2.7 , 2 );
   else if( (i == 1 && j == 1) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if( (i == 1 && j == 6) || (i == 1 && j == 7) || (i == 1 && j == 9)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 1) || (i == 2 && j == 2) || (i == 2 && j == 6)) 
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 ); 
   else if((i == 2 && j == 8))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.45 , 2 );
   else 
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
    */
  
    //MATERIAL N (MC16D) 
    /*if( (i == 0 && j == 4) ||(i == 0 && j == 8) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
    else if((i == 1 && j == 3)  )
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if((i == 1 && j == 1)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.2 , 2 );
   else if((i == 1 && j == 8)  )
   	f2->SetParameters(1, h->GetMean()  , 0.01 , -2 , 2 ); 
   else if((i == 1 && j == 6) || (i == 1 && j == 4)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   else if ((i == 2 && j == 1))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 2))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.6 , 2 ); 
   else if ((i == 2 && j == 3))
        f2->SetParameters(3.68242e+00, 9.85088e-01   , 1.17157e-01  , -4.64469e-01 , 2.95459e+00  );
   else if ((i == 2 && j == 6))
	f2->SetParameters(1 , 1.02990e+00 , 0.06  , -1. , 2  );
   else if ((i == 2 && j == 9))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -7 , 2 );
   else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   */

   //MATERIAL IBL (MC16D)
   /*if( i == 0 )
   {
	if( j == 1 || j == 3  )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
	else if ( j == 5 || j == 7  )
		f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
	else
		f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   }
   else if((i == 1 && j == 5) )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 1 && j == 9)  )
	f2->SetParameters(1.62322e+00, h->GetMean()  , 1.76737e-01 , -6.97698e-01 , 1 );
   else if ((i == 2 && j == 1) || (i == 2 && j == 2))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 3))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.8 , 2 );
  else if ((i == 2 && j == 6))
        f2->SetParameters(1, h->GetMean()  , 0.03 , -0.3 , 2 ); 
  else if ((i == 2 && j == 7))
        f2->SetParameters(1, h->GetMean()  , 0.04 , -0.3 , 2 );
   else if ((i == 2 && j == 8))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else if ((i == 2 && j == 9))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 );
   else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
   */

    //MATERIAL PP0 (MC16D)
    //if( (i == 0 && j == 0) || (i == 0 && j == 4) ||(i == 0 && j == 8) || (i == 0 && j == 9) )
      //  f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
   if(i == 0)
   {
	if( j == 1) f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
	else if ( j == 5) f2->SetParameters(1, h->GetMean()  , 0.01 , -0.1 , 2 );
	else f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );	
   }
   else if((i == 1 && j == 5)  || (i == 1 && j == 9)  )
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.4 , 2 );
   else if((i == 1 && j == 6))
	f2->SetParameters(1, h->GetMean()  , 0.01 , -0.3 , 2 ); 
   else if ((i == 2 && j == 1))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else if ((i == 2 && j == 0))
        f2->SetParameters(1, h->GetMean()  , 0.09 , -7 , 2 );
  else if ((i == 2 && j == 3))
        f2->SetParameters(1, h->GetMean()  , 0.1 , -1 , 2 ); 
  else if ((i == 2 && j == 7))
        f2->SetParameters(1, 0.99  , 0.05 , -2 , 2 );
   else if ((i == 2 && j == 8))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.45 , 2 );
   else if ((i == 2 && j == 9))
        f2->SetParameters(1, h->GetMean()  , 0.01 , -0.5 , 2 );
   else
        f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );

	


    //material A: f2->SetParLimits(0,0.,1e3);
    f2->SetParLimits(0,0.,1e6);
    f2->SetParLimits(1,0.8,1.2);
    f2->SetParLimits(2,0.01, 0.5);
    f2->SetParLimits(3,-10.,10.);
    f2->SetParLimits(4,1.,10);

    
    
    h->Fit("f2", "r");

    file << " eta " << i << " energy " << j << " " << f2->GetParameter(0) << " " << f2->GetParameter(1) << " " << f2->GetParameter(2) << " " << f2->GetParameter(3) << " " << f2->GetParameter(4) << endl; 

    
}


vector<double>  params(int i, int j, int k)
{

TString eta = Form("eta%d",i);
TString energy = Form("energy%d",j);
TString bias = Form("bias%d",k);

//cout<<"teste1"<<endl;

TString path = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/data_nominal_eta_2_3.txt";

ifstream file( path  );

	vector<double> parameters;

	std::string line;
	//std::vector<std::string> myLines;
	string str;
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
    //cout<<" SIZE " << params.size() << endl; 

    return params;
}

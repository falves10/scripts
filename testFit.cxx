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


//c->Divide(2,2);
for(int i = 0;i < 15; i++)
{
    for(int j = 0; j < 10; j++)
    {
        for(int k = 0; k < 16 ;k++)
        {
            TCanvas *c = new TCanvas();
            gStyle->SetOptFit(1111);
	       //h = (TH1F *)f_in->Get(Form("data_eta_%d_energy_%d_bias%d",i, j, k)); // data
            h = (TH1F *)f_in->Get(Form("mc_eta_%d_energy_%d_bias0",i, j));
            cout<<"eta, energy, bias " << i << " " << j << " " << k << endl; 
            cout<<h->GetMean()<<" " << h->GetRMS() << endl;
	       h->Draw();

	       TF1 *f2 = new TF1("f2","crystalball",0.9,1.3);

	       if(h->GetEntries() < 1000) continue;
	       else
	       {

		//f2->SetParameters(1, h->GetMean()  , 0.01 , -1 , 2 );
 		/* MC nominal
		if( j == 0 )
		{
			if( i  == 1 )
				f2->SetParameters(20, 0.992018, 0.0413628, -0.3, 2.72838 );
			else if ( i  == 2 ) f2->SetParameters(20, 0.992018, 0.0413628, -0.2, 2.72838 );
			else f2->SetParameters( 1e+2 , 0.992018, 0.0413628, -1, 2. );
		}
                else if ( j == 1 )
		{
			if( i  == 2 )
				f2->SetParameters(155586, 0.996261, 0.0414388, -0.3, 2.53711 );
			else
				f2->SetParameters(155586, 0.996261, 0.0414388, -0.905716, 2.53711 );
		}
		else if ( j ==  2 )  
		{
			if( i == 1 || i == 2)
				f2->SetParameters(1e2, h->GetMean(), 0.0627053, -8.74546e-01, 2.03312e+00 );
			else 
				f2->SetParameters(1e2, h->GetMean(), 0.0627053, -8.74546e-01, 2.03312e+00 );
		}    
		else if ( j == 3 )
		{
			if( i == 2 ) f2->SetParameters(1e2, h->GetMean(), 0.06, -0.7, 2.0 );
			else f2->SetParameters(1e2, h->GetMean(), 0.06, -0.9, 2.0 );
		}
		else if ( j == 4 )
		{
			if( i == 2 )
				f2->SetParameters(1e+1, h->GetMean(), 0.06, -0.3, 2.0 ); 
			else
				f2->SetParameters(1e+1, h->GetMean(), 0.06, -0.5, 2.0 );
		}
		else if ( j == 5 )	
			f2->SetParameters(1e2, h->GetMean(), 0.06, -1, 2.0 );
		else if ( j == 6 )
			f2->SetParameters(1e2, h->GetMean(), 0.06, -0.8, 2.0 );	
		else if ( j == 7 )
			f2->SetParameters(10, h->GetMean(), 0.07, -1, 2. );
		else if ( j == 8 || j == 9 )
			f2->SetParameters(1, h->GetMean(), 0.07, -0.3, 2.0 );

		*/


		if( j == 0 )
    			f2->SetParameters( 24.491, 0.996921, 0.041169, -1.02395, 2.79162 );
else if( j == 1 )
    f2->SetParameters( 71.8089, 1.00044, 0.0411136, -0.893548, 2.61719 );
else if( j == 2 )
    f2->SetParameters( 116.533, 1.00307, 0.041934, -0.800426, 2.39107 );
else if( j == 3 )
    f2->SetParameters( 119.281, 1.00593, 0.0434437, -0.727336, 2.23322 );
else if( j == 4 )
    f2->SetParameters( 88.0382, 1.0081, 0.0452797, -0.670262, 2.13681);
else if( j == 5 )
    f2->SetParameters( 50.8651, 1.01055, 0.0474722, -0.641028, 1.95129);
else if( j == 6 )
{
	if( i == 13 ) f2->SetParameters( 1, 1.01263, 0.049631, -0.618114, 1.80426);
	else f2->SetParameters( 24.2464, 1.01263, 0.049631, -0.618114, 1.80426);
}
else if( j == 7 )
{
	if( i == 8 ) f2->SetParameters( 1, 1.01458, 0.0519279, -0.610171, 1.64167);
	else f2->SetParameters( 9.89067, 1.01458, 0.0519279, -0.610171, 1.64167);

}
else if( j == 8 )
    f2->SetParameters( 3.54385, 1.01532, 0.053865, -0.581497, 1.64641);
else if( j == 9 )
    f2->SetParameters( 1.62573, 1.0182, 0.0581871, -0.592639, 1.50921);


f2->SetParLimits(0, 0.,1000);
f2->SetParLimits(1,0.8,1.2);
f2->SetParLimits(2,0.01,0.5);
f2->SetParLimits(3,-10.,10.);
f2->SetParLimits(4,1.,10);



	 	// data nominal
		/*if( j == 0 )
			f2->SetParameters(54927.8, 0.992018, 0.0413628, -1.04522, 2.72838 );	
		else if ( j == 1 )
			f2->SetParameters(155586, 0.996261, 0.0414388, -0.905716, 2.53711 ); 
		else if ( j == 2 )	
			f2->SetParameters(232059, 0.99981, 0.0427053, -0.815207, 2.28122 );
		else if ( j == 3 )
			f2->SetParameters(222901, 0.995381, 0.0438607, -0.74959, 2.10982 );
		else if ( j == 4 )
			f2->SetParameters(155332, 0.998958, 0.0460873, -0.701446, 1.99597 );
		else if ( j == 5 )
			f2->SetParameters(85932.3, 1.00334, 0.0483657, -0.660556, 1.94629);
		else if ( j == 6 )
			f2->SetParameters(39958.5, 1.00769, 0.0507526, -0.633211, 1.89242);
		else if ( j == 7 )
			f2->SetParameters(15958.2, 1.01376, 0.0551146, -0.661654, 1.66726);
		else if ( j == 8 )
			f2->SetParameters(5654.9, 1.02723, 0.0589129, -0.592765, 2.31006);
		else if ( j == 9 )
			f2->SetParameters(2751.33, 1.00918, 0.0653444, -0.700741, 1.66549);
		*/
		
		/*
		if( j  <= 7 )
			f2->SetParameters( 1.93671e+05 ,  h->GetMean()   , 0.06  , -1  , 2 );
		else
			f2->SetParameters( 1.93671e+05 ,  h->GetMean()   , 0.08  , -1  , 2 );
		*/
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

TString path = "/publicfs/atlas/atlasnew/higgs/hgg/fabiolucio/EgammCalibration/Codes_ntuples/condor/calibrated_m_eOverp/Histograms/merge_test/material_reBinning_txt/input_data17_A.txt";

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

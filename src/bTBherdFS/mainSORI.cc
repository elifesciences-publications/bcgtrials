#include <iostream>
#include "../../include/PATHS.h"
#include "../../include/bTBICBMV8.h"
#include <string.h>

// Delay Distributions

// gsl_ran_lognormal(gsl_r2,2.7085,0.7)

//gsl_ran_exponential(gsl_r2,364/2)


int reacto(gsl_rng *r, double p, int n)
{
    // Pathological case (Approx to 1)
    if(n == 0){return 1;}
    
    double probs[n-1];
    
    // Conditional probability of (i+1) reactors
    // given that we have more than 0
    double accum = 0;
    
    for(int i=0; i < n; i++)
    {
        probs[i] = gsl_ran_binomial_pdf(i+1,p,n) / (1 - pow((1 - p),n));
        accum += probs[i];
    }
    
    
    double ranx = gsl_ran_flat(r,0,1);
    accum = 0;
    
    for(int i=0; i < n;i++)
    {
        accum += probs[i];
        if( accum > ranx)
        {
            //cout << i+1 << endl;
            return((i+1));
        }
    }
    return(n);
} 


int main (int argc, char * const argv[]) 
{
    
	bool DiscoFlag = true;
	bool save_ts        = false;
	
	MTRand mrandy;
	gsl_rng * gsl_r2;
	
	gsl_rng_default_seed = mrandy.rand() * sizeof(long unsigned int);
	gsl_r2 = gsl_rng_alloc(gsl_rng_mt19937);
    
    //double max_sim_time = 10;
	
	const int no_of_PTI = 3;
	
	int PTI[no_of_PTI] = {1,2,4};
	
	int Seeders = 0;
	double Turnover = 0.0;
	// Standard and severe interpretations of skin test
	double Standard[2] = {0.0,0.0}; 
	double Repeat_Sensitivity[2] = {0.59,0.70};
	
	double Severe[2] = {0.0,0.0};
	double beta_o = 0.0;
	double beta_r = 0.0;
	double beta_i = 0.0;
	double beta_ov1 = 0.0;
	double beta_ov2 = 0.0;
	double beta_rv = 0.0;
	double beta_iv = 0.0;
	double rel1=0.0,rel2=0.0,rel3=0.0,rel4=0.0;	
	
	double q = 0.0;
    
	double p_confirm[3] = {0.0,0.0,0.0};
	double p_occult = 0.0;
	double s_slaughter = 0.1;
	
	double delta_age = 0.0;
	double shad_lay = 0.0;
	double p_shadow = 0.0;
	double pinf_move = 0.0;
	
	double xinf_PTI[4];
		
	int no_of_reps=0;
		
	double TO = 0.0;
	double TR = 0.0;
	double TV1 = 0.0;
	double TV2 = 0.0;
	double TOV1 = 0.0;
	double TOV2 = 0.0;
	double TRV = 0.0;


	int low_n,high_n,low_PTI,high_PTI;
	
	// Modifiers to testing sequence
	// Interval for final short interval test
	bool change_SIT = false;
	double second_SIT = 120.0;
	
	// Number of clear tests (at standard interpretation) to end breakdown
	int must_clear = 1;
	// Boolean variable to allow end of breakdown on clear test at severe interpretation
	bool clear_from_severe = false;
    
	char *filename;
	
	if(argc < 36)
    { 
		cout << "bTBherd <Sensitivity Standard> <Specificity Standard> <Sensitivity Severe> <Specificity Severe> <O> <R> <V1> <V2> <OV1> <OV2> <RV> <beta_o> <beta_r> <beta_i> <beta_ov1> <beta_ov2> <beta_rv> <beta_iv> <q> <Xinf PTI1> <Xinf PTI2> <Xinf PTI4> <rel1> <rel2> <rel3> <rel4> <p_confirm_I> <p_confirm_R> <p_confirm_O> <p_slaughter> <delta_age> <p_shadow> <shadow_delay> <pinf_move> <no_of_reps> <outputfile> <Fixed SIT> <final SIT> <Change Final (Bool)> <End on Severe (Bool)> <must_clear>" << endl;
		exit(EXIT_FAILURE);
    }
    else
    {
		Standard[0] = atof(argv[1]);
		Standard[1] = atof(argv[2]);
		Severe[0] = atof(argv[3]);
		Severe[1] = atof(argv[4]);
		TO = atof(argv[5]);
		TR = atof(argv[6]);
		TV1 = atof(argv[7]);
		TV2 = atof(argv[8]);
		TOV1 = atof(argv[9]);
		TOV2 = atof(argv[10]);
		TRV = atof(argv[11]);
		beta_o = atof(argv[12]);
		beta_r = atof(argv[13]);
		beta_i = atof(argv[14]);
		beta_ov1 = atof(argv[15]);
		beta_ov2 = atof(argv[16]);
		beta_rv = atof(argv[17]);
		beta_iv = atof(argv[18]);
		
		q      = atof(argv[19]);
		
		xinf_PTI[0] = atof(argv[20]);
		
		xinf_PTI[1] = atof(argv[21]);
        
		xinf_PTI[3] = atof(argv[22]);
		
		xinf_PTI[2] = 0.0;
		rel1        = atof(argv[23]);
		rel2        = atof(argv[24]);
		rel3        = atof(argv[25]);
		rel4        = atof(argv[26]);


		p_confirm[0] = atof(argv[27]);
		
		p_confirm[1] = atof(argv[28]);
		
		p_confirm[2] = atof(argv[29]);
		
		
		s_slaughter = atof(argv[30]);
		
		delta_age = atof(argv[31]);
		
		p_shadow  = atof(argv[32]);
		
		shad_lay = atof(argv[33]);
		
		pinf_move = atof(argv[34]);
		
		no_of_reps = atoi(argv[35]);
		
		filename = argv[36];
		

		low_PTI = 0;
		high_PTI = 2;
			
		low_n = 30;
		high_n = 330;
			
			
		}
		
		cout << "Sensitivity: " << Standard[0] << endl;
		cout << "Specificity: " << Standard[1] << endl;
		cout << "Sensitivity (Severe): " << Severe[0] << endl;
		cout << "Specificity (Severe): " << Severe[1] << endl;
		
		//cout << "Repeat Sensitivity: " << Repeat_Sensitivity[0] << endl;
		//cout << "Repeat Sensitivity (Severe): " << Repeat_Sensitivity[1] << endl;
        
		
		cout << "TO: " << TO << endl;
		cout << "TR: " << TR << endl;
		cout << "TV1: " << TV1 << endl;
		cout << "TV2: " << TV2 << endl;
		cout << "TOV1: " << TV1 << endl;
		cout << "TOV2: " << TV2 << endl;
		cout << "TRV: " << TRV << endl;
				
					
		cout << "beta_o: " << beta_o << endl;
		cout << "beta_r: " << beta_r << endl;
		cout << "beta_i: " << beta_i << endl;
		cout << "beta_ov1: " << beta_ov1 << endl;
		cout << "beta_ov2: " << beta_ov2 << endl;
		cout << "beta_rv: " << beta_rv << endl;
		cout << "beta_iv: " << beta_ov1 << endl;
		
		cout << "q: " << q << endl;
		cout << "Xinf PTI1: " << xinf_PTI[0] << endl;
		cout << "Xinf PTI2: " << xinf_PTI[1] << endl;
		cout << "Xinf PTI3: " << xinf_PTI[2] << endl;
		cout << "Xinf PTI4: " << xinf_PTI[3] << endl;
		cout << "Rel1: " << rel1 << endl;
		cout << "Rel2: " << rel2 << endl;
		cout << "Rel3: " << rel3 << endl;
		cout << "Rel4: " << rel4 << endl;

		cout << "Confirmation I: " << p_confirm[0] << endl;
		cout << "Confirmation R: " << p_confirm[1] << endl;
		cout << "Confirmation O: " << p_confirm[2] << endl;
				
		cout << "Slaughterhouse Sensitivity: " << s_slaughter << endl;
		cout << "Delta Age (months): " << delta_age << endl;
		cout << "p_shadow " << p_shadow << endl;
		cout << "Shadow Delay " << shad_lay << endl;
		cout << "Pinf_move " << pinf_move << endl;
						
		cout << "No. of Reps: " << no_of_reps << endl;
		
		cout << low_n << ' ' << high_n << endl;
		cout << low_PTI << ' ' << high_PTI << endl;
		//cout << SIT[0] << endl;
		
		
	ofstream breakdownFile;
	ofstream breakdownFile2;
	
	breakdownFile.open(filename);
	strcat(filename,"Fu");
	breakdownFile2.open(filename);
    
	breakdownFile << "HerdSize,PTI,Prolonged,Recurr,Break6,Break12,Break24,Forward_Trans,Reactors_at_Start,BreakTime,Slaughter,Confirmation,Confirmed,Breaklength,Reactors,Visits,Tests,Burden,Onward_Trans" << endl;  	
	breakdownFile2 << "Replicate,HerdSize,PTI,Breaklength,Break_Recurr,Onward_Trans,Reactors_at_Start,Reactors_at_VE6M,Reactors_at_VE12M,Reactors,Visits,Tests,Burden,Forward_Trans,Slaughter,Confirmed" << endl;
    
    
	//const int no_of_PTI=1;
    //int PTI[no_of_PTI] = {1};
	
    // Herds, Gamma_O, Gamma_R, Gamma_V1, Gamma_V2, Gamma_OV1, Gamma_OV2, Gamma_RV, 
    //		  Gamma_SO, Gamma_SR, Gamma_SV1, Gamma_SV2, Gamma_SOV1, Gamma_SOV2, Gamma_SRV,
    // input, output, full_save, 
    // latent_R, latent_V1, latent_V2, latent_OV1, latent_OV2, latent_RV, 
    // latent_SR, latent_SV1, latent_SV2, latent_SOV1, latent_SOV2, latent_SRV
    // verbose)
    

	//TOF = 0.1;
	//TR = 0.1;
	//beta_r = 0.0;
	//beta_i = 0.0;
	//q = 1.0;

	bTBICBM model(1,
	0,0,0,0,0,0,0,
	0,0,0,0,0,0,0,
	PATH_TO_DATA,filename,save_ts,
	TO*364.0,TR*364.0,TV1*364.0,TV2*364.0,TOV1*364.0,TOV2*364.0,TRV*364.0,
	TO*364.0*shad_lay,364.0*shad_lay*TR,TV1*364.0,TV2*364.0,TOV1*364.0,TOV2,364.0*shad_lay*TR,
	true,pinf_move, delta_age);
	
	model.set_confirmation(p_confirm);
	
	//model.set_transmission(0.0,beta_r,beta_i,0.0,0.0,beta_r,beta_i,q,xinf_PTI[0]);
	//model.set_transmission(0.0,beta_r,beta_i,0.0,0.0,beta_r,beta_i,q,xinf_PTI[0]);
	//model.set_transmission(beta_o,beta_i,beta_r,beta_ov1,beta_ov2,beta_rv,beta_iv,q,0.0);
	
	model.set_forcing(false);
	//model.set_demo(LifeExp, 50.0,1.5356);
	//model.set_demo(LifeFixed, 50.0,0,(1.0/50.0)*2000.0);
	//model.set_demo(LifeExp, 50.0,0,(1.0/50.0)*2000.0);
	//model.set_demo(LifeGamma, 50.0,3.0,(1.0/50.0)*2000.0);
	//model.set_demo(LifeExp, 50.0,0.0);
	
	double parish_testing = 0.0;
	
	//const int no_of_N = 24;
	//int N_size[no_of_N] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400};
	
	//const int no_of_N = 12;
	//int N_size[no_of_N] = {10,20,30,40,50,100,150,200,250,300,350,400};
	
	// Ensemble sums for all N, PTI
	
	int no_of_H = 7;
	
	
	int ensemble_break6U[(no_of_H*no_of_PTI)];
	int ensemble_break12U[(no_of_H*no_of_PTI)];
	int ensemble_break24U[(no_of_H*no_of_PTI)];
	int ensemble_prolongedU[(no_of_H*no_of_PTI)];
    
	int ensemble_break6C[(no_of_H*no_of_PTI)];
	int ensemble_break12C[(no_of_H*no_of_PTI)];
	int ensemble_break24C[(no_of_H*no_of_PTI)];
	int ensemble_prolongedC[(no_of_H*no_of_PTI)];
    
	int ensemble_Confirmed[(no_of_H*no_of_PTI)];
	int ensemble_Unconfirmed[(no_of_H*no_of_PTI)];
	int ensemble_ConfirmedS[(no_of_H*no_of_PTI)];
	int ensemble_UnconfirmedS[(no_of_H*no_of_PTI)];
	
	//int ensemble_Slaughter[(no_of_H*no_of_PTI)];
	
	double ensemble_ReactorsU[(no_of_H*no_of_PTI)];
	double ensemble_ReactorsC[(no_of_H*no_of_PTI)];
    
	double ensemble_BreaklengthsU[(no_of_H*no_of_PTI)];
	double ensemble_BreaklengthsC[(no_of_H*no_of_PTI)];
    
	int num_confirmed[(no_of_H*no_of_PTI)];
	int num_unconfirmed[(no_of_H*no_of_PTI)];
	
	int ensemble_VisitsU[(no_of_H*no_of_PTI)];
	int ensemble_VisitsC[(no_of_H*no_of_PTI)];
    
	int ensemble_TestsU[(no_of_H*no_of_PTI)];
	int ensemble_TestsC[(no_of_H*no_of_PTI)];
    
	double ensemble_ReactorsTotU[(no_of_H*no_of_PTI)];
	double ensemble_ReactorsTotC[(no_of_H*no_of_PTI)];
    
	int ensemble_BurdenU[(no_of_H*no_of_PTI)];
	int ensemble_BurdenC[(no_of_H*no_of_PTI)];
    
	int ensemble_Forward_TransU[(no_of_H*no_of_PTI)];
	int ensemble_Forward_TransC[(no_of_H*no_of_PTI)];
    
	int ensemble_Onward_TransU[(no_of_H*no_of_PTI)];
	int ensemble_Onward_TransC[(no_of_H*no_of_PTI)];
	
	/*
     int sum_Reactors[(no_of_H*no_of_PTI)];
     int sum_Breaklengths[(no_of_H*no_of_PTI)];
     int sum_Visits[(no_of_H*no_of_PTI)];
     int sum_Tests[(no_of_H*no_of_PTI)];
     int sum_Reactors_at_Start[(no_of_H*no_of_PTI)];
     int sum_Burden[(no_of_H*no_of_PTI)];
     int sum_Forward_Trans[(no_of_H*no_of_PTI)];
     int sum_Onward_Trans[(no_of_H*no_of_PTI)];
     */
	
	// Initialise all ensemble variables to zero
	
	for(int k=0; k < no_of_PTI; k++)
	{
		for(int n=0; n < no_of_H; n++)
		{
            
            ensemble_break6U[index2(k,n,no_of_H)] = 0;
            ensemble_break12U[index2(k,n,no_of_H)] = 0;
            ensemble_break24U[index2(k,n,no_of_H)] = 0;
            ensemble_prolongedU[index2(k,n,no_of_H)] = 0;
			
            ensemble_break6C[index2(k,n,no_of_H)]	= 0;
            ensemble_break12C[index2(k,n,no_of_H)]	= 0;
            ensemble_break24C[index2(k,n,no_of_H)]	= 0;
            ensemble_prolongedC[index2(k,n,no_of_H)]= 0;
			
            ensemble_Confirmed[index2(k,n,no_of_H)] = 0;
            ensemble_Unconfirmed[index2(k,n,no_of_H)] = 0;
            ensemble_ConfirmedS[index2(k,n,no_of_H)] = 0;
            ensemble_UnconfirmedS[index2(k,n,no_of_H)] = 0;
            
            //ensemble_Slaughter[index2(k,n,no_of_H)] = 0;
            ensemble_ReactorsU[index2(k,n,no_of_H)] = 0;
            ensemble_ReactorsC[index2(k,n,no_of_H)] = 0;
			
            ensemble_BreaklengthsU[index2(k,n,no_of_H)] = 0;
            ensemble_BreaklengthsC[index2(k,n,no_of_H)] = 0;
			
            num_confirmed[index2(k,n,no_of_H)] = 0;
            num_unconfirmed[index2(k,n,no_of_H)] = 0;
            
            ensemble_VisitsU[index2(k,n,no_of_H)] = 0;
            ensemble_VisitsC[index2(k,n,no_of_H)] = 0;
            
            ensemble_TestsU[index2(k,n,no_of_H)] = 0;
            ensemble_TestsC[index2(k,n,no_of_H)] = 0;
            
            ensemble_ReactorsTotU[index2(k,n,no_of_H)] = 0;
            ensemble_ReactorsTotC[index2(k,n,no_of_H)] = 0;
            
            ensemble_BurdenU[index2(k,n,no_of_H)] = 0;
            ensemble_BurdenC[index2(k,n,no_of_H)] = 0;
            
            ensemble_Forward_TransU[index2(k,n,no_of_H)] = 0;
            ensemble_Forward_TransC[index2(k,n,no_of_H)] = 0;
            
            ensemble_Onward_TransU[index2(k,n,no_of_H)] = 0;
            ensemble_Onward_TransC[index2(k,n,no_of_H)] = 0;
            
            /*
             sum_Visits[index2(k,n,no_of_H)] = 0;
             sum_Tests[index2(k,n,no_of_H)] = 0;
             sum_Reactors_at_Start[index2(k,n,no_of_H)] = 0;
             sum_Burden[index2(k,n,no_of_H)] = 0;
             sum_Forward_Trans[index2(k,n,no_of_H)] = 0;
             sum_Onward_Trans[index2(k,n,no_of_H)] = 0;
             */
            
		}
        
        
	}
	
	
	
	int Herd_Size = 0;
	
    model.set_susceptible_risk(1.0,rel1,rel2,rel3,rel4);
//model.set_assortative(rel1,rel2,rel3,rel3,rel3,rel4);

    model.set_test_characteristics(Standard,Severe,Standard, Standard, Standard,s_slaughter);

	for(int r=0; r < model.DarthHerdNo; r++)
	//for(int r=2366; r == 2366; r++)
	//for(int i=0;i<500;i++)
	//int r = 2366;
	{
	//int r = 2366;		
						model.DarthSelecta = r;
	
						model.set_forcing(false);
						
							
						int n = model.DarthBindex[r];
						int k = model.DarthPTI[r];
						//cout << "PTI: " << k << ' ' << PTI[k] << endl;
						//Herd_Size = 30+n*60;
						
						Herd_Size = model.DarthHerdSize[r];
						

						//cout << "Herd Size: " << Herd_Size << ' ' << model.DarthMoves[r] << endl;
						// Exponential Model
							//model.set_demo(LifeExp, 364.0/(Turnover),0.0,0.0*Turnover*Herd_Size/364.0);
						// Empirical Model
						model.set_demo(Darth, 0.0,0.0,0.0);
						
						model.set_transmission(beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_o/364.0,beta_o/364.0,beta_r/364.0,beta_i/364.0,
												beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_o/364.0,beta_o/364.0,beta_r/364.0,beta_i/364.0,
												q,xinf_PTI[PTI[k]-1],0.0,p_shadow);

						Seeders = 1;
						//cout << "Entrant " << endl;
						model.run_uk_testing(DiscoFlag,PTI[k],Herd_Size,Seeders,must_clear,clear_from_severe,change_SIT,second_SIT);
						
                         			//model.print_demo();

                        			breakdownFile2 << r << ',' << Herd_Size << ','  << PTI[k] << ',' 
						<< model.primary_breaklength << ',' 
						<< model.break_recurr  << ',' 
						<< model.onward_trans << ',' 
						<< model.Reactors_at_Start << ',' 
						<< model.Reactors_at_VE6M << ','
						<< model.Reactors_at_VE12M << ','
						<< model.Primary_Reactors << ',' 
						<< model.Total_Visits << ',' 
						<< model.Total_Tests  << ',' 
						<< model.Burden << ',' 
						<< model.forward_trans << ',' 
						<< model.slaughter_house << ',' 
						<< model.confirmed_ever
						<<  endl;
						
						
						if(model.confirmed_ever)
						{
                            //cout << primary_breaklength << ' ' << model.return_time() << ' ' << breakstart << ' ' << break6 << ' ' << break12 << ' ' << break24 << endl;
                            
							ensemble_break6C[index2(k,n,no_of_H)]  += (int) model.break6;
							ensemble_break12C[index2(k,n,no_of_H)] += (int) model.break12;
							ensemble_break24C[index2(k,n,no_of_H)] += (int) model.break24;
							ensemble_prolongedC[index2(k,n,no_of_H)] += (int) (model.primary_breaklength > 240);
							ensemble_BreaklengthsC[index2(k,n,no_of_H)]  += model.primary_breaklength;
							if(model.Reactors_at_Start != 0)
							{
                                ensemble_ReactorsC[index2(k,n,no_of_H)] += (model.Reactors_at_Start);
							}
							if(model.Primary_Reactors != 0)
							{
                                ensemble_ReactorsTotC[index2(k,n,no_of_H)] += (model.Primary_Reactors);
							}
							ensemble_VisitsC[index2(k,n,no_of_H)] += (int) model.Total_Visits;
							ensemble_TestsC[index2(k,n,no_of_H)] += (int) model.Total_Tests;
							if(model.Burden > 0)
							{
                                ensemble_BurdenC[index2(k,n,no_of_H)]++;
							}
							
							ensemble_Forward_TransC[index2(k,n,no_of_H)] += (int) model.forward_trans;
							ensemble_Onward_TransC[index2(k,n,no_of_H)] += (int) model.onward_trans;
							
							if(model.slaughter_house)
							{
                                ensemble_ConfirmedS[index2(k,n,no_of_H)]++;
							}
							else
							{
                                ensemble_Confirmed[index2(k,n,no_of_H)]++;
							}
							
							//ensemble_Slaughter[index2(k,n,no_of_H)] += (int) slaughter_house;
							//ensemble_Confirmed[index2(k,n,no_of_H)] += (int) confirmed_at_start;
							num_confirmed[index2(k,n,no_of_H)]++;
						}
						else
						{
						    //cout << "U: " <<  primary_breaklength << ' ' << model.return_time() << ' ' << breakstart << ' ' << break6 << ' ' << break12 << ' ' << break24 << endl;
                            
							ensemble_break6U[index2(k,n,no_of_H)]  += (int) model.break6;
							ensemble_break12U[index2(k,n,no_of_H)] += (int) model.break12;
							ensemble_break24U[index2(k,n,no_of_H)] += (int) model.break24;
							ensemble_prolongedU[index2(k,n,no_of_H)] += (int) (model.primary_breaklength > 240);
							ensemble_BreaklengthsU[index2(k,n,no_of_H)]  += model.primary_breaklength;
							ensemble_ReactorsU[index2(k,n,no_of_H)] += (model.Reactors_at_Start);
							
							ensemble_ReactorsTotU[index2(k,n,no_of_H)] += (model.Primary_Reactors);
							ensemble_VisitsU[index2(k,n,no_of_H)] += (int) model.Total_Visits;
							ensemble_TestsU[index2(k,n,no_of_H)] += (int) model.Total_Tests;
                            
							if(model.Burden > 0)
							{
                                ensemble_BurdenU[index2(k,n,no_of_H)]++;
							}
							ensemble_Forward_TransU[index2(k,n,no_of_H)] += (int) model.forward_trans;
							ensemble_Onward_TransU[index2(k,n,no_of_H)] += (int) model.onward_trans;
							
							if(model.slaughter_house)
							{
                                ensemble_UnconfirmedS[index2(k,n,no_of_H)]++;
							}
							else
							{
                                ensemble_Unconfirmed[index2(k,n,no_of_H)]++;
							}
							
							num_unconfirmed[index2(k,n,no_of_H)]++;
						}
						/*
                         sum_Reactors[index2(k,n,no_of_H)] += Total_Reactors;
                         sum_Breaklengths[index2(k,n,no_of_H)] += primary_breaklength;
                         sum_Visits[index2(k,n,no_of_H)] += Total_Visits;
                         sum_Tests[index2(k,n,no_of_H)] += Total_Tests;
                         sum_Reactors_at_Start[index2(k,n,no_of_H)] += Reactors_at_Start;
                         sum_Burden[index2(k,n,no_of_H)] += Burden;
                         sum_Forward_Trans[index2(k,n,no_of_H)] += forward_trans;
                         sum_Onward_Trans[index2(k,n,no_of_H)] += onward_trans;
                         */
						//}
						continue;
					}
					//Fadeout out of disease before breakdown with p=0 chance of re-introduction
					//if(breakfirstFlag && (xinf == 0.0) && policy_check[4] == 0)
					//Fadeout of disease before breakdown with low chance of re-introduction
					//if(breakfirst==0 && policy_check[4] == 0)
					/*
					 if(first_break && policy_check[4] == 0 && model.return_time() > 3*364.0)
					 {
					 // Drop r by one and loop again
					 cout << "Dropped " << thumb << endl;
					 thumb++;
					 r--;
					 break;
					 }
					 */	
				
				
				
				
				
				
			
			
			
		
	
	// Output ensemble variables
	
	for(int k=0; k < no_of_PTI; k++)
	{
		for(int n=0; n < no_of_H; n++)
		{
            
            breakdownFile	<< 30+n*60 << ','  << PTI[k] << ',' 
            << ensemble_prolongedC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ',' 
            << "NA" << ',' 
            << ensemble_break6C[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ',' 
            << ensemble_break12C[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ',' 
            << ensemble_break24C[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ',' 
            << "NA" << ',' 
            << ((ensemble_ReactorsC[index2(k,n,no_of_H)]/(double) (ensemble_Confirmed[index2(k,n,no_of_H)]))) << ',' 
            << "NA" << ','
            << (ensemble_ConfirmedS[index2(k,n,no_of_H)]) /(double) (num_confirmed[index2(k,n,no_of_H)]) << ',' 
            << num_confirmed[index2(k,n,no_of_H)]/(double) (num_confirmed[index2(k,n,no_of_H)] + num_unconfirmed[index2(k,n,no_of_H)]) << ','
            << 1 << ',' 
            << ensemble_BreaklengthsC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ',' 
            << ((ensemble_ReactorsTotC[index2(k,n,no_of_H)]/(double) ensemble_Confirmed[index2(k,n,no_of_H)])) << ','
            << ensemble_VisitsC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ','
            << ensemble_TestsC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ','
            << ensemble_BurdenC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)] << ','
            << ensemble_Onward_TransC[index2(k,n,no_of_H)]/(double) num_confirmed[index2(k,n,no_of_H)]
            << endl;
            
            
            breakdownFile	<< 30+n*60 << ','  << PTI[k] << ',' 
            << ensemble_prolongedU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ',' 
            << "NA" << ',' 
            << ensemble_break6U[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ',' 
            << ensemble_break12U[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ',' 
            << ensemble_break24U[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ','
            << "NA" << ',' 
            << (ensemble_ReactorsU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)]) << ',' 
            << "NA" << ','
            << (ensemble_UnconfirmedS[index2(k,n,no_of_H)])/(double) (num_confirmed[index2(k,n,no_of_H)]) << ',' 
            << num_confirmed[index2(k,n,no_of_H)]/(double) (num_confirmed[index2(k,n,no_of_H)] + num_unconfirmed[index2(k,n,no_of_H)]) << ','
            << 0 << ',' 
            << ensemble_BreaklengthsU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ',' 
            << (ensemble_ReactorsTotU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)]) << ','
            << ensemble_VisitsU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ','
            << ensemble_TestsU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ','
            << ensemble_BurdenU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)] << ','
            << ensemble_Onward_TransU[index2(k,n,no_of_H)]/(double) num_unconfirmed[index2(k,n,no_of_H)]
            << endl;
            
            
		}
        
        
	}
    
	
	
	breakdownFile.close();
	breakdownFile2.close();
	
	FILE * ReactorsFile;
	FILE * CReactorsFile;
	FILE * SeactorsFile;
	FILE * SlaughterFile;
	
	
	ReactorsFile = fopen("ReactorsDistro.dat","w");
	CReactorsFile = fopen("CReactorsDistro.dat","w");
	SeactorsFile = fopen("SReactorsDistro.dat","w");
	SlaughterFile = fopen("SlaughterDistro.dat","w");
	
	
	gsl_histogram_fprintf(ReactorsFile, model.AgeReactors, "%g", "%g");
	gsl_histogram_fprintf(CReactorsFile, model.AgeCReactors, "%g", "%g");
	gsl_histogram_fprintf(SeactorsFile, model.AgeSReactors, "%g", "%g");
	gsl_histogram_fprintf(SlaughterFile, model.AgeSlaughter, "%g", "%g");
	
	fclose(ReactorsFile);
	fclose(CReactorsFile);
	fclose(SeactorsFile);
	
	
	return 0;
}

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
    
	bool DiscoFlag = false;
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
	double Severe[2] = {0.0,0.0};
	double Standard_v[2] = {0.0,0.0}; 
	double Severe_v[2] = {0.0,0.0};
	double DIVA[2] = {0.0,0.0};

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
	
	double vacc_eff = 1.0;
	double vacc_lay = 1.0;
	bool doDIVA = false;

	double test_interval;
    int herd_size;

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
	
	// Vaccination Strategy
	double batch_interval = 0.0;
	bool do_vaccination = false;
	bool do_batching = false;
	bool suspend_vacc_in_break = false;
	double Vacc_Eligible_Age = 0.0;
	double vacc_target = 0.0;					
	double trial_duration = 5.0;
	bool retain_reactors=false;
	int this_herd = 0;
// Extra 6	
	
if(argc < 46)
    { 
		cout << "bTBtrial <Sensitivity Standard> <Specificity Standard> <Sensitivity Severe> <Specificity Severe> <Vaccinate Specificity Standard> <Vaccinate Specificity Severe> <DIVA Sensitivity> <DIVA Specificity>  <Vacc Eff> <Vacc Eff Latency> <doDiva> <O> <R> <V1> <V2> <OV1> <OV2> <RV> <beta_o> <beta_r> <beta_i> <beta_ov1> <beta_ov2> <beta_rv> <beta_iv> <q> <Xinf PTI1> <Xinf PTI2> <Xinf PTI4> <rel1> <rel2> <rel3> <rel4> <p_confirm_I> <p_confirm_R> <p_confirm_O> <p_slaughter> <delta_age> <p_shadow> <shadow_delay> <no_of_reps> <outputfile> <trial_duration> <test_interval> <herd_size> " << endl;
		exit(EXIT_FAILURE);
    }
    else
    {
		Standard[0] = atof(argv[1]);
		Standard[1] = atof(argv[2]);
		Severe[0] = atof(argv[3]);
		Severe[1] = atof(argv[4]);
		Standard_v[0] = atof(argv[1]);
		Standard_v[1] = atof(argv[5]);
		Severe_v[0] = atof(argv[3]);
		Severe_v[1] = atof(argv[6]);
		DIVA[0] = atof(argv[7]);
		DIVA[1] = atof(argv[8]);
		
		vacc_eff = atof(argv[9]);
		vacc_lay = atof(argv[10]);
		
		if(atoi(argv[11])==0)
		{
		 doDIVA = false;
		}
		else
		{
		 doDIVA = true;
		}
		
		TO = atof(argv[12]);
		TR = atof(argv[13]);
		TV1 = atof(argv[14]);
		TV2 = atof(argv[15]);
		TOV1 = atof(argv[16]);
		TOV2 = atof(argv[17]);
		TRV = atof(argv[18]);
		beta_o = atof(argv[19]);
		beta_r = atof(argv[20]);
		beta_i = atof(argv[21]);
		beta_ov1 = atof(argv[22]);
		beta_ov2 = atof(argv[23]);
		beta_rv = atof(argv[24]);
		beta_iv = atof(argv[25]);
		
		q      = atof(argv[26]);
		
		xinf_PTI[0] = atof(argv[27]);
		
		xinf_PTI[1] = atof(argv[28]);
        
		xinf_PTI[3] = atof(argv[29]);
		
		xinf_PTI[2] = 0.0;
		rel1        = atof(argv[30]);
		rel2        = atof(argv[31]);
		rel3        = atof(argv[32]);
		rel4        = atof(argv[33]);

		p_confirm[0] = atof(argv[34]);
		
		p_confirm[1] = atof(argv[35]);
		
		p_confirm[2] = atof(argv[36]);
		
		
		s_slaughter = atof(argv[37]);
		
		delta_age = atof(argv[38]);
		
		p_shadow  = atof(argv[39]);
		
		shad_lay = atof(argv[40]);
		
		
		no_of_reps = atoi(argv[41]);
		
		filename = argv[42];
		
		trial_duration = atof(argv[43]);
		test_interval = atoi(argv[44]);
		herd_size = atoi(argv[45]);
		
		low_PTI = 0;
		high_PTI = 2;
			
		low_n = 30;
		high_n = 330;
			
			
		}
		/*
		cout << "Sensitivity: " << Standard[0] << endl;
		cout << "Specificity: " << Standard[1] << endl;
		cout << "Sensitivity (Severe): " << Severe[0] << endl;
		cout << "Specificity (Severe): " << Severe[1] << endl;
		
		cout << "Sensitivity_v: " << Standard_v[0] << endl;
		cout << "Specificity_v: " << Standard_v[1] << endl;
		cout << "Sensitivity_v (Severe): " << Severe_v[0] << endl;
		cout << "Specificity_v (Severe): " << Severe_v[1] << endl;
		cout << "Sensitivity DIVA: " << DIVA[0] << endl;
		cout << "Specificity DIVA: " << DIVA[1] << endl;
		
		cout << "Vaccine Efficacy: " << vacc_eff << endl;
		cout << "Vaccine D'Lay: " << vacc_lay << endl;
		cout << "doDIVA: " << doDIVA << endl;
		
		//cout << "Repeat Sensitivity: " << Repeat_Sensitivity[0] << endl;
		//cout << "Repeat Sensitivity (Severe): " << Repeat_Sensitivity[1] << endl;
        
		
		cout << "TO: " << TO << endl;
		cout << "TR: " << TR << endl;
		cout << "TV1: " << TV1 << endl;
		cout << "TV2: " << TV2 << endl;
		cout << "TOV1: " << TOV1 << endl;
		cout << "TOV2: " << TOV2 << endl;
		cout << "TRV: " << TRV << endl;
				
					
		cout << "beta_o: " << beta_o << endl;
		cout << "beta_r: " << beta_r << endl;
		cout << "beta_i: " << beta_i << endl;
		cout << "beta_ov1: " << beta_ov1 << endl;
		cout << "beta_ov2: " << beta_ov2 << endl;
		cout << "beta_rv: " << beta_rv << endl;
		cout << "beta_iv: " << beta_iv << endl;
		
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
						
		cout << "No. of Reps: " << no_of_reps << endl;
		
		cout << low_n << ' ' << high_n << endl;
		cout << low_PTI << ' ' << high_PTI << endl;

		if(do_vaccination)
		{
		 cout << "Do Vaccination, eligible age: " << Vacc_Eligible_Age << endl;
		}
		if(do_batching)
		{
		cout << "Batch vaccination at interval of : " << batch_interval << endl;
		if(batch_interval < 7.0)
		{cout << "Do not use batch_interval of less than 1 week" << endl; exit(0);}
		}
		if(suspend_vacc_in_break)
		{
		cout << "Suspend vaccination during Breakdown" << endl;
		}
		cout << "Vaccination Target " << vacc_target << endl;
		cout << "Trial Duration " << trial_duration << endl;
		if(retain_reactors){cout << "Retain Reactors." << endl;}
	*/		  
    
	    

	//TOF = 0.1;
	//TR = 0.1;
	//beta_r = 0.0;
	//beta_i = 0.0;
	//q = 1.0;

		bTBICBM phase1g1(1,
        0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,
        PATH_TO_DATA,filename,save_ts,
        TO*364.0,TR*364.0,TV1*364.0,TV2*364.0,TOV1*364.0*vacc_lay,TOV2*364.0*vacc_lay,TRV*364.0*vacc_lay,
        TO*364.0*shad_lay,364.0*shad_lay,TV1*364.0,TV2*364.0,TOV1*364.0*vacc_lay,TOV2*vacc_lay,364.0*shad_lay*vacc_lay,
        false,0.0, delta_age);
        
        bTBICBM phase1g2(1,
        0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,
        PATH_TO_DATA,filename,save_ts,
        TO*364.0,TR*364.0,TV1*364.0,TV2*364.0,TOV1*364.0*vacc_lay,TOV2*364.0*vacc_lay,TRV*364.0*vacc_lay,
        TO*364.0*shad_lay,364.0*shad_lay,TV1*364.0,TV2*364.0,TOV1*364.0*vacc_lay,TOV2*vacc_lay,364.0*shad_lay*vacc_lay,
        false,0.0, delta_age);
        
	phase1g1.set_confirmation(p_confirm);
	phase1g2.set_confirmation(p_confirm);
	//phase1.set_transmission(0.0,beta_r,beta_i,0.0,0.0,beta_r,beta_i,q,xinf_PTI[0]);
	//phase1.set_transmission(0.0,beta_r,beta_i,0.0,0.0,beta_r,beta_i,q,xinf_PTI[0]);
	//phase1.set_transmission(beta_o,beta_i,beta_r,beta_ov1,beta_ov2,beta_rv,beta_iv,q,0.0);
	
	phase1g1.set_forcing(false);
	phase1g2.set_forcing(false);
	//phase1.set_demo(LifeExp, 50.0,1.5356);
	//phase1.set_demo(LifeFixed, 50.0,0,(1.0/50.0)*2000.0);
	//phase1.set_demo(LifeExp, 50.0,0,(1.0/50.0)*2000.0);
	//phase1.set_demo(LifeGamma, 50.0,3.0,(1.0/50.0)*2000.0);
	//phase1.set_demo(LifeExp, 50.0,0.0);
	
	double parish_testing = 0.0;
	
	//const int no_of_N = 24;
	//int N_size[no_of_N] = {10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400};
	
	//const int no_of_N = 12;
	//int N_size[no_of_N] = {10,20,30,40,50,100,150,200,250,300,350,400};
	
		
     phase1g1.set_susceptible_risk(1.0,rel1,rel2,rel3,rel4);
     phase1g2.set_susceptible_risk(1.0,rel1,rel2,rel3,rel4);

    //phase1.set_test_characteristics(Standard,Severe,Standard_v, Severe_v, DIVA,s_slaughter,false,doDIVA);
	
	//phase1.set_test_characteristics(Standard,Severe,Standard_v, Severe_v, DIVA,s_slaughter,false,false);
	//phase1.set_test_characteristics(Standard,Severe,Standard_v, Severe_v, DIVA,s_slaughter);
	
	phase1g1.set_test_characteristics(Standard,Severe,Standard_v, Severe_v, DIVA,s_slaughter,false,doDIVA);
	phase1g2.set_test_characteristics(Standard,Severe,Standard_v, Severe_v, DIVA,s_slaughter,false,doDIVA);

	phase1g1.DarthSelecta = 0;
	phase1g2.DarthSelecta = 0;
					
	phase1g1.set_demo(Experiment, 0.0,0.0,0.0);
	phase1g2.set_demo(Experiment, 0.0,0.0,0.0);					
	phase1g1.set_transmission(beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_ov1/364.0,beta_ov2/364.0,beta_rv/364.0,beta_iv/364.0,
												beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_ov1/364.0,beta_ov2/364.0,beta_rv/364.0,beta_iv/364.0,
												q,0.0,vacc_eff,p_shadow);
												phase1g2.set_transmission(beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_ov1/364.0,beta_ov2/364.0,beta_rv/364.0,beta_iv/364.0,
												beta_o/364.0,beta_r/364.0,beta_i/364.0,beta_ov1/364.0,beta_ov2/364.0,beta_rv/364.0,beta_iv/364.0,
												q,0.0,vacc_eff,p_shadow);
	int group_id = 0;
	int phase_id = 1;
	vector<int> policy_check;
	
	int Utot,Rtot,Wtot,Vtot;
	
	Seeders = herd_size/2;
	
	if((herd_size % 4)!=0){cout << "Group Size must be divisible by 4" << endl;}
	
	for(int r=0; r < no_of_reps; r++)
	{
	
	group_id = 0;
	phase_id = 1;
	// Start with Clean Herd // Vaccinate once at time zero
	phase1g1.initialise_herd_Experiment(herd_size);
	phase1g1.vaccinate_off_schedule();
	
	Rtot = Seeders;
	Utot = Seeders/2;
	Vtot = Seeders/2;
	Wtot = 0;
	
	
	policy_check = phase1g1.ExperimentWithDIVASinglePhase(r, group_id, phase_id,Seeders,Rtot,Utot,Vtot,Wtot, trial_duration,test_interval,false);
	
	group_id = 1;
	phase_id = 1;
	
	phase1g2.initialise_herd_Experiment(herd_size);
	phase1g2.vaccinate_off_schedule();
	
	policy_check = phase1g2.ExperimentWithDIVASinglePhase(r, group_id, phase_id,Seeders,Rtot,Utot,Vtot,Wtot, trial_duration,test_interval,false);
	
	
	
	phase1g1.RemoveSeeders();
	

	
	phase1g2.RemoveSeeders();
	
	phase_id = 2;
	group_id = 2;
	
	//cout << "Entering the funky constructor" << endl;
	
	// Vaccinated Seeders
	bTBICBM phase2g2 = bTBICBM (phase1g1, phase1g2, false);
	
	phase2g2.refresh_herd_Experiment(herd_size);
	phase2g2.vaccinate_off_schedule();
	
	Rtot = 0;
	Utot = Seeders/2;
	Vtot = Seeders/2;
	Wtot = Seeders;
	
	policy_check = phase2g2.ExperimentWithDIVASinglePhase(r, group_id, phase_id,Seeders,Rtot,Utot,Vtot,Wtot, trial_duration,test_interval,false);
	
	phase_id = 2;
	group_id = 4;
	
	Rtot = Seeders;
	Utot = Seeders/2;
	Vtot = Seeders/2;
	Wtot = 0;
	
		//cout << "Entering the second funky constructor" << endl;
	
	
	// Unvacinated Seeders	
	bTBICBM phase2g4 = bTBICBM (phase1g1, phase1g2, true); 

	phase2g4.refresh_herd_Experiment(herd_size);
	phase2g4.vaccinate_off_schedule();

policy_check = phase2g4.ExperimentWithDIVASinglePhase(r, group_id, phase_id,Seeders,Rtot,Utot,Vtot,Wtot, trial_duration,test_interval,false);

	}
	
	return 0;
}

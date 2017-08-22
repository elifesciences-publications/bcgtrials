#include <iostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include <math.h>
#include <algorithm>
#include <map>
using namespace std;


#include<time.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

#include "./mersennetwister.h"


// SORID model
const double Rstate_emp[4] = {0.0701609761,0.0001441427,0.8060078782,0.1236870030}; 

enum record {nosave,overwrite,append};

typedef record record_t;

enum cow_demo {Birth_demo, Death_demo, Off_demo};

enum slaughter_state {susceptible, occult, reactive};

enum test {wht,rht,all};

typedef test test_t;

enum epi_status_t {S,O,R,I,V1,V2,OV1,OV2,RV,IV,
				   SS,SO,SR,SI,SV1,SV2,SOV1,SOV2,SRV,SIV};

enum epi_events {Birth,Strans,V1trans,V2trans,Omove,Rmove,V1move,V2move,OV1move,Ov2move,RVmove,
				 SStrans,SV1trans,SV2trans,SOmove,SRmove,SV1move,SV2move,SOV1move,SOv2move,SRVmove};

enum test_status_t {NotReactor,Reactor};

enum demo_t {LifeExp,LifeFixed,LifeNegBin,Darth,Experiment};

enum move_t {uniform, radiation, gravity};

enum forcing_t {seasonal,constant};

const int SICCT = 0;
const int SICCT_S = 1;
const int SICCT_V1 = 2;
const int SICCT_V2 = 3;
const int DIVA = 4; 

struct cow_t {
	double Birth_time;
    double Death_time;
    double Off_time;
    double Infection_time;
    double Vaccinated_time;
	epi_status_t Epi_status;
	int Epi_stage;
	test_status_t Std_status;
	test_status_t Svr_status;
	test_status_t Diva_status;
	bool Diva_ever;	
	bool Confirmation_status;
	bool Control;
	bool Seeder;
	double rate;
	double uniq_id;
};

struct vaccination_t {
	double time;
	int herd;
	bool revaccinate;
};


bool compare_birthdate(cow_t one,cow_t two)
{
    return (one.Birth_time < two.Birth_time);
}

bool compare_firstOff(cow_t one,cow_t two)
{
	double one_off = (one.Death_time < one.Off_time) ? one.Death_time : one.Off_time;
	double two_off = (two.Death_time < two.Off_time) ? two.Death_time : two.Off_time;
    return (one_off < two_off);
}

bool compare_Off_time(cow_t one,cow_t two)
{
    return (one.Off_time <= two.Off_time);
}

const int yeardef = 364;

const double birth_seas[yeardef]={0.00170018407151270,0.00170018407151270,0.00170018407151270,0.00170018407151270,0.00170018407151270,0.00170018407151270,0.00170018407151270,0.00164370528649813,0.00164370528649813,0.00164370528649813,0.00164370528649813,0.00164370528649813,0.00164370528649813,0.00164370528649813,0.00164074569951047,0.00164074569951047,0.00164074569951047,0.00164074569951047,0.00164074569951047,0.00164074569951047,0.00164074569951047,0.00157242856654525,0.00157242856654525,0.00157242856654525,0.00157242856654525,0.00157242856654525,0.00157242856654525,0.00157242856654525,0.00167275856542702,0.00167275856542702,0.00167275856542702,0.00167275856542702,0.00167275856542702,0.00167275856542702,0.00167275856542702,0.00180105666134221,0.00180105666134221,0.00180105666134221,0.00180105666134221,0.00180105666134221,0.00180105666134221,0.00180105666134221,0.00200852370917737,0.00200852370917737,0.00200852370917737,0.00200852370917737,0.00200852370917737,0.00200852370917737,0.00200852370917737,0.00224514268884101,0.00224514268884101,0.00224514268884101,0.00224514268884101,0.00224514268884101,0.00224514268884101,0.00224514268884101,0.00273786459583707,0.00273786459583707,0.00273786459583707,0.00273786459583707,0.00273786459583707,0.00273786459583707,0.00273786459583707,0.00302538847168852,0.00302538847168852,0.00302538847168852,0.00302538847168852,0.00302538847168852,0.00302538847168852,0.00302538847168852,0.00353700240895583,0.00353700240895583,0.00353700240895583,0.00353700240895583,0.00353700240895583,0.00353700240895583,0.00353700240895583,0.00399268015215631,0.00399268015215631,0.00399268015215631,0.00399268015215631,0.00399268015215631,0.00399268015215631,0.00399268015215631,0.0043465481029812,0.0043465481029812,0.0043465481029812,0.0043465481029812,0.0043465481029812,0.0043465481029812,0.0043465481029812,0.00459574532734241,0.00459574532734241,0.00459574532734241,0.00459574532734241,0.00459574532734241,0.00459574532734241,0.00459574532734241,0.00454252208801427,0.00454252208801427,0.00454252208801427,0.00454252208801427,0.00454252208801427,0.00454252208801427,0.00454252208801427,0.00447050547131447,0.00447050547131447,0.00447050547131447,0.00447050547131447,0.00447050547131447,0.00447050547131447,0.00447050547131447,0.00422323197849525,0.00422323197849525,0.00422323197849525,0.00422323197849525,0.00422323197849525,0.00422323197849525,0.00422323197849525,0.00416537205288644,0.00416537205288644,0.00416537205288644,0.00416537205288644,0.00416537205288644,0.00416537205288644,0.00416537205288644,0.00377445993826598,0.00377445993826598,0.00377445993826598,0.00377445993826598,0.00377445993826598,0.00377445993826598,0.00377445993826598,0.00377884999229767,0.00377884999229767,0.00377884999229767,0.00377884999229767,0.00377884999229767,0.00377884999229767,0.00377884999229767,0.00335089371388163,0.00335089371388163,0.00335089371388163,0.00335089371388163,0.00335089371388163,0.00335089371388163,0.00335089371388163,0.00310623452290151,0.00310623452290151,0.00310623452290151,0.00310623452290151,0.00310623452290151,0.00310623452290151,0.00310623452290151,0.0028526472445086,0.0028526472445086,0.0028526472445086,0.0028526472445086,0.0028526472445086,0.0028526472445086,0.0028526472445086,0.00272158686740493,0.00272158686740493,0.00272158686740493,0.00272158686740493,0.00272158686740493,0.00272158686740493,0.00272158686740493,0.00271275743289174,0.00271275743289174,0.00271275743289174,0.00271275743289174,0.00271275743289174,0.00271275743289174,0.00271275743289174,0.00276410626712768,0.00276410626712768,0.00276410626712768,0.00276410626712768,0.00276410626712768,0.00276410626712768,0.00276410626712768,0.00259649499072638,0.00259649499072638,0.00259649499072638,0.00259649499072638,0.00259649499072638,0.00259649499072638,0.00259649499072638,0.00261745873188899,0.00261745873188899,0.00261745873188899,0.00261745873188899,0.00261745873188899,0.00261745873188899,0.00261745873188899,0.00260808670642806,0.00260808670642806,0.00260808670642806,0.00260808670642806,0.00260808670642806,0.00260808670642806,0.00260808670642806,0.00249429058675242,0.00249429058675242,0.00249429058675242,0.00249429058675242,0.00249429058675242,0.00249429058675242,0.00249429058675242,0.00246499067557456,0.00246499067557456,0.00246499067557456,0.00246499067557456,0.00246499067557456,0.00246499067557456,0.00246499067557456,0.00251239339382696,0.00251239339382696,0.00251239339382696,0.00251239339382696,0.00251239339382696,0.00251239339382696,0.00251239339382696,0.00254337040429783,0.00254337040429783,0.00254337040429783,0.00254337040429783,0.00254337040429783,0.00254337040429783,0.00254337040429783,0.00273673008749180,0.00273673008749180,0.00273673008749180,0.00273673008749180,0.00273673008749180,0.00273673008749180,0.00273673008749180,0.00279478731889979,0.00279478731889979,0.00279478731889979,0.00279478731889979,0.00279478731889979,0.00279478731889979,0.00279478731889979,0.00278314627674832,0.00278314627674832,0.00278314627674832,0.00278314627674832,0.00278314627674832,0.00278314627674832,0.00278314627674832,0.00290345348779681,0.00290345348779681,0.00290345348779681,0.00290345348779681,0.00290345348779681,0.00290345348779681,0.00290345348779681,0.00283903314436535,0.00283903314436535,0.00283903314436535,0.00283903314436535,0.00283903314436535,0.00283903314436535,0.00283903314436535,0.00268227368691881,0.00268227368691881,0.00268227368691881,0.00268227368691881,0.00268227368691881,0.00268227368691881,0.00268227368691881,0.00266466414434221,0.00266466414434221,0.00266466414434221,0.00266466414434221,0.00266466414434221,0.00266466414434221,0.00266466414434221,0.00250016043427796,0.00250016043427796,0.00250016043427796,0.00250016043427796,0.00250016043427796,0.00250016043427796,0.00250016043427796,0.00253498490783279,0.00253498490783279,0.00253498490783279,0.00253498490783279,0.00253498490783279,0.00253498490783279,0.00253498490783279,0.00247559586228035,0.00247559586228035,0.00247559586228035,0.00247559586228035,0.00247559586228035,0.00247559586228035,0.00247559586228035,0.00246415212592806,0.00246415212592806,0.00246415212592806,0.00246415212592806,0.00246415212592806,0.00246415212592806,0.00246415212592806,0.00239711748065750,0.00239711748065750,0.00239711748065750,0.00239711748065750,0.00239711748065750,0.00239711748065750,0.00239711748065750,0.00231513692109924,0.00231513692109924,0.00231513692109924,0.00231513692109924,0.00231513692109924,0.00231513692109924,0.00231513692109924,0.00222570806762202,0.00222570806762202,0.00222570806762202,0.00222570806762202,0.00222570806762202,0.00222570806762202,0.00222570806762202,0.00213835092503618,0.00213835092503618,0.00213835092503618,0.00213835092503618,0.00213835092503618,0.00213835092503618,0.00213835092503618,0.00190419826786225,0.00190419826786225,0.00190419826786225,0.00190419826786225,0.00190419826786225,0.00190419826786225,0.00190419826786225,0.00187277731934323,0.00187277731934323,0.00187277731934323,0.00187277731934323,0.00187277731934323,0.00187277731934323,0.00187277731934323,0.00189729256489104,0.00189729256489104,0.00189729256489104,0.00189729256489104,0.00189729256489104,0.00189729256489104,0.00189729256489104,0.00190755246644827,0.00190755246644827,0.00190755246644827,0.00190755246644827,0.00190755246644827,0.00190755246644827,0.00190755246644827};

const double confirm_AgeBreaks[21] = {200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,
							  3000, 3200, 3400, 3600, 3800, 4000, 8000};

const double confirm_by_age[21] = {0.5757576,0.5892473,0.5729814,0.4987841,0.4068989,0.3093842,0.2898807,
							  0.2675159,0.2597403,0.2579387,0.2811751,0.2803407,0.2733700,0.2766355,
							  0.3201804,0.2732026,0.3055163,0.2844037,0.2743764,0.2849462,0.2751220};

// Sensitivity of vaccinates
// As function of time from vaccination (days)

const int dim_vaccSpec = 7;

const double vaccStime[dim_vaccSpec] ={3.0*30.0,6.0*30.0,9.0*30.0,12.0*30.0,15.0*30.0,24.0*30.0,10*12*30};
const double vaccSpec[dim_vaccSpec] = {0.0,0.60,0.80,0.09,0.05,0.11,0.10};
const double vaccSpecS[dim_vaccSpec] = {0.0,0.84,0.95,0.30,0.10,0.30,0.30};

int number_testWHT(gsl_rng *r, int Herd_Size, double u)
{
	// Test all animals except those under the age of 6 weeks (i.e. practically all!)
	return(gsl_ran_binomial(r,(1-(6.0/52.0)*364.0*u/Herd_Size),Herd_Size));
}

// Sample from empirical distribution RHT

int  number_testRHT(gsl_rng *r, int Herd_Size)
{
	
	//double pval = RHT[(int) round(gsl_ran_flat(gsl_r,0,no_of_RHTs))];
	
	//return(gsl_ran_binomial(gsl_r,pval,N[i]));
	
	double pval=0.0;
	do
	{
		pval = 0.4941136064 + gsl_ran_cauchy(r,0.0932207681);
	}while(pval < 0 || pval > 1.0);
	
	return(gsl_ran_binomial(r,pval,Herd_Size));
	
}

int animals_to_test(gsl_rng *r, int Herd_Size, double u, int PTI)
{
	int animals_to_test = 0;
	
	test_t whole_herd_test;
	
	switch(PTI)
	{
			
			// Choose Routine tests to be WHT/RHT acording to proportions that initiate breakdowns in 
                        // Study Population 
                        // PTI 1: 52 RHT, 3235 WHT, SLH 781 p(WHT) = 3235/(3235+52) = 0.9841801
                    // PTI 2: 443 RHT, 478 WHT, SLH 282 p(WHT) = 443/(478+443) = 0.48
                        // PTI 4: 555 RHT, 111 WHT, SLH 268 p(WHT) = 111/(555+111) = 0.1666667
			// PTI 1
			case 1:
			if(gsl_ran_flat(r,0,1) <= 0.9841801){whole_herd_test = wht;}
			else{whole_herd_test = rht;}
			break;
			// PTI 2
			case 2:
			if(gsl_ran_flat(r,0,1) <= 0.48){whole_herd_test = wht;}
			else{whole_herd_test = rht;}
			break;
			// PTI 4
			case 4:
			if(gsl_ran_flat(r,0,1) <= 0.1666667){whole_herd_test = wht;}
			else{whole_herd_test = rht;}
			break;
	}
	
	switch(whole_herd_test)
	{
		case wht:
			animals_to_test = number_testWHT(r, Herd_Size, u); 
			break;
		case rht:
			animals_to_test = number_testRHT(r, Herd_Size);
			break;
	}
	
	return(animals_to_test);
}


int ran_zeroPois(gsl_rng *r, int n, double lambda, double p0)
{
	double ans=0;
	ans = gsl_ran_poisson(r,lambda);
	
	if(gsl_rng_uniform(r) < p0){return 0;}else{return(ans);}
}

int* loadEmpirical(int & length,const char filename[])
{
	
	ifstream inputfile(filename);
	
	int cracker=0;
	
	int *SIT;
	
	do{
		inputfile >> cracker;
		length++;
	}while(inputfile);
	length--;
	inputfile.close();
	
	inputfile.open(filename);
	
	SIT = (int *) calloc(sizeof(int),length);
	
	for(int i=0; i < length; i++)
	{
		inputfile >> SIT[i]; 
		
	}
	
	inputfile.close();
	return(SIT);
	
}

double* loadEmpiricalD(int & length,const char filename[])
{
	
	ifstream inputfile(filename);
	
	double cracker=0.0;
	
	double *SIT;
	
	do{
		inputfile >> cracker;
		length++;
	}while(inputfile);
	length--;
	inputfile.close();
	
	inputfile.open(filename);
	
	SIT = (double *) calloc(sizeof(double),length);
	
	for(int i=0; i < length; i++)
	{
		inputfile >> SIT[i]; 
		//cout << SIT[i] << endl;
	}
	
	
	return(SIT);
	
}


double** loadEmpiricalDD(int & length,const char filename[])
{
	
	ifstream inputfile(filename);
	
	double cracker=0.0;
	
	double **SIT;
	
	do{
		inputfile >> cracker >> cracker;
		length++;
	}while(inputfile);
	length--;
	inputfile.close();
	
	inputfile.open(filename);
	
	SIT = (double **) calloc(sizeof(double *),2);
	
	SIT[0] = (double *) calloc(sizeof(double),length);
	SIT[1] = (double *) calloc(sizeof(double),length);
	
 	
	for(int i=0; i < length; i++)
	{
		inputfile >> SIT[0][i] >> SIT[1][i]; 
	}
	
	return(SIT);
	
}


void mkherd(int herd[],int s,int of,int os,int r,int ds, int dof, int dos, int dr, int inf,int dinf)
{
	//cout << "Herd Structure: " << endl;
	
	for(int i =0; i < s; i++)
	{
		
		herd[i] = 0;
		//cout << herd[i] << ' ';
	}
	
	for(int i = s; i < (s+of); i++)
	{
		
		herd[i] = 1;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of); i < (s+of+os); i++)
	{
		
		herd[i] = 2;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os); i < (s+of+os+r); i++)
	{
		
		herd[i] = 3;
		//cout << herd[i] << ' ';
	}
	
	for(int i =(s+of+os+r); i < (s+of+os+r+ds); i++)
	{
		
		herd[i] = 4;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os+r+ds); i < (s+of+os+r+ds+dof); i++)
	{
		
		herd[i] = 5;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os+r+ds+dof); i < (s+of+os+r+ds+dof+dos); i++)
	{
		
		herd[i] = 6;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os+r+ds+dof+dos); i < (s+of+os+r+ds+dof+dos+dr); i++)
	{
		
		herd[i] = 7;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os+r+ds+dof+dos+dr); i < (s+of+os+r+ds+dof+dos+dr+inf); i++)
	{
		
		herd[i] = 8;
		//cout << herd[i] << ' ';
	}
	
	for(int i = (s+of+os+r+ds+dof+dos+dr+inf); i < (s+of+os+r+ds+dof+dos+dr+inf+dinf); i++)
	{
		
		herd[i] = 9;
		//cout << herd[i] << ' ';
	}
	
}


int ran_zanormal(gsl_rng *r, int n, double mu,double v, double p0)
{
	double ans=0;
	// Ensure that we test a positive number of animals less than herd size (n)
	do{
		ans = gsl_ran_gaussian(r,v);
		ans = mu + ans;
	}
	while(ans < 0.0 | ans > n);
	
	if(gsl_rng_uniform(r) < p0){return n;}else{return((int) ans);}
}




// Helper function to add number <no> to filename <s>
inline char *filename(char *s, int no)
{
	static char buf[255];
	sprintf(buf, "%s%d.dat", s, no);
	return buf;
}


// Helper functions for 2,3 matrix indexing
inline int index2(int i, int a, int y)
{
	return y*i + a;
}

inline int index3(int i, int a, int k, int y, int z)
{
	return y*z*i + z*a + k;
}

inline void deref2(int l, int y, int & i, int & a)
{
	a = l % y;
	i = (l - a)/y;
}

inline void deref3(int l, int y, int z, int & i, int & a, int & k)
{
	int m = l % (y*z);
	
	k = m % z;
	a = (m-k)/z;
	i = (l - z*a - k)/(y*z);
}


// gamma function
int gamma(int x)
{
	if(x <= 0){return 0;}
	else {return 1;}
} 


class bTBICBM {
public:
	bTBICBM(int Herds, int OrderO, int OrderR, int OrderV1, int OrderV2, int OrderOV1,int OrderOV2, int OrderRV, int OrderSO, int OrderSR, int OrderSV1, int OrderSV2, int OrderSOV1,int OrderSOV2, int OrderSRV, string input, string output,bool setfull_save,double latent_O, double latent_R, double latent_V1,double latent_V2, double latent_OV1, double latent_OV2, double latent_RV,double latent_SO, double latent_SR, double latent_SV1,double latent_SV2, double latent_SOV1, double latent_SOV2, double latent_SRV, bool verbose,double Setp_inf,double Setdelta_age, bool setdebug_save);
	bTBICBM( const bTBICBM& other);
	bTBICBM(const bTBICBM& other, const bTBICBM& one, bool pick_controls);
	void set_homogenous_mixing();
	void set_susceptible_risk(double beta1,double beta2,double beta3,double beta4,double beta5);
	void set_assortative(double beta1, double beta2, double beta3, double beta4, double beta5,double beta6);
	double run(double runterval,double time_slice,record_t record);
	bool initialise_herd(int Herd_size,int herd_index, int Marys, epi_status_t Mary_Type);
	bool initialise_herd_AllInAllOut(int Herd_size,int herd_index, int Marys, epi_status_t Mary_Type);
	bool initialise_herd_Experiment(int Herd_size);
	bool refresh_herd_Experiment(int Herd_size);
    bool set_demo(demo_t model, double Life_Expectation, double var1, double births);
    bool set_demo(demo_t model, double* Life_Expectation, double* var1, double *births, int Patch);
    void set_forcing(bool force_it);
    void set_inout(string input, string output);
	void set_time(double set_clock);
	void set_confirmation(double *set_conf);
	void set_test_characteristics(double *standard,double *severe,double *vacc1, double *vacc2, double *DIVAt,double slaughter,bool retro_slaughter,bool doDIVA);
	double return_time();
	bool disease_free(int i);
	void set_transmission(double beta_o, double beta_r,double beta_i, double beta_ov1, double beta_ov2,double beta_rv, double beta_iv,double beta_So, double beta_Sr,double beta_Si, double beta_Sov1, double beta_Sov2,double beta_Srv, double beta_Siv, double set_scaling,double xinf_set,double setVacc_Eff,double set_Pshadow);
	vector<int> Test_Protocol_Per_Animal(int this_herd,int PTI,test_t whole_herd_test,double eligible_vacc);
	vector<int> DIVA_Protocol_Per_Animal(int this_herd,int PTI,test_t whole_herd_test,double eligible_vacc);
	void print_demo();
	void Cull(int this_herd, vector<int> results,bool confirmed, bool disclose_tests,int PTI,bool retain);
	void Whole_Herd_Cull(int this_herd, vector<int> results,bool confirmed, bool disclose_tests,int PTI,bool retain);
	void RemoveSeeders();
	void PhaseII(int this_herd, vector<int> results, bool keep, int herd_size, bool control);
	void AltPhaseII(int this_herd, vector<int> results, bool keep, int herd_size, bool control);
	void AltPhaseIIB(int this_herd, vector<int> results, bool keep, int herd_size, bool control);			
	void InvPhaseII(int this_herd, vector<int> results, bool keep, int herd_size,bool control);	
	int accum_infected_cow_days();
	void cull_FP(int Mary);
	int burden(int i);
	void  run_uk_testing(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation);
	void run_until_first_test(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation);	
	void  run_ni_testing(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation);
	void  run_uk_with_Vacc(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible);	
	void  run_uk_with_DIVA(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible);	
	void  run_uk_with_Vacc_endVL(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible);	
	void  run_uk_with_DIVA_endVL(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible);	
	void  run_uk_with_VaccZero(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible);	
	void  run_uk_introduce_Vacc(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc,bool start_when_clear);	
	void  trial_uk_with_DIVA(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible, bool retain,double trial_years,double vacc_p);	
	void  DIVAtrial_uk_with_DIVA(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear,bool clear_from_severe,bool change_SIT,int second_SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible, bool retain,double trial_years,double vacc_p);	
	void  trial_fixed_with_DIVA(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int SIT,bool perfect_isolation,double batch_vacc, bool do_vaccinate,bool do_batch,bool suspend_in_break,double vacc_eligible, bool retain,double trial_years,double vacc_p);
	void ExperimentWithDIVA(int rFlag, int Herd_Size,int Seeders, double trial_days, double test_period);
	vector<int> ExperimentWithDIVASinglePhase(int rFlag, int group_id, int phase_id, int Seeders, int R, int U, int V, int W, double trial_days,double test_period,bool no_vacc);
	vector<int> ExperimentWithDIVAInversePhase(int rFlag, int group_id, int phase_id, int Seeders, int R, int U, int V, int W, double trial_days,double test_period,bool no_vacc);	
void ClearStatus(int this_herd);
	void vaccinate_off_schedule(int h);
	int CountSeeders(int this_herd);
	int InvCountSeeders(int this_herd);
	~bTBICBM();
	int reacto(double p, int n);
	void save_data();
	// Simulation status flags (model outputs)
	bool	verbose;
	double	breakstart;
	double	breakfirst;
	double	break_recurr;
	double	primary_breaklength;
	double	primary_breakstart;
	int		Reactors_at_Start;
	int		Reactors_at_VE6M;
	int		Reactors_at_VE12M;
	int		Primary_Reactors;
	int		Reactors_at_First_Test;
	int 	Total_Diva_tests;
	int		Total_Diva_negatives;
	int 	Burden;
	int 	Total_Visits;
	int 	Total_Tests;
	int 	forward_trans;
	int 	onward_trans;
	
				
	bool	break6;
	bool	break12;
	bool	break24;
	bool	confirmed;
	bool 	confirmed_ever;
	bool 	confirmed_at_start;
	bool 	slaughter_house;
	bool 	retro;
				
	bool severe;
	int severe_indicator;
				
	bool breakdown;
	// toggles at end of first breakdown
	bool breakfirstFlag;
	bool follow_up;
	bool VE6Mflag;
	bool first_test;
	bool first_break;
	bool on_first_break;
	// Endsequence ensures that we bail out as soon as recurrence has occured
	bool endsequence;
	bool slaughter_recurr;
	// Introduce first infected animal at random within parish testing interval
	double short_interval;
	test_t test_type;
	
	double * DarthMoves;
	int * DarthHerdSize;
	int * DarthBindex;
	int * DarthPTI;
	
	int DarthSelecta;
	
	int DarthHerdNo;
	
	gsl_histogram * AgeReactors;
	gsl_histogram * AgeCReactors;
	gsl_histogram * AgeSReactors;
	gsl_histogram * AgeSlaughter;
		
	int maxRangeAge;
	// Functions
	void audit_maps(const char * descriptor);
	void audit_time(const char * descriptor);
	bool IsLesioned(map<double,cow_t>::iterator iter);
	void DarthLoadHerds();
	void InitialiseAgeDistro();
	void recalculate_weights();
	void recalculate_weights1();
	void recalculate_weights(int i);
	void recalculate_rate(int i);
	void update_weight(int i, map<double,cow_t>::iterator iter);
	void do_event(int h);
	void update_next_move();	
	void calculate_totals();
	void recalculate_FOI();
	bool slaughter_house_test(cow_t cow);
	bool do_management();
	int  number_testWHT(int i);
	int  number_testRHT(int i);
	void incro_cows();
	double* Get_New_Cow(int i,bool birth);
	void update_next_special();
	void vacc_reset_schedule();
	void add_vacc(int i, double timeo, bool revacc);
	void vaccinate();
	void add_bovine(int h);
	private:
	
	
	// Random Number generators
	
	MTRand mrand;
	const gsl_rng_type * gsl_T;
	gsl_rng * gsl_r;
	
	// State Variables
    // Array of H herds with cows
	
    map<double, cow_t> *cows;
	map<double, vaccination_t> vacc_schedule;
	map<double, cow_t> *infected;
	
    // Number of exposed stages
	int     TO,TR,TV1,TV2,TOV1,TOV2,TRV;
	int     TSO,TSR,TSV1,TSV2,TSOV1,TSOV2,TSRV;
	// Rate of movement between stages (1/days) 
	double  gO,gR,gV1,gV2,gOV1,gOV2,gRV;
	double  gSO,gSR,gSV1,gSV2,gSOV1,gSOV2,gSRV;

	
	// Demography
	
	int    H;   // No. of Herds
	int    TargetHerdSize;
	demo_t Demo_model;
	move_t Movement_model;
	forcing_t Seasonal_model;
	
    double *Life_Expectation;
	double *Life_var1;
	double *Occupancy_Expectation;
	double *Occupancy_var1;
	double *Life_Births;
	int NextMove;
	double NextMove_time;
	double NextSpecial_time;
	int NextSpecial;
	bool   Special_is_move;
	// Time between vaccine catchup
	// for neonates
	double batch_time;
	
	
	// Empirical Herd Model
	
	int * Darth_cows;
	double * Darth_Age_Distro;
	int * DarthHerdOffsets;
	int * DarthHerdSamples;
	
	int DarthAgeBins;
	int DarthRows;
	
	// Transmission
	double contactO,contactR,contactI,contactOV1,contactOV2,contactRV,contactIV;    // Transmission rates (from O,R,I,OV1,OV2,RV,IV)
	double contactSO,contactSR,contactSI,contactSOV1,contactSOV2,contactSRV,contactSIV;    // Transmission rates (from O,R,I,OV1,OV2,RV,IV)
	
	//seasonal forcing
	double alpha;

	double WAIFW[5*5];

	double q; // Herd Size Scaling Factor
	double xinf; // External infection pressure (herd-size dependent)
	double Vacc_Eff;
	
	int *Stot,*Otot,*Rtot,*Itot,*V1tot,*V2tot,*OV1tot,*OV2tot,*RVtot,*IVtot;
	int *SStot,*SOtot,*SRtot,*SItot,*SV1tot,*SV2tot,*SOV1tot,*SOV2tot,*SRVtot,*SIVtot;
	
	double t;  // Time elapsed in units
	
	double dt; // time_slice for approx sim
	int bini;    // Bin index
	double x;  // Random number
	int event;    // Next event
	
	// Cumulative risk of onward transmission: Infected Cow Days
	int prob_on_trans;
	
	double told; // Time before last step 
	
	
	bool full_save; //Save ALL data?
	bool debug_save; //Save extra debug data
	
	record_t recordingstatus;
	
	// Directories from which to read parameters
	// and save data
	
	string   paramdirectory;
	string   datadirectory;

	// Variable to store amount of time in different simulation methods
	double simtime;
	// Interval to sampling state variables at
	double binterval;
	
	double * Herdrates;
	double rate;
	double * FOI;
	
	// Files for data output
	ofstream *Agefile; // Age Vectors
	ofstream *Statusfile; // Epi Status
	ofstream *Lifefile; // Life Expectation Vectors
	ofstream *disclose; // Results of Test Protocol
	ofstream *Reactorfile; // Age of Reactors
	ofstream *Confirmedfile; // Age of Confirmed Reactors
	ofstream *Slaughterfile; // Age of Slaughterhouse Reactors
	ofstream *SIndividualfile; // All end-point animals

	
	// Sampling proportions for RHT
	
	double *RHT;
	int no_of_RHTs;
	
	// Empirical Distributions of Time Between Tests

	int no_of_SIT;
	int no_of_VE6M;
	int no_of_VE12M;
	int no_of_Turns;
	int no_of_PTI1;
	int no_of_PTI2;
	int no_of_PTI4;
	
	int* SIT;
	int* VE6M;
	int* VE12M;
	int* PTI1;
	int* PTI2;
	int* PTI4;
	
	// Temp pointer to track current distribution for time between routine tests
	
	int* Current_PTI;
	int Current_PTI_Num;
	
	// Buffer last testtype for Disclosure
	
	test_t current_testtype;
	
	// Test Sensitivity and Specificity
	
	double sensitivity[5];
	double specificity[5];
	
	double  pinf_move;
	double  p_shadow;
	double  p_confirm[3];
	double  delta_age;
	
	double target_vaccination_p;
	
	// Slaughter-house detection
	double slaughter_sensitivity;
	
	
	// Attempt to negate Vaccinates with DIVA
	bool DIVAnegate;
	
	// Watermark for class
	double Cuniq_id;
	
};

void bTBICBM::vaccinate_off_schedule(int h=0)
{
int num_vaccinated = 0;
int num_revaccinated = 0;

// Vaccination of whole herd
// Code to vaccinate

map<double,cow_t>::iterator iter; 
for( iter = cows[h].begin(); iter != cows[h].end(); iter++ ) 
{
if(!(iter->second).Control)
{
switch((iter->second).Epi_status)
  	{
  	case S:
  		
  		(iter->second).Epi_status = V1;
  		Stot[h]--;
  		V1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case O:
  		(iter->second).Epi_status = OV1;
  		Otot[h]--;
  		OV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case R:
  		(iter->second).Epi_status = RV;
  		Rtot[h]--;
  		RVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case I:
  		(iter->second).Epi_status = IV;
  		Itot[h]--;
  		IVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case V2:
  		(iter->second).Epi_status = V1;
  		V2tot[h]--;
  		V1tot[h]++;
  		num_revaccinated++;
  	break;
  	case SS:
  		(iter->second).Epi_status = SV1;
  		SStot[h]--;
  		SV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	  	case SO:
  		(iter->second).Epi_status = SOV1;
  		SOtot[h]--;
  		SOV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SR:
  		(iter->second).Epi_status = SRV;
  		SRtot[h]--;
  		SRVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SI:
  		(iter->second).Epi_status = SIV;
  		SItot[h]--;
  		SIVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SV2:
  		(iter->second).Epi_status = SV1;
  		SStot[h]--;
  		SV1tot[h]++;
  		num_revaccinated++;
  	break;
  	}

}	
}
	
	

	
}


void bTBICBM::vaccinate()
{

int h = (vacc_schedule.begin())->second.herd;
int num_vaccinated = 0;
int num_revaccinated = 0;

// Vaccination of whole herd
if((vacc_schedule.begin())->second.revaccinate)
{
vacc_schedule.erase(vacc_schedule.begin());
//cout << t << ' ' << "Vaccinate whole herd" << endl;

// Code to vaccinate

map<double,cow_t>::iterator iter; 
for( iter = cows[h].begin(); iter != cows[h].end(); iter++ ) 
{
if(!(iter->second).Control)
{
switch((iter->second).Epi_status)
  	{
  	case S:
  		(iter->second).Epi_status = V1;
  		(iter->second).Vaccinated_time = t;
  		Stot[h]--;
  		V1tot[h]++;
  		num_vaccinated++;
  	break;
  	case V2:
  		(iter->second).Epi_status = V1;
  		(iter->second).Vaccinated_time = t;
  		V2tot[h]--;
  		V1tot[h]++;
  		num_revaccinated++;
  	break;
  	case SS:
  		(iter->second).Epi_status = SV1;
  		(iter->second).Vaccinated_time = t;
  		SStot[h]--;
  		SV1tot[h]++;
  		num_vaccinated++;
  	break;
  	case SV2:
  		(iter->second).Epi_status = SV1;
  		(iter->second).Vaccinated_time = t;
  		SStot[h]--;
  		SV1tot[h]++;
  		num_revaccinated++;
  	break;
  	}

}

}
	
// Reschedule for same herd, add a little noise to avoid
// collisions on list of events
add_vacc(h,t+365.0+ gsl_ran_flat(gsl_r,0,0.0001),true);
}
else // Batch catchup neonates DO NOT REVACCINATE
{
//cout << t << ' ' << "Batch Catchup" << endl;

vacc_schedule.erase(vacc_schedule.begin());
// Code to vaccinate

map<double,cow_t>::iterator iter; 
for( iter = cows[h].begin(); iter != cows[h].end(); iter++ ) 
{

if(!(iter->second).Control && (iter->second).Vaccinated_time < 0.0)
{
switch((iter->second).Epi_status)
  	{
  	case S:
  		(iter->second).Epi_status = V1;
  		Stot[h]--;
  		V1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case O:
  		(iter->second).Epi_status = OV1;
  		Otot[h]--;
  		OV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case R:
  		(iter->second).Epi_status = RV;
  		Rtot[h]--;
  		RVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case I:
  		(iter->second).Epi_status = IV;
  		Itot[h]--;
  		IVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SS:
  		(iter->second).Epi_status = SV1;
  		SStot[h]--;
  		SV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	  	case SO:
  		(iter->second).Epi_status = SOV1;
  		SOtot[h]--;
  		SOV1tot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SR:
  		(iter->second).Epi_status = SRV;
  		SRtot[h]--;
  		SRVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  	case SI:
  		(iter->second).Epi_status = SIV;
  		SItot[h]--;
  		SIVtot[h]++;
  		num_vaccinated++;
  		(iter->second).Vaccinated_time = t;
  	break;
  		
  	}

}




}	

// Reschedule for same herd 
add_vacc(0,t + batch_time + gsl_ran_flat(gsl_r,0,0.0001),false);
}


update_next_special();

//cout << t << " Vaccinated:  " << num_vaccinated
//		  << " Revaccinated: " << num_revaccinated << endl;
		  
}

void bTBICBM::set_homogenous_mixing()
{
for(int i=0;i<(5*5);i++)
{
WAIFW[i] = 1.0;
}
}

void bTBICBM::set_susceptible_risk(double beta1,double beta2,double beta3,double beta4,double beta5)
{

for(int i=0; i < 5; i++)
{

WAIFW[index2(0,i,5)] = beta1;
WAIFW[index2(1,i,5)] = beta2;
WAIFW[index2(2,i,5)] = beta3;
WAIFW[index2(3,i,5)] = beta4;
WAIFW[index2(4,i,5)] = beta5;

}

/*
for(int i=0; i < 5;i++)
{
for(int j=0; j < 5;j++)
{
cout << WAIFW[index2(i,j,5)] << ' ';
}
cout << endl;

}
*/
}

void bTBICBM::set_assortative(double beta1, double beta2, double beta3, double beta4, double beta5,double beta6)
{

for(int i=0;i<(5*5);i++)
{
WAIFW[i] = beta6;
}

WAIFW[index2(0,0,5)] = beta1;
WAIFW[index2(1,1,5)] = beta2;
WAIFW[index2(2,2,5)] = beta3;
WAIFW[index2(3,4,5)] = beta4;
WAIFW[index2(4,4,5)] = beta5;


}

void bTBICBM::vacc_reset_schedule()
{


if(vacc_schedule.size()>0)
{vacc_schedule.erase(vacc_schedule.begin(),vacc_schedule.end());}
	
update_next_special();

}

void bTBICBM::add_vacc(int h, double timeo, bool revacc)
{
if(h > H) {cout << "The look of love.\n"; exit(1);}

vaccination_t pot;

pot.time = timeo;
pot.herd = h;
pot.revaccinate = revacc;

vacc_schedule.insert(std::pair<double, vaccination_t>(timeo,pot));
//cout << t << " Perform vaccination " << revacc << " at " << timeo << " in " << h <<endl; 
update_next_special();

}


void bTBICBM::set_confirmation(double *set_conf)
{

p_confirm[0] = set_conf[0];
p_confirm[1] = set_conf[1];
p_confirm[2] = set_conf[2];

//cout << "Confirmation: " << p_confirm[0] << ' ' << p_confirm[1] << ' ' << p_confirm[2] << endl;

}

void bTBICBM::set_test_characteristics(double *standard,double *severe,double *vacc1, double *vacc2, double *DIVAt,double slaughter,bool retro_slaughter=false,bool doDIVA=false)
{

sensitivity[SICCT] = standard[0];
specificity[SICCT] = standard[1];

sensitivity[SICCT_S] = severe[0];
specificity[SICCT_S] = severe[1];

sensitivity[SICCT_V1] = vacc1[0];
specificity[SICCT_V1] = vacc1[1];

sensitivity[SICCT_V2] = vacc2[0];
specificity[SICCT_V2] = vacc2[1];

sensitivity[DIVA] = DIVAt[0];
specificity[DIVA] = DIVAt[1];

slaughter_sensitivity = slaughter;
retro = retro_slaughter;
DIVAnegate = doDIVA;
}

void bTBICBM::InitialiseAgeDistro()
{
 const int bins = 21;
 //const double rangeAgeDistro[9] = {0,500,1000,1500,2000,2500,3000,3500,8000};
 const double rangeAgeDistro[22] = {0,200,400,600,800,1000,1200,1400,1600,1800,2000,2200,2400,2600,2800,3000,3200,3400,3600,3800,4000,8000};
 
 
 //const double rangeAgeDistro[9] = {0,60,120,180,240,480,720,1200,8200};
										

if(AgeReactors != NULL)
{
gsl_histogram_reset(AgeReactors);
}
else
{
AgeReactors = gsl_histogram_alloc(bins);
gsl_histogram_set_ranges(AgeReactors,rangeAgeDistro,bins+1);
}

if(AgeCReactors != NULL)
{
gsl_histogram_reset(AgeCReactors);
}
else
{
AgeCReactors = gsl_histogram_alloc(bins);
gsl_histogram_set_ranges(AgeCReactors,rangeAgeDistro,bins+1);
}

if(AgeSReactors != NULL)
{
gsl_histogram_reset(AgeSReactors);
}
else
{
AgeSReactors = gsl_histogram_alloc(bins);
gsl_histogram_set_ranges(AgeSReactors,rangeAgeDistro,bins+1);
}

if(AgeSlaughter != NULL)
{
gsl_histogram_reset(AgeSlaughter);
}
else
{
AgeSlaughter = gsl_histogram_alloc(bins);
gsl_histogram_set_ranges(AgeSlaughter,rangeAgeDistro,bins+1);
}


}

void bTBICBM::DarthLoadHerds()
{
/*
// Load Age Distribution
ifstream inputfile("/home/andrew/include/DARTH/DarthAgeDistro.csv");
if(!inputfile){cout << "Terminal velocity" << endl;exit(0);}
*/

// Load dimensions
ifstream dimensionfile((paramdirectory + "DarthDimensions.csv").c_str());
if(!dimensionfile){cout << "Fantastic Four" << endl;exit(1);}
dimensionfile >> DarthHerdNo >> DarthRows;
dimensionfile.close();

//DarthHerdNo = 6601;
//DarthRows = 1149750;

/*
DarthAgeBins = 55;

Darth_Age_Distro = (double *) calloc(sizeof(double),DarthHerdNo*DarthAgeBins);

for(int h=0; h < DarthHerdNo; h++)
{
	for(int i=0; i < DarthAgeBins; i ++)
	{
	inputfile >> Darth_Age_Distro[index2(h,i,DarthAgeBins)];
	//cout << h << ' ' << i << ' ' << Darth_Age_Distro[index2(h,i,DarthAgeBins)] << endl;
	
	}
}

inputfile.close();

if(verbose){cout << "End Age Distributions" << endl;}
*/

//Load cow table

ifstream inputfile((paramdirectory + "DarthCowTable.csv").c_str());
if(!inputfile){cout << "Rise of the Triad" << endl;exit(0);}
DarthRows = 1149750;
DarthSelecta = 0;

Darth_cows = (int *) calloc(sizeof(int),DarthRows*4);
DarthHerdOffsets = (int *) calloc(sizeof(int),DarthHerdNo);
DarthHerdSamples = (int *) calloc(sizeof(int),DarthHerdNo);

int index,CPH,age_on,move_off,death,break_age;

int r = 0;

inputfile >> index >> move_off >> age_on >> death >> break_age;

int curr_ind = index;

Darth_cows[index2(r,0,4)] = age_on;
Darth_cows[index2(r,1,4)] = move_off;
Darth_cows[index2(r,2,4)] = death;
Darth_cows[index2(r,3,4)] = break_age;

r++;

DarthHerdOffsets[0] = 0;

for(int h=0; h < DarthHerdNo; h++)
{

	while(index == curr_ind && inputfile && r < DarthRows)
	{
		inputfile >> index >> move_off >> age_on >> death >> break_age;

		Darth_cows[index2(r,0,4)] = age_on;
		Darth_cows[index2(r,1,4)] = move_off;
		Darth_cows[index2(r,2,4)] = death;
		Darth_cows[index2(r,3,4)] = break_age;
		//cout << index << ' ' << CPH << ' ' << age_on << ' ' << move_off << ' ' << death << endl;
		r++;
	} 
	
	curr_ind = index;
	

	if(h > 0)
	{
	if(h!=(DarthHerdNo-1)) 
	{DarthHerdOffsets[h+1] = r-1;
		DarthHerdSamples[h] = DarthHerdOffsets[h+1] - DarthHerdOffsets[h];
	}
	else
	{
		DarthHerdSamples[h] = (r-1) - DarthHerdOffsets[h];
	}
	
	}
	else
	{
	DarthHerdOffsets[h+1] = r-1;
	DarthHerdSamples[h] = r-1;
	}
	
}	

inputfile.close();
/*
cout << "//COW_TABLE//" << endl;
for(int r=0;r < DarthRows;r++)
{

cout << r << ' ' << Darth_cows[index2(r,0,4)] << ' ' 
	 << Darth_cows[index2(r,1,4)] << ' ' 
	 << Darth_cows[index2(r,2,4)] << ' ' 
	 << Darth_cows[index2(r,3,4)] << endl; 

}
cout << "//COW_TABLE//" << endl;
*/

/*
cout << "// DARTH OFFSETS //" << endl;
for(int h=0;h < DarthHerdNo; h++)
{
cout << h << ' ' << DarthHerdOffsets[h] << ' ' << DarthHerdSamples[h] << endl;
}
cout << "// DARTH OFFSETS //" << endl;
*/

int booty = 0;
DarthMoves = loadEmpiricalD(booty,(paramdirectory + "DarthMoves.csv").c_str());
if(booty != DarthHerdNo)
{
cout << "Fill the moat with water " << booty << endl;exit(1);
}

booty = 0;

DarthHerdSize = loadEmpirical(booty,(paramdirectory + "DarthHerds.csv").c_str());

if(booty != DarthHerdNo)
{
cout << "A body in Sri Lanka " << booty << endl;exit(1);
}

booty = 0;

DarthBindex = loadEmpirical(booty,(paramdirectory + "DarthIndex.csv").c_str());

if(booty != DarthHerdNo)
{
cout << "Fire in Sawston " << booty << endl;exit(1);
}

booty = 0;

DarthPTI = loadEmpirical(booty,(paramdirectory + "DarthPTI.csv").c_str());

if(booty != DarthHerdNo)
{
cout << "Round the Horne " << booty << endl;exit(1);
}

if(verbose){cout << "End Import" << endl;}
}


void bTBICBM::add_bovine(int h)
{

	bool update_all = true;
	double* demo_dat;
	double key_index;

	do{

	cow_t new_cow;
	
	demo_dat = Get_New_Cow(NextMove,true);
	
	new_cow.Death_time = demo_dat[Death_demo];
	new_cow.Off_time = demo_dat[Off_demo];
	new_cow.Birth_time = demo_dat[Birth_demo];
	new_cow.Infection_time = 0;
    new_cow.Vaccinated_time = -1;
    new_cow.Control = gsl_ran_binomial (gsl_r,  1.0-target_vaccination_p, 1);
	new_cow.Seeder = false;
	
	delete [] demo_dat;

	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);

	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)

	if(gsl_ran_binomial (gsl_r, p_shadow, 1))
	{
	
	if((t-new_cow.Birth_time) > 0.0)
		{
			new_cow.Epi_status = (epi_status_t) (10+gsl_ran_binomial (gsl_r, pinf_move, 1));
			if(new_cow.Epi_status == SO)
			{
				SOtot[h]++;new_cow.Infection_time = t; // Set infection time to time onto herd
			}
			else
			{SStot[h]++;}
		}
		else
		{
			new_cow.Epi_status = SS;
			SStot[h]++;
		}
	
	
	}
	else
	{
		if((t-new_cow.Birth_time) > 0.0)
		{
			new_cow.Epi_status = (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
			if(new_cow.Epi_status == O)
			{
				Otot[h]++;new_cow.Infection_time = t; // Set infection time to time onto herd
			}
			else
			{Stot[h]++;}
		}
		else
		{
			new_cow.Epi_status = S;
			Stot[h]++;
		}
	}
	
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	
	new_cow.Confirmation_status = false;
		
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;

    // check that key does not exist in herd, if it does then resample
	if(cows[h].find(key_index)==cows[h].end())
	{
	cows[h].insert(std::pair<double, cow_t>(key_index,new_cow));
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	infected[h].insert(std::pair<double, cow_t>(new_cow.uniq_id,new_cow));
	}
	break;
	}
	else{continue;}
	}while(true);
	
	recalculate_weights(h);
	
	update_next_move();
	update_next_special();

}

bool bTBICBM::do_management()
{
// Remove top animal

	//cout << "Old Herd: " << cows[NextMove].size() << endl;

	bool update_all = true;
	bool slaughter_flag = slaughter_house_test((cows[NextMove].begin())->second);
	double* demo_dat;
	double key_index;


	int recorditas = 0;
	
	// Record age of animals going to slaughter
	
	if((t - (cows[NextMove].begin()->second).Birth_time) > maxRangeAge) 
    {
    	recorditas = maxRangeAge-1;
    }
    else
    {
        recorditas = (t - (cows[NextMove].begin()->second).Birth_time);
    }
	
	gsl_histogram_increment(AgeSlaughter,recorditas);
	
	
	switch(Seasonal_model)
	{
	case constant:
	
	if(debug_save)
	{
	SIndividualfile[NextMove] 
	<< DarthSelecta << ','
	<< 1-(cows[NextMove].begin()->second).Control << ','
	<< (cows[NextMove].begin()->second).Birth_time << ','
	<< (cows[NextMove].begin()->second).Death_time << ',' 
	<< (cows[NextMove].begin()->second).Vaccinated_time << ',' 
	<< 0 << ','
	<< slaughter_flag << endl;
	}	
	//cout << "French Cricket" << endl;
	if(slaughter_flag)
	{
	if(full_save)
	{
	Slaughterfile[NextMove] << DarthSelecta << ' ' << (t - (cows[NextMove].begin()->second).Birth_time) 
	<< ' ' << (cows[NextMove].begin()->second).Epi_status << endl;
	}

	if((t - (cows[NextMove].begin()->second).Birth_time) > maxRangeAge) 
    {
    	recorditas = maxRangeAge-1;
    }
    else
    {
        recorditas = (t - (cows[NextMove].begin()->second).Birth_time);
    }
	
	gsl_histogram_increment(AgeSReactors,recorditas);
	
	}
	
	
	switch((cows[NextMove].begin())->second.Epi_status)
	{
	case S:
	Stot[NextMove]--;
	update_all = false;
	break;
	case O:
	
	Otot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case R:

	Rtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case I:
	
	Itot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; 		exit(1);}


	break;
	case V1:
	V1tot[NextMove]--;
	update_all = false;
	break;
	case V2:
	V2tot[NextMove]--;
	update_all = false;
	break;
	case OV1:
	
	OV1tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case OV2:
	
	OV2tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case RV:
	RVtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}
	break;
	case IV:
	IVtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}
	break;
	case SS:
	SStot[NextMove]--;
	update_all = false;
	break;
	case SO:

	SOtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SR:

	SRtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SI:
	
	SItot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SV1:
	SV1tot[NextMove]--;
	update_all = false;
	break;
	case SV2:
	SV2tot[NextMove]--;
	update_all = false;
	break;
	case SOV1:
	
	SOV1tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SOV2:

	SOV2tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SRV:

	SRVtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SIV:

	SIVtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	}
	
	cows[NextMove].erase(cows[NextMove].begin());
	
	
	do{	
	demo_dat = Get_New_Cow(NextMove,true);
		
	cow_t new_cow;
	
	new_cow.Death_time = demo_dat[Death_demo];
	new_cow.Off_time = demo_dat[Off_demo];
	new_cow.Birth_time = demo_dat[Birth_demo];
	new_cow.Infection_time = 0;

	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);

	delete [] demo_dat;

	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)

	if(gsl_ran_binomial (gsl_r, p_shadow, 1))
	{
	
	if((t-new_cow.Birth_time) > 0.0)
		{
			new_cow.Epi_status = (epi_status_t) (10+gsl_ran_binomial (gsl_r, pinf_move, 1));
			if(new_cow.Epi_status == SO)
			{
				SOtot[NextMove]++;new_cow.Infection_time = t; // Set infection time to time onto herd
				
			}
			else
			{SStot[NextMove]++;}
		}
		else
		{
			new_cow.Epi_status = SS;
			SStot[NextMove]++;
		}
	
	
	}
	else
	{
		if((t-new_cow.Birth_time) > 0.0)
		{
			new_cow.Epi_status = (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
			if(new_cow.Epi_status == O)
			{
				Otot[NextMove]++;new_cow.Infection_time = t; // Set infection time to time onto herd
				
			}
			else
			{Stot[NextMove]++;}
		}
		else
		{
			new_cow.Epi_status = S;
			Stot[NextMove]++;
		}
	}
	
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	new_cow.Confirmation_status = false;
	new_cow.Vaccinated_time = -1;
    new_cow.Control = gsl_ran_binomial (gsl_r,  1.0-target_vaccination_p, 1);
   	new_cow.Seeder = false;
   	
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
	    // check that key does not exist in herd, if it does then resample
	if(cows[NextMove].find(key_index)==cows[NextMove].end())
	{
	cows[NextMove].insert(std::pair<double, cow_t>(key_index,new_cow));
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	infected[NextMove].insert(std::pair<double, cow_t>(new_cow.uniq_id,new_cow));
	}
	break;
	}
	else{continue;}
	}while(true);
	recalculate_weights(NextMove);
	
	 
	break;
	case seasonal:
		if(debug_save)
	{
	SIndividualfile[NextMove] 
	<< DarthSelecta << ','
	<< 1-(cows[NextMove].begin()->second).Control << ','
	<< (cows[NextMove].begin()->second).Birth_time << ','
	<< (cows[NextMove].begin()->second).Death_time << ',' 
	<< (cows[NextMove].begin()->second).Vaccinated_time << ',' 
	<< 0 << ','
	<< slaughter_flag << endl;
	}
	
	if(slaughter_flag){
	if(full_save)
	{
	Slaughterfile[NextMove] << DarthSelecta << ' ' << (t - (cows[NextMove].begin()->second).Birth_time) << endl;
	}
	if((t - (cows[NextMove].begin()->second).Birth_time) > maxRangeAge) 
    {
    	recorditas = maxRangeAge-1;
    }
    else
    {
        recorditas = (t - (cows[NextMove].begin()->second).Birth_time);
    }
	
	gsl_histogram_increment(AgeSReactors,recorditas);
	}

	switch((cows[NextMove].begin())->second.Epi_status)
	{
	case S:
	Stot[NextMove]--;
	update_all = false;
	break;
	case O:
	Otot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case R:
	Rtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case I:
	Itot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case V1:
	V1tot[NextMove]--;
	update_all = false;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case V2:
	V2tot[NextMove]--;
	update_all = false;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case OV1:
	OV1tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case OV2:
	OV2tot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case RV:
	RVtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case IV:
	IVtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SS:
	SStot[NextMove]--;
	update_all = false;
	break;
	case SO:
	SOtot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SR:
	SRtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SI:
	SItot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SV1:
	SV1tot[NextMove]--;
	update_all = false;
	break;
	case SV2:
	SV2tot[NextMove]--;
	update_all = false;
	break;
	case SOV1:
	SOV1tot[NextMove]--;
	if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SOV2:
	SOV2tot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SRV:
	SRVtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	case SIV:
	SIVtot[NextMove]--;
		if(infected[NextMove].erase(cows[NextMove].begin()->first)!=1){cout<<"Wendy"<<endl; exit(1);}

	break;
	}
	cows[NextMove].erase(cows[NextMove].begin());
	//cout << "Kill cow " << cows[NextMove].size() << endl;
	//if(update_all)
	{
		recalculate_weights(NextMove);
	}
	//else
	//{
	//	recalculate_rate(NextMove);
	//}
	
	break;
	}

//cout << "New Herd: " << cows[NextMove].size() << endl;

update_next_move();
update_next_special();

return(slaughter_flag);
}

void bTBICBM::update_next_special()
{

if(vacc_schedule.size() > 0 && NextMove_time > (vacc_schedule.begin())->second.time)
 {
	NextSpecial_time = (vacc_schedule.begin())->second.time;
	NextSpecial=(vacc_schedule.begin())->second.herd;
	Special_is_move = false;
 }
 else
 {
    NextSpecial_time = NextMove_time;
	NextSpecial=NextMove;
	Special_is_move = true; 
 }

}


void bTBICBM::update_next_move()
{

for(int i=0; i < H; i++)
{

 if(cows[i].size()>0)
 {
 
 if((cows[i].begin())->first <= (cows[NextMove].begin())->first)
 {
 NextMove = i;NextMove_time = (*cows[i].begin()).first;
 //if(NextMove_time < 0){ cout << "Headless Doll" << endl;}
 }
 
 }

}
//cout << t << " Nextmove: " << NextMove << " at " << NextMove_time << endl;

update_next_special();

}

void bTBICBM::print_demo()
{


map<double,cow_t>::iterator iter;
for( iter = cows[0].begin(); iter != cows[0].end(); iter++ ) 
  {
cout << DarthSelecta << ',' << (t-iter->second.Birth_time) << endl;
 }
cout << endl;
}



void bTBICBM::recalculate_FOI()
{
// Annual cohorts, everything greater than 4 years dumped into final class
 //cout << "This many infected: " << infected[h].size() << endl;

for(int h=0; h < H; h++)
{

double multipass = (1/(pow((double) cows[h].size() / 165.0,q)));

for(int i=0; i < 5; i++)
{

double lambda = 0;
map<double,cow_t>::iterator iter;
map<double,cow_t>::iterator seek;
for( iter = infected[h].begin(); iter != infected[h].end(); iter++ ) 
  {
 	// Select Age Class
   int j = floor((t-iter->second.Birth_time)/364.0);if(j > 4){j=4;}
   
   if(j < 0){cout << "BIG WARNING: " << i << ' ' << j << endl;
   			 cout << "BIG WARNING: " << t << ' ' << iter->second.Birth_time << endl;
   
   cout << "Herd: " << DarthSelecta << endl;
   
   exit(1);
   
   }
   
   // Infected map status is not updated
   // so use value from corresponding animal in cows map
   
   //cout << iter->first << endl;
   
   seek = cows[h].find(iter->first);
   if(seek == cows[h].end()){cout << "Missing infected!" << endl; exit(1);}
   
   map<double,cow_t>::iterator iter2;
   switch(seek->second.Epi_status)
   {
    case S:
    	cout << "THIS SHOULD NEVER HAPPEN!" << endl;
    	cout << "Time: " << t << endl;
    	cout.precision(15);
    	
    	cout << "Naughty member of infected: " << iter->second.uniq_id << endl;
    	cout << (iter->second).Epi_status << endl;
  		cout << (iter->second).Birth_time << endl;
    	cout << (iter->second).Death_time << endl;
    	cout << (iter->second).Off_time << endl;
   		cout << (iter->second).Infection_time << endl;
    	cout << (iter->second).Vaccinated_time << endl;
		cout << (iter->second).Epi_stage << endl;
		cout << (iter->second).Std_status << endl;
		cout << (iter->second).Svr_status << endl;
		cout << (iter->second).Diva_status << endl;
		cout << (iter->second).Confirmation_status << endl;
		cout << (iter->second).rate << endl;
    	
    	cout << "Matches: " << seek->second.uniq_id << endl;
    	cout << (seek->second).Epi_status << endl;
  		cout << (seek->second).Birth_time << endl;
    	cout << (seek->second).Death_time << endl;
    	cout << (seek->second).Off_time << endl;
   		cout << (seek->second).Infection_time << endl;
    	cout << (seek->second).Vaccinated_time << endl;
		cout << (seek->second).Epi_stage << endl;
		cout << (seek->second).Std_status << endl;
		cout << (seek->second).Svr_status << endl;
		cout << (seek->second).Diva_status << endl;
		cout << (seek->second).Confirmation_status << endl;
		cout << (seek->second).rate << endl;
    	
    	cout << "Herd Follows" << endl;
    	for( iter2 = cows[h].begin(); iter2 != cows[h].end(); iter2++ ) 
  		{
			cout << "Bovine: " << iter2->second.uniq_id << endl;
			cout << (iter2->second).Epi_status << endl;
  			cout << (iter2->second).Birth_time << endl;
    		cout << (iter2->second).Death_time << endl;
    		cout << (iter2->second).Off_time << endl;
   			cout << (iter2->second).Infection_time << endl;
    		cout << (iter2->second).Vaccinated_time << endl;
    		cout << (iter2->second).Epi_stage << endl;
			cout << (iter2->second).Std_status << endl;
			cout << (iter2->second).Svr_status << endl;
			cout << (iter2->second).Diva_status << endl;
			cout << (iter2->second).Confirmation_status << endl;
			cout << (iter2->second).rate << endl;
			
		}
    	exit(1);
  	case O:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactO);
  	break;
  	case R:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactR);
  	break;
  	case I:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactI);
  	break;
  	
  	case OV1:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactOV1);
  	break;
  	case OV2:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactOV2);
  	break;
  	case RV:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactRV);
  	break;
  	case IV:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactIV);
  	break;
  	
  	case SO:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSO);
  	break;
  	case SR:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSR);
  	break;
  	case SI:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSI);
  	break;

  	case SOV1:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSOV1);
  	break;
  	case SOV2:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSOV2);
  	break;
  	case SRV:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSRV);
  	break;
  	case SIV:
  		lambda += WAIFW[index2(i,j,5)]*(multipass*contactSIV);
  	break;
   default:
  						cout << "For once in my life: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
   }

  }

FOI[index2(h,i,5)] = lambda + xinf*multipass*WAIFW[index2(i,i,5)];
//FOI[index2(h,i,5)] = lambda + xinf;
//cout << "This FOI: " << i << ' ' << h << ' ' << lambda << ' ' <<  xinf*WAIFW[index2(i,i,5)]/(pow((double) cows[h].size(),q)) << ' ' << FOI[index2(h,i,5)] << endl;

}
}


}

void bTBICBM::recalculate_weights()
{
 double ermine=0.0;
 rate = 0.0;
 recalculate_FOI();

 for(int i=0; i < H; i++)
 {
 
 int N = cows[i].size();
 

 //FOI[i] = (1.0/(pow(N/165.0,q)))*(contactO*Otot[i] + contactR*Rtot[i] + contactI*Itot[i] + contactOV1*OV1tot[i] + contactOV2*OV2tot[i] + contactRV*RVtot[i] + contactIV*IVtot[i] 
 //							+ contactSO*SOtot[i] + contactSR*SRtot[i] + contactSI*SItot[i] + contactSOV1*SOV1tot[i] + contactSOV2*SOV2tot[i] + contactSRV*SRVtot[i] + contactSIV*SIVtot[i] ) + xinf;
 Herdrates[i] = Life_Births[i]*(1+alpha*cos(2*M_PI*t/364.0));
 //cout << "Force of Infection: " << FOI[i] << endl;
 //cout << "LifeBirths: " << Life_Births[i] << endl;

  map<double,cow_t>::iterator iter;   
  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  {
   int a_c = floor((t-iter->second.Birth_time)/364.0);
   if(a_c > 4){a_c=4;}
    
    //cout << a_c << ' ' << FOI[index2(i,a_c,5)] << endl;
    
   switch((iter->second).Epi_status)
  {
  case S:
  		
  		Herdrates[i] += ((iter->second).rate=FOI[index2(i,a_c,5)]);
  		//cout << "S : " << (iter->second).rate << endl;
  break;
  case O:
  		
  		Herdrates[i] += ((iter->second).rate=gO);
  		//cout << "O : " << (iter->second).rate << endl;
  break;
  case R:
  		
    	Herdrates[i] += ((iter->second).rate=gR);
    	//cout << "R : " << (iter->second).rate << endl;
  break;
  case I:
  		
  		(iter->second).rate = 0.0;
  		//cout << "I : " << (iter->second).rate << endl;
  break;
  case V1:
  		(iter->second).rate = gV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V1 : " <<  ' ' << Vacc_Eff*FOI[index2(i,a_c,5)] << ' ' << (iter->second).rate << endl;
  break;
  case V2:
  		(iter->second).rate=FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V2 : " << (iter->second).rate << endl;
  break;
  case OV1:
  		
  		Herdrates[i] += ((iter->second).rate=gOV1);
  		//cout << "OV1 : " << (iter->second).rate << endl;
  break;
  case OV2:
  		
  		Herdrates[i] += ((iter->second).rate=gOV2);
  		//cout << "OV2 : " << (iter->second).rate << endl;
  break;
  case RV:
  		
  		Herdrates[i] += ((iter->second).rate=gRV);
  	//cout << "RV : " << (iter->second).rate << endl;
  break;
  case IV:
  	
  		(iter->second).rate = 0.0;
  		//cout << "IV : " << (iter->second).rate << endl;
	  break;
    case SS:
  		//cout << "SS Rate" << endl;
  		Herdrates[i] += ((iter->second).rate=FOI[index2(i,a_c,5)]);
  		//cout << "SS : " << (iter->second).rate << endl;
  break;
  case SO:
  		
  		Herdrates[i] += ((iter->second).rate=gSO);
  		//cout << "SO : " << (iter->second).rate << endl;
  break;
  case SR:
  		
    	Herdrates[i] += ((iter->second).rate=gSR);
    	//cout << "SR : " << (iter->second).rate << endl;
  break;
  case SI:
  		
  		(iter->second).rate = 0.0;
  		//cout << "SI : " << (iter->second).rate << endl;
  break;
  case SV1:
  		
  		(iter->second).rate = gSV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "SV1 : " << (iter->second).rate << endl;
  break;
  case SV2:
  		(iter->second).rate=Vacc_Eff*FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "SV2 : " << (iter->second).rate << endl;
  break;
  case SOV1:
  		
  		Herdrates[i] += ((iter->second).rate=gSOV1);
  		//cout << "SOV1 : " << (iter->second).rate << endl;
  break;
  case SOV2:
  		Herdrates[i] += ((iter->second).rate=gSOV2);
  		//cout << "SOV2 : " << (iter->second).rate << endl;
  break;
  case SRV:
  		//cout << "SRV : " << (iter->second).rate << endl;
  		Herdrates[i] += ((iter->second).rate=gSRV);
  break;
  case SIV:
  		(iter->second).rate = 0.0;
//cout << "SIV : " << (iter->second).rate << endl;
  break;
  default:
  						cout << "Days and nights flying by: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  }
 
   
 }
 rate += Herdrates[i];
 //cout << i << " Herdrates: " << Herdrates[i] << " New rate: " << rate << endl;

 }
 
}


void bTBICBM::recalculate_weights(int i)
{

 int N = cows[i].size();
 recalculate_FOI();
 Herdrates[i] = Life_Births[i]*(1+alpha*cos(2*M_PI*t/364.0));
// cout << "Force of Infection: " << FOI[i] << endl;
 //cout << "LifeBirths: " << Life_Births[i] << endl;
  
  map<double,cow_t>::iterator iter;   
  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  {
    int a_c = floor((t-iter->second.Birth_time)/364.0);
   if(a_c > 4){a_c=4;}
   switch((iter->second).Epi_status)
  {
  case S:
  		Herdrates[i] += ((iter->second).rate=FOI[index2(i,a_c,5)]);
  		//cout << "S: " << Herdrates[i] << endl;
  		//cout << "S: " << FOI[index2(i,a_c,5)] << endl;
  break;
  case O:
  		Herdrates[i] += ((iter->second).rate=gO);
  		//cout << "O: " << Herdrates[i] << endl;
  break;
  case R:
    	Herdrates[i] += ((iter->second).rate=gR);
    	//cout << "R: " << Herdrates[i] << endl;
  break;
  case I:
  		(iter->second).rate = 0.0;
  		//cout << "I: " << Herdrates[i] << endl;
  break;
  case V1:
  		(iter->second).rate = gV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V1: " << Herdrates[i] << endl;
  		//cout << "V1: " << Vacc_Eff*FOI[index2(i,a_c,5)] << endl;
  break;
  case V2:
  		(iter->second).rate = FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V2: " << Herdrates[i] << endl;
  break;
  case OV1:
  		Herdrates[i] += ((iter->second).rate=gOV1);
  		//cout << "OV1: " << Herdrates[i] << endl;
  break;
  case OV2:
  		Herdrates[i] += ((iter->second).rate=gOV2);
  		//cout << "OV2: " << Herdrates[i] << endl;
  break;
  case RV:
  		Herdrates[i] += ((iter->second).rate=gRV);
  		//cout << "RV: " << Herdrates[i] << endl;
  break;
  case IV:
  		(iter->second).rate = 0.0;
  		//cout << "IV: " << Herdrates[i] << endl;
  break;
    case SS:
  		Herdrates[i] += ((iter->second).rate=FOI[index2(i,a_c,5)]);
  		//cout << "SS: " << Herdrates[i] << endl;
  break;
  case SO:
  		Herdrates[i] += ((iter->second).rate=gSO);
  		//cout << "O: " << Herdrates[i] << endl;
  break;
  case SR:
    	Herdrates[i] += ((iter->second).rate=gSR);
    	//cout << "R: " << Herdrates[i] << endl;
  break;
  case SI:
  		(iter->second).rate = 0.0;
  		//cout << "I: " << Herdrates[i] << endl;
  break;
  case SV1:
  		(iter->second).rate = gSV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V1: " << Herdrates[i] << endl;
  break;
  case SV2:
  		(iter->second).rate+=FOI[index2(i,a_c,5)];
  		Herdrates[i] += (iter->second).rate;
  		//cout << "V2: " << Herdrates[i] << endl;
  break;
  case SOV1:
  		Herdrates[i] += ((iter->second).rate=gSOV1);
  		//cout << "OV1: " << Herdrates[i] << endl;
  break;
  case SOV2:
  		Herdrates[i] += ((iter->second).rate=gSOV2);
  		//cout << "OV2: " << Herdrates[i] << endl;
  break;
  case SRV:
  		Herdrates[i] += ((iter->second).rate=gSRV);
  		//cout << "RV: " << Herdrates[i] << endl;
  break;
  case SIV:
  		(iter->second).rate = 0.0;
  		//cout << "IV: " << Herdrates[i] << endl;
  break;
  default:
  						cout << "I'm not alone anymore: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  }
   
 }
 
 rate = 0.0;
 for(int j=0; j < H; j++)
 {
 rate += Herdrates[j];
 }
 
}

void bTBICBM::update_weight(int i, map<double,cow_t>::iterator iter)
{
 int N = cows[i].size();
 int a_c = floor((t-iter->second.Birth_time)/364.0);
 if(a_c > 4){a_c=4;}
 recalculate_FOI();
 

 switch(iter->second.Epi_status)
 {
  case S:
		((iter->second).rate=FOI[index2(i,a_c,5)]);
		recalculate_rate(i);
		//recalculate_weights();
		// << "S: " << FOI[i] << endl;
  break;
  case O:
		recalculate_weights(i);
  break;
  case R:
		recalculate_weights(i);
  break;
  case I:
		recalculate_weights(i);
  break;
  case V1:
  		
		(iter->second).rate = gV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		recalculate_rate(i);
  		//recalculate_weights();
  		//cout << "V1: " << Vacc_Eff*FOI[i] << endl;
  		//cout << "V1: " << Herdrates[i] << endl;
  break;
  case V2:
  		
  		(iter->second).rate=FOI[index2(i,a_c,5)];
  		recalculate_rate(i);
  		//recalculate_weights();
  		//cout << "V2: " << FOI[i] << endl;
  		//cout << "V2: " << Herdrates[i] << endl;
  break;
  case OV1:
  		recalculate_weights(i);
  break;
  case OV2:
  		recalculate_weights(i);
  break;
  case RV:
  		recalculate_weights(i);
  break;
  case IV:
		recalculate_weights(i);
  break;
  case SS:
		((iter->second).rate=FOI[index2(i,a_c,5)]);
		recalculate_rate(i);
		//cout << "SS: " << FOI[i] << endl;
  break;
  case SO:
		recalculate_weights(i);
  break;
  case SR:
		recalculate_weights(i);
  break;
  case SI:
		recalculate_weights(i);
  break;
  case SV1:
  		
	 	(iter->second).rate = gSV1;
  		(iter->second).rate += Vacc_Eff*FOI[index2(i,a_c,5)];
  		recalculate_rate(i);
  		//cout << "SV1: " << FOI[i] << endl;
  		//cout << "V1: " << Herdrates[i] << endl;
  break;
  case SV2:
  		
  		(iter->second).rate=FOI[index2(i,a_c,5)];
  		recalculate_rate(i);
  		//cout << "SV2: " << FOI[i] << endl;
  		//cout << "V2: " << Herdrates[i] << endl;
  break;
  case SOV1:
  		recalculate_weights(i);
  break;
  case SOV2:
  		recalculate_weights(i);
  break;
  case SRV:
  		recalculate_weights(i);
  break;
  case SIV:
		recalculate_weights(i);
  break;
  default:
  						cout << "Who needs me." << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  }
}
 
void bTBICBM::recalculate_rate(int i)
{
  
  Herdrates[i] = Life_Births[i]*(1+alpha*cos(2*M_PI*t/364.0));

  map<double,cow_t>::iterator iter;   
  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  {
  		
  		Herdrates[i] += (iter->second).rate;
  }
 
  rate = 0.0;
  for(int j=0; j < H; j++)
  {
  	rate += Herdrates[j];
  }  
 
}


bTBICBM::bTBICBM(int Herds, int OrderO, int OrderR, int OrderV1, int OrderV2, int OrderOV1,int OrderOV2, int OrderRV, int OrderSO, int OrderSR, int OrderSV1, int OrderSV2, int OrderSOV1,int OrderSOV2, int OrderSRV, string input, string output,bool setfull_save,double latent_O, double latent_R, double latent_V1,double latent_V2, double latent_OV1, double latent_OV2, double latent_RV,double latent_SO, double latent_SR, double latent_SV1,double latent_SV2, double latent_SOV1, double latent_SOV2, double latent_SRV, bool varbose,double Setp_inf=0.0,double Setdelta_age=0.0, bool setdebug_save=false)
{
	verbose = varbose;
	
	recordingstatus = overwrite;
	
	full_save = setfull_save;
	debug_save = setdebug_save;
   	
	// Dimensions
	
	H = Herds;
	
	paramdirectory = input;
	datadirectory = output;
	
	// Decide type of model
	
	// EXP model
	if (OrderO == 0)
	{
		TO      = 0;
		gO      = 1.0 / latent_O;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TO      = OrderO;
		gO      = (TO / latent_O);       
	}
	
	
	// EXP model
	
	// negative latent Period corresponds to absorbing state
	if(latent_R < 0.0)
	{
	 gR = 0.0;
	 TR = 1;
	}
	else
	{
	if (OrderR == 0)
	{
		TR      = 0;
		gR      = 1.0 / latent_R;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TR      = OrderR;
		gR      = (TR / latent_R);       
	}
	}
	
	// EXP model
	if (OrderV1 == 0)
	{
		TV1      = 1;
		gV1      = 1.0 / latent_V1;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TV1      = OrderV1;
		gV1      = (TV1 / latent_V1);       
	}
	
		// EXP model
	if (OrderV2 == 0)
	{
		TV2     = 0;
		gV2      = 1.0 / latent_V2;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TV2      = OrderV2;
		gV2      = (TV2 / latent_V2);       
	}
	
		// EXP model
	if (OrderOV1 == 0)
	{
		TOV1      = 0;
		gOV1      = 1.0 / latent_OV1;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TOV1      = OrderV1;
		gOV1      = (TOV1 / latent_OV1);       
	}
	
		// EXP model
	if (OrderV2 == 0)
	{
		TOV2     = 0;
		gOV2      = 1.0 / latent_OV2;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TOV2      = OrderOV2;
		gOV2      = (TOV2 / latent_OV2);       
	}
	
	// EXP model
	
	// negative latent Period corresponds to absorbing state
	if(latent_RV < 0.0)
	{
	 gRV = 0.0;
	 TRV = 1;
	}
	else
	{
	if (OrderRV == 0)
	{
		TRV      = 0;
		gRV      = 1.0 / latent_RV;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TRV      = OrderRV;
		gRV      = (TRV / latent_RV);       
	}
	}

	// EXP model
	if (OrderO == 0)
	{
		TSO      = 0;
		gSO      = 1.0 / latent_SO;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSO      = OrderSO;
		gSO      = (TSO / latent_SO);       
	}
	
	
	// EXP model
	
	// negative latent Period corresponds to absorbing state
	if(latent_SR < 0.0)
	{
	 gSR = 0.0;
	 TSR = 1;
	}
	else
	{
	if (OrderSR == 0)
	{
		TSR      = 0;
		gSR      = 1.0 / latent_SR;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSR      = OrderSR;
		gSR      = (TSR / latent_SR);       
	}
	}
	
	// EXP model
	if (OrderSV1 == 0)
	{
		TSV1      = 1;
		gSV1      = 1.0 / latent_SV1;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSV1      = OrderSV1;
		gSV1      = (TSV1 / latent_SV1);       
	}
	
		// EXP model
	if (OrderSV2 == 0)
	{
		TSV2     = 0;
		gSV2      = 1.0 / latent_SV2;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSV2      = OrderSV2;
		gSV2      = (TSV2 / latent_SV2);       
	}
	
		// EXP model
	if (OrderSOV1 == 0)
	{
		TSOV1      = 0;
		gSOV1      = 1.0 / latent_SOV1;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSOV1      = OrderSV1;
		gSOV1      = (TSOV1 / latent_SOV1);       
	}
	
		// EXP model
	if (OrderSV2 == 0)
	{
		TSOV2     = 0;
		gSOV2      = 1.0 / latent_SOV2;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSOV2      = OrderSOV2;
		gSOV2      = (TSOV2 / latent_SOV2);       
	}
	
	// EXP model
	
	// negative latent Period corresponds to absorbing state
	if(latent_SRV < 0.0)
	{
	 gSRV = 0.0;
	 TSRV = 1;
	}
	else
	{
	if (OrderSRV == 0)
	{
		TSRV      = 0;
		gSRV      = 1.0 / latent_SRV;
	}
	// Gamma-Distributed model with OrderE stages
	else
	{
		TSRV      = OrderSRV;
		gSRV      = (TSRV / latent_SRV);       
	}
	}


	// Dimensions now immutable

    cows = new map<double, cow_t>[H];
    infected = new map<double,cow_t>[H];
	
	Demo_model = LifeExp;
	Movement_model = uniform;
	Seasonal_model = constant;
			
	Life_Expectation = new double[H];
	Life_var1 = new double[H];
	Occupancy_Expectation = new double[H];
	Occupancy_var1 = new double[H];
	
	Life_Births = new double[H];
	Herdrates = new double[H];
	pinf_move = Setp_inf;
	delta_age  = Setdelta_age;
	p_confirm[0] = 0.0;
	p_confirm[1] = 0.0;
	p_confirm[2] = 0.0;
	p_shadow = 0.0;
	
	Stot = new int[H];Otot = new int[H];Rtot = new int[H];Itot = new int[H];
	V1tot = new int[H];V2tot = new int[H];OV1tot = new int[H];OV2tot = new int[H];
	RVtot = new int[H];IVtot = new int[H];
	
	SStot = new int[H];SOtot = new int[H];SRtot = new int[H];SItot = new int[H];
	SV1tot = new int[H];SV2tot = new int[H];SOV1tot = new int[H];SOV2tot = new int[H];
	SRVtot = new int[H];SIVtot = new int[H];
	
	for(int i=0;i < H;i++)
	{
	Life_Expectation[i] = yeardef*3.0;
	Life_var1[i] = 0.0;
	Occupancy_Expectation[i] = yeardef*3.0;
	Occupancy_var1[i] = 0.0;
	Life_Births[i] = 0.0;
	Herdrates[i] = 0.0;
	Stot[i]=0;SStot[i]=0;
	Otot[i]=0;SOtot[i]=0;
	Rtot[i]=0;SRtot[i]=0;
	Itot[i]=0;SItot[i]=0;
	V1tot[i]=0;SV1tot[i]=0;
	V2tot[i]=0;SV2tot[i]=0;
	OV1tot[i]=0;SOV1tot[i]=0;
	OV2tot[i]=0;SOV2tot[i]=0;
	RVtot[i]=0;SRVtot[i]=0;
	IVtot[i]=0;SIVtot[i]=0;
	}
	
	NextMove=0;
	NextMove_time = 0.0;
	
	NextSpecial_time = -1.0;
	NextSpecial=0;
	Special_is_move = false;
	
	//Initialise all contact rates to zero
	
	contactO=0;contactR=0;contactI=0;contactOV1=0;contactOV2=0;contactRV=0;contactIV=0; 
	contactSO=0;contactSR=0;contactSI=0;contactSOV1=0;contactSOV2=0;contactSRV=0;contactSIV=0; 
	
	for(int i=0;i<(5*5);i++)
	{WAIFW[i] = 0.0;}
	
	q = 0.0;
	
	xinf = 0.0;
	
	Vacc_Eff = 0.0;
	
	FOI = new double[H*5];
	
	for(int i=0;i<(H*5);i++)
	{FOI[i] = 0.0;}
	
	Statusfile = new ofstream[H];
	Agefile = new ofstream[H];
	Lifefile = new ofstream[H];
	Reactorfile = new ofstream[H];
	Slaughterfile = new ofstream[H];
	Confirmedfile = new ofstream[H];
	SIndividualfile = new ofstream[H];
	
	if(debug_save)
	{
		if(verbose){cout << "Creating State Variable Time-Series" << endl;}
		for(int i=0; i < H; i++)
		{
			Statusfile[i].open((datadirectory + filename("StatusH",i)).c_str());
			Agefile[i].open((datadirectory + filename("AgeH",i)).c_str());
			Lifefile[i].open((datadirectory + filename("LifeH",i)).c_str());
		}
	}
	else{if(verbose){cout << "Creating the bare minimum files. \n";}}
	disclose = new ofstream[H];
	if(full_save)
	{
	for(int i=0; i < H; i++)
    {
		disclose[i].open((datadirectory + filename("DiscloseH",i)).c_str());
		disclose[i] << "Time " << "R " << "RSev " << "TrueR " << "Burden " << "Tested " << "Confirmed " << "StdS " << "StdO " << "StdR " << "StdI " << "SvrS " << "SvrO " << "SvrR " << "SvrI " << "StdV1 " << "StdV2 " << "StdOV1 " << "StdOV2 "  << "StdRV " << "StdIV " << "SvrV1 " << "SvrV2 " << "SvrOV1 " << "SvrOV2 "  << "SvrRV " << "SvrIV " << "TestS " << "TestO " << "TestR " << "TestI " << "TestV1 " << "TestV2 " << "TestOV1 " << "TestOV2 " << "TestIV " << "TestRV " << "stdDIVAtests " << "svrDIVAtests " << "stdDIVAnegate " << "svrDIVAnegate " << "Replicate " << "Removed " << "PTI " << "HerdSize " << "TestType " << "OnHerd ";
    	disclose[i] << "RU " << "RV " << "LU " << "LV " << "SU " << "IU " << "SV " << "IV" << endl;
    	
    	Reactorfile[i].open((datadirectory + filename("ReactorsH",i)).c_str());
    	Slaughterfile[i].open((datadirectory + filename("ReactorsSH",i)).c_str());
    	Confirmedfile[i].open((datadirectory + filename("ReactorsCH",i)).c_str());
    	
    	if(debug_save)
    	{
    	SIndividualfile[i].open((datadirectory + filename("Individual",i)).c_str());
    	
    	SIndividualfile[i] << "Herd, Vaccinated, Birth_time, Death_time, Vaccination_time, Reactor, VL" << endl;
    	}
    }
	}
	// Cumulative Risk of transmission to another farm
	
	prob_on_trans = 0;
	
	
	// Initialise Time
	
	t = 0.0;  // Time elapsed in days (current realisation)
	rate = 0.0;  // Current total rate of events
	
	
	bini = 0;    // Bin index
	x = 0.0;  // Random number
	event= 0;    // Next event
	
	told = 0.0; // Time before last step
	
	// Initialise some other stuff for sanity  
	
	binterval = 0.0;
		
	// Initialise gsl random number generation
	// Mersenne Twister
	
	#pragma omp critical(rnginitialise)
	{
	gsl_rng_default_seed = mrand.rand() * ULONG_MAX;
	//gsl_rng_default_seed = 10024;
	gsl_r = gsl_rng_alloc(gsl_rng_mt19937);
	Cuniq_id = gsl_rng_default_seed;
	//cout << "Uniq ID: " << Cuniq_id << endl;
	}
	
	slaughter_sensitivity = 0.0;	   
    
	//RHT = loadEmpiricalD(no_of_RHTs,"RHT.csv");
	
	// Empirical Distributions of Time Between Tests

	no_of_SIT = 0;
	no_of_VE6M = 0;
	no_of_VE12M = 0;
	no_of_PTI1  = 0;
	no_of_PTI2  = 0;
	no_of_PTI4  = 0;
	
	// Maximum values for approximate empirical distributions

	SIT = loadEmpirical(no_of_SIT,(paramdirectory + "SITtrunk.csv").c_str());
	VE6M = loadEmpirical(no_of_VE6M,(paramdirectory + "VE-6M.csv").c_str());
	VE12M = loadEmpirical(no_of_VE12M,(paramdirectory + "VE-12M.csv").c_str());
	PTI1 = loadEmpirical(no_of_PTI1,(paramdirectory + "PTI1_Times.csv").c_str());
	PTI2 = loadEmpirical(no_of_PTI2,(paramdirectory + "PTI2_Times.csv").c_str());
	PTI4 = loadEmpirical(no_of_PTI4,(paramdirectory + "PTI4_Times.csv").c_str());
	
	int* Current_PTI = NULL;
	int Current_PTI_Num = 0;

	// Simulation status flags (model outputs)
	
	breakstart=0.0;
	breakfirst=-1;
	break_recurr=0.0;
	primary_breaklength = 0.0;
	primary_breakstart  = 0.0;
	Reactors_at_Start = 0;
	Reactors_at_VE6M  = 0;
	Reactors_at_VE12M = 0;
	Primary_Reactors  = 0;
	Reactors_at_First_Test = 0;
	forward_trans = 0;
	break6=false;
	break12=false;
	break24=false;
	confirmed=false;
	confirmed_ever=false;
	confirmed_at_start = false;
	slaughter_house = false;
	retro = false;
				
	severe=false;
	severe_indicator = 1;
				
	breakdown = false;
	// toggles at end of first breakdown
	breakfirstFlag = false;
	follow_up = false;
	VE6Mflag  = false;
	first_test = true;
	first_break = true;
	on_first_break = true;
	// Endsequence ensures that we bail out as soon as recurrence has occured
	endsequence = false;
	slaughter_recurr = false;
	// Introduce first infected animal at random within parish testing interval
	//double introducebTB = mrand.rand()*test_period;
	//double introducebTB = 0.0;
	short_interval = 60;
	
	Burden			= 0;
	Total_Visits    = 0;
	Total_Tests		= 0;
	forward_trans 	= 0;
	onward_trans 	= 0;
	
	// Empirical Herd Model

	DarthHerdNo=0;
	DarthAgeBins=55;
	DarthRows=0;
	
	DarthLoadHerds();
	
	AgeReactors = NULL;
	AgeCReactors = NULL;
	AgeSReactors = NULL;
	AgeSlaughter = NULL;
	
	InitialiseAgeDistro();
	
	maxRangeAge = 7920;
	
	set_homogenous_mixing();
	
	sensitivity[0]=0.0;sensitivity[1]=0.0;sensitivity[2]=0.0;sensitivity[3]=0.0;sensitivity[4]=0.0;
	specificity[0]=0.0;specificity[1]=0.0;specificity[2]=0.0;specificity[3]=0.0;specificity[4]=0.0;
	slaughter_sensitivity = 0.0;
	
	TargetHerdSize = 0;
	
	alpha=0.0;
	//alpha=1.0;
	
	DIVAnegate=false;
	
	// Target proportion of vaccinates and controls
	target_vaccination_p = 0.0;
	
}

// Copy constructor for cloning herd

bTBICBM::bTBICBM(const bTBICBM& other)
{
	verbose = other.verbose;
	
	recordingstatus = other.recordingstatus;
	
	full_save = other.full_save;
	debug_save = other.debug_save;
   	
	// Dimensions
	
	H = other.H;
	
	paramdirectory = other.paramdirectory;
	datadirectory = other.datadirectory;
	
	TO      = other.TO;
	gO      = other.gO;
		
	gR = other.gR;
	TR = other.TR;
	
	TV1      = other.TV1;
	gV1      = other.gV1;
	
	TV2     = other.TV2;
	gV2      = other.gV2;

	TOV1      = other.TOV1;
	gOV1      = other.gOV1;
	
	TOV2     = other.TOV2;
	gOV2      = other.gOV2;

	gRV = other.gRV;
	TRV = other.TRV;

	TSO      = other.TSO;
	gSO      = other.gSO;
	
    gSR = other.gSR;
	TSR = other.TSR;
	
	TSV1      = other.TSV1;
	gSV1      = other.gSV1;
	
	TSV2     = other.TSV2;
	gSV2      = other.gSV2;
	
	TSOV1      = other.TSOV1;
	gSOV1      = other.gSOV1;
	
	TSOV2     = other.TSOV2;
	gSOV2      = other.gSOV2;
	
	gSRV = other.gSRV;
	TSRV = other.TSRV;
	
	// Dimensions now immutable

    cows = new map<double, cow_t>[H];
    infected = new map<double,cow_t>[H];
	
	// Copy cows map
	
	for(int z=0; z < H; z++)
	{
	map<double,cow_t>::iterator iter;   
  	for( iter = other.cows[z].begin(); iter != other.cows[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	cows[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	
	}
	}
	
	// Copy infected map
	
	for(int z=0; z < H; z++)
	{
	map<double,cow_t>::iterator iter;
    map<double,cow_t>::iterator seek;   
  	for( iter = other.infected[z].begin(); iter != other.infected[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	
	infected[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	
	// Check key is in cow map
	seek = cows[z].find(iter->first);
   if(seek == cows[z].end()){cout << "Dratted infected!" << endl; exit(1);}

	}
	}
	
	
	Demo_model = other.Demo_model;
	Movement_model = other.Movement_model;
	Seasonal_model = other.Seasonal_model;
			
	Life_Expectation = new double[H];
	Life_var1 = new double[H];
	Occupancy_Expectation = new double[H];
	Occupancy_var1 = new double[H];
	
	Life_Births = new double[H];
	Herdrates = new double[H];
	pinf_move = other.pinf_move;
	delta_age  = other.delta_age;
	p_confirm[0] = other.p_confirm[0];
	p_confirm[1] = other.p_confirm[1];
	p_confirm[2] = other.p_confirm[2];
	p_shadow = other.p_shadow;
	
	Stot = new int[H];Otot = new int[H];Rtot = new int[H];Itot = new int[H];
	V1tot = new int[H];V2tot = new int[H];OV1tot = new int[H];OV2tot = new int[H];
	RVtot = new int[H];IVtot = new int[H];
	
	SStot = new int[H];SOtot = new int[H];SRtot = new int[H];SItot = new int[H];
	SV1tot = new int[H];SV2tot = new int[H];SOV1tot = new int[H];SOV2tot = new int[H];
	SRVtot = new int[H];SIVtot = new int[H];
	
	for(int i=0;i < H;i++)
	{
	Life_Expectation[i] = other.Life_Expectation[i];
	Life_var1[i] = other.Life_var1[i];
	Occupancy_Expectation[i] = other.Occupancy_Expectation[i];
	Occupancy_var1[i] = other.Occupancy_var1[i];
	Life_Births[i] = other.Life_Births[i];
	Herdrates[i] = other.Herdrates[i];
	Stot[i]=other.Stot[i];SStot[i]=other.SStot[i];
	Otot[i]=other.Otot[i];SOtot[i]=other.SOtot[i];
	Rtot[i]=other.Rtot[i];SRtot[i]=other.SRtot[i];
	Itot[i]=other.Itot[i];SItot[i]=other.SItot[i];
	V1tot[i]=other.V1tot[i];SV1tot[i]=other.SV1tot[i];
	V2tot[i]=other.V2tot[i];SV2tot[i]=other.SV2tot[i];
	OV1tot[i]=other.OV1tot[i];SOV1tot[i]=other.SOV1tot[i];
	OV2tot[i]=other.OV2tot[i];SOV2tot[i]=other.SOV2tot[i];
	RVtot[i]=other.RVtot[i];SRVtot[i]=other.SRVtot[i];
	IVtot[i]=other.IVtot[i];SIVtot[i]=other.SIVtot[i];
	}
	
	NextMove=other.NextMove;
	NextMove_time = other.NextMove_time;
	
	NextSpecial_time = other.NextSpecial_time;
	NextSpecial=other.NextSpecial;
	Special_is_move = other.Special_is_move;
	
	//Initialise all contact rates
	
	contactO=other.contactO;
	contactR=other.contactR;
	contactI=other.contactI;
	contactOV1=other.contactOV1;
	contactOV2=other.contactOV2;
	contactRV=other.contactRV;
	contactIV=other.contactIV; 
	contactSO=other.contactSO;
	contactSR=other.contactSR;
	contactSI=other.contactSI;
	contactSOV1=other.contactSOV1;
	contactSOV2=other.contactSOV2;
	contactSRV=other.contactSRV;
	contactSIV=other.contactSIV; 
	
	for(int i=0;i<(5*5);i++)
	{WAIFW[i] = other.WAIFW[i];}
	
	q = other.q;
	
	xinf = other.xinf;
	
	Vacc_Eff = other.Vacc_Eff;
	
	FOI = new double[H*5];
	
	for(int i=0;i<(H*5);i++)
	{FOI[i] = other.FOI[i];}
	
	Statusfile = new ofstream[H];
	Agefile = new ofstream[H];
	Lifefile = new ofstream[H];
	Reactorfile = new ofstream[H];
	Slaughterfile = new ofstream[H];
	Confirmedfile = new ofstream[H];
	SIndividualfile = new ofstream[H];
	
	// Open no output files for copied objects
		
	disclose = new ofstream[H];
	
	// Cumulative Risk of transmission to another farm
	
	prob_on_trans = other.prob_on_trans;
	
	
	// Initialise Time
	
	t = other.t;  // Time elapsed in days (current realisation)
	rate = other.rate;  // Current total rate of events
	
	
	bini = other.bini;    // Bin index
	x = other.x;  // Random number
	event= other.event;    // Next event
	
	told = other.told; // Time before last step
	
	// Initialise some other stuff for sanity  
	
	binterval = other.binterval;
		
	// Initialise gsl random number generation
	// Mersenne Twister
	
	#pragma omp critical(rnginitialise)
	{
	gsl_rng_default_seed = mrand.rand() * ULONG_MAX;
	//gsl_rng_default_seed = 10024;
	gsl_r = gsl_rng_alloc(gsl_rng_mt19937);
	Cuniq_id = gsl_rng_default_seed;
	//cout << "Uniq ID: " << Cuniq_id << endl;
	}
	
	slaughter_sensitivity = other.slaughter_sensitivity;	   
    
	//RHT = loadEmpiricalD(no_of_RHTs,"RHT.csv");
	
	// Empirical Distributions of Time Between Tests

	no_of_SIT = 0;
	no_of_VE6M = 0;
	no_of_VE12M = 0;
	no_of_PTI1  = 0;
	no_of_PTI2  = 0;
	no_of_PTI4  = 0;
	
	// Maximum values for approximate empirical distributions

	SIT = loadEmpirical(no_of_SIT,(paramdirectory + "SITtrunk.csv").c_str());
	VE6M = loadEmpirical(no_of_VE6M,(paramdirectory + "VE-6M.csv").c_str());
	VE12M = loadEmpirical(no_of_VE12M,(paramdirectory + "VE-12M.csv").c_str());
	PTI1 = loadEmpirical(no_of_PTI1,(paramdirectory + "PTI1_Times.csv").c_str());
	PTI2 = loadEmpirical(no_of_PTI2,(paramdirectory + "PTI2_Times.csv").c_str());
	PTI4 = loadEmpirical(no_of_PTI4,(paramdirectory + "PTI4_Times.csv").c_str());
	
	int* Current_PTI = PTI1;
	int Current_PTI_Num = 0;

	// Simulation status flags (model outputs)
	
	breakstart=other.breakstart;
	breakfirst=other.breakfirst;
	break_recurr=other.break_recurr;
	primary_breaklength = other.primary_breaklength;
	primary_breakstart  = other.primary_breakstart;
	Reactors_at_Start = other.Reactors_at_Start;
	Reactors_at_VE6M  = other.Reactors_at_VE6M;
	Reactors_at_VE12M = other.Reactors_at_VE12M;
	Primary_Reactors  = other.Primary_Reactors;
	Reactors_at_First_Test = other.Reactors_at_First_Test;
	forward_trans = other.forward_trans;
	break6=other.break6;
	break12=other.break12;
	break24=other.break24;
	confirmed=other.confirmed;
	confirmed_ever=other.confirmed_ever;
	confirmed_at_start = other.confirmed_at_start;
	slaughter_house = other.slaughter_house;
	retro = other.retro;
				
	severe=other.severe;
	severe_indicator = other.severe_indicator;
				
	breakdown = other.breakdown;
	// toggles at end of first breakdown
	breakfirstFlag = other.breakfirstFlag;
	follow_up = other.follow_up;
	VE6Mflag  = other.VE6Mflag;
	first_test = other.first_test;
	first_break = other.first_break;
	on_first_break = other.on_first_break;
	// Endsequence ensures that we bail out as soon as recurrence has occured
	endsequence = other.endsequence;
	slaughter_recurr = other.slaughter_recurr;
	// Introduce first infected animal at random within parish testing interval
	//double introducebTB = mrand.rand()*test_period;
	//double introducebTB = 0.0;
	short_interval = other.short_interval;
	
	Burden			= other.Burden;
	Total_Visits    = other.Total_Visits;
	Total_Tests		= other.Total_Tests;
	forward_trans 	= other.forward_trans;
	onward_trans 	= other.onward_trans;
	
	// Empirical Herd Model

	DarthHerdNo=other.DarthHerdNo;
	DarthAgeBins=other.DarthAgeBins;
	DarthRows=other.DarthRows;
	
	DarthLoadHerds();
	
	AgeReactors = NULL;
	AgeCReactors = NULL;
	AgeSReactors = NULL;
	AgeSlaughter = NULL;
	
	InitialiseAgeDistro();
	
	maxRangeAge = 7920;
	
	sensitivity[0]=other.sensitivity[0];
	sensitivity[1]=other.sensitivity[1];
	sensitivity[2]=other.sensitivity[2];
	sensitivity[3]=other.sensitivity[3];
	sensitivity[4]=other.sensitivity[4];
	specificity[0]=other.specificity[0];
	specificity[1]=other.specificity[1];
	specificity[2]=other.specificity[2];
	specificity[3]=other.specificity[3];
	specificity[4]=other.specificity[4];
	slaughter_sensitivity = other.slaughter_sensitivity;
	
	TargetHerdSize = other.TargetHerdSize;
	
	alpha=other.alpha;
	
	DIVAnegate=other.DIVAnegate;
	
	// Target proportion of vaccinates and controls
	target_vaccination_p = other.target_vaccination_p;
	
	audit_maps("Copy Constructor");
	
	
}

// Merge constructor for transmission experiments herd

bTBICBM::bTBICBM(const bTBICBM& other, const bTBICBM& one, bool pick_controls)
{
	//Take all basic parameters from other
	//Pick relevant animals (flag matches with pick_controls) from one.
	
	//cout << "Cows Other: " << other.cows[0].size() << ' ' << "Cows One: " << one.cows[0].size() << endl;
	
	verbose = other.verbose;
	
	recordingstatus = other.recordingstatus;
	
	full_save = other.full_save;
	debug_save = other.debug_save;
   	
	// Dimensions
	
	H = other.H;
	
	paramdirectory = other.paramdirectory;
	datadirectory = other.datadirectory;
	
	TO      = other.TO;
	gO      = other.gO;
		
	gR = other.gR;
	TR = other.TR;
	
	TV1      = other.TV1;
	gV1      = other.gV1;
	
	TV2     = other.TV2;
	gV2      = other.gV2;

	TOV1      = other.TOV1;
	gOV1      = other.gOV1;
	
	TOV2     = other.TOV2;
	gOV2      = other.gOV2;

	gRV = other.gRV;
	TRV = other.TRV;

	TSO      = other.TSO;
	gSO      = other.gSO;
	
    gSR = other.gSR;
	TSR = other.TSR;
	
	TSV1      = other.TSV1;
	gSV1      = other.gSV1;
	
	TSV2     = other.TSV2;
	gSV2      = other.gSV2;
	
	TSOV1      = other.TSOV1;
	gSOV1      = other.gSOV1;
	
	TSOV2     = other.TSOV2;
	gSOV2      = other.gSOV2;
	
	gSRV = other.gSRV;
	TSRV = other.TSRV;
	
	// Dimensions now immutable

    cows = new map<double, cow_t>[H];
    infected = new map<double,cow_t>[H];
	
	// Copy cows map
	
	for(int z=0; z < H; z++)
	{
	
	map<double,cow_t>::iterator iter;   
  	for( iter = other.cows[z].begin(); iter != other.cows[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	
	// Only keep relevant vaccinated (controls=false) or unvaccinated animals 
	// (controls=true)
	if(new_cow.Control==pick_controls)
	{
	cows[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	}
	
	}
	
  
  	for( iter = one.cows[z].begin(); iter != one.cows[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	
	// Only keep relevant vaccinated (controls=false) or unvaccinated animals 
	// (controls=true)
	if(new_cow.Control==pick_controls)
	{
	cows[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	}
	
	}
	
	
	}
	
	// Copy infected map
	
	for(int z=0; z < H; z++)
	{
	map<double,cow_t>::iterator iter;
    map<double,cow_t>::iterator seek;   
  	for( iter = other.infected[z].begin(); iter != other.infected[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	// Only retain relevant animals flagged by Control
	if(new_cow.Control==pick_controls)
	{
	infected[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	// Check key is in cow map
	seek = cows[z].find(iter->first);
   if(seek == cows[z].end()){cout << "Dratted infected!" << endl; exit(1);}

	}
	
	}
	
	for( iter = one.infected[z].begin(); iter != one.infected[z].end(); iter++ ) 
	{
	
	cow_t new_cow;

	new_cow.Birth_time = (iter->second).Birth_time;
    new_cow.Death_time = (iter->second).Death_time;
    new_cow.Off_time = (iter->second).Off_time;
    new_cow.Infection_time = (iter->second).Infection_time;
    new_cow.Vaccinated_time = (iter->second).Vaccinated_time;
	new_cow.Epi_status = (iter->second).Epi_status;
	new_cow.Epi_stage = (iter->second).Epi_stage;
	new_cow.Std_status = (iter->second).Std_status;
	new_cow.Svr_status = (iter->second).Svr_status;
	new_cow.Diva_status = (iter->second).Diva_status;
	new_cow.Diva_ever = (iter->second).Diva_ever;	
	new_cow.Confirmation_status = (iter->second).Confirmation_status;
	new_cow.Control = (iter->second).Control;
	new_cow.Seeder = (iter->second).Seeder;
	new_cow.rate = (iter->second).rate;
	new_cow.uniq_id = (iter->second).uniq_id;
	
	// Only retain relevant animals flagged by Control
	if(new_cow.Control==pick_controls)
	{
	infected[z].insert(std::pair<double, cow_t>(iter->first,new_cow));
	// Check key is in cow map
	seek = cows[z].find(iter->first);
   if(seek == cows[z].end()){cout << "Dratted infected!" << endl; exit(1);}

	
	}
	
	}
	
	
	
	}
	
	
	Demo_model = other.Demo_model;
	Movement_model = other.Movement_model;
	Seasonal_model = other.Seasonal_model;
			
	Life_Expectation = new double[H];
	Life_var1 = new double[H];
	Occupancy_Expectation = new double[H];
	Occupancy_var1 = new double[H];
	
	Life_Births = new double[H];
	Herdrates = new double[H];
	pinf_move = other.pinf_move;
	delta_age  = other.delta_age;
	p_confirm[0] = other.p_confirm[0];
	p_confirm[1] = other.p_confirm[1];
	p_confirm[2] = other.p_confirm[2];
	p_shadow = other.p_shadow;
	
	Stot = new int[H];Otot = new int[H];Rtot = new int[H];Itot = new int[H];
	V1tot = new int[H];V2tot = new int[H];OV1tot = new int[H];OV2tot = new int[H];
	RVtot = new int[H];IVtot = new int[H];
	
	SStot = new int[H];SOtot = new int[H];SRtot = new int[H];SItot = new int[H];
	SV1tot = new int[H];SV2tot = new int[H];SOV1tot = new int[H];SOV2tot = new int[H];
	SRVtot = new int[H];SIVtot = new int[H];
	
	for(int i=0;i < H;i++)
	{
	Life_Expectation[i] = other.Life_Expectation[i];
	Life_var1[i] = other.Life_var1[i];
	Occupancy_Expectation[i] = other.Occupancy_Expectation[i];
	Occupancy_var1[i] = other.Occupancy_var1[i];
	Life_Births[i] = other.Life_Births[i];
	Herdrates[i] = other.Herdrates[i];
	Stot[i]=other.Stot[i];SStot[i]=other.SStot[i];
	Otot[i]=other.Otot[i];SOtot[i]=other.SOtot[i];
	Rtot[i]=other.Rtot[i];SRtot[i]=other.SRtot[i];
	Itot[i]=other.Itot[i];SItot[i]=other.SItot[i];
	V1tot[i]=other.V1tot[i];SV1tot[i]=other.SV1tot[i];
	V2tot[i]=other.V2tot[i];SV2tot[i]=other.SV2tot[i];
	OV1tot[i]=other.OV1tot[i];SOV1tot[i]=other.SOV1tot[i];
	OV2tot[i]=other.OV2tot[i];SOV2tot[i]=other.SOV2tot[i];
	RVtot[i]=other.RVtot[i];SRVtot[i]=other.SRVtot[i];
	IVtot[i]=other.IVtot[i];SIVtot[i]=other.SIVtot[i];
	}
	
	NextMove=other.NextMove;
	NextMove_time = other.NextMove_time;
	
	NextSpecial_time = other.NextSpecial_time;
	NextSpecial=other.NextSpecial;
	Special_is_move = other.Special_is_move;
	
	//Initialise all contact rates
	
	contactO=other.contactO;
	contactR=other.contactR;
	contactI=other.contactI;
	contactOV1=other.contactOV1;
	contactOV2=other.contactOV2;
	contactRV=other.contactRV;
	contactIV=other.contactIV; 
	contactSO=other.contactSO;
	contactSR=other.contactSR;
	contactSI=other.contactSI;
	contactSOV1=other.contactSOV1;
	contactSOV2=other.contactSOV2;
	contactSRV=other.contactSRV;
	contactSIV=other.contactSIV; 
	
	for(int i=0;i<(5*5);i++)
	{WAIFW[i] = other.WAIFW[i];}
	
	q = other.q;
	
	xinf = other.xinf;
	
	Vacc_Eff = other.Vacc_Eff;
	
	FOI = new double[H*5];
	
	for(int i=0;i<(H*5);i++)
	{FOI[i] = other.FOI[i];}
	
	Statusfile = new ofstream[H];
	Agefile = new ofstream[H];
	Lifefile = new ofstream[H];
	Reactorfile = new ofstream[H];
	Slaughterfile = new ofstream[H];
	Confirmedfile = new ofstream[H];
	SIndividualfile = new ofstream[H];
	
	// Open no output files for copied objects
		
	disclose = new ofstream[H];
	
	// Cumulative Risk of transmission to another farm
	
	prob_on_trans = other.prob_on_trans;
	
	
	// Initialise Time
	
	t = other.t;  // Time elapsed in days (current realisation)
	rate = other.rate;  // Current total rate of events
	
	
	bini = other.bini;    // Bin index
	x = other.x;  // Random number
	event= other.event;    // Next event
	
	told = other.told; // Time before last step
	
	// Initialise some other stuff for sanity  
	
	binterval = other.binterval;
		
	// Initialise gsl random number generation
	// Mersenne Twister
	
	#pragma omp critical(rnginitialise)
	{
	gsl_rng_default_seed = mrand.rand() * ULONG_MAX;
	//gsl_rng_default_seed = 10024;
	gsl_r = gsl_rng_alloc(gsl_rng_mt19937);
	Cuniq_id = gsl_rng_default_seed;
	//cout << "Uniq ID: " << Cuniq_id << endl;
	}
	
	slaughter_sensitivity = other.slaughter_sensitivity;	   
    
	//RHT = loadEmpiricalD(no_of_RHTs,"RHT.csv");
	
	// Empirical Distributions of Time Between Tests

	no_of_SIT = 0;
	no_of_VE6M = 0;
	no_of_VE12M = 0;
	no_of_PTI1  = 0;
	no_of_PTI2  = 0;
	no_of_PTI4  = 0;
	
	// Maximum values for approximate empirical distributions

	SIT = loadEmpirical(no_of_SIT,(paramdirectory + "SITtrunk.csv").c_str());
	VE6M = loadEmpirical(no_of_VE6M,(paramdirectory + "VE-6M.csv").c_str());
	VE12M = loadEmpirical(no_of_VE12M,(paramdirectory + "VE-12M.csv").c_str());
	PTI1 = loadEmpirical(no_of_PTI1,(paramdirectory + "PTI1_Times.csv").c_str());
	PTI2 = loadEmpirical(no_of_PTI2,(paramdirectory + "PTI2_Times.csv").c_str());
	PTI4 = loadEmpirical(no_of_PTI4,(paramdirectory + "PTI4_Times.csv").c_str());
	
	int* Current_PTI = PTI1;
	int Current_PTI_Num = 0;

	// Simulation status flags (model outputs)
	
	breakstart=other.breakstart;
	breakfirst=other.breakfirst;
	break_recurr=other.break_recurr;
	primary_breaklength = other.primary_breaklength;
	primary_breakstart  = other.primary_breakstart;
	Reactors_at_Start = other.Reactors_at_Start;
	Reactors_at_VE6M  = other.Reactors_at_VE6M;
	Reactors_at_VE12M = other.Reactors_at_VE12M;
	Primary_Reactors  = other.Primary_Reactors;
	Reactors_at_First_Test = other.Reactors_at_First_Test;
	forward_trans = other.forward_trans;
	break6=other.break6;
	break12=other.break12;
	break24=other.break24;
	confirmed=other.confirmed;
	confirmed_ever=other.confirmed_ever;
	confirmed_at_start = other.confirmed_at_start;
	slaughter_house = other.slaughter_house;
	retro = other.retro;
				
	severe=other.severe;
	severe_indicator = other.severe_indicator;
				
	breakdown = other.breakdown;
	// toggles at end of first breakdown
	breakfirstFlag = other.breakfirstFlag;
	follow_up = other.follow_up;
	VE6Mflag  = other.VE6Mflag;
	first_test = other.first_test;
	first_break = other.first_break;
	on_first_break = other.on_first_break;
	// Endsequence ensures that we bail out as soon as recurrence has occured
	endsequence = other.endsequence;
	slaughter_recurr = other.slaughter_recurr;
	// Introduce first infected animal at random within parish testing interval
	//double introducebTB = mrand.rand()*test_period;
	//double introducebTB = 0.0;
	short_interval = other.short_interval;
	
	Burden			= other.Burden;
	Total_Visits    = other.Total_Visits;
	Total_Tests		= other.Total_Tests;
	forward_trans 	= other.forward_trans;
	onward_trans 	= other.onward_trans;
	
	// Empirical Herd Model

	DarthHerdNo=other.DarthHerdNo;
	DarthAgeBins=other.DarthAgeBins;
	DarthRows=other.DarthRows;
	
	DarthLoadHerds();
	
	AgeReactors = NULL;
	AgeCReactors = NULL;
	AgeSReactors = NULL;
	AgeSlaughter = NULL;
	
	InitialiseAgeDistro();
	
	maxRangeAge = 7920;
	
	sensitivity[0]=other.sensitivity[0];
	sensitivity[1]=other.sensitivity[1];
	sensitivity[2]=other.sensitivity[2];
	sensitivity[3]=other.sensitivity[3];
	sensitivity[4]=other.sensitivity[4];
	specificity[0]=other.specificity[0];
	specificity[1]=other.specificity[1];
	specificity[2]=other.specificity[2];
	specificity[3]=other.specificity[3];
	specificity[4]=other.specificity[4];
	slaughter_sensitivity = other.slaughter_sensitivity;
	
	TargetHerdSize = other.TargetHerdSize;
	
	alpha=other.alpha;
	
	DIVAnegate=other.DIVAnegate;
	
	// Target proportion of vaccinates and controls
	target_vaccination_p = other.target_vaccination_p;
	
	// Recalculate variables dependent on map objects which have just been merged
	update_next_move();
	calculate_totals();
	
	audit_maps("Merge Constructor");
	
	
}

void bTBICBM::audit_time(const char * descriptor) 
{
 if(t < 0.0){cout << "Backwards: " << t << ' ' << descriptor;exit(1);}
}


void bTBICBM::audit_maps(const char * descriptor) 
{
int master_inf[H];

	for(int i=0; i < H; i++)
	{
	master_inf[i]=0;

  map<double,cow_t>::iterator iter;   
  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  {
	
	switch((iter->second).Epi_status)
  	{
  	case S:

  	break;
  	case O:
  		master_inf[i]++;
  	break;
  	case R:
  		master_inf[i]++;
  	break;
  	case I:
  		master_inf[i]++;
  	break;
  	case V1:
  		
  	break;
  	case V2:
  		
  	break;
  	case OV1:
  		master_inf[i]++;
  	break;
  	case OV2:
  		master_inf[i]++;
  	break;
  	case RV:
  		master_inf[i]++;
  	break;
  	case IV:
  		master_inf[i]++;
  	break;
  	  	case SS:
  		
  	break;
  	case SO:
  		master_inf[i]++;
  	break;
  	case SR:
  		master_inf[i]++;
  	break;
  	case SI:
  		master_inf[i]++;
  	break;
  	case SV1:
  		master_inf[i]++;
  	break;
  	case SV2:
  		master_inf[i]++;
  	break;
  	case SOV1:
  		master_inf[i]++;
  	break;
  	case SOV2:
  		master_inf[i]++;
  	break;
  	case SRV:
  		master_inf[i]++;
  	break;
  	case SIV:
  		master_inf[i]++;
  	break;
  	default:
  						cout << "I get a kick out of you: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  						
  	}
	
	}

	#pragma omp critical(dataout)
{
if((master_inf[i]) != infected[i].size())
{
cout << "TinTin: " << Cuniq_id << ' ' << master_inf[i] << ' ' << infected[i].size() << ' ' << cows[i].size() << ' ' << descriptor << endl;
exit(1);
}

}
}
	
}

void bTBICBM::calculate_totals()
{
//audit_maps("Biggles");
	for(int i=0; i < H; i++)
	{
	Stot[i]=0;SStot[i]=0;
	Otot[i]=0;SOtot[i]=0;
	Rtot[i]=0;SRtot[i]=0;
	Itot[i]=0;SItot[i]=0;
	V1tot[i]=0;SV1tot[i]=0;
	V2tot[i]=0;SV2tot[i]=0;
	OV1tot[i]=0;SOV1tot[i]=0;
	OV2tot[i]=0;SOV2tot[i]=0;
	RVtot[i]=0;SRVtot[i]=0;
	IVtot[i]=0;SIVtot[i]=0;

  map<double,cow_t>::iterator iter;   
  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  {
	
	switch((iter->second).Epi_status)
  	{
  	case S:
  		Stot[i]++;
  	break;
  	case O:
  		Otot[i]++;
  	break;
  	case R:
  		Rtot[i]++;
  	break;
  	case I:
  		Itot[i]++;
  	break;
  	case V1:
  		V1tot[i]++;
  	break;
  	case V2:
  		V2tot[i]++;
  	break;
  	case OV1:
  		OV1tot[i]++;
  	break;
  	case OV2:
  		OV2tot[i]++;
  	break;
  	case RV:
  		RVtot[i]++;
  	break;
  	case IV:
  		IVtot[i]++;
  	break;
  	  	case SS:
  		SStot[i]++;
  	break;
  	case SO:
  		SOtot[i]++;
  	break;
  	case SR:
  		SRtot[i]++;
  	break;
  	case SI:
  		SItot[i]++;
  	break;
  	case SV1:
  		SV1tot[i]++;
  	break;
  	case SV2:
  		SV2tot[i]++;
  	break;
  	case SOV1:
  		SOV1tot[i]++;
  	break;
  	case SOV2:
  		SOV2tot[i]++;
  	break;
  	case SRV:
  		SRVtot[i]++;
  	break;
  	case SIV:
  		SIVtot[i]++;
  	break;
  	default:
  						cout << "I get a kick out of you: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  						
  	}

	
	}


}
	//audit_maps("Benji");
/*for(int i=0; i < H; i++)
	{
	
	 cout << "States: "
	 	  << Stot[i] << ' '
		  << Otot[i] << ' '
		  << Rtot[i] << ' '
	      << Itot[i] << ' '
	      << V1tot[i] << ' '
		  << V2tot[i] << ' '
 		  << OV1tot[i] << ' '
		  << OV2tot[i] << ' '
		  << RVtot[i] << ' '
		  << IVtot[i] << ' '
		  << cows[i].size() << ' '
		  << endl;
	
	}
	

	
*/

						// #pragma omp critical(dataout)
						//{
						//cout << "Check: " << t << ' ' << DarthSelecta << ' ' << Otot[0] + Rtot[0]+ Itot[0] + OV1tot[0] + OV2tot[0] + RVtot[0] + IVtot[0] << ' ' << infected[0].size() << endl;
						//if((Otot[0] + Rtot[0]+ Itot[0] + OV1tot[0] + OV2tot[0] + RVtot[0] + IVtot[0]) != infected[0].size())
//{
//cout << "Check: " << Cuniq_id << ' ' << DarthSelecta << ' ' << Otot[0] + Rtot[0]+ Itot[0] + OV1tot[0] + OV2tot[0] + RVtot[0] + IVtot[0] << ' ' << infected[0].size() << endl;
//exit(1);
//}
//}

}

void bTBICBM::cull_FP(int Mary)
{
//audit_maps("Cull FP");

	/*
	if(S[0] >= Mary)
	{
		S[0] -= Mary;
	}
	*/
	
}


bool bTBICBM::initialise_herd(int Herd_size,int herd_index, int Marys, epi_status_t Mary_Type)
{
		//audit_maps("Borked Before");

	if(herd_index > H) {cout << "You and me, round about midnight.\n"; return false;}
	
	int these[Marys], this_index[Herd_size];
          
    for (int i = 0; i < Herd_size; i++)
    {
    	this_index[i] = (double) i;
    }
          
    gsl_ran_choose(gsl_r, these, Marys, this_index, Herd_size, sizeof (int));
     	
		
	//cout << "Marys : " << Marys << endl;
	// Empty herd if any cows are present
	if(cows[herd_index].size()>0)
	{cows[herd_index].erase(cows[herd_index].begin(),cows[herd_index].end());}
	// Empty infection list
	if(infected[herd_index].size()>0)
	{infected[herd_index].erase(infected[herd_index].begin(),infected[herd_index].end());}
	

	for(int i=0;i < Herd_size; i++)
	{
	double *noo_coo = Get_New_Cow(herd_index,false);
	double key_index = 0.0;
	
	cow_t new_cow;
	
	new_cow.Birth_time = noo_coo[Birth_demo];
	new_cow.Off_time = noo_coo[Off_demo];
	new_cow.Death_time = noo_coo[Death_demo];
	new_cow.Infection_time = 0;
	
	/*cout << (t-new_cow.Birth_time) << ' ' 
		 << (new_cow.Birth_time) << ' ' 
		 << new_cow.Off_time << ' '
		 << new_cow.Death_time << endl;
	*/
	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);
	
	delete [] noo_coo;
	
	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)
	new_cow.Epi_status = S;
	//new_cow.Epi_status = ((t-new_cow.Birth_time) > 0.0) ? S : (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
	
	//cout << "Status " << new_cow.Epi_status << endl;
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	new_cow.Confirmation_status   = false;
	 new_cow.Vaccinated_time = -1;
    new_cow.Control = gsl_ran_binomial (gsl_r,  1.0-target_vaccination_p, 1);
	new_cow.Seeder = false;

	// Is cow shadow-progressor?

	epi_status_t shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
	new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	//cout << "Shadow: " << shadow << ' ' << new_cow.Epi_status << endl;
	
	// Infect Mary's animals
	for(int k=0; k < Marys; k++)
	 { 
	 if(i==these[k]){new_cow.Epi_status = (epi_status_t) (new_cow.Epi_status + 1);
	 				
	// cout << "Infection! " << k << ' ' << new_cow.Epi_status << endl;
	 }
	 }
	 
		
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
	
	// check that key does not exist in herd, if it does then resample
	if(cows[herd_index].find(key_index)==cows[herd_index].end())
	{
	cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
	}
	else{ 
	i--;continue;}
	
	
	//std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	//ret = cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	//if (ret.second==false) {
    //std::cout << "BOVINE COLLISION" << endl;
    // Final safety check for key collision
    //exit(1);
  	
  	
  	
	
	//cows[herd_index][key_index] = new_cow;
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	ret = infected[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	if (ret.second==false) {
    std::cout << "element 'z' already existed";
  	}
  
	}
	
	}

	////audit_maps("HerdInit");

	 
	if(NextMove==herd_index)
	{
	//if(NextMove_time < 0){ cout << "Decapitated Cheerleader" << endl;}
	NextMove_time = (*cows[herd_index].begin()).first;
	}
	
	////audit_maps("NextHerd");
	update_next_move();
	////audit_maps("Update next move");
	calculate_totals();



	return true;
}


bool bTBICBM::initialise_herd_AllInAllOut(int Herd_size,int herd_index, int Marys, epi_status_t Mary_Type)
{
	// Sample all life characteristics as usual, 
	// but set move and death times such that
	// we have immortal cattle
	// For DIVA validation trials, timeframe of 1 year
	// So set to BIGNUM=100 years
		
	double BIGNUM = 364.0*100;
	
	
	if(herd_index > H) {cout << "You and me, round about midnight.\n"; return false;}
	
	int these[Marys], this_index[Herd_size];
          
    for (int i = 0; i < Herd_size; i++)
    {
    	this_index[i] = (double) i;
    }
          
    gsl_ran_choose(gsl_r, these, Marys, this_index, Herd_size, sizeof (int));
     	
		
	//cout << "Marys : " << Marys << endl;
	// Empty herd if any cows are present
	if(cows[herd_index].size()>0)
	{cows[herd_index].erase(cows[herd_index].begin(),cows[herd_index].end());}
	// Empty infection list
	if(infected[herd_index].size()>0)
	{infected[herd_index].erase(infected[herd_index].begin(),infected[herd_index].end());}
	

	for(int i=0;i < Herd_size; i++)
	{
	double *noo_coo = Get_New_Cow(herd_index,false);
	double key_index = 0.0;
	
	cow_t new_cow;
	
	new_cow.Birth_time = noo_coo[Birth_demo];
	new_cow.Off_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Death_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Infection_time = 0;
	/*cout << cows[herd_index][i].Birth_time << ' ' 
		 << cows[herd_index][i].Off_time << ' '
		 << cows[herd_index][i].Death_time << endl;*/
	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);
		 
	delete [] noo_coo;
	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)
	new_cow.Epi_status = S;
	//new_cow.Epi_status = ((t-new_cow.Birth_time) > 0.0) ? S : (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
	
	//cout << "Status " << new_cow.Epi_status << endl;
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	new_cow.Confirmation_status   = false;
	new_cow.Vaccinated_time = -1;
    new_cow.Control = gsl_ran_binomial (gsl_r,  1.0-target_vaccination_p, 1);
	new_cow.Seeder = false;
	
	
	// Is cow shadow-progressor?

	epi_status_t shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
	new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	//cout << "Shadow: " << shadow << ' ' << new_cow.Epi_status << endl;
	
	// Infect Mary's animals
	for(int k=0; k < Marys; k++)
	 { 
	 if(i==these[k]){new_cow.Epi_status = (epi_status_t) (new_cow.Epi_status + 1);
	 				
	 //cout << "Infection! " << k << ' ' << new_cow.Epi_status << endl;
	 }
	 }
	 
	
	// NOT TESTED YET! FIX FROM initialise_herd above
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
		if(cows[herd_index].find(key_index)==cows[herd_index].end())
	{
	cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
	}
	else{ 
	i--;continue;}
	
	
	//std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	//ret = cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	//if (ret.second==false) {
    //std::cout << "BOVINE COLLISION" << endl;
    // Final safety check for key collision
    //exit(1);

	
	//cows[herd_index][key_index] = new_cow;
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	ret = infected[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	if (ret.second==false) {
    std::cout << "element 'z' already existed";
  	}
  
	}
	
	}
	 
	if(NextMove==herd_index)
	{
	//if(NextMove_time < 0){ cout << "Decapitated Cheerleader" << endl;}
	NextMove_time = (*cows[herd_index].begin()).first;
	}
	
	update_next_move();
	//audit_maps("Update Next Move II");
	calculate_totals();



	return true;
}

bool bTBICBM::initialise_herd_Experiment(int Herd_size)
{
	// Trial animals kept alive to end of trial
	// For DIVA validation trials, timeframe of 1 year
	// So set to BIGNUM=1000 years
		
	int herd_index = 0;
		
	double BIGNUM = 364.0*1000;
		
	// Empty herd if any cows are present
	if(cows[herd_index].size()>0)
	{cows[herd_index].erase(cows[herd_index].begin(),cows[herd_index].end());}
	// Empty infection list
	if(infected[herd_index].size()>0)
	{infected[herd_index].erase(infected[herd_index].begin(),infected[herd_index].end());}
	

	for(int i=0;i < Herd_size; i++)
	{
	double key_index = 0.0;
	
	cow_t new_cow;
	// All animals 6 months old at beginning of experiment
	new_cow.Birth_time = -6*30;
	new_cow.Off_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Death_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Infection_time = 0;
	/*cout << cows[herd_index][i].Birth_time << ' ' 
		 << cows[herd_index][i].Off_time << ' '
		 << cows[herd_index][i].Death_time << endl;*/
	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);

	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)
	new_cow.Epi_status = S;
	//new_cow.Epi_status = ((t-new_cow.Birth_time) > 0.0) ? S : (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
	
	//cout << "Status " << new_cow.Epi_status << endl;
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	new_cow.Confirmation_status   = false;
	new_cow.Vaccinated_time = -1;
    
    
    // 1/2 Herd unvaccinated reactors, 1/4 herd unvaccinated controls
    // So want 3/4 Herd flagged as controls (not vaccinated)
    
    if( i < 3*Herd_size/4 ) 
    {
    new_cow.Control = false;
    }
    else { new_cow.Control = true;}

	// Default value, flip for reactor animals
	new_cow.Seeder = false;
	// Is cow shadow-progressor?

	epi_status_t shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
	new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	//cout << "Shadow: " << shadow << ' ' << new_cow.Epi_status << endl;
	
	// Infect Mary's animals, labelled as controls as they are not vaccinated AND seeders
	// Use post-hoc predictive proportions from SORI within-herd model 
	
	// Animals all identical 
	// So do not need to sample
	// Want Group Structure
	// Seeder animals are Reactors, unvaccinated Control = T & Seeder = T
	// Want half of incontact animals to be vaccinated, half unvaccinated

	if ( i < Herd_size/2)
	{	
	
	
	unsigned int initial_mary[4];
	gsl_ran_multinomial(gsl_r,4,1,Rstate_emp,initial_mary);
	 
	int Mary_Choose;
	for(Mary_Choose=0;Mary_Choose<4;Mary_Choose++)
	 {
	  if(initial_mary[Mary_Choose]==1){break;}
	 }
	 //cout << Mary_Choose << ' ' << initial_mary[Mary_Choose] << endl;
	 
	 new_cow.Epi_status = (epi_status_t) (Mary_Choose);
	 new_cow.Diva_ever=true;
	 new_cow.Control = true;new_cow.Seeder = true;
	 //cout << "Infection! " << k << ' ' << new_cow.Epi_status << endl;

	 }
	  				
	// NOT TESTED YET! FIX FROM initialise_herd above
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
	if(cows[herd_index].find(key_index)==cows[herd_index].end())
	{
	cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
	}
	else{ 
	i--;continue;}
	
	
	//std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	//ret = cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	//if (ret.second==false) {
    //std::cout << "BOVINE COLLISION" << endl;
    // Final safety check for key collision
    //exit(1);

	
	//cows[herd_index][key_index] = new_cow;
	if(new_cow.Epi_status != S && new_cow.Epi_status != SS)
	{
	std::pair<std::map<double, cow_t>::iterator,bool> ret;
  	ret = infected[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
  	if (ret.second==false) {
    std::cout << "element 'z' already existed";
  	}
  
	}
	
	}
	 
	if(NextMove==herd_index)
	{
	//if(NextMove_time < 0){ cout << "Decapitated Cheerleader" << endl;}
	NextMove_time = (*cows[herd_index].begin()).first;
	}
	
	update_next_move();
	//audit_maps("Update Next Move II");
	calculate_totals();


	return true;
}

bool bTBICBM::refresh_herd_Experiment(int Herd_size)
{
	// Trial animals kept alive to end of trial
	// For DIVA validation trials, timeframe of 1 year
	// So set to BIGNUM=1000 years
		
	int herd_index = 0;
		
	double BIGNUM = 364.0*1000;
		
	for(int i=0;i < Herd_size/2; i++)
	{
	double key_index = 0.0;
	
	cow_t new_cow;
	// All animals 6 months old at beginning of experiment
	new_cow.Birth_time = -6*30;
	new_cow.Off_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Death_time = BIGNUM - 1.0 + gsl_ran_flat(gsl_r,-0.25,0.25);
	new_cow.Infection_time = 0;
	/*cout << cows[herd_index][i].Birth_time << ' ' 
		 << cows[herd_index][i].Off_time << ' '
		 << cows[herd_index][i].Death_time << endl;*/
	new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);

	// If birth, animal is susceptible.
	// else chance (pinf_move) animal is infected (S=0,O=1)
	new_cow.Epi_status = S;
	//new_cow.Epi_status = ((t-new_cow.Birth_time) > 0.0) ? S : (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
	
	//cout << "Status " << new_cow.Epi_status << endl;
	new_cow.Epi_stage = 0;
	new_cow.Std_status = NotReactor;	
	new_cow.Svr_status = NotReactor;
	new_cow.Diva_status = NotReactor;
	new_cow.Diva_ever = false;
	new_cow.Confirmation_status   = false;
	new_cow.Vaccinated_time = -1;
    
    // 1/2 Herd unvaccinated reactors, 1/4 herd unvaccinated controls
    // So want 3/4 Herd flagged as controls (not vaccinated)
    
    if( i < 1*Herd_size/4 ) 
    {
    new_cow.Control = false;
    }
    else { new_cow.Control = true;}

	// Default value, flip for reactor animals
	new_cow.Seeder = false;
	// Is cow shadow-progressor?

	epi_status_t shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
	new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	//cout << "Shadow: " << shadow << ' ' << new_cow.Epi_status << endl;
		  				
	// NOT TESTED YET! FIX FROM initialise_herd above
	key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
	if(cows[herd_index].find(key_index)==cows[herd_index].end())
	{
	cows[herd_index].insert(std::pair<double, cow_t>(key_index,new_cow));
	}
	else{ 
	i--;continue;}
	

  
	}
	 
	if(NextMove==herd_index)
	{
	//if(NextMove_time < 0){ cout << "Decapitated Cheerleader" << endl;}
	NextMove_time = (*cows[herd_index].begin()).first;
	}
	
	update_next_move();
	//audit_maps("Update Next Move II");
	calculate_totals();

	return true;
}

bool bTBICBM::set_demo(demo_t model, double Life_Exp, double var1, double births=0.0)
{

	switch(model){
	case LifeExp:
		Demo_model = LifeExp;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp;
		Life_var1[i] = 0.0;
		switch(Seasonal_model)
			{
			case constant:
				{Life_Births[i] = 0.0;}
			break;
			case seasonal:
				{Life_Births[i] = births;}
			break;
		}	
	
		}
		break;
		case LifeFixed:
		Demo_model = LifeFixed;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp;
		Life_var1[i] = 0.0;
		
				switch(Seasonal_model)
			{
			case constant:
				{Life_Births[i] = 0.0;}
			break;
			case seasonal:
				{Life_Births[i] = births;}
			break;
		}	
		
		}
		break;
		case LifeNegBin:
		Demo_model = LifeNegBin;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp;
		Life_var1[i] = var1;
		
				switch(Seasonal_model)
			{
			case constant:
				{Life_Births[i] = 0.0;}
			break;
			case seasonal:
				{Life_Births[i] = births;}
			break;
		}	
		
		}
		break;
		case Darth:
		Demo_model = Darth;
		switch(Seasonal_model)
		{
		case constant:
		for(int i=0;i < H;i++)
			{Life_Births[i] = 0.0;}
		break;
		case seasonal:
		for(int i=0;i < H;i++)
			{Life_Births[i] = DarthMoves[DarthSelecta];}
		break;
		}	
		return true;
		//cout << "Hah! " << Darth << endl;
		break;
		case Experiment:
		Demo_model = Experiment;
		for(int i=0;i < H;i++)
			{Life_Births[i] = 0.0;}
		break;
		
		}

	return true;
}

double* bTBICBM::Get_New_Cow(int i,bool birth)
{

// Age, move_time, death_time 

double* demo_dat = new double[3];
demo_dat[0]=0.0;demo_dat[1]=0.0;demo_dat[2]=0.0;
double L=0.0;
double age=0.0;

//cout << "i = " << i << endl;

switch(Demo_model){
	case LifeExp:
		while(L == 0.0)
		{
	    //cout << Life_Expectation[i] << endl;
		L = gsl_ran_exponential(gsl_r,Life_Expectation[i]);
		}
		age = gsl_ran_exponential(gsl_r,Life_Expectation[i]);
		while(age > L){age = gsl_ran_exponential(gsl_r,Life_Expectation[i]);}
		
		//cout << t << ' ' << L << ' ' << age << ' ' << Life_Expectation[i] << endl;
		
		if(birth)
		{
			age = 0;
			demo_dat[Birth_demo] = t;
		}
		else
		{
			demo_dat[Birth_demo] = t-age;
		}
		
		
		demo_dat[Death_demo] = (t-age)+L;
		demo_dat[Off_demo] = (t-age)+L;	
		
		break;
	case LifeFixed:
		L = Life_Expectation[i];
		
		age = gsl_ran_flat(gsl_r,0,5*Life_Expectation[i]);	
		while(age > L){age = gsl_ran_flat(gsl_r,0,4*Life_Expectation[i]);}
		

		if(birth)
		{
			age = 0;
			demo_dat[Birth_demo] = t;
		}
		else
		{
			demo_dat[Birth_demo] = t-age;
		}
		
		demo_dat[Death_demo] = (t-age)+L;
		demo_dat[Off_demo] = (t-age)+L;	
		//cout << t << ' ' << age << ' ' << L << endl;
		break;
	case LifeNegBin:
		// Add a little random bit to avoid integer life and sort algorithm crashing
	 	L = gsl_ran_negative_binomial(gsl_r, Life_var1[i]/(Life_Expectation[i]+Life_var1[i]),Life_var1[i]) + gsl_ran_flat(gsl_r,-0.25,0.25);
		//cout << Life_Expectation[i] << ' ' << L << endl;
		
		age = gsl_ran_negative_binomial(gsl_r, Life_var1[i]/(Life_Expectation[i]+Life_var1[i]),Life_var1[i]) + gsl_ran_flat(gsl_r,-0.25,0.25);
		
		while(L < age){L = gsl_ran_negative_binomial(gsl_r, Life_var1[i]/(Life_Expectation[i]+Life_var1[i]),Life_var1[i]) + gsl_ran_flat(gsl_r,-0.25,0.25);}
		//cout << L << ' ' << age << endl;

		if(birth)
		{
			age = 0.0;
			demo_dat[Birth_demo] = t;
		}
		else
		{
			demo_dat[Birth_demo] = t-age;
		}
		
		demo_dat[Death_demo] = (t-age)+L;
		demo_dat[Off_demo] = (t-age)+L;		

		break;
	case Darth:
		//cout << DarthSelecta << ' ' << DarthHerdSamples[DarthSelecta] << ' ' 
		//     << DarthHerdOffsets[DarthSelecta] << ' ';
		int thisone = floor(gsl_rng_uniform(gsl_r)*DarthHerdSamples[DarthSelecta]);
		//cout << thisone << endl;
		
		//cout << "Base " << index2(DarthHerdOffsets[DarthSelecta],0,4) << endl;
		if(birth)
		{
			age = Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,0,4)] + gsl_ran_flat(gsl_r,0,0.25);
			
			demo_dat[Birth_demo] = t-age;
		//cout << "Birth " << age << endl;
		}
		else
		{
		age = Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,3,4)] + gsl_ran_flat(gsl_r,0,0.25);
		demo_dat[Birth_demo] = t-age;	
		
		//cout << "Init: " << age << endl;
		}
		
		
		//cout << "This One: " << thisone << " out of " << DarthHerdSamples[DarthSelecta] << ' ' << index2(DarthHerdOffsets[DarthSelecta]+thisone,1,4) << endl; 
		// Add a little random bit to avoid integer life and sort algorithm crashing
		demo_dat[Off_demo] = t + Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,1,4)] + gsl_ran_flat(gsl_r,0,0.25);
		demo_dat[Death_demo] = (t-age) + Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,2,4)] + gsl_ran_flat(gsl_r,0,0.25);
		
		//cout << age << ' ' << demo_dat[Birth_demo] << ' ' << demo_dat[Off_demo] << ' ' << demo_dat[Death_demo] << endl;
		// Adjustment to ensure simulation does not go backwards in time
		// Pathological cases are filtered out, but 
		// can still arise for animals that are slaughtered at birth...
		if(demo_dat[Death_demo] < t){demo_dat[Death_demo]=t;/*cout<<"adjust d" <<endl;*/}
		//if(demo_dat[Off_demo] < t){demo_dat[Off_demo]=t;cout<<"adjust o" <<endl;}
		
		
		
	/*	{cout << "Basic pleb" << endl;
									cout << t << ' ' << demo_dat[Birth_demo] << ' ' << demo_dat[Off_demo] << ' ' << demo_dat[Death_demo] << endl;
									
									cout << DarthSelecta << ' '
									     << Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,0,4)] << ' ' 
										 << Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,1,4)] << ' '
										 << Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,2,4)] << ' '
										 << Darth_cows[index2(DarthHerdOffsets[DarthSelecta]+thisone,3,4)] << endl;
																	}*/

		break;
		}
		


return(demo_dat);
}

bool bTBICBM::set_demo(demo_t model, double *Life_Exp, double *var1, double *births,int Patch)
{
	switch(model){
	case LifeExp:
		Demo_model = LifeExp;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp[i];
		Life_var1[i] = 0.0;
		
		switch(Seasonal_model)
		{
		case constant:
		for(int i=0;i < H;i++)
		{Life_Births[i] = 0.0;}
		break;
		case seasonal:
		for(int i=0;i < H;i++)
		{Life_Births[i] = births[i];}
		break;
		}	
		
		}
		break;
		case LifeFixed:
		Demo_model = LifeFixed;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp[i];
		Life_var1[i] = 0.0;
		
		switch(Seasonal_model)
		{
		case constant:
		for(int i=0;i < H;i++)
		{Life_Births[i] = 0.0;}
		break;
		case seasonal:
		for(int i=0;i < H;i++)
		{Life_Births[i] = births[i];}
		break;
		}	
	
		}
		break;
		case LifeNegBin:
		Demo_model = LifeNegBin;
		for(int i=0;i < H;i++)
		{
		Life_Expectation[i] = Life_Exp[i];
		Life_var1[i] = var1[i];
		}
		
		switch(Seasonal_model)
		{
		case constant:
		for(int i=0;i < H;i++)
			{Life_Births[i] = 0.0;}
		break;
		case seasonal:
		for(int i=0;i < H;i++)
			{Life_Births[i] = births[i];}
		break;
		}	
		
		break;
		case Darth:
		Demo_model = Darth;
		
		switch(Seasonal_model)
		{
		case constant:
		for(int i=0;i < H;i++)
			{Life_Births[i] = 0.0;}
		break;
		case seasonal:
		for(int i=0;i < H;i++)
			{Life_Births[i] = DarthMoves[DarthSelecta];}
		break;
		}
			
		break;
		case Experiment:
		Demo_model = Experiment;
		
		}
	
	
		
		
	return true;	
}

double bTBICBM::return_time()
{
	return t;
}

void bTBICBM::set_forcing(bool setflag)
{
	if(setflag)
	{
	Seasonal_model = seasonal;
	}
	else
	{
	Seasonal_model = constant;
	}
}

bool bTBICBM::disease_free(int i=0)
{	

for(int i=0; i < H; i++)
		{
			  map<double,cow_t>::iterator iter;   
  			  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  			  {
     			if((iter->second).Epi_status != S && (iter->second).Epi_status != V1 && (iter->second).Epi_status != V2 
     			&& (iter->second).Epi_status != SS && (iter->second).Epi_status != SV1 && (iter->second).Epi_status != SV2)
     			{return false;}
     			}
			
		}

	return(true);
	
	//return((Rtot[i] == 0 && Itot[i] == 0 && Otot[i] == 0 &&  RVtot[i] == 0 && IVtot[i] == 0 && OV1tot[i] == 0 && OV2tot[i] == 0));

}

int bTBICBM::burden(int i=0)
{
int burr=0;
for(int i=0; i < H; i++)
		{
		
		
			map<double,cow_t>::iterator iter;   
  			  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  			  {
     			
     			if((iter->second).Epi_status != S && (iter->second).Epi_status != V1 && (iter->second).Epi_status != V2
     			&& (iter->second).Epi_status != SS && (iter->second).Epi_status != SV1 && (iter->second).Epi_status != SV2)
     			{burr++;}
     			
     			}
			
		}

	return(burr);
	
/*
return((Rtot[i] + DRtot[i] + OFtot[i] +  OStot[i] + DOFtot[i] + DOStot[i] + I[i] + DI[i]));
*/
}


double bTBICBM::run(double runterval,double time_slice,record_t record)
{

	simtime = 0;
	bini = 0;
	told = t; // use to measure time from start of current simulation run
	
	binterval = time_slice;
	
	recordingstatus = record;
	// If saving date write out initial timepoint
	if(full_save)
	{save_data();}
	
	//cout << "Entrant" << endl;
	audit_maps("Entrance");
	
	calculate_totals();
	
	
	recalculate_weights();
	//cout << "Perishing" << endl;

	update_next_move();
	

	//cout << "Initial time: " << t << " Run for " << runterval << endl;

	//Reset dt at each run
	dt = 0.1;
	//cout << "Herds:\n";
	// Flag to elegantly advance when rate = 0
	bool no_event = false;
	while (simtime < runterval)
	{
		
		/*
		while(cows[0].size() < TargetHerdSize)
		{
		 add_bovine(0);
		}
		*/
		
		
		//cout << t << ' ' << cows[0].size() << endl;
		//cout << t << ' ' << rate << endl;
			//audit_maps("looper");
		calculate_totals();
		//recalculate_rate(0);
		//recalculate_weights();
		
	     //cout << "Time " << t << " Rate: " << rate << " dt: " << dt 
	     //	   << " Next Special: " << NextSpecial_time
	     //	   << " Next Move: " << NextMove_time << endl;
	     	   
		//cout << burden(0) << ' ' << cows[0].size() << endl;

		double dt = 0;
		// Select a number between 0 and 1
		// Inter-event time by sampling exponential distribution
		do{x = gsl_rng_uniform(gsl_r);}while(x == 0);
		if(x<0 || x>1.0){cout << "Ranx: " << x << endl;}
				
		// Sanity check for division by zero
		// If rate = 0 (all infected/disease extinct) set dt to 
		// jump to end time for simulation
		// no_event flag allows us to step through simulation
		// until all special events are carried out
		// (which may lead to rate becoming non-zero and 
		// returning to simulation loop)
		if(rate == 0){ dt = runterval-simtime;no_event=true;/*cout << "dt fudged: " << dt;*/}
		else{dt = - (long double) log(x) / rate;no_event=false;}
		// Select a number between 0 and 1
		// Use to choose next event
		
		//cout << t << ' ' << "Next Markov at " << (t + dt) << endl;
		//cout << "Time " << t << " Rate: " << rate << " dt: " << dt << endl;	
		// No event before next bin or next special event
		
		if((simtime + dt) > binterval*bini && (NextSpecial_time > (told + binterval*bini)))
		{
			//cout << "Save Data" << endl;
			//cout << "Time " << t;
			//cout << " Rate: " << rate;
			//cout << " dt: " << dt << endl;

			t = t + (binterval*bini-simtime);
			simtime=binterval*bini;
			incro_cows();
			if(full_save){save_data();}
			//audit_maps("Do Nothing");
			bini++;
		}
		// else if no event before next special event
		else if((t + dt) > NextSpecial_time)
		{
		//cout << "Special, next move: " << NextSpecial_time << endl;
//cout << "Time " << t << " Rate: " << rate << " dt: " << dt << endl;
		
		// Keep removing animals until next Move > current time
		if(Special_is_move)
		{
		t = NextMove_time;
		simtime = t-told;
		
		// If we detect infection at slaughterhouse halt and return 
		// Simulation time since entry
		//cout << Cuniq_id << " Manage" << endl;
		// #pragma omp critical(dataout)
			//	{
			//	cout << Cuniq_id << " Management" << endl;
			//	}
		if(do_management()){/*cout << "End Run Slaughter";*/return (simtime);}
		//audit_maps("Do Management");
		audit_time("Management");
		//cout << "Bellisimo " << endl;
		}
		else
		{
//cout << "Vaccinate at  " << t << endl;
		 vaccinate();
		 //audit_maps("Vaccinate");
		}
		
		}
		
		else // Do event
		{
			
			//cout << "Do Event" << endl;
			x = gsl_rng_uniform(gsl_r)*rate;
				
			//cout << "Time " << t << " Rate: " << rate << " dt: " << dt << endl;
			//cout << "Rate1: " << Herdrates[0] << " Rate2: " << Herdrates[1] << endl;	
			// Advance time
				
			simtime = simtime + dt;
			t       = t       + dt;
			
			if(!no_event)
			{
			double p = 0.0;
			for(int l=0; l < H; l++)
			{
				p += Herdrates[l];
				if(x < p)
					{   
						//cout << x << ' ' << l << ' ' << p << ' ' << rate << ' ' << Herdrates[l] << endl;
						
						do_event(l);
						audit_time("Event");
						//audit_maps("Do Event");
						//#pragma omp critical(dataout)
						//{
						//cout << Cuniq_id << " Event" << endl;
						//}
						break;
					}
			
			}
			}

		
		
		}
	
		
		if((simtime) > binterval * bini)
	    {

			//cout << "Marmalade" << endl;
			//cout << t << ' ' << bini << ' ' << binterval*bini << endl;
			// Total cumulative "Infection" in herd
			//cout << "Save Data" << endl;
			//cout << "Time " << t << " Rate: " << rate << " dt: " << dt << endl;	

			incro_cows();
			if(full_save){save_data();cout<<"Saving Data " << t << endl;}
			bini++;
		//audit_maps("Mark Time");
			
		}
		
	}
	//cout << "End Run Natural " << t << endl;
	//cout << "\n";	
	
	return 0;
	
}

void bTBICBM::incro_cows()
{
	
	
	for(int i=0; i < H; i++)
    {
		prob_on_trans += (Otot[i] + Rtot[i] + Itot[i] + OV1tot[i] + OV2tot[i] + RVtot[i] + IVtot[i]
						+ SOtot[i] + SRtot[i] + SItot[i] + SOV1tot[i] + SOV2tot[i] + SRVtot[i] + SIVtot[i]);
	//cout << Otot[i] << ' ' << Rtot[i] << ' ' << Itot[i]
	//	 << OV1tot[i] << ' ' << OV2tot[i] << ' ' << RVtot[i] << ' ' << IVtot[i] << endl;
	}
	//cout << "Aggrevation: " << prob_on_trans << endl;
	
}

void bTBICBM::save_data()
{

	// Record state variables
	//cout << "Agony!" << endl;
	
	
	if(debug_save)
	{
	//cout << "Bunions!" << endl;
	
		for(int i=0; i < H; i++)
		{
			int Stotal=0;int SStotal=0;
			int Ototal=0;int SOtotal=0;
			int Rtotal=0;int SRtotal=0;
			int Itotal=0;int SItotal=0;
			int V1total=0;int SV1total=0;
			int V2total=0;int SV2total=0;
			int OV1total=0;int SOV1total=0;
			int OV2total=0;int SOV2total=0;
			int RVtotal=0;int SRVtotal=0;
			int IVtotal=0;int SIVtotal=0;
			
			Agefile[i] << t; 
			Statusfile[i] << t;
			Lifefile[i] << t;
			
			Lifefile[i] << ' ' << (cows[i].size());
			
			double AgeAccum=0.0;
			
			  	map<double,cow_t>::iterator iter;   
  			  for( iter = cows[i].begin(); iter != cows[i].end(); iter++ ) 
  			  {
     			
     			
     			//Lifefile[i] << ' ' << (iter->second.Death_time-iter->second.Birth_time);
     			//Lifefile[i] << ' ' << (cows[i][j].Death_time);
     			AgeAccum += (t - (iter->second).Birth_time);
     			Agefile[i] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
     			
     			switch((iter->second).Epi_status)
     			{
     			case S:
  					Stotal++;
  				break;
  				case O:
  					Ototal++;
  				break;
  				case R:
  					Rtotal++;
  				break;
  				case I:
  					Itotal++;
  				break;
  				case V1:
  					V1total++;
  				break;
  				case V2:
  					V2total++;
  				break;
  				case OV1:
  					OV1total++;
  				break;
  				case OV2:
  					OV2total++;
  				break;
  				case RV:
  					RVtotal++;
  				break;
  				case IV:
  					IVtotal++;
  				break;
  	      		case SS:
  					SStotal++;
  				break;
  				case SO:
  					SOtotal++;
  				break;
  				case SR:
  					SRtotal++;
  				break;
  				case SI:
  					SItotal++;
  				break;
  				case SV1:
  					SV1total++;
  				break;
  				case SV2:
  					SV2total++;
  				break;
  				case SOV1:
  					SOV1total++;
  				break;
  				case SOV2:
  					SOV2total++;
  				break;
  				case SRV:
  					SRVtotal++;
  				break;
  				case SIV:
  					SIVtotal++;
  				break;
    			
     			}
     			
     			//Statusfile[i] << ' ' << (iter->second).Epi_status;
     			}
     			
			Lifefile[i] << ' ' << AgeAccum / (cows[i].size());
			Lifefile[i] << ' ' << Stotal;
			Lifefile[i] << ' ' << Ototal;
			Lifefile[i] << ' ' << Rtotal;
			Lifefile[i] << ' ' << Itotal;
			Lifefile[i] << ' ' << V1total;
			Lifefile[i] << ' ' << V2total;
			Lifefile[i] << ' ' << OV1total;
			Lifefile[i] << ' ' << OV2total;
			Lifefile[i] << ' ' << RVtotal;
			Lifefile[i] << ' ' << IVtotal;
			Lifefile[i] << ' ' << SStotal;
			Lifefile[i] << ' ' << SOtotal;
			Lifefile[i] << ' ' << SRtotal;
			Lifefile[i] << ' ' << SItotal;
			Lifefile[i] << ' ' << SV1total;
			Lifefile[i] << ' ' << SV2total;
			Lifefile[i] << ' ' << SOV1total;
			Lifefile[i] << ' ' << SOV2total;
			Lifefile[i] << ' ' << SRVtotal;
			Lifefile[i] << ' ' << SIVtotal;
			
			//Agefile[i] << endl;
			Statusfile[i] << endl;
			Lifefile[i] << endl;
			
		}
		
	}
	
}



bool bTBICBM::slaughter_house_test(cow_t cow)
{
double this_p = 0.0;

//cout << "Testing at slaughter: " << cow.Epi_status << endl;
//cout << "Retro: " << retro << " " << (cow.Death_time - cow.Off_time) << endl;
// If not a move off to slaughter, no test
// 90 day cut-off gives ~ 20% of off movements tested
if(!retro && (cow.Death_time - cow.Off_time) > 90)
{
 return false;
}
else if(cow.Epi_status == S || cow.Epi_status == SS || cow.Epi_status == V1 || cow.Epi_status == V2 || cow.Epi_status == SV1 || cow.Epi_status == SV2)
	{
	return(false);
	}
else
	{
	double rnd_rain = gsl_rng_uniform(gsl_r);
	
	int age_class=floor((t-cow.Birth_time)/200.0);
	
	if(age_class>20){age_class=20;}

	
	if(cow.Epi_status == SI || cow.Epi_status == I || cow.Epi_status == IV || cow.Epi_status == SIV)
	{
	//this_p=p_confirm[0];
	this_p=confirm_by_age[age_class];
	//cout << "I " << this_p << ' ' <<  slaughter_sensitivity << endl;
	}
	else if(cow.Epi_status == SR || cow.Epi_status == R || cow.Epi_status == RV || cow.Epi_status == SRV)
	{
	//this_p=p_confirm[1];
	this_p=confirm_by_age[age_class];
	//cout << "R " << this_p << ' ' <<  slaughter_sensitivity << endl;
	}
	else
	{
	//this_p=p_confirm[2];
	this_p=0.0;
	//cout << "O " << this_p << ' ' <<  slaughter_sensitivity << endl;
	}
		
	return((gsl_rng_uniform(gsl_r) <= this_p*slaughter_sensitivity));
	}
}
	
	
void bTBICBM::do_event(int h)
{


x = gsl_rng_uniform(gsl_r)*Herdrates[h];
				
double p = 0.0;
double y = 0.0;

//double *demo_dat = Get_New_Cow(h,true);

p = Life_Births[h];

//cout << "EVENT! Time: " << t << ' ' << h << ' ' << x << ' ' << Herdrates[h] << ' ' << Life_Births[h] << ' ' << p << endl;
		
// If Birth
if(x < p)
{
//cout << "Sanity!" << endl;
if(cows[h].size() < TargetHerdSize)
{ 

do{

double *demo_dat = Get_New_Cow(h,true);

//cout << "Birth! " << t << ' ' << h << ' ' << x << ' ' << p << endl;	
cow_t new_cow;
	
new_cow.Death_time = demo_dat[Death_demo];
new_cow.Off_time = demo_dat[Off_demo];
new_cow.Birth_time = demo_dat[Birth_demo];
new_cow.Infection_time = 0.0;
new_cow.Vaccinated_time = -1;
new_cow.Control = gsl_ran_binomial (gsl_r,  target_vaccination_p, 1);
new_cow.Seeder = false;
new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);

delete [] demo_dat;

// If birth, animal is susceptible.
// else chance (pinf_move) animal is infected (S=0,O=1)

if((t-new_cow.Birth_time) > 0.0)
{
new_cow.Epi_status = (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
if(new_cow.Epi_status != S)
{
//Truncate Infection time to time onto herd
new_cow.Infection_time = t;

}
}
else
{
new_cow.Epi_status = S;
}

// Is cow shadow-progressor?
epi_status_t shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	

new_cow.Epi_stage = 0;
new_cow.Std_status = NotReactor;	
new_cow.Svr_status = NotReactor;
new_cow.Diva_status = NotReactor;
new_cow.Diva_ever = false;
new_cow.Confirmation_status = false;

double key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;

//cout << "Death: " << new_cow.Death_time << " Off Time: "  << new_cow.Off_time << endl;
//cout << "Key Index: " << key_index << endl;

//cout << "Old herd: " << cows[h].size() << endl;

	    // check that key does not exist in herd, if it does then resample
	if(cows[h].find(key_index)==cows[h].end())
	{
	cows[h].insert(std::pair<double, cow_t>(key_index,new_cow));
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	infected[h].insert(std::pair<double, cow_t>(new_cow.uniq_id,new_cow));
	}
	switch(new_cow.Epi_status)
{
case S:
	Stot[h]++;
	update_weight(h, cows[h].insert(std::pair<double, cow_t>(key_index,new_cow)).first);
break;
case O:
	Otot[h]++;
	infected[h].insert(std::pair<double, cow_t>(key_index,new_cow));
	recalculate_weights(NextMove);
break;
case SS:
	SStot[h]++;
	update_weight(h, cows[h].insert(std::pair<double, cow_t>(key_index,new_cow)).first);
break;
case SO:
	SOtot[h]++;
	infected[h].insert(std::pair<double, cow_t>(key_index,new_cow));
	recalculate_weights(NextMove);
break;
}
	
	break;
	}
	else{continue;}
	}while(true);



//cout << "Live Birth!" << cows[h].size() << endl;

update_next_move();
}

}
else
{

// Which cow
//cout << t << ' ' << x << ' ' << p << ' ' << Herdrates[h] << endl;
map<double,cow_t>::iterator iter;   
for( iter = cows[h].begin(); iter != cows[h].end(); iter++ ) 
{
				p += (iter->second).rate;
				
				if(x < p)
					{
					 //cout << t << ' ' << x << ' ' << p << ' ' << Herdrates[h] << " TP" << endl;

					 switch((iter->second).Epi_status)
  					{
  						case S:
  							// If animal is under 20 months 
  							// switch to "fast" progression
  							//cout << "S FAST event" << endl;
  							
  							// Slow progression
  							//cout << "S event" << endl;
  							(iter->second).Epi_status = O;
  							(iter->second).Epi_stage = 0;
  							(iter->second).Infection_time = t;
  							infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							Otot[h]++;Stot[h]--;
  							update_weight(h, iter);
  							
  						break;
  						case O:
  							//cout << "O event" << endl;
  							
  							if((iter->second).Epi_stage < TO)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=R;
  							Rtot[h]++;Otot[h]--;recalculate_weights(h);
  							}
  							
  						break;
  						case R:
  							//cout << "R event" << endl;
  							if((iter->second).Epi_stage < TR)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=I;
  							Itot[h]++;Rtot[h]--;recalculate_weights(h);
  							}
  						break;
  						case I:
  							//cout << "Angry Inch: " << cows[h][l].rate << endl;
  						break;
  						case V1:
  							 //cout << "V1 event" << endl;
  							 y = gsl_rng_uniform(gsl_r)*(iter->second).rate;
  							if(y < gV1)
  							{
  							 if((iter->second).Epi_stage < TV1)
  								{(iter->second).Epi_stage++;}
  							else
  								{
  								(iter->second).Epi_stage=0;(iter->second).Epi_status=V2;
  								V1tot[h]--;V2tot[h]++;update_weight(h, iter);
  								}
  							}
  							else
  							{
  							(iter->second).Epi_status=OV1;(iter->second).Epi_stage=0;
  							OV1tot[h]++;V1tot[h]--;
  							infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							recalculate_weights(h);
  							}
  						break;
  						case V2:
  							//cout << "V2 event" << endl;			
  							(iter->second).Epi_status=OV2;(iter->second).Epi_stage=0;
  							OV2tot[h]++;V2tot[h]--;
  							infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							recalculate_weights(h);
  						break;
  						case OV1:
  							//cout << "OV1 event" << endl;
  							if((iter->second).Epi_stage < TOV1)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=RV;
  							OV1tot[h]--;RVtot[h]++;recalculate_weights(h);
  							}
  						break;
  						case OV2:
  							//cout << "OV2 event" << endl;
  							if((iter->second).Epi_stage < TOV2)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=RV;
  							OV2tot[h]--;RVtot[h]++;recalculate_weights(h);
  							}
  						break;
  						case RV:
  							//cout << "RV event" << endl;
  							if((iter->second).Epi_stage < TRV)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=IV;
  							RVtot[h]--;IVtot[h]++;recalculate_weights(h);
  							}
  						break;
  						case IV:
  							//cout << "Angry Inch: " << cows[h][l].rate << endl;
  						break;
  						case SS:
  							//cout << "SS event" << endl;
  							// If animal is under 24 months 
  							// switch to "fast" progression
  							
  							
  							if((((t - (iter->second).Birth_time))/30.0) < delta_age)
  							{
							(iter->second).Epi_status = SO;
  							(iter->second).Epi_stage = 0;
  							(iter->second).Infection_time = t;
  							SOtot[h]++;SStot[h]--;
  							//cout << "S FAST event" << endl;
  							infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							
  							update_weight(h, iter);
  							}
  							else
  							{// Slow progression
  							//cout << "S event" << endl;
  							(iter->second).Epi_status = O;
  							(iter->second).Epi_stage = 0;
  							(iter->second).Infection_time = t;
  							infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							Otot[h]++;SStot[h]--;
  							update_weight(h, iter);
  							}
  					
  						break;
  						case SO:
  							//cout << "SO event" << endl;
  							
  							if((iter->second).Epi_stage < TSO)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=SR;
  							SRtot[h]++;SOtot[h]--;recalculate_weights(h);
  							}
  							
  						break;
  						case SR:
  							//cout << "SR event" << endl;
  							if((iter->second).Epi_stage < TSR)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=SI;
  							SItot[h]++;SRtot[h]--;recalculate_weights(h);
  							}
  						break;
  						case SI:
  							//cout << "Angry Inch: " << cows[h][l].rate << endl;
  						break;
  						case SV1:
  							// cout << "SV1 event" << endl;
  							 y = gsl_rng_uniform(gsl_r)*(iter->second).rate;
  							if(y < gSV1)
  							{
  							 if((iter->second).Epi_stage < TSV1)
  								{(iter->second).Epi_stage++;}
  							else
  								{
  								(iter->second).Epi_stage=0;(iter->second).Epi_status=SV2;
  								SV1tot[h]--;SV2tot[h]++;update_weight(h, iter);
  								}
  							}
  							else
  							{
  							  // If animal is under 24 months remain in fast progression
  							  // Otherwise switch over to slow progression
  							  if((((t - (iter->second).Birth_time))/30.0) < delta_age)
  							  {
  								(iter->second).Epi_status=SOV1;(iter->second).Epi_stage=0;
  								SOV1tot[h]++;SV1tot[h]--;
  								infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  								recalculate_weights(h);
  							  }
  							  else
  							  {
  							  //Slow Progression, vaccinate
  							  (iter->second).Epi_status=OV1;(iter->second).Epi_stage=0;
  							  OV1tot[h]++;SV1tot[h]--;
  							  infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  							  recalculate_weights(h);
  							  }
  							  
  							
  							
  							}
  						break;
  						case SV2:
  							//cout << "SV2 event" << endl;
  						   // Infection
  							
  							// If animal is under 24 months remain in fast progression
  							  // Otherwise switch over to slow progression
  							  if((((t - (iter->second).Birth_time))/30.0) < delta_age)
  							  {
  								(iter->second).Epi_status=SOV2;(iter->second).Epi_stage=0;
  								SOV2tot[h]++;SV2tot[h]--;
  								infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  								recalculate_weights(h);
  							  }
  							  else
  							  {
  							   //Slow Progression, vaccinate
  							   (iter->second).Epi_status=OV2;(iter->second).Epi_stage=0;
  								OV2tot[h]++;SV2tot[h]--;
  								infected[h].insert(std::pair<double, cow_t>(iter->first, iter->second));
  								recalculate_weights(h);
  							  }
  						
  						break;
  						case SOV1:
  							//cout << "OV1 event" << endl;
  							if((iter->second).Epi_stage < TSOV1)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=SRV;
  							SOV1tot[h]--;SRVtot[h]++;recalculate_weights(h);
  							}
  						break;
  						case SOV2:
  							//cout << "OV2 event" << endl;
  							if((iter->second).Epi_stage < TSOV2)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=SRV;
  							SOV2tot[h]++;SRVtot[h]--;recalculate_weights(h);
  							}
  						break;
  						case SRV:
  							//cout << "RV event" << endl;
  							if((iter->second).Epi_stage < TSRV)
  							{(iter->second).Epi_stage++;}
  							else
  							{
  							(iter->second).Epi_stage=0;(iter->second).Epi_status=SIV;
  							SRVtot[h]--;SIVtot[h]++;recalculate_weights(h);
  							}
  						break;
  						case SIV:
  							//cout << "Angry Inch: " << cows[h][l].rate << endl;
  						break;
  						default:
  						cout << "Repeats in my ear: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(2);
  						break; 
  						}
  						
  						break;
					}

}
			
}


}

	void bTBICBM::set_time(double set_clock=0.0)
	{
		
		t = set_clock;
		told = set_clock;
		bini = 0;
		
	}
	
	
	void bTBICBM::set_inout(string input, string output)
	{
		paramdirectory = input;
		datadirectory  = output;
		    
		// Close and re-open output files

			for(int i=0; i < H; i++)
			{
				if(debug_save)
				{
				if(Statusfile[i]){Statusfile[i].close();}
				if(Agefile[i]){Agefile[i].close();}
				if(Lifefile[i]){Lifefile[i].close();}
				}
				if(disclose[i]){disclose[i].close();}
			}
		
		
			for(int i=0; i < H; i++)
			{
				if(debug_save)
				{
				Statusfile[i].open((datadirectory + filename("StatusH",i)).c_str());
				Agefile[i].open((datadirectory + filename("AgeH",i)).c_str());
				Lifefile[i].open((datadirectory + filename("LifeH",i)).c_str());
				}
				disclose[i].open((datadirectory + filename("DiscloseH",i)).c_str());
				disclose[i] << "Time " << "R " << "IR " << "TrueR " << "Burden " << "Tested " << "Confirmed " << "StdS " << "StdOF " << "StdOS " << "StdR " << "StdI " << "SvrS " << "SvrOF " << "SvrOS " << "SvrR " << "SvrI " << "StdDS " << "StdDOF " << "StdDOS " << "StdDR " << "StdDI " << "SvrDS " << "SvrDOF " << "SvrDOS " << "SvrDR " << "SvrDI " << "TestS " << "TestFO " << "TestSO " << "TestR " << "TestI " << "Replicate " << "Removed " << "PTI " << "HerdSize " << "TestType " << "OnHerd" << endl;
			}
	}
	
	
	void bTBICBM::set_transmission(double beta_o, double beta_r,double beta_i, double beta_ov1, double beta_ov2,double beta_rv, double beta_iv,double beta_So, double beta_Sr,double beta_Si, double beta_Sov1, double beta_Sov2,double beta_Srv, double beta_Siv, double set_scaling,double xinf_set,double setVacc_Eff = 0.0,double set_Pshadow=0.0)
	{
		
		
		
		contactO= beta_o;
		contactR = beta_r;
		contactI = beta_i;
		contactOV1 = beta_ov1;
		contactOV2 = beta_ov2;
		contactRV = beta_rv;
		contactIV = beta_iv;
		
		contactSO= beta_So;
		contactSR = beta_Sr;
		contactSI = beta_Si;
		contactSOV1 = beta_Sov1;
		contactSOV2 = beta_Sov2;
		contactSRV = beta_Srv;
		contactSIV = beta_Siv;
		
		q = set_scaling;
		xinf = xinf_set;
		Vacc_Eff = setVacc_Eff;
		p_shadow = set_Pshadow;
	}
	
	
	bool bTBICBM::IsLesioned(map<double,cow_t>::iterator iter)
	{
	double rnd_rain = gsl_rng_uniform(gsl_r);
	double this_p = 0.0;
	
	// Confirmation probability dependant on age only
	// Occult animals treated as susceptibles
	
	int age_class=floor((t-iter->second.Birth_time)/200.0);
	
	if(age_class>20){age_class=20;}
	
	if((iter->second).Epi_status == SI || (iter->second).Epi_status == I || (iter->second).Epi_status == IV || (iter->second).Epi_status == SIV)
	{
	//this_p=p_confirm[0];
	this_p=confirm_by_age[age_class];
	
	}
	else if((iter->second).Epi_status == SR || (iter->second).Epi_status == R || (iter->second).Epi_status == RV || (iter->second).Epi_status == SRV)
	{
    //this_p=p_confirm[1];
	this_p=confirm_by_age[age_class];
	}
	else
	{
	//this_p=p_confirm[2];
	this_p=0.0;
	}
	
	//p_confirm = 1.0-0.9*exp(-0.5*exp(-0.03*(t-(iter->second).Infection_time)));
	//cout << "Prob Lesions: " << this_p << endl;
	
	if(rnd_rain < this_p){(iter->second).Confirmation_status=true;/*cout << "Lesioned!" << endl*/;return true;}
	else
	{return false;}
			
	}
	

	vector<int> bTBICBM::Test_Protocol_Per_Animal(int this_herd,int PTI,test_t whole_herd_test,double eligible_vacc=0.0)
	{
		
	vector<int> results(41,0);
		
	
		   		
		//cout << "DIVA! " << DIVA << endl;
		
		current_testtype = whole_herd_test;
		
		int animals_to_test = 0;
		bool* elegible      = new bool[cows[this_herd].size()];
		int i=0;
		if(whole_herd_test==all)
		{
		
			animals_to_test = (cows[this_herd].size());
				
				map<double,cow_t>::iterator iter; 
  				for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
				{
				 switch((iter->second).Epi_status)
  					{
  						case S:
  								elegible[i] = true;
								i++;	
  						break;
  						case O:
								elegible[i] = true;
								i++;	
  						break;
  						case R:
								elegible[i] = true;
								i++;	
  						break;
  						case I:
  								elegible[i] = true;
								i++;	
  						break;
  						case V1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								i++;
								animals_to_test--;	
								}	
  						break;
  						case V2:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case OV1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case OV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case RV:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case IV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}		
  						break;
  						case SS:
  								elegible[i] = true;
								i++;	
  						break;
  						case SO:
  								elegible[i] = true;
								i++;	
  						break;
  						case SR:
  								elegible[i] = true;
								i++;	
  						break;
  						case SI:
  								elegible[i] = true;
								i++;	
  						break;
  						case SV1:
 								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case SV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case SOV1:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  		
  						break;
  						case SOV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case SRV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						case SIV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;	
								}	
  						break;
  						default:
  						cout << "I did it my way: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  						}
				//cout << "Eligibility: " << (iter->second).Epi_status << ' ' << elegible[i-1] << endl; 
				
				}
		}
		else
		{
			
			// Choose Routine tests to be WHT/RHT acording to proportions that initiate breakdowns in 
			// Study Population 
			// PTI 1: 52 RHT, 3235 WHT, SLH 781 p(WHT) = 3235/(3235+52) = 0.9841801
		    // PTI 2: 443 RHT, 478 WHT, SLH 282 p(WHT) = 443/(478+443) = 0.48
			// PTI 4: 555 RHT, 111 WHT, SLH 268 p(WHT) = 111/(555+111) = 0.1666667
			
			switch(PTI)
			{		
					// PTI 1
					case 1:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.9841801){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
					// PTI 2
					case 2:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.48){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
					//case 2:
					//if(gsl_ran_flat(gsl_r,0,1) <= 0.0){whole_herd_test = wht;}
					//else{whole_herd_test = rht;}
					break;
					// PTI 4
					case 4:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.1666667){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
			}
			
			
			map<double,cow_t>::iterator iter; 
			for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
			{
					
			switch(whole_herd_test)
			{
				case wht:
					
					if((t-(iter->second).Birth_time) > (6*7)){elegible[i] = true;animals_to_test++;}
     						else{ elegible[i] = false;}
				break;
				case rht:
					
     				if((t-(iter->second).Birth_time) > (2*364))
     						{elegible[i] = true;animals_to_test++;}
     					else{ elegible[i] = false;}
     			break;
			}
			// Remove eligibility for relevant vaccinates
			
			switch((iter->second).Epi_status)
  			{
 
  						case V1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								}
								else
								{
								elegible[i] = elegible[i] && false;	
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case V2:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case OV1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case OV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case RV:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						case IV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}		
  						break;
  						
  						case SV1:
 								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
							
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						case SOV1:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  		
  						break;
  						case SOV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SRV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
							
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SIV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						}
			
			i++;
			
			}
					
					
			}
			
		
	/*
	cout << "Testing " 
	<< animals_to_test/(double) cows[this_herd].size() 
	<< "% : " 
	<< animals_to_test 
	<< " out of " << cows[this_herd].size() << endl;
	*/
	
    //animals_to_test = (cows[this_herd].size());
     		
		int testS=0;
		int testO=0;
		int testOV1=0;
		int testOV2=0;
		int testV1=0;
		int testV2=0;
		int testR=0;
		int testI=0;
		int testRV=0;
		int testIV=0;
				
		int std_detected_I = 0;	    int severe_detected_I = 0;
		int std_detected_R = 0;	    int severe_detected_R = 0;
		int std_detected_O = 0;		int severe_detected_O = 0;
		int std_detected_OV1 = 0;	int severe_detected_OV1 = 0;
		int std_detected_V1 = 0;	int severe_detected_V1 = 0;
		int std_detected_OV2 = 0;	int severe_detected_OV2 = 0;
		int std_detected_V2 = 0;	int severe_detected_V2 = 0;
		int std_detected_IV = 0;	int severe_detected_IV = 0;
		int std_detected_RV = 0;	int severe_detected_RV = 0;
		
		int std_detected_S = 0;	    int severe_detected_S = 0;
		
		int vl_stddetected=0;
		
		int std_DIVAnegate=0;	int severe_DIVAnegate=0;
		int std_DIVAtests=0;	int severe_DIVAtests=0;

		i = 0;int tested=0;
		
		double rnd_rain;
		
		map<double,cow_t>::iterator iter; 
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		//cout << i << ' ' << test_ind << ' ' << daPick[test_ind] << endl;

  
		 if(!elegible[i])
		 {
		 i++;
		 }
		 else
		 {
		  tested++;i++;
		  
		double detect  = 0.0;
		double detectS = 0.0;
		double t_from_v = 0.0;
		int this_one = 0;  
		  
		  
		switch((iter->second).Epi_status)
		{
		case S:
		testS++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		if(rnd_rain <= specificity[SICCT]){std_detected_S++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= specificity[SICCT_S]){severe_detected_S++;(iter->second).Svr_status=Reactor;}
		
		break;
		case O:
		testO++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= (specificity[SICCT]))		{std_detected_O++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= (specificity[SICCT_S]))	{severe_detected_O++;(iter->second).Svr_status=Reactor;}
		
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		case R:
		testR++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_R++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_R++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case I:
		testI++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_I++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_I++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case V1:

		testV1++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		
		if(rnd_rain <= vaccSpec[this_one]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}

		// If we get reacting vaccinate, perform DIVA test to try and negate
		 rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{		
			std_DIVAtests++;
			if(rnd_rain > specificity[DIVA]){std_detected_V1--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{
			severe_DIVAtests++;
			if(rnd_rain > specificity[DIVA]){severe_detected_V1--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
			

		break;
		case V2:
		testV2++;

		rnd_rain = gsl_rng_uniform(gsl_r);
		
		//if(rnd_rain <= specificity[SICCT_V2]){std_detected_V2++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V2]){severe_detected_V2++;(iter->second).Svr_status=Reactor;}
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		
		if(rnd_rain <= vaccSpec[this_one]){std_detected_V2++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_V2++;(iter->second).Svr_status=Reactor;}

		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
			if(rnd_rain > specificity[DIVA]){std_detected_V2--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{
		severe_DIVAtests++;
		if(rnd_rain > specificity[DIVA]){severe_detected_V2--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		break;
		case OV1:
		testOV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= (specificity[SICCT_V1]))		{std_detected_OV1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= (specificity[SICCT_V1]))	{severe_detected_OV1++;(iter->second).Svr_status=Reactor;}
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		
		if(rnd_rain <= vaccSpec[this_one]){std_detected_OV1++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_OV1++;(iter->second).Svr_status=Reactor;}

		
		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_OV1--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_OV1--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		

		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case OV2:
		testOV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= (specificity[SICCT_V2]))		{std_detected_OV2++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= (specificity[SICCT_V2]))	{severe_detected_OV2++;(iter->second).Svr_status=Reactor;}
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		if(rnd_rain <= vaccSpec[this_one]){std_detected_OV2++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_OV2++;(iter->second).Svr_status=Reactor;}

		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_OV2--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_OV2--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		

		break;
		case RV:
		testRV++;
		
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;

		detect  = (vaccSpec[this_one] > sensitivity[SICCT]) ? vaccSpec[this_one] : sensitivity[SICCT];
		detectS = (vaccSpecS[this_one] > sensitivity[SICCT_S]) ? vaccSpecS[this_one] : sensitivity[SICCT_S];
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;		
		
		if(rnd_rain <= detect){std_detected_RV++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= detectS){severe_detected_RV++;(iter->second).Svr_status=Reactor;}
	
		//if(rnd_rain <= sensitivity[SICCT]){std_detected_RV++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_RV++;(iter->second).Svr_status=Reactor;}
		
		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){std_detected_RV--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){severe_detected_RV--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case IV:
		testIV++;
		
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;

		detect  = (vaccSpec[this_one] > sensitivity[SICCT]) ? vaccSpec[this_one] : sensitivity[SICCT];
		detectS = (vaccSpecS[this_one] > sensitivity[SICCT_S]) ? vaccSpecS[this_one] : sensitivity[SICCT_S];
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		rnd_rain = gsl_rng_uniform(gsl_r);
		
		if(rnd_rain <= detect){std_detected_IV++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= detectS){severe_detected_IV++;(iter->second).Svr_status=Reactor;}
				
		//if(rnd_rain <= sensitivity[SICCT]){std_detected_IV++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_IV++;(iter->second).Svr_status=Reactor;}
		
				// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){std_detected_IV--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){severe_detected_IV--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case SS:
		testS++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(rnd_rain <= specificity[SICCT]){std_detected_S++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= specificity[SICCT_S]){severe_detected_S++;(iter->second).Svr_status=Reactor;}
		
		
		break;
		case SO:
		testO++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= (specificity[SICCT]))		{std_detected_O++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= (specificity[SICCT_S]))	{severe_detected_O++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SR:
		testR++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_R++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_R++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SI:
		testI++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_I++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_I++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SV1:
		testV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		if(rnd_rain <= vaccSpec[this_one]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}

		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_V1--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_V1--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		
		break;
		case SV2:
		testV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		//if(rnd_rain <= specificity[SICCT_V2]){std_detected_V2++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V2]){severe_detected_V2++;(iter->second).Svr_status=Reactor;}

				// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		if(rnd_rain <= vaccSpec[this_one]){std_detected_V2++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_V2++;(iter->second).Svr_status=Reactor;}


		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_V2--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_V2--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		break;
		case SOV1:
		testOV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= (specificity[SICCT_V1]))		{std_detected_OV1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= (specificity[SICCT_V1]))	{severe_detected_OV1++;(iter->second).Svr_status=Reactor;}
		
				// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		if(rnd_rain <= vaccSpec[this_one]){std_detected_OV1++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_OV1++;(iter->second).Svr_status=Reactor;}

		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_OV1--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_OV1--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}

		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		
		break;
		case SOV2:
		testOV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= (specificity[SICCT_V2]))		{std_detected_OV2++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= (specificity[SICCT_V2]))	{severe_detected_OV2++;(iter->second).Svr_status=Reactor;}
		
		// Look up specificity based on time from vaccination
	
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;
		if(rnd_rain <= vaccSpec[this_one]){std_detected_OV2++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= vaccSpecS[this_one]){severe_detected_OV2++;(iter->second).Svr_status=Reactor;}

		
		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
	
			if(rnd_rain > specificity[DIVA]){std_detected_V2--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
		
			if(rnd_rain > specificity[DIVA]){severe_detected_V2--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}

		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		
		break;
		case SRV:
		testRV++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= sensitivity[SICCT]){std_detected_RV++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_RV++;(iter->second).Svr_status=Reactor;}
		
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;

		detect  = (vaccSpec[this_one] > sensitivity[SICCT]) ? vaccSpec[this_one] : sensitivity[SICCT];
		detectS = (vaccSpecS[this_one] > sensitivity[SICCT_S]) ? vaccSpecS[this_one] : sensitivity[SICCT_S];
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		if(rnd_rain <= detect){std_detected_RV++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= detectS){severe_detected_RV++;(iter->second).Svr_status=Reactor;}
		
		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){std_detected_RV--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){severe_detected_RV--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}

		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		case SIV:
		testIV++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		//if(rnd_rain <= sensitivity[SICCT]){std_detected_IV++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_IV++;(iter->second).Svr_status=Reactor;}
		
		t_from_v = (t - (iter->second).Vaccinated_time);
		// If longer than end of array, use last value
		this_one=dim_vaccSpec-1;
		// Select appropriate specificity
		for(int z=0; z < dim_vaccSpec ;z++)
		{
		 if(t_from_v < vaccStime[z]){this_one=z;break;}  
		}
		
		//if(rnd_rain <= specificity[SICCT_V1]){std_detected_V1++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= specificity[SICCT_V1]){severe_detected_V1++;(iter->second).Svr_status=Reactor;}
		//cout << "Vaccinate Specificity: " << this_one << ' ' << t_from_v << ' ' << vaccSpec[this_one] << endl;

		detect  = (vaccSpec[this_one] > sensitivity[SICCT]) ? vaccSpec[this_one] : sensitivity[SICCT];
		detectS = (vaccSpecS[this_one] > sensitivity[SICCT_S]) ? vaccSpecS[this_one] : sensitivity[SICCT_S];
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		if(rnd_rain <= detect){std_detected_IV++;(iter->second).Std_status=Reactor;}
		if(rnd_rain <= detectS){severe_detected_IV++;(iter->second).Svr_status=Reactor;}
		
		// If we get reacting vaccinate, perform DIVA test to try and negate
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(DIVAnegate & (iter->second).Std_status==Reactor)
		{	
			std_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){std_detected_IV--;(iter->second).Std_status=NotReactor;std_DIVAnegate++;}
		}
		if(DIVAnegate & (iter->second).Svr_status==Reactor)
		{	
			severe_DIVAtests++;	
			if(rnd_rain > sensitivity[DIVA]){severe_detected_IV--;(iter->second).Svr_status=NotReactor;severe_DIVAnegate++;}
		}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		default:
  						cout << "New York, New York: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
		}
		
		}
		}
		if(tested != animals_to_test){cout << "Big Star: " << tested << ' ' << animals_to_test << endl;exit(1);}
		
				
		results[0] = (int) t;
		results[1] = std_detected_S + std_detected_O + std_detected_R + std_detected_I + std_detected_V1 + std_detected_V2 + std_detected_OV1 + std_detected_OV2 + std_detected_RV + std_detected_IV;
		results[2] = severe_detected_S + severe_detected_O + severe_detected_R + severe_detected_I + severe_detected_V1 + severe_detected_V2 + severe_detected_OV1 + severe_detected_OV2 + severe_detected_RV + severe_detected_IV;
		results[3] = Rtot[this_herd] + Itot[this_herd] 
					 + SRtot[this_herd] + SItot[this_herd] 
					 + RVtot[this_herd] + IVtot[this_herd];
					 + SRVtot[this_herd] + SIVtot[this_herd]; 
		results[4] = Otot[this_herd] + Rtot[this_herd] + Itot[this_herd] 
					+ SOtot[this_herd] + SRtot[this_herd] + SItot[this_herd] 
					+ OV1tot[this_herd] + OV2tot[this_herd] + RVtot[this_herd] + IVtot[this_herd];
					+ SOV1tot[this_herd] + SOV2tot[this_herd] + SRVtot[this_herd] + SIVtot[this_herd]; 
		results[5] = animals_to_test;
		
		// Fixed probability of confirmation: only allow visible lesions to be detected in true infected reactors
		// At standard interpretation (used to flip to severe interpretation only)
		results[6] = vl_stddetected;
		
		results[7] = std_detected_S;
		results[8] = std_detected_O;
		results[9] = std_detected_R;
		results[10] = std_detected_I;
		results[11] = severe_detected_S;
		results[12] = severe_detected_O;
		results[13] = severe_detected_R;
		results[14] = severe_detected_I;
		results[15] = std_detected_V1;
		results[16] = std_detected_V2;
		results[17] = std_detected_OV1;
		results[18] = std_detected_OV2;
		results[19] = std_detected_RV;
		results[20] = std_detected_IV;
		results[21] = severe_detected_V1;
		results[22] = severe_detected_V2;
		results[23] = severe_detected_OV1;
		results[24] = severe_detected_OV2;
		results[25] = severe_detected_RV;
		results[26] = severe_detected_IV;
		results[27] = testS;
		results[28] = testO;
		results[29] = testR;
		results[30] = testI;
		results[31] = testV1;
		results[32] = testV2;
		results[33] = testOV1;
		results[34] = testOV2;
		results[35] = testIV;
		results[36] = testRV;
		results[37] = std_DIVAtests;
		results[38] = severe_DIVAtests;
		results[39] = std_DIVAnegate;
		results[40] = severe_DIVAnegate;
		
		delete [] elegible;	
		
		//audit_maps("Tested");
		
		return(results);
		
	}
	
	vector<int> bTBICBM::DIVA_Protocol_Per_Animal(int this_herd,int PTI,test_t whole_herd_test,double eligible_vacc=0.0)
	{
	
		// Count U_cases as first test positive (std) test result
		int U_Cases = 0;
		// Count V_cases as first test positive DIVA test result
		int V_Cases = 0;
		// Count R_cases as first test positive (std) test result
		int R_Cases = 0;
		// Count v_cases as first test positive DIVA test result
		int W_Cases = 0;
		
		
		//cout << "DIVA! " << DIVA << endl;

		//calculate_totals();

		current_testtype = whole_herd_test;
		
		int animals_to_test = 0;
		bool* elegible      = new bool[cows[this_herd].size()];
		int i=0;
		if(whole_herd_test==all)
		{
		
			animals_to_test = (cows[this_herd].size());
				
				map<double,cow_t>::iterator iter; 
  				for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
				{
				 switch((iter->second).Epi_status)
  					{
  						case S:
  								elegible[i] = true;
								i++;	
  						break;
  						case O:
								elegible[i] = true;
								i++;	
  						break;
  						case R:
								elegible[i] = true;
								i++;	
  						break;
  						case I:
  								elegible[i] = true;
								i++;	
  						break;
  						case V1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								i++;
								animals_to_test--;
								cout << "Nope" << endl;
								}	
  						break;
  						case V2:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case OV1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case OV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case RV:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case IV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}		
  						break;
  						case SS:
  								elegible[i] = true;
								i++;	
  						break;
  						case SO:
  								elegible[i] = true;
								i++;	
  						break;
  						case SR:
  								elegible[i] = true;
								i++;	
  						break;
  						case SI:
  								elegible[i] = true;
								i++;	
  						break;
  						case SV1:
 								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case SV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case SOV1:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  		
  						break;
  						case SOV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case SRV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						case SIV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = true;
								i++;
								}
								else
								{
								elegible[i] = false;
								animals_to_test--;
								i++;
								cout << "Nope" << endl;	
								}	
  						break;
  						default:
  						cout << "I did it my way: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
  						}
				//cout << "Eligibility: " << (iter->second).Epi_status << ' ' << elegible[i-1] << endl; 
				
				}
		}
		else
		{
			//cout << "This should never print" << endl;
			// Choose Routine tests to be WHT/RHT acording to proportions that initiate breakdowns in 
			// Study Population 
			// PTI 1: 52 RHT, 3235 WHT, SLH 781 p(WHT) = 3235/(3235+52) = 0.9841801
		    // PTI 2: 443 RHT, 478 WHT, SLH 282 p(WHT) = 443/(478+443) = 0.48
			// PTI 4: 555 RHT, 111 WHT, SLH 268 p(WHT) = 111/(555+111) = 0.1666667
			
			switch(PTI)
			{		
					// PTI 1
					case 1:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.9841801){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
					// PTI 2
					case 2:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.48){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
					//case 2:
					//if(gsl_ran_flat(gsl_r,0,1) <= 0.0){whole_herd_test = wht;}
					//else{whole_herd_test = rht;}
					break;
					// PTI 4
					case 4:
					if(gsl_ran_flat(gsl_r,0,1) <= 0.1666667){whole_herd_test = wht;}
					else{whole_herd_test = rht;}
					break;
			}
			
			
			map<double,cow_t>::iterator iter; 
			for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
			{
					
			switch(whole_herd_test)
			{
				case wht:
					
					if((t-(iter->second).Birth_time) > (6*7)){elegible[i] = true;animals_to_test++;}
     						else{ elegible[i] = false;}
				break;
				case rht:
					
     				if((t-(iter->second).Birth_time) > (2*364))
     						{elegible[i] = true;animals_to_test++;}
     					else{ elegible[i] = false;}
     			break;
			}
			// Remove eligibility for relevant vaccinates
			
			switch((iter->second).Epi_status)
  			{
 
  						case V1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								}
								else
								{
								elegible[i] = elegible[i] && false;	
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case V2:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case OV1:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case OV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case RV:
								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						case IV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}		
  						break;
  						
  						case SV1:
 								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
							
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						case SOV1:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  		
  						break;
  						case SOV2:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SRV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
							
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;
								}	
  						break;
  						case SIV:
  								if((t - (iter->second).Vaccinated_time) >= eligible_vacc)
								{
								elegible[i] = elegible[i] && true;
								
								}
								else
								{
								elegible[i] = elegible[i] && false;
								animals_to_test--;cout << "NOPE" << endl;	
								}	
  						break;
  						}
			
			i++;
			
			}
					
					
			}
			
	//cout << "End Eligibility Test" << endl;	
	/*
	cout << "Testing " 
	<< animals_to_test/(double) cows[this_herd].size() 
	<< "% : " 
	<< animals_to_test 
	<< " out of " << cows[this_herd].size() << endl;
	*/
	
    //animals_to_test = (cows[this_herd].size());
     		
		int testS=0;
		int testO=0;
		int testOV1=0;
		int testOV2=0;
		int testV1=0;
		int testV2=0;
		int testR=0;
		int testI=0;
		int testRV=0;
		int testIV=0;
				
		int std_detected_I = 0;	    int severe_detected_I = 0;
		int std_detected_R = 0;	    int severe_detected_R = 0;
		int std_detected_O = 0;		int severe_detected_O = 0;
		int std_detected_OV1 = 0;	int severe_detected_OV1 = 0;
		int std_detected_V1 = 0;	int severe_detected_V1 = 0;
		int std_detected_OV2 = 0;	int severe_detected_OV2 = 0;
		int std_detected_V2 = 0;	int severe_detected_V2 = 0;
		int std_detected_IV = 0;	int severe_detected_IV = 0;
		int std_detected_RV = 0;	int severe_detected_RV = 0;
		
		int std_detected_S = 0;	    int severe_detected_S = 0;
		
		int vl_stddetected=0;
		
		int std_DIVAnegate=0;	int severe_DIVAnegate=0;
		int std_DIVAtests=0;	int severe_DIVAtests=0;

		i = 0;int tested=0;
		
		double rnd_rain;
		
		map<double,cow_t>::iterator iter; 
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		//cout << i << ' ' << test_ind << ' ' << daPick[test_ind] << endl;

  
		 if(!elegible[i])
		 {
		 i++;
		 }
		 else
		 {
		  tested++;i++;
		  
		double detect  = 0.0;
		double detectS = 0.0;
		double t_from_v = 0.0;
		int this_one = 0;  
		  
		  
		switch((iter->second).Epi_status)
		{
		case S:
		testS++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		if(rnd_rain <= specificity[SICCT]){std_detected_S++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= specificity[SICCT_S]){severe_detected_S++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		break;
		case O:
		testO++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= (specificity[SICCT]))		{std_detected_O++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= (specificity[SICCT_S]))	{severe_detected_O++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		case R:
		testR++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_R++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_R++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case I:
		testI++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_I++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_I++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case V1:

		testV1++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		// Use DIVA test for all vaccinates
		
		if(rnd_rain <= specificity[DIVA]){std_detected_V1++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_V1++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}

		break;
		case V2:
		testV2++;

		rnd_rain = gsl_rng_uniform(gsl_r);
		
		// Use DIVA test for all vaccinates
				
		if(rnd_rain <= specificity[DIVA]){std_detected_V2++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_V2++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}

		break;		
		case OV1:
		testOV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);

	    // Use DIVA test for all vaccinates
					
		if(rnd_rain <= specificity[DIVA]){std_detected_OV1++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_OV1++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case OV2:
		testOV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		// Use DIVA test for all vaccinates
		
		if(rnd_rain <= specificity[DIVA]){std_detected_OV2++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_OV2++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		

		break;
		case RV:
		testRV++;
		
		detect  = sensitivity[DIVA];
		detectS = sensitivity[DIVA];
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;		
		
		if(rnd_rain <= detect){std_detected_RV++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= detectS){severe_detected_RV++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
	
		//if(rnd_rain <= sensitivity[SICCT]){std_detected_RV++;(iter->second).Std_status=Reactor;}
		//if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_RV++;(iter->second).Svr_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case IV:
		testIV++;
		
		//detect  = (sensitivity[DIVA]) ? vaccSpec[this_one] : sensitivity[SICCT];
		//detectS = (sensitivity[DIVA]) ? vaccSpecS[this_one] : sensitivity[SICCT_S];
		
		detect  = (sensitivity[DIVA]);
		detectS = (sensitivity[DIVA]);
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		rnd_rain = gsl_rng_uniform(gsl_r);
		
		if(rnd_rain <= detect){std_detected_IV++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= detectS){severe_detected_IV++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
						
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		break;
		case SS:
		testS++;
		rnd_rain = gsl_rng_uniform(gsl_r);
		if(rnd_rain <= specificity[SICCT]){std_detected_S++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){U_Cases++;(iter->second).Diva_ever=true;}}
		if(rnd_rain <= specificity[SICCT_S]){severe_detected_S++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		
		break;
		case SO:
		testO++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= (specificity[SICCT]))		{std_detected_O++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= (specificity[SICCT_S]))	{severe_detected_O++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SR:
		testR++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_R++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_R++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SI:
		testI++;
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= sensitivity[SICCT]){std_detected_I++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){U_Cases++;}else{R_Cases++;}}}
		if(rnd_rain <= sensitivity[SICCT_S]){severe_detected_I++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		break;
		case SV1:
		testV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		
		if(rnd_rain <= specificity[DIVA]){std_detected_V1++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_V1++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		break;
		case SV2:
		testV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
		
		if(rnd_rain <= specificity[DIVA]){std_detected_V2++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_V2++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}

		break;
		case SOV1:
		testOV1++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);

		if(rnd_rain <= specificity[DIVA]){std_detected_OV1++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_OV1++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}

		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		
		break;
		case SOV2:
		testOV2++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		if(rnd_rain <= specificity[DIVA]){std_detected_OV2++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= specificity[DIVA]){severe_detected_OV2++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}

		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}

		
		break;
		case SRV:
		testRV++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
	
		detect  = sensitivity[DIVA];
		detectS = sensitivity[DIVA];
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		if(rnd_rain <= detect){std_detected_RV++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= detectS){severe_detected_RV++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
		
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		case SIV:
		testIV++;
		
		rnd_rain = gsl_rng_uniform(gsl_r);
			
		detect  = sensitivity[DIVA];
		detectS = sensitivity[DIVA];
		
		//cout << "Detect: " << detect << ' ' << detectS << endl;

		if(rnd_rain <= detect){std_detected_IV++;(iter->second).Std_status=Reactor;(iter->second).Diva_status=Reactor;if(!(iter->second).Diva_ever){(iter->second).Diva_ever=true;if(!(iter->second).Seeder){V_Cases++;}else{W_Cases++;}}}
		if(rnd_rain <= detectS){severe_detected_IV++;(iter->second).Svr_status=Reactor;(iter->second).Diva_status=Reactor;}
				
		if(IsLesioned(iter))
		{	(iter->second).Confirmation_status = true;
			if((iter->second).Std_status==Reactor){vl_stddetected++;}
		}
		
		break;
		default:
  						cout << "New York, New York: " << endl;
  						cout << (iter->second).Epi_status << endl;
  						cout << (iter->second).Birth_time << endl;
    					cout << (iter->second).Death_time << endl;
    					cout << (iter->second).Off_time << endl;
   						cout << (iter->second).Infection_time << endl;
    					cout << (iter->second).Vaccinated_time << endl;
						cout << (iter->second).Epi_stage << endl;
						cout << (iter->second).Std_status << endl;
						cout << (iter->second).Svr_status << endl;
						cout << (iter->second).Diva_status << endl;
						cout << (iter->second).Confirmation_status << endl;
						cout << (iter->second).rate << endl;
  						exit(1);
  						break; 
		}
		
		}
		}
		if(tested != animals_to_test){cout << "Big Star: " << tested << ' ' << animals_to_test << endl;exit(1);}
		
		//cout << "End Testing " << endl;
		
		vector<int> results(45,0);
				
		results[0] = (int) t;
		results[1] = std_detected_S + std_detected_O + std_detected_R + std_detected_I + std_detected_V1 + std_detected_V2 + std_detected_OV1 + std_detected_OV2 + std_detected_RV + std_detected_IV;
		results[2] = severe_detected_S + severe_detected_O + severe_detected_R + severe_detected_I + severe_detected_V1 + severe_detected_V2 + severe_detected_OV1 + severe_detected_OV2 + severe_detected_RV + severe_detected_IV;
		results[3] = Rtot[this_herd] + Itot[this_herd] 
					 + SRtot[this_herd] + SItot[this_herd] 
					 + RVtot[this_herd] + IVtot[this_herd]
					 + SRVtot[this_herd] + SIVtot[this_herd]; 
		results[4] = Otot[this_herd] + Rtot[this_herd] + Itot[this_herd] 
					+ SOtot[this_herd] + SRtot[this_herd] + SItot[this_herd] 
					+ OV1tot[this_herd] + OV2tot[this_herd] + RVtot[this_herd] + IVtot[this_herd]
					+ SOV1tot[this_herd] + SOV2tot[this_herd] + SRVtot[this_herd] + SIVtot[this_herd]; 
		results[5] = animals_to_test;

	    // Debug Code
	    // Trap for negative total counts
	    
	    /*

		if(Rtot[this_herd] < 0 || Itot[this_herd] < 0 || Rtot[this_herd] < 0 || SRtot[this_herd] < 0 ||
		SItot[this_herd] < 0 || RVtot[this_herd] < 0 || IVtot[this_herd] < 0 || SRVtot[this_herd] < 0 || SIVtot[this_herd])
		{

		cout << "R" << ' ' << "I" << ' ' << "SR" << ' ' << "SI" << ' ' << "RV" << ' ' << "IV" << ' ' << "SRV" << ' ' << "SIV" << endl;

		cout << Rtot[this_herd] << ' ' <<  Itot[this_herd] << ' ' <<  SRtot[this_herd]  << ' ' << SItot[this_herd] 
			<< ' ' <<  RVtot[this_herd] << ' ' <<  IVtot[this_herd] << ' ' <<  SRVtot[this_herd] << ' ' <<  SIVtot[this_herd] << endl; 
		}	
		
		if(Otot[this_herd] < 0 || OV1tot[this_herd] < 0 || OV2tot[this_herd] < 0 || SOV1tot[this_herd] < 0 ||
		SOV2tot[this_herd] < 0)
		{

		cout << "O" << ' ' << "OV1" << ' ' << "OV2" << ' ' << "SOV1" << ' ' << "SOV2" << endl;

		cout << Otot[this_herd] << ' ' <<  OV1tot[this_herd] << ' ' <<  OV2tot[this_herd]  << ' ' << SOV1tot[this_herd] 
			<< ' ' <<  SOV2tot[this_herd] << endl; 
		}	
		
		*/

		// Fixed probability of confirmation: only allow visible lesions to be detected in true infected reactors
		// At standard interpretation (used to flip to severe interpretation only)
		results[6] = vl_stddetected;
		
		results[7] = std_detected_S;
		results[8] = std_detected_O;
		results[9] = std_detected_R;
		results[10] = std_detected_I;
		results[11] = severe_detected_S;
		results[12] = severe_detected_O;
		results[13] = severe_detected_R;
		results[14] = severe_detected_I;
		results[15] = std_detected_V1;
		results[16] = std_detected_V2;
		results[17] = std_detected_OV1;
		results[18] = std_detected_OV2;
		results[19] = std_detected_RV;
		results[20] = std_detected_IV;
		results[21] = severe_detected_V1;
		results[22] = severe_detected_V2;
		results[23] = severe_detected_OV1;
		results[24] = severe_detected_OV2;
		results[25] = severe_detected_RV;
		results[26] = severe_detected_IV;
		results[27] = testS;
		results[28] = testO;
		results[29] = testR;
		results[30] = testI;
		results[31] = testV1;
		results[32] = testV2;
		results[33] = testOV1;
		results[34] = testOV2;
		results[35] = testIV;
		results[36] = testRV;
		results[37] = std_DIVAtests;
		results[38] = severe_DIVAtests;
		results[39] = std_DIVAnegate;
		results[40] = severe_DIVAnegate;
		results[41] = U_Cases;
		results[42] = R_Cases;
		results[43] = V_Cases;
		results[44] = W_Cases;
		delete [] elegible;	
		
		//audit_maps("DIIVA'd");
		
	//	cout << "End Results Composition" << endl;
		
		return(results);
		
	}	
	

	int  bTBICBM::number_testWHT(int i)
	{
		// Test all animals except those under the age of 6 weeks (i.e. practically all!)
	return(gsl_ran_binomial(gsl_r,(1-(6.0/52.0)*364.0*Life_Births[i]/(cows[i].size())),cows[i].size()));
	}
	
	// Sample from empirical distribution RHT
	
	int  bTBICBM::number_testRHT(int i)
	{
		
		//double pval = RHT[(int) round(gsl_ran_flat(gsl_r,0,no_of_RHTs))];
		
		//return(gsl_ran_binomial(gsl_r,pval,N[i]));
		
		double pval=0.0;
		do
		{
			pval = 0.4941136064 + gsl_ran_cauchy(gsl_r,0.0932207681);
		}while(pval < 0 || pval > 1.0);
		
		return(gsl_ran_binomial(gsl_r,pval,cows[i].size()));
		
	}
	
	void bTBICBM::Cull(int this_herd, vector<int> results,bool confirmed, bool disclose_tests=false, int PTI=1,bool retain=false)
	{
	
		int LU = 0;
		int LV = 0;
		
		map<double,cow_t> remove_these;
		double key_index;
		epi_status_t shadow;
		//cout << "Parkway\n";
		if(results.size() < 36)
		{
			cout << "It's life, but not as we know it." << endl;
		}

		//cout << "Herd Size before testing: " << cows[this_herd].size() << endl;
		
		int Removed = 0;
		
		double L = 0.0;
		
		if(!retain)
		{
		
		if(confirmed)
		{
		map<double,cow_t>::iterator iter;   
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		 if((iter->second).Std_status == Reactor || (iter->second).Svr_status == Reactor)
		 {
			remove_these.insert(std::pair<double, cow_t>(iter->first, iter->second));
		 }
		 else // Reset lesions flag for non-reactors
		 {
		 (iter->second).Confirmation_status = false;
		 }
		}
		}
		else
		{
		map<double,cow_t>::iterator iter;
		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		
		if((iter->second).Std_status == Reactor)
		{
			remove_these.insert(std::pair<double, cow_t>(iter->first, iter->second));
		}
		else if((iter->second).Svr_status == Reactor)
		{  			
		  // Reset all severe reactors back to NR for next test
		  (iter->second).Svr_status=NotReactor;
		  (iter->second).Confirmation_status = false;
		  //cout << "RESET SEVERE " << t << endl;
		}
		else  // Reset lesion flag for non-reactors
		{
		(iter->second).Confirmation_status = false;
		}
		
		
		}
		}
		
	
		map<double,cow_t>::iterator iter;
  		for( iter = remove_these.begin(); iter != remove_these.end(); iter++ ) 
		{
		
		  Removed++;
		  double * demo_dat;
		  int recorditas = 0;
	
		  switch(Seasonal_model)
		  {
			case constant:
			/*
			for(int j=0;j < cows[this_herd].size(); j++)
			{
			 cout << cows[this_herd][j].Death_time << ' ';
			 }
			 cout << endl;
			*/
			
			
			
			//Reactorfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
			if(full_save)
			{
			Reactorfile[this_herd] << t << ' ' << (t - (iter->second).Birth_time) << endl;
				if(debug_save)
				{
				SIndividualfile[this_herd] 
				<< DarthSelecta << ','
				<< 1-(iter->second).Control << ','
				<< (iter->second).Birth_time << ','
				<< (iter->second).Death_time << ',' 
				<< (iter->second).Vaccinated_time << ',' 
				<< 1 << ','
				<< (iter->second).Confirmation_status << endl;
				}
			}
			
		  	if((t - (iter->second).Birth_time) > maxRangeAge)
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			//cout << "Record!" << ' ' << gsl_histogram_sum(AgeReactors) << endl;
			gsl_histogram_increment(AgeReactors,recorditas);
			//cout << gsl_histogram_sum(AgeReactors) << endl;

			if((iter->second).Confirmation_status)
			{
			
			
				if(full_save)
			{
			Confirmedfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
			}
			
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			
			gsl_histogram_increment(AgeCReactors,recorditas);
			// Bookkeep vaccinated and unvaccinated with lesions
			if((iter->second).Epi_status == S || (iter->second).Epi_status == SS ||
			 (iter->second).Epi_status == O  || (iter->second).Epi_status == SO ||
			 (iter->second).Epi_status == R  || (iter->second).Epi_status == SR ||
			(iter->second).Epi_status == I  || (iter->second).Epi_status == SI )
			{LU++;}
			else{LV++;}
			
			
			}
			
			switch(iter->second.Epi_status)
			{
			case S:
			Stot[this_herd]--;
			break;
			case O:
			Otot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case R:
			Rtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case I:
			Itot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case V1:
			V1tot[this_herd]--;
			break;
			case V2:
			V2tot[this_herd]--;
			break;
			case OV1:
			OV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case OV2:
			OV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case RV:
			RVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case IV:
			IVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			
			case SS:
			SStot[this_herd]--;
			break;
			case SO:
			
			SOtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
		
			break;
			case SR:
			
			SRtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}

			break;
			case SI:
			SItot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SV1:
			SV1tot[this_herd]--;
			break;
			case SV2:
			SV2tot[this_herd]--;
			break;
			case SOV1:
			SOV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SOV2:
			SOV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SRV:
			SRVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SIV:
			SIVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			}
			
			cows[this_herd].erase(iter->first);
		
			do{
			demo_dat = Get_New_Cow(NextMove,true);
			cow_t new_cow;
	
			new_cow.Death_time = demo_dat[Death_demo];
			new_cow.Off_time = demo_dat[Off_demo];
			new_cow.Birth_time = demo_dat[Birth_demo];
	
			new_cow.uniq_id = gsl_ran_flat(gsl_r,0,3e8);
	
			delete [] demo_dat;
	
			// If birth, animal is susceptible.
			// else chance (pinf_move) animal is infected (S=0,O=1)
			new_cow.Epi_status = ((t-new_cow.Birth_time) > 0.0) ? S : (epi_status_t) gsl_ran_binomial (gsl_r, pinf_move, 1);
			if(new_cow.Epi_status==S)
			{
			 new_cow.Infection_time = 0;
			}
			 else
			{
			 new_cow.Infection_time = t; // Set infection time to time on herd
			}
			
			//new_cow.Epi_stage = 0;
			new_cow.Std_status = NotReactor;	
			new_cow.Svr_status = NotReactor;
			new_cow.Diva_status = NotReactor;
				new_cow.Diva_ever = false;
			new_cow.Confirmation_status = false;
			new_cow.Vaccinated_time = -1;
    		new_cow.Control = gsl_ran_binomial (gsl_r,  1.0-target_vaccination_p, 1);
    		new_cow.Seeder = false;
    		
			// Is cow shadow-progressor?

			shadow = (epi_status_t) gsl_ran_binomial (gsl_r, p_shadow, 1);
			new_cow.Epi_status = (epi_status_t) ((int) new_cow.Epi_status + 10*shadow);
	
			key_index = (new_cow.Death_time < new_cow.Off_time) ? new_cow.Death_time : new_cow.Off_time;
	
			    // check that key does not exist in herd, if it does then resample
	if(cows[this_herd].find(key_index)==cows[this_herd].end())
	{
	cows[this_herd].insert(std::pair<double, cow_t>(key_index,new_cow));
	if(new_cow.Epi_status == O || new_cow.Epi_status == SO)
	{
	infected[this_herd].insert(std::pair<double, cow_t>(new_cow.uniq_id,new_cow));
	}
	break;
	}
	else{continue;}
	}while(true);
	
			break;
			case seasonal:
			if(full_save)
			{
			//cout << "This should not print" << endl;
			//Reactorfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
			Reactorfile[this_herd] << t << ' ' << (t - (iter->second).Birth_time) << endl;
				if(debug_save)
				{
				SIndividualfile[this_herd] 
				<< DarthSelecta << ','
				<< 1-(iter->second).Control << ','
				<< (iter->second).Birth_time << ','
				<< (iter->second).Death_time << ',' 
				<< (iter->second).Vaccinated_time << ',' 
				<< 1 << ','
				<< (iter->second).Confirmation_status << endl;
				}
			}
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			//cout << "Record!" << ' ' << gsl_histogram_sum(AgeReactors) << endl;

			gsl_histogram_increment(AgeReactors,recorditas);
			
			//cout << "Record!" << ' ' << gsl_histogram_sum(AgeReactors) << endl;


			if((iter->second).Confirmation_status)
			{
				if(full_save)
	{
			Confirmedfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
	}		
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			
			
			// Bookkeep vaccinated and unvaccinated with lesions
			if((iter->second).Epi_status == S || (iter->second).Epi_status == SS ||
			 (iter->second).Epi_status == O  || (iter->second).Epi_status == SO ||
			 (iter->second).Epi_status == R  || (iter->second).Epi_status == SR ||
			(iter->second).Epi_status == I  || (iter->second).Epi_status == SI )
			{LU++;}
			else{LV++;}
			
			gsl_histogram_increment(AgeCReactors,recorditas);
			}
			
			
			
			switch(iter->second.Epi_status)
			{
			case S:
			Stot[this_herd]--;
			break;
			case O:
			Otot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}

			break;
			case R:
			Rtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case I:
			Itot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case V1:
			V1tot[this_herd]--;
			
			break;
			case V2:
			V2tot[this_herd]--;
			
			break;
			case OV1:
			OV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case OV2:
			OV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case RV:
			RVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case IV:
			IVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SS:
			SStot[this_herd]--;

			break;
			case SO:
			SOtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SR:
			SRtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SI:
			SItot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SV1:
			SV1tot[this_herd]--;
			break;
			case SV2:
			SV2tot[this_herd]--;
			break;
			case SOV1:
			SOV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SOV2:
			SOV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SRV:
			SRVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SIV:
			SIVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			}
			
			cows[this_herd].erase(iter->first);
			//cout << "CULL COW " << t << " " << cows[this_herd].size() << endl;
			break;
			}	

		
		  
		 }
		 // Clean up
		 if(remove_these.size()>0)
		{remove_these.erase(remove_these.begin(),remove_these.end());}
		 
		 update_next_move();	
		 
		 }
				
		// Record Testing History
		if(disclose_tests && full_save)
		{
			for(int i=0; i < 41; i++)
			{
				disclose[this_herd] << results[i] << ' ';
			}
			disclose[this_herd] << 0 << ' ' << Removed << ' ' << PTI << ' ' << DarthSelecta << ' ' << current_testtype << ' ' << cows[this_herd].size() << ' ';
			
			// Output RU, RV, LU, LV, SU, IU, SV, IV
			disclose[this_herd] << (results[7]+results[8]+results[9]+results[10]) << ' '
								<< (results[15]+results[16]+results[17]+results[18]+results[19]+results[20]) << ' '
			 					<< LU << ' '
			 					<< LV << ' '
								 << (Stot[this_herd] + SStot[this_herd]) << ' '
			 			 
				 << (Otot[this_herd] + Rtot[this_herd] + Itot[this_herd] 
				+ SOtot[this_herd] + SRtot[this_herd] + SItot[this_herd]) << ' '
				
			 << (V1tot[this_herd] + V2tot[this_herd] + SV1tot[this_herd] + SV2tot[this_herd]) << ' '
				
			 << (OV1tot[this_herd] + OV2tot[this_herd] + RVtot[this_herd] + IVtot[this_herd]
					+ SOV1tot[this_herd] + SOV2tot[this_herd] + SRVtot[this_herd] + SIVtot[this_herd]) << endl;
					
					
			 
			
		//cows[this_herd].size()
		}
	
	
				  
		
	}
	
	void bTBICBM::Whole_Herd_Cull(int this_herd, vector<int> results,bool confirmed, bool disclose_tests=false, int PTI=1,bool retain=false)
	{
	
		int LU = 0;
		int LV = 0;
		
		map<double,cow_t> remove_these;
		double key_index;
		epi_status_t shadow;
		//cout << "Parkway\n";
		if(results.size() < 36)
		{
			cout << "It's life, but not as we know it." << endl;
		}

		//cout << "Herd Size before testing: " << cows[this_herd].size() << endl;
		
		int Removed = 0;
		
		double L = 0.0;
		
		if(!retain)
		{
		
		if(confirmed)
		{
		map<double,cow_t>::iterator iter;   
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		 if((iter->second).Std_status == Reactor || (iter->second).Svr_status == Reactor)
		 {
			remove_these.insert(std::pair<double, cow_t>(iter->first, iter->second));
		 }
		 else // Reset lesions flag for non-reactors
		 {
		 (iter->second).Confirmation_status = false;
		 }
		}
		}
		else
		{
		map<double,cow_t>::iterator iter;
		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		
		if((iter->second).Std_status == Reactor)
		{
			remove_these.insert(std::pair<double, cow_t>(iter->first, iter->second));
		}
		else if((iter->second).Svr_status == Reactor)
		{  			
		  // Reset all severe reactors back to NR for next test
		  (iter->second).Svr_status=NotReactor;
		  (iter->second).Confirmation_status = false;
		  //cout << "RESET SEVERE " << t << endl;
		}
		else  // Reset lesion flag for non-reactors
		{
		(iter->second).Confirmation_status = false;
		}
		
		
		}
		}
		
	
		map<double,cow_t>::iterator iter;
  		for( iter = remove_these.begin(); iter != remove_these.end(); iter++ ) 
		{
		
		  Removed++;
		  
		  int recorditas = 0;
	
		  switch(Seasonal_model)
		  {
			case constant:
	
			if(full_save)
			{		
		
			//Reactorfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
			Reactorfile[this_herd] << t << ' ' << (t - (iter->second).Birth_time) << endl;
			
				if(debug_save)
				{
				SIndividualfile[this_herd] 
				<< DarthSelecta << ','
				<< 1-(iter->second).Control << ','
				<< (iter->second).Birth_time << ','
				<< (iter->second).Death_time << ',' 
				<< (iter->second).Vaccinated_time << ',' 
				<< 1 << ','
				<< (iter->second).Confirmation_status << endl;
				}
	}		
			
		  	if((t - (iter->second).Birth_time) > maxRangeAge)
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			gsl_histogram_increment(AgeReactors,recorditas);

			if((iter->second).Confirmation_status)
			{
			
			
				if(full_save)
	{
			Confirmedfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
	}		
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			
			gsl_histogram_increment(AgeCReactors,recorditas);
			// Bookkeep vaccinated and unvaccinated with lesions
			if((iter->second).Epi_status == S || (iter->second).Epi_status == SS ||
			 (iter->second).Epi_status == O  || (iter->second).Epi_status == SO ||
			 (iter->second).Epi_status == R  || (iter->second).Epi_status == SR ||
			(iter->second).Epi_status == I  || (iter->second).Epi_status == SI )
			{LU++;}
			else{LV++;}
			
			
			}
			
			switch(iter->second.Epi_status)
			{
			case S:
			Stot[this_herd]--;
			break;
			case O:
			Otot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}

			break;
			case R:
			Rtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case I:
			Itot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case V1:
			V1tot[this_herd]--;
			break;
			case V2:
			V2tot[this_herd]--;
			break;
			case OV1:
			OV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case OV2:
			OV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case RV:
			RVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case IV:
			IVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			
			case SS:
			SStot[this_herd]--;
			break;
			case SO:
			
			SOtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
		
			break;
			case SR:
			
			SRtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}

			break;
			case SI:
			SItot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SV1:
			SV1tot[this_herd]--;
			break;
			case SV2:
			SV2tot[this_herd]--;
			break;
			case SOV1:
			SOV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SOV2:
			SOV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SRV:
			SRVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SIV:
			SIVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			}
			
			cows[this_herd].erase(iter->first);
		
			break;
			case seasonal:
			
					if(full_save)
	{
			
			Reactorfile[this_herd] << t << ' ' << (t - (iter->second).Birth_time) << endl;
				if(debug_save)
				{
				SIndividualfile[this_herd] 
				<< DarthSelecta << ','
				<< 1-(iter->second).Control << ','
				<< (iter->second).Birth_time << ','
				<< (iter->second).Death_time << ',' 
				<< (iter->second).Vaccinated_time << ',' 
				<< 1 << ','
				<< (iter->second).Confirmation_status << endl;
			}
	}		
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			//cout << "Record!" << ' ' << gsl_histogram_sum(AgeReactors) << endl;

			gsl_histogram_increment(AgeReactors,recorditas);
			
			//cout << "Record!" << ' ' << gsl_histogram_sum(AgeReactors) << endl;


			if((iter->second).Confirmation_status)
			{
				if(full_save)
	{
			Confirmedfile[this_herd] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) << endl;
	}		
			if((t - (iter->second).Birth_time) > maxRangeAge) 
   			 {
    		recorditas = maxRangeAge-1;
   			 }
   			 else
    		{
        	recorditas = (t - (iter->second).Birth_time);
    		}
			
			
			// Bookkeep vaccinated and unvaccinated with lesions
			if((iter->second).Epi_status == S || (iter->second).Epi_status == SS ||
			 (iter->second).Epi_status == O  || (iter->second).Epi_status == SO ||
			 (iter->second).Epi_status == R  || (iter->second).Epi_status == SR ||
			(iter->second).Epi_status == I  || (iter->second).Epi_status == SI )
			{LU++;}
			else{LV++;}
			
			gsl_histogram_increment(AgeCReactors,recorditas);
			}
			
			
			
			switch(iter->second.Epi_status)
			{
			case S:
			Stot[this_herd]--;
			break;
			case O:
			Otot[this_herd]--;
						if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
break;
			case R:
			Rtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case I:
			Itot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case V1:
			V1tot[this_herd]--;
			
			break;
			case V2:
			V2tot[this_herd]--;
			
			break;
			case OV1:
			OV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case OV2:
			OV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case RV:
			RVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case IV:
			IVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SS:
			SStot[this_herd]--;

			break;
			case SO:
			SOtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SR:
			SRtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SI:
			SItot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SV1:
			SV1tot[this_herd]--;
			
			break;
			case SV2:
			SV2tot[this_herd]--;
			
			break;
			case SOV1:
			SOV1tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SOV2:
			SOV2tot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SRV:
			SRVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			case SIV:
			SIVtot[this_herd]--;
			if(infected[this_herd].erase(iter->first) == 0){cout << "Failed to match!"; exit(1);}
			break;
			}
			
			cows[this_herd].erase(iter->first);
			//cout << "CULL COW " << t << " " << cows[this_herd].size() << endl;
			break;
			}	

		
		  
		 }
		 // Clean up
		 if(remove_these.size()>0)
		{remove_these.erase(remove_these.begin(),remove_these.end());}
		 
		 update_next_move();	
		 
		 }
				
		// Record Testing History
		if(disclose_tests && full_save)
		{
			for(int i=0; i < 41; i++)
			{
				disclose[this_herd] << results[i] << ' ';
			}
			disclose[this_herd] << 0 << ' ' << Removed << ' ' << PTI << ' ' << DarthSelecta << ' ' << current_testtype << ' ' << cows[this_herd].size() << ' ';
			
			// Output RU, RV, LU, LV, SU, IU, SV, IV
			disclose[this_herd] << (results[7]+results[8]+results[9]+results[10]) << ' '
								<< (results[15]+results[16]+results[17]+results[18]+results[19]+results[20]) << ' '
			 					<< LU << ' '
			 					<< LV << ' '
								 << (Stot[this_herd] + SStot[this_herd]) << ' '
			 			 
				 << (Otot[this_herd] + Rtot[this_herd] + Itot[this_herd] 
				+ SOtot[this_herd] + SRtot[this_herd] + SItot[this_herd]) << ' '
				
			 << (V1tot[this_herd] + V2tot[this_herd] + SV1tot[this_herd] + SV2tot[this_herd]) << ' '
				
			 << (OV1tot[this_herd] + OV2tot[this_herd] + RVtot[this_herd] + IVtot[this_herd]
					+ SOV1tot[this_herd] + SOV2tot[this_herd] + SRVtot[this_herd] + SIVtot[this_herd]) << endl;
					
					
			 
			
	
		}
		
		// Kill everything else
		double recorditas = 0.0;
		map<double,cow_t>::iterator iter; 
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		bool slaughter_flag = slaughter_house_test(iter->second);
			if(debug_save)
	{
		SIndividualfile[0] 
		<< DarthSelecta << ','
		<< 1-(iter->second).Control << ','
		<< (iter->second).Birth_time << ','
		<< (iter->second).Death_time << ',' 
		<< (iter->second).Vaccinated_time << ',' 
		<< 0 << ','
		<< slaughter_flag << endl;
	}
		
	//cout << "French Cricket" << endl;
	if(slaughter_flag)
	{
	if(full_save)
	{
	Slaughterfile[0] << DarthSelecta << ' ' << (t - (iter->second).Birth_time) 
	<< ' ' << (iter->second).Epi_status << endl;
	}

	if((t - (iter->second).Birth_time) > maxRangeAge) 
    {
    	recorditas = maxRangeAge-1;
    }
    else
    {
        recorditas = (t - (iter->second).Birth_time);
    }
	
	gsl_histogram_increment(AgeSReactors,recorditas);
	
	}
		

	}
		
	//audit_maps("Cull");			  
		
	}

void bTBICBM::ClearStatus(int this_herd)
	{

	int herd_index = this_herd;

		map<double,cow_t>::iterator iter;   
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
	
	(iter->second).Std_status = NotReactor;	
	(iter->second).Svr_status = NotReactor;
	(iter->second).Diva_status = NotReactor;
	//	(iter->second).Diva_ever DO NOT CLEAR DIVA EVER STATUS!
	(iter->second).Confirmation_status   = false;
		
		}
	}
	
	int bTBICBM::CountSeeders(int this_herd)
	{

	int herd_index = this_herd;
    int cnt=0;
    
	map<double,cow_t>::iterator iter;   
  	for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
	{
	
	if((iter->second).Seeder)
	{cnt++;}
		
	}
	return(cnt);
	}
	
	int bTBICBM::InvCountSeeders(int this_herd)
	{

	int herd_index = this_herd;
    int cnt=0;
    
	map<double,cow_t>::iterator iter;   
  	for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
	{
	
	if(!(iter->second).Diva_ever && (iter->second).Seeder)
	{cnt++;}
		
	}
	return(cnt);
	}
	
	
void bTBICBM::RemoveSeeders()
	{

	double BIGNUM = 364.0*1000;

	int herd_index = 0;
	int this_herd = 0;
	

	
	//clear infected list
	if(infected[herd_index].size()>0)
	{infected[herd_index].erase(infected[herd_index].begin(),infected[herd_index].end());}
	

	  // Retain all test positives with Control flag=keep from previous phase
		// Remove all other animals from group
		
		map<double,cow_t> remove_these;
	
		map<double,cow_t>::iterator iter;   
  		for( iter = cows[this_herd].begin(); iter != cows[this_herd].end(); iter++ ) 
		{
		//cout << "Iterate" << endl;
		
		 // If seeder animal or DIVA negative, flag to remove from group
		 if( ((iter->second).Seeder))
		 {
		
		 remove_these.insert(std::pair<double, cow_t>(iter->first, iter->second));
		 }
		else
		 {
		
		// Toggle flag to indicate animal is now a seeder animal
  		(iter->second).Seeder = true;
		
		 // Rebuild infectious list
		 if((iter->second).Epi_status != S & (iter->second).Epi_status != V1 & (iter->second).Epi_status != V2 & (iter->second).Epi_status != SS & (iter->second).Epi_status != SV1 & (iter->second).Epi_status != SV2)
		{
	

		// Add to infected list
		std::pair<std::map<double, cow_t>::iterator,bool> ret;
  		ret = infected[herd_index].insert(std::pair<double, cow_t>((iter->first),(iter->second)));

  		//cout << "Seeding" << endl;
  		if (ret.second==false) {
    	std::cout << "element 'z' already existed";
  		}

		}
  	
		}
		
    	// Remove animals from Group	
		
		}
		
		for( iter = remove_these.begin(); iter != remove_these.end(); iter++ ) 
		{
		cows[this_herd].erase(iter->first);
		}
	}	
	


	
	int bTBICBM::accum_infected_cow_days()
	{
		
		int rval = prob_on_trans;
		prob_on_trans = 0;
		
		return( rval);
		
	}
	
	// gamma function
	int gamma(int x);
	
	// Helper functions for 2,3 matrix indexing
	inline void deref3(int l, int y, int z, int & i, int & a, int & k);
	inline void deref2(int l, int y, int & i, int & a);
	inline int index3(int i, int a, int k, int y, int z);
	inline int index2(int i, int a, int y);
	
	// Helper function to add number <no> to filename <s>
	inline char *filename(char *s, int no);


int bTBICBM::reacto(double p, int n)
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
	
	
	double ranx = gsl_ran_flat(gsl_r,0,1);
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


	void  bTBICBM::run_uk_testing(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear = 1,bool clear_from_severe = false,bool change_SIT=false,int second_SIT=60, bool perfect_isolation=false)
	{
	
	TargetHerdSize = Herd_Size;
	
	int this_PTI = set_this_PTI;
	
	int clear_cnt = 0;
				
	vector<int> policy_check;			
	// Set current routine testing distribution
		switch(this_PTI)
		{
		// PTI 1
		case 1:
		Current_PTI = PTI1;
		Current_PTI_Num = no_of_PTI1;
		break;
		// PTI 2
		case 2:
		Current_PTI = PTI2;
		Current_PTI_Num = no_of_PTI2;
		break;
		// PTI 4
		case 4:
		Current_PTI = PTI4;
		Current_PTI_Num = no_of_PTI4;
		break;
		}
	
		//int test_period=this_PTI*364.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
				
		int test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
		
				
		int this_seed = 0;
				
		test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
				
		set_time();
		
		//cout << model.return_time() << endl;
				
		breakstart=0.0;
		breakfirst=-1;
		break_recurr=0.0;
		primary_breaklength = 0.0;
		primary_breakstart  = 0.0;
		Reactors_at_Start = 0;
		Reactors_at_VE6M  = 0;
		Reactors_at_VE12M = 0;
		Primary_Reactors  = 0;
		Reactors_at_First_Test = 0;
		forward_trans = 0;
		break6=false;
		break12=false;
		break24=false;
		confirmed=false;
		confirmed_ever=false;
		confirmed_at_start = false;
		slaughter_house = false;
		
		Burden			= 0;
		Total_Visits    = 0;
		Total_Tests		= 0;
		Total_Diva_tests = 0;
		Total_Diva_negatives = 0;
		
		
		forward_trans 	= 0;
		onward_trans 	= 0;
	
				
		severe=false;
		severe_indicator = 1;
				
		breakdown = false;
				
		// toggles at end of first breakdown
		breakfirstFlag = false;
		follow_up = false;
		VE6Mflag  = false;
		on_first_break = true;
		// Endsequence ensures that we bail out as soon as recurrence has occured
		endsequence = false;
		slaughter_recurr = false;
		// Introduce first infected animal at random within parish testing interval
		//double introducebTB = gsl_rng_uniform(gsl_r)*test_period;
				
		//double introducebTB = 0.0;
				
		short_interval = 60;
				
		//cout << "Testing Interval: " << short_interval << endl;
		//model.set_transmission(0.0,beta/364.0,0.0,1,0.0,0.0);
		first_test = true;
				
		// Toggles at start of first breakdown
		first_break = true;
	
		//cout << "Creeper \n";
		
		if(this_PTI==1){test_type=wht;}else{test_type=rht;}
		
		// Start with Clean Herd
		initialise_herd(Herd_Size,0,0,O);
		//set_time();
		//run(10*364,1.0,nosave);

		while(true)
		{

	   //cout << (Itot[0]+Otot[0]+OV1tot[0] + OV2tot[0] + SOV1tot[0] + SOV2tot[0] +Rtot[0]+SRtot[0]+SItot[0]+SOtot[0]+SItot[0]) << ' ' << infected[0].size() << ' ' << endl;

		// Can improve efficiency of simulation
		// by simulating proportion of 
		// false positive breakdowns through
		// relative rate of infectious pressure 
		// and standard specificity

		// Adds fixed difference to Entropy score
		// Look into again when you have time...
		//cout << first_break << endl;
		// Simulate to first (false positive) breakdown
		// or first infectious import			
		
		//cout << "First Break: " << first_break << ' ' << disease_free() << endl;
		// Restock herd to target size if depleted
		if(cows[0].size() < 1){initialise_herd(Herd_Size,0,0,O);/*run(10*364,1.0,nosave)*/;}

		
		if(first_break && disease_free())
		{
		
		// Calculate time to next infectious import
		double timport = - (long double) log(gsl_rng_uniform(gsl_r)) / (xinf*Herd_Size);
	    while(isinf(timport))
		{
			cout << "-Inf" << endl;exit(0);
		}
						
		int no_of_tests = 0;
						
		//cout << "Time to first import "	<< timport << " out of " << test_period << endl;						
		while(timport > test_period)
		{
			no_of_tests++;
			timport = timport - test_period;
			test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
		}
			
		//cout << "Number of Tests: " << no_of_tests << endl;
						
		int last_test = 0;
		int animal_tests = 0;
		
		double births = 0.0;
		switch(Seasonal_model)
	{
	case constant:
	//cout << "Constant " << endl;
	
	births = Herd_Size*(1.0/Life_Expectation[0]);
	break;
	case seasonal:
	//cout << "Seasonal " << endl;
	births = Life_Births[0];
	break;
	}	

		//cout << "Herd Size " << Herd_Size << " Births " << births << endl;
	
		for(int i=0; i < no_of_tests; i++)
		{
			// Exponential approximation
			animal_tests += (last_test = animals_to_test(gsl_r,Herd_Size,births,this_PTI));
			//cout << "Testsamps: " << i << ' ' << animal_tests << endl;
		}				
						
		//cout << "Timport: " << timport << ' ' << no_of_tests << ' ' << disease_free() << endl;
						
						
		double prob_fp = 1.0 - pow((1 - specificity[SICCT]),(double) animal_tests);
		//cout << animal_tests << ' ' << prob_fp << endl;
						if(gsl_rng_uniform(gsl_r) < prob_fp)
						{
							// False positive breakdown
							
							// Start with Clean Herd
							initialise_herd(Herd_Size,0,0,O);

							breakdown = true; follow_up = false; breakstart=return_time();	
							clear_cnt = 0;
							//cout << "FP! " << endl;
							//cout << 0 << ' ';
							confirmed=false;severe=false;severe_indicator = 1;
							
							//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
							//cout << "Skido: Short Interval: " << pickero << ' ' << no_of_SIT << ' ' << short_interval << endl;
							
							test_period=short_interval;
							test_type=all;
							
							// Movement restrictions //
							//set_transmission(0.0,beta/364.0,0.0,1,0.0,0.0,0.0);
							//cout << "Breakdown Start: " << return_time() << ' ' << breakfirstFlag << endl;
							
							first_break=false;Reactors_at_Start = reacto(specificity[SICCT], last_test);
							//cout << "Reactors at Start " << Reactors_at_Start << ' ' << last_test << endl;
							cull_FP(Reactors_at_Start);
							primary_breakstart = breakstart;
							Primary_Reactors+=Reactors_at_Start;
							//cout << "FP: " <<  prob_fp << ' ' << endl;
							
						}
						else
						{
							// Infectious import has happened before false positive
							test_period = timport;
							initialise_herd(Herd_Size,0,1,O);
							//cout << 1 << ' ';
							//cout << "Seeding: " << prob_fp << endl;
							//cout << "Test_Period: " << test_period << endl;
						}
					}
					
					
		double kerplunk = 0.0;
		double trigger_time=0.0;
	  
	// Simulate next period, test and cull
						//cout << "Test Period: " << test_period << endl;
						
						if((trigger_time = run(test_period,1.0,nosave)) != 0.0)
						{
							//Treat slaughterhouse incidents as confirmed breakdowns 
							confirmed = true;severe=true;severe_indicator = 2;
				
							// Slaughter-house Breakdown
							
							if(!breakdown)
							{
							
								breakdown = true; follow_up = false; 
								
							//	cout << "Trigger Time: " << trigger_time << ' ' << endl;
								breakstart=t;	
								
								clear_cnt = 0;
								
								short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
								//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							//	cout << "Short Interval Start: " << short_interval << endl;
								
								test_period=short_interval;
								test_type=all;
								
								// Movement restrictions //
								//set_transmission(0.0,beta/364.0,0.0,beta_K,K,0.0,0.0,0.0);
								//cout << "Slaughterhouse Breakdown Start: " << return_time() << endl;
								
								if(first_break){first_break=false;Reactors_at_Start = 0; Primary_Reactors++;slaughter_house = true; confirmed_ever = true;primary_breakstart = breakstart;}
								else{	
									// If recurrence //
									kerplunk = t - breakfirst;
									endsequence = true;
									slaughter_recurr = true;
									if(kerplunk <= 30*6){break6=true;break12=true;break24=true;break_recurr = kerplunk;}
									if(kerplunk <= 30*12){break12=true;break24=true;break_recurr = kerplunk;}
									if(kerplunk <= 30*24){break24=true;break_recurr = kerplunk;}
								}
								
							
							// Slaughter-house case triggers whole herd test (immediately) at severe interpretation
							//Fixed Confirmation Model
							policy_check = Test_Protocol_Per_Animal(0,this_PTI,test_type);
		
							
							
							//Empirical Confirmation Model
							//policy_check = Test_Protocol_Per_Animal(0,Sensitivity[0],Sensitivity[1],Specificity[0],Specificity[1],p_occult,p_occult,conf_modelP[index2(k,n,no_of_H)],this_PTI,test_type,Repeat_Sensitivity[0],Repeat_Sensitivity[1],conf_modelTheta[index2(k,n,no_of_H)]);
							
							
							Cull(0,policy_check,true,DiscoFlag);
							if(on_first_break)
							{Primary_Reactors += policy_check[severe_indicator];
							Total_Visits++;
							Total_Tests += policy_check[5];
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
							
							}
							
							}
							else
							{// Need to run model until next test is scheduled
							// bookkeeping slaughterhouse animals as they arise
							
							double running_time = 0.0;
						
							do{
							if((trigger_time = run(test_period-running_time,1.0,nosave)) != 0.0)
							{
							//cout << trigger_time << ' ' << running_time << ' ' << test_period << ' ' << Primary_Reactors << endl;
							//Switch to confirmed if not already so
							confirmed = true;severe=true;severe_indicator = 2;
							running_time = running_time + trigger_time;
							if(on_first_break){Primary_Reactors++;}
							} 
							else
							{
							 running_time = test_period;
							}
							}while(running_time < test_period); 
							
							// Slaughter-house cases ensure scheduled SIT at severe interpretation
							//Fixed Confirmation Model
							policy_check = Test_Protocol_Per_Animal(0,this_PTI,test_type);
							
							//Empirical Confirmation Model
							//policy_check = Test_Protocol_Per_Animal(0,Sensitivity[0],Sensitivity[1],Specificity[0],Specificity[1],p_occult,p_occult,conf_modelP[index2(k,n,no_of_H)],this_PTI,test_type,Repeat_Sensitivity[0],Repeat_Sensitivity[1],conf_modelTheta[index2(k,n,no_of_H)]);
						
							Cull(0,policy_check,true,DiscoFlag);
							if(on_first_break)
							{Primary_Reactors += policy_check[severe_indicator];
							Total_Visits++;
							Total_Tests += policy_check[5];
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
							
							
							}
							
							
							
							
							}
						}
						else
						{

					
							first_test=false;
							policy_check = Test_Protocol_Per_Animal(0,this_PTI,test_type);
							
							
							// If breakdown is confirmed remove reactors
							// at severe interpretation
							if(policy_check[6] > 0)
							{ 
								//cout << "Confirm" << endl;
                                confirmed=true;severe=true;severe_indicator = 2;
                            	if(on_first_break){confirmed_ever = true;}
                                Cull(0,policy_check,true,DiscoFlag);
                                Reactors_at_First_Test = policy_check[severe_indicator];
                                if(on_first_break)
                                {Primary_Reactors += policy_check[severe_indicator];
                                Total_Visits++;
							Total_Tests += policy_check[5];
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
                                   
                                }
							}
							else
							{
								if(!breakdown)
								{
								//cout << "Don't Confirm" << endl;
                                confirmed=false;severe=false;severe_indicator = 1;
                                }
                                Cull(0,policy_check,false,DiscoFlag);
                            	if(on_first_break)
                            	{Primary_Reactors += policy_check[severe_indicator];
                            	Total_Visits++;
							Total_Tests += policy_check[5];
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
                            	
                            	}
							
							}
								

						}
	
		
					
		int pval = accum_infected_cow_days(); 
		if(!isinf(pval)){ forward_trans += pval;/*cout << pval << endl;*/}
		if(!isinf(pval) && !breakdown){ onward_trans += pval;/*cout << pval << endl;*/}
					
		// Which test next?
					
		if(breakdown)
		{
						// Clear Test
						if(policy_check[severe_indicator] == 0)
						{
							if(confirmed && !severe)
							{
								
								// End movement restrictions
								//set_transmission(0.0,beta/364.0,0.0,1,xinf,0.0,0.0);
								clear_cnt++;
								if(clear_cnt >= must_clear)
								{
									breakdown = false; 
								
									test_period = VE6M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE6M-1))];
									
									//cout << "Scheduling VE6M" << endl;
									follow_up = true;
								
									VE6Mflag = true;
								
									test_type=all;
								
									if(!breakfirstFlag)
									{breakfirstFlag=true;breakfirst=return_time();primary_breaklength = return_time()-breakstart;confirmed_at_start = confirmed;primary_breakstart = breakstart;}
								
								severe = false;confirmed=false;clear_cnt = 0;
								if(on_first_break){on_first_break = false;}
								// Burden, all missed infection, remaining at end of breakdown
                                Burden =  burden(0);// This is correct and logs burden at end of primary breakdown, as simulation ends with recurrence.
                                //cout << "End Confirmed Break: " << return_time() << ' ' << (return_time()-breakstart) << endl;
                                // Perfect isolation sets external infectious pressure (xinf) to 0 at end of breakdown
                                if(perfect_isolation){cout << "ISOLATE!" << endl;xinf=0.0;}
                                
                                }
                                else
                                {
                                
                                // If we are on the final short interval test, use alternative duration if specified
                                    
                                    if(change_SIT & (clear_cnt >= 1))
                                    {
                                        short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                    
                                 //   cout << "Cleared " << clear_cnt << " tests, short interval of " << short_interval << endl;
                                    
                                    test_period = short_interval;
                                    test_type=all;
								}
                                
								
							}
							else if(!confirmed)
							{
								clear_cnt++;
								
								// If not confirmed, must only clear 1 test
								// if we have cleared required number of SIT, end movement restrictions
								if(clear_cnt >= must_clear)
								{
								
									breakdown = false; 
								
									test_period = VE6M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE6M-1))];
									//test_period = 6*30.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
									//cout << "Scheduling VE6M" << endl;
									follow_up = true;
									VE6Mflag  = true;
								
									test_type=all;
								
									// End movement restrictions
									//set_transmission(0.0,beta/364.0,0.0,1,xinf,0.0,0.0);
								
									if(!breakfirstFlag)
									{breakfirstFlag=true;breakfirst=return_time();primary_breaklength = return_time()-breakstart;primary_breakstart = breakstart;}
								
								   // cout << "End unconfirmed break: " << return_time() << ' ' << (return_time()-breakstart) << endl;		
									if(on_first_break){on_first_break = false;}
									Burden =  burden(0);// This is correct and logs burden at end of primary breakdown, as simulation ends with recurrence.
                                   
                                    // Perfect isolation sets externam infectious pressure (xinf) to 0 at end of breakdown
                                    if(perfect_isolation){xinf = 0.0;}
                                
									clear_cnt=0;
								}
								else
								{
								
								// If we are on the final short interval test, use alternative duration if specified
                                    
                                    if(change_SIT & (clear_cnt >= 1))
                                    {
                                        short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                   //cout << "Unconfirmed Breakdown, cleared " << clear_cnt << " tests, short interval of " << short_interval << endl;
                                    
                                    test_period = short_interval;
                                    test_type=all;
								
								}	
									
									
							}	
							else if(confirmed && severe)
							{
								// Coming off a confirmed breakdown
								// dropping down to standard interpretation
                                clear_cnt++;
								// If we can clear from severe interpretation
								// End Breakdown
								if(clear_from_severe && (clear_cnt >= must_clear))
								{
                               
                               }
                               else
                               {
                               
								// Coming off severe interpretation
                                    // If we can clear from severe, (must_clear-1) tests to go
                                    if(!clear_from_severe)
                                    {
                                        clear_cnt = 0;
                                    }
                                    
                                    // otherwise must clear must_clear tests at standard interperetation
                                  
                                    severe = false;severe_indicator = 1;
                                   // cout << "Confirmed Break: Drop severe: " << clear_cnt << endl;
                                    
                                    // Have already cleared a test by definition to get here
                                    // so switch to secondary frequency here if specified
                                    
                                    if(change_SIT)
                                    {
                                        short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                    
                                    //short_interval = fixed_SIT + gsl_ran_exponential(gsl_r2,364/2);
                                    
                                    test_period = short_interval;
                                    test_type=all;
                                    //cout << "New Short Interval: " << short_interval << endl;
		
				   //cout << "Confirmed Break: Drop severe: " << return_time() << ' ' << (return_time()-breakstart) << endl;
								
							}
							
						}
						}
						// Failed Test
						else
						{
							if(confirmed & !severe){severe=true;severe_indicator = 2;}
							short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
							//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							test_period=short_interval;
							test_type=all;
							clear_cnt=0;
							//cout << "Short Interval: " << short_interval << endl;
							
							//cout << "Failed Test: " << (return_time()-breakstart) << ' ' << breakfirstFlag << ' ' << policy_check[4] << ' ' << policy_check[1] << ' ' << policy_check[2] << ' ' << policy_check[5] << ' ' << Herd_Size  << endl;
						}
					}
		else
		{
						// Starting Breakdown
						// Breakdowns can only start under *standard* interpretation so check standard reactors here //
						// Avoids logical inconsistency if severe_reactors < standard_reactors... //
						if(policy_check[1] > 0) 
						{
							breakdown = true; follow_up = false; breakstart=return_time();	
							clear_cnt = 0;
							//int pickero = (int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1));
							
							short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
							//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							
							//cout << "Skido: Short Interval: " << ' ' << no_of_SIT << ' ' << short_interval << endl;
							
							
							test_period=short_interval;
							test_type=all;
							// Movement restrictions //
							//set_transmission(0.0,beta/364.0,0.0,1,0.0,0.0,0.0);
							//cout << "Breakdown Startled: " << return_time() << ' ' << breakfirstFlag << endl;
							
							if(first_break){first_break=false;Primary_Reactors = policy_check[severe_indicator];Reactors_at_Start = policy_check[severe_indicator];if(confirmed){confirmed_ever=true;} primary_breakstart = breakstart;slaughter_house = false;}
							else{	
								// If recurrence //
								kerplunk = return_time() - breakfirst;
								endsequence = true;
								if(VE6Mflag)
								{
									Reactors_at_VE6M =  policy_check[severe_indicator];
									
								}
								else
								{
									Reactors_at_VE12M =  policy_check[severe_indicator];
								}
								
								if(kerplunk <= 30*6){break6=true;break12=true;break24=true;break_recurr = kerplunk;}
								if(kerplunk <= 30*12){break12=true;break24=true;break_recurr = kerplunk;}
								if(kerplunk <= 30*24){break24=true;break_recurr = kerplunk;}
							}
							
							
						}
						// Moving through compulsory follow_up tests at 30*6 days and 30*12 days
						else if(follow_up)
						{
							if(VE6Mflag)
							{
								VE6Mflag = false;
								test_type=all;
								test_period = VE12M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE12M-1))];
								//test_period = 12*30.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
								Reactors_at_VE6M = 0;
								//cout << "Schedule VE-12M" << endl;
							}
							else
							{
								follow_up = false; 
								//test_period = this_PTI*364.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
								test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
								if(this_PTI==1){test_type=wht;}else{test_type=rht;}
								Reactors_at_VE12M = 0;
								endsequence = true;
								//cout << "Schedule Endgame" << endl;
							}
							
						}
						// Clear test, schedule next routine test
						else
						{
						
						test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
						//cout << "Next scheduled Test: " << test_period << ' ' << "Disease Status: " << disease_free() << endl; 
					}

					}		
					// Run for maximum of 20 years	
					//if((policy_check[4] == 0 && !breakdown))
					
					//if((policy_check[4] == 0 && !breakdown) || return_time() >= 364*20)
					// If we are 24 months past end of first breakdown or recurrence has occured - bail out.
					// or 
					// If first breakdown has lasted longer than 5 years - bail out
					// 
					//if((breakfirstFlag && (return_time()-breakfirst) >= 30*24) || (!first_break && !breakfirstFlag && (return_time()-breakstart) >= 364.0*5) )
					
					
					
					
					//|| ( (!first_break) && !breakfirstFlag && (return_time()-breakstart) >= 5*364))
					// or we are 36 months past the end of the first breakdown - bail out
					// this should be accounted for by endsequence
					//|| ( (!first_break) && breakfirstFlag && (return_time()-(primary_breakstart+primary_breaklength)) >= 3*364) )
					// If bailed out due to recurrence at follow up test or cleared VE12M
					// or If we are 5 years (or 200 reactors) past start of first breakdown and haven't cleared - bail out.  
					if( endsequence || ( (!first_break) && !breakfirstFlag && ((return_time()-breakstart) >= 100*364 || Primary_Reactors > 1000)))
					{
					if(on_first_break){Burden =  burden(0);primary_breaklength = return_time()-breakstart;}
					if(full_save){save_data();};/*cout<<"Endsequence " << first_break << ' ' << breakfirstFlag << ' ' << (return_time()-breakstart) << endl;*/break;}
	
	}
					
	
	
	}





void  bTBICBM::trial_uk_with_DIVA(bool DiscoFlag,int set_this_PTI, int Herd_Size,int Seeders,int must_clear = 1,bool clear_from_severe = false,bool change_SIT=false,int second_SIT=60, bool perfect_isolation=false,double batch_vacc=7.0, bool do_vaccinate=true,bool do_batch=false,bool suspend_in_break=false,double eligible_vacc = 0.0,bool retain = false, double trial_years = 5.0,double vacc_p = 1.0)
	{
	batch_time = batch_vacc;

	TargetHerdSize = Herd_Size;
	
	target_vaccination_p = vacc_p;
	
	int this_PTI = set_this_PTI;
	
	int clear_cnt = 0;
				
	vector<int> policy_check;			
	// Set current routine testing distribution
		switch(this_PTI)
		{
		// PTI 1
		case 1:
		Current_PTI = PTI1;
		Current_PTI_Num = no_of_PTI1;
		break;
		// PTI 2
		case 2:
		Current_PTI = PTI2;
		Current_PTI_Num = no_of_PTI2;
		break;
		// PTI 4
		case 4:
		Current_PTI = PTI4;
		Current_PTI_Num = no_of_PTI4;
		break;
		}
	
		//int test_period=this_PTI*364.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
				
		//int test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
		// Initial testing period of 6 months
		int routine_interval = 180;
		int test_period = routine_interval;
		
				
		int this_seed = 0;
				
		//test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
	
		set_time();
		
		//cout << model.return_time() << endl;
				
		breakstart=0.0;
		breakfirst=-1;
		break_recurr=0.0;
		primary_breaklength = 0.0;
		primary_breakstart  = 0.0;
		Reactors_at_Start = 0;
		Reactors_at_VE6M  = 0;
		Reactors_at_VE12M = 0;
		Primary_Reactors  = 0;
		Reactors_at_First_Test = 0;
		forward_trans = 0;
		break6=false;
		break12=false;
		break24=false;
		confirmed=false;
		confirmed_ever=false;
		confirmed_at_start = false;
		slaughter_house = false;
		Total_Diva_tests = 0;
		Total_Diva_negatives = 0;
		
		Burden			= 0;
		Total_Visits    = 0;
		Total_Tests		= 0;
		
		forward_trans 	= 0;
		onward_trans 	= 0;
			
		severe=false;
		severe_indicator = 1;
				
		breakdown = false;
				
		// toggles at end of first breakdown
		breakfirstFlag = false;
		follow_up = false;
		VE6Mflag  = false;
		on_first_break = true;
		// Endsequence ensures that we bail out as soon as recurrence has occured
		endsequence = false;
		slaughter_recurr = false;
		// Introduce first infected animal at random within parish testing interval
		//double introducebTB = gsl_rng_uniform(gsl_r)*test_period;
				
		//double introducebTB = 0.0;
				
		short_interval = 60;
				
		//cout << "Testing Interval: " << short_interval << endl;
		//model.set_transmission(0.0,beta/364.0,0.0,1,0.0,0.0);
		first_test = true;
				
		// Toggles at start of first breakdown
		first_break = true;
	
		//cout << "Creeper \n";
		test_type=all;
		//if(this_PTI==1){test_type=wht;}else{test_type=rht;}
		
		// Start with Clean Herd
		initialise_herd(Herd_Size,0,0,O);
		//vaccinate_off_schedule();
		//set_time();
		//run(10*364,1.0,nosave);

	   // Initialise vaccination schedule
	   // Annual revaccination with batch_time
	   // catchup of neonates
	   vacc_reset_schedule();
	   
	   if(do_vaccinate)
	   {
	   add_vacc(0,t + gsl_ran_flat(gsl_r,0,0.0001),true);
	
	   }
		if(do_batch)
		{
		   add_vacc(0,t+batch_time + gsl_ran_flat(gsl_r,0,0.0001),false);
		}
	   
	   
		while(true)
		{

	   //cout << (Itot[0]+Otot[0]+OV1tot[0] + OV2tot[0] + SOV1tot[0] + SOV2tot[0] +Rtot[0]+SRtot[0]+SItot[0]+SOtot[0]+SItot[0]) << ' ' << infected[0].size() << ' ' << endl;

		// Can improve efficiency of simulation
		// by simulating proportion of 
		// false positive breakdowns through
		// relative rate of infectious pressure 
		// and standard specificity

		// Adds fixed difference to Entropy score
		// Look into again when you have time...
		//cout << first_break << endl;
		// Simulate to first (false positive) breakdown
		// or first infectious import			
		
		//cout << "First Break: " << first_break << ' ' << disease_free() << endl;
		// Restock herd to target size if depleted
		//if(cows[0].size() < 1){initialise_herd(Herd_Size,0,0,O);/*run(10*364,1.0,nosave)*/;}

		
		//if(first_break && disease_free())
	//	{
		// Reseed herd
	//	initialise_herd(Herd_Size,0,1,O);
	//	}			
			
		double kerplunk = 0.0;
		double trigger_time=0.0;
	  
	// Simulate next period, test and cull
						//cout << "Test Period: " << test_period << endl;
						
						if((trigger_time = run(test_period,1.0,nosave)) != 0.0)
						{
							//Treat slaughterhouse incidents as confirmed breakdowns 
							confirmed = true;severe=true;severe_indicator = 2;
							
							// Slaughter-house Breakdown
							
							if(!breakdown)
							{
							
								breakdown = true; follow_up = false; 
								
							//	cout << "Trigger Time: " << trigger_time << ' ' << endl;
								breakstart=t;	
								
								clear_cnt = 0;
								
								//short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
								//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							//	cout << "Short Interval Start: " << short_interval << endl;
								
								test_period=short_interval;
								test_type=all;
								
								// Suspend Vaccination
	  							 if(do_vaccinate && suspend_in_break)
	  							 {
	  						    	vacc_reset_schedule();
	  							 }
								
								// Movement restrictions //
								//set_transmission(0.0,beta/364.0,0.0,beta_K,K,0.0,0.0,0.0);
								//cout << "Slaughterhouse Breakdown Start: " << return_time() << endl;
								
								if(first_break){first_break=false;Reactors_at_Start = 0; Primary_Reactors++;slaughter_house = true; confirmed_ever = true;primary_breakstart = breakstart;}
								else{	
									// If recurrence //
									kerplunk = t - breakfirst;
									endsequence = true;
									slaughter_recurr = true;
									if(kerplunk <= 30*6){break6=true;break12=true;break24=true;break_recurr = kerplunk;}
									if(kerplunk <= 30*12){break12=true;break24=true;break_recurr = kerplunk;}
									if(kerplunk <= 30*24){break24=true;break_recurr = kerplunk;}
								}
								
							
							// Slaughter-house case triggers whole herd test (immediately) at severe interpretation
							//Fixed Confirmation Model
							policy_check = DIVA_Protocol_Per_Animal(0,this_PTI,test_type,eligible_vacc);
			
							//Empirical Confirmation Model
							//policy_check = DIVA_Protocol_Per_Animal(0,Sensitivity[0],Sensitivity[1],Specificity[0],Specificity[1],p_occult,p_occult,conf_modelP[index2(k,n,no_of_H)],this_PTI,test_type,Repeat_Sensitivity[0],Repeat_Sensitivity[1],conf_modelTheta[index2(k,n,no_of_H)]);
							
						
							
							Cull(0,policy_check,true,DiscoFlag,this_PTI,retain);
							// For trial simulations, want total number of animal tests over trial
							Total_Tests += policy_check[5];
							if(on_first_break && breakdown)
							{Primary_Reactors += policy_check[severe_indicator];
							Total_Visits++;
							
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
							
							}
							
							}
							else
							{// Need to run model until next test is scheduled
							// bookkeeping slaughterhouse animals as they arise
							//cout << "Pinko Eyo" << endl;
							double running_time = 0.0;
						
							do{
							if((trigger_time = run(test_period-running_time,1.0,nosave)) != 0.0)
							{
							//cout << trigger_time << ' ' << running_time << ' ' << test_period << ' ' << Primary_Reactors << endl;
							//Switch to confirmed if not already so
							confirmed = true;severe=true;severe_indicator = 2;
							running_time = running_time + trigger_time;
							// For trial simulations, want total number of animal tests over trial
							Total_Tests += policy_check[5];
							if(on_first_break  && breakdown)
							{Primary_Reactors++;
							Total_Visits++;
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
							
							}
							} 
							else
							{
							 running_time = test_period;
							}
							}while(running_time < test_period); 
							
							// Slaughter-house cases ensure scheduled SIT at severe interpretation
							//Fixed Confirmation Model
							policy_check = DIVA_Protocol_Per_Animal(0,this_PTI,test_type,eligible_vacc);
							//Empirical Confirmation Model
							//policy_check = DIVA_Protocol_Per_Animal(0,Sensitivity[0],Sensitivity[1],Specificity[0],Specificity[1],p_occult,p_occult,conf_modelP[index2(k,n,no_of_H)],this_PTI,test_type,Repeat_Sensitivity[0],Repeat_Sensitivity[1],conf_modelTheta[index2(k,n,no_of_H)]);
							
							Cull(0,policy_check,true,DiscoFlag,this_PTI,retain);
							// For trial simulations, want total number of animal tests over trial
							Total_Tests += policy_check[5];
							if(on_first_break  && breakdown)
							{Primary_Reactors += policy_check[severe_indicator];
							Total_Visits++;
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
							
							}
							
							
							
							
							}
						}
						else
						{

							first_test=false;
							policy_check = DIVA_Protocol_Per_Animal(0,this_PTI,test_type,eligible_vacc);
							// For trial simulations, want total number of animal tests over trial
							Total_Tests += policy_check[5];
							// If breakdown is confirmed remove reactors
							// at severe interpretation
							
							if(policy_check[6] > 0)
							{ 
								//cout << "Confirm" << endl;
                                confirmed=true;severe=true;severe_indicator = 2;
                            	if(on_first_break  && breakdown){confirmed_ever = true;}
                                Cull(0,policy_check,true,DiscoFlag,this_PTI,retain);
                                Reactors_at_First_Test = policy_check[severe_indicator];
                                if(on_first_break  && breakdown){Primary_Reactors += policy_check[severe_indicator];
                                Total_Visits++;
		
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
                                
                                }
							}
							else
							{
								if(!breakdown)
								{
								//cout << "Don't Confirm" << endl;
                                confirmed=false;severe=false;severe_indicator = 1;
                                }
                                Cull(0,policy_check,false,DiscoFlag,this_PTI,retain);
                            	if(on_first_break  && breakdown){Primary_Reactors += policy_check[severe_indicator];
                            	Total_Visits++;
							
							if(severe_indicator==2)
							{
							Total_Diva_tests += policy_check[38];
							Total_Diva_negatives += policy_check[40];
							}
							else
							{
							Total_Diva_tests += policy_check[37];
							Total_Diva_negatives += policy_check[39];
							}
                            	}
							
							}
								

						}
	
		
					
		int pval = accum_infected_cow_days(); 
		if(!isinf(pval)){ forward_trans += pval;/*cout << pval << endl;*/}
		if(!isinf(pval) && !breakdown){ onward_trans += pval;/*cout << pval << endl;*/}
					
		// Which test next?
					
		if(breakdown)
		{
						// Clear Test
						if(policy_check[severe_indicator] == 0)
						{
							if(confirmed && !severe)
							{
								
								// End movement restrictions
								//set_transmission(0.0,beta/364.0,0.0,1,xinf,0.0,0.0);
								clear_cnt++;
								if(clear_cnt >= must_clear)
								{
									breakdown = false; 
								
									//test_period = VE6M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE6M-1))];
									test_period = 6*30;
									
									//test_period = 6*30.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
									//cout << "Scheduling VE6M" << endl;
									follow_up = true;
								
									VE6Mflag = true;
								
									test_type=all;
								
									if(!breakfirstFlag)
									{breakfirstFlag=true;breakfirst=return_time();primary_breaklength = return_time()-breakstart;confirmed_at_start = confirmed;primary_breakstart = breakstart;}
								
								severe = false;confirmed=false;clear_cnt = 0;
								if(on_first_break){on_first_break = false;}
								// Burden, all missed infection, remaining at end of breakdown
                                Burden =  burden(0);// This is correct and logs burden at end of primary breakdown, as simulation ends with recurrence.
                                //cout << "End Confirmed Break: " << return_time() << ' ' << (return_time()-breakstart) << endl;
                                // Perfect isolation sets external infectious pressure (xinf) to 0 at end of breakdown
                                if(perfect_isolation){cout << "ISOLATE!" << endl;xinf=0.0;}
                                // Restart Vaccination
	  							 if(do_vaccinate && suspend_in_break)
	  							 {
	  							 add_vacc(0,t + gsl_ran_flat(gsl_r,0,0.0001),true);
	   							}
	   							if(do_batch && suspend_in_break)
								{
		  						 add_vacc(0,t+batch_time + gsl_ran_flat(gsl_r,0,0.0001),false);
								}
	   								
                                }
                                else
                                {
                                
                                // If we are on the final short interval test, use alternative duration if specified
                                    
                                    if(change_SIT & (clear_cnt >= 1))
                                    {
                                        //short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        //short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                    
                                 //   cout << "Cleared " << clear_cnt << " tests, short interval of " << short_interval << endl;
                                    
                                    test_period = short_interval;
                                    test_type=all;
								}
                                
								
							}
							else if(!confirmed)
							{
								clear_cnt++;
								
								// If not confirmed, must only clear 1 test
								// if we have cleared required number of SIT, end movement restrictions
								if(clear_cnt >= must_clear)
								{
								
									breakdown = false; 
								
									//test_period = VE6M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE6M-1))];
									test_period = 6*30;
									
									//test_period = 6*30.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
									//cout << "Scheduling VE6M" << endl;
									follow_up = true;
									VE6Mflag  = true;
								
									test_type=all;
								
									// End movement restrictions
									//set_transmission(0.0,beta/364.0,0.0,1,xinf,0.0,0.0);
								
									if(!breakfirstFlag)
									{breakfirstFlag=true;breakfirst=return_time();primary_breaklength = return_time()-breakstart;primary_breakstart = breakstart;}
								
								   // cout << "End unconfirmed break: " << return_time() << ' ' << (return_time()-breakstart) << endl;		
									if(on_first_break){on_first_break = false;}
									Burden =  burden(0);// This is correct and logs burden at end of primary breakdown, as simulation ends with recurrence.
                                   
                                    // Perfect isolation sets externam infectious pressure (xinf) to 0 at end of breakdown
                                    if(perfect_isolation){xinf = 0.0;}
                                
									clear_cnt=0;
								}
								else
								{
								
								// If we are on the final short interval test, use alternative duration if specified
                                    
                                    if(change_SIT & (clear_cnt >= 1))
                                    {
                                        //short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        //short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                   //cout << "Unconfirmed Breakdown, cleared " << clear_cnt << " tests, short interval of " << short_interval << endl;
                                    
                                    test_period = short_interval;
                                    test_type=all;
								
								}	
									
									
							}	
							else if(confirmed && severe)
							{
								// Coming off a confirmed breakdown
								// dropping down to standard interpretation
                                clear_cnt++;
								// If we can clear from severe interpretation
								// End Breakdown
								if(clear_from_severe && (clear_cnt >= must_clear))
								{
                               
                               }
                               else
                               {
                               
								// Coming off severe interpretation
                                    // If we can clear from severe, (must_clear-1) tests to go
                                    if(!clear_from_severe)
                                    {
                                        clear_cnt = 0;
                                    }
                                    
                                    // otherwise must clear must_clear tests at standard interperetation
                                  
                                    severe = false;severe_indicator = 1;
                                   // cout << "Confirmed Break: Drop severe: " << clear_cnt << endl;
                                    
                                    // Have already cleared a test by definition to get here
                                    // so switch to secondary frequency here if specified
                                    
                                    if(change_SIT)
                                    {
                                        //short_interval = second_SIT;
                                    }
                                    else
                                    {
                                        //short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
                                    }
                                    
                                    //short_interval = fixed_SIT + gsl_ran_exponential(gsl_r2,364/2);
                                    
                                    test_period = short_interval;
                                    test_type=all;
                                    //cout << "New Short Interval: " << short_interval << endl;
		
				   //cout << "Confirmed Break: Drop severe: " << return_time() << ' ' << (return_time()-breakstart) << endl;
								
							}
							
						}
						}
						// Failed Test
						else
						{
							if(confirmed & !severe){severe=true;severe_indicator = 2;}
							//short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
							//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							test_period=short_interval;
							test_type=all;
							clear_cnt=0;
							//cout << "Short Interval: " << short_interval << endl;
							
							//cout << "Failed Test: " << (return_time()-breakstart) << ' ' << breakfirstFlag << ' ' << policy_check[4] << ' ' << policy_check[1] << ' ' << policy_check[2] << ' ' << policy_check[5] << ' ' << Herd_Size  << endl;
						}
					}
		else
		{
						// Starting Breakdown
						// Breakdowns can only start under *standard* interpretation so check standard reactors here //
						// Avoids logical inconsistency if severe_reactors < standard_reactors... //
						if(policy_check[1] > 0) 
						{
							breakdown = true; follow_up = false; breakstart=return_time();	
							clear_cnt = 0;
							//int pickero = (int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1));
							
							//short_interval = SIT[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_SIT-1))];
							//short_interval = 60.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
							
							//cout << "Skido: Short Interval: " << ' ' << no_of_SIT << ' ' << short_interval << endl;
							
							// Suspend Vaccination during breakdown?
	  					    if(do_vaccinate && suspend_in_break)
	   							{
	      							vacc_reset_schedule();
	   							}
							
							test_period=short_interval;
							test_type=all;
							// Movement restrictions //
							//set_transmission(0.0,beta/364.0,0.0,1,0.0,0.0,0.0);
							//cout << "Breakdown Startled: " << return_time() << ' ' << breakfirstFlag << endl;
							
							if(first_break){first_break=false;Primary_Reactors = policy_check[severe_indicator];Reactors_at_Start = policy_check[severe_indicator];if(confirmed){confirmed_ever=true;} primary_breakstart = breakstart;slaughter_house = false;}
							else{	
								// If recurrence //
								kerplunk = return_time() - breakfirst;
								endsequence = true;
								if(VE6Mflag)
								{
									Reactors_at_VE6M =  policy_check[severe_indicator];
									
								}
								else
								{
									Reactors_at_VE12M =  policy_check[severe_indicator];
								}
								
								if(kerplunk <= 30*6){break6=true;break12=true;break24=true;break_recurr = kerplunk;}
								if(kerplunk <= 30*12){break12=true;break24=true;break_recurr = kerplunk;}
								if(kerplunk <= 30*24){break24=true;break_recurr = kerplunk;}
							}
							
							
						}
						// Moving through compulsory follow_up tests at 30*6 days and 30*12 days
						else if(follow_up)
						{
							if(VE6Mflag)
							{
								VE6Mflag = false;
								test_type=all;
								//test_period = VE12M[(int) floor(gsl_rng_uniform(gsl_r)*(no_of_VE12M-1))];
								test_period = 6*30;
								//test_period = 6*30.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
								Reactors_at_VE6M = 0;
								//cout << "Schedule VE-12M" << endl;
							}
							else
							{
								follow_up = false; 
								//test_period = this_PTI*364.0 + gsl_ran_lognormal(gsl_r,2.7085,0.7);
								test_period = routine_interval;
								//if(this_PTI==1){test_type=wht;}else{test_type=rht;}
								test_type=all;
								Reactors_at_VE12M = 0;
							}
							
						}
						// Clear test, schedule next routine test
						else
						{
						
						//test_period = Current_PTI[(int) floor(gsl_rng_uniform(gsl_r)*(Current_PTI_Num-1))];
						test_period = routine_interval;
						//cout << "Next scheduled Test: " << test_period << ' ' << "Disease Status: " << disease_free() << endl; 
					}

					}		
					// Run for maximum of 20 years	
					//if((policy_check[4] == 0 && !breakdown))
					
					//if((policy_check[4] == 0 && !breakdown) || return_time() >= 364*20)
					// If we are 24 months past end of first breakdown or recurrence has occured - bail out.
					// or 
					// If first breakdown has lasted longer than 5 years - bail out
					// 
					//if((breakfirstFlag && (return_time()-breakfirst) >= 30*24) || (!first_break && !breakfirstFlag && (return_time()-breakstart) >= 364.0*5) )
					
					
					
					
					//|| ( (!first_break) && !breakfirstFlag && (return_time()-breakstart) >= 5*364))
					// or we are 36 months past the end of the first breakdown - bail out
					// this should be accounted for by endsequence
					//|| ( (!first_break) && breakfirstFlag && (return_time()-(primary_breakstart+primary_breaklength)) >= 3*364) )
					// If bailed out due to recurrence at follow up test or cleared VE12M
					// or If we are 5 years (or 200 reactors) past start of first breakdown and haven't cleared - bail out.  // Is this too short?
					if(return_time() >= trial_years*364)
					{
					
					if(on_first_break && breakdown){Burden =  burden(0);primary_breaklength = return_time()-breakstart;}
					if(on_first_break && !breakdown){Burden =  burden(0);}
					if(full_save){save_data();};/*cout<<"Endsequence " << first_break << ' ' << breakfirstFlag << ' ' << (return_time()-breakstart) << endl;*/break;}
	
	}				
	
	
	}  
	

// Prototype experiment, vaccinated in-contact animals	

// Simulate single phase only, do not initialise group
vector<int> bTBICBM::ExperimentWithDIVASinglePhase(int rFlag, int group_id, int phase_id, int Seeders, int R, int U, int V, int W, double trial_days = 364.0,double test_period=60,bool no_vacc=true)
	{
	int t_points = (int) trial_days/test_period;
	
   int Unegative	= 0;
   int Upositive	= 0;
   int Utot			= U;
   int U_cases		= 0;
   int UVL			= 0;
     
   int Rpositive	= 0;
   int Rnegative	= 0;
   int R_tot		= R;
   int R_cases		= 0;
   int RVL			= 0;
	   
   int Vpositive	= 0;
   int Vnegative	= 0;
   int Vtot			= V;
   int V_cases		= 0;
   int VVL			= 0;

   int Wpositive	= 0;
   int Wnegative	= 0;
   int Wtot			= W;
   int W_cases		= 0;
   int WVL			= 0;
	
	// Count positive and negative U,R,V,W (necessary for Phase II)
	
	vector<int> policy_check;			
	map<double,cow_t>::iterator iter; 
	
	 for( iter = cows[0].begin(); iter != cows[0].end(); iter++ ) 
		{
			if((iter->second).Seeder && (iter->second).Vaccinated_time>=0)
			{if((iter->second).Diva_ever){Wpositive++;}else{Wnegative++;}}
						
			if((iter->second).Seeder && (iter->second).Vaccinated_time<0)
			{if((iter->second).Diva_ever){Rpositive++;}else{Rnegative++;}}
			
			if(!(iter->second).Seeder && (iter->second).Vaccinated_time<0)
			{if((iter->second).Diva_ever){Upositive++;}else{Unegative++;}}
			
			if(!(iter->second).Seeder && (iter->second).Vaccinated_time>=0)
			{if((iter->second).Diva_ever){Vpositive++;}else{Vnegative++;}}
			
	    }	
	
	
   int Audit_Seeder = 0;
   int Audit_Control = 0;
   int Audit_Vaccinated = 0;	
	

	// Initialise vaccination schedule
	// Annual revaccination with batch_time
	// catchup of neonates
	// vacc_reset_schedule();

	for(int t=0; t <= t_points; t++)
	{
		//Revaccinated at annual intervals
	   if(t%6==0 && t!=0 && !no_vacc){vaccinate_off_schedule();}
	   if(t!=0){run(test_period,1.0,nosave);}
	   ClearStatus(0);
	   policy_check = DIVA_Protocol_Per_Animal(0,1,all,0.0);
	   
 	   U_cases = policy_check[41];
 	   R_cases = policy_check[42];
 	   V_cases = policy_check[43];
 	   W_cases = policy_check[44];
 	    	    	   	   
	   Upositive += U_cases;
	   Rpositive += R_cases;
	   Vpositive += V_cases;
	   Wpositive += W_cases;
	   
	   Unegative -= U_cases;
	   Rnegative -= R_cases;
	   Vnegative -= V_cases;
	   Wnegative -= W_cases;
	         
	   UVL = 0;
	   VVL = 0;
	   RVL = 0;
	   WVL = 0;

	   Audit_Seeder = 0;
   	   Audit_Control = 0;
   	   Audit_Vaccinated = 0;	

	   // PM for lesions at end of experiment
	   for( iter = cows[0].begin(); iter != cows[0].end(); iter++ ) 
		{
			if((iter->second).Control){Audit_Control++;}
			if((iter->second).Seeder){Audit_Seeder++;}
			if((iter->second).Vaccinated_time>=0){Audit_Vaccinated++;}
						
			bool slaughter_flag = slaughter_house_test(iter->second);
			if(slaughter_flag)
			{	if((iter->second).Control)
					{
					if((iter->second).Seeder){RVL++;}else{UVL++;}
					}
					else
					{
					if((iter->second).Seeder){WVL++;}else{VVL++;}
					}
	    	}
	    }
	    
	   //Whole_Herd_Cull(0,policy_check,true,DiscoFlag);
 	   
 	   //Vnegative = Vnegative - policy_check[42];
 	   
 	   cout << rFlag << ' ' << group_id << ' ' << phase_id << ' ' << t << ' ' 
 	   << R_tot << ' ' << Rpositive << ' ' << Rnegative << ' ' << R_cases << ' ' << RVL << ' ' 
 	   << Vtot << ' ' << Vpositive << ' ' << Vnegative << ' ' << V_cases << ' ' << VVL << ' '
 	   << Utot << ' ' << Upositive << ' ' << Unegative << ' ' << U_cases << ' ' << UVL << ' '
	   << Wtot << ' ' << Wpositive << ' ' << Wnegative << ' ' << W_cases << ' ' << WVL << ' ' << cows[0].size() << ' ' << Audit_Seeder << ' ' << Audit_Control << ' ' << Audit_Vaccinated << endl;
	
 	}

	   return(policy_check);
	}


		
	bTBICBM::~bTBICBM()
	{

		delete [] Life_Expectation;
		delete [] Life_var1;
		delete [] Occupancy_Expectation;
		delete [] Occupancy_var1;

		delete [] Life_Births;
		delete [] FOI;

		delete [] cows;
		delete [] infected;

		delete [] SIT;
		delete [] VE6M;
		delete [] VE12M;
		delete [] PTI1;
		delete [] PTI2;
		delete [] PTI4;

		delete [] Stot;	delete [] SStot;
		delete [] Otot; delete [] SOtot;
		delete [] Rtot; delete [] SRtot;
		delete [] Itot; delete [] SItot;
		delete [] V1tot; delete [] SV1tot;
		delete [] V2tot; delete [] SV2tot;
		delete [] OV1tot; delete [] SOV1tot;
		delete [] OV2tot; delete [] SOV2tot;
		delete [] RVtot; delete [] SRVtot;
		delete [] IVtot; delete [] SIVtot;

		delete [] Darth_Age_Distro;
		delete [] Darth_cows;
		delete [] DarthHerdOffsets;
		delete [] DarthHerdSamples;
		delete [] DarthHerdSize;
		delete [] DarthBindex;
		delete [] DarthPTI;
		
		gsl_histogram_free(AgeReactors);
		gsl_histogram_free(AgeCReactors);
		gsl_histogram_free(AgeSReactors);
		gsl_histogram_free(AgeSlaughter);
		
		for(int i = 0; i < H; i++)
		{
			
			if(Agefile[i])
			{
			Agefile[i].close();	
			}
			
			if(Statusfile[i])
			{
			Statusfile[i].close();
			}
			
			if(Lifefile[i])
			{
			Lifefile[i].close();
			}
		
		}

				
		for(int i = 0; i < H; i++)
		{
				
			if(disclose[i])
			{
				disclose[i].close();
			}
				if(Reactorfile[i])
			{
				Reactorfile[i].close();
			}
			if(Confirmedfile[i])
			{
				Confirmedfile[i].close();
			}
			
				if(Slaughterfile[i])
			{
				Slaughterfile[i].close();
			}
			
			if(SIndividualfile[i])
			{
				SIndividualfile[i].close();
			}
			
			
		}
		
		//delete [] disclose;
		
		
	}
	
	

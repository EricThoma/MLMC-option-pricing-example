#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <time.h>


//make sim_data struct instead
struct eq_data {
	double sigma; // volatility
    double mu; // drift
    double T;
};

void run_MC_sim(gsl_rng * seed, struct eq_data* params);
void run_MLMC_sim(gsl_rng * seed, struct eq_data* params);

double estimator_MLMC(gsl_rng * seed, double epsilon, struct eq_data * params);
double estimator_MC(gsl_rng * seed, int K, struct eq_data * params, int N);
double estimator_exact_MC_EU(gsl_rng * seed, int K, struct eq_data * params);

void gen_gaussian_vec(gsl_rng * seed, int N, double var, double* res);

double payoff_euro(int N, double* path, double strike);

void solve_EM(double * gauss_vec, double* path, struct eq_data * params, int N);
double eval_Y_diff(gsl_rng * seed, struct eq_data * params, double* path1, double* path2, double* gauss_vec1, double* gauss_vec2, int N, int factor);

int main (void) {
	struct eq_data* params;
	double sigma = 1.0;
	double mu = 0.5;

	int N = 500;
	double T = 1.0;
	gsl_rng * seed;
	const gsl_rng_type * TYPE;
	gsl_rng_env_setup();
	TYPE = gsl_rng_default;
	seed = gsl_rng_alloc (TYPE);

	params = (struct eq_data*)malloc(sizeof(struct eq_data));
	params->sigma = sigma;
	params->mu = mu;
	params->T = T;

	run_MC_sim(seed, params);


	return 0;
}

void run_MC_sim(gsl_rng * seed, struct eq_data* params){
	printf("Pricing option...\n");

	int ntrials = 1000;
	int kcur[ntrials];
	int ncur[ntrials];
	int msecs[ntrials];
	double vals[ntrials];

	kcur[0] = 16*16*2;
	ncur[0] = 16*2*2;
	FILE* output = (FILE*)fopen("./data_MC.dat", "w");

	clock_t start, diff;

	for(int i = 0; i < ntrials; i++){
		start = clock();
		double total = 0.0;
		for(int j = 0; j < 40; j++){
			total += fabs(estimator_MC(seed, kcur[i], params, ncur[i]) - 0.8873); // 0.8873 is true val
		}
		vals[i] = total/40.0;
		diff = clock() - start;
    	int msec = diff * 1000 / CLOCKS_PER_SEC;
    	msecs[i] = msec;

		printf("%i %i %i %f %i\n", i, kcur[i], ncur[i], vals[i], msec);

		fwrite(&(kcur[i]), sizeof(int), 1, output);
		fwrite(&(ncur[i]), sizeof(int), 1, output);
		fwrite(&(vals[i]), sizeof(double), 1, output);
		fwrite(&(msecs[i]), sizeof(int), 1, output);

		if(msec > 240*1000){ // cutoff at 240 seconds per sim
			break;
		}

		if(i+1 < ntrials){
			kcur[i+1] = kcur[i]*2;//ceil((double)kcur[i]*2);
			ncur[i+1] = ceil((double)ncur[i]*sqrt(2));
		}


	}
}


// warning: writes to file
void run_MLMC_sim(gsl_rng * seed, struct eq_data* params){

	printf("Pricing option...\n");


	int ntrials = 1000;
	double epsilon[ntrials];
	int msecs[ntrials];
	double vals[ntrials];
	epsilon[0] = 0.02;


	FILE* output = (FILE*)fopen("./data_MLMC.dat", "w");

	clock_t start, diff;

     for(int i = 0; i < ntrials; i++){
    	printf("i = %i / %i\n", i, ntrials);
    	printf("eps = %f\n", epsilon[i]);
    	start = clock();
    	double total = 0.0;
		for(int j = 0; j < 40; j++){
			total += fabs(estimator_MLMC(seed, epsilon[i], params) - 0.8873);
		}
    	double val = total/40.0;
    	
    	vals[i] = val;
    	diff = clock() - start;
    	int msec = diff * 1000 / CLOCKS_PER_SEC;
    	msecs[i] = msec;
    	printf("val= %f\n", val);
    	printf("Time taken %d milliseconds\n", msec);
    	fwrite(&(epsilon[i]), sizeof(double), 1, output);
		fwrite(&(msecs[i]), sizeof(int), 1, output);
		fwrite(&(vals[i]), sizeof(double), 1, output);

		if(i +1 < ntrials){
			epsilon[i+1] = epsilon[i]/(sqrt(2));
		}

		if(msec > 240*1000){
			break;
		}
    }
}

// use analytical solution and MC method
// specialized to European option with strike 1 and initial val X(0)=1.0
double estimator_exact_MC_EU(gsl_rng * seed, int K, struct eq_data * params){
	double total = 0.0;
	double val;
	double sigsquare = (params->sigma * params->sigma)/2.0;
	for(int i = 0; i < K; i++){
		val = exp((params->sigma) * gsl_ran_gaussian(seed, params->T) 
			       + (params->mu - sigsquare)*params->T);
		total += fmax(val - 1.0, 0.0); // Euro option with strike 1
	}

	return total / ((double)K);
}

//try to estimate within varepsilon-ish using MLMC
double estimator_MLMC(gsl_rng * seed, double epsilon, struct eq_data * params){
	int L = 1;


	//finicky parameter
	L = ceil(fabs(log2(epsilon*epsilon))/2.0);

	double lambda = 5*((double)L)/(epsilon*epsilon);

	int* ncurs = (int*)malloc(sizeof(int)*L);


	ncurs[0] = 64;
	for(int i = 1; i < L; i++){
		ncurs[i] = ncurs[i-1]*2;
	}


	int* kcurs = (int*)malloc(sizeof(int)*L);
	for(int i = 0; i < L; i++){
		kcurs[i] = (int)(lambda/((double)ncurs[i])); // lambda * estimated cost

		if(kcurs[i] == 0){
			L = i-1;
			break;

		}
	}



	double * Ys = (double*)malloc(sizeof(int)*L);

	for(int ell = 0; ell < L; ell++){
		Ys[ell] = 0;
	}

	fflush(stdout);

	Ys[0] = estimator_MC(seed, kcurs[0], params, ncurs[0]);




	for(int ell = 1; ell < L; ell++){
		double* paths1 = (double*)malloc(sizeof(double)*ncurs[ell]+1);
		double* paths2 = (double*)malloc(sizeof(double)*ncurs[ell-1]+1);
		double* gauss_vec1 = (double*)malloc(sizeof(double)*ncurs[ell]);
		double* gauss_vec2 = (double*)malloc(sizeof(double)*ncurs[ell-1]);

		fflush(stdout);
		for(int i = 0; i < kcurs[ell]; i++){

			double ydiff = eval_Y_diff(seed, params, 
							paths1, paths2, gauss_vec1, gauss_vec2, ncurs[ell-1], 2);

			
			Ys[ell] += ydiff;

			
		}
		Ys[ell] = Ys[ell] / ((double)kcurs[ell]);
		free(paths1);
		free(paths2);
		free(gauss_vec1);
		free(gauss_vec2);
	}

	double res = 0.0;
	for(int ell = 0; ell < L; ell++){
		res += Ys[ell];
	}

	free(Ys);

	free(ncurs);
	free(kcurs);

	return res;
}


// K trials
double estimator_MC(gsl_rng * seed, int K, struct eq_data * params, int N){
	double total = 0.0;
	double val;

	double* gauss_vec = (double *) malloc(sizeof(double)*N);
	double* path = (double *) malloc(sizeof(double)*N + 1);

	for(int i = 0; i < K; i++){
		gen_gaussian_vec(seed, N, (params->T)/((double)N), gauss_vec);
		solve_EM(gauss_vec, path, params, N);
		val = payoff_euro(N, path, 1.0); //strike at 1.0
		total += val;
	}

	return total/((double)K);
}

//
double payoff_euro(int N, double* path, double strike){
	// utility is max(x - s,0) where s is strike price and x is spot price
	return fmax(path[N] - strike, 0.0);
}



// result is pre-allocated to size N
void gen_gaussian_vec(gsl_rng * seed, int N, double var, double* res){
	for(int i = 0; i < N; i++){
		res[i] = sqrt(var)*gsl_ran_gaussian(seed, 1);
	}
	return;
}

void solve_EM(double * gauss_vec, double* path, struct eq_data * params, int N){
	double step = (params->T)/((double)(N));
	// set initial data
	path[0] = 1.0;
	double mudt = (params->mu) * step;
	
	for(int i=0; i < N; i++){
		path[i+1] = path[i] + path[i]*(mudt + (params->sigma) * gauss_vec[i]);
	}

	return;
}

// internal sum part of MLMC
double eval_Y_diff(gsl_rng * seed, struct eq_data * params, double* path1, double* path2, double* gauss_vec1, double* gauss_vec2, int N, int factor){

	// pre-generate all used gaussians at the finest resolution
	gen_gaussian_vec(seed, N*factor, (params->T)/((double)N*factor), gauss_vec1);

	// reduce resolution
	for(int i = 0; i < N; i++){
		gauss_vec2[i] = 0;
		for(int j = 0; j < factor; j++){
			gauss_vec2[i] += gauss_vec1[i*factor + j];
		}
	}

	// now run both simulations
	solve_EM(gauss_vec1, path1, params, N*factor);
	solve_EM(gauss_vec2, path2, params, N);

	double y1 = payoff_euro(N*factor, path1, 1.0);
	double y2 = payoff_euro(N, path2, 1.0);

	return y1 - y2;
}


#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

#define N_ENER_CONV  4096  // number of bins for the convolution, not that it needs to follow 2^N because of the FFT
#define EMIN_RELXILL 0.00035  // minimal energy of the convolution (in keV)
#define EMAX_RELXILL 2000.0

void rebin_spectrum(double* ener, double* flu, int nbins, double* ener0, double* flu0, int nbins0);
void get_log_grid(double* ener, int n_ener, double emin, double emax);
// double get_local_spec(double **xill_spec, double theta, int i);
void get_local_spec(double **xill_spec, double *local_spec, double ctheta, int nbins);
double trapz(double* ener, double* photar);

double Incl[] = {18.194874, 31.788328, 41.409622, 49.458397, 56.632984, 63.256313, 69.51268, 75.522484, 81.37307, 87.13402};
int xill_n_incl = 10, xill_n_ener = 2999;
const int number_of_bins = 5000;
int len_count = 0;
double iobs;

int main(int argc, char *argv[]) {

	double iobs_deg;
	double lxi;
	double max_energy = 300;  // in keV
	double alpha = -3.0;
	double rstep;

	// for (int i = 0; i < 10; i++) {
	// 	Incl[i] = 0.95 - (double)i * 0.1;
	// 	// printf(" %lf\n", Incl[i]);
	// }
    char ph_file[300], xillver_file[300], out_file[300];

    // Input parameters: photons_data_filename, xill_spec_filename, final_filename, incl, logxi, alpha	

    strcpy(ph_file, argv[1]);
    strcpy(xillver_file, argv[2]);
    strcpy(out_file, argv[3]);
    iobs_deg = atof(argv[4]);
    lxi = atof(argv[5]);
    alpha = atof(argv[6]);
    rstep = atof(argv[7]);
	iobs = iobs_deg * M_PI / 180.0;


    /*------------Reading xillver data file------------*/

	double **xill_spec, *xill_energy;

	xill_spec = (double**) malloc(xill_n_incl*sizeof(double*));
	for (int ii = 0; ii < xill_n_incl; ii++) {
		xill_spec[ii] = (double *) malloc(xill_n_ener * sizeof(double));
	}
	xill_energy = (double *) malloc((xill_n_ener + 1) * sizeof(double));


	FILE *xillver_data = fopen(xillver_file, "r");

    if (xillver_data == NULL){
        printf("Error Reading xillver data file!\n");
        exit (2);
    }
    
    for (int i = 0; i < xill_n_incl ; i++) {
    	for (int j = 0; j < xill_n_ener; j++) {
	        fscanf(xillver_data, "%lf %lf\n", &xill_energy[j], &xill_spec[i][j]);
	        // if (i == 0 && j < 10)
	        //     printf("%lf %lf\n", xill_energy[j], xill_spec[i][j]);
	        // xill_spec[i][j] *= 0.5*cos(iobs_deg*M_PI/180);  // comment out
	    }
    }
    xill_energy[xill_n_ener] = xill_energy[xill_n_ener-1] + (xill_energy[xill_n_ener-1] - xill_energy[xill_n_ener-2]);

    fclose(xillver_data);

    /*------------End of reading xillver data file------------*/


    // double E_obs[number_of_bins+1], N_obs[number_of_bins], local_numflux[number_of_bins];
    double *E_obs, *N_obs, *local_numflux, *n_local_numflux;
    E_obs = (double *) malloc((number_of_bins + 1) * sizeof(double));
    N_obs = (double *) malloc((number_of_bins) * sizeof(double));
    local_numflux = (double *) malloc((number_of_bins) * sizeof(double));
    n_local_numflux = (double *) malloc((number_of_bins) * sizeof(double));
    
    for (int i = 0; i < number_of_bins; i++)
    {
		N_obs[i] = 0.0;
		E_obs[i] = pow(10, (i*(log10(max_energy) + 1)/4999 - 1));	
        local_numflux[i] = 0.0;
    }
    E_obs[number_of_bins] = E_obs[number_of_bins-1] + (E_obs[number_of_bins-1] - E_obs[number_of_bins-2]) / 2.0;

	// double **local_spec = (double**) malloc(xill_n_incl*sizeof(double*));
	double *local_spec = (double *) malloc(number_of_bins * sizeof(double));
	double *xill_spec_interpolated = (double *) malloc(xill_n_ener * sizeof(double));


    /*----------------- Reading photon data file and convolving with xillver spectra ---------------*/

    FILE *ph_data;
    ph_data = fopen(ph_file, "r");
    if (ph_data == NULL){
        printf("Error Reading photons data File!\n");
        exit (1);
    }
	
	int end_of_file = 0;
	double xobs, yobs, r_disk, gfactor, cosne;
	int index;

    double r, g, rsquare, elocal, llocal, pp, qq, ctheta;
    // double rstep = 1.0001;
    double rstep2 = (rstep - 1)/rstep;
    int k;
    int i;


	while (!end_of_file) {
		fscanf(ph_data, "%d %lf %lf %lf %lf %lf\n", &index, &xobs, &yobs, &r_disk, &gfactor, &cosne);

		r = r_disk;
		g = gfactor;
		ctheta = cosne;
		// ctheta = cos(iobs);
		rsquare = xobs*xobs + yobs*yobs;
		get_local_spec(xill_spec, xill_spec_interpolated, ctheta, xill_n_ener);
		rebin_spectrum(E_obs, local_spec, number_of_bins, xill_energy, xill_spec_interpolated, xill_n_ener);

	    /*--------conv start------------*/
		for (int j = 0; j < number_of_bins; j++)
		{
			elocal = E_obs[j];
			// llocal = local_spec[j] / pow(10,lxi) / E_obs[j];
			llocal = local_spec[j] / E_obs[j];
			// printf(" %d, %d, %lf\n", i, j, llocal);

			local_numflux[j] += llocal;
			n_local_numflux[j] += 1.0;
			pp = g*elocal;
			qq = (g*g)*pow(r, alpha);

			if(pp >= 0.1 && pp<= max_energy)
			{
				k = floor(4999*(log10(pp) + 1)/(log10(max_energy) + 1));
				N_obs[k] = N_obs[k] + llocal*qq*rsquare*rstep2;
			}
		}

		if (len_count < 10)
			printf(" i=%d r=%lf, g=%lf, cosne=%lf, r2=%lf\n", len_count, r, g, ctheta, rsquare);

		if (feof(ph_data)) end_of_file = 1;
		len_count++;
	}


    printf(" Number of photons in the photon data file = %d\n", len_count);
    fclose(ph_data);
    
    printf("Finished importing the photon data and convolution!\n");


	for (int j = 0; j < number_of_bins; j++) {
		local_numflux[j] /= n_local_numflux[j];
	}

	double Q0, Q1;
	Q0 = trapz(E_obs, local_numflux);
	Q1 = trapz(E_obs, N_obs);

	printf("Q0 = %lf,\tQ1 = %lf\n", Q0, Q1);

	/*------------writing the full spectrum------------------*/

	// printf(" 3\n");
	FILE *foutput;
	foutput = fopen(out_file, "w");

	for (int i = 0; i < number_of_bins - 1; i++)
	{
	    fprintf(foutput, "%.10lf %.10lf\n", E_obs[i], N_obs[i]*Q0/Q1);
	}
	fclose(foutput);

	free(E_obs);
	free(N_obs);
	free(local_numflux);
	free(n_local_numflux);
	for (int ii = 0; ii < xill_n_incl; ii++) {
		free(xill_spec[ii]);
	}
	free(xill_spec);
	free(xill_energy);
	free(local_spec);
	free(xill_spec_interpolated);

	return 0;
}




/** rebin spectrum to a given energy grid
 *  length of ener is nbins+1       **/

void rebin_spectrum(double* ener, double* flu, int nbins, double* ener0, double* flu0, int nbins0){

	int ii; int jj;
	int imin = 0;
	int imax = 0;

	for (ii=0; ii<nbins; ii++){

		flu[ii] = 0.0;

		/* check of the bin is outside the given energy range */
		if ( (ener0[0] <= ener[ii+1]) && (ener0[nbins0] >= ener[ii]) ){

			/* need to make sure we are in the correct bin */
			while ( ener0[imin]<=ener[ii] && imin<=nbins0){
				imin++;
			}
			// need to set it back, as we just crossed to the next bin
			if (imin>0){
				imin--;
			}
			while ( (ener0[imax]<=ener[ii+1] && imax<nbins0)){
				imax++;
			}
			if (imax>0){
				imax--;
			}

			double elo = ener[ii];
			double ehi = ener[ii+1];
			if (elo < ener0[imin]) elo=ener0[imin];
			if (ehi > ener0[imax+1]) ehi=ener0[imax+1];

			if (imax==imin){
				flu[ii] = (ehi-elo) / (ener0[imin+1] - ener0[imin]) * flu0[imin];
			} else {

				double dmin=(ener0[imin+1]-elo)/(ener0[imin+1]-ener0[imin]);
				double dmax=(ehi-ener0[imax])/(ener0[imax+1]-ener0[imax]);

				flu[ii] += flu0[imin]*dmin + flu0[imax]*dmax;

				for (jj=imin+1; jj <= imax-1; jj++) {
					flu[ii] += flu0[jj];
				}

			}

		}

	}
	// printf(" flu0 %d %lf %lf %lf\n", nbins-1, flu[nbins-1], flu[nbins-2], flu[nbins-3]);
	flu[nbins - 1] = flu[nbins-2] - 0.9*(flu[nbins-3] - flu[nbins-2]);
	// printf(" flu1 %d %lf %lf %lf\n", nbins-1, flu[nbins-1], flu[nbins-2], flu[nbins-3]);

}

// double get_local_spec(double **xill_spec, double ctheta, int i) {
// 	int ind;
// 	double c2 = ctheta;
// 	if (ctheta < 0.05) {
// 		ind = 8;
// 		ctheta = 0.05;
// 	}
// 	else if (ctheta > 0.95) {
// 		ctheta = 0.95;
// 		ind = 0;
// 	}
// 	else {
// 		// int n, n2;
// 		// n = (int)((ctheta - 0.05) / 0.1) + 1;
// 		// n2 = 10 - ((int)((ctheta - 0.05) / 0.1) + 1);
// 		ind = 10 - ((int)((ctheta - 0.05) / 0.1) + 1) - 1;
// 	}



// 	double frac = (acos(ctheta) - acos(Incl[ind])) / (acos(Incl[ind + 1]) - acos(Incl[ind]));
// 	assert(frac <= 1);
// 	// if (frac > 1 || frac < 0)
// 	// 	printf(" frac=%lf, ind=%d %lf %lf %lf %lf %lf\n", frac, ind, ctheta, c2, acos(ctheta), acos(Incl[ind]), acos(Incl[ind+1]));
// 	assert(frac >= 0);
// 	return xill_spec[ind][i] + frac * (xill_spec[ind + 1][i] - xill_spec[ind][i]);

// }

// double Incl[] = {18.194874, 31.788328, 41.409622, 49.458397, 56.632984, 63.256313, 69.51268, 75.522484, 81.37307, 87.13402};

void get_local_spec(double **xill_spec, double *local_spec, double ctheta, int nbins) {
	double theta = fabs(acos(fabs(ctheta))) * 180.0 / M_PI;
	int ind = 0;
	if (theta <= 18.194874) {
		theta = 18.194874;
		ind = 0;
	}
	else if (theta >= 87.13402) {
		theta = 87.13402;
		ind = xill_n_incl - 2;
	}
	else {
		while(theta > Incl[ind])
			ind++;
		ind--;
	}

	assert(theta >= Incl[ind]);
	assert(theta <= Incl[ind+1]);
	double frac = (theta - Incl[ind]) / (Incl[ind+1] - Incl[ind]);

	for(int i = 0; i < nbins; i++)
		local_spec[i] = (xill_spec[ind][i] + frac * (xill_spec[ind + 1][i] - xill_spec[ind][i]));
		// local_spec[i] = (xill_spec[ind][i] + frac * (xill_spec[ind + 1][i] - xill_spec[ind][i])) * 0.5*cos(iobs);

	// return xill_spec[ind][i] + frac * (xill_spec[ind + 1][i] - xill_spec[ind][i]);

}

/* get a logarithmic grid from emin to emax with n_ener bins  */
void get_log_grid(double* ener, int n_ener, double emin, double emax){
	int ii;
	for (ii=0; ii<n_ener; ii++){
		ener[ii] = 1.0*ii / (n_ener-1) * ( log(emax) - log(emin)) + log(emin);
		ener[ii] = exp(ener[ii]);
	}
}



double trapz(double* ener, double* photar)
{
    int i;
    double Q;
    Q = 0.0;
    for(i = 0; i < number_of_bins - 1; i++)
    {
        Q = Q + photar[i]*(ener[i + 1] - ener[i]);
    }
    
    return Q;
}

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "../tempo2.h"
#include "/usr/include/fitsio.h"
#include <time.h>

using namespace std;   		/* Is this required for a plugin ? Yes, for Linux */

void extra_delays_fermi(pulsar *psr,int npsr);
void clock_corrections_fermi(pulsar *psr,int npsr);
void ephemeris_routines_fermi(pulsar *psr,int npsr);
void formBatsAll_fermi(pulsar *psr,int npsr);

static char random_letter(int is_cap);
static char random_number();
static void random_string(int length, char *str);

// met2mjd: conversion of Mission Elapsed Time in TT to MJD
double met2mjd(double met)
{
	double mjd_ref = 51910.0007428703703703703;	
	double mjd = mjd_ref + met / 86400.;
	
	return mjd;
}

// mjd2met : conversion of MJD in TT to Mission Elapsed Time
double mjd2met(double mjd)
{
	double mjd_ref = 51910.0007428703703703703;
	double met = 86400. * (mjd - mjd_ref);
	
	return met;
}

double inner_product(double vect_x[], double vect_y[])
{
  	return vect_x[0]*vect_y[0] + vect_x[1]*vect_y[1] + vect_x[2]*vect_y[2];
}

void outer_product(double vect_x[], double vect_y[], double vect_z[])
{
  	vect_z[0] = vect_x[1]*vect_y[2] - vect_x[2]*vect_y[1];
  	vect_z[1] = vect_x[2]*vect_y[0] - vect_x[0]*vect_y[2];
  	vect_z[2] = vect_x[0]*vect_y[1] - vect_x[1]*vect_y[0];
}

extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
	// Provide basic information
	printf("\n");
	printf("------------------------------------------\n");
	printf("Output interface:    photons\n");
	printf("Author:              Anne Archibald\n");
	printf("                     based on the fermi plugin v5.2\n");
	printf("                     by Lucas Guillemot\n");
	printf("Updated:             ?? November 2010\n");
	printf("Version:             0.0\n");
	printf("------------------------------------------\n");
	printf("\n");

	int i,j,k,l,event = 0,event2;

	char parFile[1][MAX_FILELEN];
	char timFile[1][MAX_FILELEN];
	char temptim[9];
	char gr[256]="/xs";					// default graphical output is x-window

	srand((unsigned)time(NULL));
	temptim[8] = '\0';
	random_string(8,temptim);
	strcpy(timFile[0],temptim);
	strcat(timFile[0],".tim");

	int nbins  = 20;					// default number of bins for the output phase histogram
	int Hbins  = 20;					// default number of bins for the H-test vs time plot		

	char FT1[MAX_FILELEN];
	char output[MAX_FILELEN];
	char output_pos_file[MAX_FILELEN];
	char phasecol[32];
	strcpy(phasecol,"PULSE_PHASE");
	char command[128];
	
    char error_buffer[128];
    char history[128];
	int  par_file      = 0;
	int  FT1_file	   = 0;
	
	int  ophase        = 0;
	int  graph	       = 0;
	int  output_file   = 0;
	int  output_pos    = 0;
	int  phase_replace = 0;

	double intpart;

	FILE *outputf;
	FILE *output_posf;
	FILE *temp_tim;

	/* ------------------------------------------------- //
	// fits file definitions
	// ------------------------------------------------- */

	fitsfile *ft1;
	fitsfile *ft2;

	char colname[FLEN_VALUE];
	int open_status = 0, status = 0, status2 = 0, hdupos, ncols_FT1,anynul;
	long nrows_FT1;
	double nulval = 0.;

	int FT1_time_col,FT1_phase_col;
	int nrows2, rows_status, rows_left;
	int max_rows = 10000;

	/* ------------------------------------------------- //
	// Time and satellite position definitions
	// ------------------------------------------------- */

	double temptime;
	double minFT1time = 999999999., maxFT1time = 0.;
	
	double time_MET_TT[max_rows], time_MJD_TT;
	double obs_earth[max_rows][3];

	double sctime1 = 0., sctime2 = 0.;
	double scposn1[3] = {0., 0., 0.};
	double scposn2[3] = {0., 0., 0.};
	double intposn[3] = {0., 0., 0.};

	double fract = 0.;
	double length1, length2, length12, intlength;
    double vector12[3], vectprod_out[3];
	double inttheta, factor_cos, factor_sin;

	longdouble lasttime, tzrmjd_bary;
	double lastpos[3];
	
	double tpb;

	/* ------------------------------------------------- //
	// Additional definitions
	// ------------------------------------------------- */

	psr[0].correctTroposphere = 0;

	/* ------------------------------------------------- //
	// Command-line arguments
	// ------------------------------------------------- */

	for (i=2;i<argc;i++)
	{
		if (strcmp(argv[i],"-f")==0)
		{
			par_file = 1;
			strcpy(parFile[0],argv[i+1]); 
		}
		else if (strcmp(argv[i],"-ft1")==0)
		{
			FT1_file = 1;
			strcpy(FT1,argv[i+1]);
		}
		else if (strcmp(argv[i],"-output")==0)
		{
			output_file = 1;
			strcpy(output,argv[i+1]);
		}
		else if (strcmp(argv[i],"-phase")==0)
		{
			phase_replace = 1;
		}
		else if (strcmp(argv[i],"-ophase")==0)
		{
			ophase = 1;
		}
		else if (strcmp(argv[i],"-col")==0)
		{
			strcpy(phasecol,argv[i+1]);
		}	
		else if (strcmp(argv[i],"-h")==0)
		{
			printf("\n TEMPO2 fermi plugin\n");
			printf("======================\n");
			printf("\n USAGE: \n\t tempo2 -gr photons -ft1 FT1.fits -f par.par\n");
			printf("\n Command line options:\n");
			printf("\t -output XXX: writes event times and phases in file XXX\n");
			printf("\t -phase: stores phases in the FT1 by the ones calculated by TEMPO2\n");
			printf("\t -ophase: will calculate orbital phases instead of pulse phases. Changes default column name to ORBITAL_PHASE.\n");
			printf("\t -col XXX: phases will be stored in column XXX. Default is PULSE_PHASE\n");
			printf("\t -h: this help.\n");
			printf("===============================================");
			printf("===============================================\n");
			exit(0);			
		}
    }
    
	if (!FT1_file)
	{
		printf("No input FT1 file !\n");
		printf("At least it should look like :\n");
		printf("\t tempo2 -gr photons -ft1 FT1.fits -f par.par\n");
		printf("More options are available. If you need help, please type:\n");
		printf("\t tempo2 -gr photons -h\n");
		exit(0);
	}
	if (!par_file)
	{
		printf("No ephemeris file !\n");
		printf("At least it should look like :\n");
		printf("\t tempo2 -gr photons -ft1 FT1.fits -f par.par\n");
		printf("More options are available. If you need help, please type:\n");
		printf("\t tempo2 -gr photons -h\n");
		exit(0);
	}
	if (output_file)
	{
		outputf = fopen(output,"w+");
	}
	if (output_pos)
	{
		output_posf = fopen(output_pos_file,"w+");
	}
	if (ophase)
	{
		if (strcmp(phasecol,"PULSE_PHASE") == 0) strcpy(phasecol,"ORBITAL_PHASE");
	}

	/* ------------------------------------------------- //
	// FT1 file
	// ------------------------------------------------- */

        // FIXME: make sure it's TIMESYS==TDB and ??==SOLARSYSTEM
  	if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
    {
		fits_movabs_hdu(ft1,2,NULL,&status);

		fits_get_num_rows(ft1, &nrows_FT1, &status);
		fits_get_num_cols(ft1, &ncols_FT1, &status);

		fits_get_colname(ft1,CASESEN,"TIME",colname,&FT1_time_col,&status);
		
		//
		
		for (i=1;i<=nrows_FT1;i++)
		{
			fits_read_col_dbl(ft1,FT1_time_col,i,1,1,nulval,&temptime,&anynul,&status);
			if (temptime < minFT1time) minFT1time = temptime;
			if (temptime > maxFT1time) maxFT1time = temptime;
		}
		
		//

		if (phase_replace)
		{
            fits_get_colname(ft1,CASESEN,phasecol,colname,&FT1_phase_col,&status);
            
            if (status != 0)
            {
                    fits_insert_col(ft1,ncols_FT1 + 1,phasecol,"1D", &status2);
                            
                    if (status2 != 0)
                    {
                            fits_get_errstatus( status2, error_buffer);
                            printf( "fits_insert_col: %s\n", error_buffer);
                            exit(-1);
                    }
            }
                        
			status = 0;
			fits_get_colname(ft1,CASESEN,phasecol,colname,&FT1_phase_col,&status);
    	}
	}
	
	if (open_status != 0)
	{
		printf("Can't find %s !\n",FT1);
		exit(-1);
	}

	fits_close_file(ft1, &status);

	printf("First photon date in FT1: %lf MET (s)\n",minFT1time);
	printf(" Last photon date in FT1: %lf MET (s)\n\n",maxFT1time);
	
	// The FT1 file is cut into small tim files with # rows <= 10000

	rows_status = 1;
	rows_left	= nrows_FT1 - rows_status + 1;
	nrows2 		= min(rows_left,max_rows);

	float *phase;
	float *times;
        
        phase  = (float *)calloc(max_rows,sizeof(float));
        times  = (float *)calloc(max_rows,sizeof(float));
	
	float tmin   = 100000., tmax   = -100000.;
	
	
	/* ------------------------------------------------- //
	// Barycentric TZRMJD
	// ------------------------------------------------- */
	
	// A temporary file is created. It is first used to get the barycentered TZR
	temp_tim = fopen(timFile[0],"w+");
	fprintf(temp_tim,"FORMAT 1\n");
	fprintf(temp_tim," photons 0.0 %.12lf 0.00000 BAT\n",1.);
	fclose(temp_tim);

	// Load the arrival times
	readParfile(psr,parFile,timFile,*npsr);
	readTimfile(psr,timFile,*npsr);
	
	if (ophase && (psr[0].param[param_pb].paramSet[0] == 0))
	{
		printf("Error: no binary parameters found in %s !\n",parFile[0]);
		sprintf(command,"rm -f %s",timFile[0]);
		system(command);
		exit(-1);
	}
	
	if (psr->tzrsite[0] == '@')
	{
		tzrmjd_bary = psr->param[param_tzrmjd].val[0];
	}
	else
	{
		psr->obsn[0].sat = psr->param[param_tzrmjd].val[0];
		psr->obsn[0].freq = psr->param[param_tzrfrq].val[0];
		strcpy(psr->obsn[0].telID, psr->tzrsite); 
		psr->obsn[0].deleted = 0;
		psr->obsn[0].nFlags = 0;
		psr->obsn[0].delayCorr = 1;
		psr->obsn[0].clockCorr = 1;
		
		// Form barycentric TZR
		formBatsAll(psr,*npsr);
		tzrmjd_bary = psr->obsn[0].bat;
	}
	
	
	

	/* ------------------------------------------------- //
	// Main loop
	// ------------------------------------------------- */	

	while (rows_left > 0)
	{
		printf("Treating events # %d to %d... \n",rows_status,rows_status + nrows2 - 1);

		// Acquisition of the non-barycentered TOAs in MET TT
		j = 0;

		if (!fits_open_file(&ft1,FT1, READONLY, &status))
		{
			fits_movabs_hdu(ft1,2,NULL,&status);

			for (i=rows_status;i<rows_status + nrows2;i++)
			{
				fits_read_col_dbl(ft1,FT1_time_col,i,1,1,nulval,&time_MET_TT[j],&anynul,&status);
				j++;
			}
		}

		fits_close_file(ft1, &status);


		// A temporary file is created. TOAs in MJD TT are stored is this file
		temp_tim = fopen(timFile[0],"w+");
		fprintf(temp_tim,"FORMAT 1\n");

		for (i=0;i<nrows2;i++)
		{
			time_MJD_TT = met2mjd(time_MET_TT[i]);		
			fprintf(temp_tim," photons 0.0 %.12lf 0.00000 BAT\n",time_MJD_TT);
		}

		fclose(temp_tim);

		// Load the arrival times
		readParfile(psr,parFile,timFile,*npsr); 
		readTimfile(psr,timFile,*npsr);
		
		


		/* ------------------------------------------------- //
		// Delays, corrections and positions
		// ------------------------------------------------- */
		for (i=0;i<nrows2;i++)
		{
			psr->obsn[i].deleted = 0;
			psr->obsn[i].nFlags = 0;
			psr->obsn[i].delayCorr = 0;	// Don't correct delays for event i
			psr->obsn[i].clockCorr = 0;	// Don't make clock correction TT -> TDB

			// Position replacements
			for (k=0;k<3;k++) 
			{
				psr[0].obsn[i].observatory_earth[k] = 0;
			}
		}


		/* ------------------------------------------------- //
		// Calculation of the event phases - step 1
                // Step 1 is for all but the last photon
		// ------------------------------------------------- */

		// keep track of the last TOA before shifting the others
		lasttime = psr->obsn[psr[0].nobs-1].sat;
		for (k=0;k<3;k++) lastpos[k] = psr->obsn[psr[0].nobs-1].observatory_earth[k];

		for (i=psr[0].nobs-1;i>0;i--)
		{
			psr->obsn[i].sat = psr->obsn[i-1].sat;

			for (k=0;k<3;k++) psr[0].obsn[i].observatory_earth[k] = psr[0].obsn[i-1].observatory_earth[k];
		}

	
		/* ------------------------------------------------- //
		// Stick in a fake obs to get the reference phase
		// ------------------------------------------------- */
		psr->obsn[0].sat = tzrmjd_bary;
		psr->obsn[0].freq = 0.;
		psr->obsn[0].deleted = 0;
		psr->obsn[0].nFlags = 0;
		psr->obsn[0].delayCorr = 0;
		psr->obsn[0].clockCorr = 0;
                printf("%20.15Lf\t%20.15Lf\n", psr->obsn[0].sat, psr->obsn[0].bat);
                printf("%20.15Lf\t%20.15Lf\n", psr->obsn[1].sat, psr->obsn[1].bat);
		
		// ------------------------------------------------- //
		// Form barycentric arrival times - step 1
		// ------------------------------------------------- //
		formBatsAll(psr,*npsr);
                printf("%20.15Lf\t%20.15Lf\n", psr->obsn[0].sat, psr->obsn[0].bat);
                printf("%20.15Lf\t%20.15Lf\n", psr->obsn[1].sat, psr->obsn[1].bat);

		// ------------------------------------------------- //
		// Calculate event phases - step 1
		// ------------------------------------------------- //		
		formResiduals(psr,*npsr,0.0);

		for (i=1;i<nrows2;i++)
		{
                        if (ophase == 0)
                        {
                                phase[event] = modf(psr[0].obsn[i].phase,&intpart);			
                                if (phase[event] < 0.) phase[event]++;
                        }
                        else
                        {
                                if(psr[0].param[param_tasc].paramSet[0])
                                {
                                        tpb = (psr[0].obsn[i].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
                                }
                                else if(psr[0].param[param_t0].paramSet[0])
                                {
                                        tpb = (psr[0].obsn[i].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]);
                                }
                                
                                phase[event] = modf(tpb+1000000.0,&intpart);;
                                if (phase[event] < 0.) phase[event]++;
                        }

                        times[event] = psr[0].obsn[i].bat;
                        if (times[event] > tmax) 	tmax = times[event];
                        if (times[event] < tmin)	tmin = times[event];
				
			if (output_file)
			{
				fprintf(outputf,"%d\t",event + rows_status);
				fprintf(outputf,"%20.15Lf %12.10le\n",psr[0].obsn[i].bat,phase[event]);
			}
	
			event++;
		}		

		// ------------------------------------------------- //
		// Calculation of the event phases - step 2
                // Step 2 repeats the above code for the last photon
		// ------------------------------------------------- //
		psr[0].obsn[1].sat = lasttime;
		for (k=0;k<3;k++) psr[0].obsn[1].observatory_earth[k] = lastpos[k];
	
		// ------------------------------------------------- //
		// Form barycentric arrival times - step 2
		// ------------------------------------------------- //
		formBatsAll_fermi(psr,*npsr);

		// ------------------------------------------------- //
		// Calculate event phases - step 2
		// ------------------------------------------------- //		
		formResiduals(psr,*npsr,0.0);
	
                if (ophase == 0)
                {
                        phase[event] = modf(psr[0].obsn[1].phase,&intpart);
                        if (phase[event] < 0.) phase[event]++;
                }
                else
                {
                        if(psr[0].param[param_tasc].paramSet[0])
                        {
                                tpb = (psr[0].obsn[1].bat-psr[0].param[param_tasc].val[0])/(psr[0].param[param_pb].val[0]);
                        }
                        else if(psr[0].param[param_t0].paramSet[0])
                        {
                                tpb = (psr[0].obsn[1].bat-psr[0].param[param_t0].val[0])/(psr[0].param[param_pb].val[0]);
                        }
                        
                        phase[event] = modf(tpb+1000000.0,&intpart);;
                        if (phase[event] < 0.) phase[event]++;
                }
        
                times[event] = psr[0].obsn[1].bat;		
                if (times[event] > tmax) 	tmax = times[event];
                if (times[event] < tmin)	tmin = times[event];

		if (output_file)
		{
			fprintf(outputf,"%d\t",event + rows_status);
			fprintf(outputf,"%20.15Lf %12.10le\n",psr[0].obsn[1].bat,phase[event]);
		}

		// ------------------------------------------------- //
		// Put new phases into the input FT1 file
		// ------------------------------------------------- //
		if (phase_replace)
		{
  			if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
    			{
				fits_movabs_hdu(ft1,2,NULL,&status);
				
				for (event2=rows_status;event2<rows_status+nrows2;event2++)
				{
					if (graph == 0)	i = event2 - rows_status;
					else i = event2 - 1;
				
					fits_write_col_flt(ft1,FT1_phase_col,event2,1,1,&phase[i],&status);
				}
			}

			fits_close_file(ft1, &status);		
		}

		/* ------------------------------------------------- //
		// End of the loop
		// ------------------------------------------------- */	
		rows_status	+= nrows2;
		rows_left	= nrows_FT1 - rows_status + 1;
		nrows2 		= min(rows_left,max_rows);
		
                event = 0;
	}
        
    // ------------------------------------------------- //
	// Add a little bit of history to the header
	// ------------------------------------------------- //
	if (phase_replace)
	{
	        sprintf(history,"Pulse phases calculated with the TEMPO2 photons plugin using ephemeris %s",parFile[0]);
	        
	        if (!fits_open_file(&ft1,FT1, READWRITE, &open_status))
	        {
	                status = 0;
	                fits_movabs_hdu(ft1,2,NULL,&status);
	                fits_write_history(ft1,history,&status);
	                
	                if (status != 0)
	                {
	                        fits_get_errstatus(status,error_buffer);
	                        printf( "fits_insert_col: %s\n", error_buffer);
	                        exit(-1);
	                }
	        }
	        
	        fits_close_file(ft1, &status);
	}

	if (output_file) fclose(outputf);
	if (output_pos)  fclose(output_posf);

	printf("Done with %s\n",psr[0].name);

	// ------------------------------------------------- //
	// Clean temporary tim file
	// ------------------------------------------------- //	

	sprintf(command,"rm -f %s",timFile[0]);
	system(command);
	
	return 0;
}





/* ------------------------------------------------- //
// Barycentering tools
// ------------------------------------------------- */

void clock_corrections_fermi(pulsar *psr,int npsr)
{
	// toa2utc(psr,npsr);                      // UTC(Observatory) -> UTC(NIST) 
	// tai2ut1(psr,npsr);                      // TAI -> UT1			 
	tt2tb(psr,npsr);                        // Rough estimate of TT-TB (+-2.2 microsec) 
}

void ephemeris_routines_fermi(pulsar *psr,int npsr)
{
	vectorPulsar(psr,npsr);                 // Form a vector pointing at the pulsar 
	readEphemeris(psr,npsr,0);              // Read the ephemeris 
	// get_obsCoord(psr,npsr);                 // Get Coordinate of observatory relative to Earth's centre 
	tt2tb(psr,npsr);                        // Observatory/time-dependent part of TT-TB 
	readEphemeris(psr,npsr,0);              // Re-evaluate ephemeris with correct TB 
}

void extra_delays_fermi(pulsar *psr,int npsr)
{
	calculate_bclt(psr,npsr);               // Calculate bclt
}

void formBatsAll_fermi(pulsar *psr,int npsr)
{
	clock_corrections_fermi(psr,npsr);            	// Clock corrections  ... 
	
	//	printf("Reading ephemeris routines...\n");
	ephemeris_routines_fermi(psr,npsr);           	// Ephemeris routines ... 
	
	//	printf("Reading extra delays...\n");
	extra_delays_fermi(psr,npsr);                 	// Other time delays  ... 
	
	//	printf("Forming barycentric arrival times...\n");
	formBats(psr,npsr);                     	// Form Barycentric arrival times 
	
	//	printf("Evaluating secular motion...\n");
	secularMotion(psr,npsr);
}

//


static char random_letter(int is_cap)
{
   int letter = (int)(26 * (rand() / (RAND_MAX + 1.0)));
   return((char)((is_cap == 1) ? (letter + 65) : (letter + 97)));
}

static char random_number()
{
   int number = (int)(10 * (rand() / (RAND_MAX + 1.0)));
   return((char)(number + 48));
}

static void random_string(int length, char *str)
{
	int i;
	int char_type;
	
	for(i = 0; i < length; i++)
	{
		char_type = (int)(3 * (rand() / (RAND_MAX + 1.0)));
	
		switch(char_type)
		{
			case 0:
				str[i] = random_letter(0);
				break;
			case 1:
				str[i] = random_letter(1);
				break;
			case 2:
				str[i] = random_number();
				break;
			default:
				str[i] = random_number();
				break;
		}
	}  
}


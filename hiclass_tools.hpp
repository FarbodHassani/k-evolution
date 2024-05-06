//////////////////////////
// class_tools.hpp
//////////////////////////
//
// interface to linear Boltzmann code CLASS
//
// Author (k-evolution): Farbod Hassani (Université de Genève & Universitetet i Oslo)
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: Jan 2023
//
//////////////////////////

#ifndef CLASS_TOOLS_HEADER
#define CLASS_TOOLS_HEADER

#if defined(HAVE_HICLASS)

#include <gsl/gsl_spline.h>
#include "parser.hpp"

using namespace std;
using namespace LATfield2;

//////////////////////////
// initializeCLASSstructures
//////////////////////////
// Description:
//   initializes CLASS structures containing interpolation tables for various transfer functions
//
// Arguments:
//   sim               simulation metadata structure
//   ic                settings for IC generation
//   cosmo             cosmological parameter structure
//   class_background  CLASS structure that will contain the background
//   class_perturbs    CLASS structure that will contain perturbations
//   class_spectra     CLASS structure that will contain spectra
//   output_value      CLASS parameter value specifying the output (optional)
//
// Returns:
//
//////////////////////////

void initializeCLASSstructures(metadata & sim, icsettings & ic, cosmology & cosmo, background & class_background,  thermo & class_thermo, perturbs & class_perturbs, parameter * params = NULL, int numparam = 0, const char * output_value = "dTk, vTk, mPk")
{
	precision class_precision;
	transfers class_transfers;
  primordial class_primordial;
	nonlinear class_nonlinear;
  spectra class_spectra;
  lensing class_lensing;
  output class_output;
	file_content class_filecontent;
	ErrorMsg class_errmsg;
	char filename[] = "initializeCLASSstructures";
	char tmp[8 * PARAM_MAX_LENGTH];
	double recfast_z_initial;
	double perturb_sampling_stepsize;
	int recfast_Nz0;
	int i;
  int num_entries;
  #ifdef HAVE_HICLASS_BG
  if (cosmo.gravity_model==0)
    num_entries = 29; // If we have EFT theory - propto_omega we have to read 30 class params
  else if (cosmo.gravity_model==1)
  {
    num_entries = 27;
  }
  else
  {
    if(parallel.isRoot())  cout << " \033[1;31m ERROR: The k-essence gravity is not specified correctly!  \033[0m"<< endl;
    parallel.abortForce();
  }
  #else
    num_entries = 29;
  #endif
#ifdef CLASS_K_PER_DECADE_FOR_PK
	int k_per_decade_for_pk;
	if (numparam == 0 || !parseParameter(params, numparam, "k_per_decade_for_pk", k_per_decade_for_pk))
	{
		num_entries += 1;
		k_per_decade_for_pk = CLASS_K_PER_DECADE_FOR_PK;
	}
#endif
	if (cosmo.num_ncdm > 0) num_entries += 3;
	if (parallel.isRoot()) num_entries += 8;

	if ((1.015 * ic.z_ic + 0.01) > 9999.)
	{
		if (numparam == 0 || !parseParameter(params, numparam, "recfast_z_initial", recfast_z_initial))
		{
			num_entries += 1;
			recfast_z_initial = ic.z_ic * 1.02;
		}
		if (numparam == 0 || !parseParameter(params, numparam, "recfast_Nz0", recfast_Nz0))
		{
			num_entries += 1;
			recfast_Nz0 = 2 * (int) ceil(ic.z_ic * 1.02);
		}
	}
	if (numparam == 0 || !parseParameter(params, numparam, "perturb_sampling_stepsize", perturb_sampling_stepsize))
	{
		num_entries += 1;
		perturb_sampling_stepsize = 0.01;
	}

	num_entries += numparam;

	parser_init(&class_filecontent, num_entries, filename, class_errmsg);

	for (i = 0; i < num_entries; i++)
		class_filecontent.read[i] = _FALSE_;

	i = 0;

	sprintf(class_filecontent.name[i], "root");
	sprintf(class_filecontent.value[i++], "%s%s_class", sim.output_path, sim.basename_generic);

	sprintf(class_filecontent.name[i], "k_pivot");
	sprintf(class_filecontent.value[i++], "%e", ic.k_pivot);

	sprintf(class_filecontent.name[i], "A_s");
	sprintf(class_filecontent.value[i++], "%e", ic.A_s);

	sprintf(class_filecontent.name[i], "n_s");
	sprintf(class_filecontent.value[i++], "%f", ic.n_s);

	sprintf(class_filecontent.name[i], "z_pk");
	if (ic.z_ic > sim.z_in)
		sprintf(class_filecontent.value[i++], "%f, %f, 0", 1.015 * ic.z_ic + 0.01, sim.z_in);
	else
		sprintf(class_filecontent.value[i++], "%f, 0", 1.015 * ic.z_ic + 0.01);

	sprintf(class_filecontent.name[i], "output");
	sprintf(class_filecontent.value[i++], "%s", output_value);

  sprintf(class_filecontent.name[i], "gauge");
  sprintf(class_filecontent.value[i++], "synchronous");

  sprintf(class_filecontent.name[i], "extra metric transfer functions");
  sprintf(class_filecontent.value[i++], "y");


	sprintf(class_filecontent.name[i], "P_k_ini type");
	sprintf(class_filecontent.value[i++], "analytic_Pk");

    sprintf(class_filecontent.name[i], "P_k_max_h/Mpc");
	sprintf(class_filecontent.value[i++], "%f", 2. * M_PI * (double) sim.numpts / sim.boxsize);

	sprintf(class_filecontent.name[i], "h");
	sprintf(class_filecontent.value[i++], "%f", cosmo.h);

	sprintf(class_filecontent.name[i], "Omega_cdm");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_cdm);

	sprintf(class_filecontent.name[i], "Omega_b");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_b);

	sprintf(class_filecontent.name[i], "Omega_g");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_g);

	sprintf(class_filecontent.name[i], "Omega_ur");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_ur);

	sprintf(class_filecontent.name[i], "Omega_fld");
	sprintf(class_filecontent.value[i++], "%e", 0.);
  //
	// sprintf(class_filecontent.name[i], "w0_fld");
	// sprintf(class_filecontent.value[i++], "%g", cosmo.w0_fld);
  //
	// sprintf(class_filecontent.name[i], "wa_fld");
	// sprintf(class_filecontent.value[i++], "%g", cosmo.wa_fld);
  //
	// sprintf(class_filecontent.name[i], "cs2_fld");
	// sprintf(class_filecontent.value[i++], "%g", cosmo.cs2_fld);

	sprintf(class_filecontent.name[i], "N_ncdm");
	sprintf(class_filecontent.value[i++], "%d", cosmo.num_ncdm);

  sprintf(class_filecontent.name[i], "write parameters");
  sprintf(class_filecontent.value[i++], "yeap");

  sprintf(class_filecontent.name[i], "write background");
  sprintf(class_filecontent.value[i++], "yeap");

  sprintf(class_filecontent.name[i], "format");
  sprintf(class_filecontent.value[i++], "class");

  sprintf(class_filecontent.name[i], "headers");
  sprintf(class_filecontent.value[i++], "yeap");

	sprintf(class_filecontent.name[i], "perturb_sampling_stepsize");
	sprintf(class_filecontent.value[i++], "%g", perturb_sampling_stepsize);

  sprintf(class_filecontent.name[i], "output_background_smg");
  sprintf(class_filecontent.value[i++],"%d", 4);

  sprintf(class_filecontent.name[i], "skip_stability_tests_smg");
  sprintf(class_filecontent.value[i++],"yes");

  sprintf(class_filecontent.name[i], "Omega_Lambda");
  sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_Lambda);

  sprintf(class_filecontent.name[i], "Omega_smg");
  sprintf(class_filecontent.value[i++], "%d", -1);

  #ifdef HAVE_HICLASS_BG
  if (cosmo.gravity_model==0)
  {
    sprintf(class_filecontent.name[i], "gravity_model");
    sprintf(class_filecontent.value[i++], "propto_omega");

    sprintf(class_filecontent.name[i], "expansion_model");
    sprintf(class_filecontent.value[i++], "wowa");

    sprintf(class_filecontent.name[i], "expansion_smg");
    sprintf(class_filecontent.value[i++],"%e, %e, %e",cosmo.Omega_kessence, cosmo.w_kessence, cosmo.w_a_kessence);

    sprintf(class_filecontent.name[i], "parameters_smg");
    sprintf(class_filecontent.value[i++],"%e, %e, %e, %e, %e", 3.0 * (1+cosmo.w_kessence)/cosmo.cs2_kessence, 0.0, 0.0, 0.0, 1.0);
  }
  else if (cosmo.gravity_model==1)
  {
    sprintf(class_filecontent.name[i], "gravity_model");
    sprintf(class_filecontent.value[i++], "k_essence_power");

    sprintf(class_filecontent.name[i], "parameters_smg");
    sprintf(class_filecontent.value[i++],"%e, %e, %e, %e, %e, %e",cosmo.Xt, cosmo.g0, cosmo.g2, cosmo.g4, cosmo.phi_i, cosmo.X_i);
  }
  #else
  sprintf(class_filecontent.name[i], "gravity_model");
  sprintf(class_filecontent.value[i++], "propto_omega");

  sprintf(class_filecontent.name[i], "expansion_model");
  sprintf(class_filecontent.value[i++], "wowa");

  sprintf(class_filecontent.name[i], "expansion_smg");
  sprintf(class_filecontent.value[i++],"%e, %e, %e",cosmo.Omega_kessence, cosmo.w_kessence, 0.0);

  sprintf(class_filecontent.name[i], "parameters_smg");
  sprintf(class_filecontent.value[i++],"%e, %e, %e, %e, %e", 3.0 * (1+cosmo.w_kessence)/cosmo.cs2_kessence, 0.0, 0.0, 0.0, 1.0);
  #endif
  //

  // sprintf(class_filecontent.name[i], "tuning_dxdy_guess_smg");
  // sprintf(class_filecontent.value[i++], "%e", 1.e-7);
  //
  // sprintf(class_filecontent.name[i], "tuning_index_smg");
  // sprintf(class_filecontent.value[i++], "%f", 1.);

#ifdef CLASS_K_PER_DECADE_FOR_PK
	sprintf(class_filecontent.name[i], "k_per_decade_for_pk");
	sprintf(class_filecontent.value[i++], "%d", k_per_decade_for_pk);
#endif

	if ((1.015 * ic.z_ic + 0.01) > 9999.)
	{
		sprintf(class_filecontent.name[i], "recfast_z_initial");
		sprintf(class_filecontent.value[i++], "%e", recfast_z_initial);

		sprintf(class_filecontent.name[i], "recfast_Nz0");
		sprintf(class_filecontent.value[i++], "%d", recfast_Nz0);
	}

	if (parallel.isRoot())
	{
		sprintf(class_filecontent.name[i], "background_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "perturbations_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "spectra_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "thermodynamics_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "transfer_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "primordial_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "nonlinear_verbose");
		sprintf(class_filecontent.value[i++], "1");

    sprintf(class_filecontent.name[i], "input_verbose");
    sprintf(class_filecontent.value[i++], "1");
	}

	while (numparam > 0)
	{
		numparam--;
		if (!params[numparam].used)
		{
			sprintf(class_filecontent.name[i], "%s", params[numparam].name);
			sprintf(class_filecontent.value[i++], "%s", params[numparam].value);
		}
	}
	if (cosmo.num_ncdm > 0)
	{
		sprintf(class_filecontent.name[num_entries-3], "m_ncdm");
		sprintf(class_filecontent.value[num_entries-3], "%f", cosmo.m_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-3], cosmo.m_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-3], "%s", tmp);
		}

		sprintf(class_filecontent.name[num_entries-2], "T_ncdm");
		sprintf(class_filecontent.value[num_entries-2], "%f", cosmo.T_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-2], cosmo.T_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-2], "%s", tmp);
		}

		sprintf(class_filecontent.name[num_entries-1], "deg_ncdm");
		sprintf(class_filecontent.value[num_entries-1], "%f", cosmo.deg_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-1], cosmo.deg_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-1], "%s", tmp);
		}
	}

	COUT << " gevolution is calling hiclass..." << endl << endl;

  if (input_init(&class_filecontent, &class_precision, &class_background, &class_thermo, &class_perturbs, &class_transfers, &class_primordial, &class_spectra, &class_nonlinear, &class_lensing, &class_output, class_errmsg) == _FAILURE_)
	{
		COUT << " error: calling input_init from hiclass library failed!" << endl << " following error message was passed: " << class_errmsg << endl;
		parallel.abortForce();
	}

	parser_free(&class_filecontent);

	if (background_init(&class_precision, &class_background) == _FAILURE_)
	{
		COUT << " error: calling background_init from hiclass library failed!" << endl << " following error message was passed: " << class_background.error_message << endl;
		parallel.abortForce();
	}

	if (thermodynamics_init(&class_precision, &class_background, &class_thermo) == _FAILURE_)
	{
		COUT << " error: calling thermodynamics_init from hiclass library failed!" << endl << " following error message was passed: " << class_thermo.error_message << endl;
		parallel.abortForce();
	}

	if (perturb_init(&class_precision, &class_background, &class_thermo, &class_perturbs) == _FAILURE_)
	{
		COUT << " error: calling perturb_init from hiclass library failed!" << endl << " following error message was passed: " << class_perturbs.error_message << endl;
		parallel.abortForce();
	}

	if (primordial_init(&class_precision, &class_perturbs, &class_primordial) == _FAILURE_)
	{
		COUT << " error: calling primordial_init from hiclass library failed!" << endl << " following error message was passed: " << class_primordial.error_message << endl;
		parallel.abortForce();
	}

	if (nonlinear_init(&class_precision, &class_background, &class_thermo, &class_perturbs, &class_primordial, &class_nonlinear) == _FAILURE_)
	{
		COUT << " error: calling nonlinear_init from hiclass library failed!" << endl << " following error message was passed: " << class_nonlinear.error_message << endl;
		parallel.abortForce();
	}

	if (transfer_init(&class_precision, &class_background, &class_thermo, &class_perturbs, &class_nonlinear, &class_transfers) == _FAILURE_)
	{
		COUT << " error: calling transfer_init from hiclass library failed!" << endl << " following error message was passed: " << class_transfers.error_message << endl;
		parallel.abortForce();
	}

  if (spectra_init(&class_precision, &class_background, &class_perturbs, &class_primordial, &class_nonlinear, &class_transfers, &class_spectra) == _FAILURE_)
	{
		COUT << " error: calling spectra_init from hiclass library failed!" << endl << " following error message was passed: " << class_spectra.error_message << endl;
		parallel.abortForce();
	}

	// Now, free unneeded structures

	if (lensing_free(&class_lensing) == _FAILURE_)
	{
		COUT << " error: calling lensing_free from hiclass library failed!" << endl << " following error message was passed: " << class_lensing.error_message << endl;
		parallel.abortForce();
	}

	if (transfer_free(&class_transfers) == _FAILURE_)
	{
		COUT << " error: calling transfer_free from hiclass library failed!" << endl << " following error message was passed: " << class_transfers.error_message << endl;
		parallel.abortForce();
	}

	if (nonlinear_free(&class_nonlinear) == _FAILURE_)
	{
		COUT << " error: calling nonlinear_free from hiclass library failed!" << endl << " following error message was passed: " << class_nonlinear.error_message << endl;
		parallel.abortForce();
	}

	if (primordial_free(&class_primordial) == _FAILURE_)
	{
		COUT << " error: calling primordial_free from hiclass library failed!" << endl << " following error message was passed: " << class_primordial.error_message << endl;
		parallel.abortForce();
	}

	COUT << endl << " hiclass structures initialized successfully." << endl;
}


//////////////////////////
// freeCLASSstructures
//////////////////////////
// Description:
//   frees CLASS structures containing interpolation tables for various transfer functions
//
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_thremo      CLASS structure that contains thermodynamics
//   class_perturbs    CLASS structure that contains perturbations
//
// Returns:
//
//////////////////////////

void freeCLASSstructures(background & class_background, thermo & class_thermo, perturbs & class_perturbs)
{
	if (perturb_free(&class_perturbs) == _FAILURE_)
	{
		COUT << " error: calling perturb_free from CLASS library failed!" << endl << " following error message was passed: " << class_perturbs.error_message << endl;
		parallel.abortForce();
	}

	if (thermodynamics_free(&class_thermo) == _FAILURE_)
	{
		COUT << " error: calling thermodynamics_free from CLASS library failed!" << endl << " following error message was passed: " << class_thermo.error_message << endl;
		parallel.abortForce();
	}

	if (background_free(&class_background) == _FAILURE_)
	{
		COUT << " error: calling background_free from CLASS library failed!" << endl << " following error message was passed: " << class_background.error_message << endl;
		parallel.abortForce();
	}

}


//////////////////////////
// loadTransferFunctions
//////////////////////////
// Description:
//   loads a set of tabulated transfer functions from some precomputed CLASS structures
//
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   tk_delta          will point to the gsl_spline which holds the tabulated
//                     transfer function for delta (memory will be allocated)
//   tk_theta          will point to the gsl_spline which holds the tabulated
//                     transfer function for theta (memory will be allocated)
//   qname             string containing the name of the component (e.g. "cdm"); if no string is
//                     specified, the transfer functions for phi and psi are returned instead!
//   boxsize           comoving box size (in the same units as used in the CLASS output)
//   z                 redshift at which the transfer functions are to be obtained
//   h                 conversion factor between 1/Mpc and h/Mpc (theta is in units of 1/Mpc)
//
// Returns:
//
//////////////////////////
void loadTransferFunctions(background & class_background, perturbs & class_perturbs, gsl_spline * & tk_delta, gsl_spline * & tk_theta, const char * qname, const double boxsize, const double z, double h,  double Hconf_class, double Omega_m, double Omega_rad, double Omega_mg, double w_mg)
{
  int cols = 0, dcol = -1, tcol = -1, kcol = -1, phicol = -1, psicol= -1, etacol= -1, h_primecol= -1, eta_primecol= -1;
  double alpha, alpha_prime;
	double * k;
	double * tk_d;
	double * tk_t;
	double * data;
	char coltitles[_MAXTITLESTRINGLENGTH_] = {0};
	char dname[32];
	char tname[32];
	char kname[16];
	char * ptr;
  // hiclass bg inputs
  double a = 1./(1.+z);
  perturb_output_titles(&class_background, &class_perturbs, class_format, coltitles);
  if (strncmp(qname,"vx",strlen("vx")) == 0)
	{
		sprintf(dname, "vx_smg");
		sprintf(tname, "vx_prime_smg");
  }
	else if (qname != NULL)
	{
		sprintf(dname, "d_%s", qname);
		sprintf(tname, "t_%s", qname);
  }
  else if (qname == NULL)
	{
    COUT << " error in loadTransferFunctions (HAVE_CLASS)! why phi/psi columns are requested?! The radiation perturbation part should be implemented and tested!" << endl;
		parallel.abortForce();
	}
	sprintf(kname, "k (h/Mpc)");

	ptr = strtok(coltitles, _DELIMITER_);
	while (ptr != NULL)
	{
    if (strncmp(ptr, dname, strlen(dname)) == 0) dcol = cols;
    else if (strncmp(ptr, tname, strlen(tname)) == 0) tcol = cols;
    else if (strncmp(ptr, kname, strlen(kname)) == 0) kcol = cols;
    // gauge transformation
    else if (strncmp(ptr, "phi", strlen("phi")) == 0) phicol = cols;
    else if (strncmp(ptr, "psi", strlen("psi")) == 0) psicol = cols;
    else if (strncmp(ptr, "eta_prime", strlen("eta_prime")) == 0) eta_primecol = cols;
    else if (strncmp(ptr, "eta", strlen("eta")) == 0) etacol = cols;
    else if (strncmp(ptr, "h_prime", strlen("h_prime")) == 0) h_primecol = cols;
    cols++;
      ptr = strtok(NULL, _DELIMITER_);
  }

	if (dcol < 0 || (tcol < 0 && strncmp(qname,"cdm",strlen("cdm")) != 0 ) || kcol < 0 || (qname != NULL && (phicol < 0 || psicol < 0 || etacol < 0 || h_primecol < 0 || eta_primecol < 0) ) )
	{
		COUT << " error in loadTransferFunctions (HAVE_CLASS)! Unable to identify requested columns!" << endl;
		parallel.abortForce();
	}

	data = (double *) malloc(sizeof(double) * cols*class_perturbs.k_size[class_perturbs.index_md_scalars]);
	k = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_d = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_t = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);

	perturb_output_data(&class_background, &class_perturbs, class_format, z, cols, data);

	for (int i = 0; i < class_perturbs.k_size[class_perturbs.index_md_scalars]; i++)
	{
    k[i] = data[i*cols + kcol] * boxsize;
    alpha = (data[i*cols + h_primecol] + 6.0*data[i*cols + eta_primecol])/(2.0*data[i*cols + kcol]*data[i*cols + kcol]*h*h);// Note that in hiclass k and time derivative do not have the same unit, there is an h difference which is corrected. prime is [1/Mpc] while k[h/Mpc] and k*h give in 1/Mpc
    alpha_prime = data[i*cols + psicol] + data[i*cols + phicol] - data[i*cols + etacol];
    if (strncmp(qname,"vx",strlen("vx")) == 0)
     {
      tk_d[i] = -data[i*cols + dcol] + alpha; // gauge correction NOTE that v_x is in [1/Mpc] and pi(Newt)[hiclass unit] = V_x +alpha where pi(sync) =  V_x
      tk_t[i] = -data[i*cols + tcol] + alpha_prime;// gauge correction where pi'(Newt)[hiclass unit==1] = V'_x == pi'
      tk_t[i] += tk_d[i] * Hconf_class  - data[i*cols + psicol];  // zeta  = pi'(conformal_Newtonian) + H(conf)*pi - psi
     }

    else if (strncmp(qname,"cdm",strlen("cdm")) == 0) // tk_t cdm is 0 in sync gauge!
    {
      tk_d[i] = data[i*cols + dcol] - alpha * 3. * Hconf_class; // gauge correction; [alpha]=Mpc and H[1/Mpc]!
      tk_t[i] = (alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;// gauge correction + converting in gevolution!
      // boxsize/h --> unit in Mpc
      // k *h --> 1/Mpc
      // alpha is in 1/Mpc
      // unit as tk_t in hiclass is 1/Mpc and alpha = (h_prime[1/Mpc] + 6.*eta_prime[1/Mpc])/(2.*k[h/Mpc]*h*k[h/Mpc]*h); Results in  [alpha] = Mpc and [alpha k*h * k *h] = Mpc * 1/Mpc^2 = 1/Mpc
      // Also note that tk_t[i] is in 1/Mpc in class
      // boxsize[Mpc/h] / h give boxsize in Mpc. so at the end we have tk_t in gevolution without unit.
    }
    else if (strncmp(qname,"b",strlen("b")) == 0)
    {
      tk_d[i] = data[i*cols + dcol] - alpha * 3. * Hconf_class;
      tk_t[i] = (data[i*cols + tcol] + alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;
    }
    // else if (strncmp(qname,"fld",strlen("fld")) != 0)
    // {
    //   tk_d[i] = data[i*cols + dcol] - alpha * 3. * (1+ w_fld) * Hconf_class;
    //   tk_t[i] = (data[i*cols + tcol] + alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;
    // }
    else if (strncmp(qname,"g",strlen("g")) == 0)
    {
      tk_d[i] = data[i*cols + dcol] - alpha * 4. * Hconf_class;
      tk_t[i] = (data[i*cols + tcol] + alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;
    }
    else if (strncmp(qname,"ur",strlen("ur")) == 0)
    {
      tk_d[i] = data[i*cols + dcol] - alpha * 4. * Hconf_class;
      tk_t[i] = (data[i*cols + tcol] + alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;
    }
    else if (strncmp(qname,"tot",strlen("tot")) == 0)
    {
      tk_d[i] = data[i*cols + dcol] - alpha * 3 * (Omega_mg * (1. + w_mg) + Omega_m + Omega_rad * (1.+1./3.)) * Hconf_class;
      // Omega_fld*(1. + w_fld)  fld part is not included!
      tk_t[i] = (data[i*cols + tcol] + alpha * data[i*cols + kcol] * h * data[i*cols + kcol] * h) * boxsize/h;
    }
		if (i > 0)
		{
			if (k[i] < k[i-1])
			{
				COUT << " error in loadTransferFunctions (HAVE_CLASS)! k-values are not strictly ordered." << endl;
				parallel.abortForce();
			}
		}
	}

	free(data);

	tk_delta = gsl_spline_alloc(gsl_interp_cspline, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_theta = gsl_spline_alloc(gsl_interp_cspline, class_perturbs.k_size[class_perturbs.index_md_scalars]);

	gsl_spline_init(tk_delta, k, tk_d, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	gsl_spline_init(tk_theta, k, tk_t, class_perturbs.k_size[class_perturbs.index_md_scalars]);

	free(k);
	free(tk_d);
	free(tk_t);
}

#ifdef HAVE_HICLASS_BG
//////////////////////////
// Load background functions
//////////////////////////
// Description:
//   loads a set of tabulated background parameters from some precomputed CLASS structures
//
// Arguments:
//   class_background  CLASS structure that contains the background
//   bg_data           .......
//   qname             string containing the name of the component (e.g. "H [1/Mpc]"); if no string is
//   boxsize           comoving box size (in the same units as used in the CLASS output)
//   z                 redshift at which the transfer functions are to be obtained
//   h                 conversion factor between 1/Mpc and h/Mpc (theta is in units of 1/Mpc)
//
// Returns:
//
//////////////////////////

void loadBGFunctions(background & class_background, gsl_spline * & bg_data, const char * qname, double z_in)
{
	int cols = 0, bgcol = -1, bgcol2 = -1, zcol = -1;
	char coltitles[_MAXTITLESTRINGLENGTH_] = {0};
	char zname[16];
	double * a;
	double * bg;
	double * data;
	char * ptr;
  int num_points=30;
	int bg_size=0;
  double dz=0.005/num_points;
  double z1, z2;
	// Get the names of the background variables as a single string
	background_output_titles(&class_background, coltitles);

	// Get the redshift variable
	sprintf(zname, "z");
	ptr = strtok(coltitles, _DELIMITER_);
	while (ptr != NULL)
	{
  		if ((strncmp(ptr, qname, strlen(qname)) == 0) && (strlen(qname)== strlen(ptr))) bgcol = cols;
  		if (strncmp(ptr, zname, strlen(zname)) == 0) zcol = cols;
      cols++;
      ptr = strtok(NULL, _DELIMITER_);
	}

	if (bgcol < 0 || zcol < 0)
	{
		COUT << " error in loadBGFunctions (HAVE_CLASS)! Unable to identify requested columns!" << endl;
		parallel.abortForce();
	}

	data = (double *) malloc(sizeof(double) * cols * class_background.bt_size);
	if(!data)
  {
    COUT << " error in loadBGFunctions (HAVE_HICLASS_BG)! Unable to allocate memory!" << endl;
    parallel.abortForce();
  }

	background_output_data(&class_background, cols, data);

  a = (double *) malloc(sizeof(double) * (class_background.bt_size-bg_size+num_points));
  bg = (double *) malloc(sizeof(double) * (class_background.bt_size-bg_size+num_points));
	if(!a || !bg)
  {
    COUT << " error in loadBGFunctions (HAVE_HICLASS_BG)! Unable to allocate memory!" << endl;
    parallel.abortForce();
  }

  for (int i = bg_size; i < class_background.bt_size; i++)
	{
		a[i-bg_size] = 1. / (1. + data[i*cols + zcol]);
		bg[i-bg_size] = data[i*cols + bgcol];
		if (i > bg_size)
		{
			if (a[i-bg_size] < a[i-bg_size-1])
			{
				COUT << " error in loadBGFunctions (HAVE_CLASS)! a-values are not strictly ordered." << endl;
				parallel.abortForce();
			}
		}
	}
  for (int i=0;i<num_points;i++)
  {
    z1 = 1./a[class_background.bt_size-bg_size+i-1] -1.;
    z2 = 1./a[class_background.bt_size-bg_size+i-2] -1.;
    a[class_background.bt_size-bg_size+i] =  a[class_background.bt_size-bg_size+i-1] + 1./(1.+dz);
    bg[class_background.bt_size-bg_size+i] = bg[class_background.bt_size-bg_size+i-1] - dz *(bg[class_background.bt_size-bg_size+i-1] -  bg[class_background.bt_size-bg_size+i-2])/(z1 - z2);
  }
	free(data);

  bg_data = gsl_spline_alloc(gsl_interp_cspline, class_background.bt_size-bg_size+num_points);

	gsl_spline_init(bg_data, a, bg, class_background.bt_size-bg_size+num_points);

	free(a);
	free(bg);
}
#endif

#endif

#endif

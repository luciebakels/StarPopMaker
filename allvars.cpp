#include "allvars.h"

extern int ThisTask	{0};
extern int NTask	{0};

namespace constants {
	extern const double G_cm3_gi_si2			{ 6.672599e-8 };
	extern const double G_kpc_km2_Msi_si2		{ 4.3e-6 };
	extern const double rhocrit_Ms_kpci3		{ 2.7755e2 };
	extern const double rhocrit_g_cmi3			{ 1.9790e-29 };
	extern const double kpc_to_cm				{ 3.085678e21 };
	extern const double s_to_yr					{ 2.893777e-8 };
	extern const double Msun_kg					{ 1.9891e30 };
	extern const double Msun_g 					{ 1.9891e33 };
	extern const double Lsun_erg_s				{ 3.846e33 };
	extern const double Rsun_cm					{ 6.9551e10 };
	extern const double sigmaSB_erg_cm2_si_K4	{ 5.6704e-5 };
	extern const double atomicMass_kg			{ 1.660539e-27 };
	extern const double h_eV_s					{ 4.135668e-15 };
	extern const double h_erg_s					{ 1.054573e-27 };
	extern const double c_cm_si					{ 2.997925e10 };
	extern const double kB_eV_Ki				{ 8.613303e-5 };
	extern const double kB_erg_Ki				{ 1.380649e-16 };
	extern const double erg_to_eV				{ 6.2415e11 };
	extern const double accretionEfficiency		{ 1.0 };
	extern const double radiativeEfficiency		{ 1.0 };
	extern const double adiabaticIndex			{ 1.666667 };
	extern const double hydrogen_g				{ 1.673724e-24 };
	extern const double lum_cyg_X1				{ 3.828e37 };
	extern const double T_cyg_X1				{ 1.0e7 };
	extern const double nH_fraction				{ (12.0/27.0) };
}

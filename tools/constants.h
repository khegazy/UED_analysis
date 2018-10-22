#ifndef CONSTANTS_H
#define CONSTANTS_H


// Conversion to a.u.
// website http://physics.nist.gov/cgi-bin/cuu/CCValue?threv|ShowSecond=hartree&First=auefld
#define eV_to_au 0.0367493		// hartree/eV
#define icm_to_au 4.55634e-6		// hartree/imc
#define m_to_au (1/5.291772108e-11) 	// bohr/m
#define angs_to_au (m_to_au*1e-10) 	// bohr/angs
#define s_to_au (1/2.418884326e-17)	// jiffy/s
#define fs_to_au (s_to_au*1e-15)	// jiffy/fs
#define Vpm_to_au (1/5.142206707e11)	// (hartree/(e*bohr))/(V/m)

// Atom Constants
#define Mel 9.1093829e-31		// kg
#define CHARGE_EL = 1.60217662e-19	// Coulombs
#define A_BOHR 5.291772108e-11 		// bohr/m
#define ALPHA 0.0072973525664		// fine structure constant

// Energy Constants
#define H_EV 4.135667662e-15		// eV*s
#define HBAR_EV 6.582119514e-16		// eV*s
#define H_J 6.626070040e-34		// J*s
#define HBAR_J 1.054571800e-34		// J*s
#define EVpHRT 27.2113845 		// eV/hartree

// Other Conversions
#define eV_to_j CHARGE_EL
#define icm_to_eV 1.23981e-4            // eV/cm

// Other Constants
#define PI 3.14159265358979312
#define C_SI 299792458			// m/s
#define C_AU (1/ALPHA)			// m/s
#define C_CMpS 29979245800		// cm/s
#define EPS0_SI 8.854187817e-12		// F/m
#define EPS0_AU (1/(4.0*PI))		// a.u.
#define KBOLTZ_SI 1.38064852e-23	// J/K 
#define KBOLTZ_EV 8.6173324e-5		// eV/K 
#define KBOLTZ_AU (KBOLTZ_EV*eV_to_au)	// hartree/K 
#define GOLD 1.618034
#define TINY 1.0e-20
#define NANVAL 1.234567e-10

/*
#define nmPau 0.05291772108 // nm/bohr
#define icm 8065.54445 // icm/eV
#define hc 1239.84190604789 // 197.326968*2*pi eV nm
#define C 2.99792458e10 // cm/s
#define MC2 931.494e6 // eV
#define hbar 6.58211915e-16
#define kb 8.617343e-5 // eV / K
#define fsPau .02418884326505 // fs / au
#define icmPau (icm * Eh) // icm/hartree
#define amuPau (9.1093826e-31/1.66053886e-27) // unified AMU / au_mass
#define e2P4pieps0 14.3996445 //eV Ang
#define auPe2P4pieps0 (e2P4pieps0 /Eh/a0)  // hartree borh
#define aufor10PW 0.5336 // Atomic unit of field for 10^16 W cm^-1
#define auenergy 183.631526 // 4pi eps_0 * E_au^2 in units of  eV ang^-3
#define mbarn 1e-24 // Mb/cm^2
#define auIntensity 3.55e16 // W/cm^2
#define aupolarizability (a0*a0*a0) // ang^3
*/

#endif


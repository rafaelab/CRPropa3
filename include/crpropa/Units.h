#ifndef CRPROPA_UNITS_H
#define CRPROPA_UNITS_H

#include <cmath>
#include <numbers>

namespace crpropa {

/**
 * \addtogroup PhysicsDefinitions
 * @{
 */

/**
 * @file
 * @brief Definition of SI base units and constants
 *
 * Definition of SI base units and constants used elsewhere in the code
 * Based on:
 * - CODATA recommended values of the fundamental physical constants: 2006
 *   doi:10.1103/RevModPhys.80.633
 * - IAU 2012 Resolution B2, IAU 2015 Resolution B2
 *   https://www.iau.org/administration/resolutions/
 */

// SI units
inline constexpr double meter = 1;
inline constexpr double second = 1;
inline constexpr double kilogram = 1;
inline constexpr double ampere = 1;
inline constexpr double mol = 1;
inline constexpr double kelvin = 1;

// derived units
inline constexpr double newton = 1 * kilogram * meter / second / second;
inline constexpr double pascal = 1 * newton / meter / meter;
inline constexpr double joule = 1 * newton * meter;
inline constexpr double tesla = 1 * newton / ampere / meter;
inline constexpr double volt = 1 * kilogram * meter * meter / ampere / second / second / second;
inline constexpr double coulomb = 1 * ampere * second;
inline constexpr double hertz = 1 / second;
inline constexpr double rad = 1;
inline constexpr double deg = std::numbers::pi / 180.;

// SI Prefixes
inline constexpr double yocto = 1e-24;
inline constexpr double zepto = 1e-21;
inline constexpr double atto = 1e-18;
inline constexpr double femto = 1e-15;
inline constexpr double pico = 1e-12;
inline constexpr double nano = 1e-9;
inline constexpr double micro = 1e-6;
inline constexpr double milli = 1e-3;

inline constexpr double kilo = 1e3;
inline constexpr double mega = 1e6;
inline constexpr double giga = 1e9;
inline constexpr double tera = 1e12;
inline constexpr double peta = 1e15;
inline constexpr double exa = 1e18;
inline constexpr double zetta = 1e21;
inline constexpr double yotta = 1e24;


// physical constants
inline constexpr double eplus = 1.602176487e-19 * ampere * second;
inline constexpr double c_light = 2.99792458e8 * meter / second;
inline constexpr double c_squared = c_light * c_light;
inline constexpr double amu = 1.660538921e-27 * kilogram;
inline constexpr double mass_proton = 1.67262158e-27 * kilogram;
inline constexpr double mass_neutron = 1.67492735e-27 * kilogram;
inline constexpr double mass_electron = 9.10938291e-31 * kilogram;
inline constexpr double h_planck = 6.62606957e-34 * joule * second;
inline constexpr double k_boltzmann = 1.3806488e-23 * joule / kelvin;
inline constexpr double mu0 = 4 * std::numbers::pi * 1e-7 * newton / ampere / ampere;
inline constexpr double epsilon0 = 1.0 / mu0 / c_squared * ampere * second / volt / meter;
inline constexpr double alpha_finestructure = eplus * eplus / 2. / epsilon0 / h_planck / c_light;
inline constexpr double radius_electron = eplus * eplus / 4. / std::numbers::pi / epsilon0 / mass_electron / c_squared;
inline constexpr double sigma_thomson = 8. * std::numbers::pi / 3. * radius_electron * radius_electron;

// gauss
inline constexpr double gauss = 1e-4 * tesla;
inline constexpr double microgauss = 1e-6 * gauss;
inline constexpr double nanogauss = 1e-9 * gauss;
inline constexpr double muG = microgauss;
inline constexpr double nG = nanogauss;

inline constexpr double erg = 1e-7 * joule;

// electron volt
inline constexpr double electronvolt = eplus * volt;
inline constexpr double kiloelectronvolt = 1e3 * electronvolt;
inline constexpr double megaelectronvolt = 1e6 * electronvolt;
inline constexpr double gigaelectronvolt = 1e9 * electronvolt;
inline constexpr double teraelectronvolt = 1e12 * electronvolt;
inline constexpr double petaelectronvolt = 1e15 * electronvolt;
inline constexpr double exaelectronvolt = 1e18 * electronvolt;
inline constexpr double eV = electronvolt;
inline constexpr double keV = kiloelectronvolt;
inline constexpr double MeV = megaelectronvolt;
inline constexpr double GeV = gigaelectronvolt;
inline constexpr double TeV = teraelectronvolt;
inline constexpr double PeV = petaelectronvolt;
inline constexpr double EeV = exaelectronvolt;

inline constexpr double barn = 1e-28 * meter * meter;

// astronomical distances
inline constexpr double au = 149597870700 * meter;
inline constexpr double ly = 365.25 * 24 * 3600 * second * c_light;
inline constexpr double parsec = 648000 / std::numbers::pi * au;
inline constexpr double kiloparsec = 1e3 * parsec;
inline constexpr double megaparsec = 1e6 * parsec;
inline constexpr double gigaparsec = 1e9 * parsec;
inline constexpr double pc = parsec;
inline constexpr double kpc = kiloparsec;
inline constexpr double Mpc = megaparsec;
inline constexpr double Gpc = gigaparsec;

// meter
inline constexpr double kilometer = 1000 * meter;
inline constexpr double centimeter = 0.01 * meter;
inline constexpr double km = kilometer;
inline constexpr double cm = centimeter;

// second
inline constexpr double nanosecond = 1e-9 * second;
inline constexpr double microsecond = 1e-6 * second;
inline constexpr double millisecond = 1e-3 * second;
inline constexpr double minute = 60 * second;
inline constexpr double hour = 3600 * second;
inline constexpr double ns = nanosecond;
inline constexpr double mus = microsecond;
inline constexpr double ms = millisecond;
inline constexpr double sec = second;

// volume
inline constexpr double ccm = cm * cm * cm; // bad name; change to cm3 in the future
inline constexpr double cm3 = cm * cm * cm;

// other derived constants
inline constexpr double mec2 = mass_electron * c_squared;
inline constexpr double h_dirac = h_planck / (2 * std::numbers::pi);


/** @}*/

} // namespace crpropa

#endif // CRPROPA_UNITS_H

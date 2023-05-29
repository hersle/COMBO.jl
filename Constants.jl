module Constants

using Unitful
using UnitfulAstro
import PhysicalConstants

export c, ħ, h, kB, G, α
export me, mH, σT
export km, pc, Mpc, Gpc
export yr, Gyr
export EH1ion, EHe1ion, EHe2ion

# fundamental constants
const c  = PhysicalConstants.CODATA2018.c_0 / u"m/s"
const ħ  = PhysicalConstants.CODATA2018.ħ / u"J*s"
const h  = 2*π * ħ
const kB = PhysicalConstants.CODATA2018.k_B / u"J/K"
const G  = PhysicalConstants.CODATA2018.G / u"m^3/kg/s^2"
const α  = PhysicalConstants.CODATA2018.α # fine structure constant, ≈ 1/137

# masses and Thomson scattering cross section
const me = PhysicalConstants.CODATA2018.m_e / u"kg"
const mH = PhysicalConstants.CODATA2018.m_p / u"kg"
const σT = PhysicalConstants.CODATA2018.σ_e / u"m^2"

# distances
const km = 1e3
const pc = ustrip(uconvert(u"m", 1u"pc")) # pc / m
const Mpc = 1e6 * pc
const Gpc = 1e9 * pc

# times
const yr = 365.25 * 24 * 60 * 60
const Gyr = 1e9 * yr

# ionization energies (from https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page))
const eV = PhysicalConstants.CODATA2018.e / u"C" # J
const EH1ion  = 13.59844 * eV
const EHe1ion = 24.58738 * eV
const EHe2ion = 54.41776 * eV

end

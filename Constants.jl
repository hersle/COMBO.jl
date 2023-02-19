module Constants

using Unitful
using UnitfulAstro
import PhysicalConstants

export c, ħ, kB, G
export km, pc, Mpc, Gpc
export yr, Gyr

const c  = PhysicalConstants.CODATA2018.c_0 / u"m/s"
const ħ  = PhysicalConstants.CODATA2018.ħ / u"J*s"
const kB = PhysicalConstants.CODATA2018.k_B / u"J/K"
const G  = PhysicalConstants.CODATA2018.G / u"m^3/kg/s^2"

# distances
const km = 1e3
const pc = ustrip(uconvert(u"m", 1u"pc")) # pc / m
const Mpc = 1e6 * pc
const Gpc = 1e9 * pc

# times
const yr = 365.25 * 24 * 60 * 60
const Gyr = 1e9 * yr

end

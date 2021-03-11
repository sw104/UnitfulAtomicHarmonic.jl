module UnitfulAtomicHarmonic

import Unitful;
using Unitful: @unit, Dimension, Dimensions, dimension, Units, NoDims, NoUnits, uconvert, ustrip;

export ahunit, ahuconvert, ahustrip, initahu;

"""
    initahu(ω::Number, m::Number)
    initahu(ω::Unitful.Quantity, m::Unitful.Quantity)

Finishes the units initialisation by setting the harmonic frequency base unit
value to `ω` and the mass base unit to `m`. If neither `ω` of `m` have units,
`ω` is assumed to be in units of `rad s⁻¹` and `m` is assumed to be in units of
the atomic mass unit `u`.

Note that conversions from `ω` to SI units should use derived units based on
`s_ahu` account for the factor of `2π` in the conversion.
"""
function initahu(ωh::Unitful.Quantity, mh::Unitful.Quantity)
  # Define base units.
  @eval @unit ħ_u "ħ"   ReducedPlanckConstant   Unitful.ħ                         false;
  @eval @unit ω   "ω"   HarmonicFrequency uconvert(Unitful.rad*Unitful.s^-1, $ωh) false;
  @eval @unit mₕ  "mₕ"  HarmonicMass      uconvert(Unitful.kg, $mh)               false;
  @eval @unit aₕ  "aₕ"  HarmonicLength    (1ħ_u / (1mₕ * 1ω))^(1/2)               false;

  # Define convenience units for handling radians properly in conversions.
  @eval @unit s_ahu   "s"   AngularSecond Unitful.s/(2π*Unitful.rad)  true;
  @eval @unit Hz_ahu  "Hz"  AngularHertz  1/s_ahu                     true;

  # Simplify energy units.
  unit = :(ħ_u*ω);
  @eval ahunit(::typeof(dimension($unit))) = $unit

  # Only register symbols with Unitful when they are all defined.
  Unitful.register(UnitfulAtomicHarmonic);
end
# Assume any dimensionless numbers given are in angular frequency and atomic
# mass units.
initahu(ω::Number, m::Number) = initahu(ω*Unitful.rad*Unitful.s^-1, m*Unitful.u);

"""
  ahunit(x::Unitful.Quantity)
  ahunit(x::Unitful.Units)
  ahunit(x::Unitful.Dimensions)

Returns the appropriate atomic harmonic unit (a `Unitful.Units` object) for the
dimension of `x`.
"""
ahunit(x) = ahunit(dimension(x));

# Harmonic units for each re-defined dimension.
ahunit(x::Dimension{:Length})  = aₕ^x.power;
ahunit(x::Dimension{:Mass})    = mₕ^x.power;
ahunit(x::Dimension{:Time})    = ω^(-x.power);

ahunit(::Dimension{D}) where D = throw(ArgumentError(string(
                                   "no atomic harmonic unit defined for the ",
                                   "dimension $D."
                                 )));

# Define ahunit() when multiple dimensions are passed.
@generated ahunit(::Dimensions{N}) where N = prod(ahunit, N);

# Define ahunit() for dimensionless quantities.
ahunit(::typeof(NoDims)) = NoUnits;

"""
    ahuconvert(x::Unitful.Quantity)

Convert a quantity to appropriate atomic harmonic units.
"""
ahuconvert(x) = uconvert(ahunit(x), x);

"""
    ahuconvert(u::Unitful.Units, x::Number)

Associate the passed units `u` with the passed number `x` which is in atomic
harmonic units.
"""
ahuconvert(u::Units, x::Number) = uconvert(u, x*ahunit(u));

"""
    ahustrip(x::Unitful.Quantity)

Return the value of `x` without its' units.
"""
ahustrip(x) = ustrip(ahuconvert(x));

end

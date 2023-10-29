module Harmonics
using SIMIlluminationPatterns
using SIMIlluminationPatterns: IlluminationPattern
# using SIMIlluminationPatterns.IlluminationPatterns: ParameterEstimator, PE
# using StructuredIlluminationMicroscopy.Utils

using Distributed

abstract type ParameterEstimator{IP<:IlluminationPattern} end
const PE = ParameterEstimator

include("harmonic/peak_frequency.jl")
include("harmonic/autocorrelation_frequency.jl")
include("harmonic/peak_phase_shift.jl")
include("harmonic/cross_correlation_phase_shift.jl")
include("harmonic/IP_cross_correlation_phase_offset.jl")
include("harmonic/cross_correlation_modulation.jl")

export estimate

export AutoCorrelationThroughInputInterpolation, ACTII
export AutoCorrelationThroughOutputInterpolation, ACTOI
export PeakPhaseShift, PPS
export CrossCorrelationPhaseShift, CCPS
export IPCrossCorrelationPhaseOffset, IPCCPO
export CrossCorrelationModulation, CCM

end

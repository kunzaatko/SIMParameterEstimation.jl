module SIMParameterEstimation

using SIMIlluminationPatterns
using SIMIlluminationPatterns: IlluminationPattern
using IterTools
# using SIMIlluminationPatterns.IlluminationPatterns: ParameterEstimator, PE
# using StructuredIlluminationMicroscopy.Utils

using Distributed

abstract type ParameterEstimator{IP<:IlluminationPattern} end
const PE = ParameterEstimator

include("utils.jl")
include("xcorr.jl")
# include("estimators/peak_frequency.jl")
# include("estimators/autocorrelation_frequency.jl")
# include("estimators/peak_phase_shift.jl")
# include("estimators/cross_correlation_phase_shift.jl")
# include("estimators/IP_cross_correlation_phase_offset.jl")
# include("estimators/cross_correlation_modulation.jl")

# export AutoCorrelationThroughInputInterpolation, ACTII
# export AutoCorrelationThroughOutputInterpolation, ACTOI
# export PeakPhaseShift, PPS
# export CrossCorrelationPhaseShift, CCPS
# export IPCrossCorrelationPhaseOffset, IPCCPO
# export CrossCorrelationModulation, CCM

# TODO: The plan for this packages is to instead of creating a million types have an estimator type that holds the
# parameter names that it estimates and "requirements" for the estimation (perhaps) preparation functions and the
# symbolic name that I think of of the estimator. This means that I can create the estimator by a constructor that gets
# the symbol and everything else is determined based on the functions that are overloaded on the symbol (using
# `Val{:name}`). Then a reconstructions can be done by creating a sequence of estimators that are to be applied and then
# maybe some preparation functions yet to be determined, how this will work. The parameters should be determined by some
# type `SIMParams` or similar. <02-09-24> 

include("estimators/frequency_shift/fair_sim.jl")

include("base_overloads.jl")

# export estimate

end

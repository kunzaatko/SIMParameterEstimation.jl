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

include("estimators/frequency_shift/fair_sim.jl")

# export estimate

end

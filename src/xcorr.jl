using ImageFiltering: AbstractBorder, Fill, imfilter!, padarray
using IntervalSets, IntegralArrays, Statistics, ShiftedArrays

# TODO: I really need to make unit tests for this. Or I will end up fucking up the real cross-correlation while
# implementing complex cross-correlation <22-08-24> 
# TODO: raw correlation should really be called covariance when we demean the signal. Similarly the non-demeaned version
# is the second moment <22-08-24> 

struct Window{N}
    inds::CartesianIndices{N}
end
Window(A::AbstractArray) = Window(CartesianIndices(A))
Base.Indices(w::Window) = Tuple(first(d):last(d) for d in w.inds.indices)
Base.broadcasted(::typeof(+), w::Window{N}, ind::CartesianIndex{N}) where {N} = Window(w.inds .+ ind)
Base.broadcasted(::typeof(-), w::Window{N}, ind::CartesianIndex{N}) where {N} = Window(w.inds .- ind)
Base.intersect(w::Window{N}, ws::Window{N}...) where {N} = Window(intersect(w.inds, getfield.(ws, :inds)...))
Base.length(w::Window{N}) where {N} = length(w.inds)

Base.@propagate_inbounds function Base.getindex(A::IntegralArray{T,N}, w::Window{N}) where {T,N}
    A[Tuple(ClosedInterval(first(d), last(d)) for d in w.inds.indices)...]
end

abstract type SummingSet{N} end
struct WindowOverlap{N} <: SummingSet{N}
    signal_win::Window{N}
    template_win::Window{N}
end
function WindowOverlap(signal::AbstractArray{TS,N}, template::AbstractArray{TT,N}) where {N,TS,TT}
    return WindowOverlap(Window(signal), Window(template))
end

# PERF: Is this OK to leave it an abstract type instead of specializing over the array types? <15-08-24> 
struct TrueOverlap{N} <: SummingSet{N}
    # NOTE: It can not be a BitArray, since we need to accept OffsetArrays <02-09-24> 
    signal_supp::AbstractArray{Bool,N}
    template_supp::AbstractArray{Bool,N}
end
# TrueOverlap(signal_supp::, template_supp::BitArray) = TrueOverlap(signal_supp, template_supp)
TrueOverlap(signal::AbstractArray{ST,N}, template::AbstractArray{TT,N}) where {ST,TT,N} = TrueOverlap(signal .!= zero(ST), template .!= zero(TT))

struct Global{N} <: SummingSet{N} end

# TODO: `means!` without the count and `means` could be defined as a single method on the abstract types. A similar
# thing can be done with `overlap_length` <15-08-24> 
overlap_length!(out::AbstractArray, A_win::Window, B_win::Window) = map!(ind -> length(intersect(A_win, B_win .+ ind)), out, A_win.inds)
overlap_length(A_win::Window, B_win::Window) = overlap_length!(similar(A_win.inds, Int), A_win, B_win)
overlap_length!(out::AbstractArray, set::WindowOverlap) = overlap_length!(out, set.signal_win, set.template_win)
overlap_length(set::WindowOverlap) = overlap_length!(similar(set.signal_win.inds, Float32), set)

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    set::WindowOverlap
)
    counts = overlap_length(set.signal_win, set.template_win)
    return means!(signal_mean, template_mean, signal, template, counts, set)
end

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    counts::AbstractArray{Int},
    set::WindowOverlap
) where {ST,TT}
    signal_p = padarray(ST, signal, Fill(zero(ST), Base.Indices(set.template_win)))
    template_p = padarray(TT, template, Fill(zero(TT), Base.Indices(set.signal_win)))
    signal_iA = IntegralArray(signal_p)
    template_iA = IntegralArray(template_p)

    @inbounds @simd for uv in CartesianIndices(signal)
        signal_mean[uv] = signal_iA[set.template_win.+uv]
        template_mean[uv] = template_iA[set.signal_win.-uv]
    end
    @inbounds signal_mean ./= counts
    @inbounds template_mean ./= counts
    return signal_mean, template_mean
end

means(signal::AbstractArray, template::AbstractArray, set::WindowOverlap) =
    means!(similar(signal), similar(signal), signal, template, set)

# FIX: This is repetition that can be avoided by using a sum function with the `set` argument. Then the `means!` and
# `powers!` functions can be written using this function. <15-08-24> 
function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::WindowOverlap
) where {ST,TT}
    signal_energy_p = padarray(ST, signal .^ 2, Fill(zero(ST), Base.Indices(set.template_win)))
    template_energy_p = padarray(TT, template .^ 2, Fill(zero(TT), Base.Indices(set.signal_win)))
    signal_iA = IntegralArray(signal_energy_p)
    template_iA = IntegralArray(template_energy_p)

    @inbounds @simd for uv in CartesianIndices(signal)
        signal_energy[uv] = signal_iA[set.template_win.+uv]
        template_energy[uv] = template_iA[set.signal_win.-uv]
    end
    return signal_energy, template_energy
end

energies(signal::AbstractArray, template::AbstractArray, set::WindowOverlap) =
    energies!(similar(signal), similar(signal), signal, template, set)

# NOTE: There are quite large numerical errors that lead to the `Int` Fourier transform to give `InexactError`s. This is
# not ideal but it should work correctly as the numerical errors are never as large to tip the rounding. <15-08-24> 
overlap_length!(out::AbstractArray, A_supp::AbstractArray{Bool}, B_supp::AbstractArray{Bool}) = round.(Int, imfilter!(out, A_supp, B_supp, Fill(0)))
overlap_length(A_supp::AbstractArray{Bool}, B_supp::AbstractArray{Bool}) = overlap_length!(similar(A_supp, Float32), A_supp, B_supp)
overlap_length!(out::AbstractArray, set::TrueOverlap) = overlap_length!(out, set.signal_supp, set.template_supp)
overlap_length(set::TrueOverlap) = overlap_length!(similar(set.signal_supp, Float32), set)

# FIX: Should this return the mean of the adjoint array? <21-08-24> 
function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    set::TrueOverlap
)
    counts = overlap_length(set.signal_supp, set.template_supp)
    return means!(signal_mean, template_mean, signal, template, counts, set)
end

function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    counts::AbstractArray{Int},
    set::TrueOverlap
) where {ST,TT}
    imfilter!(signal_mean, signal, ST.(set.template_supp), Fill(zero(ST)))
    imfilter!(template_mean, TT.(set.signal_supp), template, Fill(zero(TT)))
    @inbounds signal_mean ./= counts
    @inbounds template_mean ./= counts
    return signal_mean, template_mean
end

means(signal::AbstractArray, template::AbstractArray, set::TrueOverlap) =
    means!(similar(signal), similar(signal), signal, template, set)

abstract type AbstractTrueOverlapSumAlgorithm end
struct Direct <: AbstractTrueOverlapSumAlgorithm end
struct Transform <: AbstractTrueOverlapSumAlgorithm end

function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::TrueOverlap,
    ::Transform
) where {ST,TT}
    # TODO: Test whether this template_supp doesn't have to be reflected before, i.e. whether this is the correlation
    # rather than convolution as it is supposed to be. <19-08-24> 
    imfilter!(signal_energy, signal, set.template_supp, Fill(zero(ST)))
    imfilter!(template_energy, TT.(set.signal_supp), template, Fill(zero(TT)))
    return signal_energy, template_energy
end

# PERF: This can be approached in multiple ways, we can directly sum the intersection of the `signal_supp` and
# `template_supp` or a similar method to the window approach can be used. This would lead to better performance if the
# areas under the intersection are large, hence the border is much smaller in terms of pixel counts.  <15-08-24> 
function energies!(
    signal_energy::AbstractArray,
    template_energy::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray{TT},
    set::TrueOverlap,
    ::Direct
) where {ST,TT}
    signal_supp_p = padarray(Bool, set.signal_supp, Fill(false, set.template_supp))
    template_supp_p = padarray(Bool, set.template_supp, Fill(false, set.signal_supp))
    signal_energy_p = padarray(ST, signal .^ 2, Fill(zero(ST), set.template_supp))
    template_energy_p = padarray(TT, template .^ 2, Fill(zero(ST), set.signal_supp))
    @inbounds @simd for uv in CartesianIndices(signal)
        signal_energy[uv] = sum(signal_energy_p[signal_supp_p.*ShiftedArray(template_supp_p, Tuple(uv); default=false)])
        template_energy[uv] = sum(template_energy_p[signal_supp_p.*ShiftedArray(template_supp_p, Tuple(-uv); default=false)])
    end
    return signal_energy, template_energy
end

energies!(signal_energy::AbstractArray, template_energy::AbstractArray, signal::AbstractArray, template::AbstractArray, set::TrueOverlap) =
    energies!(signal_energy, template_energy, signal, template, set, Transform())
energies(signal::AbstractArray, template::AbstractArray, set::TrueOverlap) =
    energies!(similar(signal), similar(signal), signal, template, set)

overlap_length!(out::AbstractArray, ::Global) = fill!(out, length(out))
function means!(
    signal_mean::AbstractArray,
    template_mean::AbstractArray,
    signal::AbstractArray,
    template::AbstractArray,
    ::Global
)
    signal_mean .= mean(signal)
    template_mean .= mean(template)
    return signal_mean, template_mean
end

# TODO: It doesn't make any sense to calculate the `xcorr` for shifts where there is no overlap <22-08-24> 
# TODO: There should also be a `range` argument, which determines on which indices the cross-correlation is done. Even
# better, there should be an indices/shifts argument that determines where to calculate the xcorr <10-08-24> 
function xcorr!(
    out::AbstractArray,
    signal::AbstractArray{ST},
    template::AbstractArray;
    normset::SummingSet=TrueOverlap(signal, template)
) where {ST}
    counts = overlap_length(normset)
    signal_mean, template_mean = means(signal, template, normset)

    imfilter!(out, signal, template, Fill(zero(ST)))
    out .-= (counts .* signal_mean .* template_mean) # Demeaning

    signal_energy, template_energy = energies(signal, template, normset)
    signal_denom, template_denom = signal_energy .- (counts .* signal_mean .^ 2), template_energy .- (counts .* template_mean .^ 2)
    norm_factor!(signal_denom, signal_denom)
    norm_factor!(template_denom, template_denom)
    denom = signal_denom .* template_denom

    out ./= denom # Normalization
    return out
end
xcorr(signal::AbstractArray, template::AbstractArray; vargs...) =
    xcorr!(similar(signal), signal, template; vargs...)

function norm_factor!(out, denom::AbstractArray{ST}) where {ST<:Real}
    out .= sqrt.(denom)
    return out
end
function norm_factor!(out, denom::AbstractArray{ST}) where {ST<:Complex}
    # FIX: Is this correct? <22-08-24> 
    out .= abs.(denom)
    return out
end

include("xcorr_docs.jl")

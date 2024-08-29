# ```jldoctest xcorr; setup = :(using OffsetArrays; using SIMParameterEstimation: Window)
"""
    Window(arr::AbstractArray)

Window of indices of the array `arr`, representing a Cartesian product of the contained index intervals. 

Useful for computing the overlap of two signals during cross-correlation.


```julia-repl
julia> signal = [0 0 0 0 0 0 0 0 0;
                 0 0 0 1 1 1 0 0 0;
                 0 0 1 1 1 1 1 0 0;
                 0 1 1 1 1 1 1 1 0;
                 0 1 1 1 1 1 1 1 0;
                 0 0 1 1 1 1 1 0 0;
                 0 0 0 1 1 1 0 0 0;
                 0 0 0 0 0 0 0 0 0];

julia> signal = OffsetArrays.Origin(-4,-4)(signal)
8×9 OffsetArray(::Matrix{Int64}, -4:3, -4:4) with eltype Int64 with indices -4:3×-4:4:
 0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  0  0  0
 0  0  1  1  1  1  1  0  0
 0  1  1  1  1  1  1  1  0
 0  1  1  1  1  1  1  1  0
 0  0  1  1  1  1  1  0  0
 0  0  0  1  1  1  0  0  0
 0  0  0  0  0  0  0  0  0

julia> w_signal = Window(signal)
Window{2}(CartesianIndices((OffsetArrays.IdOffsetRange(values=-4:3, indices=-4:3), OffsetArrays.IdOffsetRange(values=-4:4, indices=-4:4))))

julia> template = [0 1 0;
                   1 1 1;
                   0 1 0];

julia> template = OffsetArrays.centered(template)
3×3 OffsetArray(::Matrix{Int64}, -1:1, -1:1) with eltype Int64 with indices -1:1×-1:1:
 0  1  0
 1  1  1
 0  1  0

julia> w_template = Window(template)
Window{2}(CartesianIndices((OffsetArrays.IdOffsetRange(values=-1:1, indices=-1:1), OffsetArrays.IdOffsetRange(values=-1:1, indices=-1:1))))

julia> @assert length(w_template) == 9

julia> intersect(w_signal, w_template .+ CartesianIndex((3, 3)))
Window{2}(CartesianIndices((2:3, 2:4)))
```
"""
Window


"""
     mean!(out::AbstractArray, f::AbstractArray, inds::SummingSet)

 Compute the mean of `f` for each of its index with the given indice set `inds`.

# Examples 

 """
mean!(out::AbstractArray, f::AbstractArray, inds::SummingSet)
# ```jldoctest xcorr
#
# ```


"""
    mean!(out::AbstractArray, A::AbstractArray, set::SummingSet)

Compute the mean of array `A` for each shift  with the specified `set` that depends on the
template for the cross-correlation.

__Warning__: When the template does not contain the origin (`CartesianIndex(0,...)`) in 
    its axes, it may return `NaN`s wherever the shifted template does not overlap with 
    the signal.
"""
mean!

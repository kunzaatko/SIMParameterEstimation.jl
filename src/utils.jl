# TODO: Documentation <29-08-24> 
function mixin_matrix(
    T::Type{<:Real},
    ϕ::NTuple{N,<:Real},
    μ::NTuple{N,<:Real}=ntuple(_ -> 1, Val(N))
) where {N}
    M = ones(Complex{T}, N, 3)
    M[:, 2] = collect(@. μ / 2 * exp(im * ϕ))
    M[:, 3] = collect(@. μ / 2 * exp(-im * ϕ))

    return M
end
# NOTE: There are unbound type parameters here, but it does not matter really <29-08-24> 
mixin_matrix(ϕ::NTuple{N,TP}, μ::NTuple{N,TM}=ntuple(_ -> 1, Val(N))) where {N,TP<:Real,TM<:Real} = mixin_matrix(promote_type(TP, TM), ϕ, μ)
mixin_matrix(ϕ_0::Real, μ::NTuple{N,T}) where {N,T<:Real} = mixin_matrix(ϕ_0 .+ Tuple(LinRange(0, 2π, N + 1)[begin:(end-1)]), μ)
mixin_matrix(ϕ_0::Real, N::Int) = mixin_matrix(ϕ_0, ntuple(_ -> 1, Val(N)))

function separation_matrix(
    args...
)
    # PERF: Could be faster and more precise, if the matrix was created from the analytical inversion
    # FIX: This does not work for a non-square matrix. It needs to be handled differently for N != 3. Perhaps it will
    # lead to more variants of components -1 and 1 and then they will be averaged?? <29-08-24> 
    return inv(mixin_matrix(args...))
end

# TODO: There should be a mutating version. Some times we only want the components and do not need to keep them. Mixin
# should be symmetric to this. <29-08-24> 
# TODO: Docs. Should include the fact that it is a 4 dim array by default <29-08-24> 
function mixin_components(f_imgs::AbstractArray{<:Number,3}, M::AbstractMatrix{<:Complex})
    nphases = size(M, 1)
    @assert mod(size(f_imgs, 3), nphases) == 0 """The number of images supplied must be a multiple of the number of phases.
    The complete stack of images should be of `size(f_imgs, 3)` == `norientations`×`nphases`"""
    comps = stack(IterTools.partition(eachslice(f_imgs, dims=3), nphases, nphases)) do single_orientation
        stack(row -> sum(single_orientation .* row), eachrow(M))
    end
    return comps
end
mixin_components(f_imgs::AbstractArray{<:Number,4}, M::AbstractMatrix{<:Complex}) = mixin_components(reshape(f_imgs, size(f_imgs, 1), size(f_imgs, 2), :), M)

# TODO: There should be a mutating version. Some times we only want the components and do not need to keep them. Mixin
# should be symmetric to this. <29-08-24> 
# TODO: Docs. Should include the fact that it is a 4 dim array by default <29-08-24> 
function separate_components(f_imgs::AbstractArray{<:Number,3}, M_inv::AbstractMatrix{<:Complex})
    # NOTE: This currently only works if there are 3 phases... <29-08-24> 
    nphases = size(M_inv, 1)
    @assert mod(size(f_imgs, 3), nphases) == 0 """The number of images supplied must be a multiple of the number of phases.
    The complete stack of images should be of `size(f_imgs, 3)` == `norientations`×`nphases`"""
    comps = stack(IterTools.partition(eachslice(f_imgs, dims=3), nphases, nphases)) do single_orientation
        stack(row -> sum(single_orientation .* row), eachrow(M_inv))
    end
    return comps
end
separate_components(f_imgs::AbstractArray{<:Number,4}, M_inv::AbstractMatrix{<:Complex}) = separate_components(reshape(f_imgs, size(f_imgs, 1), size(f_imgs, 2), :), M_inv)

include("utils_docs.jl")

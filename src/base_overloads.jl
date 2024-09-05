for T in (:Window, :WindowOverlap, :TrueOverlap)
    # NOTE: Taken from Distributions.jl <kunzaatko> 
    for func in (:(==), :isequal, :isapprox)
        @eval function Base.$func(thing_1::$T, thing_2::$T; kwargs...)
            for f in fieldnames($T)
                isdefined(thing_1, f) && isdefined(thing_2, f) || return false
                # perform equivalence check to support types that have no defined equality, such
                # as `missing`
                getfield(thing_1, f) === getfield(thing_2, f) || $func(getfield(thing_1, f), getfield(thing_2, f); kwargs...) || return false
            end

            return true
        end
    end

    # NOTE: Taken from Distributions.jl <kunzaatko> 
    @eval function Base.hash(thing::$T, h::UInt)
        hashed = hash(IlluminationPattern, h)
        hashed = hash(nameof($T), hashed)

        for f in fieldnames($T)
            hashed = hash(getfield(thing, f), hashed)
        end

        return hashed
    end
end

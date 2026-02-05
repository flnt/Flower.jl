  """
            safefill(T; sources=(;), defaults=nothing)

        Constructs a struct `T` (like `Numerical`) using field values from one or more
        sources (`NamedTuple`, `struct`, or `Dict`). If a field is missing, uses the
        default from `T()`.
        """
        function safefill(T; sources=(;), defaults=nothing)
            defaults = isnothing(defaults) ? T() : defaults
            fnames = fieldnames(T)
            kwargs = Dict{Symbol, Any}()

            # Try to find each field in sources
            for f in fnames
                found = false
                for src in sources
                    if hasproperty(src, f)
                        kwargs[f] = getproperty(src, f)
                        found = true
                        break
                    elseif src isa AbstractDict && haskey(src, f)
                        kwargs[f] = src[f]
                        found = true
                        break
                    end
                end
                if !found
                    kwargs[f] = getproperty(defaults, f)
                end
            end

            return T(; kwargs...)   # ✅ keyword construction (works with @with_kw)
        end

        function safefill_with_aliases(::Type{Numerical{T,D}}, sim, phys, io,aliases) where {T<:Real, D<:Integer}
            default = Numerical{T,D}()  # construct default parametric instance
          
            args = Dict{Symbol,Any}()

            for field in fieldnames(Numerical{T,D})
                srcsym = haskey(aliases, field) ? aliases[field] : field

                val = nothing
                if hasproperty(sim, srcsym)
                    val = getproperty(sim, srcsym)
                elseif hasproperty(phys, srcsym)
                    val = getproperty(phys, srcsym)
                elseif hasproperty(io, srcsym)
                    val = getproperty(io, srcsym)
                end

                args[field] = val === nothing ? getproperty(default, field) : val
            end

            # Use keyword or positional constructor depending on your struct definition
            return Numerical{T,D}(; args...)
        end

        function safefill_with_aliases_and_extra(::Type{Numerical{T,D}}, sim, phys, io, aliases, extra) where {T<:Real, D<:Integer}
            default = Numerical{T,D}()
            args = Dict{Symbol,Any}()

            for field in fieldnames(Numerical{T,D})
                srcsym = get(aliases, field, field)

                val =
                    hasproperty(sim,  srcsym)  ? getproperty(sim, srcsym)  :
                    hasproperty(phys, srcsym)  ? getproperty(phys, srcsym) :
                    hasproperty(io,   srcsym)  ? getproperty(io, srcsym)   :
                    haskey(extra, srcsym)      ? extra[srcsym]             :
                    nothing

                args[field] = val === nothing ? getproperty(default, field) : val
            end

            return Numerical{T,D}(; args...)
        end

        function safefill_with_aliases_and_extra_already_init(::Type{Numerical{T,D}}, default,sim, phys, io, aliases, extra) where {T<:Real, D<:Integer}
            # default = Numerical{T,D}()
            args = Dict{Symbol,Any}()

            for field in fieldnames(Numerical{T,D})
                srcsym = get(aliases, field, field)

                val =
                    hasproperty(sim,  srcsym)  ? getproperty(sim, srcsym)  :
                    hasproperty(phys, srcsym)  ? getproperty(phys, srcsym) :
                    hasproperty(io,   srcsym)  ? getproperty(io, srcsym)   :
                    haskey(extra, srcsym)      ? extra[srcsym]             :
                    nothing

                args[field] = val === nothing ? getproperty(default, field) : val
            end

            return Numerical{T,D}(; args...)
        end
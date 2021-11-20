module SigmaShift

using MITgcmTools, Dierckx, Interpolations

export sigma, vars2sigma1, sigma1column, var2sigmacolumn, sigma1grid
export mixinversions!, dedup!

"""
    function vars2sigma1(vars,p,sig1grid,γ,spline_order)
    map variables from regular 3D grid onto sigma1 surfaces
# Arguments
- `vars::Dict{String,Array{T,3}}}`: dict of 3d arrays
- `p::Array{T,1}` : vertical profile of standard pressures
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline
# Output
- `varsσ::Dict{String,Array{T,3}`: dict of 3d arrays of variables on sigma1 surfaces
"""
function vars2sigma1(vars::Dict{String,Array{T,3}},pressure::Vector{T},σ₁grid::Vector{T},splorder::Integer) where T<:AbstractFloat

    # is there a problem if pressure is not the same type as the input vars?
    # could introduce two parametric types to function definition above
    
    # θ and S must exist

    # a list of potential temperature names that are recognized
    θkeys = ("THETA","θ","theta","T")
    Skeys = ("SALT","S","Sp")

    θname = string()

    nx = Integer(); ny = Integer(); nz = Integer()
    
    for k in θkeys
        if haskey(vars,k)
            nx,ny,nz = size(vars[k])
            θname = k
        end
    end
    θname == "" && error("Potential temperature missing")

    Sname = string()
    for k in Skeys
        if haskey(vars,k)
            nx,ny,nz = size(vars[k])
            Sname = k
        end
    end
    Sname == "" && error("Salinity missing")

    # loop over faces
    nσ = length(σ₁grid)

    # vcol = Dict with profile/column data
    # pre-allocate each key in vars
    vcol = Dict{String,Array{T,1}}() # vars in a column
    varsσ = Dict{String,Array{T,3}}()

    for (key, value) in vars
        vcol[key] = fill(convert(T,NaN),nz)
        varsσ[key] = fill(convert(T,NaN),(nx,ny,nσ))
    end

    # allocate standard pressure by hand.
    varsσ["p"] = fill(convert(T,NaN),(nx,ny,nσ))

    for xx = 1:nx
        for yy = 1:ny
            for (vcolname, vcolval) in vars
                # vcol = Dict with profile/column data
                vcol[vcolname] = vcolval[xx,yy,:]
            end

            # also need to filter dry values and to change zz
            # Consider using `isdry` function and dryval in future.
            nw = count(notnanorzero,vcol[θname]) # number of wet points in column
            nwS = count(notnanorzero,vcol[Sname])
            if nw != nwS
                error("T,S zeroes inconsistent")
            end

            if nw > 3
                # incurs error if splorder > number of points in column
                # if nw > splorder #need >=n+1 points to do order-n interpolation
                σ₁=sigma1column(vcol[θname][1:nw],vcol[Sname][1:nw],pressure[1:nw])

                for (vckey,vcval) in vcol
                    varσ = var2sigmacolumn(σ₁,vcval[1:nw],σ₁grid,splorder)
                    [varsσ[vckey][xx,yy,ss] = convert(T,varσ[ss]) for ss = 1:nσ]
                end

                # do standard pressure by hand.
                pσ = var2sigmacolumn(σ₁,pressure[1:nw],σ₁grid,splorder)
                [varsσ["p"][xx,yy,ss] = convert(T,pσ[ss]) for ss = 1:nσ]

            end
        end
    end
    return varsσ
end

"""
    function sigmacolumn(θ,S,p,p0)
    σ for a water column
# Arguments
- `θz::Vector{T}`: potential temperature
- `Sz::Vector{T}`: practical salinity
- `pz::Vector{T}`: vertical profile of standard pressures
- `p0::Integer`: reference pressure
# Output
- `σ::Vector{T}`:  sigma for wet points in column
"""
function sigmacolumn(θz::Vector{T},Sz::Vector{T},pz::Vector{T2},p0::Integer)::Vector{T} where T<:AbstractFloat where T2<:AbstractFloat
    nz = length(θz)
    σ = similar(θz)
    #p0T = convert(T,p0)
    
    # hard coded for sigma1
    σa,σb,σc = MITgcmTools.SeaWaterDensity(θz,Sz,pz,p0)
    [σ[zz] = convert(T,σc[zz]) .- 1000.0 for zz = 1:nz]
    return σ
end

"""
    function sigma1column(θ,S,p)
    σ₁ for a water column
# Arguments
- `θz::Array{Float,1}}`: potential temperature
- `Sz::Array{Float,1}}`: practical salinity
- `pz::Array{Float,1}`: vertical profile of standard pressures
# Output
- `σ₁`:  sigma-1 for wet points in column
"""
sigma1column(θz,Sz,pz) = sigmacolumn(θz,Sz,pz,1000)

notnanorzero(z) = !iszero(z) && !isnan(z)

""" function sigma1grid()
    Standard choice of sigma1 surfaces
"""
function sigma1grid()
    σ₁grida = 24:0.05:31
    σ₁gridb = 31.02:0.02:33
    σ₁grid = vcat(σ₁grida,σ₁gridb)
    σ₁grid = σ₁grid[1:3:end]
    return σ₁grid
end


"""
    function var2sigmacolumn(σ,v,σ₁grid,splorder)
    map θ,S, p onto σ₁ surfaces for a water column
# Arguments
- `σ`::Array{Float,1}}`: sigma values of input variable
- `v::Array{Float,1}}`: variable of interest
- `sig1`: σ surface values
- `splorder`: 1-5, order of spline
# Output
- `θonσ`: variable on sig1 sigma surfaces
"""
function var2sigmacolumn(σorig::Vector{T},v,σgrid,splorder::Integer) where T<:AbstractFloat where T2<:AbstractFloat 
    # choose a univariate spline with s = magic number
    #θspl = Spline1D(σ₁,θz;k=splorder,s=length(σ₁))

    linearinterp = false

    σ = copy(σorig) # make sure sigma-1 doesn't mutate and pass back

    nσout = length(σgrid)
    θonσ = fill(convert(T,NaN),nσout)

    # some velocity variables may include NaN's
    # add this kludge to make NaN velocity equal to zero
    # a useful approximation to no-flow boundary condition?
    replace!(v,convert(T,NaN)=>zero(T))
    
    # eliminate homogeneities
    dedup!(σ,v)

    # mix inversions
    mixinversions!(σ,v)

    nσin = length(σ)

    # 1) no inversions or homogeneity (this constraint relaxed now because too many profiles thrown out)
    # 2) range of sig1, 3) no extrapolation
    # if sum(diff(σ₁).<0)==0 && count(minimum(σ₁).<=σ₁grid.<=maximum(σ₁)) > 0
    # added nσin > 1 to fix error
    if nσin > 1 && count(minimum(σ) .<= σgrid .<= maximum(σ)) > 0

        # eliminate any extrapolation
        sgood = findall(minimum(σ).<=σgrid.<= maximum(σ))

        if nσin > splorder
            θspl = Spline1D(σ,v;k=splorder)
            for ss in sgood
                θonσ[ss] = θspl(σgrid[ss])
            end

            # check for spline instability

            # this version didn't work if length(sgood)==1
            #if maximum(θonσ[sgood]) - minimum(θonσ[sgood]) > 1.00 * (maximum(v) - minimum(v))

            # give some leeway with "0.1"
            if (maximum(θonσ[sgood]) > maximum(v) + 0.1*(maximum(v)-minimum(v)) ||
                minimum(θonσ[sgood]) < minimum(v) - 0.1*(maximum(v)-minimum(v)))
                linearinterp = true
                #println("unstable spline")
            end
        else # spline interp
            linearinterp = true
        end 

        if linearinterp
            #println("doing linear interp")
#            println(size(σ),size(v))
            interp_linear = LinearInterpolation(σ, v)
            for ss in sgood
                θonσ[ss] = interp_linear(σgrid[ss])
            end
        end # linearinterp

    else
        #println("not doing any calcs")
    end # any good points?

    return θonσ
end

"""
function mixinversions!(a,b)
"""
function mixinversions!(a,b)
    while sum(diff(a).<=0) > 0
        length(a) == 1 ? da = 1. : da = diff(a)
        ii = findfirst(isnotpositive,da)
        a[ii] = (a[ii] + a[ii+1])/2
        b[ii] = (b[ii] + b[ii+1])/2
        deleteat!(a,ii+1)
        deleteat!(b,ii+1)
    end
end

"""
function dedup!(a,b)
"""
function dedup!(a,b)
    length(a) == 1 ? da = 1. : da = diff(a)
    while count(iszero,da) > 0
        dedupfirst!(a,b)
        length(a) ==1 ? da = 1. : da = diff(a)
    end
end

function dedupfirst!(a,b)
    counter = 1
    da = diff(a) # requires length of 2 or more
    ii = findfirst(iszero,da)
    while iszero(da[ii])
        b[ii+1] += b[ii]
        deleteat!(a,ii)
        deleteat!(b,ii)
        deleteat!(da,ii)
        counter += 1
        (length(da) < ii) ? break : nothing
    end
    b[ii] /= counter
end

isnotpositive(x) = (abs(x) == -x)

end

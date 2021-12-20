module IsopycnalSurfaces

using Dierckx, Interpolations

export sigma0column, sigma1column, sigma2column,
 vars2sigma1,  var2sigmacolumn, sigma1grid,
 mixinversions!, dedup!, density, sigmacolumn

"""
    function vars2sigma1(vars,p,sig1grid,γ,spline_order)
    map variables from regular 3D grid onto sigma1 surfaces
# Arguments
- `vars::Dict{String,Array{T,3}}}`: dict of 3d arrays
- `p::Array{T,1}` : vertical profile of standard pressures
- `sig1grid`: σ₁ surface values
- `γ`: grid description needed for preallocation
- `splorder`: 1-5, order of spline, optional keyword argument, default=3
- `linearinterp`: optional keyword logical argument, default=false
# Output
- `varsσ::Dict{String,Array{T,3}`: dict of 3d arrays of variables on sigma1 surfaces
"""
function vars2sigma1(vars::Dict{String,Array{T,3}},pressure::Vector{T},σ₁grid::Vector{T};splorder=3,linearinterp=false) where T<:AbstractFloat

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
                    varσ = var2sigmacolumn(σ₁,vcval[1:nw],σ₁grid,splorder=splorder,linearinterp=linearinterp)
                    [varsσ[vckey][xx,yy,ss] = convert(T,varσ[ss]) for ss = 1:nσ]
                end

                # do standard pressure by hand.
                pσ = var2sigmacolumn(σ₁,pressure[1:nw],σ₁grid,splorder=splorder,linearinterp=linearinterp)
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
function sigmacolumn(θz::Vector{T},Sz::Vector{T},pz::Vector{T2},p0::Integer, eos="EOS80")::Vector{T} where T<:AbstractFloat where T2<:AbstractFloat
    nz = length(θz)
    σ = similar(θz)
    
    # choose EOS method, added by Ray Dec 09 2021
    if eos == "EOS80" # unesco, saunders et al., 1980
        σ = density.(Sz,θz,p0) .- 1000.
    elseif eos == "JMD95" # Jackett McDougall 1995, JAOT
        σa,σb,σc = SeaWaterDensity(θz,Sz,pz,p0)
        [σ[zz] = convert(T,σc[zz]) .- 1000.0 for zz = 1:nz]
    else
        error("The entered EOS is not supported currently, please try the supported one, like EOS94")
    end

    return σ
end

"""
    function sigma0column(θ,S,p)
    σ₀ for a water column
    Untested for a mix of float values

# Arguments
- `θz::Vector{T}`: potential temperature
- `Sz::Vector{T}`: practical salinity
- `pz::Vector{T}`: vertical profile of standard pressures
# Output
- `σ₀`:  sigma-0 for wet points in column
"""
sigma0column(θz,Sz,pz,eos=missing) = ismissing(eos) ? sigmacolumn(θz,Sz,pz,0) : sigmacolumn(θz,Sz,pz,0,eos)

"""
    function sigma1column(θ,S,p)
    σ₁ for a water column
    Untested for a mix of float values

# Arguments
- `θz::Vector{T}`: potential temperature
- `Sz::Vector{T}`: practical salinity
- `pz::Vector{T}`: vertical profile of standard pressures
# Output
- `σ₁`:  sigma-1 for wet points in column
"""
sigma1column(θz,Sz,pz,eos=missing) = ismissing(eos) ? sigmacolumn(θz,Sz,pz,1000) : sigmacolumn(θz,Sz,pz,1000,eos)

"""
    function sigma2column(θ,S,p)
    σ₂ for a water column
    Untested for a mix of float values

# Arguments
- `θz::Vector{T}`: potential temperature
- `Sz::Vector{T}`: practical salinity
- `pz::Vector{T}`: vertical profile of standard pressures
# Output
- `σ₂`:  sigma-2 for wet points in column
"""
# revised by Ray, Dec 09 2021
sigma2column(θz,Sz,pz,eos=missing) = ismissing(eos) ? sigmacolumn(θz,Sz,pz,2000) : sigmacolumn(θz,Sz,pz,2000,eos)

"""
   function notnanorzero

     true is argument is not a NaN nor zero
"""
notnanorzero(z) = !iszero(z) && !isnan(z)

""" function sigma1grid()
    Standard (from Susan Wijffels, WHOI) choice of sigma1 surfaces
# Arguments
- `z`: value
# Output
- `σ₁grid`: list (vector) of σ₁ values
"""
function sigma1grid()
    σ₁grida = 24:0.05:31
    σ₁gridb = 31.02:0.02:33
    σ₁grid = vcat(σ₁grida,σ₁gridb)
    σ₁grid = σ₁grid[1:3:end]
    return σ₁grid
end

"""
    function var2sigmacolumn(σ,v,σ₁grid;splorder,linearinterp)
    map θ,S, p onto σ₁ surfaces for a water column
# Arguments
- `σ`::Array{Float,1}}`: sigma values of input variable
- `v::Array{Float,1}}`: variable of interest
- `sig1`: σ surface values
- `splorder`: optional argument of 1-5, order of spline, default=3
- `linearinterp`: optional argument, true to force linear interpolation, default=false
# Output
- `θonσ`: variable on sig1 sigma surfaces
"""
function var2sigmacolumn(σorig::Vector{T},v,σgrid; splorder=3,linearinterp=false) where T<:AbstractFloat where T2<:AbstractFloat 
    # choose a univariate spline with s = magic number
    #θspl = Spline1D(σ₁,θz;k=splorder,s=length(σ₁))

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

For values of `a` that are not increasing, 
an inversion is defined.
Average these values of a until values of `a` are sorted.
Do the same averaging on the accompanying vector `b`.

# Arguments
- `a::Vector{T}`: a density variable
- `b::Vector{T}`: an accompanying tracer variable

Arguments are mutated by the function. 
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

Remove values of `a` that are duplicates.
Remove values of `b` that have the same location as the duplicates in `a`.
The length of `a` and `b` should be identical before and after invoking this function. 

# Arguments
- `a::Vector{T}`: a density variable
- `b::Vector{T}`: an accompanying tracer variable

Arguments are mutated by the function. 
"""
function dedup!(a,b)
    length(a) == 1 ? da = 1. : da = diff(a)
    while count(iszero,da) > 0
        dedupfirst!(a,b)
        length(a) ==1 ? da = 1. : da = diff(a)
    end
end

"""
function dedup!(a,b)

Remove the first duplicate of `a`.
Remove the value of `b` that is located at the same entry as the first duplicate in `a`.
The length of `a` and `b` should be identical before and after invoking this function. 

# Arguments
- `a::Vector{T}`: a density variable
- `b::Vector{T}`: an accompanying tracer variable

Arguments are mutated by the function. 
"""
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


"""
SeaWaterDensity(Θ,Σ,Π,Π0) from MITgcmTools.jl/PhysicalOceanography.jl, From Gael Forget

Compute potential density (ρP), in situ density (ρI), and density
referenced to PREF (Π0 in decibars) from potential temperature (Θ in °C),
salinity (Σ in psu) and pressure (Π in decibars) according to the
UNESCO / Jackett & McDougall 1994 equation of state.

Credits: code based on a Matlab implementation by B. Ferron

Reference: https://www.jodc.go.jp/info/ioc_doc/UNESCO_tech/059832eb.pdf

Check value: ρI = `1041.83267kg/m^3` for Θ=`3°Celcius`, Σ=`35psu`, Π=`3000dbar`
```
(ρP,ρI,ρR) = SeaWaterDensity(3.,35.5,3000.)
isapprox(ρI,1041.83267, rtol=1e-6)
```
"""
function SeaWaterDensity(Θ,Σ,Π,Π0=missing)

   #square root salinity
   sqrtΣ= sqrt.(Σ)
   #compute density pure water at atm pressure
   ZR1= ((((6.536332E-9*Θ .-1.120083E-6).*Θ .+1.001685E-4).*Θ
   .-9.095290E-3).*Θ .+6.793952E-2).*Θ .+999.842594
   #seawater density atm pressure
   ZR2= (((5.3875E-9*Θ .-8.2467E-7).*Θ .+7.6438E-5).*Θ
   .-4.0899E-3).*Θ .+0.824493
   ZR3= (-1.6546E-6*Θ .+1.0227E-4).*Θ .-5.72466E-3
   ZR4= 4.8314E-4

   #potential density (referenced to the surface)
   ρP= (ZR4*Σ + ZR3.*sqrtΣ + ZR2).*Σ + ZR1

   #add the compression terms
   ZE = (-3.508914E-8*Θ .-1.248266E-8).*Θ .-2.595994E-6
   ZBW= ( 1.296821E-6*Θ .-5.782165E-9).*Θ .+1.045941E-4
   ZB = ZBW + ZE .* Σ

   ZD = -2.042967E-2
   ZC = (-7.267926E-5*Θ .+2.598241E-3).*Θ .+0.1571896
   ZAW= ((5.939910E-6*Θ .+2.512549E-3).*Θ .-0.1028859).*Θ .-4.721788
   ZA = ( ZD*sqrtΣ + ZC).*Σ + ZAW

   ZB1= (-0.1909078*Θ .+7.390729).*Θ .-55.87545
   ZA1= ((2.326469E-3*Θ .+1.553190).*Θ .-65.00517).*Θ .+1044.077
   ZKW= (((-1.361629E-4*Θ .-1.852732E-2).*Θ .-30.41638).*Θ
   .+2098.925).*Θ .+190925.6
   ZK0= (ZB1.*sqrtΣ + ZA1).*Σ + ZKW

   #in situ density
   ρI = ρP ./ (1.0 .-Π./(ZK0-Π.*(ZA-Π.*ZB)))

   #density referenced to level Π0
   if !ismissing(Π0)
      ρR = ρP ./ (1.0 .-Π0./(ZK0-Π0.*(ZA-Π0.*ZB)))
   else
      ρR = ρP
   end

   return ρP,ρI,ρR
end



"""
density(S,T,p) and dependent functions from PhysOcean.jl/EOS80.jl, From Alexander Barth

Compute the density of sea-water (kg/m³) at the salinity `S` (psu, PSS-78), temperature `T` (degree Celsius, ITS-90) and pressure `p` (decibar) using the UNESCO 1983 polynomial.

Reference: Fofonoff, N.P.; Millard, R.C. (1983). Algorithms for computation of fundamental properties of seawater. UNESCO Technical Papers in Marine Science, No. 44. UNESCO: Paris. 53 pp.
http://web.archive.org/web/20170103000527/http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

```
Check value: ρI = `1041.87651kg/m^3` for Θ=`3°Celcius`, Σ=`35.5psu`, Π=`3000dbar`
ρI = density(35.5, 3, 3000)
```
added by Ray, Dec 09, 2021
"""
function density(S,T,p)
    ρ = density0(S,T)

    if (p == 0)
        return ρ
    end

    K = secant_bulk_modulus(S,T,p)
    # convert decibars to bars
    p = p/10
    return ρ / (1 - p/K)
end

"""
    temperature68(T)
Convert temperature `T` from ITS-90 scale to the IPTS-68 scale following Saunders, 1990.
Saunders, P.M. 1990, The International Temperature Scale of 1990, ITS-90. No.10, p.10.
https://web.archive.org/web/20170304194831/http://webapp1.dlib.indiana.edu/virtual_disk_library/index.cgi/4955867/FID474/wocedocs/newsltr/news10/news10.pdf
"""
temperature68(T) = 1.00024 * T

"""
    density_reference_pure_water(T)
density of pure water at the temperature `T` (degree Celsius, ITS-90)
"""
function density_reference_pure_water(T)

    t = temperature68(T)

    # page 21, equation 14 of
    # http://web.archive.org/web/20170103000527/http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    a0 = 999.842594 # why [-28.263737]
    a1 = 6.793952e-2
    a2 = -9.095290e-3
    a3 = 1.001685e-4
    a4 = -1.120083e-6
    a5 = 6.536332e-9

    ρ_w = a0 + (a1 + (a2 + (a3 + (a4 + a5 * t) * t) * t) * t) * t
    return ρ_w
end

function density0(S,T)
    t = temperature68(T)
    # page 21, equation (13) of
    # http://web.archive.org/web/20170103000527/http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf

    b0 = 8.24493e-1
    b1 = -4.0899e-3
    b2 = 7.6438e-5
    b3 = -8.2467e-7
    b4 = 5.3875e-9
    c0 = -5.72466e-3
    c1 = 1.0227e-4
    c2 = -1.6546e-6
    d0 = 4.8314e-4

    ρ = density_reference_pure_water(T) + ((b0 + (b1 + (b2 + (b3 + b4 * t) * t) * t) * t) + (c0 + (c1 + c2 * t) * t) * sqrt(S) + d0 * S) * S;
    return ρ
end
"""
    secant_bulk_modulus(S,T,p)
Compute the secant bulk modulus of sea-water (bars) at the salinity `S` (psu, PSS-78), temperature `T` (degree Celsius, ITS-90) and pressure `p` (decibar) using the UNESCO polynomial 1983.
Fofonoff, N.P.; Millard, R.C. (1983). Algorithms for computation of fundamental properties of seawater. UNESCO Technical Papers in Marine Science, No. 44. UNESCO: Paris. 53 pp.
http://web.archive.org/web/20170103000527/http://unesdoc.unesco.org/images/0005/000598/059832eb.pdf
"""
function secant_bulk_modulus(S,T,p)
    # convert decibars to bars
    p = p/10

    t = temperature68(T)

    # page 18, equation (19)
    e0 = +19652.21 # [-1930.06]
    e1 = +148.4206
    e2 = -2.327105
    e3 = +1.360477E-2
    e4 = -5.155288E-5
    Kw = e0 + (e1  + (e2 + (e3  + e4 * t) * t) * t) * t

    # page 18, equation (16)
    # probably typo f3 vs f2
    f0 = +54.6746;    g0 = +7.944E-2
    f1 = -0.603459;   g1 = +1.6483E-2
    f2 = +1.09987E-2; g2 = -5.3009E-4
    f3 = -6.1670E-5

    K0 = Kw + ((f0 + (f1 + (f2  + f3 * t) * t) * t) + (g0 + (g1 + g2 * t) * t) * sqrt(S)) * S

    if (p == 0)
        return K0
    end

    # page 19
    h0 = +3.239908 # [-0.1194975]
    h1 = +1.43713E-3
    h2 = +1.16092E-4
    h3 = -5.77905E-7
    Aw = h0 + (h1 + (h2 + h3 * t) * t) * t


    k0 = +8.50935E-5 # [+ 3.47718E-5]
    k1 = -6.12293E-6
    k2 = +5.2787E-8
    Bw = k0 + (k1 + k2 * t) * t

    # page 18, equation (17)
    i0 = +2.2838E-3; j0 = +1.91075E-4
    i1 = -1.0981E-5
    i2 = -1.6078E-6
    A = Aw + ((i0 + (i1 + i2 * t) * t) + j0 * sqrt(S)) * S

    # page 18, equation (18)
    m0 = -9.9348E-7
    m1 = +2.0816E-8
    m2 = +9.1697E-10
    B = Bw + (m0 + (m1 + m2 * t) * t) * S

    K = K0 + (A + B * p) * p

    return K
end

end

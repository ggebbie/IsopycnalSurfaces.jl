using SigmaShift
using Test, Revise

@testset "SigmaShift.jl" begin
    #############################
    # Test dedup! function
    na = 10
    for ndupes = 1:na-1
        println(ndupes)
        a = sort(randn(na))
        b = randn(na)
        # make duplicates at end
        counter = 0
        while counter < ndupes  
            a[end-counter-1] = a[end-counter]
            counter += 1
        end
        dedup!(a,b)
        println(a)
        println(b)
        println(@test issorted(a))
    end

    ################################
    # Idealized mapping onto sigma1
    pz = collect(0.:500.:4000.) # pressure levels
    nz = length(pz)

    θz = collect(range(20,stop=10,length=nz))
    Sz = collect(range(36,stop=35,length=nz))
    ztest = sort(rand(2:8,2))
    σ₁true = sigma1column(θz[ztest],Sz[ztest],pz[ztest])

    sig1grid = range(minimum(σ₁true),stop=maximum(σ₁true),length=20)
    splorder = 3
    σ₁=sigma1column(θz,Sz,pz)
    sgood = findall(minimum(σ₁) .<= sig1grid .<= maximum(σ₁))
    pσ = var2sigmacolumn(σ₁,pz,sig1grid[sgood],splorder)

    @test isapprox(pσ[begin],pz[ztest[begin]])
    @test isapprox(pσ[end],pz[ztest[end]])

    ################################

    # GibbsSeaWater
    SA = gsw_sa_from_sp.(Sz[ztest],pz[ztest],90,45 )   # Absolute S from (pracital S, pressure, lon, lat)
    CT = gsw_ct_from_pt.(SA ,θz[ztest])                # Conservative T from (Absolute S, sigma0)
    σ₁Gibbs = gsw_rho.(SA, CT,1000.) .- 1000.         



    # MOVE THIS TEST TO ECCOTOUR.JL
    # Test the mapping onto sigma1.#
    # expt = "test"
    # diagpath = pwd()
    # path_out = pwd()
    
    # γ = setupLLCgrid("grid/"))
    # nf = length(γ.fSize)

    # # get standard levels of MITgcm
    # z = depthlevels(γ)
    # pstdz = pressurelevels(z)
    # p₀ = 1000.0 ; # dbar

    # # sig1 value of interestn
    # sig1grid = 30.0;

    # TSroot = "state_3d_set1" # 1: θ, 2: S
    # splorder = 100 # spline order

    # # first filter for state_3d_set1
    # filelist = searchdir(diagpath,TSroot)

    # # second filter for "data"
    # datafilelist  = filter(x -> occursin("data",x),filelist)

    # filelist2 = searchdir(diagpath,RProot) 
    # datafilelist2  = filter(x -> occursin("data",x),filelist)

    # # make an output directory for each expteriment
    # !isdir(path_out) ? mkdir(path_out) : nothing;
    # nt = length(datafilelist)
    
    # global tt = 0
    # for datafile in datafilelist
    #     tt += 1

    #     #print timestamp
    #     year,month = timestamp_monthly_v4r4(tt)

    #     # eliminate suffix
    #     fileroot = rstrip(datafile,['.','d','a','t','a'])
    #     fileroot2 = RProot*fileroot[14:end] # a better way?
    #     fileroots = (fileroot,fileroot2)
    
    #     # Read from filelist, map to sigma-1, write to file
    #     mdsio2sigma1(diagpath,path_out,fileroots,γ,pstdz,sig1grid,splorder)
    # end
end

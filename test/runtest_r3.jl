using Revise
using IsopycnalSurfaces, Test

@testset "IsopycnalSurfaces.jl" begin
    # test input from a 3D array (failed Dec 09 2021)
    @testset "3d_array" begin

        # currently empty
        pz = collect(0.:500.:4000.) # pressure levels
        nz = length(pz)

        θz = collect(range(20,stop=10,length=nz))
        Sz = collect(range(36,stop=35,length=nz))

        ztest = sort(rand(2:8,2))
        σ₁true = sigma1column(θz[ztest],Sz[ztest],pz[ztest])

        # put θ,S,p into 3D array
        nx = 10; ny = 10;
        θ = Array{Float64,3}(undef,nx,ny,nz)
        S = Array{Float64,3}(undef,nx,ny,nz)
        p = Array{Float64,3}(undef,nx,ny,nz)
        
        [θ[i,j,k] = θz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]
        [S[i,j,k] = Sz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]
        [p[i,j,k] = pz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]

        σ₁grid = collect(range(minimum(σ₁true),stop=maximum(σ₁true),length=20))

        vars = Dict("θ" => θ, "Sp" => S)

        @testset "3d_array_spline" begin
            splorder = 3
            varsσ = vars2sigma1(vars,pz,σ₁grid,splorder)
            xx = rand(1:nx); yy = rand(1:ny)
            @test isapprox(varsσ["p"][xx,yy,begin],pz[ztest[begin]])
            @test isapprox(varsσ["p"][xx,yy,end],pz[ztest[end]])
        end

        @testset "3d_array_linear" begin
            splorder = 3
            linearinterp = true
            varsσ = vars2sigma1(vars,pz,σ₁grid,splorder,linearinterp)
            xx = rand(1:nx); yy = rand(1:ny)
            @test isapprox(varsσ["p"][xx,yy,begin],pz[ztest[begin]])
            @test isapprox(varsσ["p"][xx,yy,end],pz[ztest[end]])
        end
    end
end
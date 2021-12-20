using Revise
using IsopycnalSurfaces, Test

@testset "IsopycnalSurfaces.jl" begin
    #############################
    # Test dedup! function
    @testset "column" begin
        ################################
        # Idealized mapping onto sigma1
        pz = collect(0.:500.:4000.) # pressure levels
        nz = length(pz)

        θz = collect(range(20,stop=10,length=nz))
        Sz = collect(range(36,stop=35,length=nz))
        ztest = sort(rand(2:8,2))
        σ₁true = sigma1column(θz[ztest],Sz[ztest],pz[ztest])

        σ₁grid = range(minimum(σ₁true),stop=maximum(σ₁true),length=20)
        σ₁=sigma1column(θz,Sz,pz)

        @testset "column_spline" begin
            splorder = 3
            sgood = findall(minimum(σ₁) .<= σ₁grid .<= maximum(σ₁))
            pσ = var2sigmacolumn(σ₁,pz,σ₁grid[sgood],splorder)

            @test isapprox(pσ[begin],pz[ztest[begin]])
            @test isapprox(pσ[end],pz[ztest[end]])
        end

        @testset "column_linear" begin
            splorder = 100
            sgood = findall(minimum(σ₁) .<= σ₁grid .<= maximum(σ₁))
            pσ = var2sigmacolumn(σ₁,pz,σ₁grid[sgood],splorder)

            @test isapprox(pσ[begin],pz[ztest[begin]])
            @test isapprox(pσ[end],pz[ztest[end]])
        end

    end
end
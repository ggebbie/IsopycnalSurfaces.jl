#using Revise
using IsopycnalSurfaces, Test

@testset "IsopycnalSurfaces.jl" begin
    #############################
    # Test dedup! function
    @testset "deduplication" begin
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
    end

    @testset "point" begin

        # some examples of different interfaces.
        
        # classic definition for sigma1
        θ = Vector{Float64}(undef,1); θ[1] = 15.0
        S = Vector{Float64}(undef,1); S[1] = 35.0
        p = Vector{Float64}(undef,1); p[1] = 2000.0

        # classic EOS takes 3 classic variables.
        σ₁A = sigma1column(θ,S,p,"EOS80")
        σ₁B = sigma1column(θ,S,p,"JMD95")

        # Thermodynamic EOS will also take classic variables. 
        σ₁C = sigma1column(θ,S,p,"TEOS10")

        # Argument list can be packed into a Dict with symbols.
        vars = Dict(:θ => θ, :S => S, :p => p)
        #σ₁A2= sigmacolumn(vars,p₀=1000,eos="EOS80") # fails
        #σ₁B2= sigmacolumn(vars,p₀=1000,eos="JMD95") # fails
        σ₁C2= sigmacolumn(vars,p₀=1000,eos="TEOS10")

        # A more general version with strings instead of symbols
        vars = Dict("θ" => θ, "S" => S, "p" => p) 
        σ₁D= sigmacolumn(vars,p₀=1000,eos="TEOS10")

        # smaller than 5% error?
        @test abs(σ₁C[1] - σ₁B[1])/(abs(σ₁C[1]) + abs(σ₁B[1])) < 0.05
        @test abs(σ₁A[1] - σ₁B[1])/(abs(σ₁A[1]) + abs(σ₁B[1])) < 0.05
        @test abs(σ₁D[1] - σ₁B[1])/(abs(σ₁D[1]) + abs(σ₁B[1])) < 0.05

        # can one supply a depth index instead of pressure?
        vars = Dict("θ" => θ, "S" => S, "z" => p) 
        σ₁E= sigmacolumn(vars,p₀=1000,eos="TEOS10")
        
    end
    
        
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
            sgood = findall(minimum(σ₁) .<= σ₁grid .<= maximum(σ₁))
            pσ = var2sigmacolumn(σ₁,pz,σ₁grid[sgood],splorder=3)
            @test isapprox(pσ[begin],pz[ztest[begin]])
            @test isapprox(pσ[end],pz[ztest[end]])

        end

        @testset "column_linear" begin
            sgood = findall(minimum(σ₁) .<= σ₁grid .<= maximum(σ₁))
            pσ = var2sigmacolumn(σ₁,pz,σ₁grid[sgood],linearinterp=true)

            @test isapprox(pσ[begin],pz[ztest[begin]])
            @test isapprox(pσ[end],pz[ztest[end]])
        end

        @testset "EOS" begin
            # just check that this runs
            # would be better to show some similarity in results
            varscol = Dict(:θ => θz, :S => Sz, :p => pz) 
            namescol = parsevars(varscol)

            println("Using the EOS from the Thermodynamic Equation of State 2010")
            display(sigmacolumn(varscol,p₀=2000,eos="TEOS10"))
            println("Using the EOS from MITgcmTools.jl:")
            display(sigma2column(θz,Sz,pz,"JMD95")) # using the EOS from MITgcmTools.jl
            println("")
            println("Default, using the EOS from PhysOcean.jl/EOS80.jl:")
            display(sigma2column(θz,Sz,pz))  # default, using the EOS from PhysOcean.jl/EOS80.jl

            zz = rand(1:nz)
            
            # EOS80 versus JMD95
            @test abs(sigma2column(θz,Sz,pz,"JMD95")[zz] - sigma2column(θz,Sz,pz,"EOS80")[zz])/(sigma2column(θz,Sz,pz,"JMD95")[zz] + sigma2column(θz,Sz,pz,"EOS80")[zz]) < 0.005

            # EOS80 versus TEOS10
            @test abs(sigmacolumn(varscol,p₀=2000,eos="TEOS10")[zz] - sigma2column(θz,Sz,pz,"EOS80")[zz])/(sigmacolumn(varscol,p₀=2000,eos="TEOS10")[zz] + sigma2column(θz,Sz,pz,"EOS80")[zz]) < 0.005
            
        end # EOS

        # test input from a 3D array
        @testset "3d_array" begin
            
            eoses = ("EOS80","JMD95")
            for eos in eoses
                ztest = sort(rand(2:8,2))
                σ₁true = sigma1column(θz[ztest],Sz[ztest],pz[ztest],eos)

                # put θ,S,p into 3D array
                nx = 10; ny = 10;
                θ = Array{Float64,3}(undef,nx,ny,nz)
                S = Array{Float64,3}(undef,nx,ny,nz)
                p = Array{Float64,3}(undef,nx,ny,nz)
                
                [θ[i,j,k] = θz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]
                [S[i,j,k] = Sz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]
                [p[i,j,k] = pz[k] for i = 1:nx for j = 1:ny for k = 1:nz ]

                σ₁grid = collect(range(minimum(σ₁true),stop=maximum(σ₁true),length=20))

                vars_strings = Dict("θ" => θ, "Sₚ" => S)
                vars_symbols = Dict(:θ => θ, :Sₚ => S, :p => p)

                # pick a random location to sample
                xx = rand(1:nx); yy = rand(1:ny)

                @testset "spline" begin

                    # test when input dictionary uses strings
                    varsσ = vars2sigma1(vars_strings,σ₁grid,p=pz,splorder=3,eos=eos)
                    @test isapprox(varsσ["p"][xx,yy,begin],pz[ztest[begin]])
                    @test isapprox(varsσ["p"][xx,yy,end],pz[ztest[end]])

                    # test when input dictionary uses symbols
                    varsσ = vars2sigma1(vars_symbols,σ₁grid,p=pz,splorder=3,eos=eos)
                    @test isapprox(varsσ[:p][xx,yy,begin],pz[ztest[begin]])
                    @test isapprox(varsσ[:p][xx,yy,end],pz[ztest[end]])

                end

                @testset "linear" begin
                    #splorder = 3
                    varsσ = vars2sigma1(vars_symbols,σ₁grid,p=pz,linearinterp=true,eos=eos)
                    xx = rand(1:nx); yy = rand(1:ny)
                    @test isapprox(varsσ[:p][xx,yy,begin],pz[ztest[begin]])
                    @test isapprox(varsσ[:p][xx,yy,end],pz[ztest[end]])

                end

            end
        end
    end # column
end

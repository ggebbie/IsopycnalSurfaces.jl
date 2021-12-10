using Revise
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
end
using Revise
using IsopycnalSurfaces, Test

pz = collect(0.:500.:4000.) # pressure levels
nz = length(pz)
θz = collect(range(20,stop=10,length=nz))
Sz = collect(range(36,stop=35,length=nz))

#println("Default, using the EOS from PhysOcean.jl/EOS80.jl:")
#display(sigma2column(θz,Sz,pz))  # default, using the EOS from PhysOcean.jl/EOS80.jl

println("")
println("Using the EOS from MITgcmTools.jl:")
display(sigma2column(θz,Sz,pz,"pressure","JMD95")) # using the EOS from MITgcmTools.jl

println("")
println("Using the EOS from GibbsSeaWater.jl:")
display(sigma2column(θz,Sz,pz,"pressure","Gibbs")) # using the EOS from GibbsSeaWater.jl

println("")
println("Using the EOS from PhysOcean.jl.jl:")
display(sigma2column(θz,Sz,pz,"pressure","EOS80")) # using the EOS from PhysOcean.jl

println("")
display((sigma2column(θz,Sz,pz,"pressure","Gibbs") .- sigma2column(θz,Sz,pz,"pressure","JMD95")) ./ sigma2column(θz,Sz,pz,"pressure","JMD95"))
# realtive difference between Gibbs and JMD95 is around 0.015~0.02%
display((sigma2column(θz,Sz,pz,"pressure","EOS80") .- sigma2column(θz,Sz,pz,"pressure","JMD95")) ./ sigma2column(θz,Sz,pz,"pressure","JMD95"))
# realtive difference between EOS80 and JMD95 is around 0.15~0.3%



println("")
println("Default: pressure, GibbsSeaWater")
display(sigmacolumn(θz,Sz,pz,2000))

using Revise
using IsopycnalSurfaces, Test

pz = collect(0.:500.:4000.) # pressure levels
nz = length(pz)
θz = collect(range(20,stop=10,length=nz))
Sz = collect(range(36,stop=35,length=nz))


println("")
display(sigmacolumn(θz,Sz,pz,2000,"pressure","JMD95"))

println("")
display(sigmacolumn(θz,Sz,pz,2000,"depth","JMD95"))

println("")
println("Default: pressure coordinate and GibbsSeaWater")
display(sigma2column(θz,Sz,pz))

println("")
println("")
println("Test EOS in MITgcmTools.jl, depth or pressure coordinate affects the results of in-situ density:")
display(SeaWaterDensity(3.,35.5,505., 2000.) .- 1000. )
display(SeaWaterDensity(3.,35.5,500., 2000.) .- 1000. )

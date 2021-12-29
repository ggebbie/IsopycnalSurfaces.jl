using Revise
using IsopycnalSurfaces, Test

pz = collect(0.:500.:4000.) # pressure levels
nz = length(pz)
θz = collect(range(20,stop=10,length=nz))
Sz = collect(range(36,stop=35,length=nz))


println("Using the EOS from MITgcmTools.jl:")
display(sigma2column(θz,Sz,pz,"JMD95")) # using the EOS from MITgcmTools.jl
println("")
println("Default, using the EOS from PhysOcean.jl/EOS80.jl:")
display(sigma2column(θz,Sz,pz))  # default, using the EOS from PhysOcean.jl/EOS80.jl


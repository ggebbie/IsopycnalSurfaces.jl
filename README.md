# SigmaShift

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/SigmaShift.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/SigmaShift.jl/dev)
[![Build Status](https://github.com/ggebbie/SigmaShift.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/SigmaShift.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/SigmaShift.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/SigmaShift.jl)

## Installation

pkg> add SigmaShift


## Usage (load module)

include("./SigmaShift.jl/src/SigmaShift.jl")

## Test1
pz = collect(0.:500.:4000.) # pressure levels
nz = length(pz)

θz = collect(range(20,stop=10,length=nz))
Sz = collect(range(36,stop=35,length=nz))
ztest = sort(rand(2:8,2))
σ₁true = SigmaShift.sigma1column(θz[ztest],Sz[ztest],pz[ztest])

Results: 
	30.041106141810587
    30.906751163764284

sig1grid = range(minimum(σ₁true),stop=maximum(σ₁true),length=20)
splorder = 3
σ₁=SigmaShift.sigma1column(θz,Sz,pz)
sgood = findall(minimum(σ₁) .<= sig1grid .<= maximum(σ₁))
pσ = SigmaShift.var2sigmacolumn(σ₁,pz,sig1grid[sgood],splorder)

Results:
  500.0
  594.1585568815079
  689.2920330984479
  785.4352311816727
  882.6229536620141
  980.8900030703251
 1080.2727053085036
 1180.8184014302499
 1282.5792758100033
 1385.6075337303162
 1489.9553804737648
 1595.6775607975524
 1702.842371037887
 1811.5226095980352
 1921.7910778786677
 2033.7207448928955
 2147.397582546121
 2262.929623162734
 2380.427038421197
 2499.9999999999995
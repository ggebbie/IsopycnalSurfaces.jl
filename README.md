# IsopycnalSurfaces.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ggebbie.github.io/IsopycnalSurfaces.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ggebbie.github.io/IsopycnalSurfaces.jl/dev)
[![Build Status](https://github.com/ggebbie/IsopycnalSurfaces.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ggebbie/IsopycnalSurfaces.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ggebbie/IsopycnalSurfaces.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ggebbie/IsopycnalSurfaces.jl)

* IsopycnalSurfaces.jl was originally started by G Jake Gebbie as a part of the ECCOtour.jl package (https://github.com/ggebbie/ECCOtour.jl)

* Core functionality: Take a profile (water column) of temperature and salinity on pressure or depth coordinates and transfer onto a vertical coordinate of density. 

* Other goals: Permit other properties besides temperature and salinity to be mapped onto depth. Allow 3D or 4D input fields. Allow various equations of state (not implemented yet). 

* This Julia package is in early development and breaking changes are expected.

* See the function list in the documentation linked through the badge above

* After setting up the Julia environment (instructions below), check that all tests pass via the following shell command in the repository base directory:
`julia --project=@. test/runtests.jl`

* This project was originally called *SigmaShift.jl*.

* Also see Greg Wagner's package Isopycnal.jl (https://github.com/glwagner/Isopycnal.jl).

# Requirements

Compatibility with the built-in tests requires Julia 1.6+. 

# Installation (from jrayshi)
pkg> add IsopycnalSurfaces (not working)

Optional:
cd ./IsopycnalSurfaces.jl\
pkg> activate .

# Setting up project environment

Details about setting up a Julia environment are available at https://github.com/ggebbie/ECCOtour.jl#readme .

# Data files

The Julia code is designed to download input files from Google Drive (using GoogleDrive.jl) and to place them in the `data` directory. 

# Functions

Available functions are listed in the documentation at https://ggebbie.github.io/IsopycnalSurfaces.jl/dev/ .

# How this Julia package was started

This package was generated using PkgTemplates.jl. 

Steps: 
1. Use PkgTemplates to make git repo.

Run the following Julia code

`using Revise, PkgTemplates`

`t = Template(; 
    user="ggebbie",
    dir="~/projects",
    authors="G Jake Gebbie",
    julia=v"1.6",
    plugins=[
        License(; name="MIT"),
        Git(; manifest=true, ssh=true),
        GitHubActions(; x86=false),
        Codecov(),
        Documenter{GitHubActions}(),
        Develop(),
    ],
             )`

`t("IsopycnalSurfaces.jl")`

2. Make a new empty repository on GitHub.\
	
3. Then push this existing repository from the command line:
    `git remote add origin git@github.com:ggebbie/IsopycnalSurfaces.jl.git`\
    `git branch -M main`\
    `git push -u origin main`

4. Use Documenter.jl and DocumenterTools to automatically deploy documentation following: https://m3g.github.io/JuliaNotes.jl/stable/publish_docs/ .









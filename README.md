[![Open in Visual Studio Code](https://classroom.github.com/assets/open-in-vscode-718a45dd9cf7e7f842a935f5ebbe5719a5e09af4491e668f4dbf3b35d5cca122.svg)](https://classroom.github.com/online_ide?assignment_repo_id=11858194&assignment_repo_type=AssignmentRepo)
# Astro 528 [Class Project](https://psuastro528.github.io/Fall2023/project/)

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://PsuAstro528.github.io/project-template/stable)

GitHub Actions : [![Build Status](https://github.com/PsuAstro528/project-template/workflows/CI/badge.svg)](https://github.com/PsuAstro528/project-template/actions?query=workflow%3ACI+branch%3Amain)


##  Overview
In this project, I begin by simulating the March 20, 2015 solar eclipse from Gottingen, Germany via the sun's observed velocity with respect to the observer. I query JPL horizons using SPICE for required body parameters. The surface of the sun is gridded to allow for determinaiton of which cells are blocked by the moon during the eclipse and more accurate solar velocities (astrophysical reasons lie with convective blueshift from granulation). The limb darkening weighted velocity is determined for each visual cell on the solar surface grid along the line of sight to the observer. Once my project was confirmed to be working using the 2015 eclipse as a benchmark via Reiners 2015, I updated the code to simulate the October 14th 2023 eclipse from NEID and EXPRES. 
To collect (1) RM curve (2) change in relative intensity (3) eclipse movie run:
using MyProject
MyProject.get_kernels()
MyProject.low_loop(num_lats, num_lons) for express
MyProject.kitt_loop(num_lats, num_lons) for neid
MyProject.gottingen_loop(num_lats, num_lons) for gottingen
(plots done in python under plots folder)

## Parallel Code V1
This project will parallelize over the grid size to evaluate performance and accuracy of recovered velocity with respect to grid size. For now my parallel code is over multiple cores using a shared memory system making use of Threads.@threads, ThreadsX.collect, and ThreadsX.map. My most inner functions in coordinates.jl, moon.jl, and velocity.jl have been duplicated to have a parallel version. The main script that computes the RVs is the compute_rv function in epoch_computations.jl - this has been duplicated in epoch_computations_pa. So compute_rv in epoch_computations.jl is my serial code and compute_rv_pa in epoch_computations_pa.jl is my parallel code. I determine the time taken for each one over a range of grid sizes (200-375). To get this run parallel_loop() and serial_loop(). The results are found in gridvstime.png and as you can see the parallel code has a shorter compute time than the serial code as the grid size increases YAY.  

## Parallel Code V2
Version one of my parallel code used multi-threading and is described above. For version two, I use multi-processing. This parallel code can be found in parallel_v2.jl - where I have rewritten the computation compute_rv for a single patch so I can distribute individual patches to the processors where the computation is then done; this is done in compute_rv_pa2. Then the function parallel_v2 computes the grid for a given resolution and distributes individual patches so the components of the velocity can be done on separate processors. For this parallel code I use a single timestamp during the Gottingen 2015 eclipse and I have confirmed that the restructured computation returns the expected velocity. For this submission, I have determine strong scaling for two problem sizes: 50x100 and 250x500 grid and have also determined the weak scaling for 1-12 processors. To run, the following is done in the julia terminal: 
using Distributed
addprocs(12)
@everywhere using MyProject
@everywhere MyProject.get_kernels()
MyProject.parallel_v2()
This will save the run time for each case into a jld2 file which I then plot using python. Figure scaling_v2.png under test shows the strong scaling for both problem sizes and weak scaling. 
*For final project submission, I will both increases the number of processors and benchmark additional problem sizes. I will also benchmark the time loop throughout all timestamps using different number of processors. 


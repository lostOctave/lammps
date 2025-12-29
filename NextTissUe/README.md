![image](https://github.com/user-attachments/assets/c28c14ec-5e6b-44fc-85a4-b1c8a65ab0fa)


# NextTissUe
A new approach to modeling epithelial cells and tissues, by A. Pasupalak, W. Zu and M. Pica Ciamarra. 
See [here](https://arxiv.org/abs/2409.16128) for an introduction.

We provide NextTissUe as an additional package to Lammps, a Large-scale Atomic/Molecular Massively Parallel Simulator, and thank its developers for their work.
Basic knowledge of Lammps is recommended to run NextTissUe.

Follow the steps below to install LAMMPS and NextTissUe.

1) Download and install Lammps. Our code has been built and tested on a Linux architecture and version 2 of Lammps Aug-2023. Link to download https://www.lammps.org/download.html
2) Extract the zip file USER-AREABOND.zip into your Lammps src folder.
3) Modify the file 'YOUR_PATH_TO_LAMMPS/lammps-2Aug2023/cmake/CMakeLists.txt' by adding include 'USER-AREABOND' in "set(STANDARD_PACKAGES ....)".
5) Compile Lammps using packages molecule, extra-pair, and USER-AREABOND, following the instructions from https://docs.lammps.org/Build_cmake.html

# How to run NextTissuE
We provide an example Lammps input script that sets up an initial configuration and evolves the system, "in.demo".
Follow the instructions at https://docs.lammps.org/Run_basics.html to run it.

The script includes comments to describe the different performed steps.
It relies on the executables "template.x" to create a template ring-polymer model (a molecule template, in Lammps), and "random_propelling_directions.x" to initialise the self-propelling directions.

To proceed:

1) Download all the files.
2) Compile "template.cpp" into an executable with the name "template.x".
4) Compile "random_propelling_directions.cpp" into an executable with the name "random_propelling_directions.x".
5) Run the Lammps executable file on the "in.demo" input script. 


# charge_transport_package
The charge transport package assists to perform charge transport simulations in organic solar cells and 
amorphous donor-acceptor blend materials. The goal here is to compute the mobility of electrons and holes 
in a multiscale model via kinetic Monte Carlo simulations using bimolecular charge transfer rates based on
the Marcus theory using Marcus, Jortner or Weiss-Dorsey rates.  
Therefore, accurate hopping rates between neighbouring molecules are mandatory, depending on the internal and 
outer-sphere reorganization energies, site-energy differences, charge transfer integrals, and external driving forces.  
One can combine intermolecular and intramolecular charge transfer in amorphous materials.  
The hopping rates are then fed into the kinetic Monte Carlo (kMC) simulations, which allow to monitor the charge dynamics 
in the system as well as to get ensemble averages of the occupation probability, charge carrier mobility, 
electronic fluxes and pathways of minimal local resistance. Hence, the multiscale model links the microscopic structure 
to macroscopic observables.


The kMC simulations can be performed using three different programms.
A) ctp_run / kmc_run (VOTCA, C++, fast),
B) kmc_multigeo_V4 (Python Jim B., on-the-fly modus)
C) charge_transport_package, (Fortran Tobias K., a priori modus)

The options for different calculations are set using
A) options_VOTCA.xml
B) options are set at the beginning of the kmc_multigeo_V4 file
C) options_charge_transport_package.xml

run_kmc_votca_multiple_V3.sh  helps to handle multiple kMC simulations.
The simulation data for A) and C) is stored in state.sql files.


A) build_order_VOTCA.sh provides the basic workflow in VOTCA.

C) multi_morphologies_VOTCA_prepare_kMC.sh helps to prepare the system and sample multiple frames.
   It includes the initial steps for a single morphology.gro file from GROMACS MM simulations. 
   
C) xtp_createStatefile_fromtxt_t_koch08.py is called by the charge_transport_package and helps to generate state.sql files.
run_kmc_votca_multiple_V3.sh  helps to handle multiple kMC simulations.

A short summary how to start each calculator is provided in kMC_Anleitung_Kurzfassung.pdf.
Note that many programmes expect the folder structure given here. 

 =================================
 |   charge transport package    |
 =================================
If you want to use charge_transport_package for DIPRO calculations you need to 
1) Delete the ! in before !CALL DSYEV in  charge_transport_package.f90 
and compile with LAPACK and BLAS package:
2)  gfortran charge_transport_package.f90  -o charge_transport_package    /usr/lib/lapack/liblapack.so.3  /usr/lib/libblas/libblas.so.3 

 charge_transport_package
 Possible calculators: 
 reduce_kmc_xyz
 calc_pathway
 Dijkstra_kmcNetwork
 Dijkstra_newNetwork
 Dijkstra_loadNetwork
 KMC_FOR_MULTIPLE_CHARGES
 KMC_FOR_MULTIPLE_CHARGES_xml      use votca optionsfile for kmc run options
 gro_to_VOTCA
 rates_to_VOTCA_sqlfile
 calc_current_I
 calc_current_I_from_sql
 calc_current_I_xml
 DIPRO
 DIPRO_pair
 DIPRO_xyz
 DIPRO_nMO
 lambda_out_dE_out
 lambda_out
 lambda_in
 lambda_in_oniom
 rate_calculator
 data_analysis



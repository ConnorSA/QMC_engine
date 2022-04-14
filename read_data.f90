!> @brief Subroutines to read files and transfer user inputs into fortran
MODULE read_data

  USE shared_constants

  IMPLICIT NONE

  CONTAINS

    !> @brief Read in data from a text file
    !>
    !> This module reads txt, csv files line by line.
    !> The implementation is basic but will read each line properly
    !> as it was written specifically for the parameters needed from the user.
    SUBROUTINE readtxt(filename, n_electrons_in, n_atoms_in, bond_length, n_trials, &
       n_MCMC_steps, n_omp_threads, use_biased_optimiser, n_optimiser_steps, &
       basis_type_in, proton_number_int_in, &
       search_seed, n_basis_functions_per_atom_in, n_Jastrow_in, &
       fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in, plot_distance, &
       plot_points)

       ! User Inputs required - the user must define all of these
       ! Will obtain these from either txt file or a gui that generates the txt file

       !> Filename of input text file
       CHARACTER(LEN=20) :: filename
       !> Number of electrons
       integer :: n_electrons_in
       !> Number of atoms
       integer :: n_atoms_in
       !> Bond length for 2 atoms
       real(dp) :: bond_length
       !> number of trials in the hypercube search ( for each step of the biased optimiser)
       integer :: n_trials
       !> Total number of steps for each MCMC run
       integer :: n_MCMC_steps
       !> Number of OMP Threads
       integer :: n_omp_threads
        !> Flag for use of biased optimiser: true uses biased optimiser, false just uses hypercube search
       logical :: use_biased_optimiser
       !> Number of steps in the optimser
       integer :: n_optimiser_steps
       ! Optional user Inputs - these may be defined by the user, but have good default values below
       !> Number of linear terms in the single-electron wavefunction per atom
       integer :: n_basis_functions_per_atom_in
       !> Seed for the latin hypercube
       integer :: search_seed
        !> Number of dofs in the Jastrow (interaction) term
       integer :: n_Jastrow_in
       !> Lengthscale of the finite difference code
       real(dp) :: fd_length_in
       !> Inverse lengthscale of nuclear-electron interaction
       real(dp) :: Jastrow_b_length_in
       !> Inverse lengthscale of electron-electron interaction
       real(dp) :: Jastrow_d_length_in
       !> Proton number of atoms
       integer :: proton_number_int_in
       !> integer code for type of basis. Codes defined in shared_constants.f90
       integer :: basis_type_in

       !> Amount of distance away from atoms to plot
       real(dp) :: plot_distance
       !> Number of points along each axis to plot
       integer :: plot_points

       !> Integer value to determine if biased optimiser should be used
       integer :: bias_opt_0_1





      OPEN(UNIT = 7, FILE = filename)

      ! Read in lines containing only one value
      READ(7,*) n_electrons_in
      READ(7,*) n_atoms_in
      READ(7,*) bond_length
      READ(7,*) n_trials
      READ(7,*) n_MCMC_steps
      READ(7,*) bias_opt_0_1
      READ(7,*) n_optimiser_steps
      READ(7,*) basis_type_in
      READ(7,*) proton_number_int_in
      READ(7,*) n_omp_threads
      READ(7,*) search_seed
      READ(7,*) n_basis_functions_per_atom_in
      READ(7,*) n_Jastrow_in
      READ(7,*) fd_length_in
      READ(7,*) Jastrow_b_length_in
      READ(7,*) Jastrow_d_length_in
      READ(7,*) plot_distance
      READ(7,*) plot_points

      CLOSE(7)

      if (bias_opt_0_1 == 1) then
        use_biased_optimiser = .true.
      else if ( bias_opt_0_1 == 0 ) then
        use_biased_optimiser = .false.
      end if


    ! Print to see if data was properly read into fortran
      ! PRINT *, n_electrons_in
      ! PRINT *, n_atoms_in
      ! PRINT *, bond_length
      ! PRINT *, n_trials
      ! PRINT *, n_MCMC_steps
      ! PRINT *, n_basis_functions_per_atom_in
      ! PRINT *, search_seed
      ! PRINT *, n_Jastrow_in
      ! PRINT *, fd_length_in
      ! PRINT *, Jastrow_b_length_in
      ! PRINT *, Jastrow_d_length_in

    END SUBROUTINE readtxt

END MODULE read_data

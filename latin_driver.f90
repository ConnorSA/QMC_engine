!> @brief Main driver for simulation using the latin_hypercube search
!!
!! This code computes the ground state energy and corresponding trial wavefunction
!! for a given atom and electron configuration.
!! It saves the resulting wavefunction or electron density in netcdf.
!!
program main_driver
  use shared_constants
  use read_data
  use basis_functions
  use param_search
  use mcmc
  use calculations
  use write_file
  !$ use omp_lib
  implicit none

  !> Input text filename
  CHARACTER(LEN=20) :: input_filename = "init_params.txt"
  ! These parameters are set directly from the input script
  !> Number of electrons
  integer :: n_electrons_in
  !> Number of atoms
  integer :: n_atoms_in
  !> Bond length for 2 atoms
  real(dp) :: bond_length
  !> number of trials in the hypercube search
  integer :: n_trials
   !> Total number of steps for each MCMC run
  integer :: n_MCMC_steps
   !> Flag for use of biased optimiser: true uses biased optimiser, false just uses hypercube search
  logical :: use_biased_optimiser
  !> Number of steps in the optimser
  integer :: n_optimiser_steps
  ! Optional user Inputs - these may be defined by the user, but have good default values
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
  !>integer code for type of basis. In shared_constants.f90
  integer :: basis_type_in

  ! Following assigned based on input script parameters
  !> Coordinates of the atoms. Shape (3,n_atoms)
  real(dp), allocatable, dimension(:,:) :: atom_coords_in
  real(dp), allocatable ,dimension(:):: proton_numbers_in
  ! MCMC variables
  !> MCMC initial coordinates
  real(dp), allocatable, dimension(:) :: x_0
  integer ::  n_burned, thinning_interval
  integer ::  e_code
  real(dp) :: s, E, norm_coeff
  real(dp), allocatable, dimension(:,:) :: mcmc_run

  !MCMC Adaption parameters
  integer :: mcmc_adapt_steps
  real(dp) :: mcmc_adapt_s_0
  real(dp) :: mcmc_adapt_s_max
  real(dp) :: mcmc_adapt_s_min
  real(dp) :: mcmc_adapt_memory
  integer :: mcmc_adapt_interval

  ! Internal variables
  real(dp), dimension(:, :), allocatable :: trials ! array of trials
  integer :: i, loop ! Loop variable
  real(dp) :: cpu_start_time, cpu_finish_time ! CPU timing
  integer :: real_start_time,real_finish_time, rate ! Real timing
  integer :: n_omp_threads
  ! Output Variables
  integer :: n_save_points
  real(dp), dimension(:), allocatable :: trial_energies
  real(dp) :: box
  real(dp), dimension(:), allocatable :: density_result
  integer :: n_density_integration_points
  ! Start timers
  call cpu_time(cpu_start_time)
  call system_clock(real_start_time,rate)

  ! Read in text file to get user inputs
     call readtxt(input_filename, n_electrons_in, n_atoms_in, bond_length, n_trials, &
     n_MCMC_steps, n_omp_threads, use_biased_optimiser, n_optimiser_steps, &
     basis_type_in, proton_number_int_in, &
     search_seed, n_basis_functions_per_atom_in, n_Jastrow_in, &
     fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in, box, &
     n_save_points)
  ! Uses this to test that inputs were read in properly
  if (.true.) then
    print *, "testing input parameters:"
    PRINT *, "n_electrons_in=",n_electrons_in
    PRINT *, "n_atoms_in=",n_atoms_in
    PRINT *, "bond_length=",bond_length
    PRINT *, "n_trials=",n_trials
    PRINT *, "n_MCMC_steps=",n_MCMC_steps
    PRINT *, "n_omp_threads=",n_omp_threads
    print*, "use_biased_optimiser=",use_biased_optimiser
    print*, "n_optimiser_steps=",n_optimiser_steps
    print*, "basis_type_in=",basis_type_in
    print*, "proton_number_int_in=",proton_number_int_in
    PRINT *, "search_seed=",search_seed
    PRINT *, "n_basis_functions_per_atom_in=",n_basis_functions_per_atom_in
    PRINT *, "n_Jastrow_in=",n_Jastrow_in
    PRINT *, "fd_length_in=",fd_length_in
    PRINT *, "Jastrow_b_length_in=",Jastrow_b_length_in
    PRINT *, "Jastrow_d_length_in=",Jastrow_d_length_in
    print*, "box=",box
    print*, "n_save_points=",n_save_points
  end if
  !As a debug, can set parameters manually
  if (.false.) then
  n_electrons_in = 2
  n_atoms_in =2
   bond_length = 1.5
    n_trials = 50
     n_MCMC_steps = 100000
     search_seed = 200
      n_basis_functions_per_atom_in = 1
        n_Jastrow_in = 7
     fd_length_in = 0.01
      Jastrow_b_length_in = 1.0
       Jastrow_d_length_in = 1.0
       use_biased_optimiser = .false.
  n_optimiser_steps = 10
  basis_type_in = slater_1s_code
  n_save_points = 20
  box = 5.0
  n_omp_threads = OMP_GET_MAX_THREADS() !Change this if you want less threads
  end if

  ! Biased optimiser isn't yet implemented
  if (use_biased_optimiser) then
    print*, "Sorry, Biased Optimiser not implemented yet"
  end if

  ! Additional parameters
  ! MCMC adapt parameters
  mcmc_adapt_steps = 10000
  mcmc_adapt_s_0 = 0.1_dp
  mcmc_adapt_s_max = 0.5_dp
  mcmc_adapt_s_min = 0.01_dp
  mcmc_adapt_memory = 500.0_dp
  mcmc_adapt_interval = 100

  ! Allocate atom positions and proton numbers
  ! If this is changed the the box in calc.f90 will need changing
  allocate( atom_coords_in(3,n_atoms_in) )
  if (n_atoms_in == 1) then
   atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp] ! 1 atom at origin
  else if (n_atoms_in == 2) then
   atom_coords_in(:,1)=[bond_length/2.0_dp,0.0_dp,0.0_dp] !2 atoms evenly placed about origin
   atom_coords_in(:,2)=[-bond_length/2.0_dp,0.0_dp,0.0_dp]
  end if
  allocate( proton_numbers_in(n_atoms_in) )
  proton_numbers_in = real(proton_number_int_in,dp)

  ! Initialisation of basis and trials
  call omp_set_num_threads(n_omp_threads)
  call initialise_basis(n_electrons_in, n_basis_functions_per_atom_in, n_atoms_in, atom_coords_in,&
                        n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in,&
                        proton_numbers_in,basis_type_in)

  ! Initialise trials
  allocate(trials(n_trials, number_dofs))
  call  latin_hypercube(trials,dof_bounds, n_trials, search_seed)

  allocate(trial_energies(size(trials,1)))

  ! MCMC setup
  n_burned = n_MCMC_steps/10 ! Integer division
  thinning_interval = n_MCMC_steps/10000
  allocate(x_0(n_space_dims))
  x_0(1:3) = atom_coords_in(:,1) + 0.2_dp ! Displace starting electrons from atoms
  if (n_electrons_in == 2) then
    x_0(4:6) = atom_coords_in(:,1) - 0.2_dp
  end if
  norm_coeff = real(thinning_interval,dp)/real(n_MCMC_steps-n_burned,dp)
  allocate(mcmc_run((n_MCMC_steps-n_burned)/thinning_interval+1,n_space_dims))

  ! Main loop over trials in parameter space!ordered
  !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
  do i=1, size(trials,1)
    ! Run MCMC to generate coordinates in space
    call mcmc_adapt(s, log_density, x_0, mcmc_adapt_steps, mcmc_adapt_s_0, e_code, mcmc_adapt_s_max,&
    mcmc_adapt_s_min, mcmc_adapt_memory, mcmc_adapt_interval, trials(i,:),n_space_dims)
    call mcmc_sample(mcmc_run, log_density, x_0, n_MCMC_steps, n_burned, thinning_interval, s,e_code,&
     trials(i,:),n_space_dims)

    ! Compute energy by summing over positions
    E=0.0_dp
    do loop=1,(n_MCMC_steps-n_burned)/thinning_interval
      E = E + reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
    end do

    trial_energies(i) = norm_coeff*E

    print*, "completed trial" ,i," out of ", size(trials,1)
    print*, "energy = ",trial_energies(i)
  end do


  call find_best_params(trial_energies,trials)
  call param_wall_check(dof_bounds)

  print*, "Search Complete"

  print*,"best trial=", best_trial
  print*,"min energy=", trial_energies(best_trial),"Hartree Energy"
  print*,"min energy=", electronvolt*trial_energies(best_trial),"eV"
  print*,"best parameters=", best_params

  print *, "Computing output density"
  n_density_integration_points = 10000
  density_result  = calc(best_params, n_save_points, box, n_electrons_in , n_density_integration_points)
  print*, "Computed output density"
  call result_netcdf(best_params, density_result, n_electrons_in , n_atoms_in)
  call xy_grid(n_save_points, box)





  ! deallocation
  deallocate(trials,trial_energies,atom_coords_in,x_0,mcmc_run,proton_numbers_in,density_result)
  call deinitialise_basis

  ! Report Timings
  call cpu_time(cpu_finish_time)
  call system_clock(real_finish_time,rate)
  print '("Total cpu runtime = ",f7.3," seconds.")',cpu_finish_time-cpu_start_time
  print'("Total real runtime = ",f7.3," seconds.")',real(real_finish_time-real_start_time,kind=dp)/real(rate,kind=dp)

end program main_driver

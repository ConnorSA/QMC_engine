!> @brief Main driver for simulation
!!
!! More description to be added
!!
!> @param n_electrons_in Number of electrons
program main_driver
  use read_data
  use basis_functions
  use param_search
  use mcmc
  use Biased_Optim

  !use calculations
  !use write_file
  !$ use omp_lib
  implicit none

  ! User Inputs required - the user must define all of these
  ! Will obtain these from either txt file or a gui that generates the txt file
  !> Filename of input text file
  CHARACTER(LEN=20) :: input_filename = "init_params.txt"
  !> Number of electrons
  integer :: n_electrons_in
  !> Number of atoms
  integer :: n_atoms_in
  !> Bond length for 2 atoms
  real(dp) :: bond_length
  !> Coordinates of the atoms. Shape (3,n_atoms)
  real(dp), allocatable, dimension(:,:) :: atom_coords_in
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
  !> Proton numbers of atoms
  real(dp), allocatable ,dimension(:):: proton_numbers_in
  !> integer code for type of basis. Codes defined in shared_constants.f90
  integer :: basis_type_in

  !> Amount of distance away from atoms to plot
  real(dp) :: plot_distance
  !> Number of points along each axis to plot
  integer :: plot_points

  ! MCMC variables
  !> MCMC initial coordinates
  real(dp), allocatable, dimension(:) :: x_0
  integer ::  n_burned, thinning_interval
  integer ::  e_code
  real(dp) :: s, E, norm_coeff
  real(dp), allocatable, dimension(:,:) :: mcmc_run

  !Optimiser variables
  real(dp) :: ker_var_in, ker_lengthscale_in, constant_mean_prior_in, optim_rate_para_in
  real(dp), dimension(:, :), allocatable :: trials_temp
  integer :: j
  integer :: n_seed
  integer, dimension(:), allocatable:: seed
  real(dp), dimension(:), allocatable :: min_gp_params
  real(dp) :: min_gp_E
  !MCMC Adaption parameters
  integer :: mcmc_adapt_steps
  real(dp) :: mcmc_adapt_s_0
  real(dp) :: mcmc_adapt_s_max
  real(dp) :: mcmc_adapt_s_min
  real(dp) :: mcmc_adapt_memory
  integer :: mcmc_adapt_interval
  ! Outputs

  real(dp), dimension(:), allocatable :: trial_energies
  real(dp), dimension(:), allocatable :: output_best_parameters
  ! Internal variables
  real(dp), dimension(:, :), allocatable :: trials ! array of trials
  integer :: i, loop ! Loop variables
  real(dp) :: cpu_start_time, cpu_finish_time ! CPU timing
  integer :: real_start_time,real_finish_time, rate ! Real timing

  ! Calculations and NetCDF
  ! 1 atom variables
  character(LEN=30) :: filename_result = 'atom_1_result.nc4'
  real(dp), dimension(:), allocatable :: wave_1d, wave_2d, wave_3d

  ! 2 atom variables
  real(dp), dimension(:), allocatable :: ele_2atom, wave_2atom
  character(LEN=30) :: filename_result_2 = 'atom_2_result.nc4'

  !> Start timers
  call cpu_time(cpu_start_time)
  call system_clock(real_start_time,rate)

  ! Read in text file to get user inputs
  call readtxt(filename, n_electrons_in, n_atoms_in, bond_length, n_trials, &
     n_MCMC_steps, n_omp_threads, use_biased_optimiser, n_optimiser_steps, &
     basis_type_in, proton_numbers_in, &
     search_seed, n_basis_functions_per_atom_in, n_Jastrow_in, &
     fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in, plot_distance, &
     plot_points)



  ! Uses this to test that inputs were read in properly
     if (.true.) then
      print *, "testing input parameters:"
  PRINT *, "n_electrons_in=",n_electrons_in
  PRINT *, "n_atoms_in=",n_atoms_in
  PRINT *, "bond_length=",bond_length
  PRINT *, "n_trials=",n_trials
  PRINT *, "n_MCMC_steps=",n_MCMC_steps
  PRINT *, "search_seed=",search_seed
  PRINT *, "n_basis_functions_per_atom_in=",n_basis_functions_per_atom_in
  PRINT *, "n_Jastrow_in=",n_Jastrow_in
  PRINT *, "fd_length_in=",fd_length_in
  PRINT *, "Jastrow_b_length_in=",Jastrow_b_length_in
  PRINT *, "Jastrow_d_length_in=",Jastrow_d_length_in
     end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! REPLACE FOLLOWING WITH INPUT CODE !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! These are now all read in or set through readtxt subroutine
    ! THESE STILL NEED ADDING TO INPUT CODE
     ! use_biased_optimiser = .false.
     ! n_optimiser_steps = 10
     ! allocate( proton_numbers_in(n_atoms_in) )
     ! proton_numbers_in = 1.0_dp

     ! basis_type_in = slater_1s_code

     ! MCMC adapt parameters
    mcmc_adapt_steps = 10000
    mcmc_adapt_s_0 = 0.1_dp
    mcmc_adapt_s_max = 0.4_dp
    mcmc_adapt_s_min = 0.03_dp
    mcmc_adapt_memory = 500.0_dp
    mcmc_adapt_interval = 100

    ! Biased Optimiser Init parameters
    ker_var_in = 1.0_dp
    ker_lengthscale_in = 1.0_dp
    constant_mean_prior_in = 0.0_dp
    optim_rate_para_in = 0.5_dp



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! END OF INPUT CODE!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Allocate atom positions
  allocate( atom_coords_in(3,n_atoms_in) )
  if (n_atoms_in == 1) then
   atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp] ! 1 atom at origin
  else if (n_atoms_in == 2) then
   atom_coords_in(:,1)=[bond_length/2.0_dp,0.0_dp,0.0_dp] !2 atoms evenly placed about origin
   atom_coords_in(:,2)=[-bond_length/2.0_dp,0.0_dp,0.0_dp]
  end if

  ! Initialisation of basis and trials

  call omp_set_num_threads(n_omp_threads)

  call initialise_basis(n_electrons_in, n_basis_functions_per_atom_in, n_atoms_in, atom_coords_in,&
                        n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in,&
                        proton_numbers_in,basis_type_in)

  !!!MCMC SETUP!!!!
  allocate(x_0(3*n_electrons_in))
  x_0(1:3) = atom_coords_in(:,1) + 0.2_dp ! Displace starting electrons from atoms
  if (n_electrons_in == 2) then
    x_0(4:6) = atom_coords_in(:,2) + 0.2_dp
  end if
  n_burned = n_MCMC_steps/10 ! Integer division
  thinning_interval = n_MCMC_steps/10000
  norm_coeff = real(thinning_interval,dp)/real(n_MCMC_steps-n_burned,dp)
  allocate(mcmc_run((n_MCMC_steps-n_burned)/thinning_interval+1,n_space_dims))

  print*, "Setup Complete"
  !!!!!!!!!!! SEARCH BEGINS HERE!!!!!!!!!

  allocate(trials(n_trials, number_dofs))
  call  latin_hypercube(trials,dof_bounds, n_trials, search_seed)

  allocate(trial_energies(n_trials))

  ! Main loop over trials in parameter space
  !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
  do i=1, n_trials
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
  print*, "Initial Hypecube Search Done"
  !!!!!!!!! NO OPTIMISER!!!!!!!!!!!!
  if(.not.(use_biased_optimiser)) then

    call find_best_params(trial_energies,trials)
    call param_wall_check(dof_bounds)
    output_best_parameters = best_params
    print*,"best trial=", best_trial
    print*,"min energy=", trial_energies(best_trial),"Hartree Energy"
    print*,"min energy=", electronvolt*trial_energies(best_trial),"eV"
    print*,"best parameters=", best_params
  end if

  !!!!!!!!!OPTIMISER!!!!!!!!!!!!!!!
  if(use_biased_optimiser) then
    ! THESE HARDCODED CONSTANTS NEED CHANGING TO VARIABLES
    call Bi_Op_init(trials, trial_energies, n_trials, number_dofs,&
      ker_var_in, ker_lengthscale_in, constant_mean_prior_in, optim_rate_para_in,&
      n_trials,n_trials, n_optimiser_steps)

    print*, "Optimiser Initiated"

    allocate(trials_temp(n_trials, number_dofs))

    do j=1,n_optimiser_steps
      !seeding is needed for stoch_grad, will crash without
      call random_seed(size=n_seed)
      allocate(seed(n_seed))
      call random_seed(get=seed)
      trials_temp = Bi_Op_step(trials, trial_energies, n_trials, n_trials, seed, j)
      print*,"Bi_Op stepped", j, "times"
      trials=trials_temp
      !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
      do i=1, n_trials

          !print *, trials(i,:)
          ! Run MCMC to generate coordinates in space
        call mcmc_adapt(s, log_density, x_0, mcmc_adapt_steps, mcmc_adapt_s_0, e_code, mcmc_adapt_s_max,&
        mcmc_adapt_s_min, mcmc_adapt_memory, mcmc_adapt_interval, trials(i,:),n_space_dims)
        call mcmc_sample(mcmc_run, log_density, x_0, n_MCMC_steps, n_burned, thinning_interval, s,e_code,&
         trials(i,:),n_space_dims)
          ! Compute energy by summing over positions
          E=0.0_dp
          do loop=1,(n_MCMC_steps-n_burned)/thinning_interval
          !print*, "mcmcrun=",mcmc_run(loop,:)
          !print*, "trials=",trials(i,:)
          !print*, "Ham=",reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
          E = E + reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
          end do

          trial_energies(i) = norm_coeff*E
      end do
      call stoch_grad_exit()
      deallocate(seed)
    end do
    print*, "Completed Optimiser trials"
    allocate(min_gp_params(size(dof_bounds,1)))

    call gp_min_stored(min_gp_params,min_gp_E)
    call find_best_params(trial_energies,trials)
    if (trial_energies(best_trial(1)) .le. min_gp_E) then
        print*, 'smallest vaule was found in final batch'
        print*,"best trial=", best_trial
        print*,"min energy=", trial_energies(best_trial),"Hartree Energy"
        print*,"min energy=", electronvolt*trial_energies(best_trial),"eV"
        print*,"best parameters=", best_params
        output_best_parameters = best_params
    else
        print*, 'smallest vaule was found in a previous batch'
        print*,"min energy=", min_gp_E,"Hartree Energy"
        print*,"min energy=", electronvolt*min_gp_E,"eV"
        print*,"best parameters=", min_gp_params
        output_best_parameters = min_gp_params
    end if
    deallocate(trials_temp,min_gp_params)
  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! OUTPUT CODE GOES HERE!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 !!!!!NOW USE output_best_parameters FOR THE OUTPUT PLOTTING ETC!!!!!
  ! !1 atom results
  !  if (n_atoms_in == 1) then
  !    wave_1d = wave1d(best_params(1:2), -10,10)
  !    wave_2d = wave2d(best_params(1:2), -10,10)
  !    wave_3d = wave3d(best_params(1:2), -10,10)

  !    ! Writing wavefunction results to NetCDF file
  !    call atom1_netcdf(filename_result, best_params, wave_1d, wave_2d, wave_3d)

  !    ! Creating coordinate text files
  !    call xyz_coords1d(-10,10)
  !    call xyz_coords2d(-10,10)
  !    call xyz_coords3d(-10,10)
  !  end if

  ! ! 2 atom results
  ! ! Use xyz3d.txt for the coordinates, these are only calculated in 3D for now.

  !  if (n_atoms_in == 2) then
  !    ! Electron density
  !    ele_2atom = ele_den_2atom(best_params, -10,10, n_MCMC_steps)
  !    ! Fixed wavefunction
  !    wave_2atom  = wave_2atoms(best_params, -10, 10)

  ! !   ! Netcdf
  !    call atom2_netcdf(filename_result_2, best_params, wave_2atom, ele_2atom)
  !    call xyz_coords3d(-10,10)
  !  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!! OUTPUT CODE FINISH   !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! deallocation
  deallocate(trials,trial_energies,atom_coords_in,x_0,mcmc_run,proton_numbers_in)
  call deinitialise_basis

  ! Report Timings
  print*, "Output Complete"
  call cpu_time(cpu_finish_time)
  call system_clock(real_finish_time,rate)
  print '("Total cpu runtime = ",f7.3," seconds.")',cpu_finish_time-cpu_start_time
  print'("Total real runtime = ",f7.3," seconds.")',real(real_finish_time-real_start_time,kind=dp)/real(rate,kind=dp)

end program main_driver

!> @brief Driver for testng the biased optimization routines

module log_rho_mod
use shared_constants
use basis_functions
implicit none

contains

function log_rho(x,dof) result(retval)
    implicit none
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: dof
    real(dp) :: retval


    retval = 2*log(abs(wave_function(x,dof)))

end function log_rho

end module log_rho_mod

program main
    use Biased_Optim
    use shared_constants
    use stoch_grad
    use shared_constants
    use gp_surrogate
    use basis_functions
    use param_search
    use log_rho_mod
    use mcmc
    !$ use omp_lib
    implicit none

    !actually define this somehow

    integer :: n_threads
    real(dp), dimension(3,1) :: atom_coords_in
    real(dp), dimension(:, :), allocatable :: trials, trials_temp
    real(dp), dimension(:), allocatable :: trial_energies !one batch
    real(dp) :: min_gp_E !minimum of stored data in gp
    real(dp), dimension(:), allocatable :: min_gp_params !params at min_gp_E

    integer :: i,j,n_batches

    integer, parameter :: space_dim = 6

    !MCMC variables
    real(dp), dimension(space_dim) :: x_0=0.1_dp
    integer, parameter :: n_steps=1000000, n_burned=100000, thinning_interval=1000
    integer :: n_seed, e_code, loop
    real(dp) :: s, a_ave, E, norm_coeff, start_time, finish_time
    real(dp), dimension((n_steps-n_burned)/thinning_interval+1,space_dim) :: mcmc_run
    integer, dimension(:), allocatable:: seed

    n_threads = OMP_GET_MAX_THREADS()
    print*, n_threads

    n_batches = 100

    call omp_set_num_threads(n_threads)
    call cpu_time(start_time)


    ! Initialisation of basis and trials
    atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp]
    !atom_coords_in(:,2)=[-1.5_dp,0.0_dp,0.0_dp]
    call initialise_basis(2,1,1,atom_coords_in)

    allocate(trials(n_threads, number_dofs))

    call latin_hypercube(trials,dof_bounds, n_threads, 123) !Now takes (dof_bounds, trials, seed_in) as input. Same is true for random_search_grid.
    !print*, trials
    allocate(trials_temp(size(trials,1),size(trials,2)))
    allocate(trial_energies(size(trials,1)))
    !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
    do i=1, n_threads

        !print *, trials(i,:)
        ! Run MCMC to generate coordinates in space
        call mcmc_adapt(s, log_rho, x_0, 10000, 0.1_dp, e_code, 0.4_dp, 0.03_dp, 500.0_dp, 100, trials(i,:),space_dim)
        ! print*,"s"
        call mcmc_sample(mcmc_run, log_rho, x_0, n_steps, n_burned, thinning_interval, s,e_code, trials(i,:),space_dim,a_ave)
        !print*, "sample"
        ! Compute energy by summing over positions
        E=0.0_dp
        do loop=1,(n_steps-n_burned)/thinning_interval
        !print*, "mcmcrun=",mcmc_run(loop,:)
        !print*, "trials=",trials(i,:)
        !print*, "Ham=",reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
        E = E + reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
        end do
        norm_coeff = real(thinning_interval,dp)/real(n_steps-n_burned,dp)
        trial_energies(i) = norm_coeff*E
    end do

    print*, "Latin hcube done"

    !note the constants set the prior, might want to be user inputs
    call Bi_Op_init(trials, trial_energies, n_threads, size(dof_bounds,1),&
    1.0_dp, 1.0_dp, 0.0_dp, 0.5_dp, n_threads,n_threads, n_batches)

    print*,"inited Bi_Op"

   !call gp_debug_print_state()

    !the main loop, either set number of iterations or converge to tolerence
    do j=1,n_batches
        !seeding is needed for stoch_grad, will crash without
        call random_seed(size=n_seed)
        allocate(seed(n_seed))
        call random_seed(get=seed)
        trials_temp = Bi_Op_step(trials, trial_energies, n_threads, n_threads, seed, j)
        print*,"Bi_Op stepped", j, "times"
        trials=trials_temp
        !$OMP parallel do default(shared) private(i,s,e_code,mcmc_run, E,loop)
        do i=1, n_threads

            !print *, trials(i,:)
            ! Run MCMC to generate coordinates in space
            call mcmc_adapt(s, log_rho, x_0, 10000, 0.1_dp, e_code, 0.4_dp, 0.03_dp, 500.0_dp, 100, trials(i,:),space_dim)
            !print*,"s"
            call mcmc_sample(mcmc_run, log_rho, x_0, n_steps, n_burned, thinning_interval, s,e_code, trials(i,:),space_dim,a_ave)
            !print*, "sample"
            ! Compute energy by summing over positions
            E=0.0_dp
            do loop=1,(n_steps-n_burned)/thinning_interval
            !print*, "mcmcrun=",mcmc_run(loop,:)
            !print*, "trials=",trials(i,:)
            !print*, "Ham=",reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
            E = E + reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
            end do
            norm_coeff = real(thinning_interval,dp)/real(n_steps-n_burned,dp)
            trial_energies(i) = norm_coeff*E
        end do
        call stoch_grad_exit()
        deallocate(seed)
    end do


    allocate(min_gp_params(size(dof_bounds,1)))

    call gp_min_stored(min_gp_params,min_gp_E)
    call find_best_params(trial_energies,trials)
    if (trial_energies(best_trial(1)) .le. min_gp_E) then
        print*, 'smallest vaule was found in final batch'
        print*,"best trial=", best_trial
        print*,"min energy=", trial_energies(best_trial),"Hartree Energy"
        print*,"min energy=", electronvolt*trial_energies(best_trial),"eV"
        print*,"best parameters=", best_params
    else
        print*, 'smallest vaule was found in a previous batch'
        print*,"min energy=", min_gp_E,"Hartree Energy"
        print*,"min energy=", electronvolt*min_gp_E,"eV"
        print*,"best parameters=", min_gp_params
    end if

    call param_wall_check(dof_bounds)
    !could be important for figuring out scalling.
    call cpu_time(finish_time)
    print*, "Total runtime = ",finish_time-start_time


    call gp_debug_print_state()
    call gp_exit()
    ! deallocation
    deallocate(trials,trial_energies)
    call deinitialise_basis

    call initialise_basis(2,1,1,atom_coords_in)
            !print *, trials(i,:)
            ! Run MCMC to generate coordinates in space
    call mcmc_adapt(s, log_rho, x_0, 10000, 0.1_dp, e_code, 0.4_dp, 0.03_dp, 500.0_dp, 100, min_gp_params,space_dim)
    !print*,"s"
    call mcmc_sample(mcmc_run, log_rho, x_0, n_steps, n_burned, thinning_interval, s,e_code,min_gp_params,space_dim,a_ave)
    !print*, "sample"
    ! Compute energy by summing over positions
    E=0.0_dp
    do loop=1,(n_steps-n_burned)/thinning_interval
    !print*, "mcmcrun=",mcmc_run(loop,:)
    !print*, "trials=",trials(i,:)
    !print*, "Ham=",reduced_hamiltonian(mcmc_run(loop,:),trials(i,:))
    E = E + reduced_hamiltonian(mcmc_run(loop,:),min_gp_params)
    end do
    norm_coeff = real(thinning_interval,dp)/real(n_steps-n_burned,dp)
    print*, norm_coeff*E
    call deinitialise_basis


end program main

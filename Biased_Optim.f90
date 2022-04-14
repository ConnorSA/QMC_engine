!> @brief Biased optimization subroutines and functions

module Biased_Optim
    use shared_constants
    use stoch_grad
    use gp_surrogate
    !$ use omp_lib
    implicit none

    !takes batch of data, returns next set of points to explore
    !algorithm from Parallel Bayesian Global Optimization of Expensive Functions Wang et al. 2019

    real(dp) :: current_best_E, constant_mean_value, gamma
    integer :: n_dof, no_samples
    logical :: GP_uptodate = .False.

    interface Bi_Op_init
        module procedure Bi_Op_init_constant_mean, Bi_Op_init_arb_mean
    end interface

    private current_best_E, constant_mean_value, n_dof, gamma, no_samples, gaussian, zero_func, constant_mean
contains

    function gaussian(density_dimension) result(x_prop)
        implicit none
        integer, intent(in) :: density_dimension
        real(dp), dimension(2*((density_dimension+1)/2)) :: u !uniform random sample,
        !note the dimension is the nearest even number above the given dimension
        !integer division will autmatically round towards 0
        real(dp), dimension(density_dimension) :: xi !normal random sample
        real(dp), dimension(density_dimension) :: x_prop
        integer :: i
        !proposes, uses same dist as described in log_prop_distinteger, dimension(:), allocatable:: seed
        !if proposal is changed, might need updated

        !implementing box muller to generate a normal
        call random_number(u)
        do i=1,density_dimension/2
            xi(i) = sqrt(-2.0_dp*log(u(i)))*cos(2.0_dp*pi*u(i*2))
            xi(i*2) = sqrt(-2.0_dp*log(u(i)))*sin(2.0_dp*pi*u(i*2))
        end do
        if (modulo(density_dimension,2).ne.0) then
            xi(density_dimension) = sqrt(-2.0_dp*log(u(density_dimension)))*cos(2.0_dp*pi*u(density_dimension+1))
        end if

        !basic GRW
        x_prop = xi

    end function gaussian

    function zero_func(x,dim) result(out)
        !a 0 function
        implicit none
        real(dp), dimension(:), intent(in) :: x !param space point
        integer, intent(in) :: dim
        real(dp) :: out

        out=0.0_dp

    end function  zero_func

    function constant_mean(x) result(out)
        !a constant mean prior
        implicit none
        real(dp), dimension(:), intent(in) :: x !param space point
        real(dp) :: out

        out=constant_mean_value

    end function  constant_mean

    !a tempoary function, repersents a check for restart files
    function find_restart_file() result(retval)
        implicit none
        logical :: retval

        retval=.False.

    end function find_restart_file

    !initialises based on a first batch of data (taken from a latin hypercube)
    subroutine Bi_Op_init_constant_mean(param_init_data, energy_init_data, n_data_in, n_dof_in,&
        ker_var, ker_lengthscale, constant_mean_prior, optim_rate_para, optim_no_samples,n_threads,n_loops_to_do)
        implicit none
        integer, intent(inout) :: n_loops_to_do
        integer, intent(in) :: n_data_in, n_dof_in, optim_no_samples,n_threads
        real(dp), dimension(n_data_in, n_dof_in), intent(in) :: param_init_data
        real(dp), dimension(n_data_in), intent(in) :: energy_init_data
        real(dp), intent(in) :: ker_var, ker_lengthscale, constant_mean_prior, optim_rate_para
        integer gp_n_data, gp_n_dof, n_cycles
        real(dp) :: kernel_var, kernal_inv_length
        real(dp), dimension(:,:), allocatable :: param_data
        real(dp), dimension(:), allocatable :: E_data
        real(dp), dimension(:,:), allocatable :: param_pres, param_cov
        real(dp), dimension(:), allocatable :: data_mean

        if (find_restart_file()) then
            call read_restart_file_sizes(gp_n_data, gp_n_dof)
            allocate(param_data(gp_n_data, gp_n_dof))
            allocate(E_data(gp_n_data))
            allocate(param_pres(gp_n_data, gp_n_data))
            allocate(param_cov(gp_n_data, gp_n_data))
            allocate(data_mean(gp_n_data))
            call read_restart_file_data(n_cycles, no_samples, current_best_E,&
            constant_mean_value, kernel_var, gamma, kernal_inv_length,&
            gp_n_data, E_data, gp_n_data, data_mean, gp_n_data, gp_n_data, param_pres, &
            gp_n_data, gp_n_data, param_pres, gp_n_data, gp_n_dof, param_data)
            n_dof = gp_n_dof

            call gp_restart(gp_n_data, gp_n_dof,kernel_var, kernal_inv_length,&
            param_data, E_data,param_pres, param_cov, data_mean, constant_mean,&
             zero_func, n_threads)
            GP_uptodate = .False.
            current_best_E = minval(E_data)
            n_loops_to_do=n_loops_to_do-(n_cycles-1)
        else
            constant_mean_value = constant_mean_prior
            gamma = optim_rate_para
            no_samples=optim_no_samples
            n_dof = n_dof_in

            call gp_init(constant_mean, zero_func, ker_var, ker_lengthscale, param_init_data, energy_init_data, n_data_in,&
             n_dof_in,n_threads)
            GP_uptodate = .False.

            current_best_E = minval(energy_init_data)
        end if

    end subroutine Bi_Op_init_constant_mean

    subroutine Bi_Op_init_arb_mean(param_init_data,  energy_init_data, n_data_in, n_dof_in, ker_var,&
        ker_lengthscale, mean_prior_func, mean_prior_dx, optim_rate_para,&
         optim_no_samples,n_threads, n_loops_to_do)
       implicit none
       integer, intent(inout) :: n_loops_to_do
       integer, intent(in) :: n_data_in, n_dof_in, optim_no_samples,n_threads
       real(dp), dimension(n_data_in, n_dof_in), intent(in) :: param_init_data
       real(dp), dimension(n_data_in), intent(in) :: energy_init_data
       real(dp), intent(in) :: ker_var, ker_lengthscale,  optim_rate_para
       procedure(mean_func_interface) :: mean_prior_func
       procedure(mean_func_interface_dx) :: mean_prior_dx

       n_dof = n_dof_in
       gamma = optim_rate_para
       no_samples=optim_no_samples

       call gp_init(mean_prior_func, mean_prior_dx, ker_var, ker_lengthscale, param_init_data, energy_init_data,&
        n_data_in, n_dof_in, n_threads)
       GP_uptodate = .False.

       current_best_E = minval(energy_init_data)

    end subroutine Bi_Op_init_arb_mean

    !do one batch of optim
    function Bi_Op_step(param_update_data,  energy_update_data, n_new_data, threads, seed, n_cycles) result(best_locs)
        implicit none
        integer, intent(in) :: n_new_data, n_cycles
        integer, intent(in) :: threads !number of threads/restarts goes here
        integer, dimension(:), intent(IN):: seed
        real(dp), dimension(n_new_data, n_dof), intent(in) :: param_update_data
        real(dp), dimension(n_new_data), intent(in) :: energy_update_data
        real(dp), dimension(threads,threads, n_dof) ::  per_thread_best_loc
        real(dp), dimension(0:threads) ::  m,z
        real(dp), dimension(0:threads,0:threads) ::  c
        real(dp), dimension(threads) ::  per_thread_best_qei
        real(dp), dimension(threads,n_dof) :: best_locs
        real(dp) :: best_E_in_batch
        integer :: i, j, best_restart, n_seed, gp_n_data, gp_n_dof
        real(dp) :: kernel_var, kernal_inv_length
        real(dp), dimension(:,:), allocatable :: param_data
        real(dp), dimension(:), allocatable :: E_data
        real(dp), dimension(:,:), allocatable :: param_pres, param_cov
        real(dp), dimension(:), allocatable :: data_mean

        if (GP_uptodate) then
            call gp_update(param_update_data,  energy_update_data, n_new_data, 'F')
        end if

        !updating best point
        best_E_in_batch = minval(energy_update_data)
        if (best_E_in_batch < current_best_E) current_best_E=best_E_in_batch
        print*, current_best_E

        call stoch_grad_init(threads,current_best_E,seed)

        ! !doing the restart file thing
        ! call gp_return_size_data(gp_n_data, gp_n_dof)
        ! allocate(param_data(gp_n_data, gp_n_dof))
        ! allocate(E_data(gp_n_data))
        ! allocate(param_pres(gp_n_data, gp_n_data))
        ! allocate(param_cov(gp_n_data, gp_n_data))
        ! allocate(data_mean(gp_n_data))
        ! call gp_return_state_data(kernel_var, kernal_inv_length,  param_data, E_data,&
        ! param_pres, param_cov, data_mean)
        ! !assuming calling order of (scalars(rank 0),..., vector_dim, vector(rank 1),..., matrix_dim_1,matrix_dim2, matrix(rank 2))
        ! call write_restart_file(gp_n_data, gp_n_dof,n_cycles,no_samples, current_best_E,&
        ! constant_mean_value, gamma, kernel_var, kernal_inv_length, &
        ! gp_n_data, E_data, gp_n_data, data_mean, gp_n_data, gp_n_data, param_pres, &
        ! gp_n_data, gp_n_data, param_pres, gp_n_data, gp_n_dof, param_data)


        !doing grad_accent, can be parrallel
        !$OMP parallel do default(shared) private(i,m,c,z)
        do i=1,threads
            !the "no change" state
            z=0.0_dp
            m=0.0_dp
            call grad_accent(N_points(i,:,:),threads,no_samples,gamma,per_thread_best_loc(i,:,:))
            m(1:threads) = current_best_E-gp_mu_post(per_thread_best_loc(i,:,:), threads)
            call get_cholesky(c, gp_k_post_same_x(per_thread_best_loc(i,:,:), threads))
            c= -c
            !c = -get_cholesky(gp_k_post_same_x(per_thread_best_loc(i,:,:), threads))
            do j=1,no_samples
                z(1:threads) = gaussian(threads)
                per_thread_best_qei(i) =  maxval(m + matmul(c, z))
            end do
        end do


        best_restart = maxloc(per_thread_best_qei,dim=1)

        best_locs = per_thread_best_loc(best_restart,:,:)

        GP_uptodate = .True.

    end function Bi_Op_step

end module Biased_Optim

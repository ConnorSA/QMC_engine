 !> @brief Functions/subroutines associated with building/initialisation the parameter search space
 !> and finding the best paramters from a given MCMC run.
 module param_search
  use iso_fortran_env
  use shared_constants
  implicit none


  ! Public Variables
  real(dp), dimension(:), allocatable, protected :: best_params
  integer, dimension(1), protected :: best_trial
  contains
  !> @brief random_search_grid returns MxN grid of test points for M free parameters and N trials distributed
  !> @param trials is a 2D array containing the parameter set for each trial used in MCMC.
  !> @param param_bounds is a 2D array containing the upper and lower bounds for each parameter in the wavefunciton.
  !> @param n_trials is the number of trials to be tested in MCMC. This defines the binning of the Latin hypercube.
  !> @param seed_in is an integer that defines the seed being used to generate the latin Latin hypercube.
  subroutine random_search_grid(trials,param_bounds,n_trials, seed_in)
    real(dp), dimension(:,:), allocatable ,intent(in) :: param_bounds !length of this will be the "N"
    real(dp), dimension(:, :),intent(inout) :: trials ! (x,y) x is each trial, y is each param
    integer, intent(in) :: n_trials, seed_in!this is the "M"
    real(dp), dimension(:),allocatable :: rand_num
    integer, dimension(2) :: arr_shape
    integer :: i, n
    integer,allocatable :: seed(:)
    !initialising a certain seed is done as shown below.
    call random_seed(size=n)
    allocate(seed(n))
    seed = seed_in   !set seed to that read in from input.
    call random_seed(put=seed)


    ! arr_shape(1) contains the number of variables, param_bounds(x,y) where x is for each param and y for limits.
    arr_shape = shape(param_bounds)

    allocate(rand_num(n_trials))

    do i = 1,arr_shape(1)
    call random_number(rand_num)
    trials(:, i) = rand_num*(param_bounds(i,1) - param_bounds(i,2)) + param_bounds(i,2)
    end do

  end subroutine random_search_grid




!> @brief latin_hypercube returns MxN grid of test points for M free parameters and N trials on a latin hypercube for given bounds.
!resource for scrambling arrays:http://fortranwiki.org/fortran/show/scramble
!> @param trials is a 2D array containing the parameter set for each trial used in MCMC.
!> @param param_bounds is a 2D array containing the upper and lower bounds for each parameter in the wavefunciton.
!> @param n_trials is the number of trials to be tested in MCMC. This defines the binning of the Latin hypercube.
!> @param seed_in is an integer that defines the seed being used to generate the latin Latin hypercube.
  subroutine latin_hypercube(trials,param_bounds,n_trials, seed_in)
    real(dp), dimension(:,:), intent(in) :: param_bounds !length of this will be the "N"
    real(dp), dimension(:, :),intent(inout) :: trials ! (x,y) x is each trial, y is each param
    integer, intent(in) :: n_trials, seed_in !n_trials is the "M" and is specified as an input, seed_in is an integer used to set the seed.
    real(dp), dimension(:),allocatable ::  zero_to_one
    real(dp) :: rand_num, temp_val
    integer, dimension(2) :: arr_shape
    integer :: i,j, k, h, n
    integer,allocatable :: seed(:)
    !setting a certain seed is done as shown below.
    call random_seed(size=n)
    allocate(seed(n))
    seed = seed_in    !set seed value to that read in from input.
    call random_seed(put=seed) !actually set the random seed for repeatability


    ! arr_shape(1) contains the number of variables, param_bounds(x,y) where x is for each param and y for limits.
    arr_shape = shape(param_bounds)

    allocate(zero_to_one(n_trials))
    do i = 1,n_trials
      zero_to_one(i) = (real(i, dp)-1.0_dp)/(n_trials-1.0_dp) !extra bit needed to make sure we are binning right in the param search space!
    end do

    ! initialize latin cube in order.
    do i = 1,arr_shape(1)
      trials(:, i) = zero_to_one*(param_bounds(i,1) - param_bounds(i,2)) + param_bounds(i,2)
    end do

    !scramble the array for each param.
    do h = 1, 2 !sufficient scrambling.
      do i = 1,arr_shape(1) !loop for each parameter.
        do j = 1, n_trials !scramble over this loop
          call random_number(rand_num)
          k = 1 + FLOOR((n_trials+1-1)*rand_num) !this is how we are randomly selecting an index to swap.
          temp_val=trials(k,i) !store temp value to swap
          trials(k,i) = trials(j,i) !swap points in the grid
          trials(j,i) = temp_val !swap back the temp value. simple as.
        end do
      end do
    end do

    deallocate(zero_to_one,seed)
  end subroutine latin_hypercube



  !> @brief find_best_params returns the 1D array that contain the current best set of parameters found from a MCMC run.
  !> the output is stored in the global variable best_params.
  !> @param trial_energies is a 1D array of the energies found for each parameter set.
  !> @param trials is a 2D array containing the parameter set for each trial used in MCMC.

  subroutine find_best_params(trial_energies, trials)
    real(dp), dimension(:, :), allocatable, intent(in) :: trials
    real(dp), dimension(:), allocatable, intent(in) :: trial_energies
    integer, dimension(1) :: arr_shape
    arr_shape = shape(trials(1,:))
    if(allocated(best_params)) deallocate(best_params)
    allocate(best_params(arr_shape(1)))
    best_trial = minloc(trial_energies) !index corresponding to the most negative, "best", energy.
    best_params = trials(best_trial(1),:)
  end subroutine find_best_params

  !> @brief param_wall_check evaluates if the best parameter set is close to the parameter bounds used to generate
  !> the test points. The tolerence is set at 5%.
  !> @param param_bounds is a 2D array containing the upper and lower bounds for each parameter in the wavefunciton.
  subroutine param_wall_check(param_bounds)
    real(dp), dimension(:,:), allocatable ,intent(in) :: param_bounds
    real(dp), parameter :: tol = 0.05_dp !percentage tol as determined by the min and max bounds.
    real(dp) :: bound_gap
    integer :: i

    do i = 1, size(best_params)
      bound_gap = abs(param_bounds(i,1) - param_bounds(i,2))
      if ( abs(best_params(i) - param_bounds(i,1))/bound_gap < tol &
      .or. abs(best_params(i) - param_bounds(i,2))/bound_gap < tol ) then
        print*, "parameter ", i,  " is close to wall"
        !for later: might want to instead store i in and array of potential problem parameters.
        !As oppsed to just printing out the problem.
      end if
    end do
  end subroutine param_wall_check

end module

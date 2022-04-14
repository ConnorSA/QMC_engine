!> @brief Subroutines for obtaining a gradient estimation through Cholesky decomposition
module gradient_estimator
  use iso_fortran_env
  use shared_constants
  use gp_surrogate
  use priors
  implicit none

  !public variables to be used between subroutines and functions (in this module only).
  !F_ij is the Workspace, L_ij is the Cholesky decomposition of the covariance matrix.
  !gradient_matrix contains the derivatives for each test point (threads) and each parameter (q by n_dof)
  real(dp), dimension(:,:), allocatable :: F_ij, L_ij, gradient_matrix
  real(dp) :: check_max

  contains
  !> @brief performs Cholesky decomposition without destroying the old array. also with 0 at i=0 padding.
  !> @param cholesky_decomp is where the result is stored.
  !> @param matrix_in 2D matrix to be decomposed.
  subroutine get_cholesky(cholesky_decomp,matrix_in)
  real(dp), dimension(:,:) :: matrix_in, cholesky_decomp
  real(dp), dimension(:,:), allocatable :: get_cholesky_temp
  integer :: info
  integer, dimension(2) :: shape_arr
  get_cholesky_temp = matrix_in
  shape_arr = shape(matrix_in)
  call DPOTRF( 'L', shape_arr(1), get_cholesky_temp, shape_arr(1),info) !this is the LAPACK routine for getting Cholesky 'L' = lower triangle. info has an error state
  cholesky_decomp = 0.0_dp
  cholesky_decomp(1:shape_arr(1), 1:shape_arr(1)) = -get_cholesky_temp
  end subroutine get_cholesky



  !> @brief builds the "m_vector" as described in: arXiv:1602.05149v4  [stat.ML]  5 May 2019
  !> contains the difference between the mean and best result in a vector.
  !> @param x is 2D array containing the param_inputs for each thread.
  !> @param previous_best is a real number containing the best energy from the prior trial.
  function m_vector(x, previous_best)
  real(dp), dimension(:), allocatable :: m_vector, pbest_arr
  real(dp), dimension(:,:) , intent(in):: x
  real(dp), intent(in) :: previous_best
  integer, dimension(2) :: shape_arr
  shape_arr = shape(x)
  allocate(m_vector(0:shape_arr(1)))
  allocate(pbest_arr(shape_arr(1)))
  pbest_arr = previous_best !need an array to compare each thread/trial/point in param space.
  m_vector(0) = 0.0_dp
  m_vector(1:shape_arr(1))=pbest_arr - gp_mu_post(x, shape_arr(1))
  end function m_vector



  !> @brief build the "z_vector" as described in: arXiv:1602.05149v4  [stat.ML]  5 May 2019
  !> @param z_in is a 1D array of normally distributed numbers.
  function z_vector(z_in)
  real(dp), dimension(:), allocatable :: z_vector
  real(dp), dimension(:), intent(in) :: z_in
  integer :: size_z
  size_z = size(z_in)
  allocate(z_vector(0:size_z))
  z_vector(0) = 0
  z_vector(1:size_z) = z_in
  end function z_vector


  !> @brief finds the location in the vector that corresponds to the maximum of "m + CZ"
  !> @param m_vec as returned from the m_vector function.
  !> @param z_vec as returned from the z_vector function.
  !> @param c_mat corresponds to the Cholesky decomposition of the covariance matrix.
  function find_max_thread(m_vec, z_vec, c_mat)
    real(dp), dimension(:,:) , intent(in):: c_mat
    real(dp), dimension(:), intent(in) :: m_vec, z_vec
    integer, dimension(1) :: find_max_thread
    integer, dimension(2) :: shape_arr
    shape_arr = shape(c_mat)
    find_max_thread = maxloc(m_vec + matmul(c_mat, z_vec)) - 1
    !-1 is shift the indexing back to 0 from 1.
  end function find_max_thread



  !> @brief check if the maximum found in the find_max_thread is a unique value.
  !> @param m_vec as returned from the m_vector function.
  !> @param z_vec as returned from the z_vector function.
  !> @param c_mat corresponds to the Cholesky decomposition of the covariance matrix.
  !> @param find_max_thread integer containing the location of the maximum of "m + CZ"
  subroutine check_max_is_unique(m_vec, z_vec, c_mat, find_max_thread)
    integer, dimension(1), intent(in) :: find_max_thread
    real(dp), dimension(:,:) , intent(in):: c_mat
    real(dp), dimension(:), intent(in) :: m_vec, z_vec
    real(dp), dimension(:), allocatable :: f_vector
    integer :: i
    !make sure indexing starts from 0, so is consistent.
    allocate(f_vector(0:( size(z_vec) -1 )))
    !by default set = 1, so it does nothing when multiplied with something.
    check_max = 1.0_dp
    f_vector = m_vec + matmul(c_mat, z_vec)
    !loop over each element of the calculated vector and check if it equals the max of that vector
    !but in more than one location.
    do i = 0, size(f_vector) -1
      if ( (f_vector(i) .EQ. f_vector(find_max_thread(1))) .AND. (i .NE. find_max_thread(1)) ) then
        check_max=0.0_dp
      end if
    end do
    deallocate(f_vector)
  end subroutine check_max_is_unique


  !> @brief initialises the workspace F and Cholesky decomposition of the covariance matrix.
  !> @param x is 2D array containing the param_inputs for each thread.
  !> @param z_in is a 1D array of normally distributed numbers.
  !> @param previous_best is a real number containing the best energy from the prior trial.
  subroutine init_diff_F(x,z_in, previous_best)
  real(dp), dimension(:,:) , intent(in):: x
  real(dp), intent(in) :: previous_best
  real(dp), dimension(:), intent(in) :: z_in
  real(dp), dimension(:), allocatable :: m_vec, z_vec
  real(dp), dimension(:,:), allocatable :: c_mat, init_F_ij
  integer, dimension(2) :: shape_arr
  integer, dimension(1) :: loc_max_thread
  integer :: i, j
  shape_arr = shape(x)
  allocate(m_vec(0:shape_arr(1)))
  allocate(z_vec(0:shape_arr(1)))
  allocate(c_mat(0:shape_arr(1), 0:shape_arr(1)  ))
  allocate(init_F_ij(0:shape_arr(1), 0:shape_arr(1)))

  if (allocated(F_ij)) then
    deallocate(F_ij)
  end if
  allocate(F_ij(1:shape_arr(1), 1:shape_arr(1)))

  if (allocated(L_ij)) then
    deallocate(L_ij)
  end if
  allocate(L_ij(1:shape_arr(1), 1:shape_arr(1)))
  init_F_ij = 0.0_dp
  m_vec = m_vector(x, previous_best)
  z_vec = z_vector(z_in)

  call get_cholesky(c_mat, gp_k_post(x, shape_arr(1))  )
  loc_max_thread = find_max_thread(m_vec, z_vec, c_mat)
  !print*, loc_max_thread
  !This is the analytical result of diffing the max(f) w.r.t each L(ij)
  init_F_ij( loc_max_thread(1) , 1:loc_max_thread(1)  ) = -z_vec(1:loc_max_thread(1))
  !putting array back into correct size (removing index 0, padding)
  F_ij =   init_F_ij(1:shape_arr(1) , 1:shape_arr(1))
  !put cholesky in same size array, removing 0 index padding.
  L_ij = c_mat(1:shape_arr(1) , 1:shape_arr(1))


  call check_max_is_unique(m_vec, z_vec, c_mat, loc_max_thread)
  deallocate(c_mat, z_vec, m_vec, init_F_ij)
  end subroutine init_diff_F


  !> @brief performs the backward difference differentiation on the workspace F
  subroutine bkwd_diff_F()
  integer :: i, j, k, N
  integer, dimension(2) :: shape_arr
  shape_arr = shape(F_ij)
  N = shape_arr(1)

  do k=N,1, -1
    if (abs(L_ij(k,k)) > 0.0_dp) then
    !row operations
    do j = k+1, N
      do i = j, N
        F_ij(i,k) = F_ij(i,k) - F_ij(i,j)*L_ij(j,k)
        F_ij(j,k) = F_ij(j,k) - F_ij(i,j)*L_ij(i,k)
      end do
    end do
    !lead column
    do j = k+1, n
      F_ij(j,k) = F_ij(j,k)/L_ij(k,k)
      F_ij(k,k) = F_ij(k,k) - L_ij(j,k)* F_ij(j,k)
    end do
    !pivot
    F_ij(k,k) = 0.5_dp*F_ij(k,k)/L_ij(k,k)
    end if
  end do
  end subroutine bkwd_diff_F



  !> @brief combines the gradient contributions from the workspace F and those associated with the m_vector and covariance matrix.
  !> @param x is 2D array containing the param_inputs for each thread.
  !> @param n_params is the total number of parameters being evaluated
  subroutine combine_F_grad_K(x,n_params)
  real(dp), dimension(:,:) :: x
  real(dp), dimension(:,:), allocatable :: grad_mu
  real(dp), dimension(:,:,:,:), allocatable:: post_grad_in
  integer, dimension(2) :: shape_arr
  integer, dimension(4) :: shape_arr_4
  integer :: i, j, k, l, n_params


  shape_arr = shape(x)

  if allocated(post_grad_in) then
    deallocate(post_grad_in)
  end if
  allocate(post_grad_in (shape_arr(1), shape_arr(1), shape_arr(1), n_params))
  if allocated(grad_mu) then
    deallocate(grad_mu)
  end if
  allocate(grad_mu( shape_arr(1), n_params ))
  call gp_post_k_grad(x, shape_arr(1), post_grad_in)
  shape_arr_4 = shape(post_grad_in)

  !print*, shape_arr_4(:)
  if allocated(gradient_matrix) then
    deallocate(gradient_matrix)
  end if
  allocate(gradient_matrix (shape_arr_4(3), shape_arr_4(4))  )
  gradient_matrix = 0.0_dp
  do k = 1, shape_arr_4(3)
    do l = 1, shape_arr_4(4)
      do i = 1, shape_arr_4(1)
        do j = 1, shape_arr_4(2)
          gradient_matrix(k,l) = gradient_matrix(k,l) + F_ij(i,j)*post_grad_in(i,j,k,l)
        end do
      end do
    end do
  end do
  !gradient in param space associated with each mean.
  grad_mu = gp_mu_post_grad(x,shape_arr(1))

  gradient_matrix = (gradient_matrix - grad_mu)


  !deallocate(F_ij, L_ij)
  end subroutine combine_F_grad_K


  !> @brief runs the gradient estimator associated routines from one block.
  !> @param x is 2D array containing the param_inputs for each thread.
  !> @param n_params is the total number of parameters being evaluated.
  !> @param z_in is a 1D array of normally distributed numbers.
  !> @param previous_best is a real number containing the best energy from the prior trial.
  function run_gradient_estimator(x, z_in, n_params, previous_best)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), allocatable :: run_gradient_estimator
    real(dp), dimension(:), intent(in) :: z_in
    integer, intent(in) :: n_params
    real(dp), intent(in) :: previous_best

    call init_diff_F(x,z_in, previous_best)
    call bkwd_diff_F()
    call combine_F_grad_K(x, n_params)
    if allocated(run_gradient_estimator) then
      deallocate(run_gradient_estimator)
    end if
    allocate(run_gradient_estimator, mold=gradient_matrix)
    run_gradient_estimator = gradient_matrix*check_max
    !deallocate(gradient_matrix)
  end function run_gradient_estimator

end module gradient_estimator

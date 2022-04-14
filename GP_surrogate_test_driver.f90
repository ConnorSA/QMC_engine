!> @brief Driver for testing the Gaussain process surrogate

program main
    use priors
    use shared_constants
    use gp_surrogate
    implicit none
    real(dp) :: ker_var=0.5_dp, ker_length=1.0_dp
    integer, parameter :: n_data_in=2, n_dof_in=2
    real(dp), dimension(n_data_in,n_dof_in) :: param_data_in
    real(dp), dimension(n_data_in) :: E_data_in
    real(dp), dimension(1,n_dof_in) :: test_point = reshape((/1.0_dp, 1.0_dp/), (/1,n_dof_in/)), &
    more_data = reshape((/0.5_dp, 1.0_dp/), (/1,n_dof_in/))
    real(dp), dimension(1) :: more_E_data = 1.5_dp
    real(dp), dimension(2,2,2,n_dof_in) :: post_grad
    integer :: i,j

    do i=1,n_data_in
        do j=1,n_dof_in
            param_data_in(i,j) = real(i*j)
        end do
        E_data_in(i) = real(i)*0.5
    end do

    call gp_init(prior_mean, prior_mean_dx,  ker_var, ker_length, param_data_in, E_data_in, n_data_in, n_dof_in)
    print*,'init'
    print*, gp_k_post(test_point, 1)
    print*, gp_mu_post(test_point, 1)
    call gp_update(more_data, more_E_data, 1, 'U')
    print*,'update'
    print*, gp_k_post(test_point, 1)
    print*, gp_mu_post(test_point, 1)
    print*, 'testing with 1st data'
    print*, gp_k_post(param_data_in,n_data_in)
    print*, gp_mu_post(param_data_in,n_data_in)
    call gp_post_k_grad(param_data_in, n_data_in, post_grad)
    print*,'k_gradded'
    print*,post_grad
    print*,'mu_grad'
    print*,gp_mu_post_grad(param_data_in,n_data_in)


end program main

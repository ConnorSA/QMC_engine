!> @brief Gaussian process surrogate submodules and functions

module gp_surrogate
    use shared_constants
    !$ use omp_lib
    implicit none
    ! current state to be stored as private variariables/subprograms

    logical, protected :: init = .false.
    integer, protected :: n_data, n_dof,n_threads
    real(dp), protected :: kernel_var, kernal_inv_length
    real(dp), dimension(:,:), allocatable, protected :: param_data
    real(dp), dimension(:), allocatable, protected :: E_data
    real(dp), dimension(:,:), allocatable, protected :: param_pres, param_cov
    real(dp), dimension(:), allocatable, protected :: data_mean !predicted mean of data (by prior)

    procedure (mean_func_interface),pointer :: mu_prior => null()
    procedure (mean_func_interface_dx),pointer :: mu_prior_dx => null()
    procedure (cov_kernal_interface),pointer :: k_prior => null()
    procedure (cov_kernal_dx_1_interface),pointer :: k_prior_dx_1_i => null()
    procedure (cov_kernal_xx_dx_interface),pointer :: k_prior_xx_dx_i => null()
    procedure (cov_kernal_xx_dx_interface),pointer :: gp_k_post_xx_d_x_i => null()


    abstract interface
        function mean_func_interface(x) result(out)
            use shared_constants
            implicit none
            real(dp), dimension(:), intent(in) :: x !param space point
            real(dp) :: out
        end function  mean_func_interface
    end interface

    abstract interface
    function mean_func_interface_dx(x,dim) result(out)
        use shared_constants
        implicit none
        integer, intent(in) :: dim !dimension along which differentiation occurs
        real(dp), dimension(:), intent(in) :: x !param space point
        real(dp) :: out
    end function  mean_func_interface_dx
    end interface

    !note this means prior covariance needs only to be pointwise evaluatable
    abstract interface
        function cov_kernal_interface(x_1,x_2) result(out)
            use shared_constants
            implicit none
            real(dp), dimension(:), intent(in) :: x_1 !param space point
            real(dp), dimension(:), intent(in) :: x_2 !param space point
            real(dp) :: out
        end function  cov_kernal_interface
    end interface

    !interfaces for derivatives of covariances
    abstract interface
        function cov_kernal_dx_1_interface(x_1,x_2,dim) result(out)
            use shared_constants
            implicit none
            integer, intent(in) :: dim
            real(dp), dimension(:), intent(in) :: x_1 !param space point
            real(dp), dimension(:), intent(in) :: x_2 !param space point
            real(dp) :: out
        end function  cov_kernal_dx_1_interface
    end interface

    abstract interface
        function cov_kernal_xx_dx_interface(x,dim) result(out)
            use shared_constants
            implicit none
            integer, intent(in) :: dim
            real(dp), dimension(:), intent(in) :: x !param space point
            real(dp) :: out
        end function  cov_kernal_xx_dx_interface
    end interface

    !use overloaded functions unless otherwise needed
    !>@brief function for postieror covairance kernal
    !>use in format x_1, x_2, x1_dim, x2_dim for cov(x_1,x_2)
    !>use in format x, x_dim for cov(x,x)
    interface gp_k_post
        module procedure gp_k_post_diff_x,gp_k_post_same_x
    end interface

    !>@brief general initialisation function, see gp_init_gausscov
    interface gp_init
        module procedure gp_init_arbcov, gp_init_gausscov
    end interface

    private mu_prior, k_prior, n_data, n_dof, param_data, E_data, param_pres, param_cov, data_mean, init, gaussian_kernel, &
     kernel_var, kernal_inv_length, gaussian_kernel_dx_1_i, k_prior_dx_1_i, k_prior_xx_dx_i, zero_func, gp_k_post_d_x_1_i, &
     gp_k_post_xx_d_x_i, svd_inverse, gp_update_full_svd, gp_update_woodbury_block_update,n_threads
contains

        !>@brief a gaussian kernal with variance kernel_var, and lengthsacle 1/kernal_inv_length, both module variables
    function gaussian_kernel(x_1,x_2) result(out)
        implicit none
        real(dp), dimension(:), intent(in) :: x_1 !>param space point
        real(dp), dimension(:), intent(in) :: x_2 !>param space point
        real(dp) :: out

        out =abs(kernel_var*exp(-norm2(x_1-x_2)**2*0.5*kernal_inv_length**2))

    end function  gaussian_kernel

    !>@brief derivative of gaussian_kernel wrt dimension dim, of x_1
    function gaussian_kernel_dx_1_i(x_1,x_2,dim) result(out)
        implicit none
        real(dp), dimension(:), intent(in) :: x_1 !>param space point
        real(dp), dimension(:), intent(in) :: x_2 !>param space point
        integer, intent(in) :: dim !>dimension along which derivative occurs
        real(dp) :: out

        out = -kernal_inv_length*(x_1(dim)-x_2(dim))*gaussian_kernel(x_1,x_2)

    end function  gaussian_kernel_dx_1_i

    !>@brief returns 0, used for some derivatives
    function zero_func(x,dim) result(out)
        use shared_constants
        implicit none
        real(dp), dimension(:), intent(in) :: x !param space point
        integer, intent(in) :: dim
        real(dp) :: out

        out=0.0_dp

    end function  zero_func

    !>@brief intialises with a gaussian covariance
    subroutine gp_init_gausscov(mean_prior, mean_prior_dx, ker_var, ker_length, param_data_in, E_data_in, n_data_in,&
        n_dof_in,n_threads_in)
        implicit none
        procedure(mean_func_interface) :: mean_prior !>mean_pior, following mean_func_interface form
        procedure(mean_func_interface_dx) :: mean_prior_dx !>derivative of mean_pior, following mean_func_interface_dx form
        integer, intent(in) :: n_data_in, n_dof_in,n_threads_in  !>number of data points in param_data_in, number of dimensions of data, number of threads assigned
        real(dp), intent(in) :: ker_var, ker_length !>descriptors for rbf/gaussian cov kernal
        real(dp), dimension(n_data_in,n_dof_in), intent(in) :: param_data_in !>initialation data parameter space locations
        real(dp), dimension(n_data_in), intent(in) :: E_data_in !>initialation data energies
        integer :: i,j

        n_threads=n_threads_in

        !setting prior parameters
        kernel_var = ker_var
        kernal_inv_length = 1.0_dp/ker_length
        !setting priors
        mu_prior => mean_prior
        k_prior => gaussian_kernel

        !derivatives of priors, needed for SGA step later
        mu_prior_dx => mean_prior_dx
        k_prior_dx_1_i => gaussian_kernel_dx_1_i
        k_prior_xx_dx_i => zero_func
        gp_k_post_xx_d_x_i => zero_func !prior 0=>post 0, this saves having to evaluate (or write, until support is expanded) code to get the 0

        !setting data
        n_data = n_data_in
        n_dof = n_dof_in
        allocate(param_data(n_data,n_dof))
        allocate(E_data(n_data))
        E_data = E_data_in
        param_data = param_data_in


        allocate(param_cov(n_data,n_data))
        allocate(param_pres(n_data,n_data))

        !finds covarience for new data
        do i=1,n_data
            do j=1,i
                !is symmteric, this reduces calls
                param_cov(i,j) = k_prior(param_data(i,:), param_data(j,:))
                param_cov(j,i) = param_cov(i,j)
            end do
        end do

        !compute inverse of cov (precsion)
        param_pres = param_cov
        call svd_inverse(param_pres,  n_data)

        allocate(data_mean(n_data))
        !updating the prior mean list
        do i=1,n_data
            data_mean(i) = mu_prior(param_data(i,:))
        end do

        init = .True.

    end subroutine gp_init_gausscov

    !>do not use, don't have a function for the dervatives of the kernal/mean
    subroutine gp_init_arbcov(mean_prior,  cov_prior, param_data_in, E_data_in, n_data_in, n_dof_in)
        implicit none
        procedure(mean_func_interface) :: mean_prior
        procedure(cov_kernal_interface) :: cov_prior
        integer, intent(in) :: n_data_in, n_dof_in
        real(dp), dimension(n_data_in,n_dof_in), intent(in) :: param_data_in
        real(dp), dimension(n_data_in), intent(in) :: E_data_in
        !work and ipiv are needed for lapack call, have no useful info
        integer :: i,j

        !setting priors
        mu_prior => mean_prior
        k_prior => cov_prior

        !setting data
        n_data = n_data_in
        n_dof = n_dof_in
        allocate(param_data(n_data,n_dof))
        allocate(E_data(n_data))
        E_data = E_data_in
        param_data = param_data_in


        allocate(param_cov(n_data,n_data))
        allocate(param_pres(n_data,n_data))

        !finds covarience for new data
        do i=1,n_data
            do j=1,i
                !is symmteric, this reduces calls
                param_cov(i,j) = k_prior(param_data(i,:), param_data(j,:))
                param_cov(j,i) = param_cov(i,j)
            end do
        end do

        !compute inverse of cov (precsion)
        param_pres = param_cov
        call svd_inverse(param_pres,  n_data)

        allocate(data_mean(n_data))
        !updating the prior mean list
        do i=1,n_data
            data_mean(i) = mu_prior(param_data(i,:))
        end do

        init = .True.

    end subroutine gp_init_arbcov

    !>@brief posterior mean
    function gp_mu_post(x,x_dim) result(out)
        implicit none
        integer, intent(in) :: x_dim !>dimesnsions of input
        real(dp), dimension(x_dim,n_dof), intent(in) :: x !>param space points
        real(dp), dimension(x_dim) :: out
        real(dp), dimension(x_dim) :: x_means
        real(dp), dimension(x_dim,n_data) :: x_cov_data
        integer :: i,j

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        do i=1,x_dim
            x_means(i) = mu_prior(x(i,:))
            do j=1,n_data
                x_cov_data(i,j) = k_prior(x(i,:), param_data(j,:))
            end do
        end do

       out = x_means + matmul(matmul(x_cov_data,param_pres),(E_data-data_mean))

    end function gp_mu_post

    !>@brief derivative of the posterior mean, wrt dimension dim
    function gp_mu_post_dx(x, dim) result(out)
        implicit none
        !only gets evaluated at one point
        integer, intent(in) :: dim
        real(dp), dimension(n_dof), intent(in) :: x !>param space point
        real(dp) :: out
        real(dp), dimension(1,n_data) :: x_cov_data_dx
        real(dp), dimension(1) :: prod_temp
        integer :: j

        do j=1,n_data
            x_cov_data_dx(1,j) = k_prior_dx_1_i(x,param_data(j,:),dim)
        end do

        prod_temp = matmul(x_cov_data_dx,matmul(param_pres,(E_data-data_mean)))

        out = mu_prior_dx(x, dim) + prod_temp(1)

    end function gp_mu_post_dx

    !>@brief derivative of the posterior mean, needed for stoch_grad
    function gp_mu_post_grad(x, x_dim, only) result(out)
        implicit none
        !row i contain the gradient of mu wrt data_i
        integer, intent(in) :: x_dim
        integer, intent(in), optional :: only !>if used will only evaluate for this data point, all others are set to 0
        real(dp), dimension(x_dim,n_dof), intent(in) :: x !param space point
        real(dp), dimension(x_dim,n_dof) :: out
        integer :: i,j

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        out=0

        if (present(only)) then
            do j=1,n_dof
                out(only,j) = gp_mu_post_dx(x(only,:),j)
            end do
        else
            do i=1,x_dim
                do j=1,n_dof
                    out(i,j) = gp_mu_post_dx(x(i,:),j)
                end do
            end do
        end if

    end function gp_mu_post_grad

    !>@brief specific function for cov(x,x), see gp_k_post
    function gp_k_post_same_x(x,x_dim) result(out)
        implicit none
        integer, intent(in) :: x_dim
        real(dp), dimension(x_dim,n_dof), intent(in) :: x !param space point
        real(dp), dimension(x_dim,x_dim) :: out, x_cov_x
        real(dp), dimension(x_dim,n_data) :: x_cov_data
        real(dp), dimension(n_data,x_dim) :: data_cov_x
        integer :: i,j

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        do i=1,x_dim
            do j=1,n_data
                data_cov_x(j,i) = k_prior(param_data(j,:), x(i,:))
            end do
        end do

        x_cov_data = transpose(data_cov_x)
        do i=1,x_dim
            do j=1,i
                !is symmteric, this reduces calls
                x_cov_x(i,j) = k_prior(x(i,:), x(j,:))
                x_cov_x(j,i) = x_cov_x(i,j)
            end do
        end do

       out = x_cov_x + matmul(matmul(x_cov_data,param_pres),data_cov_x)

    end function gp_k_post_same_x

    !>@brief specific function for cov(x_1,x_2), see gp_k_post
    function gp_k_post_diff_x(x1,x2,x1_dim,x2_dim) result(out)
        implicit none
        integer, intent(in) :: x1_dim, x2_dim
        real(dp), dimension(x1_dim,n_dof), intent(in) :: x1 !param space point
        real(dp), dimension(x2_dim,n_dof), intent(in) :: x2 !param space point
        real(dp), dimension(x1_dim,x2_dim) :: out, x1_cov_x2
        real(dp), dimension(n_data,x2_dim) :: data_cov_x2
        real(dp), dimension(x1_dim,n_data) :: x1_cov_data
        integer :: i,j

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        do i=1,x2_dim
            do j=1,n_data
                data_cov_x2(j,i) = k_prior(x2(i,:), param_data(j,:))
            end do
        end do

        do i=1,x1_dim
            do j=1,n_data
                x1_cov_data(i,j) = k_prior(param_data(j,:),x1(i,:))
            end do
        end do

        do i=1,x1_dim
            do j=1,x2_dim
                x1_cov_x2(i,j) = k_prior(x1(i,:), x2(j,:))
            end do
        end do

       out = x1_cov_x2 + matmul(matmul(x1_cov_data,param_pres),data_cov_x2)

    end function gp_k_post_diff_x

    !>@brief derivative of posterior covariance wrt x_1, dimension dim
    function gp_k_post_d_x_1_i(x1,x2,dim) result(out)
        implicit none
        !only gets evaluated at one point
        integer, intent(in) :: dim
        real(dp), dimension(n_dof), intent(in) :: x1 !param space point
        real(dp), dimension(n_dof), intent(in) :: x2 !param space point
        real(dp) :: out
        real(dp), dimension(n_data,1) :: data_cov_x2
        real(dp), dimension(1,n_data) :: x1_cov_data_dx_1_i
        real(dp), dimension(1,1) :: prod_temp
        integer :: j


        do j=1,n_data
            data_cov_x2(j,1) = k_prior(x2, param_data(j,:))
        end do

        do j=1,n_data
            x1_cov_data_dx_1_i(1,j) = k_prior_dx_1_i(x1,param_data(j,:),dim)
        end do

        prod_temp = matmul(matmul(x1_cov_data_dx_1_i,param_pres),data_cov_x2)

        out = k_prior_dx_1_i(x1,x2,dim) + prod_temp(1,1)

    end function gp_k_post_d_x_1_i

    !>@brief grad of the posterior covariance, evaluated at X of size x_dim,n_dof, needed for stoch_grad
    subroutine gp_post_k_grad(X, x_dim, out)
        !triple nested for loops are slow,try not to call often
        implicit none
        integer,intent(in) :: x_dim !<size of X
        real(dp), dimension(x_dim,n_dof), intent(in) :: X !<point to evaluate at
        real(dp), dimension(x_dim,x_dim,x_dim,n_dof),intent(out) ::  out !< the output matrix, size (x_dim,x_dim,x_dim,n_dof)
        integer :: i,j,l

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        do i=1,x_dim
            do j=1,x_dim
                if (i/=j) then
                    do l=1,n_dof
                        out(i,j,j,l) = gp_k_post_d_x_1_i(X(j,:),X(i,:),l)
                        out(i,j,i,l) = gp_k_post_d_x_1_i(X(i,:),X(j,:),l)
                    end do
                else
                    do l=1,n_dof
                        out(i,j,i,l) = gp_k_post_xx_d_x_i(X(i,:),l)
                    end do
                end if
            end do
        end do

    end subroutine gp_post_k_grad

    !>@brief update routine for adding data to gp
    subroutine gp_update(param_data_in_top, E_data_in_top, n_top, Algo_choice)
        integer, intent(in) :: n_top !>number of data points added
        real(dp), dimension(n_top,n_dof), intent(in) :: param_data_in_top !>nXdof array of points in param space
        real(dp), dimension(n_top), intent(in) :: E_data_in_top !>the corresponding energy
        character(len=1), intent(in) :: Algo_choice !>select inversion algorithm for finding presicion matrix, F for full inverse (from SVD), U for updating inverse (based on block inversion and woodbury matrix update formula)

        select case (Algo_choice)
        case ("F")
            call gp_update_full_svd(param_data_in_top, E_data_in_top, n_top)
        case ("U")
            call gp_update_woodbury_block_update(param_data_in_top, E_data_in_top, n_top)
        case default
            print*, "please make a valid algorithm choice"
        end select

    end subroutine gp_update

    !>@brief see gp_update
    subroutine gp_update_full_svd(param_data_in, E_data_in, n)
        integer, intent(in) :: n !number of data points added
        real(dp), dimension(n,n_dof), intent(in) :: param_data_in !nXdof array of points in param space
        real(dp), dimension(n), intent(in) :: E_data_in !the corresponding energy
        real(dp), dimension(n+n_data,n+n_data) :: data_cov, data_pres
        integer :: i,j
        real(dp), dimension(n+n_data,n_dof) :: param_data_temp
        real(dp), dimension(n+n_data) :: E_data_temp, data_mean_temp


        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        !temmpareraly hold existing data
        E_data_temp(1:n_data)=E_data
        param_data_temp(1:n_data,:)=param_data
        E_data_temp(n_data+1:n+n_data) = E_data_in
        param_data_temp(n_data+1:n+n_data,:) = param_data_in

        data_cov(1:n_data,1:n_data) = param_cov

        !finds covarience for new data
        do i=n_data+1,n_data+n
            do j=1,i
                !is symmteric, this reduces calls
                data_cov(i,j) = k_prior(param_data_temp(i,:), param_data_temp(j,:))
                data_cov(j,i) = data_cov(i,j)
            end do
        end do

        !compute inverse of cov (precsion)
        data_pres = data_cov


        ! do i = 1,n+n_data
        !     data_pres(i,i) = data_pres(i,i) + 10.0_dp**(-4.0_dp)
        ! end do
        call mkl_set_num_threads(n_threads)
        call svd_inverse(data_pres,n+n_data)
        call mkl_set_num_threads(1)
        !updating the prior mean list
        data_mean_temp(1:n_data) = data_mean
        do i=n_data+1,n+n_data
            data_mean_temp(i) = mu_prior(param_data_temp(i,:))
        end do

        !only changing stored values after lapack has confirmed success
        data_mean=0
        deallocate(data_mean)
        allocate(data_mean(n+n_data))
        data_mean=data_mean_temp
        E_data=0
        param_data=0
        deallocate(param_data,E_data)
        allocate(param_data(n+n_data,n_dof))
        allocate(E_data(n+n_data))
        param_data = param_data_temp
        E_data = E_data_temp
        param_cov=0
        param_pres=0
        deallocate(param_cov,param_pres)
        allocate(param_cov(n+n_data,n+n_data))
        allocate(param_pres(n+n_data,n+n_data))
        !stores coveriance and precsion
        param_cov = data_cov
        param_pres = data_pres
        n_data = n_data+n

    end subroutine gp_update_full_svd

    !>@brief see gp_update
    subroutine gp_update_woodbury_block_update(param_data_in, E_data_in, n)
        ! an update scheme based on block inversion and woodbury updates
        ! is less stable than the whole matrix SVD approach
        ! might be an idea for slow decaying kernals though
        integer, intent(in) :: n !number of data points added
        real(dp), dimension(n,n_dof), intent(in) :: param_data_in !nXdof array of points in param space
        real(dp), dimension(n), intent(in) :: E_data_in !the corresponding energy
        real(dp), dimension(n+n_data,n+n_data) :: data_cov, data_pres
        real(dp), dimension(n,n) :: new_data_cov, new_data_pres, inv_block_new
        real(dp), dimension(n_data,n_data) :: old_data_cov, old_data_pres, inv_block_old
        real(dp), dimension(n_data,n) :: old_new_data_cov
        integer :: i,j
        real(dp), dimension(n+n_data,n_dof) :: param_data_temp
        real(dp), dimension(n+n_data) :: E_data_temp, data_mean_temp

        if (.not. init) then
            print*,'not intialised'
            stop
        end if

        !temmpareraly hold existing data
        E_data_temp(1:n_data)=E_data
        param_data_temp(1:n_data,:)=param_data
        E_data_temp(n_data+1:n+n_data) = E_data_in
        param_data_temp(n_data+1:n+n_data,:) = param_data_in

        data_cov(1:n_data,1:n_data) = param_cov
        !finds covarience for new data
        do i=n_data+1,n_data+n
            do j=1,i
                !is symmteric, this reduces calls
                data_cov(i,j) = k_prior(param_data_temp(i,:), param_data_temp(j,:))
                data_cov(j,i) = data_cov(i,j)
            end do
        end do

        !finds the blocks of the matricies
        old_data_cov = param_cov
        old_data_pres = param_pres
        new_data_cov = data_cov(n_data+1:n+n_data,n_data+1:n+n_data)
        old_new_data_cov = data_cov(1:n_data,n_data+1:n_data+n)

        !finding the block inverse
        ! new_data_pres = new_data_cov
        ! call svd_inverse(new_data_pres,n)

        !setting the diagonal blocks of the larger matrix
        inv_block_new = new_data_cov - matmul(transpose(old_new_data_cov),matmul(old_data_pres,old_new_data_cov))

        ! do i = 1,n
        !     inv_block_new(i,i) = inv_block_new(i,i) + 10.0_dp**(-4.0_dp)
        ! end do
        call svd_inverse(inv_block_new,n)
        !via woodbury
        inv_block_old = matmul(old_new_data_cov,matmul(inv_block_new,transpose(old_new_data_cov)))
        inv_block_old = old_data_pres + matmul(old_data_pres,matmul(inv_block_old,old_data_pres))

        !combining
        data_pres(1:n_data,1:n_data) = inv_block_old
        data_pres(n_data+1:n_data+n,n_data+1:n_data+n) = inv_block_new
        data_pres(1:n_data,n_data+1:n_data+n) = -matmul(old_data_pres,matmul(old_new_data_cov,inv_block_new))
        data_pres(n_data+1:n_data+n,1:n_data) = -matmul(inv_block_new,matmul(transpose(old_new_data_cov),old_data_pres))

        !updating the prior mean list
        data_mean_temp(1:n_data) = data_mean
        do i=n_data+1,n+n_data
            data_mean_temp(i) = mu_prior(param_data_temp(i,:))
        end do

        !only changing stored values after lapack has confirmed success
        data_mean=0
        deallocate(data_mean)
        allocate(data_mean(n+n_data))
        data_mean=data_mean_temp
        E_data=0
        param_data=0
        deallocate(param_data,E_data)
        allocate(param_data(n+n_data,n_dof))
        allocate(E_data(n+n_data))
        param_data = param_data_temp
        E_data = E_data_temp
        param_cov=0
        param_pres=0
        deallocate(param_cov,param_pres)
        allocate(param_cov(n+n_data,n+n_data))
        allocate(param_pres(n+n_data,n+n_data))
        !stores coveriance and precsion
        param_cov = data_cov
        param_pres = data_pres
        n_data = n_data+n

        !cleaning up some allocated variables
    end subroutine gp_update_woodbury_block_update

    !>@brief prints all scalar stae variables, size of all array state variables, and min and max stored energy
    subroutine gp_debug_print_state()
    implicit none

    print*, "init=", init
    print*, "n_data=", n_data
    print*, "n_dof=", n_dof
    print*, "kernel_var=", kernel_var
    print*, "kernal_inv_length=", kernal_inv_length
    print*, "shape of param_data=", shape(param_data)
    print*, "shape of E_data=", shape(E_data)
    print*, "shape of param_pres=", shape(param_pres)
    print*, "shape of param_cov=", shape(param_cov)
    print*, "shape of data_mean=", shape(data_mean)
    print*, "min energy", minval(E_data)
    print*, "max energy", maxval(E_data)

    end subroutine gp_debug_print_state

    !>@brief deallocates state variables
    subroutine gp_exit()
    implicit none

        !deallocating state variables
        deallocate(param_data,E_data,param_pres, param_cov,data_mean)

        !deallocating pointers
        nullify(mu_prior,mu_prior_dx ,k_prior ,k_prior_dx_1_i ,k_prior_xx_dx_i ,gp_k_post_xx_d_x_i)

    end subroutine gp_exit

    !>@brief internal function for computing matrix inverses based on SVD
    subroutine svd_inverse(data_pres,  mat_size)
        implicit none
        integer, intent(in) :: mat_size
        real(dp), dimension(mat_size,mat_size), intent(inout) :: data_pres
        real(dp), dimension(mat_size,mat_size) :: U,VT,S_inv, data_temp
        real(dp), dimension(mat_size) :: S
        integer, dimension(8*mat_size) :: iwork
        real(dp), dimension(:,:), allocatable :: work
        real(dp) :: max_s, min_neg_eig
        logical :: bad_eigvals, not_psd
        integer :: info, Lwork, i
        min_neg_eig=0.0_dp
        bad_eigvals = .False.
        not_psd = .False.
        allocate(work(1,1))
        call dgesdd('A',mat_size,mat_size,data_temp,mat_size,S,U,mat_size,VT,mat_size,work,-1,iwork,info)
        if (info.ne.0) then
            print*, 'Lapack error (SVD query, dgesdd), info=',info
            stop
        end if
        data_temp=data_pres !resetting
        Lwork = int(work(1,1))
        deallocate(work)
        allocate(work(Lwork,Lwork))
        call dgesdd('A',mat_size,mat_size,data_temp,mat_size,S,U,mat_size,VT,mat_size,work,Lwork,iwork,info)
        if (info.ne.0) then
            print*, 'Lapack error (SVD, dgesdd), info=',info
            stop
        end if
        max_s=maxval(s)
        s_inv = 0
        do i=1,mat_size
            if (abs(S(i))<epsilon(max_s)*max_s*(mat_size)) then
                S_inv(i,i)=0
                bad_eigvals = .True.
            else if (S(i) <0) then
                not_psd = .True.
                S_inv(i,i) = 0!1/s(i)
                min_neg_eig = min(S(i),min_neg_eig)
            else
                S_inv(i,i) = 1/s(i)
            end if
        end do
        if (bad_eigvals) print*, 'Warning, Covariance matrix close to singular'
        if (not_psd) then
            print*,'Matrix not PSD, with min eigenvalue:', min_neg_eig
        end if
        data_pres = matmul(transpose(VT),matmul(S_inv,transpose(U)))

        deallocate(work)

    end subroutine svd_inverse

    !>@brief for finding minimum of the energy and associated parameter space point, from stored data
    subroutine gp_min_stored(min_params,  min_E)
        implicit none
        real(dp), dimension(n_dof), intent(out) :: min_params !> the parameter space point for the min energy
        real(dp), intent(out) :: min_E !> min stored energy
        integer :: minloc_E

        minloc_E = minloc(E_data,1)

        min_E = E_data(minloc_E)
        min_params = param_data(minloc_E,:)

    end subroutine gp_min_stored

    !>@brief returns the integers that control the size of the state variables
    subroutine gp_return_size_data(n_data_out,  n_dof_out)
    implicit none
    integer,intent(out) :: n_data_out,  n_dof_out

    n_data_out = n_data
    n_dof_out = n_dof

    end subroutine gp_return_size_data

    !>for acessing a copy of the state variables, they are returned to the intent(out) parameter with the corresponding name
    subroutine gp_return_state_data(kernel_var_out, kernal_inv_length_out,  param_data_out, E_data_out,&
         param_pres_out, param_cov_out, data_mean_out)
        implicit none
        real(dp), intent(out)  :: kernel_var_out, kernal_inv_length_out
        real(dp), dimension(n_data,n_dof), intent(out):: param_data_out
        real(dp), dimension(n_data), intent(out) :: E_data_out
        real(dp), dimension(n_data,n_data), intent(out) :: param_pres_out, param_cov_out
        real(dp), dimension(n_data), intent(out) :: data_mean_out

        kernel_var_out = kernel_var
        kernal_inv_length_out = kernal_inv_length
        param_data_out = param_data
        E_data_out = E_data
        param_pres_out = param_pres
        param_cov_out = param_cov
        data_mean_out = data_mean

    end subroutine gp_return_state_data

    !>@brief sets a state from inputed data, intended for use with restart files. Variables must correspond to named counterpart
    subroutine gp_restart(n_data_in, n_dof_in,kernel_var_in, kernal_inv_length_in,&
        param_data_in, E_data_in,param_pres_in, param_cov_in, data_mean_in,&
         mean_prior,mean_prior_dx,n_threads_in)
        implicit none
        integer, intent(in) :: n_data_in, n_dof_in,n_threads_in
        real(dp), intent(in) :: kernel_var_in, kernal_inv_length_in
        real(dp), dimension(n_data_in,n_dof_in), intent(in):: param_data_in
        real(dp), dimension(n_data_in), intent(in) :: E_data_in
        real(dp), dimension(n_data_in,n_data_in), intent(in) :: param_pres_in, param_cov_in
        real(dp), dimension(n_data_in), intent(in) :: data_mean_in
        procedure(mean_func_interface) :: mean_prior
        procedure(mean_func_interface_dx) :: mean_prior_dx

        n_data = n_data_in
        n_dof = n_dof_in
        n_threads = n_threads_in

        kernel_var=kernel_var_in
        kernal_inv_length=kernal_inv_length_in
        allocate(param_data(n_data,n_dof))
        param_data=param_data_in
        allocate(E_data(n_data))
        E_data=E_data_in
        allocate(param_pres(n_data,n_data))
        param_pres=param_pres_in
        allocate(param_cov(n_data,n_data))
        param_cov=param_cov_in
        allocate(data_mean(n_data))
        data_mean=data_mean_in


        mu_prior => mean_prior
        k_prior => gaussian_kernel

        !derivatives of priors, needed for SGA step later
        mu_prior_dx => mean_prior_dx
        k_prior_dx_1_i => gaussian_kernel_dx_1_i
        k_prior_xx_dx_i => zero_func
        gp_k_post_xx_d_x_i => zero_func

        init = .True.

    end subroutine gp_restart

end module gp_surrogate

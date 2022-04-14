!> @brief Optimization modules
module stoch_grad
    use shared_constants
    use param_search
    use gradient_estimator
    use basis_functions

    implicit none
    !>performs the stochastic gradient ascent for a given gradient function

    integer :: n_dof,N_restarts
    real(dp) :: best_from_last_run
    real(dp), DIMENSION(:,:,:), ALLOCATABLE :: N_points

    PRIVATE n_dof,gaussian,best_from_last_run
contains

!> @brief Initalises module variables and sets up starts points for the restarts of the gradient assent
!> @param N sets number of restarts
!> @param best_from_last_run_in gives the best energy from the last set of points
!> @param seed random seed in
subroutine stoch_grad_init(N,best_from_last_run_in,seed)
    implicit none

    INTEGER,INTENT(IN) :: N
    real(dp), intent(in) :: best_from_last_run_in
    integer, dimension(:), INTENT(IN):: seed
    integer :: i

    n_dof = size(dof_bounds,1)
    N_restarts = N
    ALLOCATE(N_points(N_restarts,N_restarts,n_dof))

    best_from_last_run=best_from_last_run_in

    do i=1,N_restarts
        call latin_hypercube(N_points(i,:,:),dof_bounds, N_restarts ,seed(1))
    end do

end subroutine stoch_grad_init

!>deallocates state variables
subroutine stoch_grad_exit()
    implicit none

    deALLOCATE(N_points)

end subroutine stoch_grad_exit

!> Performs the gradient assent for sequence length
!> X_0 is the intial point to start the assent at
!> sequence_length is how long the sequence to average over is
!> No_samlples is the number of samples taken in the gradient estimator,
!> gamma is parameter to tune for how much to move the next point long the estimated gradient
!> out is the output and used to compute the sum on the go

subroutine grad_accent(X_0,sequence_length,no_samples,gamma,out)

   implicit none
   real(dp), dimension(:,:),intent(in) :: X_0
   real(dp), dimension(:,:), intent(out) :: out
   real(dp), dimension(:,:), allocatable :: X_t,X_t_prev
   real(dp) :: gamma
   INTEGER::i,sequence_length,no_samples

   ALLOCATE(X_t(size(X_0,1),size(X_0,2)))
   ALLOCATE(X_t_prev(size(X_0,1),size(X_0,2)))


   X_t_prev=X_0
   out = X_t_prev
   do i = 1,sequence_length

    X_t= nearest_point(dof_bounds,(X_t_prev + gamma*G_t(X_t_prev,no_samples,shape(X_t))))
    out = out + X_t !>as result is taken as average, computing on the go uses less memory
    X_t_prev = X_t

   end do

   out = out/real(sequence_length+1)

end subroutine grad_accent

!> @brief Returns a standard gaussian mean = 0, var = 1  of given dimension
!> @param density_dimension is the size of gaussian you want
function gaussian(density_dimension) result(x_prop)
    implicit none
    integer, intent(in) :: density_dimension
    real(dp), dimension(2*((density_dimension+1)/2)) :: u !>uniform random sample,
    !>note the dimension is the nearest even number above the given dimension
    !>integer division will autmatically round towards 0
    real(dp), dimension(density_dimension) :: xi !>normal random sample
    real(dp), dimension(density_dimension) :: x_prop
    integer :: i

    !>implementing box muller to generate a normal
    call random_number(u)
    do i=1,density_dimension/2
        xi(i) = sqrt(-2.0_dp*log(u(i)))*cos(2.0_dp*pi*u(i*2))
        xi(i*2) = sqrt(-2.0_dp*log(u(i)))*sin(2.0_dp*pi*u(i*2))
    end do
    if (modulo(density_dimension,2).ne.0) then
        xi(density_dimension) = sqrt(-2.0_dp*log(u(density_dimension)))*cos(2.0_dp*pi*u(density_dimension+1))
    end if

    !>basic GRW
    x_prop = xi

end function gaussian

!> @brief Return an estimate for the gradient which is an average over no_samples
!> @param X initial input
!> @param no_samples how many samples you want to use from the gradient estimator
!> @param X_shape defines the shape of the output
function G_t (X,no_samples,X_shape)

    implicit none
    INTEGER:: no_samples
    real(dp), dimension(:,:), intent(in) :: X
    integer, dimension(2), intent(in) :: X_shape
    real(dp), DIMENSION(X_shape(1),X_shape(2)) :: G_t
    INTEGER :: i

    G_t=X

    do i = 1,no_samples

        G_t = G_t + run_gradient_estimator(G_t,gaussian(size(X,1)),n_dof,best_from_last_run)

    end do

    G_t = G_t/no_samples

end function G_t

!> @brief Returns a new point back to the search space if it leaves
!> @param X new suggested point
!> @param bounds the degree of freedom bounds defined on input
function nearest_point(bounds,X) result(X_Out)

    implicit none
    real(dp), dimension(:,:), intent(in) :: bounds
    real(dp), dimension(:,:),INTENT(IN):: X
    real(dp), dimension(:,:),ALLOCATABLE:: X_Out
    INTEGER :: i,j

    ALLOCATE(X_Out(size(X,1),size(X,2)))

    do i = 1,size(X,1)

        do j = 1 ,size(X,2)

            if (X(i,j)<bounds(j,1)) then
                X_Out(i,j)=bounds(j,1)

            else if (X(i,j)>bounds(j,2)) then
                X_Out(i,j)=bounds(j,2)
            else
                X_Out(i,j) = X(i,j)
            end if

        end do



    end do
end function nearest_point
end module stoch_grad

!> @brief Functions for obtaining the prior mean and its derivative 
module priors
    use shared_constants
    implicit none

contains
    !> @brief Gives the prior mean
    function prior_mean(x) result(out)
        implicit none
        real(dp), dimension(:), intent(in) :: x !param space point
        real(dp) :: out

        out = 1.0_dp

    end function prior_mean

    !> @brief Derivative of the prior mean
    function prior_mean_dx(x,dim) result(out)
        implicit none
        integer, intent(in) :: dim
        real(dp), dimension(:), intent(in) :: x !param space point
        real(dp) :: out

        out = 0.0_dp

    end function prior_mean_dx

end module priors

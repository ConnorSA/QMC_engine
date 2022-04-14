!> @brief Function completing calculations of either the wavefunction or electron density (depending on the system) at a set of coordinates [x,y,0] on a grid of equally spaced points on a xy grid.

MODULE calculations

  USE ISO_FORTRAN_ENV
  !> Uses the modules \ref component_functions , \ref basis_functions and \ref electron_density_functions
  USE component_functions
  USE basis_functions
  USE electron_density_functions

  IMPLICIT NONE

  CONTAINS

  !> @brief calc function returns the wavefunction or electron density evaluated on a 2D xy plane  (taking a slice) depending on the number of electrons in the system.
  !> Sets the z coordinate to 0 as only evaluating on a xy plane.
  FUNCTION calc(dof, points, box, n_ele, n_MCMC_steps)

    IMPLICIT NONE
    INTEGER, PARAMETER :: dp=kind(1.0d0)
    !> @param points is the user defined number of evaluation points for the plotting grid.
    !> @param n_ele is the number of electrons in the system
    !> @param n_MCMC_steps is the number of steps needed for the MCMC used for evaluating the electron density for the 2 electron case.
    INTEGER, INTENT(IN) :: points, n_ele, n_MCMC_steps
    !> @param box is the user defined size of the domain of the system wanting to be looked at.
    real(dp),intent(in) :: box
    !> @param dof is the optimal degrees of freedom
    REAL(dp), DIMENSION(:), INTENT(IN) :: dof
    !> Loop Variables
    INTEGER :: array_len, track, i,j
    REAL(dp) :: a,b
    real(dp),dimension(3,2) :: box_in
    !> @param wave_norm_squared is the result returned from \ref wave_function_normalisation
    real(dp) :: wave_norm_squared
    !> @param calc Ouput array which contains the results of the function
    REAL(dp), DIMENSION(:), ALLOCATABLE :: calc

    do i=1,3
      box_in(i,:) = [-box,box]
    end do

    !> This allocates memory to the array based on the number of grid points defined by the user.
    array_len = (abs(points) + abs(points) + 1)**2
    ALLOCATE(calc(array_len))

    !> Uses \ref wave_function_normalisation which returns the wave_norm_squared
    wave_norm_squared = wave_function_normalisation(box_in,dof,n_MCMC_steps)
    track = 0

    !> This is for the 1 electron case, will just calculate and returns the wavefunction. Uses the \ref wave_function function.
    IF (n_ele .eq. 1) THEN

      DO i = -points, points
        DO j = -points, points
	  !> Track puts the result in the right place in the resulting output calc array.
          track = track + 1
          !> Calculating the positions for input and gives coordinates between -box/2 to box/2
          a =  (real(i) / (abs(points) + abs(points)))*box
          b =  (real(j) / (abs(points) + abs(points)))*box
          !> Returns the result of the wavefunction evaluated at one coordinate point using \ref wave_function function.
          calc(track) = wave_function([a,b,0.0d0], dof)
        END DO
      END DO
       !> Divides the result by the square root of the results from the \ref wave_norm_squared result defined earlier.
      calc = calc/sqrt(wave_norm_squared)

    !> For the 2 electron case, will calculate and return the electron density. Uses the \ref electron_density function
    ELSE IF (n_ele .eq. 2) THEN

      DO i = -points, points
        DO j = -points, points
          !> Calculating the positions for input
          track = track + 1
          a =  (real(i) / (abs(points) + abs(points)))*box
          b =  (real(j) / (abs(points) + abs(points)))*box
          !> Returns the result of the electron density evaluated at one coordinate point using \ref electron_density function.
          calc(track) = electron_density([a,b,0.0d0],box_in, dof, n_MCMC_steps )
        END DO
      END DO
      !> Divides the result by the results from the \ref wave_norm_squared result defined earlier.
      calc = calc/wave_norm_squared

    !> Returns a error print statement if the number of electrons is not defined.
    ELSE
      PRINT*, 'Incorrect number of electrons in output calculation'
    END IF

  END FUNCTION calc


END MODULE calculations

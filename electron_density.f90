!> @brief Electron density function
module electron_density_functions
  use basis_functions

  !$ use omp_lib
  implicit none
  contains
  !> @brief Computes electron density for 2 electron simulations
  !> Uses basic Monte Carlo Integration
  function electron_density(fixed_position,integral_bounds,dof_coefficients,n_MC_points,seed_in)
    real(dp) :: electron_density
    !> @param fixed_position Space position of density
    real(dp), dimension(3), intent(in) :: fixed_position
    !> @param integral_bounds bounds of integration box
    real(dp), dimension(3,2), intent(in) :: integral_bounds
    !> @param dof_coefficients Current dof parameters
    real(dp), dimension(:), intent(in) :: dof_coefficients
    !> @param n_MC_points Number of Monte Carlo Points
    integer, intent(in) :: n_MC_points
    !> @param seed_in input seed
    integer, intent(in),optional :: seed_in
    integer, allocatable :: seed(:)
    integer :: i ! loop variables
    integer :: n
    real(dp) :: x,y,z,Lx,Ly,Lz
    real(dp), dimension(6) :: position_total
    if (present(seed_in)) then
      call random_seed(size=n)
      allocate(seed(n))
      seed = seed_in
      call random_seed(put=seed)
    end if

    do i=1,3
    if ((integral_bounds(i,2)-integral_bounds(i,1))<= 0.1) then
      print*, "Error, box too small or bounds incorrect"
    end if
    end do
    position_total(1:3) = 0.0_dp
    position_total(4:6) = fixed_position
    Lx = integral_bounds(1,2)-integral_bounds(1,1)
    Ly = integral_bounds(2,2)-integral_bounds(2,1)
    Lz = integral_bounds(3,2)-integral_bounds(3,1)
    !$OMP parallel do default(shared) private(i,x,y,z) firstprivate(position_total) reduction(+:electron_density)
    do i=1,n_MC_points
      call random_number(x)
      call random_number(y)
      call random_number(z)
      x = Lx*x +integral_bounds(1,1)
      y = Ly*y +integral_bounds(2,1)
      z = Lz*z +integral_bounds(3,1)
      position_total(4:6) = [x,y,z]
      electron_density = electron_density+wave_function(position_total,dof_coefficients)**2
    end do
    electron_density = 2*Lx*Ly*Lz*electron_density/n_MC_points
    if (allocated(seed)) deallocate(seed)
  end function electron_density

  !> @brief Computes the integral of the wavefunction squared
  !> To normalise the wavefunction for outputing
  !> Uses basic Monte Carlo Integration
  function wave_function_normalisation(integral_bounds,dof_coefficients,n_MC_points,seed_in)
    real(dp) :: wave_function_normalisation
    !> bounds of integration box
    real(dp), dimension(3,2), intent(in) :: integral_bounds
    !> Current dof parameters
    real(dp), dimension(:), intent(in) :: dof_coefficients
    !> Number of Monte Carlo Points
    integer, intent(in) :: n_MC_points
    !> input seed
    integer, intent(in),optional :: seed_in
    integer, allocatable :: seed(:)
    integer :: i ! loop variables
    integer :: n
    real(dp) :: x,y,z,Lx,Ly,Lz
    real(dp), dimension(:),allocatable :: position_total
    if (present(seed_in)) then
      call random_seed(size=n)
      allocate(seed(n))
      seed = seed_in
      call random_seed(put=seed)
    end if

    do i=1,3
    if ((integral_bounds(i,2)-integral_bounds(i,1))<= 0.1) then
      print*, "Error, box too small or bounds incorrect"
    end if
    end do
    allocate(position_total(n_space_dims))
    Lx = integral_bounds(1,2)-integral_bounds(1,1)
    Ly = integral_bounds(2,2)-integral_bounds(2,1)
    Lz = integral_bounds(3,2)-integral_bounds(3,1)
    !$OMP parallel do default(shared) private(i,x,y,z) firstprivate(position_total) reduction(+:wave_function_normalisation)
    do i=1,n_MC_points
      call random_number(x)
      call random_number(y)
      call random_number(z)
      x = Lx*x +integral_bounds(1,1)
      y = Ly*y +integral_bounds(2,1)
      z = Lz*z +integral_bounds(3,1)
      position_total(1:3) = [x,y,z]
      if (n_space_dims == 6) then
        call random_number(x)
        call random_number(y)
        call random_number(z)
        x = Lx*x +integral_bounds(1,1)
        y = Ly*y +integral_bounds(2,1)
        z = Lz*z +integral_bounds(3,1)
        position_total(4:6) = [x,y,z]
      end if
      wave_function_normalisation = wave_function_normalisation&
        +wave_function(position_total,dof_coefficients)**2
    end do
    wave_function_normalisation = Lx*Ly*Lz*wave_function_normalisation/n_MC_points
    if (n_space_dims == 6) then
      wave_function_normalisation = Lx*Ly*Lz*wave_function_normalisation
    end if
    if (allocated(seed)) deallocate(seed)
    deallocate(position_total)
  end function wave_function_normalisation




end module electron_density_functions

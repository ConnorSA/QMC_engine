!> @brief Basic subfunctions
!> This module contains basic functions used in basis_functions.f90
module component_functions
  use shared_constants
  implicit none

  contains

  !> @brief Gaussian distribution centred at the origin
  function centered_gaussian(position, alpha)
    implicit none
    real(dp) :: centered_gaussian
    !> Space Coordinate
    real(dp), dimension(3), intent(in) :: position 
    !> Length scale
    real(dp), intent(in) :: alpha
    centered_gaussian = (2.0_dp * alpha /pi)**(0.75_dp)*exp(-alpha * dot_product(position,position))
  end function centered_gaussian

 !> @brief Analytic Laplacian of Gaussian distribution centred at the origin
  function centered_gaussian_laplacian(position, alpha)
    implicit none
    real(dp) :: centered_gaussian_laplacian
    !> Space Coordinate
    real(dp), dimension(3), intent(in) :: position
    !> Length scale
    real(dp), intent(in) :: alpha

    centered_gaussian_laplacian = &
    (-6.0_dp*alpha + 4.0_dp*alpha**2*dot_product(position,position))*centered_gaussian(position,alpha)
  end function centered_gaussian_laplacian

  !> @brief Slater-1s distribution centred at the origin
  function centered_slater_1s(position, zeta)
    implicit none
    real(dp) :: centered_slater_1s
    !> Space Coordinate
    real(dp), dimension(3), intent(in) :: position
    !> Length scale
    real(dp),intent(in) :: zeta

    centered_slater_1s = ( zeta**3 /pi)**(0.5_dp)*exp(-zeta * norm2(position))
  end function centered_slater_1s

  !> @brief Analytic Laplacian of Slater-1s distribution centred at the origin
  function centered_slater_1s_laplacian(position, zeta)

    implicit none
    real(dp) :: centered_slater_1s_laplacian
    !> Space Coordinate
    real(dp), dimension(3), intent(in) :: position
    !> Length scale
    real(dp), intent(in) :: zeta

    centered_slater_1s_laplacian = &
    (-2.0_dp*zeta/norm2(position) + zeta**2)*centered_slater_1s(position,zeta)
  end function centered_slater_1s_laplacian

  !> Functions that construct the Jastrow interaction function
  !> Notation here follows Schmidt and Moskowitz 1990, refered to as SM90

  !> @brief Jastrow subfuction: correlation range
  !> $\bar{r}$ function in Schmidt and Moskowitz 1990
  function correlation_range(r, d)
    implicit none
    real(dp) :: correlation_range
    !> distance
    real(dp), intent(in) :: r
    !> inverse length scale
    real(dp), intent(in) :: d

    correlation_range = d * r/(1+d*r)
  end function correlation_range

  !> @brief Jastrow subfuction: basic subterm
  !> Subterm of the correlation function, one for each Jastrow dof per atom
  !> This is the term in the k sum in Schmidt and Moskowitz 1990
  function correlation_subterm(r_12, r_I1, r_I2, m, n, o)
    implicit none
    real(dp) :: correlation_subterm
    !> Distance between electrons
    real(dp), intent(in) :: r_12
    !> Distance from atom to electron 1
    real(dp), intent(in) :: r_I1
    !> Distance from atom to electron 2
    real(dp), intent(in) :: r_I2
    !> Triple of integer parameters. This is a row of \ref mno_parameters
    integer,intent(in) :: m, n, o

    real(dp) :: delta ! coefficient to match Schmidt and Moskowitz 1990

    ! Value of delta, following Schmidt and Moskowitz 1990
    if (m==n) then
      delta = 0.5_dp
    else
      delta = 1.0_dp
    end if

    correlation_subterm = delta * (r_I1**m*r_I2**n+r_I1**n*r_I2**m)*r_12**o
  end function correlation_subterm

  !> @brief Jastrow subfuction: atom subterm
  !> Correlation term for each atom. U_I12 from Schmidt and Moskowitz 1990
  function correlation_atom_term(atom_coord,electron_coords,mno_parameters,c,b,d)

    implicit none
    real(dp) :: correlation_atom_term
    !> Space coordinates of atom
    real(dp), dimension(3), intent(in) :: atom_coord
    !> Space positions of electrons
    real(dp), dimension(6), intent(in) :: electron_coords
    !> Jastrow interger parameters
    integer, dimension(:,:), intent(in) :: mno_parameters
    !> Jastrow dofs
    real(dp), dimension(:), intent(in) :: c
    !> Inverse lengthscale of nuclear-electron interaction
    real(dp), intent(in) :: b
    !> Inverse lengthscale of electron-electron interaction
    real(dp), intent(in) :: d

    ! Internal variables
    real(dp) :: r_12 ! Distance between electrons
    real(dp) :: r_I1 ! Distance from atom to electron 1
    real(dp) :: r_I2 ! Distance from atom to electron 2
    integer :: k ! Loop variable

    ! Compute distances
    r_12 = correlation_range(norm2(electron_coords(1:3)-electron_coords(4:6)),d)
    r_I1 = correlation_range(norm2(atom_coord-electron_coords(1:3)),b)
    r_I2 = correlation_range(norm2(atom_coord-electron_coords(4:6)),b)

    ! k sum from Schmidt and Moskowitz 1990 to compute U_I12
    correlation_atom_term = 0.0_dp
    do k = 1, size(mno_parameters,2)
      correlation_atom_term = correlation_atom_term + &
      c(k)*correlation_subterm(r_12,r_I1,r_I2,mno_parameters(1,k),mno_parameters(2,k)&
        ,mno_parameters(3,k))
    end do

  end function correlation_atom_term

  !> @brief Jastrow correlation fuction
  !> Function F from Schmidt and Moskowitz 1990
  function correlation_function(atom_coords,electron_coords,mno_parameters,c,b,d)
    real(dp) :: correlation_function
    !> Coordinates of atoms
    real(dp), dimension(:,:), intent(in) :: atom_coords
    !> Coordinates of electrons
    real(dp), dimension(6), intent(in) :: electron_coords
     !> Jastrow interger parameters
    integer, dimension(:,:), intent(in) :: mno_parameters
    !> Jastrow dofs
    real(dp),dimension(:), intent(in) :: c
    !> Inverse lengthscale of nuclear-electron interaction
    real(dp), intent(in) :: b
    !> Inverse lengthscale of electron-electron interaction
    real(dp), intent(in) :: d

    integer :: i ! Loop variable

    correlation_function = 0.0_dp
    do i = 1, size(atom_coords,2) ! sum over atoms
      correlation_function = correlation_function +&
       correlation_atom_term(atom_coords(:,i),electron_coords,mno_parameters,c,b,d)
    end do
    ! Compute exponential
    correlation_function = exp(correlation_function)
  end function correlation_function
end module component_functions

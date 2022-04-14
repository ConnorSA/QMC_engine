!> @brief Wavefunction, hamiltonian and basis set choice.
!> This module contains functions that define the hamiltonian that specifies the problem and 
!! the basis set being used. 
!> Makes significant use of module variables defined in an initialisation
!! subroutine. This allocates 2 function pointers using \ref wave_function_interface
!> - \ref wave_function 
!> - \ref reduced_hamiltonian
!> @param position a real(dp) array of electron coordinates
!> @param dof_coefficients a real(dp) array of degree of freedom (dof) parameters
!> These are the only functions that should be used outside of the module.
module basis_functions
  use shared_constants
  use component_functions
  implicit none
  
  ! Main functions are pointers assigned by the initialise routine
  ! Main program should only use wave_function and reduced_hamiltonian
  ! wave_function is the quantum wavefunction of the system
  ! reduced_hamiltonian is the hamiltonian divided by the wavefunction $(H\psi)/\psi$
  procedure(wave_function_interface), pointer :: wave_function
  procedure(wave_function_interface), pointer :: reduced_hamiltonian
  
  ! This is a "template" used for the procedure pointers
  abstract interface 
    function wave_function_interface(position, dof_coefficients)
      use shared_constants
      real(dp), dimension(:), intent(in) :: position
      real(dp), dimension(:), intent(in) :: dof_coefficients
      real(dp) :: wave_function_interface
    end function  wave_function_interface 
  end interface

  ! Module Variables. Protected, so can be read but not written by external code. 
  ! All set by the initialisation routine

  logical, protected :: initialised = .false. ! Flag set by initialisation routine, checked by other routines

  ! Variables directly set as inputs to the initialisation routine
  integer, protected :: n_electrons ! Number of electrons
  integer, protected :: n_basis_functions_per_atom ! Number of linear terms in the single-electron wavefunction per atom
  integer, protected :: n_atoms ! Number of atoms
  real(dp), protected, allocatable, dimension(:,:) :: atom_coords ! Coordinates of the atoms. Shape (3,n_atoms)
  integer, protected :: n_Jastrow_dofs ! Number of dofs in the Jastrow (interaction) term
  real(dp), protected :: fd_h ! Lengthscale of the finite difference code
  real(dp), protected :: b_length ! Inverse lengthscale of nuclear-electron interaction
  real(dp), protected :: d_length ! Inverse lengthscale of electron-electron interaction
  real(dp),dimension(:),allocatable, protected :: proton_numbers ! Proton Numbers of atoms
  integer, protected :: basis_type ! integer code for type of basis. In shared_constants module
  ! Variables used outside of this module
  integer, protected :: number_dofs ! Total number of degrees of freedom (dofs) in the wavefunction
  integer, protected :: n_space_dims ! Number of space dimensions, 3*n_electrons, for convenience
  real(dp), allocatable, dimension(:,:), protected :: dof_bounds ! Default bounds for the possible dof values

  ! Variables used internally
  integer, protected :: n_dofs_per_atom ! Number of dofs per atom, for convenience
  integer, protected :: n_dofs_no_Jastrow ! Number of dofs without Jastrow term, for convenience
  integer, protected, allocatable, dimension(:,:) :: mno_parameters ! Parameters that specify the particular choice of Jastrow term
  
  ! individual electron wavefunction for 2 electron code
  procedure(wave_function_interface), pointer :: wave_function_single
  contains
    
  !> @brief Initialisation routine. 
  !> Must be run at the start of the main program, after user input. This defines 
  !! the problem and selects the basis being used. It allocates module variables and procedure pointers. 
  subroutine initialise_basis(n_electrons_in, n_basis_functions_per_atom_in, n_atoms_in, atom_coords_in,&
                              n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in, &
                              proton_numbers_in,basis_type_in)
    implicit none
    ! Required Parameters
    !> Number of electrons
    integer, intent(in) :: n_electrons_in 
    !> Number of linear terms in the single-electron wavefunction per atom
    integer, intent(in) :: n_basis_functions_per_atom_in 
    !> Number of atoms
    integer, intent(in) :: n_atoms_in 
    !> Coordinates of the atoms. Shape (3,n_atoms)
    real(dp), dimension(:,:), intent(in) :: atom_coords_in 
    ! Optional Paramters
    !> Number of dofs in the Jastrow (interaction) term
    integer, optional, intent(in) :: n_Jastrow_in 
    !> Lengthscale of the finite difference code
    real(dp), optional, intent(in) :: fd_length_in 
     !> Inverse lengthscale of nuclear-electron interaction
    real(dp), optional, intent(in) :: Jastrow_b_length_in
     !> Inverse lengthscale of electron-electron interaction
    real(dp), optional, intent(in) :: Jastrow_d_length_in
    !> Proton Numbers of atoms
    real(dp), dimension(:), optional, intent(in) :: proton_numbers_in 
    !> integer code for type of basis. Codes are listed in shared_constants.f90
    integer, optional, intent(in) :: basis_type_in 
    ! Internal Variables
    integer :: i ! Loop variable
    integer :: n_dofs_per_basis_function ! number of dofs per linear term of the wavefunction. Always 2 for slater and Gaussian

    ! Check not already initialised
    if (initialised) then
      print *, "Error, cannot initialise basis twice"
      stop
    end if
    
    if (n_atoms_in <= 0) then
      print*,"Error, number of atoms must be positive"
      stop
    end if
    if ((n_atoms_in .ne. 1).and.(n_atoms_in .ne. 2)) then
      print*,"Error, only 1 or 2 atoms supported"
      stop
    end if
    n_atoms = n_atoms_in

    ! Check and allocate the atom coordinates
    if (size(atom_coords_in,2) .ne. n_atoms_in ) then
      print*,"Error, wrong number of atoms"
      stop
    end if
    if (size(atom_coords_in,1) .ne. 3 ) then
      print*,"Error, wrong dimension of atom coordinates"
      stop
    end if
    allocate(atom_coords(3,size(atom_coords_in,2)))
    atom_coords = atom_coords_in

    ! Check and assign other required inputs
    if (n_basis_functions_per_atom_in <= 0) then
      print*,"Error, number of basis functions per atom must be positive"
      stop
    end if
    n_basis_functions_per_atom = n_basis_functions_per_atom_in

    if ((n_electrons_in .ne. 1).and.(n_electrons_in .ne. 2)) then
      print*,"Error, number of electrons must be 1 or 2" 
      stop
    end if
    n_electrons = n_electrons_in

    n_space_dims = 3*n_electrons
    
    

    ! Check and assign optional paramters
    if (present(n_Jastrow_in)) then
      n_Jastrow_dofs = n_Jastrow_in
    else
      n_Jastrow_dofs = 7 ! Default Value
    end if

    if (present(fd_length_in)) then
      if (fd_length_in>0.00009) then
        fd_h = fd_length_in
      else 
        print*, "Error, finite difference length must be at least 0.0001"
        stop
      end if
    else
      fd_h = 0.01_dp ! Default Value
    end if

    if (present(Jastrow_b_length_in)) then
      if (Jastrow_b_length_in>0.05) then
        b_length = Jastrow_b_length_in
      else 
        print*, "Error, Jastrow b length must be positive"
        stop
      end if
    else
      b_length = 1.0_dp ! Default value
    end if

    if (present(Jastrow_d_length_in)) then
      if (Jastrow_d_length_in>0.05) then
        d_length = Jastrow_d_length_in
      else 
        print*, "Error, Jastrow d length must be positive"
        stop
      end if
    else
      d_length = 1.0_dp ! Default value
    end if
    allocate(proton_numbers(n_atoms))
    if (present(proton_numbers_in)) then
      if (size(proton_numbers_in)==n_atoms) then
        proton_numbers = proton_numbers_in
      else 
        print*, "Error, mismatch in number of atoms and proton numbers"
        stop
      end if
    else
      proton_numbers = 1.0_dp ! Default value is hydrogen
    end if

    if(present(basis_type_in)) then
      basis_type = basis_type_in
    else
      basis_type = slater_1s_code !Default basis set
    end if

    ! Assign individual electron wavefunction with input basis set
    select case(basis_type)
    case (slater_1s_code)
      wave_function_single => wave_function_slater_1s
      reduced_hamiltonian => reduced_hamiltonian_slater_1s
    case (sto_3g_code)
      wave_function_single => wave_function_sto3g
      reduced_hamiltonian => reduced_hamiltonian_sto3g
    case default 
      print *, "Basis type not recognised"
      stop
    end select
    
    ! Set number of dofs per basis function. This may change for future choices of basis sets
    n_dofs_per_basis_function = 2
    n_dofs_per_atom = n_dofs_per_basis_function * n_basis_functions_per_atom 
    n_dofs_no_Jastrow = n_dofs_per_atom * n_atoms 

    ! Assignments for one electron, slater type orbital
    if (n_electrons_in .eq. 1) then
      
      number_dofs = n_dofs_per_atom * n_atoms
      
      allocate(dof_bounds(number_dofs,2))

      ! Assign dof bounds for single electron part
      ! Ordering here is for 2 dofs per basis function, in pairs, first linear dof then width dof
      do i = 1, n_basis_functions_per_atom*n_atoms ! Loop over pairs of dofs
          dof_bounds(2*i-1,:)=[-1.5_dp,1.5_dp] ! Linear coeffs defaults -1.5 to 1.5
          dof_bounds(2*i,:)=[0.1_dp,1.5_dp]   ! width positive, up to 1.5 
      end do

      ! Assign wave_function pointer. Ham pointer already assigned
      wave_function => wave_function_single
    end if

    ! Assignments for 2 electrons, slater type orbital
    if (n_electrons_in .eq. 2) then
      
      ! Allocate the mno parameters that determine the form of Jastrow term
      call mno_allocate(n_Jastrow_dofs)

      number_dofs = n_dofs_per_atom * n_atoms + n_Jastrow_dofs
      
      allocate(dof_bounds(number_dofs,2))

      ! Assign dof bounds
      ! Ordering here is for 2 dofs per basis function, first linear cooef then width coeff
      do i = 1,n_basis_functions_per_atom*n_atoms ! Loop over pairs of dofs
          dof_bounds(2*i-1,:)=[-1.5_dp,1.5_dp] ! Linear coeffs bounds -1.5 to 1.5
          dof_bounds(2*i,:)=[0.1_dp,1.5_dp]   ! width positive, up to 1.5 
      end do
      if (n_Jastrow_dofs.ne.0) then
        do i = n_dofs_per_atom*n_atoms+1, number_dofs
          dof_bounds(i,:) = [-1.5_dp,1.5_dp] ! Jastrow dofs bounds -1.5 to 1.5
        end do
      end if
      ! Assign the procedure pointers
      wave_function => wave_function_2_electrons
      reduced_hamiltonian => reduced_hamiltonian_2_electrons 
    end if

    ! Set initialised flag, which is tested in the routines
    initialised = .true.
  end subroutine initialise_basis

  !> @brief Dinitialisation routine.
  !> MUST be run at the end of main program to deallocate module variables.
  subroutine deinitialise_basis
    implicit none
    if (initialised) deallocate(dof_bounds,atom_coords,proton_numbers)
    if (allocated(mno_parameters)) deallocate(mno_parameters)
    initialised = .false.
  end subroutine deinitialise_basis

  !> @brief Single Electron wavefunction: slater 1s basis.
  !> Uses \ref centered_slater_1s in component_functions.f90
  function wave_function_slater_1s(position,dof_coefficients)
    implicit none
    real(dp) :: wave_function_slater_1s
    !> Space coordinate of electron
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    integer :: i, j ! loop variables
    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.3) then
      print *, "Error, wrong space dimension, got", size(position), "wanted 3"
      stop
    end if
    if (size(dof_coefficients).ne.n_dofs_no_Jastrow) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", n_dofs_no_Jastrow
      stop
    end if
    
    ! Wavefunction is sum over linear terms with 2 dofs each - a linear coefficient and a lengthscale
    wave_function_slater_1s = 0.0_dp 
    do j = 0, n_atoms-1 ! loop over atoms
      do i= 1, n_basis_functions_per_atom ! loop over linear terms in wavefunction
        wave_function_slater_1s = wave_function_slater_1s + &
          dof_coefficients( 2*i-1 + j*n_dofs_per_atom )& ! linear coefficient dof
            *centered_slater_1s(position-atom_coords(:,j+1), & ! position relative to atom
              dof_coefficients( 2*i + j*n_dofs_per_atom  ) ) ! lengthscale dof
      end do
    end do
   
  end function wave_function_slater_1s

  !> @brief Single Electron Reduced Hamiltonian: slater 1s basis.
  !> Reduced Hamiltonian $(H\psi)/\psi$ for 1 electron, with slater 1s type orbital
  !! Uses \ref centered_slater_1s_laplacian in component_functions.f90
  !> Used only for 1 electron problems.
  function reduced_hamiltonian_slater_1s(position,dof_coefficients)
    implicit none
    real(dp) :: reduced_hamiltonian_slater_1s
    !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    integer :: i, j ! Loop variables

    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.n_space_dims) then
      print *, "Error, wrong space dimension, got", size(position), "wanted", n_space_dims
      stop
    end if
    if (size(dof_coefficients).ne.number_dofs) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", number_dofs
      stop
    end if
  
    reduced_hamiltonian_slater_1s = 0.0_dp
    ! Kinetic energy term $-0.5\nabla^2\psi$
    ! This wavefunction is linear, and it is easy to compute this analytically
    do j = 0, n_atoms-1 !> loop over atoms
      do i= 1,n_basis_functions_per_atom !> loop over linear terms in wavefunction
        reduced_hamiltonian_slater_1s = reduced_hamiltonian_slater_1s &
          -0.5_dp * dof_coefficients( 2*i-1 +j*n_dofs_per_atom) & ! linear coefficient dof
            *centered_slater_1s_laplacian(position-atom_coords(:,j+1), & ! position relative to atom
              dof_coefficients( 2*i +j*n_dofs_per_atom ) ) ! lengthscale dof
      end do
    end do
    ! Divide by the wavefunction to get reduced laplacian $(-0.5\nabla^2\psi)/\psi$
    reduced_hamiltonian_slater_1s = reduced_hamiltonian_slater_1s/wave_function_slater_1s(position,dof_coefficients)

    ! Potential energy, already reduced, so no wavefunction in this
    do j = 1, n_atoms ! loop over atoms
    reduced_hamiltonian_slater_1s = reduced_hamiltonian_slater_1s &
      -proton_numbers(j) /norm2(position-atom_coords(:,j)) ! -1/r, r distance from electron to atom
    end do

    if (n_atoms ==2) then
      reduced_hamiltonian_slater_1s = reduced_hamiltonian_slater_1s &
      +proton_numbers(1)*proton_numbers(2) /norm2(atom_coords(:,1)-atom_coords(:,2))
    end if
    
  end function reduced_hamiltonian_slater_1s

  !> @brief Finite Difference Reduced Laplacian
  !> Computes finite difference approximation to the reduced laplacian $(\nabla^2\psi)/\psi$.
  !! Uses that $f''(x)/f(x) =((f(x+h)+f(x-h))/f(x)-2)/(h^2)+O(h^2)$ (*)
  function discrete_Laplacian_reduced(position, h, dofs)
    
    implicit none
    real(dp) :: discrete_Laplacian_reduced
    !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position
    !> Finite difference lengthscale 
    real(dp), intent(in) :: h 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dofs 
    integer :: i ! Loop variable
    real(dp), dimension(:,:), allocatable :: delta ! temporary array for vector displacements
    
    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.n_space_dims) then
      print *, "Error, wrong space dimension, got", size(position), "wanted", n_space_dims
      stop
    end if
    if (size(dofs).ne.number_dofs) then
      print *, "Error, wrong number of dofs, got", size(dofs), "wanted", number_dofs
      stop
    end if
    
    ! Create matrix of cartesian basis vectors of length h. Set as 0 initially, h added during main loop
    allocate(delta(n_space_dims,n_space_dims))
    delta = 0.0_dp ! Matrix assignment
    discrete_Laplacian_reduced = 0.0_dp
    do i = 1, n_space_dims ! Loop over space dimensions
      delta(i,i) = h 
      discrete_Laplacian_reduced = discrete_Laplacian_reduced &
      ! Following (*), sum of (f(x+h)+f(x-h)) in all space dimensions
      + wave_function(position + delta(:,i), dofs) + wave_function(position - delta(:,i), dofs)
    end do 
    ! Following (*), divide by $f(x)$, and subtract 2 for each space dimension
    discrete_Laplacian_reduced = discrete_Laplacian_reduced/wave_function(position ,dofs) - 2*n_space_dims
    ! Following (*), divide by h^2
    discrete_Laplacian_reduced = discrete_Laplacian_reduced/(h**2)

    ! deallocate temporary array
    deallocate(delta)
  end function discrete_Laplacian_reduced

  !> @brief Wavefunction for 2 electrons.
  !> Uses subprocedure \ref correlation_function in component_functions.f90
  !> Uses the function pointer \ref wave_function_single 
  function wave_function_2_electrons(position, dof_coefficients)
    implicit none
    real(dp) :: wave_function_2_electrons
    !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.n_space_dims) then
      print *, "Error, wrong space dimension, got", size(position), "wanted", n_space_dims
      stop
    end if
    if (size(dof_coefficients).ne.number_dofs) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", number_dofs
      stop
    end if

    ! Slater determinant (trivial for 2 electron spin up/down pair)
    ! position(1:3) is coords of first electron, (4:6) is the second. dof(:n_no_jastrow) are the relevant dofs
    wave_function_2_electrons = wave_function_single(position(1:3),dof_coefficients(:n_dofs_no_Jastrow))&
    * wave_function_single(position(4:6),dof_coefficients(:n_dofs_no_Jastrow))
    
    ! Multiply by correlation function. This is the Jastrow term. dof(n_dofs_no_Jastrow+1:) are relevant terms
    if (n_Jastrow_dofs.ne.0) then
    wave_function_2_electrons = &
      correlation_function(atom_coords,position,mno_parameters,dof_coefficients(n_dofs_no_Jastrow+1:),&
      b_length,d_length)*wave_function_2_electrons
    end if
  end function wave_function_2_electrons

  !> @brief Reduced Hamiltonian for 2 electrons.
  !> Computes the reduced Hamiltonian $(H\psi)/\psi$ for 2 electrons.
  !! Uses \ref discrete_Laplacian_reduced in this module.
  function reduced_hamiltonian_2_electrons(position,dof_coefficients)
    implicit none
    real(dp) :: reduced_hamiltonian_2_electrons
     !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    integer :: j ! Loop variable
    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.n_space_dims) then
      print *, "Error, wrong space dimension, got", size(position), "wanted", n_space_dims
      stop
    end if
    if (size(dof_coefficients).ne.number_dofs) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", number_dofs
      stop
    end if

    ! Reduced Kinetic energy and electron-electron reduced potential energy
    ! $-0.5\nabla^2\psi$ +1/r, r distance between electrons
    reduced_hamiltonian_2_electrons = -0.5_dp*discrete_Laplacian_reduced(position, fd_h, dof_coefficients) &
      + 1.0_dp /norm2(position(1:3)-position(4:6))

    ! Reduced electron-nuclear potential energy, -1/r, r distance from electron to atom
    do j = 1, n_atoms ! Loop over atoms
      reduced_hamiltonian_2_electrons  = reduced_hamiltonian_2_electrons  &
          -proton_numbers(j) /norm2(position(1:3)-atom_coords(:,j))&
          -proton_numbers(j) /norm2(position(4:6)-atom_coords(:,j))
    end do

    if (n_atoms ==2) then
      reduced_hamiltonian_2_electrons  = reduced_hamiltonian_2_electrons  &
      +proton_numbers(1)*proton_numbers(2) /norm2(atom_coords(:,1)-atom_coords(:,2))
    end if
  end function reduced_hamiltonian_2_electrons

  !> @brief Allocate parameters for Jastrow factor.
  !> Allocates the paramters that determine type of Jastrow function. Sizes 0,3,7,9 supported.
  !! Functions are from Schmidt and Moskowitz 1990
  subroutine mno_allocate(n_terms)  
    implicit none
    !> Number of terms in Jastrow correlation function
    integer,intent(in) :: n_terms 
    select case (n_terms)
      case (0)
        n_Jastrow_dofs = 0
      case (3)
        allocate(mno_parameters(3,3))
        mno_parameters(:,1) = [0,0,1]
        mno_parameters(:,2) = [0,0,2]
        mno_parameters(:,3) = [2,0,0]
        n_Jastrow_dofs = 3
      case (7)
        allocate(mno_parameters(3,7))
        mno_parameters(:,1) = [0,0,1]
        mno_parameters(:,2) = [0,0,2]
        mno_parameters(:,3) = [0,0,3]
        mno_parameters(:,4) = [0,0,4]
        mno_parameters(:,5) = [2,0,0]
        mno_parameters(:,6) = [3,0,0]
        mno_parameters(:,7) = [4,0,0]
        n_Jastrow_dofs = 7
      case (9)
        allocate(mno_parameters(3,9))
        mno_parameters(:,1) = [0,0,1]
        mno_parameters(:,2) = [0,0,2]
        mno_parameters(:,3) = [0,0,3]
        mno_parameters(:,4) = [0,0,4]
        mno_parameters(:,5) = [2,0,0]
        mno_parameters(:,6) = [3,0,0]
        mno_parameters(:,7) = [4,0,0]
        mno_parameters(:,8) = [2,2,0]
        mno_parameters(:,9) = [2,0,2]
        n_Jastrow_dofs = 9
      case default 
        print *, "Jastrow term of size:", n_terms, "not supported"
        stop
    end select
  end subroutine mno_allocate

  !> @brief Log of the Probability Density.
  !> For the MCMC integration.
  function log_density(position,dof_coefficients) 
    implicit none
    real(dp) :: log_density
    !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    
    log_density = 2*log(abs(wave_function(position,dof_coefficients)))

  end function log_density

  !> @brief Single Electron wavefunction: gaussian sto3g basis.
  !> Uses \ref centered_gaussian in component_functions.f90
  function wave_function_sto3g(position,dof_coefficients)
    implicit none
    real(dp) :: wave_function_sto3g
    !> Space coordinate of electron
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 

    integer :: i, j ! loop variables
    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.3) then
      print *, "Error, wrong space dimension, got", size(position), "wanted 3"
      stop
    end if
    if (size(dof_coefficients).ne.n_dofs_no_Jastrow) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", n_dofs_no_Jastrow
      stop
    end if
    
    ! Wavefunction is sum over linear terms with 2 dofs each - a linear coefficient and a lengthscale
    wave_function_sto3g = 0.0_dp 
    do j = 0, n_atoms-1 ! loop over atoms
      do i= 1, n_basis_functions_per_atom ! loop over linear terms in wavefunction
        wave_function_sto3g = wave_function_sto3g + &
          dof_coefficients( 2*i-1 + j*n_dofs_per_atom )& ! linear coefficient dof
            *centered_gaussian(position-atom_coords(:,j+1), & ! position relative to atom
              dof_coefficients( 2*i + j*n_dofs_per_atom  ) ) ! lengthscale dof
      end do
    end do
   
  end function wave_function_sto3g

  !> @brief Single Electron Reduced Hamiltonian: gaussian sto3g basis.
  !> Reduced Hamiltonian $(H\psi)/\psi$ for 1 electron, with basic gaussian sto3g type orbital
  !! Uses \ref centered_gaussian_laplacian in component_functions.f90
  !> Used only for 1 electron problems.
  function reduced_hamiltonian_sto3g(position,dof_coefficients)
    implicit none
    real(dp) :: reduced_hamiltonian_sto3g
    !> Space coordinate of electrons
    real(dp), dimension(:), intent(in) :: position 
    !> Values of the dofs
    real(dp), dimension(:), intent(in) :: dof_coefficients 
    integer :: i, j ! Loop variables

    ! Initialise and bounds tests
    if (.not.(initialised)) then
      print *, "Error, basis not initialised"
      stop
    end if
    if (size(position).ne.n_space_dims) then
      print *, "Error, wrong space dimension, got", size(position), "wanted", n_space_dims
      stop
    end if
    if (size(dof_coefficients).ne.number_dofs) then
      print *, "Error, wrong number of dofs, got", size(dof_coefficients), "wanted", number_dofs
      stop
    end if
  
    reduced_hamiltonian_sto3g = 0.0_dp
    ! Kinetic energy term $-0.5\nabla^2\psi$
    ! This wavefunction is linear, and it is easy to compute this analytically
    do j = 0, n_atoms-1 ! loop over atoms
      do i= 1,n_basis_functions_per_atom ! loop over linear terms in wavefunction
        reduced_hamiltonian_sto3g = reduced_hamiltonian_sto3g &
          -0.5_dp * dof_coefficients( 2*i-1 +j*n_dofs_per_atom) & ! linear coefficient dof
            *centered_gaussian_laplacian(position-atom_coords(:,j+1), & ! position relative to atom
              dof_coefficients( 2*i +j*n_dofs_per_atom ) ) ! lengthscale dof
      end do
    end do
    ! Divide by the wavefunction to get reduced laplacian $(-0.5\nabla^2\psi)/\psi$
    reduced_hamiltonian_sto3g = reduced_hamiltonian_sto3g/wave_function_sto3g(position,dof_coefficients)

    ! Potential energy, already reduced, so no wavefunction in this
    do j = 0, n_atoms-1 ! loop over atoms
      reduced_hamiltonian_sto3g = reduced_hamiltonian_sto3g &
      -1.0_dp /norm2(position-atom_coords(:,j+1)) !-1/r, r distance from electron to atom
    end do
    
    if (n_atoms ==2) then
      reduced_hamiltonian_sto3g = reduced_hamiltonian_sto3g &
      +proton_numbers(1)*proton_numbers(2) /norm2(atom_coords(:,1)-atom_coords(:,2))
    end if

  end function reduced_hamiltonian_sto3g

end module basis_functions
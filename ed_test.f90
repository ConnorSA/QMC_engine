!> @brief Driver for testing electron density function

program main_driver
  use basis_functions


  use electron_density_functions

  implicit none

  ! User Inputs required - the user must define all of these
  integer :: n_electrons_in ! Number of electrons
  integer :: n_atoms_in ! Number of atoms
  real(dp) :: bond_length ! Bond length for 2 atoms
  real(dp), allocatable, dimension(:,:) :: atom_coords_in ! Coordinates of the atoms. Shape (3,n_atoms)
  integer :: n_trials ! number of trials in the hypercube search
  integer :: n_MCMC_steps ! Total number of steps for each MCMC run

  ! Optional user Inputs - these may be defined by the user, but have good default values below
  integer :: n_basis_functions_per_atom_in ! Number of linear terms in the single-electron wavefunction per atom
  integer :: search_seed ! Seed for the latin hypercube
  integer :: n_Jastrow_in ! Number of dofs in the Jastrow (interaction) term
  real(dp) :: fd_length_in ! Lengthscale of the finite difference code
  real(dp) :: Jastrow_b_length_in ! Inverse lengthscale of nuclear-electron interaction
  real(dp) :: Jastrow_d_length_in ! Inverse lengthscale of electron-electron interaction
  real(dp), dimension(3,2) :: integral_bounds


  real(dp),dimension(7) :: dof


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!! REPLACE FOLLOWING WITH INPUT CODE !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Assign Required Inputs
  n_electrons_in = 2
  n_atoms_in = 2

  ! Can either take in bond length or the actual coordinates.

  bond_length = 2.0_dp

  allocate( atom_coords_in(3,n_atoms_in) )
  if (n_atoms_in == 1) then
    atom_coords_in(:,1)=[0.0_dp,0.0_dp,0.0_dp] ! 1 atom at origin
  else if (n_atoms_in == 2) then
    atom_coords_in(:,1)=[bond_length/2.0_dp,0.0_dp,0.0_dp] !2 atoms evenly placed about origin
    atom_coords_in(:,2)=[-bond_length/2.0_dp,0.0_dp,0.0_dp]
  end if

  n_trials = 10 ! This is just a test value - should be much larger, and scale with number of dofs
  n_MCMC_steps = 10000000 ! This is just a test value, could be larger, try 10^7


  ! Assign optional arguments - these are sensible default values
  search_seed = 132
  n_basis_functions_per_atom_in = 1 ! For slater 1 is sensible. Increase for Gaussians if we implement them
  n_Jastrow_in = 3 ! This is low. Would like to run with 7, but that might be too high for testing
  fd_length_in = 0.1_dp ! Not sure on this
  Jastrow_b_length_in = 1.0_dp ! This is fine, but could experiment with this
  Jastrow_d_length_in = 1.0_dp ! This is fine, but could experiment with this

  integral_bounds(:,1) = -4.0
  integral_bounds(:,2) = 4.0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! END OF INPUT CODE!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! Initialisation of basis and trials

  call initialise_basis(n_electrons_in, n_basis_functions_per_atom_in, n_atoms_in, atom_coords_in,&
                        n_Jastrow_in, fd_length_in, Jastrow_b_length_in, Jastrow_d_length_in)




  dof =[0.11111111111111116_dp ,      0.60000000000000009_dp    ,&
     0.33333333333333337_dp   ,    0.90000000000000002_dp   ,    0.83333333333333337_dp   ,  &
       0.16666666666666674_dp   ,   -0.50000000000000000_dp]
  print *, electron_density([2.0_dp,0.1_dp,0.1_dp],integral_bounds,dof,n_MCMC_steps)

  ! deallocation
  deallocate(atom_coords_in)
  call deinitialise_basis



end program main_driver

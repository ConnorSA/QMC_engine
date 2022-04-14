!> @brief Contains all the subroutines related to writing results to file.
!> This module contains all the subroutines needed for all the output result files which will be used for plotting and restarting.

MODULE write_file

  USE ISO_FORTRAN_ENV
  USE netcdf

  IMPLICIT NONE

  CONTAINS

  !> @brief Main results writing to NetCDF file.
  !> The result_netcdf subroutine outputs the main NetCDF result file containing the optimal degrees of freedom, either the electron density or wavefunction and the number of electrons and nuclei used.
  !> Saves the file as results.nc4 which is used
  !> Has error statements printed out at each section in order to aid with debugging.

  SUBROUTINE result_netcdf(dof, ele, num_ele, num_nuc)

    INTEGER, PARAMETER :: dp=kind(1.0d0)
    !> @param filename which outputs a NetCDF file names results.nc4
    CHARACTER(LEN=100) :: filename = 'results.nc4'
    !> @param dof the optimal degrees of freedom
    REAL(dp), DIMENSION(:), INTENT(IN) :: dof
    !> @param ele an rank 1 array containing either the wavefunction or electron density results evaluated on a grid op coordinates (for 1 or 2 electron system respectively).
    REAL(dp), DIMENSION(:), INTENT(IN) :: ele
    !> @param num_ele the number of electrons used in the calculation
    !> @param num_nuc the number of nuclei used in the calculation
    INTEGER, INTENT(IN) :: num_ele, num_nuc

    !> Dimension names for the arrays.
    CHARACTER(LEN=100), DIMENSION(1) :: dim_dof = (/"Optimal_DOF"/)
    CHARACTER(LEN=100), DIMENSION(1) :: dim_ele = (/"Electron_Density"/)

    !> ID's
    ! File ID's
    INTEGER :: ierr, file_id
    ! Variable ID's
    INTEGER :: var_dof,  var_ele, var_num_ele, var_num_nuc
    ! Dimension ID's for the arrays
    INTEGER, DIMENSION(1) :: did_dof, did_ele
    ! Implicit types for the sizes for the arrays
    INTEGER, DIMENSION(1) :: size_dof, size_ele

    ! Create the file
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! Variable sizes
    size_dof = SHAPE(dof)
    size_ele = SHAPE(ele)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Inputting the data

    ! Number of electrons
    ierr = nf90_def_var(file_id, "Num_of_Electrons", NF90_DOUBLE, var_num_ele)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! Number of nuclei
    ierr = nf90_def_var(file_id, "Num_of_Nuclei", NF90_DOUBLE, var_num_nuc)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! DOF
    ierr = nf90_def_dim(file_id, dim_dof(1), size_dof(1), did_dof(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "Optimal_DOF", NF90_DOUBLE, did_dof, var_dof)

    ! Electron Density
    ierr = nf90_def_dim(file_id, dim_ele(1), size_ele(1), did_ele(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "Electron_Density", NF90_DOUBLE, did_ele, var_ele)


    ! Metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    !> Writing the data to file

    ! Number of electrons
    ierr = nf90_put_var(file_id, var_num_ele, num_ele)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! Number of nuclei
    ierr = nf90_put_var(file_id, var_num_nuc, num_nuc)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! DOF
    ierr = nf90_put_var(file_id, var_dof, dof)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! Electron Density
    ierr = nf90_put_var(file_id, var_ele, ele)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    !> Closing the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !Finishing statement if writing file is successful.
    PRINT *, "Success in writing the NETCDF file: ", filename

  END SUBROUTINE result_netcdf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Output results for bond lengths and corresponding energy used for plotting.
  !> Will only be used for the 2 nuclei case.

  ! Output file for the energies vs bond length plotting
  SUBROUTINE energies_netcdf(energies, bondlen)

    INTEGER, PARAMETER :: dp=kind(1.0d0)
    ! Inputting variables
    !> @param filename outputs a NetCDF file named bond_length.nc4
    CHARACTER(LEN=100) :: filename = 'bond_length.nc4'
    REAL(dp), DIMENSION(:), INTENT(IN) :: energies, bondlen

    !> Dimension names
    CHARACTER(LEN=100), DIMENSION(1) :: dim_e = (/"Energy"/)
    CHARACTER(LEN=100), DIMENSION(1) :: dim_b = (/"Bond_Lengths"/)

    !> ID's
    ! File IDs
    INTEGER :: ierr, file_id
    ! Variable ID's
    INTEGER :: var_e, var_b
    ! Dimension ID's for the arrays
    INTEGER, DIMENSION(1) :: did_e, did_b
    ! Implicit types for the array sizes
    INTEGER, DIMENSION(1) :: size_e, size_b

    ! Create the file
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ! Variable sizes
    size_e = SHAPE(energies)
    size_b = SHAPE(bondlen)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Inputting the data
    ! Energies
    ierr = nf90_def_dim(file_id, dim_e(1), size_e(1), did_e(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "Energy", NF90_DOUBLE, did_e, var_e)

    ! Bond length
    ierr = nf90_def_dim(file_id, dim_b(1), size_b(1), did_b(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "Bond_length", NF90_DOUBLE, did_b, var_b)

      ! Metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Writing data to file
    ! Energy
    ierr = nf90_put_var(file_id, var_e, energies)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_b, bondlen)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Closing the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !Finishing statement if writing file is successful.
    PRINT *, "Success in writing the NETCDF file: ", filename

  END SUBROUTINE energies_netcdf

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Writes the xyz coordinates used in the calculation of either the wavefunction or electron density to a text file names xyz.txt.
  !> This only produces a equally spaced square grid of points in the x-y plane, with the z coordinate being set to zero.
  !> These coordinates are used for the 2D contour plot in the resulting output plot.

  ! Outputs 2D coordinates to a text file
  SUBROUTINE xy_grid(points, box_size)

    INTEGER, PARAMETER :: dp=kind(1.0d0)
    !> @param box_size is the user defined parameter which defines the size of the domain in each axis.
    REAL(dp), INTENT(IN) :: box_size
    !> @param points is the user defined
    INTEGER, INTENT(IN) :: points
    ! Loop variables
    INTEGER :: i,j
    ! Resulting position values in x-y axis
    REAL(dp) :: a, b ! a is x, b is y
    INTEGER, PARAMETER :: file_no=3

    ! Opens an empty text file
    OPEN (unit=file_no,file="xyz.txt",action="write")
    DO i = -points, points
      DO j = -points, points
        ! Scales the positions to be between -box_size and +box_size
        a =  (real(i) / (abs(points) + abs(points)))*box_size
        b =  (real(j) / (abs(points) + abs(points)))*box_size
        ! Writes the coordinates to file, the z coordinate is set to zero as only using this for 2D contour plot which requires x and y coordinates to change.
        write (file_no,*) a,b,0
      END DO
    END DO
    ! Closes the file after writing
    CLOSE (file_no)
    PRINT*, 'Success in writing: xyz.txt'
  END SUBROUTINE xy_grid

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Creates a NetCDF file which contains information needed for the restart of the \ref Bi_Op_step function.

  !> This file will be created in the event that the \ref Bi_Op_step function crashes.


  ! The subroutine which created the NetCDF file containing the required information needed for
  SUBROUTINE write_restart_file(gp_n_data, gp_n_dof,n_cycles,no_samples, current_best_E,&
         constant_mean_value, gamma, kernel_var, kernal_inv_length, &
         E_data, data_mean, param_pres, param_cov, param_data)


    INTEGER, PARAMETER :: dp=kind(1.0d0)
    ! Declaring the input variables
    !> @param filename outputs a NetCDF file names restart.nc4
    CHARACTER(LEN=100) :: filename = 'restart.nc4'

    !> File ID's and loop variables
    INTEGER :: ierr, file_id, i

    !> Defining everything for the scalar variables
    !> @param gp_n_data is the data used for the GP surrogate
    !> @param gp_n_dof is the degrees of freedom used for the GP surrogate
    !> @param n_cycles is the number of cycles
    !> @param no_samples is the number of samples
    !> @param current_best_E is the current best energy found so far during the run.
    !> @param constant_mean_value is the current mean value used by the optimiser
    !> @param gamma is the current gamma value used by the optimiser
    !> @param kernel_var is the current kernal variance used by the optimiser
    !> @param kernal_inv_length is the current kernal lengthscale used by the optimiser
    REAL(dp) :: gp_n_data, gp_n_dof, n_cycles, no_samples, current_best_E
    REAL(dp) :: constant_mean_value, gamma, kernel_var, kernal_inv_length

    !> Variable IDs for the scalar variables
    INTEGER :: var_gp_n_data, var_gp_n_dof, var_n_cycles, var_no_samples
    INTEGER :: var_current_best_E, var_constant_mean_value, var_gamma
    INTEGER :: var_kernel_var, var_kernal_inv_length

    !> Defining everything for the array variables
    !> Rank 1 arrays
    !> @param E_data is a rank 1 array containing the energy data
    !> @param data_mean is a rank 1 array containing the current mean data used by the optimiser
    REAL(dp), DIMENSION(:), INTENT(IN) :: E_data
    REAL(dp), DIMENSION(:), INTENT(IN) :: data_mean

    !> Rank 2 arrays
    !> @param param_pres rank 2 array containing
    !> @param_cov rank 2 array containing the covariance matrix
    !> @param param_data rank 1 array containing the stored parameter data used by the optimiser

    REAL(dp), DIMENSION(:,:), INTENT(IN) :: param_pres
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: param_cov
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: param_data

    !> Dimension Names - named the same as the input variables
    CHARACTER(LEN=100), DIMENSION(1) :: dim_e_data = (/"E_Data"/)
    CHARACTER(LEN=100), DIMENSION(1) :: dim_data_mean = (/"data_mean"/)

    CHARACTER(LEN=100), DIMENSION(2) :: dim_param_pres = (/"param_pres","param_pres"/)
    CHARACTER(LEN=100), DIMENSION(2) :: dim_param_cov = (/"param_cov","param_cov"/)
    CHARACTER(LEN=100), DIMENSION(2) :: dim_param_data = (/"param_data","param_data"/)

    !> ID's
    ! Variable IDs
    INTEGER :: var_e_data, var_data_mean, var_param_pres
    INTEGER :: var_param_cov, var_param_data

    ! Dimension IDs for the arrays
    INTEGER, DIMENSION(1) :: did_e_data, did_data_mean
    INTEGER, DIMENSION(2) :: did_param_cov, did_param_data, did_param_pres

    ! Sizes for the arrays
    INTEGER, DIMENSION(1) :: size_e_data, size_data_mean
    INTEGER, DIMENSION(2) :: size_param_cov, size_param_data, size_param_pres


    !> Create the file
    ierr = nf90_create(filename, NF90_CLOBBER, file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    ! Variable sizes
    size_e_data = SHAPE(E_data)
    size_data_mean = SHAPE(data_mean)
    size_param_pres = SHAPE(param_pres)
    size_param_cov = SHAPE(param_cov)
    size_param_data = SHAPE(param_data)

    !> Defining the Scalar Variables

    ierr = nf90_def_var(file_id, "gp_n_data", NF90_DOUBLE, var_gp_n_data)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "gp_n_dof", NF90_DOUBLE, var_gp_n_dof)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "n_cycles", NF90_DOUBLE, var_n_cycles)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "no_samples", NF90_DOUBLE, var_no_samples)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "current_best_e", NF90_DOUBLE, var_current_best_E)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "constant_mean_value", NF90_DOUBLE, var_constant_mean_value)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "gamma", NF90_DOUBLE, var_gamma)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "kernel_var", NF90_DOUBLE, var_kernel_var)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_def_var(file_id, "kernal_inv_length", NF90_DOUBLE, var_kernal_inv_length)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    !> Inputting in the Array variables and defining the dimensions

    !> Rank 1 arrays first

    ierr = nf90_def_dim(file_id, dim_e_data(1), size_e_data(1), did_e_data(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "E_data", NF90_DOUBLE, did_e_data, var_e_data)

    ierr = nf90_def_dim(file_id, dim_data_mean(1), size_data_mean(1), did_data_mean(1))
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    ierr = nf90_def_var(file_id, "data_mean", NF90_DOUBLE, did_data_mean, var_data_mean)

    !> Rank 2 arrays

    DO i = 1,2
      ierr = nf90_def_var(file_id, dim_param_pres(i), size_param_pres(i), did_param_pres(i))
      IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    END DO
    ierr = nf90_def_var(file_id, "param_pres", NF90_DOUBLE, did_param_pres, var_param_pres)

    DO i = 1,2
      ierr = nf90_def_var(file_id, dim_param_cov(i), size_param_cov(i), did_param_cov(i))
      IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    END DO
    ierr = nf90_def_var(file_id, "param_cov", NF90_DOUBLE, did_param_cov, var_param_cov)

    DO i = 1,2
      ierr = nf90_def_var(file_id, dim_param_data(i), size_param_data(i), did_param_data(i))
      IF (ierr /= nf90_noerr) CALL handle_err(ierr)
    END DO
    ierr = nf90_def_var(file_id, "param_data", NF90_DOUBLE, did_param_data, var_param_data)


    ! Metadata
    ierr = nf90_enddef(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    !> Writing the array data to file

    !> Scalar Values

    ierr = nf90_put_var(file_id, var_gp_n_data, gp_n_data)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_gp_n_dof, gp_n_dof)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_n_cycles, n_cycles)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_no_samples, no_samples)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_current_best_E, current_best_E)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_constant_mean_value, constant_mean_value)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_gamma, gamma)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_kernel_var, kernel_var)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_kernal_inv_length, kernal_inv_length)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Rank 1 arrays

    ierr = nf90_put_var(file_id, var_e_data, E_data)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_data_mean, data_mean)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    !> Rank 2 arrays

    ierr = nf90_put_var(file_id, var_param_pres, param_pres)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_param_cov, param_cov)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    ierr = nf90_put_var(file_id, var_param_data, param_data)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)


    !> Closing the file
    ierr = nf90_close(file_id)
    IF (ierr /= nf90_noerr) CALL handle_err(ierr)

    PRINT *, "Success in writing file: ", filename

  END SUBROUTINE write_restart_file

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief Handles the NetCDF errors and print statements.
  !> This subroutine is only used and called within this module.

  ! This subroutine handles the errors for the writing netcdf file
  SUBROUTINE handle_err(ierr)

    INTEGER, INTENT(IN) :: ierr

    IF (ierr /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(ierr))
      ! Stops the code if an error is found
      STOP "Stopped"
    END IF

  END SUBROUTINE handle_err

END MODULE write_file

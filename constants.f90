!> @brief Definitions of shared constants used throughout the software
module shared_constants
  use iso_fortran_env

  implicit none
  integer, parameter :: dp = real64
  real(dp),parameter :: pi = 3.141592653589793238_dp
  real(dp),parameter :: electronvolt = 27.211386245988_dp

  ! Basis type codes:
  integer, parameter :: slater_1s_code = 100
  integer, parameter :: sto_3g_code = 200
end module

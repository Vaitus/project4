! Basic module used for common definitions of number kinds

module base

    use, intrinsic :: iso_fortran_env

    implicit none

    ! Default real type
    integer, parameter :: R8 = real64
    ! Long integer
    integer, parameter :: I8 = int64

end module base

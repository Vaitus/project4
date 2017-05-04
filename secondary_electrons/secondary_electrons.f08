! Generates positions of secondary electron emissions from the target
module secondary_electrons

    use base
    use rng

    private
    public :: setupTarget, getSE

    real(R8), parameter :: pi = 2d0 * asin(1d0)

    real(R8) :: x0, y0, z0  ! Position of the target center
    real(R8) :: tr  ! Target radius
    real(R8) :: eze  ! Exponent of erosion zone profile

    contains

    ! Sets target location and its diameter
    subroutine setupTarget(x, y, z, diameter, ezexp)
        implicit none
        real(R8), intent(in) :: x, y, z, diameter
        real(R8), intent(in), optional :: ezexp

        call initRN()

        x0 = x
        y0 = y
        z0 = z
        tr = diameter / 2d0

        if (present(ezexp)) then
            eze = ezexp
        else
            eze = 2.5d0
        end if
    end subroutine setupTarget

    ! Gets positions of a new secondary electron
    subroutine getSE(x, y, z)
        implicit none
        real(R8), intent(out) :: x, y, z
        real(R8) :: af, ar, rn

        ! Random position in phi coordinate
        af = 2.0 * pi * getRN1()

        ! Random position in r coordinate
        do
            ar = tr * getRN1()
            rn = getRN()
            if (rn < ezp(ar)) exit
        end do

        ! Calculate position in x, y, z
        x = ar * sin(af) + x0
        y = ar * cos(af) + y0
        z = z0
    end subroutine getSE

    ! Erosion zone profile function in form ezp = sin(pi r^e / tr^e)
    function ezp(r) result(res)
        implicit none
        real(R8), intent(in) :: r
        real(R8) :: res
        res = sin((pi * r**eze) / tr**eze)
    end function ezp

end module secondary_electrons

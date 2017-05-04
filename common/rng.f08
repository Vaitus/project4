! Module for random number generation
module rng

    use base
    use mt19937_64, only: init_genrand64, genrand64_real1, genrand64_real2

    implicit none

    contains

    ! Initialize random number generator according to current cpu time
    subroutine initRN()
        integer :: clk

        call system_clock(count=clk)
        call init_genrand64(int(clk, I8))
    end subroutine initRN

    ! Random number in the interval <0, 1)
    function getRN() result(r)
        real(R8) :: r
        r = genrand64_real1()
    end function getRN

    ! Random number in the interval <0, 1>
    function getRN1() result(r)
        real(R8) :: r
        r = genrand64_real2()
    end function getRN1

end module rng

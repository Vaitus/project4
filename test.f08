! Program for module testing
program test

    use base
    use cross_section
    use secondary_electrons

    implicit none

    type(cs_t) :: cs
    integer :: i
    real :: rn
    real(R8) :: E, r, x, y, z

! ****************************************************************************

    print *, "Test of cross_section module!"
    cs = loadCS("ics_ar.csv")

    print *, "E < Ei: E = 10eV, CS = ", getCS(cs, 1d1)
    open(100, file='outCS.txt', status='replace')
    do i = 1, 250
        call random_number(rn)
        E = cs%minE + (cs%maxE - cs%minE) * real(rn, R8)
        r = getCS(cs, E)
        write(100, *) E, ";", r
        print *, "E = ", E, "eV, CS = ", r
    end do
    close(100)

    call freeCS(cs)

! ****************************************************************************

    print *, "Test of SE generation!"
    call setupTarget(0d0, 0d0, 0d0, 100d0)

    open(100, file='outSE.txt', status='replace')
        do i = 1, 100000
            call getSE(x, y, z)
            write(100, *) x, ";", y, ";", z
            print *, "x = ", x, "; y = ", y, "; z = ", z
        end do
    close(100)

end program test

! Module for interpolation of cross-section
module cross_section

    use base

    implicit none

    ! Type for storing data about cross-section
    type cs_t
        real(R8), dimension(:,:), allocatable :: tbl  ! Table of cross-section values
        real(R8) :: minE, maxE  ! Minimum and maximum energy
        real(R8) :: a, b  ! Coefs. for determination of index
        integer :: N  ! Number of rows in tbl
    end type cs_t

    contains

        ! Loads cross-section data from file
        function loadCS(fname) result(cs)
            implicit none
            character(len=*), intent(in) :: fname
            type(cs_t) :: cs
            real(R8) :: x, y
            integer :: ios, num, i

            open(100, file=fname, status='old', iostat=ios)
            if (ios == 0) then
                read(100, *) num
                allocate(cs%tbl(num, 2))
                do i = 1, num
                    read(100, *) x, y
                    cs%tbl(i, 1) = x
                    cs%tbl(i, 2) = y
                end do
                close(100)
            else
                print *, "I can't read the cross-section from file ", fname, "!"
                return
            end if

            cs%minE = cs%tbl(1, 1)
            cs%maxE = cs%tbl(num, 1)
            cs%N = num
            cs%b = (real(num, 8) - 1) / (cs%maxE - cs%minE)
            cs%a = 1d0 - cs%minE * cs%b
        end function loadCS

        ! Gets interpolated value of cross-section for given energy
        function getCS(cs, E) result(r)
            implicit none
            type(cs_t), intent(inout) :: cs
            real(R8), intent(in) :: E
            real(R8) :: r
            real(R8) :: a, b
            integer(4) :: i

            if ((E < cs%minE) .or. (E > cs%maxE)) then
                r = 0.0
            else
                i = int(cs%a + cs%b * E)
                if (i == cs%N) i = i - 1  ! When E=maxE, i must be lowered
                b = (cs%tbl(i + 1, 2) - cs%tbl(i, 2)) / (cs%tbl(i + 1, 1) - cs%tbl(i, 1))
                a = cs%tbl(i, 2) - b * cs%tbl(i, 1)
                r = a + b * E
            end if
        end function getCS

        ! Destroys allocated memory for given cross-section
        subroutine freeCS(cs)
            implicit none
            type(cs_t), intent(inout) :: cs
            deallocate(cs%tbl)
        end subroutine freeCS

end module cross_section

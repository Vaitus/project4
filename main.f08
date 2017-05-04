
program main
    use leapfrog
    use mcc
    use stack
    use mcc

    !Type for a particle attributes
    type electron_t
            real(R8) :: x(3)
            real(R8) :: v(3)
            real(R8) :: timeInCell
            integer(I8) :: cellId
    end type electron_t

    real(R8), parameter :: m_e = 9.10938356e-31, q_e = -1.6021766208e-19
    real(R8), parameter :: m_ar = 39.944 * 1.660539040e-27 - m_e, q_ar = -q_e

    integer, parameter :: N = 400 !number of cells in a mesh
    integer, parameter :: N_shift = N/2 !shift for storing particles in array
    integer, parameter :: STEPS = 10000!00 !number of steps of simulation
    real(R8) :: dt = 0.000001 !delta time
    integer(I8), parameter :: numberOfElectrons = 5000

    type(electron_t) :: electron(numberOfElectrons) !array of electrons
    type(stack_var) :: stack_e !stack for free spots after used electrons
    real(R8) :: cellTime(N**3) !time of electrons in a cell given by electrons' position
    integer(I8) :: lastArraySize = 0 !Size of array of created electrons

    type(cs_t) :: cs !type for saving the cross-section
    cs = loadCS("ics_ar.csv") !loading interpolation of cross-section

    call setupTarget(0d0, 0d0, 0d0, 100d0) !setuping a target
    call simulate(q_e, m_e, dt, 'elektron.txt', STEPS) !calling a simulation

        !debugging stuff
        do i = 1, size(cellTime)
            if (cellTime(i) /= 0) then
                print *, cellTime(i), i
            end if
        end do

    contains

        !Subroutine for simulating of position of a particle
        subroutine simulate(q, m, dt, name, steps)
            implicit none

            real(R8), intent(in) :: dt !cast jak velke kroky
            real(R8) :: B(3), E(3) !pole pro mag. indukci a E. pole
            real(R8), intent(in) :: q, m  !naboj a hmotnost
            integer, intent(in) :: steps

            integer :: indexForWrtToFile
            real(R8), dimension(:), allocatable  :: positionsOfSecondaryElectrons

            character(len=*), intent(in) :: name
            integer :: i, j, k
            real(R8) :: qm
            qm = q/m

            !electron(1)%x = startingPositionForX()
            !electron(1)%v = startingVelocity()
            !B = getB(electron(1)%x)
            !E = getE(electron(1)%x)
            !call updateVelocity(electron(1)%x, electron(1)%v, E, B, qm, -0.5*dt) !posunuti castice o 0.5 "dozadu"

            indexForWrtToFile = 1

            call creationOfElectrons(startingPositionForX(), qm, E, B)

            do i = 1, steps
                do j = 1, lastArraySize
                    if(electron(j)%cellId == -1) cycle !checking if electron exists

                    !collisions
                    !call collisionMCC(getCS(cs, energy(electron(j)%v, m_e)), electron(j)%v, dt)

                    if (.not. allocated(positionsOfSecondaryElectrons)) allocate(positionsOfSecondaryElectrons(5000))
                    E = getE(electron(j)%x)
                    B = getB(electron(j)%x)

                    call updateVelocity(electron(j)%x, electron(j)%v, E, B, qm, dt)
                    call particlePush(electron(j)%x, electron(j)%v, dt)

                    if (modulo(i,10) == 0 .and. .not. isOutOfCell(electron(j)%x)) then
                        positionsOfSecondaryElectrons(indexForWrtToFile) = electron(j)%x(1)
                        positionsOfSecondaryElectrons(indexForWrtToFile+1) = electron(j)%x(2)
                        positionsOfSecondaryElectrons(indexForWrtToFile+2) = electron(j)%x(3)
                        indexForWrtToFile = indexForWrtToFile + 3
                    end if

                    print *, electron(j)%x, j

                    if ((size(positionsOfSecondaryElectrons)-indexForWrtToFile) <= 3) then
                        call writePositionsToFile(positionsOfSecondaryElectrons, indexForWrtToFile)
                        indexForWrtToFile = 1
                        deallocate(positionsOfSecondaryElectrons)
                    end if

                    if (isOutOfCell(electron(j)%x)) then
                        electron(j)%cellId = -1
                        call push(stack_e, j)
                        call writePositionsToFile(positionsOfSecondaryElectrons,indexForWrtToFile)
                        indexForWrtToFile = 1
                        deallocate(positionsOfSecondaryElectrons)

                        !if out making new electrons
                        !call creationOfElectrons(startingPositionForX(), qm, E, B)
                    else
                        if (electron(j)%cellId /= getCellId(electron(j)%x)) then
                            cellTime(electron(j)%cellId) = cellTime(electron(j)%cellId) + electron(j)%timeInCell
                            electron(j)%timeInCell = 0.0
                            electron(j)%cellId = getCellId(electron(j)%x)
                        else
                            electron(j)%timeInCell = electron(j)%timeInCell + dt
                        end if
                    end if
                end do
            end do

            do k = 1, lastArraySize
                !po simulaci pridani zbytku
                if (electron(k)%cellId /= -1) then
                    cellTime(electron(k)%cellId) = cellTime(electron(k)%cellId) + electron(k)%timeInCell
                    if (allocated(positionsOfSecondaryElectrons)) then
                        call writePositionsToFile(positionsOfSecondaryElectrons, indexForWrtToFile)
                        deallocate(positionsOfSecondaryElectrons)
                    end if
                end if
            end do

        end subroutine simulate

        subroutine creationOfElectrons(x, qm, E, B)
            implicit none

            real(R8), intent(in) :: x(3)
            real(R8), intent(in) :: qm
            real(R8), intent(out) :: B(3), E(3)
            integer(I8) :: indexForElectron

            if(stack_e%size > numberOfElectrons/2) then
                indexForElectron = pop(stack_e)
            else
                lastArraySize = lastArraySize + 1
                indexForElectron = lastArraySize
            end if

            electron(indexForElectron)%x = x
            electron(indexForElectron)%v = startingVelocity()
            electron(indexForElectron)%cellId = getCellId(x)
            electron(indexForElectron)%timeInCell = 0.0
            E = getE(electron(indexForElectron)%x)
            B = getB(electron(indexForElectron)%x)
            call updateVelocity(electron(indexForElectron)%x, electron(indexForElectron)%v, E, B, qm, -0.5*dt) !posunuti castice o 0.5 "dozadu"
        end subroutine

        !Subroutine for writing array of stored positions of a particle to a file
        subroutine writePositionsToFile(array, indexL)
            implicit none

            real(R8), dimension(:), intent(in)  :: array
            integer, intent(in) :: indexL
            integer(I8) :: i

            logical :: exist

            inquire(file="test.txt", exist=exist)
            if (exist) then
                open(12, file="test.txt", status="old", position="append", action="write")
            else
                open(12, file="test.txt", status="new", action="write")
            end if
                do i = 1, indexL, 3
                    write(12, *) array(i),";", array(i+1),";", array(i+2)
                end do
            close(12)

        end subroutine

        !Function for getting cellId of a particle for its position
        function getCellId(x) result(out)
            real(R8), intent(in) :: x(3)
            real(R8) :: out
            integer(R8) :: i,j,k

            !k = floor(((x(3)))*2)
            !j = floor(((x(2))+100)*2)
            !i = floor(((x(1))+100)*2)
            k = x(3)/0.5
            j = x(2)/0.5 + N_shift
            i = x(1)/0.5 + N_shift
            !k*N**2 + j*N + i + 1
            out = ((k*N) + j)*N + i + 1
        end function getCellId

        !Function for checking if particle is out of cell
        function isOutOfCell(x) result(out)
            real(R8), intent(in) :: x(3)
            logical :: out

            out = x(1) <= real(-100, R8) .OR. x(1) >= real(100, R8) .OR. &
                  x(2) <= real(-100, R8) .OR. x(2) >= real(100, R8) .OR. &
                  x(3) < real(0, R8) .OR. x(3) >= real(200, R8)
        end function isOutOfCell



    !real(kind=p16) :: dt
    !character :: s
    !real :: Bmag, tc
    !real, parameter :: PI = 3.141592653589793238462643383279
    !real, parameter :: TMUL(9) = (/ 0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0 /)
    !Bmag = sqrt(dot_product(Bindukce((/0.0, 0.0, 0.0/)+0.0_p16), Bindukce((/0.0, 0.0, 0.0/)+0.0_p16)))

    ! Cas gyrace pro elektron
    !tc = 2 * PI * m_e / (q_e * Bmag)

   ! do i = 1,9
   !     dt = TMUL(i) * tc
   !     write (s, '(I1)') i
   !     call simulate(q_e, m_e, dt, 'elektron'//s//'.txt', STEPS)
   ! end do

    ! Cas gyrace pro argonovy iont
    !tc = 2 * PI * m_ar / (q_e * Bmag)

    !do i = 1,9
    !    dt = TMUL(i) * tc
    !    write (s, '(I1)') i
    !    call simulate(q_ar, m_ar, dt, 'argon'//s//'.txt', STEPS)
    !end do

end program main

module mcc
    use leapfrog
    use mt_64

    type node
            real(kind=p16) :: position(3)
            real(kind=p16) :: Efield(3),Bfield(3)
    end type node

    type particle
            real(kind=p16) :: x(3)
            real(kind=p16) :: v(3)
            real(kind=p16) :: E(3),B(3)
    end type particle

    contains

    subroutine collisionMCC(dt)
        implicit none
        real(kind=p16), intent(in) :: dt !cast jak velke skoky
        real(kind=p16) :: nn !neutral density
        real(kind=p16) :: g !velocity
        real(kind=p16) :: a !cross section
        real(kind=p16) :: sigma ! cross section
        real(kind=p16) :: probability !

        nn = 1
        g = 1
        a = 1 !temporal
        sigma = a*a !temp
        probability = 1 - exp(-nn*sigma*g*dt)

        if (genrand64_real1()<probability) then
            call performCollision()
        end if
    end subroutine collisionMCC

    subroutine performCollision()
        implicit none
    end subroutine performCollision

end module mcc

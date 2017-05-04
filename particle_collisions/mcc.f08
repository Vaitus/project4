! Module for Monte Carlo Collisions
module mcc
    use base
    use rng
    use cross_section

    implicit none

    contains

        subroutine collisionMCC(sigma, v, dt)
            implicit none

            real(R8), intent(in) :: dt !delta time
            real(R8) :: nn !neutral density
            real(R8), intent(in) :: v(3) !velocity
            real(R8) :: g !speed
            real(R8), intent(in) :: sigma !cross section
            real(R8) :: probability !

            nn = 1
            g = sqrt(dot_product(v,v))
            !sigma = a !temp
            probability = 1 - exp(-nn*sigma*g*dt)

            if (getRN1()<probability) then
                call performCollision()
            end if
        end subroutine collisionMCC

        subroutine performCollision()
            implicit none

        end subroutine performCollision

end module mcc

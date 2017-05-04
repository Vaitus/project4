module leapfrog
    use base
    use cross_section
    use secondary_electrons

    implicit none

    !real, parameter :: permitivityOfFreeSpace = 8.854187817e-12
    !real, parameter :: stability_e = 0.1/(sqrt((1*q_e**2)/(m_e*permitivityOfFreeSpace)))
    !real, parameter :: stability_ar = 0.1/(sqrt((17*q_e**2)/(m_e*permitivityOfFreeSpace)))

    contains

        !vypocet energie
        function energy(v, m) result(energ)
            implicit none
            real(R8), intent(in) :: v(3), m
            real(R8) :: energ
            real, parameter :: c_light = 299792458.0

            energ = 0.5 * m * dot_product(v,v)
            !energ = sqrt( dot_product((v),(v))*m*m*c_light**2 + (m*c_light**2)**2 )
            !energ = ((1/(sqrt(1-(dot_product((v),(v))/c_light**2))))-1)*m*c_light*c_light

        end function energy

        !nastaveni startovni pozice
        function startingPositionForX() result(out)
            implicit none
            real(R8) :: out(3)
            real(R8) :: x1,y1,z1
            call getSE(x1, y1, z1)

            out = (/x1,y1,z1/)
        end function startingPositionForX

        !nastaveni startovni rychlosti
        function startingVelocity() result(out)
            implicit none
            real(R8) :: out(3)

            out = (/0.0,0.0,0.0002/)
        end function startingVelocity

        !vektorovej soucin
        function cross_product(array1, array2) result(out)
            implicit none
            real(R8), intent(in) :: array1(3), array2(3)
            real(R8) :: out(3)

            out(1) = array1(2) * array2(3) - array1(3)*array2(2)
            out(2) = -array1(1) * array2(3) + array1(3)*array2(1)
            out(3) = array1(1) * array2(2) - array1(2)*array2(1)
        end function cross_product

        !posunuti castice o rychlost
        subroutine particlePush(x, v, dt)
            real(R8), intent(inout) :: x(3)
            real(R8), intent(in) :: v(3)
            real(R8), intent(in) :: dt

            x = x + v * dt
        end subroutine particlePush

        !upraveni rychlosti
        subroutine updateVelocity(x, v, E, B, qm, dt)
            real(R8), intent(in) :: x(3), E(3), B(3)
            real(R8), intent(inout) :: v(3)
            real(R8) :: v_minus(3), v_prime(3), v_plus(3)
            real(R8) :: t(3), s(3), t_norma
            real(R8) :: v_minus_crossproduct_t(3), v_prime_cross_s(3)
            real(R8), intent(in) :: dt
            real(R8), intent(in)  :: qm

            !vektor t uprav
            t = qm * B * 0.5 * dt

            !norma vektoru t
            t_norma = t(1) * t(1) + t(2) * t(2) + t(3) * t(3)

            ! 2t/1+t^2
            s = 2 * (t / (1.0 + t_norma))

            !v_minus
            v_minus = v + qm * E * 0.5 * dt

            !v_prime
            v_minus_crossproduct_t = cross_product(v_minus, t)
            v_prime = v_minus + v_minus_crossproduct_t

            !v_plus
            v_prime_cross_s = cross_product(v_prime, s)
            v_plus = v_minus + v_prime_cross_s

            !v^n+1/2
            v = v_plus + qm * E * 0.5 * dt
        end subroutine updateVelocity

        !nastaveni magneticke indukce
        function getB(x) result(out)
            implicit none
            real(R8), intent(in) :: x(3)
            real(R8) :: out(3)

            out = (/0.0,0.0,0.0/)
        end function getB

        !elektro. pole
        function getE(x) result(out)
            implicit none
            real(R8), intent(in) :: x(3)
            real(R8) :: out(3)

            out = (/0.0,0.0,-0.0002/)
        end function getE

end module leapfrog


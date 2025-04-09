module math_functions
    use iso_fortran_env, only: int32, real64
contains

    ! given the the three catersian points of a triangule
    ! this function returns all interior angles and the side lengths 
    subroutine cosine_law(A, B, C, response)
        implicit none
        real (real64), intent(in) :: A(3), B(3), C(3) !triangle sides
        real (real64), intent(inout) :: response(6)
        real (real64) :: AB, BC, AC, cos_a, cos_b, cos_c

        !Respose will return the triangule measurements AB,AC,BC, and the angles CAB, ABC, BCA

        AB = norm2(B - A)
        AC = norm2(C - A)
        BC = norm2(B - C)

        response(1) = AB
        response(2) = AC
        response(3) = BC

        cos_a = (AB**2 + AC**2 - BC**2) / (2d0 * AC * AB)
        cos_b = (AB**2 + BC**2 - AC**2) / (2d0 * BC * AB)
        cos_c = (AC**2 + BC**2 - AB**2) / (2d0 * BC * AC)

        response(4) = cos_a
        response(5) = cos_b
        response(6) = cos_c

        return
    end subroutine cosine_law

    subroutine molecular_rotation_zyz(arr, array_size, alpha, beta, gamma, rotarr)
        implicit none
        integer (int32), intent(in) :: array_size

        real (real64), intent(in) :: arr(3, array_size)
        real (real64), intent(inout) :: rotarr(3, array_size)
        real (real64), intent(in) :: alpha, gamma, beta

        real (real64) :: a, b, g, ca, cb, cg, sa, sb, sg, U_rot(3, 3)

        a = alpha
        g = gamma
        b = beta

        ca = dcos(a)
        cb = dcos(b)
        cg = dcos(g)

        sa = dsin(a)
        sb = dsin(b)
        sg = dsin(g)


        ! construct rotation matrix

        U_rot = 0d0
        U_rot(1, 1) = ca * cb * cg - sa * sg
        U_rot(2, 2) = ca * cg - cb * sa * sg
        U_rot(3, 3) = cb
        U_rot(1, 2) = -cg * sa - ca * cb * sg
        U_rot(1, 3) = ca * sb
        U_rot(2, 1) = ca * sg + cb * cg * sa
        U_rot(2, 3) = sa * sb
        U_rot(3, 1) = -cg * sb
        U_rot(3, 2) = sb * sg

        call rotate_molecule_v2 (array_size, arr, rotarr, U_rot)

    end subroutine molecular_rotation_zyz

    subroutine interatomic_distance(arr, array_size, distance)
        implicit none
        integer (int32), intent(in) :: array_size
        integer (int32) :: count, i, j
        real (real64), intent(in) :: arr(3, array_size)
        real (real64), intent(inout) :: distance(array_size * (array_size - 1) / 2)
        real (real64) :: xij, yij, zij

        count = 0;




        do i = 1, array_size - 1
            do j = i + 1, array_size

                count = count + 1;
                xij = (arr(1, i) - arr(1, j))**2
                yij = (arr(2, i) - arr(2, j))**2
                zij = (arr(3, i) - arr(3, j))**2
                distance(count) = DSQRT(xij + yij + zij)

            enddo
        enddo

    end subroutine interatomic_distance

    integer (int32) function cmp_function(a, b)
        integer (int32) a, b
        if (a < b) compar = -1
        if (a == b) compar = 0
        if (a > b) compar = 1
        return
    end

    subroutine compared_distance(distance1, distance2, array_size, result)
        implicit none

        integer (int32), intent(in) :: array_size
        real (real64), intent(in) :: distance1(array_size), distance2(array_size)
        real (real64) :: D1(array_size), D2(array_size)
        real (real64), intent(out) :: result(2)
        real (real64) :: sr1, sr2
        integer (int32) :: i

        result = 0d0

        D1 = distance1
        D2 = distance2
        call sort_pick(D1, array_size)
        call sort_pick(D2, array_size)

        do i = 1, array_size
            result(1) = result(1) + (D1(i) - D2(i))**2
        enddo
        result(1) = DSQRT(result(1))

        !Global difference

        sr1 = 0d0
        sr2 = 0d0
        do i = 1, array_size

            sr1 = sr1 + D1(i)
            sr2 = sr2 + D2(i)
        enddo

        result(2) = DABS(sr1 - sr2)

    end subroutine compared_distance

    subroutine sort_pick(arr, array_size)

        implicit none
        integer (int32), intent(in) :: array_size
        !Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced on output by its sorted rearrangement. 
        real (real64), dimension(:), intent(inout) :: arr(array_size)
        integer (int32) :: i, j
        real (real64) :: a

        do j = 2, array_size
            a = arr(j)
            do i = j - 1, 1, -1 !Pick out each element in turn. Look for the place to insert it.
                if (arr(i) <= a) exit
                arr(i + 1) = arr(i)
            end do
            arr(i + 1) = a !Insert it.
        end do
    end subroutine sort_pick

    subroutine print_vector(vec, array_size, str)
        implicit none
        integer (int32), intent(in) :: array_size
        real (real64), intent(in) :: vec(3, array_size)
        character(len = 20, kind = 1), intent(in) :: str
        integer (int32) :: i

        write(*, *) "******* Vector Printing: ", " ************"
        write(*, *) "******* ", str, " ************"
        do i = 1, array_size
            write(*, *)vec(:, i)
        enddo

        write(*, *) "******* end Printing ************"

    end subroutine print_vector

    subroutine find_euler_angles(natom, cart_ref1, cart_mat1, masses, alpha, beta, gamma)

        implicit none

        integer (int32), intent(in) :: natom
        real (real64), intent(out) :: alpha, beta, gamma
        real (real64), intent(in) :: masses(natom), cart_ref1(3, natom), cart_mat1(3, natom)
        integer (int32) :: ierr, i
        real (real64) :: U_rot(3, 3), quat(4), threshold, U33
        !!! ----------------------------------------------------------------------------------
        !!! 4) once the geometry is in the proper frame, find the corresponding Euler angles:
        !!! ----------------------------------------------------------------------------------

        call qtrfit(natom, cart_ref1, cart_mat1, masses, quat, U_rot, ierr)

        U33 = U_rot(3, 3)

        if(U33>1d0)then
            U33 = 1d0
        elseif(U33<-1d0)then
            U33 = -1d0
        endif

        ! solve for Euler angles
        threshold = 1d0 - 1d-12!0.9993d0

        beta = dacos(U33)

        if(U33>threshold)then
            gamma = datan2(-U_rot(1, 2), U_rot(1, 1))
            alpha = 0d0
        elseif(U33 <-1d0 * threshold)then
            gamma = datan2(U_rot(1, 2), -U_rot(1, 1))
            alpha = 0d0
        else
            alpha = datan2(U_rot(2, 3), U_rot(1, 3))
            gamma = datan2(U_rot(3, 2), -U_rot(3, 1))
        endif

        return
    end subroutine find_euler_angles

    subroutine center_of_mass(cart, cm, mass, natom, natom2)
        implicit none
        integer (int32) :: k, kp, natom, natom2
        real (real64) :: mass(natom), cart(natom * 3), mtot, cm(3)
        mtot = 0d0
        do k = natom - natom2 + 1, natom
            mtot = mtot + mass(k)
        enddo
        cm = 0d0
        do k = natom - natom2 + 1, natom
            do kp = 1, 3
                cm(kp) = cm(kp) + cart((k - 1) * 3 + kp) * mass(k)
            enddo
        enddo
        cm = cm / mtot
        return
    end subroutine center_of_mass

    subroutine center_of_mass_v2(cart, cm, mass, natoms)
        implicit none
        integer (int32) :: k, kp, natoms
        real (real64) :: mass(natoms), cart(natoms * 3), mtot, cm(3)
        mtot = 0d0
        do k = 1, natoms
            mtot = mtot + mass(k)
        enddo
        cm = 0d0
        do k = 1, natoms
            do kp = 1, 3
                cm(kp) = cm(kp) + cart((k - 1) * 3 + kp) * mass(k)
            enddo
        enddo
        cm = cm / mtot
        return
    end subroutine center_of_mass_v2

    subroutine remove_center_of_mass(cart, mass, natom, natom1)
        implicit none
        integer (int32) :: k, kp, natom, natom1
        real (real64) :: mass(natom), cart(natom * 3), mtot, cmass1(3)
        mtot = 0d0
        do k = 1, natom1
            mtot = mtot + mass(k)
        enddo
        cmass1 = 0d0
        do k = 1, natom1
            do kp = 1, 3
                cmass1(kp) = cmass1(kp) + cart((k - 1) * 3 + kp) * mass(k)
            enddo
        enddo
        cmass1 = cmass1 / mtot

        do k = 1, natom
            do kp = 1, 3
                cart((k - 1) * 3 + kp) = cart((k - 1) * 3 + kp) - cmass1(kp)
            enddo
        enddo
        return
    end subroutine remove_center_of_mass

    subroutine vec_to_mat2(cart_perms, cart_mat, natom)
        implicit none
        integer (int32) :: k, kp, natom
        real (real64) :: cart_perms(3 * natom), cart_mat(3, natom)
        do k = 1, natom
            do kp = 1, 3
                cart_mat(kp, k) = cart_perms((k - 1) * 3 + kp)
            enddo
        enddo
        return
    end subroutine vec_to_mat2

    subroutine mat_to_vec2(cart_mat, cart_perms, natom)
        implicit none
        integer (int32) :: k, kp, natom
        real (real64) :: cart_perms(3 * natom), cart_mat(3, natom)
        do k = 1, natom
            do kp = 1, 3
                cart_perms((k - 1) * 3 + kp) = cart_mat(kp, k)
            enddo
        enddo
        return
    end subroutine mat_to_vec2


    subroutine rotate_molecule_v2 (array_size, x, molrot, u)
        implicit none

        integer (int32), intent(in) :: array_size
        real (real64), intent(in) :: x(3, array_size)
        real (real64), intent(out) :: molrot(3, array_size)
        real (real64), intent(in) :: u(3, 3)
        real (real64) :: mol(3, array_size), urot(3, 3)
        integer (int32) :: i

        urot = u
        mol = x

        DO i = 1, array_size
            molrot(1, i) = urot(1, 1) * mol(1, i) + urot(1, 2) * mol(2, i) + urot(1, 3) * mol(3, i)
            molrot(2, i) = urot(2, 1) * mol(1, i) + urot(2, 2) * mol(2, i) + urot(2, 3) * mol(3, i)
            molrot(3, i) = urot(3, 1) * mol(1, i) + urot(3, 2) * mol(2, i) + urot(3, 3) * mol(3, i)

        enddo

        return
    end subroutine  rotate_molecule_v2


end module math_functions
    
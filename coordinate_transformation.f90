module coordinateTransf
    use iso_fortran_env, only: int32, real64
contains


    ! convert: Internal coordinates in zyz--> BiSpherical for 5D ONLY
    subroutine zyz_to_bispherical(internal, internal0_zyz, internal0_bisph)

        !************** subroutine definition
        ! subroutine zyz --> BiSpherical (convert: Internal coordinates in zyz--> BiSpherical for 5D)
        ! internal0_zyz(1) = R
        ! internal0_zyz(2) = cos_b1
        ! internal0_zyz(3) = cos_b2
        ! internal0_zyz(4) = alpha
        ! internal0_zyz(5) = gamma

        ! internal0_bisph(1) = R
        ! internal0_bisph(2) = theta1
        ! internal0_bisph(3) = theta2
        ! internal0_bisph(4) = phi1
        ! internal0_bisph(5) = phi2
        !********************************************************
        implicit none
        real(real64), intent(in) :: internal0_zyz(5), internal(5)
        real(real64), intent(out) :: internal0_bisph(5)

        real(real64) :: term1, term2, threshold, internal0p_zyz(5)
        real(real64), parameter :: PII = dacos(-1d0)

        threshold = 1d0 - 1d-12!0.9993d0
        internal0p_zyz = internal0_zyz

        internal0_bisph(1) = internal0_zyz(1) ! R

        if(internal0p_zyz(2)>1d0)then
            internal0p_zyz(2) = 1d0
        elseif(internal0p_zyz(2)<-1d0)then
            internal0p_zyz(2) = -1d0
        endif

        if(internal0p_zyz(3)>1d0)then
            internal0p_zyz(3) = 1d0
        elseif(internal0p_zyz(3)<-1d0)then
            internal0p_zyz(3) = -1d0
        endif


        ! Note : There is a discontinuity in Phi =0 for cos_b1 =0
        ! Taylor serie for cos b1
        ! Analyse b1 = 0 and b2 = 0

        if(dabs(internal0p_zyz(2))>threshold .and. dabs(internal0p_zyz(3))>threshold)then

            if(internal0p_zyz(2)>threshold)then
                internal0_bisph(2) = 0d0
                internal0_bisph(4) = internal(4)
            elseif(internal0p_zyz(2) <-1d0 * threshold)then
                internal0_bisph(2) = PII
                internal0_bisph(4) = internal(4)
            else
                internal0_bisph(2) = dacos(internal0p_zyz(2))
                internal0_bisph(4) = PII - internal0p_zyz(5)
            end if
            if(internal0p_zyz(3)>threshold)then
                internal0_bisph(3) = 0d0
                internal0_bisph(5) = internal(5)
            elseif(internal0p_zyz(3) <-1d0 * threshold)then
                internal0_bisph(3) = PII
                internal0_bisph(5) = internal(5)
            end if

        else

            if(internal0p_zyz(2)>threshold)then
                internal0_bisph(2) = 0d0
                internal0_bisph(3) = dacos(internal0p_zyz(3))
                internal0_bisph(4) = internal(4)
                internal0_bisph(5) = -1d0 * (internal0p_zyz(4) + internal0p_zyz(5))
            elseif(internal0p_zyz(2) <-1d0 * threshold)then
                internal0_bisph(2) = PII
                internal0_bisph(3) = dacos(-1d0 * internal0p_zyz(3))
                internal0_bisph(4) = internal(4)
                internal0_bisph(5) = internal0p_zyz(4) + internal0p_zyz(5) - PII
            elseif(internal0p_zyz(3) >threshold)then
                internal0_bisph(2) = dacos(internal0p_zyz(2))
                internal0_bisph(3) = 0d0
                internal0_bisph(4) = PII - internal0p_zyz(5)
                internal0_bisph(5) = internal(5)
            elseif(internal0p_zyz(3) <-1d0 * threshold)then
                internal0_bisph(2) = dacos(internal0p_zyz(2))
                internal0_bisph(3) = PII
                internal0_bisph(4) = PII - internal0p_zyz(5)
                internal0_bisph(5) = internal(5)
            else

                internal0_bisph(2) = dacos(internal0p_zyz(2))
                internal0_bisph(3) = dacos(internal0p_zyz(2) * internal0p_zyz(3) + dcos(internal0p_zyz(4)) &
                        * dsqrt((1d0 - internal0p_zyz(2)**2) * (1d0 - internal0p_zyz(3)**2)))

                internal0_bisph(4) = PII - internal0p_zyz(5)

                term1 = dcos(internal0p_zyz(5)) * dsin(internal0p_zyz(4)) * dsqrt(1d0 - internal0p_zyz(3)**2) &
                        + dsin(internal0p_zyz(5)) * (dcos(internal0p_zyz(4)) * internal0p_zyz(2) * &
                        dsqrt(1d0 - internal0p_zyz(3)**2) - dsqrt(1d0 - internal0p_zyz(2)**2) * internal0p_zyz(3))

                term2 = dsqrt(1d0 - internal0p_zyz(3)**2) * (dcos(internal0p_zyz(4)) * &
                        internal0p_zyz(2) * dcos(internal0p_zyz(5)) - dsin(internal0p_zyz(4)) * dsin(internal0p_zyz(5))) &
                        - internal0p_zyz(3) * dcos(internal0p_zyz(5)) * dsqrt(1d0 - internal0p_zyz(2)**2)

                internal0_bisph(5) = datan2(-1d0 * term1, term2)
            end if
        end if

        !Output Alpha Check

        if (internal0_bisph(4) > PII)then
            internal0_bisph(4) = internal0_bisph(4) - 2d0 * PII
        end if
        if (internal0_bisph(4) < -PII)then
            internal0_bisph(4) = internal0_bisph(4) + 2d0 * PII
        end if

    end subroutine zyz_to_bispherical

    ! transformation from Autosurf internal angles coordinates to zyz angles
    subroutine internal0_to_angles_zyz(internal0, xdim, angles)

        implicit none
        integer, intent(in) :: xdim
        real(real64), intent(in) :: internal0(xdim)
        real(real64), intent(inout) :: angles(6)

        real(real64) :: PII

        PII = dacos(-1d0)

        if (xdim==2)then
            angles(1) = 0d0
            angles(2) = dacos(internal0(2))
            angles(3) = 0d0
            angles(4) = 0d0
            angles(5) = 0d0
            angles(6) = 0d0
        elseif (xdim==3)then
            angles(1) = 0d0
            angles(2) = dacos(internal0(2))
            angles(3) = PII - internal0(3)
            angles(4) = 0d0
            angles(5) = 0d0
            angles(6) = 0d0
        elseif (xdim==4)then
            angles(1) = internal0(4)
            angles(2) = dacos(internal0(2))
            angles(3) = 0d0
            angles(4) = 0d0
            angles(5) = dacos(internal0(3))
            angles(6) = 0d0
        elseif (xdim==5)then
            angles(1) = internal0(4)
            angles(2) = dacos(internal0(2))
            angles(3) = internal0(5)
            angles(4) = 0d0
            angles(5) = dacos(internal0(3))
            angles(6) = 0d0
        elseif (xdim==6)then
            angles(1) = internal0(4)
            angles(2) = dacos(internal0(2))
            angles(3) = internal0(5)
            angles(4) = 0d0
            angles(5) = dacos(internal0(3))
            angles(6) = internal0(6)
        endif

    end subroutine internal0_to_angles_zyz

    subroutine euler_angles_to_bispherical(inter, R_zyz, internal0)
        ! *** internal coordinates ***
        ! ----------------------------

        !************** subroutine definition
        ! subroutine zyz --> BiSpherical (convert: Internal coordinates in zyz--> BiSpherical for 5D)

        ! Input : R and Euler Angles in radians

        ! Output : R and cos_theta1,cos_theta2, phi_1,phi_2

        !********************************************************
        implicit none
        real(real64), intent(in) :: R_zyz(7), inter(5)
        real(real64), intent(out) :: internal0(5)

        real(real64) :: alpha1, beta1, gamma1, alpha2, beta2, gamma2
        real(real64) :: term1, term2, threshold, alpha
        real(real64), parameter :: PII = dacos(-1d0)

        internal0(1) = R_zyz(1)
        alpha1 = R_zyz(2)
        beta1 = R_zyz(3)
        gamma1 = R_zyz(4)
        alpha2 = R_zyz(5)
        beta2 = R_zyz(6)
        gamma2 = R_zyz(7)

        alpha = alpha1 - alpha2
        if(alpha>PII)then
            alpha = alpha - 2d0 * PII
        endif
        if(alpha<-PII)then
            alpha = alpha + 2d0 * PII
        endif

        threshold = 1d-8
        internal0(2) = beta1

        if (dabs(beta1)<threshold  .or. dabs(beta1 - PII)< threshold)then

            internal0(4) = inter(4)

            if (dabs(beta1)<threshold)then
                internal0(2) = 0d0
                internal0(3) = beta2
            else
                internal0(2) = PII
                internal0(3) = PII - beta2
            end if

            if (dabs(beta2)<threshold .or. dabs(beta2 - PII)<threshold)then
                internal0(5) = inter(5)
            else
                if (dabs(beta1)<threshold)then
                    internal0(5) = -1d0 * (alpha + gamma1)
                else
                    internal0(5) = alpha + gamma1 - PII
                end if
            end if

        else
            if (dabs(beta1)<2d-6  .or. dabs(beta1 - PII)< 2d-6)then
                internal0(4) = inter(4)
            else
                internal0(4) = PII - gamma1
            end if

            if (dabs(beta2)<threshold .or. dabs(beta2 - PII)<threshold)then
                internal0(5) = inter(5)
                if (dabs(beta2)<threshold)then
                    internal0(3) = 0d0
                else
                    internal0(3) = PII
                end if
            else
                internal0(3) = dacos(dcos(beta1) * dcos(beta2) + dcos(alpha) * dsin(beta1) * dsin(beta2))

                term1 = dcos(gamma1) * dsin(alpha) * dsin(beta2) + dsin(gamma1) * (dcos(alpha) * dcos(beta1) * dsin(beta2) &
                        - dsin(beta1) * dcos(beta2))

                term2 = dsin(beta2) * (dcos(alpha) * dcos(beta1) * dcos(gamma1) - dsin(alpha) * dsin(gamma1)) &
                        - dcos(beta2) * dcos(gamma1) * dsin(beta1)

                internal0(5) = datan2(-1d0 * term1, term2)
            end if

        end if




        !Output Alpha Check

        if (internal0(4) > PII)then
            internal0(4) = internal0(4) - 2d0 * PII
        end if
        if (internal0(4) < -PII)then
            internal0(4) = internal0(4) + 2d0 * PII
        end if

        return
    end subroutine euler_angles_to_bispherical

    subroutine euler_angles_to_autosurf(xdim, R_zyz, internal0)
        implicit none
        Integer, intent(in) :: xdim
        real(real64), intent(in) :: R_zyz(7)
        real(real64), intent(out) :: internal0(xdim)
        real(real64) :: PII, alpha1, beta1, gamma1, alpha2, beta2, gamma2

        PII = dacos(-1d0)
        internal0(1) = R_zyz(1)
        alpha1 = R_zyz(2)
        beta1 = R_zyz(3)
        gamma1 = R_zyz(4)
        alpha2 = R_zyz(5)
        beta2 = R_zyz(6)
        gamma2 = R_zyz(7)
        ! *** internal coordinates ***
        ! ----------------------------
        if(xdim==2)then
            internal0(2) = dcos(beta1)
        elseif(xdim==3)then
            internal0(2) = dcos(beta1)
            ! ----------------------------
            internal0(3) = PII - gamma1
        elseif(xdim==4)then
            internal0(2) = dcos(beta1)
            ! ----------------------------
            internal0(3) = dcos(beta2)
            ! ----------------------------
            internal0(4) = alpha1 - alpha2
            if(internal0(4)>PII)then
                internal0(4) = internal0(4) - 2d0 * PII
            endif
            if(internal0(4)<-PII)then
                internal0(4) = internal0(4) + 2d0 * PII
            endif
        elseif(xdim==5)then
            internal0(5) = gamma1
            ! ----------------------------
            internal0(2) = dcos(beta1)
            ! ----------------------------
            internal0(3) = dcos(beta2)
            ! ----------------------------
            internal0(4) = alpha1 - alpha2
            if(internal0(4)>PII)then
                internal0(4) = internal0(4) - 2d0 * PII
            endif
            if(internal0(4)<-PII)then
                internal0(4) = internal0(4) + 2d0 * PII
            endif
        elseif(xdim==6)then
            internal0(5) = gamma1
            ! ----------------------------
            internal0(6) = gamma2
            ! ----------------------------
            internal0(2) = dcos(beta1)
            ! ----------------------------
            internal0(3) = dcos(beta2)
            ! ----------------------------
            internal0(4) = alpha1 - alpha2
            if(internal0(4)>PII)then
                internal0(4) = internal0(4) - 2d0 * PII
            endif
            if(internal0(4)<-PII)then
                internal0(4) = internal0(4) + 2d0 * PII
            endif
        endif

        return
    end subroutine euler_angles_to_autosurf

    subroutine int_to_cart(internal, internal_length, mass, ref1, ref2, xdim, natom1, natom2, input_coordinates, cart)

        implicit none

        integer, intent(in) :: internal_length, natom1, natom2, xdim
        real(real64), intent(in) :: internal(internal_length), mass(natom1 + natom2), ref1(natom1 * 3), ref2(natom2 * 3)
        integer :: intfunc
        real(real64), allocatable, intent(out) :: cart(:)
        real(real64) :: ref1_temp0(natom1 * 3), ref2_temp0(natom2 * 3), cart_temp((natom1 + natom2) * 3)
        character(*), intent(in) :: input_coordinates

        allocate(cart(3 * (natom1 + natom2)))

        !************** subroutine definition
        ! Place your custom function here!
        !********************************************************
        ! make temporary copies of the input data
        ref1_temp0 = ref1! reference vector for frag1 (original frame)
        ref2_temp0 = ref2! reference vector for frag2 (original frame)

        if (input_coordinates == "User")then
            call int_to_cart_user(internal, mass, cart_temp)
        elseif (input_coordinates == "Spherical")then
            ! Only for 3D: (for CH3CN-He system: internal coordinates taken as spherical coords.)
            call int_to_cart_Spherical(internal, cart_temp, mass, natom1, natom2, ref1_temp0)
        elseif (input_coordinates == "Autosurf") then
            call int_to_cart_zyz(internal, xdim, cart_temp, mass, natom1, natom2, ref1_temp0, ref2_temp0)
            ! Only if internal is already in cartesian
        elseif (input_coordinates == "Cartesian")then

            cart_temp = internal
        endif

        cart = cart_temp

    end subroutine int_to_cart


    subroutine zyz_to_output(internal, internal0, internal0_length, R_zyz, ref1, ref2&
            , xdim, natom1, natom2, output_coordinates)

        implicit none

        integer, intent(in) :: natom1, natom2, xdim, internal0_length
        real(real64), intent(in) :: R_zyz(7), ref1(natom1 * 3), ref2(natom2 * 3), internal(xdim)
        real(real64), intent(inout) :: internal0(internal0_length)
        real(real64) :: ref1_temp0(natom1 * 3), ref2_temp0(natom2 * 3)
        character(*), intent(in) :: output_coordinates





        !************** subroutine definition
        ! Place your custom function here!
        !********************************************************
        ! make temporary copies of the input data
        ref1_temp0 = ref1! reference vector for frag1 (original frame)
        ref2_temp0 = ref2! reference vector for frag2 (original frame)

        if (output_coordinates=="Autosurf")then
            call euler_angles_to_autosurf(xdim, R_zyz, internal0)
        elseif (output_coordinates=="BiSpherical")then
            call euler_angles_to_bispherical(internal, R_zyz, internal0)
        elseif (output_coordinates == "Cartesian")then
        elseif (output_coordinates=="User")then
        elseif (output_coordinates== "Spherical")then
        end if

    end subroutine zyz_to_output

    subroutine int_to_cart_zyz(internal, xdim, cart, mass, natom1, natom2, ref1_0, ref2_0)

        !************** subroutine definition 
        ! subroutine int_to_cart_zyz (convert: Internal coordinates in zyz--> Cartesian coordinates)
        ! This function considers local axis parallels to ref1_0 and ref2_0 (systems known by Autosurf)
        ! and global axis with origin in Center Mass of A with z-axis passing through Center Mass of B

        ! * the second angle beta  (rotation around Y axis) is given by cos_beta
        ! * the otehr two, alpha and gamma are given in radians
        !********************************************************
        use math_functions, only : remove_center_of_mass, vec_to_mat2, molecular_rotation_zyz, mat_to_vec2

        implicit none
        integer, intent(in) :: natom1, natom2, xdim
        real(real64), intent(in) :: internal(xdim), mass(natom1 + natom2), ref1_0(natom1 * 3), ref2_0(natom2 * 3)
        real(real64), intent(out) :: cart((natom1 + natom2) * 3)
        integer :: i, natom
        real(real64) :: cm(3)
        real(real64) :: ref1(natom1 * 3), rotarr1(3, natom1), ref_mat1(3, natom1), arr1(natom1 * 3)
        real(real64) :: ref2(natom2 * 3), rotarr2(3, natom2), R_arr(3, natom2), ref_mat2(3, natom2), arr2(natom2 * 3)
        real(real64) :: R, angles(6), internal_temp(xdim)

        internal_temp = internal
        natom = natom1 + natom2
        ref1 = ref1_0
        ref2 = ref2_0

        R = internal_temp(1)

        if (xdim==2)then
            angles(1) = 0d0
            angles(2) = dacos(internal_temp(2))
            angles(3) = 0d0
            angles(4) = 0d0
            angles(5) = 0d0
            angles(6) = 0d0
        elseif (xdim==3)then
            angles(1) = 0d0
            angles(2) = dacos(internal_temp(2))
            angles(3) = internal_temp(3)
            angles(4) = 0d0
            angles(5) = 0d0
            angles(6) = 0d0
        elseif (xdim==4)then
            angles(1) = internal_temp(4)
            angles(2) = dacos(internal_temp(2))
            angles(3) = 0d0
            angles(4) = 0d0
            angles(5) = dacos(internal_temp(3))
            angles(6) = 0d0
        elseif (xdim==5)then
            angles(1) = internal_temp(4)
            angles(2) = dacos(internal_temp(2))
            angles(3) = internal_temp(5)
            angles(4) = 0d0
            angles(5) = dacos(internal_temp(3))
            angles(6) = 0d0
        elseif (xdim==6)then
            angles(1) = internal_temp(4)
            angles(2) = dacos(internal_temp(2))
            angles(3) = internal_temp(5)
            angles(4) = 0d0
            angles(5) = dacos(internal_temp(3))
            angles(6) = internal_temp(6)
        endif



        !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of
        !!!    reference (origin at the new_CM of frag. 1)



        ! set new CM of fragment 1 at the origin
        call remove_center_of_mass(ref1, mass(1:natom1), natom1, natom1)
        call remove_center_of_mass(ref2, mass(natom1 + 1:natom1 + natom2), natom2, natom2)

        call vec_to_mat2(ref1, ref_mat1, natom1)
        call vec_to_mat2(ref2, ref_mat2, natom2)

        call molecular_rotation_zyz(ref_mat1, natom1, angles(1), angles(2), angles(3), rotarr1)
        call molecular_rotation_zyz(ref_mat2, natom2, angles(4), angles(5), angles(6), rotarr2)

        call mat_to_vec2(rotarr1, arr1, natom1)
        R_arr(1:2, :) = rotarr2(1:2, :)
        R_arr(3, :) = rotarr2(3, :) + R
        call mat_to_vec2(R_arr, arr2, natom2)

        do i = 1, 3 * natom1
            cart(i) = arr1(i)! Cartesian coordinates for the atoms in fragment 1
        enddo
        do i = 1, 3 * natom2
            cart(3 * natom1 + i) = arr2(i)! Cartesian coordinates for the atoms in fragment 2
        enddo

    end subroutine int_to_cart_zyz

    subroutine int_to_cart_Spherical(internal, cart, mass, natom1, natom2, ref1_0)
        use math_functions, only : remove_center_of_mass
        implicit none
        integer, intent(in) :: natom1, natom2
        real(real64), intent(in) :: internal(3), mass(natom1 + natom2), ref1_0(natom1 * 3)

        integer :: i, natom
        real(real64) :: cm(3)
        real(real64) :: ref1(natom1 * 3)
        real(real64) :: ref2(natom2 * 3)
        real(real64) :: cart((natom1 + natom2) * 3)
        real(real64) :: sin_theta

        natom = natom1 + natom2
        ref1 = ref1_0

        !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of
        !!!    reference (origin at the new_CM of frag. 1)

        ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
        ! (for CH3CN-He system: internal coordinates taken as spherical coords.)

        ! set new CM of fragment 1 at the origin
        call remove_center_of_mass(ref1, mass(1:natom1), natom1, natom1)
        do i = 1, 3 * natom1
            cart(i) = ref1(i)! Cartesian coordinates for the atoms in fragment 1
        enddo
        sin_theta = dsqrt(1.d0 - internal(2)**2)
        cm(1) = internal(1) * sin_theta * dcos(internal(3))
        cm(2) = internal(1) * sin_theta * dsin(internal(3))
        cm(3) = internal(1) * internal(2)
        ref2 = cm
        do i = 1, 3
            cart(3 * natom1 + i) = ref2(i)! Cartesian coordinates for the atoms in fragment 2
        enddo

    end subroutine int_to_cart_Spherical

    ! int_to_cart_user: Example for H2O-HCN
    subroutine int_to_cart_user(internal, mass, cart)
        implicit none
        integer, parameter :: xdim = 5, natom1 = 3, natom2 = 3
        real(real64), parameter :: bohr2ang = 0.529177249d0
        real(real64), intent(in) :: internal(xdim)
        real(real64), intent(in) :: mass(natom1 + natom2)
        real(real64), intent(out) :: cart((natom1 + natom2) * 3)
        real(real64) :: rOH, thetaOH, Xcm, Ycm, Zcm, R, theta1, theta2, phi1, phi2
        real(real64) :: r1, r2, rH0, rC1, rN1
        real(real64) :: xO, yO, zO, xH1, yH1, zH1, xH2, yH2, zH2
        real(real64) :: xH0, yH0, zH0, xC1, yC1, zC1, xN1, yN1, zN1
        real(real64) :: mH1, mH2, mH3, mN, mO, mC, mtA, mtB
        real(real64) :: RCMx, RCMy, RCMz
        real(real64), parameter :: PII = dacos(-1d0)

        ! Masses
        ! -------------------
        ! H2O

        mO = mass(1)!15.99491461956d0
        mH1 = mass(2)!1.00782503207d0
        mH2 = mass(3)!1.00782503207d0
        ! HCN
        mH3 = mass(4)!1.00782503207d0
        mC = mass(5)!12.0d0
        mN = mass(6)!14.0030740048d0
        ! -------------------
        mtA = mH1 + mH2 + mO
        mtB = mN + mC + mH3

        ! H2O parameters
        rOH = 1.8361d0              ! in bohr
        rOH = rOH * bohr2ang          ! convert to Ang.
        thetaOH = 104.69d0 / 2.0d0    ! degrees
        thetaOH = thetaOH * PII / 180d0 ! convert to radians
        ! HCN parameters
        r1 = 2.0286d0 * bohr2ang      ! equilibre CH
        r2 = 2.1874d0 * bohr2ang      ! equilibre CN
        ! distance from each atom to the HCN center-of-mass
        rH0 = ((r1 + r2) * mN + r1 * mC) / (mtB)
        rC1 = (r2 * mN - r1 * mH3) / (mtB)
        rN1 = ((r1 + r2) * mH3 + r2 * mC) / (mtB)

        ! H2O Cartesian coordinates
        ! ----------------------------
        ! (H2O is fixed at the origin)
        xO = 0d0
        yO = 0d0
        zO = (mH1 + mH2) * rOH * dcos(thetaOH) / (mtA)
        xH1 = -1d0 * rOH * dsin(thetaOH)
        yH1 = 0d0
        zH1 = zO - rOH * dcos(thetaOH)
        xH2 = rOH * dsin(thetaOH)
        yH2 = 0d0
        zH2 = zH1

        RCMx = (mO * xO + xH1 * mH1 + xH2 * mH2) / mtA
        RCMy = (mO * yO + yH1 * mH1 + yH2 * mH2) / mtA
        RCMz = (mO * zO + zH1 * mH1 + zH2 * mH2) / mtA

        cart(1) = xO - RCMx
        cart(2) = yO - RCMy
        cart(3) = zO - RCMz
        cart(4) = xH1 - RCMx
        cart(5) = yH1 - RCMy
        cart(6) = zH1 - RCMz
        cart(7) = xH2 - RCMx
        cart(8) = yH2 - RCMy
        cart(9) = zH2 - RCMz
        !write(6,*)cart(1:3)
        !write(6,*)cart(4:6)
        !write(6,*)cart(7:9)

        ! internal coordinates
        R = internal(1)
        theta1 = internal(2)
        theta2 = internal(3)
        phi1 = internal(4)
        phi2 = internal(5)

        ! HCN Cartesian coordinates
        ! -------------------------
        Xcm = R * dsin(theta1) * dcos(phi1)
        Ycm = R * dsin(theta1) * dsin(phi1)
        Zcm = R * dcos(theta1)
        xH0 = Xcm - 1.0d0 * rH0 * dsin(theta2) * dcos(phi2)
        yH0 = Ycm - 1.0d0 * rH0 * dsin(theta2) * dsin(phi2)
        zH0 = Zcm - 1.0d0 * rH0 * dcos(theta2)
        xC1 = Xcm - 1.0d0 * rC1 * dsin(theta2) * dcos(phi2)
        yC1 = Ycm - 1.0d0 * rC1 * dsin(theta2) * dsin(phi2)
        zC1 = Zcm - 1.0d0 * rC1 * dcos(theta2)
        xN1 = Xcm + 1.0d0 * rN1 * dsin(theta2) * dcos(phi2)
        yN1 = Ycm + 1.0d0 * rN1 * dsin(theta2) * dsin(phi2)
        zN1 = Zcm + 1.0d0 * rN1 * dcos(theta2)
        cart(10) = xH0
        cart(11) = yH0
        cart(12) = zH0
        cart(13) = xC1
        cart(14) = yC1
        cart(15) = zC1
        cart(16) = xN1
        cart(17) = yN1
        cart(18) = zN1
        ! write(6,*)
        ! write(6,*)cart(3*natom1+1:3*natom1+3*natom2)
        !pause

        return

    end subroutine int_to_cart_user

end module coordinateTransf
    
module testing_auxiliar_functions
    use mod_types, only : int32, real64
contains

    subroutine check_cartesian_frames(  ptos, &
                                        cart_model, &
                                        cart, &
                                        natom1, &
                                        natom2, &
                                        counterCase, &
                                        test_failed, &
                                        maxerr, &
                                        err_tolerance, &
                                        xdim, &
                                        internal0, &
                                        file_output_number)
        use math_functions,only: cosine_law
        implicit none
        integer (int32), intent (inout) :: counterCase, test_failed
        integer (int32), intent (in) :: xdim, ptos(3), natom1, natom2
        real (real64), intent (in) :: err_tolerance, internal0(xdim)
        real (real64), intent (in) :: cart_model((natom1 + natom2) * 3), cart((natom1 + natom2) * 3)
        real (real64), intent (inout) :: maxerr
        integer (int32) :: i, j, k
        real (real64) :: err_test, modelA(6), modelB(6), A(6), B(6)
        integer (int32), optional :: file_output_number
        counterCase = counterCase + 1

        i = ptos(1)
        j = ptos(2)
        k = ptos(3)

        call cosine_law( cart_model(1 + 3 * (i - 1):3 * i), &
                        cart_model(1 + 3 * (j - 1):3 * j), &
                        cart_model(1 + 3 * (k - 1):3 * k), &
                        modelA)
        call cosine_law( cart_model(3 * natom1 + 1 + 3 * (i - 1):3 * natom1 + 3 * i), &
                        cart_model(3 * natom1 + 1 + 3 * (j - 1):3 * natom1 + 3 * j), &
                        cart_model(3 * natom1 + 1 + 3 * (k - 1):3 * natom1 + 3 * k), &
                        modelB)

        call cosine_law( cart(1 + 3 * (i - 1):3 * i), cart(1 + 3 * (j - 1):3 * j), cart(1 + 3 * (k - 1):3 * k), A)
        call cosine_law( cart(3 * natom1 + 1 + 3 * (i - 1):3 * natom1 + 3 * i),&
                        cart(3 * natom1 + 1 + 3 * (j - 1):3 * natom1 + 3 * j),&
                        cart(3 * natom1 + 1 + 3 * (k - 1):3 * natom1 + 3 * k), B)

        err_test = max(Norm2(A - modelA), Norm2(B - modelB))

        if (maxerr <  err_test) then
            maxerr = err_test
        endif

        if(err_test > err_tolerance) then

            test_failed = test_failed + 1

            if (present(file_output_number))then
                write(file_output_number, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(file_output_number, *)"Cartesians Rigid Molec Checker            : ", counterCase
                write(file_output_number, *)"Error           : ", err_test
                write(file_output_number, *)"Norm2 Molec A: ", Norm2(A - modelA), "Norm2 Molec B: ", Norm2(B - modelB), &
                        "Norm2 AB : ", Norm2(A - B)

                write(file_output_number, *)"Model Molec A: ", modelA
                write(file_output_number, *)"Molec A: ", A
                write(file_output_number, *)"Model Molec B: ", modelB
                write(file_output_number, *)"Molec B: ", B
                write(file_output_number, *)"Dimension :", xdim
                write(file_output_number, *) internal0
                write(file_output_number, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(file_output_number, *)

            else
                write(*, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(*, *)"Cartesians Rigid Molec Checker           : ", counterCase
                write(*, *)"Error           : ", err_test
                write(*, *)"Norm2 Molec A: ", Norm2(A - modelA), "Norm2 Molec B: ", Norm2(B - modelB), &
                        "Norm2 AB : ", Norm2(A - B)

                write(*, *)"Model Molec A: ", modelA
                write(*, *)"Molec A: ", A
                write(*, *)"Model Molec B: ", modelB
                write(*, *)"Molec B: ", B
                write(*, *)"Dimension :", xdim
                write(*, *) internal0
                write(*, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(*, *)
            end if

        end if

    end subroutine
    !Reporting Error

    subroutine report_cartesian_error(counterCase, test_failed, maxerr, td, err_tolerance, natoms, xdim, cart, internal0&
            , file_output_number)

        implicit none
        integer (int32), intent (inout) :: counterCase, test_failed
        integer (int32), intent (in) :: natoms, xdim
        real (real64), intent (in) :: td(4), err_tolerance
        real (real64), intent (in) ::  cart(natoms * 3), internal0(xdim)
        real (real64), intent (inout) :: maxerr
        integer (int32), optional :: file_output_number
        
        integer (int32) :: j
        real (real64) :: err_test
        counterCase = counterCase + 1

        err_test = 0;
        err_test = sum(DABS(td))

        if (maxerr < err_test) then
            maxerr = err_test
        endif

        if (err_test > err_tolerance)Then
            test_failed = test_failed + 1

            if (present(file_output_number))then
                write(file_output_number, *)
                write(file_output_number, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(file_output_number, *)"Cartesians-Autosurf  InterAtomic-Distance Test          : ", counterCase
                write(file_output_number, *)"Error           : ", err_test, td
                do j = 1, natoms
                    write(file_output_number, *)"Atom ", j, " :", cart(1 + (j - 1) * 3:j * 3)
                    write (file_output_number, *)"----------------------------------------------"
                end do

                write(file_output_number, *)"Dimension :", xdim
                write(file_output_number, *) internal0
                write(file_output_number, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(file_output_number, *)

            else

                write(*, *)
                write(*, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(*, *)"Cartesians-Autosurf  InterAtomic-Distance Test           : ", counterCase
                write(*, *)"Error           : ", err_test, td
                do j = 1, natoms
                    write(*, *)"Atom ", j, " :", cart(1 + (j - 1) * 3:j * 3)
                    write (*, *)"----------------------------------------------"
                end do

                write(*, *)"Dimension :", xdim
                write(*, *) internal0
                write(*, *)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%", '\n'
                write(*, *)
            end if
        end if
    end subroutine report_cartesian_error


    subroutine test_cart_to_autosurf(system_path, data_path, file_output_number)
        !use math_functions,only: check
        use helper_functions,only : get_pes_coordinates

        implicit none
        character(*), intent (in) :: system_path, data_path
        integer (int32), optional :: file_output_number

        real (real64) :: cart(18), cart_model(18), energies(4)
        character(len = 1) :: Atom_label
        real (real64) :: td(4), err_test, err_tolerance, maxerr_1, maxerr_2
        integer (int32) :: stat2, stat3, i, natoms, int_to_cart_func, xdim, ptos(3)
        integer (int32) :: failedTest_1, counterCase_1, failedTest_2, counterCase_2
        real (real64) :: response1_model(6), response2_model(6), response1(6), response2(6)
        integer (int32) :: stat
        real (real64), allocatable :: internal0(:)

        xdim = 6

        open(unit = 100, file = data_path, status = 'old', action = 'read')

        failedTest_1 = 0
        counterCase_1 = 0
        failedTest_2 = 0
        counterCase_2 = 0
        maxerr_1 = -1d0
        maxerr_2 = -1d0
        err_tolerance = 1d-7

        ptos = (/1, 2, 3 /)
        i = 0

        do i = 1, 1000

            read(100, *, iostat = stat)natoms
            if (stat /= 0) exit

            read(100, *)energies(1:4)
            read(100, *)Atom_label, cart(1:3)
            read(100, *)Atom_label, cart(4:6)
            read(100, *)Atom_label, cart(7:9)
            read(100, *)Atom_label, cart(10:12)
            read(100, *)Atom_label, cart(13:15)
            read(100, *)Atom_label, cart(16:18)

            !-1 is for cartesian input

            if (i==1)then
                cart_model = cart
            end if

            call get_pes_coordinates(internal0, cart, size(cart), 6, "Cartesian", "Autosurf", system_path, errors = td)
            !  internal,internalLength,internal0_,xdim,inputCoord,outputCoord,filePath,&
            !                                     testArr_Errors,newflag


            !Reporting Error

            call report_cartesian_error(counterCase_1, failedTest_1, maxerr_1, td, err_tolerance, natoms, xdim, cart, internal0, &
                    file_output_number)
            call check_cartesian_frames(ptos, cart_model, cart, 3, 3, counterCase_2, failedTest_2, maxerr_2&
                    , err_tolerance, xdim, internal0, file_output_number)
            !end Error Testing

        enddo

        if (present(file_output_number))then
            write(file_output_number, *)" "
            write(file_output_number, *) "%%%%%%%%%%%%%%%    Cartesian Test    %%%%%%%%%%%%%%%"
            write(file_output_number, *) "Number of failed Cartesian coordinates: ", failedTest_1, " out of ", &
                    counterCase_1, ", maxErr: ", maxerr_1, ' '
            write(file_output_number, *) "Number of failed check_cartesian_frames: ", failedTest_2, " out of ", &
                    counterCase_2, ", maxErr: ", maxerr_2, ' '
            write(file_output_number, *)" "
        endif

        call print_test_results("Cartesians-Autosurf  InterAtomic-Distance Test", maxerr_1, "Cartesian", "Autosurf", &
                counterCase_1, counterCase_1 - failedTest_1)
        call print_test_results("Cartesians Rigid Molecule Checker                ", maxerr_2, "Cartesian", "Autosurf", &
                counterCase_2, counterCase_2 - failedTest_2)

        close(100)


    end subroutine test_cart_to_autosurf


    subroutine print_test_results(testName, maxErr, inputCoord, outputCoord, ntest, success_test)

        implicit none
        character(*), intent (in) :: testName, inputCoord, outputCoord
        integer (int32), intent (in) :: ntest, success_test
        character (len = 25) :: f_res1, f_res2
        character (len = 7) :: str, str1
        integer (int32) :: failed, ptg
        real (real64), intent (in) :: maxErr

        f_res1 = achar(27) // '[32m TEST PASSED ' // achar(27) // '[0m'
        f_res2 = achar(27) // '[31m TEST FAILED :(' // achar(27) // '[0m'

        ptg = FLOOR(success_test * 100d0 / ntest)
        failed = ntest - success_test
        write (str1, '(I4,A)') ptg, "%"

        write(*, *)
        if (ntest > success_test)then

            write(*, '(A, I5, A, I1, A, I5, A)') f_res2 // "(" // achar(27) // '[31m' // str1 // achar(27) // '[0m)' // " ----  "&
                    // achar(27) // '[35m' // testName // achar(27) // '[0m)'
        else
            write(*, *) f_res1 // "(" // achar(27) // '[32m 100% ' // achar(27) // '[0m)' // " ----  "&
                    // achar(27) // '[35m' // testName // achar(27) // '[0m)'
        end if
        write(*, '(A, I8, A, I8, A, I8)') "                                Num. Test: ", ntest, "     ---> success: "&
                , success_test, " / failed : ", failed
        write(*, '(A, F21.15)') "                                Error Max: ", maxErr
        write(*, *) "                               Input/Output Coordinates: " // achar(27) // '[33m' // inputCoord // " / "&
                // outputCoord // achar(27) // '[0m)'
        write(*, *)

    end subroutine print_test_results


    subroutine generate_random_data(internal, &
                                    szi, &
                                    rms, &
                                    new_mass, &
                                    mass, &
                                    mass0, &
                                    natom1, &
                                    natom2, &
                                    ref1_0, &
                                    ref2_0, &
                                    xdim, &
                                    inputCoord, &
                                    ntest)
        use helper_functions, only:
        use coordinateTransf,only: int_to_cart

        implicit none
        integer (int32), intent (in) :: xdim, rms, natom1, natom2, ntest, szi ! rm = random mass
        real (real64), intent (out) :: internal(szi, ntest)
        real (real64), intent (out) :: new_mass(natom1 + natom2, ntest)
        real (real64) :: PII, threshold, maxdist, passingThrough, iter(xdim), rmin, rmax
        real (real64), intent (in) :: mass(natom1 + natom2), mass0(natom1 + natom2), ref1_0(natom1 * 3), ref2_0(natom2 * 3)
        character(*), intent (in) :: inputCoord
        real (real64), allocatable :: rand_cases(:, :), cart(:)
        real (real64) :: ref1_temp0(natom1 * 3), ref2_temp0(natom2 * 3), nm(natom1 + natom2, ntest), Autosurf(xdim, ntest)

        integer (int32) :: i, natom, sz, nseeds

        natom = natom1 + natom2
        PII = acos(-1d0)
        sz = rms * (natom) + xdim
        rmin = 5d0
        rmax = 15d0
        threshold = 0.9993d0

        if (rms==1)then

            sz = rms * (natom) + xdim
        else
            sz = xdim
        end if
        ! make temporary copies of the input data
        ref1_temp0 = ref1_0! reference vector for frag1 (original frame)
        ref2_temp0 = ref2_0! reference vector for frag2 (original frame)

        allocate(rand_cases(sz, ntest))

        call random_seed()
        call random_number(rand_cases)

        nseeds = ntest + 0

        call MassGenerator(rms, rand_cases, sz, nseeds, xdim, natom, mass, mass0, nm)

        new_mass = nm

        do i = 1, ntest
            !write(*,*)(sz-xdim)+1,sz,rand_cases(sz-xdim+1:sz,i)
            !iter = rand_cases(sz-xdim+1:sz,i)
            Autosurf(1, i) = rand_cases(sz - xdim + 1, i) * 10d0 + 5d0
            Autosurf(2, i) = (rand_cases(sz - xdim + 2, i) * 2d0 - 1d0) * threshold
            if (xdim > 2) then
                if (xdim == 3) then
                    Autosurf(3, i) = rand_cases(sz - xdim + 3, i) * 2d0 * PII
                else
                    Autosurf(3, i) = (rand_cases(sz - xdim + 3, i) * 2d0 - 1d0) * threshold
                end if
            end if
            if (xdim>3)then
                Autosurf(4, i) = rand_cases(sz - xdim + 4, i) * 2d0 * PII
            end if
            if (xdim>4)then
                Autosurf(5, i) = rand_cases(sz - xdim + 5, i) * 2d0 * PII
            end if
            if (xdim>5)then
                Autosurf(6, i) = rand_cases(sz - xdim + 6, i) * 2d0 * PII
            end if
        end do

        if (inputCoord=="Autosurf")then
            internal = Autosurf
        elseif (inputCoord == "Cartesian")then
            do i = 1, ntest
                call int_to_cart( Autosurf(:, i), &
                     xdim, &
                     new_mass(:, i), &
                     ref1_temp0, &
                     ref2_temp0, &
                     xdim, &
                     natom1, &
                     natom2, &
                     "Autosurf", &
                     cart)
                internal(:, i) = cart
            end do
        endif

        deallocate(rand_cases)

    end subroutine generate_random_data

    subroutine MassGenerator(rm, seedArr, sz1, sz2, xdim, natom, mass, mass0, nmass)
        implicit none
        integer (int32), intent (in) :: sz1, sz2, natom, rm, xdim
        real (real64), intent (in) :: seedArr(sz1, sz2), mass(natom), mass0(natom)
        real (real64), intent (inout) :: nmass(natom, sz2)
        real (real64) :: rand(natom), fact_mass
        integer (int32) :: i, j

        fact_mass = 0.2d0 ! maximun percentage of mass changing

        do i = 1, sz2

            if (rm==1)then

                rand = seedArr(1:natom, i) * 2d0 - 1d0
                do j = 1, natom
                    nmass(j, i) = mass0(j) * (1d0 + rm * rand(j) * fact_mass)
                end do
            elseif (rm==0)then

                nmass(1:natom, i) = mass0(1:natom)
            elseif (rm==-1)then

                nmass(:, i) = mass
            end if
        end do

    end subroutine MassGenerator


    subroutine generate_random_data2(internal, rm, mass, mass0, natom1, natom2, ref1_0, ref2_0, xdim, inputCoord, ntest)
        use helper_functions

        implicit none
        integer (int32), intent (in) :: xdim, rm, natom1, natom2, ntest ! rm = random mass
        real (real64), allocatable, intent (out) :: internal(:, :)
        real (real64) :: PII, threshold, maxdist, passingThrough
        real (real64), intent (in) :: mass(natom1 + natom2), mass0(natom1 + natom2), ref1_0(natom1 * 3), ref2_0(natom2 * 3)
        character(*), intent (in) :: inputCoord
        real (real64), allocatable :: rand_cases(:, :)

        allocate(rand_cases(rm * (natom1 + natom2) + xdim, ntest + 2000))

        if (inputCoord == "Cartesian")then
            allocate(internal(3 * (natom1 + natom2), ntest + 2000))
        else
            allocate(internal(xdim, ntest + 2000))
        end if

        PII = acos(-1d0)

        call random_seed()
        call random_number(rand_cases)

    end subroutine generate_random_data2


    subroutine interatomic_distance_test(filename, xdim, inputCoord, outputCoord, test_failed, file_output_number)
        use helper_functions,only:get_pes_coordinates

        implicit none
        integer (int32) :: i, k, nc
        real (real64), allocatable :: internal(:), internal0(:)
        real (real64) :: PII, threshold, maxdist, passingThrough
        real (real64), allocatable :: rand_cases(:, :)
        real (real64) :: Max_test_dist, testArr_Errors(4), err
        integer (int32) :: counterCase, test_number
        character(*), intent (in) :: filename, inputCoord, outputCoord
        integer (int32), intent (out) :: test_failed
        integer (int32), intent (in) :: xdim
        integer (int32), optional :: file_output_number
        character(len = 50) :: testName

        test_number = 1000

        allocate(rand_cases(xdim, test_number + 2000))
        PII = acos(-1d0)
        Max_test_dist = 0
        err = 10d0**(-8)
        counterCase = 0

        call random_seed()
        call random_number(rand_cases)

        allocate(internal(xdim))!,internal0(xdim))

        threshold = 0.9993d0
        nc = 0

        counterCase = 0
        test_failed = 0
        Max_test_dist = 0

        do while (counterCase < test_number)
            nc = nc + 1
            passingThrough = 0d0

            internal(1) = rand_cases(1, nc) * 10d0 + 5d0
            internal(2) = (rand_cases(2, nc) * 2d0 - 1d0) * threshold
            if (xdim > 2) then
                if (xdim == 3) then
                    internal(3) = rand_cases(3, nc) * 2d0 * PII
                else
                    internal(3) = (rand_cases(3, nc) * 2d0 - 1d0) * threshold
                end if
            end if
            if (xdim>3)then
                internal(4) = rand_cases(4, nc) * 2d0 * PII
            end if
            if (xdim>4)then
                internal(5) = rand_cases(5, nc) * 2d0 * PII
            end if
            if (xdim>5)then
                internal(6) = rand_cases(6, nc) * 2d0 * PII
            end if

            call get_pes_coordinates(internal0, internal, size(internal), xdim, inputCoord, outputCoord, filename, testArr_Errors)

            passingThrough = 0d0
            if (dabs(internal0(2))< threshold) then
                passingThrough = 1d0
                if (size(internal0)>3)then
                    if (dabs(internal0(3)) < threshold)then
                        passingThrough = 1d0
                    else
                        passingThrough = 0d0
                    end if
                end if

            endif

            if (passingThrough>0d0) then

                counterCase = counterCase + 1
                maxdist = MAXVAL(testArr_Errors);
                if (maxdist > Max_test_dist) then
                    Max_test_dist = maxdist
                endif

                if (maxdist > err) then
                    write(file_output_number, *)"--------------------------------------- "
                    write(file_output_number, *)"Intermolecular Distances Test : ", testArr_Errors
                    write(file_output_number, *)"Internal : ", internal
                    write(file_output_number, *)"Internal0 : ", internal0
                    write(file_output_number, *)"----------------------------------------"
                    test_failed = test_failed + 1

                endif
            endif
            !end Testing Section
        enddo

        !testName = "Interatomic Distances("xdim//"D), ifun: "//internalFunction
        write(testName, '(A(I1)A)') "Interatomic Distances(", xdim, "D)"
        testName = trim(testName);

        call print_test_results(testName, Max_test_dist, inputCoord, outputCoord, counterCase, counterCase - test_failed)
        write(file_output_number, *)
        write(file_output_number, *)"TEST NAME: ", testName
        write(file_output_number, *)"System Features : Dimension/ Number of Atoms", xdim
        write(file_output_number, *)"Max_test_dist : ", Max_test_dist
        write(file_output_number, *)"Input/Output Coordinates: ", inputCoord, " / ", outputCoord
        write(file_output_number, *)"Number of Tests", counterCase
        write(file_output_number, *)"Failed Tests: ", test_failed, " out of", test_number
        write(file_output_number, *)

    end subroutine interatomic_distance_test


    subroutine interatomic_distance_test_v2(filename, xdim, inputCoord, outputCoord, test_failed, rm, file_output_number)

        use helper_functions,only: read_file,convert_isotopic_coordinates
        use math_functions!,only: genera
        use coordinateTransf,only:zyz_to_output,int_to_cart

        implicit none

        character(*), intent (in) :: filename, inputCoord, outputCoord
        integer (int32), intent (out) :: test_failed
        integer (int32), intent (in) :: xdim
        integer (int32), intent (in) :: rm! random mass switcher  rm =1 random mass, rm =0 new mass = mass0, rm =-1 new_mass=mass of file
        integer (int32), optional :: file_output_number

        integer (int32) :: natom1, natom2, natom
        real (real64), allocatable :: ref1_0(:), ref2_0(:), new_mass(:, :), nmass(:), mass(:), mass0(:), &
                ref1_temp0(:), ref2_temp0(:), cart(:)
        real (real64) :: R_ZYZ(7)
        integer (int32) :: xdim_file, internalLength
        integer (int32) :: i, k, nc
        real (real64), allocatable :: internal(:, :), internal0(:)
        real (real64) :: PII, threshold, maxdist, passingThrough
        real (real64) :: Max_test_dist, testArr_Errors(4), err
        character(len = 50) :: testName
        character(len = 20) :: massStatement
        real (real64) :: td(4), internal_i(xdim)
        integer (int32) :: counterCase, ntest, internal0_length, szi
        real (real64), allocatable :: internal0_(:), new_mass_(:, :)

        natom = natom1 + natom2
        threshold = 0.9993d0
        nc = 0
        PII = acos(-1d0)
        err = 10d0**(-8)
        counterCase = 0
        test_failed = 0
        Max_test_dist = 0

        call read_file(filename, mass, mass0, natom1, natom2, ref1_0, ref2_0, xdim_file)

        ntest = 1000000

        allocate(new_mass(natom1 + natom2, ntest), new_mass_(natom1 + natom2, ntest))

        if (inputCoord == "Cartesian")then
            szi = 3 * (natom1 + natom2)
        else
            szi = xdim
        endif

        allocate(internal(szi, ntest))

        call generate_random_data(  internal, &
                                    szi, &
                                    rm, &
                                    new_mass, &
                                    mass, &
                                    mass0, &
                                    natom1, &
                                    natom2, &
                                    ref1_0, &
                                    ref2_0, &
                                    xdim, &
                                    inputCoord, &
                                    ntest)

        allocate(ref1_temp0(natom1 * 3), ref2_temp0(natom2 * 3), nmass(natom))

        ! make temporary copies of the input data
        ref1_temp0 = ref1_0! reference vector for frag1 (original frame)
        ref2_temp0 = ref2_0! reference vector for frag2 (original frame)

        ! make sure the original Cartesian frame of reference is at the CM of each fragment
        if(natom1>1)call remove_center_of_mass(ref1_temp0, mass0(1:natom1), natom1, natom1)
        if(natom2>1)call remove_center_of_mass(ref2_temp0, mass0(natom1 + 1:natom), natom2, natom2)

        if (outputCoord == "Cartesian")then
            internal0_length = 3 * (natom1 + natom2)
        else
            internal0_length = xdim
        endif

        allocate(internal0_(internal0_length), internal0(internal0_length))

        do while (nc < ntest)

            nc = nc + 1
            passingThrough = 0d0
            internal_i = internal(:, nc)
            nmass = new_mass(:, nc)

            call int_to_cart(internal_i, szi, nmass, ref1_temp0, ref2_temp0, xdim, natom1, natom2, inputCoord, cart)

            call convert_isotopic_coordinates(  cart, &
                                                R_ZYZ, &
                                                nmass, &
                                                mass0, &
                                                natom1, &
                                                natom2, &
                                                ref1_temp0, &
                                                ref2_temp0, &
                                                xdim, &
                                                td, &
                                                1)

            call zyz_to_output( internal_i, &
                                internal0_, &
                                internal0_length, &
                                R_ZYZ, &
                                ref1_temp0, &
                                ref2_temp0, &
                                xdim, &
                                natom1, &
                                natom2, &
                                outputCoord)

            internal0_ = internal0

            passingThrough = 0d0

            if (outputCoord =="Cartesian")then
            else

                if (dabs(internal0(2))< threshold) then
                    passingThrough = 1d0
                    if (size(internal0)>3)then
                        if (dabs(internal0(3)) < threshold)then
                            passingThrough = 1d0
                        else
                            passingThrough = 0d0
                        end if
                    end if

                endif

                if (passingThrough>0d0) then

                    natom = natom1 + natom2

                    call read_file(filename, mass, mass0, natom1, natom2, ref1_0, ref2_0, xdim_file)


                    counterCase = counterCase + 1
                    maxdist = MAXVAL(td);
                    if (maxdist > Max_test_dist) then
                        Max_test_dist = maxdist
                    endif

                    if (maxdist > err) then
                        write(file_output_number, *)"--------------------------------------- "
                        write(file_output_number, *)testName, td
                        write(file_output_number, *)"Internal : ", internal
                        write(file_output_number, *)"Internal0 : ", internal0
                        write(file_output_number, *)"----------------------------------------"
                        test_failed = test_failed + 1

                    endif
                endif
                !end Testing Section

            endif

        end do

        if (rm==-1)then
            massStatement = "-File Masses"
        elseif(rm==0)then
            massStatement = "-Same Masses"
        else
            massStatement = "-Random Masses"
        endif

        write(testName, '(A(I1)A)') "Interatomic Distances(", xdim, "D)" // massStatement
        testName = trim(testName);
        call print_test_results(testName, Max_test_dist, inputCoord, outputCoord, counterCase, counterCase - test_failed)

        write(file_output_number, *)
        write(file_output_number, *)"TEST NAME: ", testName
        write(file_output_number, *)"System Features : Dimension/ Number of Atoms", xdim
        write(file_output_number, *)"Max_test_dist : ", Max_test_dist
        write(file_output_number, *)"Input/Output Coordinates: ", inputCoord, " / ", outputCoord
        write(file_output_number, *)"Number of Tests", counterCase
        write(file_output_number, *)"Failed Tests: ", test_failed, " out of", ntest
        write(file_output_number, *)

        deallocate(ref1_0, ref2_0, new_mass, mass, mass0, ref1_temp0, ref2_temp0, cart, internal0_, nmass)


    end subroutine interatomic_distance_test_v2


     subroutine axis_change_test(filename,inputCoord,test_failed,file_output_number)
         use helper_functions,only : convert_isotopic_coordinates
         use math_functions, only : remove_center_of_mass, cosine_law,center_of_mass_v2
         use coordinateTransf, only: int_to_cart
         implicit none
         character(*,kind=1),intent (in)::filename
         character(*,kind=1),intent (in):: inputCoord
         integer (int32),intent (out)::test_failed
         integer (int32),intent (in)::file_output_number
         real (real64), allocatable ::ref1_temp0(:),ref2_temp0(:),cart(:)
         real (real64)::cmA(3),cmA0(3),cmB(3),cmB0(3)

         integer (int32) :: natom1,natom2,natom,xdim,i,k,nc,szi
         real (real64), allocatable :: ref1(:),ref2(:),mass(:),mass0(:),internal(:),internal0(:)
         real (real64) :: threshold,maxdist,passingThrough
         real (real64), allocatable :: rand_cases(:,:),MassRand(:)
         real (real64) :: Max_test_dist,test_dist(5),err,dcmA(3),dcmB(3)
         integer (int32) :: counterCase, test_number
         real(real64),parameter:: PII=acos(-1d0)

         real (real64):: R0,R_,cb1,cb10,cb2,cb20,R0_test,cb1_test,cb10_test,cb2_test,cb20_test
         real (real64):: R_error,cb10_error,cb1_error,cb2_error,cb20_error,O(3),O0(3),B1(3)
         real (real64)::response1(6),response2(6),response3(6)

         test_number=10000
         open(unit=10,file=filename,status='old',action='read')

         read(10,*) xdim
         read(10,*) natom1! number of atoms in fragment1
         read(10,*) natom2! number of atoms in fragment2

         allocate(rand_cases(xdim,test_number + 5000),MassRand(test_number + 5000))
         
         Max_test_dist =0d0
         err = 10d0**(-8)
         counterCase = 0

         call random_seed()
         call random_number(rand_cases)
         call random_number(MassRand)
         natom=natom1+natom2 ! total number of atoms

         allocate(ref1(3*natom1),ref2(3*natom2),mass(natom),mass0(natom),internal(xdim),internal0(xdim))
         allocate(ref1_temp0(3*natom1),ref2_temp0(3*natom2),cart(natom*3))

         if (inputCoord == "Cartesian")then
             szi = 3 * (natom1 + natom2)
         else
             szi = xdim
         endif
         
         if(natom1>1)then
         do i=1,natom1
             read(10,*)mass0(i),(ref1(k),k=i*3-2,i*3)
         enddo
         elseif(natom1==1)then
         read(10,*)mass0(1)
         ref1=0.d0
         endif

         if(natom2>1)then
         do i=1,natom2
             read(10,*)mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
         enddo
         elseif(natom2==1)then
         read(10,*)mass0(natom1+1)
         ref2=0.d0
         endif
         close(10)

         mass=mass0


             threshold = 0.9993d0
             nc = 0

             counterCase = 0
             test_failed = 0
             Max_test_dist = 0



             do while (counterCase < test_number)

                 nc=nc+1
                 passingThrough =0d0

                 mass(1)=mass0(1)+MassRand(nc)*10d0

                 internal(1)=rand_cases(1,nc)*10d0 + 5d0
                 internal(2)=(rand_cases(2,nc)*2d0-1d0)*threshold

                 if (xdim > 2) then
                     if (xdim == 3) then
                         internal(3)=rand_cases(3,nc)*2d0*PII
                     else
                         internal(3)=(rand_cases(3,nc)*2d0-1d0)*threshold
                     end if
                 end if
                 if (xdim>3)then
                     internal(4)=rand_cases(4,nc)*2d0*PII
                 end if
                 if (xdim>4)then
                     internal(5)=rand_cases(5,nc)*2d0*PII
                 end if
                 if (xdim>5)then
                     internal(6)=rand_cases(6,nc)*2d0*PII
                 end if





                call convert_isotopic_coordinates(internal,&
                                                internal0,&
                                                mass,&
                                                mass0,&
                                                natom1,&
                                                natom2,&
                                                ref1,&
                                                ref2,&
                                                xdim,&
                                                test_dist,&
                                                1)



                passingThrough = 0d0
                 if (dabs(internal0(2))< threshold  .and. test_dist(5)>-1 ) then
                     passingThrough = 1d0
                     if (size(internal0)>3)then
                         if (dabs(internal0(3)) < threshold)then
                             passingThrough = 1d0
                         else
                             passingThrough = 0d0
                         end if
                     end if


                 endif


                 if (passingThrough>0d0  ) then






                     maxdist = MAXVAL(test_dist);
                     if (maxdist > Max_test_dist ) then
                         Max_test_dist = maxdist
                     endif

                     if (maxdist < err  ) then
                         counterCase = counterCase + 1

                         R_error = 0d0
                         cb1_error = 0d0
                         cb10_error = 0d0
                         cb2_error = 0d0
                         cb20_error = 0d0

                                     ! make temporary copies of the input data
                             ref1_temp0=ref1! reference vector for frag1 (original frame)
                             ref2_temp0=ref2! reference vector for frag2 (original frame)

                             ! make sure the original Cartesian frame of reference is at the CM of each fragment
                             if(natom1>1)call remove_center_of_mass(ref1_temp0,mass0(1:natom1),natom1,natom1)
                             if(natom2>1)call remove_center_of_mass(ref2_temp0,mass0(natom1+1:natom),natom2,natom2)

                             !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of
                             !!!    reference (origin at the new_CM of frag. 1)

                             ref1=ref1_temp0
                             ref2=ref2_temp0
                             ! (for CH3CN-He system: internal coordinates taken as spherical coords.)

                             call int_to_cart(internal, szi, mass, ref1_temp0, ref2_temp0, xdim, natom1, natom2, inputCoord, cart)

                             call center_of_mass_v2(cart(1:3*natom1),cmA,mass(1:natom1),natom1)
                             call center_of_mass_v2(cart(1:3*natom1),cmA0,mass0(1:natom1),natom1)
                             call center_of_mass_v2(cart(1+3*natom1:3*(natom1+natom2)),cmB,mass(1+natom1:natom1+natom2),natom2)
                             call center_of_mass_v2(cart(1+3*natom1:3*(natom1+natom2)),cmB0,mass0(1+natom1:natom1+natom2),natom2)


                             O  =   cmA
                             O0 =   cmA0


                             call cosine_law(O,O0,cmB,response1)




                             R0  = internal0(1)
                             R_ = internal(1)
                             cb10 = internal0(2)
                             cb1 = internal(2)


                             R0_test = response1(3)

                             cb1_test = response1(4)
                             cb10_test = response1(5)


                             R_error = dabs(R0-R0_test)
                             cb1_error = dabs(dabs(cb1)-dabs(cb1_test))
                             cb10_error = dabs(dabs(cb10)-dabs(cb10_test))



                             if (xdim > 3)then

                                     cb20 = internal0(3)
                                     cb2 = internal(3)
                                     B1 =   cart(1+3*natom1:3+3*natom1)

                                     if (NORM2(B1-cmB)<err)then
                                         B1 =   cart(4+3*natom1:6+3*natom1)
                                     end if

                                     call cosine_law(O0,cmB,B1,response2)
                                     cb20_test =-1d0*response2(5)
                                     cb20_error = dabs(dabs(cb20)-dabs(cb20_test))

                                     !write(*,*)"cb20 comparation: ",cb20,cb20_test,"error: ", cb20_error

                                     call cosine_law(O,cmB,B1,response3)
                                     cb2_test =-1d0*response3(5)
                                     cb2_error = dabs(dabs(cb2)-dabs(cb2_test))

                                     !write(*,*)"cb2 comparation: ",cb2,cb2_test,"error: ", cb2_error
                             end if

                         if (R_error > err .or. cb1_error > err .or. cb10_error > err &
                         .or. cb2_error > err .or. cb20_error > err) then
                                     write(file_output_number,*)"--------------------------------------- "
                                     write(file_output_number,*)"TestDistance :",test_dist
                                     write(file_output_number,*)"R comparation: ", R0,R0_test,"error: ", R_error
                                     write(file_output_number,*)"cb1 comparation: ", cb1,cb1_test,"error: ", cb1_error
                                     write(file_output_number,*)"cb10 comparation: ", cb10,cb10_test,"error: ", cb10_error
                                     write(file_output_number,*)"cb20 comparation: ",cb20,cb20_test,"error: ", cb20_error
                                     write(file_output_number,*)"cb2 comparation: ",cb2,cb2_test,"error: ", cb2_error
                                     write(file_output_number,*)"Error Tests R,cb1,cb10,cb2,cb20: ", &
                                     R_error,cb1_error,cb10_error,cb2_error,cb20_error
                                     write(file_output_number,*)"Internal : ", internal
                                     write(file_output_number,*)"Internal0 : ",internal0
                                     write(file_output_number,*)"----------------------------------------"
                                     test_failed = test_failed + 1

                         endif
                     endif
                 endif

                 !end Testing Section
         enddo
             write(file_output_number,*)
             write(file_output_number,*)"TEST NAME: Testing_AxisChanges"
             write(file_output_number,*)"System Features : Dimension/ Number of Atoms", xdim,"/",natom1,natom2
             write(file_output_number,*)"Int2Cart = "//inputCoord
             write(file_output_number,*)"Number of Tests", counterCase
             write(file_output_number,*)"Failed Tests: ", test_failed, " out of",test_number
             write(file_output_number,*)


     end subroutine axis_change_test


end module testing_auxiliar_functions
    
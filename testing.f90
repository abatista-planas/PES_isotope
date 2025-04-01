program comprenhensive_unit_test

    use mod_types, only : int32, real64
    use helper_functions, only : file_checking
    use testing_auxiliar_functions, only : test_cart_to_autosurf, interatomic_distance_test,performance

    implicit none
    integer (int32), parameter ::file_out_number=12,test_number = 100000
    integer (int32):: date_time (8),failed_test
    character (len = 10, kind = 1),dimension (3) ::big_ben


    ! Create a report in './Testing_Output.txt'
    call date_and_time (big_ben (1), big_ben (2), big_ben (3), date_time)

    call file_checking('./Testing_Output.txt', file_out_number)
    rewind(file_out_number)

    write(file_out_number, *)"******************************************************************************"
    write(file_out_number, *)  "Test Day and Time Record"
    write(file_out_number, *)  "Month / Day / Year: ", date_time(2), "/", date_time(3), "/", date_time(1)
    write(file_out_number, *)  "Hr    / Min / Sec : ", date_time(5), ":", date_time(6), ":", date_time(7)

    call test_cart_to_autosurf('./test/Cartesian/input_H2O.dat', "./test/Cartesian/2b.xyz", file_out_number)

    call interatomic_distance_test("./test/2DCase/input.dat", 2, "Autosurf", "Autosurf", failed_test, file_out_number)
    call interatomic_distance_test("./test/3DCase/input.dat", 3, "Autosurf", "Autosurf", failed_test, file_out_number)
    call interatomic_distance_test("./test/4DCase/input.dat", 4, "Autosurf", "Autosurf", failed_test, file_out_number)
    call interatomic_distance_test("./test/5DCase/input.dat", 5, "Autosurf", "Autosurf", failed_test, file_out_number)
    call interatomic_distance_test("./test/6DCase/input.dat", 6, "Autosurf", "Autosurf", failed_test, file_out_number)

    !Performance
    call performance(ntest=1000000)


    close(file_out_number)


end program comprenhensive_unit_test

    

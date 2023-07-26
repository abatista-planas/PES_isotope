module testingFunc
  
    contains

    SUBROUTINE Check_Cartesian_Frames(ptos,cart_model,cart,natom1,natom2,counterCase,test_failed,maxerr, &
                                        err_tolerance,XDIM,internal0,fileOutputNumber)
        use mathFunc
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     XDIM,ptos(3),natom1,natom2
        real*8,INTENT(IN) ::      err_tolerance,internal0(XDIM),cart_model((natom1+natom2)*3),cart((natom1+natom2)*3)
        real*8,INTENT(INOUT) ::   maxerr
        integer:: i,j,k
        real*8::err_test,modelA(6),modelB(6),A(6),B(6)
        integer,optional:: fileOutputNumber
        counterCase = counterCase + 1

        i =   ptos(1)
        j =   ptos(2)
        k =   ptos(3)

        Call CosineLaw(cart_model(1+3*(i-1):3*i),cart_model(1+3*(j-1):3*j),cart_model(1+3*(k-1):3*k),modelA)
        Call CosineLaw(cart_model(3*natom1+1+3*(i-1):3*natom1+3*i),cart_model(3*natom1+1+3*(j-1):3*natom1+3*j)&
                        ,cart_model(3*natom1+1+3*(k-1):3*natom1+3*k),modelB)

        Call CosineLaw(cart(1+3*(i-1):3*i),cart(1+3*(j-1):3*j),cart(1+3*(k-1):3*k),A)
        Call CosineLaw(cart(3*natom1+1+3*(i-1):3*natom1+3*i),cart(3*natom1+1+3*(j-1):3*natom1+3*j)&
                        ,cart(3*natom1+1+3*(k-1):3*natom1+3*k),B) 


        err_test = max(Norm2(A-modelA),Norm2(B-modelB))

                if (maxerr <  err_test) then 
                    maxerr =  err_test
                endif               
        
            
              if(err_test > err_tolerance) then
                 
                  test_failed = test_failed + 1

                 if (present(fileOutputNumber))then
                  write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(fileOutputNumber,*)"Cartesians Rigid Molec Checker            : ", counterCase
                  write(fileOutputNumber,*)"Error           : ", err_test
                  write(fileOutputNumber,*)"Norm2 Molec A: ",Norm2(A-modelA),"Norm2 Molec B: ",Norm2(B-modelB),&
                            "Norm2 AB : ",Norm2(A-B)
                            
                  write(fileOutputNumber,*)"Model Molec A: ",modelA
                  write(fileOutputNumber,*)"Molec A: ",A
                  write(fileOutputNumber,*)"Model Molec B: ",modelB
                  write(fileOutputNumber,*)"Molec B: ",B
                  write(fileOutputNumber,*)"Dimension :",XDIM 
                  write(fileOutputNumber,*) internal0
                  write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(fileOutputNumber,*)
                 
                 else
                  write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(*,*)"Cartesians Rigid Molec Checker           : ", counterCase
                  write(*,*)"Error           : ", err_test
                  write(*,*)"Norm2 Molec A: ",Norm2(A-modelA),"Norm2 Molec B: ",Norm2(B-modelB),&
                            "Norm2 AB : ",Norm2(A-B)
                            
                  write(*,*)"Model Molec A: ",modelA
                  write(*,*)"Molec A: ",A
                  write(*,*)"Model Molec B: ",modelB
                  write(*,*)"Molec B: ",B
                  write(*,*)"Dimension :",XDIM 
                  write(*,*) internal0
                  write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(*,*)
                 end if
                  
                            
              end if
              
         
    END SUBROUTINE
    !Reporting Error

    SUBROUTINE Report_Cartesian_Error(counterCase,test_failed,maxerr,td,err_tolerance,natoms,XDIM,cart,internal0&
                                      ,fileOutputNumber)
        
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     natoms,XDIM
        real*8,INTENT(IN) ::      td(4),err_tolerance,cart(natoms*3),internal0(XDIM)
        real*8,INTENT(INOUT) ::   maxerr
        integer,optional:: fileOutputNumber
        integer:: j
        real*8::err_test
          counterCase = counterCase + 1
          
          err_test = 0;
          err_test = sum(DABS(td))
          
          if (maxerr < err_test ) then 
            maxerr = err_test 
          endif
          
          if (err_test > err_tolerance)Then
            test_failed = test_failed + 1

            if (present(fileOutputNumber))then
              write(fileOutputNumber,*)
              write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(fileOutputNumber,*)"Cartesians-Autosurf  InterAtomic-Distance Test          : ", counterCase
              write(fileOutputNumber,*)"Error           : ", err_test,td
              do j=1,natoms
                  write(fileOutputNumber,*)"Atom ",j," :", cart(1+(j-1)*3:j*3)
                  write (fileOutputNumber,*)"----------------------------------------------"
              end do

              write(fileOutputNumber,*)"Dimension :",XDIM 
              write(fileOutputNumber,*) internal0
              write(fileOutputNumber,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(fileOutputNumber,*)

            else

              write(*,*)
              write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(*,*)"Cartesians-Autosurf  InterAtomic-Distance Test           : ", counterCase
              write(*,*)"Error           : ", err_test,td
              do j=1,natoms
                  write(*,*)"Atom ",j," :", cart(1+(j-1)*3:j*3)
                  write (*,*)"----------------------------------------------"
              end do

              write(*,*)"Dimension :",XDIM 
              write(*,*) internal0
              write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
              write(*,*)
            end if
          end if 
    END SUBROUTINE Report_Cartesian_Error

    
SUBROUTINE Test_Cart_to_Autosurf(systPath,dataPath,fileOutputNumber)
  use mathFunc
  use helperFunc

  IMPLICIT NONE
  character(*),INTENT(IN) :: systPath,dataPath
  integer,optional:: fileOutputNumber

  real*8::cart(18),cart_model(18),energies(4)
  character(len=1)::Atom_label
  real*8 :: internal0(6),td(4),err_test,err_tolerance,maxerr_1,maxerr_2
  integer :: stat2,stat3,i,natoms,int_to_cart_func,XDIM,ptos(3)
  integer :: failedTest_1,counterCase_1,failedTest_2,counterCase_2
  real*8::response1_model(6),response2_model(6),response1(6),response2(6)
  integer :: stat


  XDIM=6
  
  
  open(unit=100,file=dataPath,status='old',action='read')
    ! OPEN(unit=200, FILE="./output/Paesani/coord_H2O_Dymer.txt",ACTION='write',IOSTAT=stat2,STATUS='OLD')
    ! if (stat2 .ne. 0) then
    !         open(unit=200,file="./output/Paesani/coord_H2O_Dymer.txt",status='new',action='write')
    ! end if

    ! OPEN(unit=300, FILE="./output/Paesani/coord_H2O_Dymer_filtered.txt",ACTION='write',IOSTAT=stat3,STATUS='OLD')
    ! if (stat3 .ne. 0) then
    !         open(unit=300,file="./output/Paesani/coord_H2O_Dymer_filtered.txt",status='new',action='write')
    ! end if

    failedTest_1 = 0
    counterCase_1 = 0
    failedTest_2 = 0
    counterCase_2 = 0
    maxerr_1 = -1d0
    maxerr_2 = -1d0
    err_tolerance=1d-7
    int_to_cart_func = -1
    ptos = (/1,2,3 /)
    i=0

       do !i=1,2510!42508  !42508
         
         i=i+1
         
  
       
          read(100,*,iostat=stat)natoms
          if (stat /= 0 .or. i>2512) exit

          read(100,*)energies(1:4)
          read(100,*)Atom_label,cart(1:3)
          read(100,*)Atom_label,cart(4:6)
          read(100,*)Atom_label,cart(7:9)
          read(100,*)Atom_label,cart(10:12)
          read(100,*)Atom_label,cart(13:15)
          read(100,*)Atom_label,cart(16:18)

           !-1 is for cartesian input

         if (i==1)then
              cart_model=cart
         end if

   
         Call Get_ISOTOP_COORDINATES(cart,size(cart),internal0,6, int_to_cart_func ,systPath,testArr_Errors = td)



!Reporting Error

          
          Call Report_Cartesian_Error(counterCase_1,failedTest_1,maxerr_1,td,err_tolerance,natoms,XDIM,cart,internal0,&
          fileOutputNumber)
          call Check_Cartesian_Frames(ptos,cart_model,cart,3,3,counterCase_2,failedTest_2,maxerr_2&
                                     ,err_tolerance,XDIM,internal0,fileOutputNumber)
!End Error Testing 
          
        !   write(200, *) i , internal0, energies(2)
        !   if (internal0(1)>5) then 
        !    write(300, *) i , internal0, energies(2)
          
        !   endif


   
      enddo
      
        if (present(fileOutputNumber))then
                write(fileOutputNumber,*)" "
                write(fileOutputNumber,*) "%%%%%%%%%%%%%%%    Cartesian Test    %%%%%%%%%%%%%%%"
                write(fileOutputNumber,*) "Number of failed Cartesian coordinates: ",failedTest_1, " out of ",&
                                         counterCase_1, ", maxErr: ",maxerr_1,' ' 
                write(fileOutputNumber,*) "Number of failed Check_Cartesian_Frames: ",failedTest_2, " out of ",&
                                         counterCase_2,", maxErr: ",maxerr_2 ,' '
                write(fileOutputNumber,*)" "
        endif




    call Print_Test_Results("Cartesians-Autosurf  InterAtomic-Distance Test",counterCase_1,counterCase_1 - failedTest_1)
    call Print_Test_Results("Cartesians Rigid Molec Checker                ",counterCase_2,counterCase_2 - failedTest_2)

      
   close(100)
  !  close(200)
  !  close(300)

End SUBROUTINE Test_Cart_to_Autosurf

  
SUBROUTINE Print_Test_Results(testName,ntest,success_test)

    IMPLICIT NONE
    character(*),INTENT(IN) :: testName
    integer,INTENT(IN)::ntest,success_test
    character (len =25 ) ::f_res1,f_res2 ,str
    character (len =25) :: str1
    integer::failed,ptg




    f_res1 = achar(27)//'[32m TEST PASSED ! '//achar(27)//'[0m'
    f_res2 = achar(27)//'[31m TEST FAILED :('//achar(27)//'[0m'

    ptg = FLOOR(success_test*100d0/ntest)
    failed = ntest-success_test
    write (str1, *) ptg
    !str = adjustl(str)
     str = trim(str1)

    if (ntest > success_test)then
 
      write(*, '(A, I5, A, I1, A, I5, A)') f_res2//" ----  "//testName//"            success: "&
                                       ,success_test,"  and failed : ",failed, " -----> out of " ,ntest &
                                       ,"("//achar(27)//'[31m'//str//'% '//achar(27)//'[0m)'
    else
      write(*,  '(A, I5, A, I1, A, I5, A)') f_res1//" ----  "//testName//"            success: "&
                                       ,success_test," / failed : ",failed, " -----> out of " ,ntest,&
                                       "("//achar(27)//'[32m          100            % '//achar(27)//'[0m)'
    end if



end SUBROUTINE Print_Test_Results










end module testingFunc
    
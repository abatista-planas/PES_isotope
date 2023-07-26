module testingFunc
  
    contains

    SUBROUTINE Check_Cartesian_Frames(ptos,cart_model,cart,natom1,natom2,counterCase,test_failed,maxerr, &
                                        err_tolerance,XDIM,internal0)
        use mathFunc
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     XDIM,ptos(3),natom1,natom2
        real*8,INTENT(IN) ::      err_tolerance,internal0(XDIM),cart_model((natom1+natom2)*3),cart((natom1+natom2)*3)
        real*8,INTENT(INOUT) ::   maxerr
        integer:: i,j,k
        real*8::err_test,modelA(6),modelB(6),A(6),B(6)

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

                 
                  write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
                  write(*,*)"Cartesian coordinates           : ", counterCase
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
              
         
    END SUBROUTINE
    !Reporting Error

    SUBROUTINE Report_Cartesian_Error(counterCase,test_failed,maxerr,td,err_tolerance,natoms,XDIM,cart,internal0)
        
        IMPLICIT NONE
        integer,INTENT(INOUT) ::  counterCase,test_failed
        integer,INTENT(IN) ::     natoms,XDIM
        real*8,INTENT(IN) ::      td(4),err_tolerance,cart(natoms*3),internal0(XDIM)
        real*8,INTENT(INOUT) ::   maxerr
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
            write(*,*)
            write(*,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%",'\n'
            write(*,*)"Cartesian coordinates           : ", counterCase
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
    END SUBROUTINE Report_Cartesian_Error

    












end module testingFunc
    
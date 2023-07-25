module testingFunc
  
    contains

    SUBROUTINE MolecularDistanceTest(arr1,ref1,ref2,N1,N2,internal0,XDIM,Distance)
        use mathFunc
        use coordinateTransf
        IMPLICIT NONE
        integer,INTENT(IN) ::  N1,N2,XDIM
        real*8,INTENT(IN) ::  internal0(XDIM),arr1(3,N1+N2),ref1(3,N1),ref2(3,N2)
        real*8,INTENT(INOut) ::  Distance(2)
        real*8 ::  d1((N1+N2)*(N1+N2-1)/2),d0((N1+N2)*(N1+N2-1)/2),ref1_temp(3,N1),ref2_temp(3,N2)
        real*8:: R,rotarr1(3,N1),rotarr2(3,N2),arr0(3,N1+N2),pii,angles(6)
        real*8:: a1,b1,c1,a2,b2,c2,internal0_temp(XDIM)
        
  
        pii=dacos(-1d0) 
        R=internal0(1)
  
        internal0_temp = internal0
        ref1_temp=ref1
        ref2_temp=ref2
  
  
        call Internal0_to_anglesZYZ(internal0_temp,XDIM,angles)
  
  
            
            a1 = angles(1)
            b1 = angles(2)
            c1 = angles(3)
            a2 = angles(4)
            b2 = angles(5)
            c2 = angles(6)
  
            ! write(*,*)"a1,b1,c1 ",a1,b1,c1
            ! write(*,*)"a2,b2,c2 ",a2,b2,c2
  
  
            call MolecularRotation_ZYZ(ref1_temp,N1,a1,b1,c1,rotarr1)
            call MolecularRotation_ZYZ(ref2_temp,N2,a2,b2,c2,rotarr2)
        
        
            arr0(1:3,1:N1) = rotarr1
            arr0(1:2,N1+1:N2) = rotarr2(1:2,1:N2)
            arr0(3,N1+1:N2) = rotarr2(3,1:N2)+R
  
  
        call InterAtomicDistance(arr1,(N1+N2),d1)
        call InterAtomicDistance(arr0,(N1+N2),d0)
  
    
  
        call ComparedDistance(d1,d0,(N1+N2)*((N1+N2)-1)/2,Distance)
        
    
    END SUBROUTINE MolecularDistanceTest


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
    
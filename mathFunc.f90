module mathFunc
  
    contains
  
    ! given the the three catersian points of a triangule
    ! this function returns all interior angles and the side lengths 
    subroutine CosineLaw(A,B,C,response)
        IMPLICIT NONE
        real*8,INTENT(IN) ::A(3),B(3),C(3)
        real*8,INTENT(INOUT) ::response(6)
        real*8 ::AB,BC,AC,cos_a,cos_b,cos_c
    
        !Respose will return the triangule measurements AB,AC,BC, and the angles CAB, ABC, BCA
    
        AB = NORM2(B-A)
        AC = NORM2(C-A)
        BC = NORM2(B-C)
    
    
        response(1) = AB
        response(2) = AC
        response(3) = BC
    
    
        cos_a = (AB**2 + AC**2 - BC**2)/(2d0*AC*AB)
        cos_b = (AB**2 + BC**2 - AC**2)/(2d0*BC*AB)
        cos_c = (AC**2 + BC**2 - AB**2)/(2d0*BC*AC)
    
        response(4) = cos_a
        response(5) = cos_b
        response(6) = cos_c
    
        return
    end subroutine CosineLaw
  
    SUBROUTINE MolecularDistanceTest(arr1,ref1,ref2,N1,N2,internal0,XDIM,Distance)
        
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
  
    SUBROUTINE MolecularRotation_ZYZ(arr,N,alpha,beta,gamma,rotarr)
        IMPLICIT NONE
        integer,INTENT(IN) ::  N
        
        real*8,INTENT(IN) ::  arr(3,N)
        real*8,INTENT(INOut) ::  rotarr(3,N)
        real*8,INTENT(IN):: alpha,gamma,beta
    
        real*8::a,b,g,ca,cb,cg,sa,sb,sg,U_rot(3,3)
  
    
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
  
        
        U_rot=0d0
        U_rot(1,1)=ca*cb*cg - sa*sg
        U_rot(2,2)=ca*cg - cb*sa*sg
        U_rot(3,3)=cb
        U_rot(1,2)=-cg*sa-ca*cb*sg
        U_rot(1,3)=ca*sb
        U_rot(2,1)=ca*sg+cb*cg*sa
        U_rot(2,3)=sa*sb
        U_rot(3,1)=-cg*sb
        U_rot(3,2)=sb*sg
  
        
        call rotmol2 (N, arr, rotarr, U_rot)
        
        
        
    END SUBROUTINE MolecularRotation_ZYZ
  
    SUBROUTINE InterAtomicDistance(arr,N,Distance)
      IMPLICIT NONE
      integer,INTENT(IN) ::  N
      Integer::count,i,j
      real*8,INTENT(IN) ::  arr(3,N)
      real*8,INTENT(INOut) ::  Distance(N*(N-1)/2)
      real*8 :: pii,xij,yij,zij
      
      pii=dacos(-1d0) 
  
      count=0;
  
  
  
  
      do i=1,N-1
          do j=i+1,N
  
              count = count + 1;
              xij = (arr(1,i)-arr(1,j))**2
              yij = (arr(2,i)-arr(2,j))**2
              zij = (arr(3,i)-arr(3,j))**2
              Distance(count) = DSQRT(xij + yij + zij)
          
  
          enddo
      enddo
  
  
  
    END SUBROUTINE InterAtomicDistance
  
    integer*2 function cmp_function( a, b )
        INTEGER*4 a, b
        if ( a .lt. b ) compar = -1
        if ( a .eq. b ) compar = 0
        if ( a .gt. b ) compar = 1
        return
    end
  
    SUBROUTINE ComparedDistance(Distance1,Distance2,N,result)
        IMPLICIT NONE
    
        integer,INTENT(IN) ::  N
        real*8,INTENT(IN) ::  Distance1(N),Distance2(N)
        real*8 ::  D1(N),D2(N)
        real*8,INTENT(OUT)::result(2)
        real*8::sr1,sr2
        integer :: i
  
        
        result=0d0
  
        D1 = Distance1
        D2 = Distance2
        call sort_pick(D1,N) 
        call sort_pick(D2,N) 
        
        do i = 1,N
            result(1) = result(1) + (D1(i)-D2(i))**2
        enddo
        result(1) = DSQRT(result(1))
        !Global difference
  
  
        sr1 = 0d0
        sr2 = 0d0
        do i = 1,N
        
            sr1  = sr1 + D1(i)
            sr2  = sr2 + D2(i)
        enddo
  
        result(2)=DABS(sr1-sr2)
  
  
    END SUBROUTINE ComparedDistance
  
    SUBROUTINE sort_pick(arr,N) 
        
        IMPLICIT NONE 
        INTEGER, INTENT(IN) :: N
        REAL*8, DIMENSION(:), INTENT(INOUT) :: arr(N) !Sorts an array arr into ascending numerical order, by straight insertion. arr is replaced on output by its sorted rearrangement. 
        INTEGER :: i,j 
        REAL*8 :: a 
        
  
        
        do j=2,N 
            a=arr(j) 
            do i=j-1,1,-1 !Pick out each element in turn. Look for the place to insert it. 
                if (arr(i) <= a) exit 
                arr(i+1)=arr(i) 
            end do 
            arr(i+1)=a !Insert it. 
        end do 
    END SUBROUTINE sort_pick
  
    SUBROUTINE Print_Vector(vec,N,str)
        IMPLICIT NONE
        integer,INTENT(IN) ::  N
        real*8,INTENT(IN) ::  vec(3,N)
        character(len = 20) ,INTENT(IN) ::  str
        integer :: i
  
        write(*,*) "******* Vector Printing: "," ************"
        write(*,*) "******* ",str," ************"
        do i = 1,N
            write(*,*)vec(:,i)
        enddo
  
        write(*,*) "******* End Printing ************"
  
    END SUBROUTINE Print_Vector
  
    subroutine Find_EulerAngles(natom,cart_ref1,cart_mat1,masses,alpha,beta,gamma)
  
        IMPLICIT NONE
  
        integer,INTENT(IN) :: natom
        real*8,INTENT(OUT):: alpha,beta,gamma
        real*8,INTENT(IN):: masses(natom),cart_ref1(3,natom),cart_mat1(3,natom)
        integer :: ierr
        real*8 :: U_rot(3,3),quat(4),threshold,U33
      !!! ----------------------------------------------------------------------------------
      !!! 4) once the geometry is in the proper frame, find the corresponding Euler angles:
      !!! ----------------------------------------------------------------------------------
  
      call qtrfit(natom,cart_ref1,cart_mat1,masses,quat,U_rot,ierr)
  
      U33 = U_rot(3,3)
  
      if(U33>1d0)then
        U33=1d0
      elseif(U33<-1d0)then
        U33=-1d0
      endif
  
      ! solve for Euler angles
      threshold = 1d0 - 1d-16!0.9993d0

  
      
      beta=dacos(U33)
      
      if(U33>threshold)then
        gamma =datan2(-U_rot(1,2),U_rot(1,1))
        alpha =0d0
      elseif( U33 <-1d0*threshold)then
        gamma=datan2(U_rot(1,2),-U_rot(1,1))
        alpha=0d0
      else
        alpha=datan2(U_rot(2,3),U_rot(1,3))
        gamma=datan2(U_rot(3,2),-U_rot(3,1))
      endif
  
  
  
      return
    end subroutine Find_EulerAngles
  
    subroutine cmass(cart,cm,mass,natom,natom2)
      IMPLICIT NONE
      integer :: k,kp,natom,natom2
      real*8 :: mass(natom),cart(natom*3),mtot,cm(3)
      mtot=0d0
      do k=natom-natom2+1,natom
        mtot=mtot+mass(k)
      enddo
      cm=0d0
      do k=natom-natom2+1,natom
        do kp=1,3
            cm(kp)=cm(kp)+cart((k-1)*3+kp)*mass(k)
        enddo
      enddo
      cm=cm/mtot
      return
    end subroutine cmass

    subroutine cmass2(cart,cm,mass,natoms)
      IMPLICIT NONE
      integer :: k,kp,natoms
      real*8 :: mass(natoms),cart(natoms*3),mtot,cm(3)
      mtot=0d0
      do k=1,natoms
        mtot=mtot+mass(k)
      enddo
      cm=0d0
      do k=1,natoms
        do kp=1,3
            cm(kp)=cm(kp)+cart((k-1)*3+kp)*mass(k)
        enddo
      enddo
      cm=cm/mtot
      return
    end subroutine cmass2

    subroutine rm_cmass(cart,mass,natom,natom1)
      IMPLICIT NONE
      integer :: k,kp,natom,natom1
      real*8 :: mass(natom),cart(natom*3),mtot,cmass1(3)
      mtot=0d0
      do k=1,natom1
        mtot=mtot+mass(k)
      enddo
      cmass1=0d0
      do k=1,natom1
        do kp=1,3
            cmass1(kp)=cmass1(kp)+cart((k-1)*3+kp)*mass(k)
        enddo
      enddo
      cmass1=cmass1/mtot
  
      do k=1,natom
        do kp=1,3
            cart((k-1)*3+kp)=cart((k-1)*3+kp)-cmass1(kp)      
        enddo
      enddo
      return
    end subroutine rm_cmass

    subroutine vec_to_mat2(cart_perms,cart_mat,natom)
      IMPLICIT NONE
      integer :: k,kp,natom
      real*8 :: cart_perms(3*natom),cart_mat(3,natom)
      do k=1,natom
        do kp=1,3
            cart_mat(kp,k)=cart_perms((k-1)*3+kp)
        enddo
      enddo
      return
    end subroutine vec_to_mat2
  
    subroutine mat_to_vec2(cart_mat,cart_perms,natom)
      IMPLICIT NONE
      integer :: k,kp,natom
      real*8 :: cart_perms(3*natom),cart_mat(3,natom)
      do k=1,natom
        do kp=1,3
            cart_perms((k-1)*3+kp)=cart_mat(kp,k)
        enddo
      enddo
      return
    end subroutine mat_to_vec2
  
  
    SUBROUTINE rotmol2 (n, x, molrot, u)
      IMPLICIT None
  
      INTEGER,INTENT(IN) ::n
      real*8,INTENT(IN):: x(3, n)
      real*8,INTENT(Out):: molrot(3, n)
      real*8,INTENT(IN):: u(3, 3)
      real*8:: mol(3, n),urot(3, 3)
      INTEGER:: i
  
      urot = u
      mol  = x
  
      DO i = 1, n
        molrot(1,i) = urot(1, 1) * mol(1, i) + urot(1, 2) * mol(2, i) + urot(1, 3) * mol(3, i)
        molrot(2,i) = urot(2, 1) * mol(1, i) + urot(2, 2) * mol(2, i) + urot(2, 3) * mol(3, i)
        molrot(3,i) = urot(3, 1) * mol(1, i) + urot(3, 2) * mol(2, i) + urot(3, 3) * mol(3, i)
  
      enddo
  
      RETURN
    end subroutine  rotmol2

  
end module mathFunc
    
module helperFunc
  
    contains

    SUBROUTINE Read_File(filename,mass,mass0,natom1,natom2,ref1,ref2,XDIM,i0Type)
      IMPLICIT NONE
      integer,INTENT(INOUT) :: natom1,natom2,XDIM
      real*8, allocatable,INTENT(INOUT) :: ref1(:),ref2(:),mass0(:),mass(:)
      integer :: natom,i,k
      character(len=25) :: i0Type ! it defines Internal0 coodinate system in the output [options: "BiSpherical","Autosurf"]
    
     
      character(*),INTENT(IN)::filename
    
      open(unit=10,file=filename,status='old',action='read')
    
      read(10,*) XDIM
      read(10,*) natom1! number of atoms in fragment1
      read(10,*) natom2! number of atoms in fragment2

      natom=natom1+natom2 ! total number of atoms
    
      read(10,*) i0Type! Internal0 coodinate system in the output [options: "BiSpherical","Autosurf"]

      
    
      allocate(ref1(3*natom1),ref2(3*natom2),mass0(natom),mass(natom))
    
    
    
      if(natom1>1)then
      do i=1,natom1
          read(10,*)mass(i),mass0(i),(ref1(k),k=i*3-2,i*3)
      enddo
      elseif(natom1==1)then
      read(10,*)mass(1),mass0(1)
      ref1=0.d0
      endif
    
      if(natom2>1)then
      do i=1,natom2
          read(10,*)mass(natom1+i),mass0(natom1+i),(ref2(k),k=i*3-2,i*3)
      enddo
      elseif(natom2==1)then
      read(10,*)mass(natom1+1),mass0(natom1+1)
      ref2=0.d0
      endif
      close(10)
    END SUBROUTINE  Read_File

  
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

    SUBROUTINE convert_isotopic_coordinates(inter,inter0,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM &
                                            ,iFun,t_dist,doTest,inter0_Sys)
       
        IMPLICIT NONE
  
        integer,INTENT(IN) :: natom1,natom2,XDIM
        real*8,INTENT(IN) ::  inter(XDIM),mass(natom1+natom2),mass0(natom1+natom2),ref1_0(natom1*3),ref2_0(natom2*3)
        real*8,INTENT(OUT) :: inter0(XDIM)
        integer ::natom,intFunc
        real*8 :: ref1_temp0(natom1*3),ref2_temp0(natom2*3),cm(3)
        real*8 :: ref1(natom1*3),ref2(natom1*3),cart_ref1(3,natom1),cart_mat1(3,natom1)
        real*8 :: cart_ref2(3,natom2),cart_mat2(3,natom2),cart_frag2(natom2*3)
        real*8 :: cart((natom1+natom2)*3),cart_mat(3,natom1+natom2)
        real*8 :: U_rot(3,3),inter0_test(XDIM)
        real*8 :: gamma1,gamma2,beta1,beta2,alpha1,alpha2,theta,phi,pii,cos_theta,sin_theta
        real*8 :: test_cart_i(3,natom1+natom2),test_cart_f(3,natom1+natom2),distance_arr_test(2),td(5)
        real*8,INTENT(OUT) :: t_dist(5)
        integer,INTENT(IN) :: iFun,doTest ! doTest = -1 will not do the distance test, any other value will 
        character(*) :: inter0_Sys ! it defines what EulerAngle_to_FinalInternal0 has to be defined 
  
       pii=dacos(-1d0) 
       natom=natom1+natom2
       td = 0
       intFunc = iFun
   
       ! make temporary copies of the input data
       ref1_temp0=ref1_0! reference vector for frag1 (original frame)
       ref2_temp0=ref2_0! reference vector for frag2 (original frame)
   
       ! make sure the original Cartesian frame of reference is at the CM of each fragment 
       if(natom1>1)call rm_cmass(ref1_temp0,mass0(1:natom1),natom1,natom1)
       if(natom2>1)call rm_cmass(ref2_temp0,mass0(natom1+1:natom),natom2,natom2)  
   
       !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
       !!!    reference (origin at the new_CM of frag. 1)
   
       ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
        ref1=ref1_temp0
        ref2=ref2_temp0
       ! (for CH3CN-He system: internal coordinates taken as spherical coords.)
   
  
        call Int_to_Cart(inter,mass,ref1,ref2,XDIM,natom1,natom2,intFunc,cart)
  
  
  
       call vec_to_mat2(cart,test_cart_i,natom) ! this is for testing purposes
       !write(*,*)"test_cart_i ",test_cart_i(:,7)
       !call Print_Vector(test_cart_i,7,"test_cart_i")
  
  
          !call dcmass(cart,mass,mass0,natom1,natom2)
  
         !!! 2) transfer the origin of the reference frame to the original CM of frag1 
          call rm_cmass(cart,mass0,natom1+natom2,natom1)
   
          !write(*,*)"cart",cart(3*natom1+1:3*(natom1+natom2))
  
          !!! 3) align Z and Z' frames
          ! find new relative position of the original CM of frag2 (with respect to the original CM of frag1)
  
          call cmass(cart,cm,mass0,natom1+natom2,natom2)
    
          !write(*,*)"cm",cm
  
  
          ! find spherical coordinates of the c.m. of frag2
          inter0(1) = dsqrt(cm(1)**2+cm(2)**2+cm(3)**2)! "R"! Norm2
  
   
  
          theta = dacos(cm(3)/inter0(1))! "TH"
          cos_theta = cm(3)/inter0(1)
          sin_theta = DSQRT(1d0 - cos_theta**2)
  
          if(dabs(dcos(theta))>1d0-1d-12)then! "PHI"
            td(5) = -1
            phi=0d0
          else
            phi=datan2(cm(2),cm(1))
          endif
  
  
          !write(*,*)"Urot : ",internal0(1) ,theta,phi
  
       !write(*,*)"internal0",internal,internal0
  
        ! construct rotation matrix
  
       U_rot=0d0
       U_rot(1,1)=dcos(phi)*cos_theta
       U_rot(2,2)=dcos(phi)
       U_rot(3,3)=cos_theta
       U_rot(1,2)=dsin(phi)*cos_theta
       U_rot(1,3)=-1d0*sin_theta
       U_rot(2,1)=-1d0*dsin(phi)
       U_rot(3,1)=dcos(phi)*sin_theta
       U_rot(3,2)=dsin(phi)*sin_theta
  
  
  
       ! transform vector "cart[1:3*natom]" into a matrix: "cart_mat[1:3,1:natom]"
       call vec_to_mat2(cart,cart_mat,natom1+natom2)
       ! rotate the initial geometry to align (make "Z" axis contain the c.m. of both frags!)
       call rotmol2(natom1+natom2,cart_mat,cart_mat,U_rot)! Now "cart_mat" is in the proper Frame
  
       !write(*,*)"cart_mat",cart_mat(:,7)
  
       ! transform back matrix "cart_mat[1:3,1:natom]" into a vector: "cart[1:3*natom]"
       call mat_to_vec2(cart_mat,cart,natom1+natom2)! Now "cart" is in the proper/original Frame: where the 
       ! PES was initially fitted, with its origin at the CM of frag1 and the Z-axis containing both CMs
  
       call vec_to_mat2(cart,test_cart_f,natom) ! testing 
       !call Print_Vector(test_cart_f,natom,"test_cart_f")
  
  
       !!! ----------------------------------------------------------------------------------
       !!! 4) once the geometry is in the proper frame, find the corresponding Euler angles:
       !!! ----------------------------------------------------------------------------------
  
       ! Fragment 1 
       ! -----------
       ! cart_mat1 = new (rotated/aligned) coordinates of frag 1 
       cart_mat1(1:3,1:natom1)=cart_mat(1:3,1:natom1)
       ! transform vector "ref1_temp0[1:3*natom1]" (original reference for frag1) into a matrix: "cart_ref1[1:3,1:natom1]"
       call vec_to_mat2(ref1_temp0,cart_ref1,natom1)
       ! find rotation matrix (U_rot) that transforms (rotates) coordinates "cart_ref2" into "cart_mat2"
  
  
  
       call Find_EulerAngles(natom1,cart_ref1,cart_mat1,mass0(1:natom1),alpha1, beta1,gamma1)
   
  
  
       ! Fragment 2 
       ! -----------
       cart_frag2=cart(3*natom1+1:3*(natom1+natom2))
       ! set center of mass of frag2 at origin
       call rm_cmass(cart_frag2,mass0(natom1+1:natom1+natom2),natom2,natom2)
       ! transform vector "cart_frag2[1:3*natom2]" (new --rotated/aligned-- coord. of frag 1) into a matrix: "cart_mat2[1:3,1:natom2]"
       call vec_to_mat2(cart_frag2,cart_mat2,natom2)
       ! transform vector "ref2_temp0[1:3*natom2]" (original reference for frag2) into a matrix: "cart_ref2[1:3,1:natom2]"
       call vec_to_mat2(ref2_temp0,cart_ref2,natom2)
       ! find rotation matrix (U_rot2) that transforms (rotates) coordinates "cart_ref2" into "cart_mat2"
  
       call Find_EulerAngles(natom2,cart_ref2,cart_mat2,mass0(natom1+1:natom1+natom2),alpha2, beta2,gamma2)
  
  
       ! *** internal coordinates ***
       ! ----------------------------
  
       
  
     
  
       if (inter0_Sys=="Autosurf")Then
        call EulerAngles_2_Autosurf(XDIM,inter0(1),alpha1,beta1,gamma1,alpha2,beta2,gamma2,inter0)
       elseif (inter0_Sys=="BiSpherical")Then
        call EulerAngles_2_BiSpherical(inter,inter0(1),alpha1,beta1,gamma1,alpha2,beta2,gamma2,inter0)
       end if
  
  
       !write(*,*)"internal0 : ",internal0
  
      !!! Testing Unit
  
       ! This function test if the cart coodinates given by Int_to_Cart in step 1
       ! has the same interatomic distances than the bimolecular system in the internal coordiantes internal0
       ! This test calculate an array d1 and d0 of all interatomic distance of all atoms taken in pairs for both 
       ! cartesian configurations   and print the sum of abs(d1-d0)
   
       if(doTest>0)then
              call EulerAngles_2_Autosurf(XDIM,inter0(1),alpha1,beta1,gamma1,alpha2,beta2,gamma2,inter0_test)
              call  MolecularDistanceTest(test_cart_f,cart_ref1,cart_ref2,natom1,natom2,inter0_test,XDIM,distance_arr_test)
              td(1:2) = distance_arr_test
              call  MolecularDistanceTest(test_cart_i,cart_ref1,cart_ref2,natom1,natom2,inter0_test,XDIM,distance_arr_test)
              td(3:4) = distance_arr_test
  
              t_dist=td
        endif
  
  
       !!! End Testing Unit
  
    END SUBROUTINE convert_isotopic_coordinates

    SUBROUTINE Get_ISOTOP_COORDINATES(internal,internal0,XDIM,ifun,filePath)
     
      
        IMPLICIT NONE
        integer,INTENT(IN) :: XDIM,ifun
        character(*),INTENT(IN) :: filePath
        real*8, INTENT(IN):: internal(XDIM)
        real*8, INTENT(OUT):: internal0(XDIM)
        integer :: natom1,natom2
        real*8, allocatable :: ref1_0(:),ref2_0(:),mass(:),mass0(:)
        real*8 ::  int0_AS(XDIM),int0_BiSph(XDIM),td(5)
        Integer :: ifun_temp,Xdim_file
        integer :: initflag
        save initflag
        character(len=25) :: i0Type ! it defines Internal0 coodinate system in the output
        data initflag /1/
        save mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file,i0Type
        
        ! i0Type options: "BiSpherical","Autosurf","Cartesian"
        
      
        IF(initflag==1)THEN! initialize 
         Call Read_File(filePath,mass,mass0,natom1,natom2,ref1_0,ref2_0,Xdim_file,i0Type)
         initflag=2  
        ENDIF
      
        i0Type = trim(i0Type)
        ifun_temp =ifun
      
      
       Call convert_isotopic_coordinates(internal,internal0,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM,&
                                          ifun_temp,td,-1,i0Type)
      
      
      
    END SUBROUTINE Get_ISOTOP_COORDINATES
  
end module helperFunc
    
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

  

    SUBROUTINE convert_isotopic_coordinates(inter,inter0,mass,mass0,natom1,natom2,ref1_0,ref2_0,XDIM &
                                            ,iFun,t_dist,doTest,inter0_Sys)
       
        use mathFunc
        use coordinateTransf
        use testingFunc , only: MolecularDistanceTest

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
        character(*),INTENT(IN) :: inter0_Sys ! it defines what EulerAngle_to_FinalInternal0 has to be defined 
  
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
    
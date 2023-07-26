module coordinateTransf
    
    contains


    ! convert: Internal coordinates in ZYZ--> BiSpherical for 5D ONLY
    SUBROUTINE ZYZ_to_BiSpherical(internal,internal0_ZYZ,internal0_BiSph)

      !************** Subroutine definition 
          ! subroutine ZYZ --> BiSpherical (convert: Internal coordinates in ZYZ--> BiSpherical for 5D)  
          ! internal0_ZYZ(1) = R
          ! internal0_ZYZ(2) = cos_b1
          ! internal0_ZYZ(3) = cos_b2
          ! internal0_ZYZ(4) = alpha
          ! internal0_ZYZ(5) = gamma

          ! internal0_BiSph(1) = R
          ! internal0_BiSph(2) = theta1
          ! internal0_BiSph(3) = theta2
          ! internal0_BiSph(4) = phi1
          ! internal0_BiSph(5) = phi2
      !********************************************************
          IMPLICIT NONE
          real*8,INTENT(IN) :: internal0_ZYZ(5) ,internal(5)     
          real*8,INTENT(OUT)::internal0_BiSph(5) 
          real*8:: pii,term1,term2,threshold,internal0p_ZYZ(5)  
          
          !write(*,*)"internal0_ZYZ",internal0_ZYZ
          pii=dacos(-1d0)
          threshold = 1d0 - 1d-12!0.9993d0
          internal0p_ZYZ = internal0_ZYZ

          internal0_BiSph(1) =  internal0_ZYZ(1) ! R

          if(internal0p_ZYZ(2)>1d0)then
            internal0p_ZYZ(2) = 1d0
          elseif(internal0p_ZYZ(2)<-1d0)then
            internal0p_ZYZ(2) = -1d0
          endif

          if(internal0p_ZYZ(3)>1d0)then
            internal0p_ZYZ(3) = 1d0
          elseif(internal0p_ZYZ(3)<-1d0)then
            internal0p_ZYZ(3) = -1d0
          endif


     ! Note : There is a discontinuity in Phi =0 for cos_b1 =0 
       ! Taylor serie for cos b1 
       ! Analyse b1 = 0 and b2 = 0


            
  
          if(DABS(internal0p_ZYZ(2))>threshold .and. DABS(internal0p_ZYZ(3))>threshold)then

              if(internal0p_ZYZ(2)>threshold)then
                internal0_BiSph(2) = 0d0
                internal0_BiSph(4) = internal(4)
              elseif( internal0p_ZYZ(2) <-1d0*threshold)then
                internal0_BiSph(2) = pii
                internal0_BiSph(4) = internal(4)
              else
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(4) = pii- internal0p_ZYZ(5)
              end if
              if(internal0p_ZYZ(3)>threshold)then
                internal0_BiSph(3) = 0d0
                internal0_BiSph(5) = internal(5)
              elseif( internal0p_ZYZ(3) <-1d0*threshold)then
                internal0_BiSph(3) = pii
                internal0_BiSph(5) = internal(5)
              end if

          else

              if(internal0p_ZYZ(2)>threshold)then
                internal0_BiSph(2) = 0d0
                internal0_BiSph(3) = DACOS(internal0p_ZYZ(3))
                internal0_BiSph(4) = internal(4)
                internal0_BiSph(5) = -1d0*(internal0p_ZYZ(4)+internal0p_ZYZ(5))
              elseif( internal0p_ZYZ(2) <-1d0*threshold)then
                internal0_BiSph(2) = pii
                internal0_BiSph(3) = DACOS(-1d0*internal0p_ZYZ(3))
                internal0_BiSph(4) = internal(4)
                internal0_BiSph(5) = internal0p_ZYZ(4)+internal0p_ZYZ(5)-pii
              elseif( internal0p_ZYZ(3) >threshold)then
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(3) = 0d0
                internal0_BiSph(4) =  pii- internal0p_ZYZ(5)
                internal0_BiSph(5) = internal(5)
              elseif( internal0p_ZYZ(3) <-1d0*threshold)then
                internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                internal0_BiSph(3) = pii
                internal0_BiSph(4) = pii- internal0p_ZYZ(5)
                internal0_BiSph(5) = internal(5) 
              else

                  internal0_BiSph(2) = DACOS(internal0p_ZYZ(2))
                  internal0_BiSph(3) = DACOS(internal0p_ZYZ(2)*internal0p_ZYZ(3)+dcos(internal0p_ZYZ(4)) &
                                        *DSQRT((1d0-internal0p_ZYZ(2)**2)*(1d0-internal0p_ZYZ(3)**2)))

                  internal0_BiSph(4) = pii- internal0p_ZYZ(5)
              
                  term1 = DCos(internal0p_ZYZ(5))*DSIN(internal0p_ZYZ(4))*DSQRT(1d0-internal0p_ZYZ(3)**2) &
                          + DSIN(internal0p_ZYZ(5))*(DCOS(internal0p_ZYZ(4))*internal0p_ZYZ(2)*DSQRT(1d0-internal0p_ZYZ(3)**2) &
                          - DSQRT(1d0-internal0p_ZYZ(2)**2)*internal0p_ZYZ(3) )
                          
              
              
                  term2 = DSQRT(1d0-internal0p_ZYZ(3)**2)*(DCOS(internal0p_ZYZ(4))* & 
                          internal0p_ZYZ(2)*DCos(internal0p_ZYZ(5)) - DSIN(internal0p_ZYZ(4))*DSIN(internal0p_ZYZ(5))) &
                          - internal0p_ZYZ(3)*DCos(internal0p_ZYZ(5))*DSQRT(1d0-internal0p_ZYZ(2)**2)
              
                  internal0_BiSph(5) = datan2(-1d0*term1,term2)
              end if
          end if
          
          
       
     

          

          !Output Alpha Check
         
          if (internal0_BiSph(4) > pii)then
            internal0_BiSph(4) = internal0_BiSph(4) -  2d0*pii 
          end if 
          if (internal0_BiSph(4) < -pii)then
            internal0_BiSph(4) = internal0_BiSph(4) + 2d0*pii 
          end if 
         
    END SUBROUTINE ZYZ_to_BiSpherical

    ! transformation from Autosurf internal angles coordinates to ZYZ angles
    SUBROUTINE Internal0_to_anglesZYZ(internal0,XDIM,angles)
        
        IMPLICIT NONE
        integer,INTENT(IN) ::  XDIM
        real*8,INTENT(IN) ::  internal0(XDIM)
        real*8,INTENT(INOut) ::  angles(6)
        
        real*8::pii
  
        pii = dacos(-1d0)
  
        
    
  
        if (XDIM==2)then
            angles(1)=0d0
            angles(2)=dacos(internal0(2))
            angles(3)=0d0
            angles(4)=0d0
            angles(5)=0d0
            angles(6)=0d0
        elseif (XDIM==3)then
            angles(1)=0d0
            angles(2)=dacos(internal0(2))
            angles(3)=pii - internal0(3)
            angles(4)=0d0
            angles(5)=0d0
            angles(6)=0d0
        elseif (XDIM==4)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=0d0
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=0d0     
        elseif (XDIM==5)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=internal0(5)
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=0d0 
        elseif (XDIM==6)then
            angles(1)=internal0(4)
            angles(2)=dacos(internal0(2))
            angles(3)=internal0(5)
            angles(4)=0d0
            angles(5)=dacos(internal0(3))
            angles(6)=internal0(6)                             
        endif
  
    
    END SUBROUTINE Internal0_to_anglesZYZ

    subroutine EulerAngles_2_BiSpherical(inter,R_ZYZ,internal0)
            ! *** internal coordinates ***
      ! ----------------------------
      
      !************** Subroutine definition 
          ! subroutine ZYZ --> BiSpherical (convert: Internal coordinates in ZYZ--> BiSpherical for 5D)  

          ! Input : R and Euler Angles in radians

          ! Output : R and cos_theta1,cos_theta2, phi_1,phi_2
          
      !********************************************************
      IMPLICIT NONE
      real*8,INTENT(IN):: R_ZYZ(7),inter(5)
      real*8,INTENT(OUT):: internal0(5)
      real*8::pii  ,alpha1,beta1,gamma1,alpha2,beta2,gamma2
      real*8:: term1,term2,threshold,alpha

      pii=dacos(-1d0) 
      
      internal0(1) = R_ZYZ(1)
      alpha1= R_ZYZ(2)
      beta1= R_ZYZ(3)
      gamma1= R_ZYZ(4)
      alpha2= R_ZYZ(5)
      beta2= R_ZYZ(6)
      gamma2= R_ZYZ(7)

      
      alpha=alpha1-alpha2
      if(alpha>pii)then
        alpha=alpha-2d0*pii
      endif
      if(alpha<-pii)then
        alpha=alpha+2d0*pii
      endif
      
      
      !write(*,*)"internal0_ZYZ",internal0_ZYZ
  
      threshold = 1d-8


      internal0(2) = beta1
      

      If (DABS(beta1)<threshold  .or. DABS(beta1-pii)< threshold)Then
          
          
          internal0(4) = inter(4)

          If (DABS(beta1)<threshold)Then
            internal0(2) = 0d0
            internal0(3) = beta2
          else
            internal0(2) = pii
            internal0(3) = pii - beta2
          end if


          If (DABS(beta2)<threshold .or. DABS(beta2-pii)<threshold)Then        
            internal0(5) = inter(5)
          else 
            If (DABS(beta1)<threshold)Then
              internal0(5) = -1d0*(alpha + gamma1)
            else
              internal0(5) = alpha + gamma1-pii  
            end if 
          end if

      else
        If (DABS(beta1)<2d-6  .or. DABS(beta1-pii)< 2d-6)Then
          internal0(4) = inter(4)
        else
          internal0(4) = pii- gamma1
        end if 
        

        
          
          If (DABS(beta2)<threshold .or. DABS(beta2-pii)<threshold)Then   
            internal0(5) = inter(5)     
            If (DABS(beta2)<threshold)Then
              internal0(3) = 0d0
            else
              internal0(3) = pii
            end if
          else 
            internal0(3) = DACOS(DCOS(beta1)*DCOS(beta2)+dcos(alpha)*DSIN(beta1)*DSIN(beta2))


            term1 = DCos(gamma1)*DSIN(alpha)*DSIN(beta2)+ DSIN(gamma1)*(DCOS(alpha)*Dcos(beta1)*DSIN(beta2) &
            - DSIN(beta1)*DCOS(beta2) )
            


            term2 = DSIN(beta2)*(DCOS(alpha)*DCos(beta1)*DCos(gamma1) - DSIN(alpha)*DSIN(gamma1)) &
                    - DCOS(beta2)*DCos(gamma1)*DSIN(beta1)

            internal0(5) = datan2(-1d0*term1,term2)
          end if

      end if


      

      !Output Alpha Check
     
      if (internal0(4) > pii)then
        internal0(4) = internal0(4) -  2d0*pii 
      end if 
      if (internal0(4) < -pii)then
        internal0(4) = internal0(4) + 2d0*pii 
      end if 

    
  
      return
    end subroutine EulerAngles_2_BiSpherical

    subroutine EulerAngles_2_Autosurf(XDIM,R_ZYZ,internal0)
      IMPLICIT NONE
      Integer,INTENT(IN)::XDIM
      real*8,INTENT(IN):: R_ZYZ(7)
      real*8,INTENT(OUT):: internal0(XDIM)
      real*8::pii  ,alpha1,beta1,gamma1,alpha2,beta2,gamma2
  
      pii=dacos(-1d0) 
      internal0(1) = R_ZYZ(1)
      alpha1= R_ZYZ(2)
      beta1= R_ZYZ(3)
      gamma1= R_ZYZ(4)
      alpha2= R_ZYZ(5)
      beta2= R_ZYZ(6)
      gamma2= R_ZYZ(7)
      ! *** internal coordinates ***
      ! ----------------------------
      if(XDIM==2)then
        internal0(2)=dcos(beta1)
      elseif(XDIM==3)then
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=pii-gamma1
      elseif(XDIM==4)then
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      elseif(XDIM==5)then
        internal0(5)=gamma1
        ! ----------------------------
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      elseif(XDIM==6)then
        internal0(5)=gamma1
        ! ----------------------------
        internal0(6)=gamma2
        ! ----------------------------
        internal0(2)=dcos(beta1)
        ! ----------------------------
        internal0(3)=dcos(beta2)
        ! ----------------------------
        internal0(4)=alpha1-alpha2
        if(internal0(4)>pii)then
          internal0(4)=internal0(4)-2d0*pii
        endif
        if(internal0(4)<-pii)then
          internal0(4)=internal0(4)+2d0*pii
        endif
      endif
  
    
  
      return
    end subroutine EulerAngles_2_Autosurf
  
    SUBROUTINE Int_to_Cart(internal,internal_length,mass,ref1,ref2,XDIM,natom1,natom2,internalFunction,cart)
          
          IMPLICIT NONE

          integer,INTENT(IN) :: internal_length,natom1,natom2,XDIM,internalFunction
          real*8,INTENT(IN) :: internal(internal_length),mass(natom1+natom2),ref1(natom1*3),ref2(natom2*3)
          integer :: intfunc
          real*8,allocatable,INTENT(OUT)  :: cart(:)
          real*8::ref1_temp0(natom1*3),ref2_temp0(natom2*3),cart_temp((natom1+natom2)*3)
          

          allocate(cart(3*(natom1+natom2)))
          intfunc = internalFunction
          !************** Subroutine definition 
              ! Place your custom function here! 
          !********************************************************
          ! make temporary copies of the input data
          ref1_temp0 = ref1! reference vector for frag1 (original frame)
          ref2_temp0 = ref2! reference vector for frag2 (original frame)

          if ( intfunc == 0)then
              call Int_to_Cart_User(internal,mass,cart_temp)
            elseif ( intfunc == 1)then
            ! Only for 3D: (for CH3CN-He system: internal coordinates taken as spherical coords.)
            call Int_to_Cart_Spherical(internal,cart_temp,mass,natom1,natom2,ref1_temp0)
            elseif ( intfunc == 2)then
            call Int_to_Cart_ZYZ(internal,XDIM,cart_temp,mass,natom1,natom2, ref1_temp0, ref2_temp0)
            ! Only if internal is already in cartesian
            elseif ( intfunc == -1)then
           
                    cart_temp = internal 
          endif

          cart = cart_temp
              
    END SUBROUTINE Int_to_Cart

    SUBROUTINE Int_to_Cart_ZYZ(internal,XDIM,cart,mass,natom1,natom2,ref1_0,ref2_0)

        !************** Subroutine definition 
            ! subroutine Int_to_Cart_ZYZ (convert: Internal coordinates in ZYZ--> Cartesian coordinates)  
            ! This function considers local axis parallels to ref1_0 and ref2_0 (systems known by Autosurf) 
            ! and global axis with origin in Center Mass of A with z-axis passing through Center Mass of B

            ! * the second angle beta  (rotation around Y axis) is given by cos_beta
            ! * the otehr two, alpha and gamma are given in radians
        !********************************************************
          use mathFunc
          IMPLICIT NONE
          integer,INTENT(IN) :: natom1,natom2,XDIM
          real*8,INTENT(IN) :: internal(XDIM),mass(natom1+natom2),ref1_0(natom1*3),ref2_0(natom2*3)
          real*8,INTENT(OUT)::cart((natom1+natom2)*3)
          integer :: i,natom
          real*8 :: cm(3)
          real*8 :: ref1(natom1*3),rotarr1(3,natom1),ref_mat1(3,natom1),arr1(natom1*3)
          real*8 :: ref2(natom2*3),rotarr2(3,natom2),R_arr(3,natom2),ref_mat2(3,natom2),arr2(natom2*3)
          real*8 :: R,angles(6),internal_temp(XDIM)
          
          internal_temp=internal
          natom=natom1+natom2
          ref1=ref1_0
          ref2=ref2_0

          R = internal_temp(1)

          if (XDIM==2)then
              angles(1)=0d0
              angles(2)=dacos(internal_temp(2))
              angles(3)=0d0
              angles(4)=0d0
              angles(5)=0d0
              angles(6)=0d0
          elseif (XDIM==3)then
              angles(1)=0d0
              angles(2)=dacos(internal_temp(2))
              angles(3)=internal_temp(3)
              angles(4)=0d0
              angles(5)=0d0
              angles(6)=0d0
          elseif (XDIM==4)then
              angles(1)=internal_temp(4)
              angles(2)=dacos(internal_temp(2))
              angles(3)=0d0
              angles(4)=0d0
              angles(5)=dacos(internal_temp(3))
              angles(6)=0d0     
          elseif (XDIM==5)then
              angles(1)=internal_temp(4)
              angles(2)=dacos(internal_temp(2))
              angles(3)=internal_temp(5)
              angles(4)=0d0
              angles(5)=dacos(internal_temp(3))
              angles(6)=0d0 
          elseif (XDIM==6)then
              angles(1)=internal_temp(4)
              angles(2)=dacos(internal_temp(2))
              angles(3)=internal_temp(5)
              angles(4)=0d0
              angles(5)=dacos(internal_temp(3))
              angles(6)=internal_temp(6)                                
          endif

          

          !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
          !!!    reference (origin at the new_CM of frag. 1)
          
          
          
          ! set new CM of fragment 1 at the origin
          call rm_cmass(ref1,mass(1:natom1),natom1,natom1) 
          call rm_cmass(ref2,mass(natom1+1:natom1+natom2),natom2,natom2) 



          call vec_to_mat2(ref1,ref_mat1,natom1)
          call vec_to_mat2(ref2,ref_mat2,natom2)

          call MolecularRotation_ZYZ(ref_mat1,natom1,angles(1),angles(2),angles(3),rotarr1)
          call MolecularRotation_ZYZ(ref_mat2,natom2,angles(4),angles(5),angles(6),rotarr2)

          

          call mat_to_vec2(rotarr1,arr1,natom1)
          R_arr(1:2,:) = rotarr2(1:2,:)
          R_arr(3,:) = rotarr2(3,:)+R
          call mat_to_vec2(R_arr,arr2,natom2)
          

          do i=1,3*natom1
              cart(i)=arr1(i)! Cartesian coordinates for the atoms in fragment 1
          enddo
          do i=1,3*natom2
              cart(3*natom1+i)=arr2(i)! Cartesian coordinates for the atoms in fragment 2
          enddo

        
    END SUBROUTINE Int_to_Cart_ZYZ

    SUBROUTINE Int_to_Cart_Spherical(internal,cart,mass,natom1,natom2,ref1_0)
          use mathFunc
          IMPLICIT NONE
          integer,INTENT(IN) :: natom1,natom2
          real*8,INTENT(IN) :: internal(3),mass(natom1+natom2),ref1_0(natom1*3)

          integer :: i,natom
          real*8 :: cm(3)
          real*8 :: ref1(natom1*3)
          real*8 :: ref2(natom2*3)
          real*8 :: cart((natom1+natom2)*3)
          real*8 :: sin_theta


          natom=natom1+natom2
          ref1=ref1_0

          !!! 1) find the Cartesian coordinates of all atoms in the NEW frame of 
          !!!    reference (origin at the new_CM of frag. 1)
          
          ! subroutine INT_Cart (convert: Internal coordinates --> Cartesian coordinates)  INT_Cart(internal,cart)
          ! (for CH3CN-He system: internal coordinates taken as spherical coords.)

          ! set new CM of fragment 1 at the origin
          call rm_cmass(ref1,mass(1:natom1),natom1,natom1) 
          do i=1,3*natom1
          cart(i)=ref1(i)! Cartesian coordinates for the atoms in fragment 1
          enddo
          sin_theta=dsqrt(1.d0-internal(2)**2)
          cm(1)=internal(1)*sin_theta*dcos(internal(3))
          cm(2)=internal(1)*sin_theta*dsin(internal(3))
          cm(3)=internal(1)*internal(2)
          ref2=cm
          do i=1,3
          cart(3*natom1+i)=ref2(i)! Cartesian coordinates for the atoms in fragment 2
          enddo




    END SUBROUTINE Int_to_Cart_Spherical 

    ! Int_to_Cart_User: Example for H2O-HCN
    subroutine Int_to_Cart_User(internal,mass,cart)
        implicit none
        integer,parameter :: XDIM=5,natom1=3,natom2=3
        real*8,parameter :: bohr2ang=0.529177249d0
        real*8,INTENT(IN) :: internal(XDIM)
        real*8,INTENT(IN) :: mass(natom1+natom2)
        real*8,INTENT(OUT) :: cart((natom1+natom2)*3)
        real*8 :: pii,rOH,thetaOH,Xcm,Ycm,Zcm,R,theta1,theta2,phi1,phi2
        real*8 :: r1,r2,rH0,rC1,rN1
        real*8 :: xO,yO,zO,xH1,yH1,zH1,xH2,yH2,zH2
        real*8 :: xH0,yH0,zH0,xC1,yC1,zC1,xN1,yN1,zN1
        real*8 :: mH1,mH2,mH3,mN,mO,mC,mtA,mtB
        real*8 :: RCMx,RCMy,RCMz
        
        pii=dacos(-1d0) 

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
        mtA = mH1+mH2+mO
        mtB = mN+mC+mH3

        ! H2O parameters
        rOH=1.8361d0              ! in bohr
        rOH=rOH*bohr2ang          ! convert to Ang.
        thetaOH=104.69d0/2.0d0    ! degrees
        thetaOH=thetaOH*pii/180d0 ! convert to radians
        ! HCN parameters
        r1=2.0286d0*bohr2ang      ! equilibre CH
        r2=2.1874d0*bohr2ang      ! equilibre CN  
        ! distance from each atom to the HCN center-of-mass
        rH0=((r1+r2)*mN+r1*mC)/(mtB) 
        rC1=(r2*mN - r1*mH3)/(mtB)
        rN1=((r1+r2)*mH3 + r2*mC)/(mtB)

        ! H2O Cartesian coordinates
        ! ----------------------------
        ! (H2O is fixed at the origin)
        xO=0d0 
        yO=0d0
        zO=(mH1+mH2)*rOH*dcos(thetaOH)/(mtA)
        xH1=-1d0*rOH*dsin(thetaOH)
        yH1=0d0
        zH1=zO-rOH*dcos(thetaOH)
        xH2=rOH*dsin(thetaOH)
        yH2=0d0
        zH2=zH1
        
        RCMx = (mO*xO + xH1*mH1 + xH2*mH2)/mtA
        RCMy = (mO*yO + yH1*mH1 + yH2*mH2)/mtA
        RCMz = (mO*zO + zH1*mH1 + zH2*mH2)/mtA
        
        
        cart(1)=xO  - RCMx
        cart(2)=yO  - RCMy
        cart(3)=zO  - RCMz
        cart(4)=xH1 - RCMx
        cart(5)=yH1 - RCMy
        cart(6)=zH1 - RCMz
        cart(7)=xH2 - RCMx
        cart(8)=yH2 - RCMy
        cart(9)=zH2 - RCMz
        !write(6,*)cart(1:3)
        !write(6,*)cart(4:6)
        !write(6,*)cart(7:9)
        
        ! internal coordinates
        R=internal(1)
        theta1=internal(2)
        theta2=internal(3)
        phi1=internal(4)
        phi2=internal(5)

        ! HCN Cartesian coordinates
        ! -------------------------
        Xcm=R*dsin(theta1)*dcos(phi1) 
        Ycm=R*dsin(theta1)*dsin(phi1) 
        Zcm=R*dcos(theta1) 
        xH0=Xcm-1.0d0*rH0*dsin(theta2)*dcos(phi2)
        yH0=Ycm-1.0d0*rH0*dsin(theta2)*dsin(phi2) 
        zH0=Zcm-1.0d0*rH0*dcos(theta2)
        xC1=Xcm-1.0d0*rC1*dsin(theta2)*dcos(phi2)
        yC1=Ycm-1.0d0*rC1*dsin(theta2)*dsin(phi2)
        zC1=Zcm-1.0d0*rC1*dcos(theta2)
        xN1=Xcm+1.0d0*rN1*dsin(theta2)*dcos(phi2)
        yN1=Ycm+1.0d0*rN1*dsin(theta2)*dsin(phi2)
        zN1=Zcm+1.0d0*rN1*dcos(theta2)
        cart(10)=xH0
        cart(11)=yH0
        cart(12)=zH0
        cart(13)=xC1
        cart(14)=yC1
        cart(15)=zC1
        cart(16)=xN1
        cart(17)=yN1
        cart(18)=zN1
        ! write(6,*)
        ! write(6,*)cart(3*natom1+1:3*natom1+3*natom2)
        !pause

        return

    end subroutine Int_to_Cart_User

end module coordinateTransf
    
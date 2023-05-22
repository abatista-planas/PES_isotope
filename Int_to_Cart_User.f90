! Units: angles in radians, distance in Angstroms
!
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
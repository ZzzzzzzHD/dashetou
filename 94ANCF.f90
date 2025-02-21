    SubRoutine ANCF_INIT
    Use Variable,Only:Eta_L,NSBN_ALE
    Use CONSTANT,Only:Re_Start,Irota
    Use ANCF
    Implicit None

    Integer::I

    AFP%NNODE  = AFP%NE +1
    AFP%L      = Eta_L
    AFP%LE     = AFP%L / AFP%NE

    Call CALSBE2SBN_Index

    ALLOCATE(AF_MASS(AFP%NNode*AFP%NumD,AFP%NNode*AFP%NumD))
    AF_MASS = 0.0D0
    Call ANCF_MAKE_MASS_MATRIX
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Allocate(AF_Fex(AFP%NNode*AFP%NumD))
    AF_Fex = 0.0D0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Allocate(AF_FeL(AFP%NNode*AFP%NumD))
    AF_FeL = 0.0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    NGQ1 = 3
    Allocate(GQ3xi(NGQ1),GQ3wei(NGQ1))
    Call ANCF_Gauss_Jacob(GQ3xi,GQ3wei,NGQ1)

    NGQ2 = 5
    Allocate(GQ5xi(NGQ2),GQ5wei(NGQ2))
    Call ANCF_Gauss_Jacob(GQ5xi,GQ5wei,NGQ2)


    If(Re_Start ==0)Then
        Allocate(AF_e(AFP%NNode*AFP%NumD),AF_ed(AFP%NNode*AFP%NumD),AF_edd(AFP%NNode*AFP%NumD))
        Allocate(AF_e0(AFP%NNode*AFP%NumD),AF_ed0(AFP%NNode*AFP%NumD),AF_edd0(AFP%NNode*AFP%NumD))
        AF_e        = 0.0
        AF_ed       = 0.0
        AF_edd      = 0.0
        AF_e0       = 0.0
        AF_ed0      = 0.0
        AF_edd0     = 0.0
        Do I = 1,AFP%NNode
            AF_e0(AFP%NumD*I-AFP%NumD+1) = (I-1) * AFP%LE ! X Coordinate, Y Initial equals to 0
            AF_e0(AFP%NumD*I-AFP%Dim+1) =  1.0
        End Do
        AF_e = AF_e0

        If((Irota == 5102).Or.(Irota == 520 ))Then
            Open(5200,file='CylinderY.dat')
            Write(5200,'(A50)')' VARIABLES = "T" "DIS_Y_Cylinder"'
            Open(5201,file='CylinderForceY.dat')
            Write(5201,'(A50)')' VARIABLES = "T" "Force_Y_Cylinder"'
        End If

        IF((IROTA == 510 ).Or.(IROTA == 5102 ).OR.(IROTA==5101).OR.(IROTA==520))THEN
        Open(5101,file='PlateX.dat')
        Open(5102,file='PlateY.dat')

        Write(5101,'(A50)')' VARIABLES = "T" "DIS_X"'
        Write(5102,'(A50)')' VARIABLES = "T" "DIS_Y"'

        Open(5181,file='PlateOsilation.dat')
        Open(5182,file='Force_X_OnPlate.dat')
        Open(5183,file='Force_Y_OnPlate.dat')
        Open(58845,file='Tail_Force.txt')
        END IF
    Else If(Re_Start==1)Then

        Allocate(AF_e(AFP%NNode*AFP%NumD),AF_ed(AFP%NNode*AFP%NumD),AF_edd(AFP%NNode*AFP%NumD))
        Allocate(AF_e0(AFP%NNode*AFP%NumD),AF_ed0(AFP%NNode*AFP%NumD),AF_edd0(AFP%NNode*AFP%NumD))
        AF_e        = 0.0
        AF_ed       = 0.0
        AF_edd      = 0.0
        AF_e0       = 0.0
        AF_ed0      = 0.0
        AF_edd0     = 0.0
        IF((IROTA == 510 ).Or.(IROTA == 5102 ).OR.(IROTA==5101).OR.(IROTA==520))THEN
        Open(5101,file='PlateX.dat',POSITION='APPEND')
        Open(5102,file='PlateY.dat',POSITION='APPEND')
        Open(5181,file='PlateOsilation.dat',POSITION='APPEND')
        Open(5182,file='Force_X_OnPlate.dat',POSITION='APPEND')
        Open(5183,file='Force_Y_OnPlate.dat',POSITION='APPEND')
        Open(58845,file='Tail_Force.txt',POSITION='APPEND')
        END IF
        
        If((Irota == 5102).Or.(Irota ==520 ))Then
            Open(5200,file='CylinderY.dat',POSITION='APPEND')
            Open(5201,file='CylinderForceY.dat',POSITION='APPEND')
        End If

    End If

    End

    SubRoutine ANCFMAIN
    Use Variable
    Use CONSTANT
    Use ANCF
    Use Lapack95
    Use Blas95
    Implicit None

    Integer::iter,maxIt,RateF2S,I
    Real(8)::h,tol,beta,gam,inv_beta,maxtol,AF_Q(AFP%NNode*AFP%NumD)
    Real(8)::R(AFP%NNode*AFP%NumD),AF_Fex_Sub(AFP%NNode*AFP%NumD),AF_Fex_Delta(AFP%NNode*AFP%NumD),Number,RateF2S_Real
    Real(8)::del(AFP%NNode*AFP%NumD),AF_e_n_1(AFP%NNode*AFP%NumD),AF_e_n_0(AFP%NNode*AFP%NumD)
    Real(8)::Mqddot(1,1),AF_edd_Matrix(AFP%NNode*AFP%NumD,1)

    Real(8)::PI

    PI = DACOS(-1.0D0)

    R           = 0.0D0
    AF_Fex_Sub  = 0.0D0
    del         = 0.0D0

    RateF2S = 20
    RateF2S_Real = 20.0D0
    Number = 1.0D0

    h           = DT / RateF2S_Real
    tol         = 1e-7
    maxIt       = 1000
    beta        = 0.3025D0
    gam         = 0.60D0
    inv_beta    = 1.0D0 / (beta * h**2.0D0)

    AF_Fex_Sub = AF_Fex
    Call ANCF_FORCE_Ex
    AF_Fex_Delta = (AF_Fex - AF_Fex_Sub )* Number / RateF2S_Real

    MASSR   = AFP%RHO

    Do I = 1,RateF2S

        AF_Fex_Sub = AF_Fex_Sub +  AF_Fex_Delta

        Call ANCF_FORCE_EL
        AF_Q = AF_Fex_Sub - AF_Fel
        Call SYSVX(AF_MASS,AF_Q,AF_edd0)

        AF_e = AF_e0
        AF_ed = AF_ed0
        AF_edd = AF_edd0

        Do iter = 1,maxIt
            AF_e  = AF_e0 + h*AF_ed0 + 0.50D0*h**2.0D0*((1.0D0-2.0D0*beta)*AF_edd0 + 2.0D0*beta*AF_edd)
            AF_ed = AF_ed0 + h*((1.0D0-gam)*AF_edd0 + gam*AF_edd)

            Call ANCF_FORCE_EL
            AF_Q = AF_Fex_Sub - AF_Fel

            Call Gemv(AF_MASS,AF_edd,R)
            R = R - AF_Q
            Call SYSVX(AF_MASS,R,del)
            AF_edd = AF_edd - del

            maxtol = maxval(Dabs(del))
            If(maxtol<=Tol)Then
                AF_e0 = AF_e
                AF_ed0 = AF_ed
                AF_edd0 = AF_edd
                Exit
            End If

            If(iter==maxIt) Then
                Print*,'iter Num reach maxIt, Maxtol : ',maxtol
                Pause
                Stop
            End If
        End do

        Call ANCF_FORCE_EL
        AF_edd_Matrix(:,1) = AF_edd(:)
        Call Gemm(AFMass2,AF_edd_Matrix,Mqddot)
        
        AFFex2 = 0.0D0
        !AFFel2 = 0.0D0
        !Mqddot = 0.0D0

        DXYK =(2.0D0*PI*FRE_Y)**2.0D0+8.0D0*PI*FRE_Y*DAMPING_Y/h+4.0D0/h/h
        DXYF = ForceYFlexi / ( MASSR * PI/4.0D0)  &
            & + (AFFex2 - AFFel2- Mqddot(1,1))* (AFP%RHO * AFP%A) /(MASSR * PI/4.0D0) &
            & + 4.0D0 * PI * FRE_Y * DAMPING_Y * ( 2.0D0 * DIS0_Y / h + VEL0_Y ) &
            & + 4.0D0 * DIS0_Y / h / h + 4.0D0 * VEL0_Y / h + ACC0_Y

        DIS_Y = DXYF / DXYK
        VEL_Y = 2.0D0 / h*(DIS_Y-DIS0_Y)-VEL0_Y
        ACC_Y = 4.0D0 / h / h * ( DIS_Y - DIS0_Y )- 4.0D0 / h * VEL0_Y - ACC0_Y

        AF_e(2)    = DIS_Y
        AF_ed(2)   = VEL_Y
        AF_edd(2)  = ACC_Y

        AF_e0(2)   = AF_e(2)
        AF_ed0(2)  = AF_ed(2)
        AF_edd0(2) = AF_edd(2)

        DIS0_Y = DIS_Y
        VEL0_Y = VEL_Y
        ACC0_Y = ACC_Y

    End Do
    Write(*,'(A10,I12,10X,A10,E12.4,5X,A10,F12.6)') 'ITER NUM =',iter, 'MAXTOL =',maxtol,'TAIL DIS =',AF_e(AFP%NNode*AFP%NumD-4)


    IF(MOD(NSTEP,5)==0)THEN
        WRITE(5101,'(4E18.10)') T,AF_e(AFP%NNode*AFP%NumD-5)
        WRITE(5102,'(4E18.10)') T,AF_e(AFP%NNode*AFP%NumD-4)
        Write(5200,'(4E18.10)') T,AF_e(2)
    END IF

    End SubRoutine ANCFMAIN

    !
    !
    !

    SubRoutine ANCF_Make_Mass_Matrix
    Use ANCF
    Use CONSTANT,Only:IROTA

    Implicit None

    Integer::I,J,K,II,JJ,KK
    Real(8)::Me(AFP%NumE,AFP%NumE),Temp,PI
    Integer::NGQ,ml,mr,il,ir                      !Node of GaussQuadrature
    Real(8),Allocatable::GQxi(:),GQwei(:)         !Local coordinate and weight of GaussQuadrature point
    Real(8)::Lo_xi,Lo_S(4),Lo_SI(AFP%DIM,AFP%NumE),Lo_SI_Tran(AFP%NumE,AFP%DIM)
    Character(12)::FileName

    !Ele Mass Matrix
    PI=DACOS(-1.0D0)
    NGQ = 4
    Allocate(GQxi(NGQ),GQwei(NGQ))
    Call ANCF_Gauss_Jacob(GQxi,GQwei,NGQ)

    ME = 0.0D0              ! Element Mass Matrix
    Do I = 1,NGQ
        Lo_xi = ( 1 + GQxi(I)) / 2.0 ! * AFP%LE will be divided in shape fun, so here just omit it
        Call ANCF_Shape_Fun(Lo_S,Lo_SI,Lo_xi,0)
        Lo_SI_Tran = Transpose(Lo_SI)
        Me = Me +  Matmul(Lo_SI_Tran, Lo_SI) * GQwei(I)
    End Do

    !Me = AFP%RHO *AFP%A * Me * AFP%Le *1.0D0/2.0D0            	! This is The Element Matrix, 1/2 from Gauss Integrand

    Me = Me * AFP%Le *1.0D0/2.0D0  								! This is Nondimensional Form

    IF(AFP%DIM == 2 )Then
        Do I = 0,AFP%NE -1
            AF_MASS(4*I+1:4*I+8,4*I+1:4*I+8) = AF_MASS(4*I+1:4*I+8,4*I+1:4*I+8) + Me
        End Do
    ElseIf(AFP%Dim == 3) Then
        Do I = 0,AFP%NE -1
            AF_MASS(6*I+1:6*I+12,6*I+1:6*I+12) = AF_MASS(6*I+1:6*I+12,6*I+1:6*I+12) + Me
        End Do
    End If


    ml = AFP%fixed(1);
    mr = AFP%fixed(2);

    If ((ml > 0).And.(IROTA /= 520).And.(IROTA /= 5102))Then
        AF_MASS(1:AFP%NumD,:) = 0.0
        AF_MASS(:,1:AFP%NumD) = 0.0
        Forall(I=1:AFP%NumD) AF_MASS(I,I) = 1.0

    Else IF (IROTA == 520)Then
        Temp = AF_MASS(2,2)
        AF_MASS(1,:) = 0.0D0
        AF_MASS(:,1) = 0.0D0
        AF_MASS(3:AFP%NumD,:) = 0.0D0
        AF_MASS(:,3:AFP%NumD) = 0.0D0
        Forall(I=1:AFP%NumD) AF_MASS(I,I) = 1.0
        !AF_MASS(2,2) = Temp  + (PI * 1.0D0**2.0D0)/(4.0D0 * AFP%A)
        !AF_MASS(2,2) = Temp * (Eta_T *1.0D0 / (PI * 1.0D0**2.0D0/4.0D0)) +  1.0D0
        !AF_MASS(2,2) = 1.0D0
        AF_MASS(2,2) = Temp + (PI * 1.0D0**2.0D0)/(4.0D0 * AFP%A)

    Else IF(IROTA == 5102)Then
        ALLOCATE(AFMass2(1,AFP%NNode*AFP%NumD))
        AFMass2(1,:)  = AF_MASS(2,:)
        AF_MASS(1:AFP%NumD,:) = 0.0
        AF_MASS(:,1:AFP%NumD) = 0.0
        Forall(I=1:AFP%NumD) AF_MASS(I,I) = 1.0
    End If

    !If (mr > 0) Then
    !    ir = AFP%NNODE*6;
    !    AF_MASS(ir+1:ir+mr,:) = 0.0
    !    AF_MASS(:,ir+1:ir+mr) = 0.0
    !    AF_MASS(ir+1:ir+mr,ir+1:ir+mr) = 1.0
    !End IF

    End SubRoutine ANCF_MAKE_MASS_MATRIX

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SubRoutine ANCF_FORCE_Ex
    Use ANCF
    Use Local
    Use blas95
    USE VARIABLE,ONLY: NPBE,PBE,NSBE,SBE,IPE,PN,X,Y,NODEELE,U,V,WEI,T,&
        & XPDKSI11,YPDITA22,XPDITA21,YPDKSI12,FPDX,FPDKSI,FPDITA,FPDY,TOR,M1,M2,M3,M4,&
        & DT,VEL0_ANG, VEL_ANG,&
        & WORK,WORK_TOR,WORK_F_Y,WORK_F_X,CL,CD,TOR0,CL0,CD0,&
        & VEL0_X,VEL_X,VEL0_Y,VEL_Y,DIS_X,DIS0_X,DIS_Y,DIS0_Y,DIS_ANG,DIS0_ANG,&
        & FORCE_X_CENTER_0,&
        & CL_CUT,ETA_R,ETA_T,ETA_L,X_ALE,Y_ALE,Nsbn_ORIG,ETA_T
    USE CONSTANT,ONLY: VISCOSITY,DENSITY,UUXX00,IROTA,LENGTH_LMM,MASSR,LL0,IROTA,NSTEP

    Implicit None
    !!!Local!!!
    INTEGER:: I,J,JJ,K,IE,NL,NK,IP,NODE1,NODE2,NET(4),MaxIndex
    REAL(8):: P1,P2,PM,X1,X2,Y1,Y2,DX,DY,DS,FNX,FNY,FTX,FTY,FD,Q,DVDX,DUDY,UT,VT,TAU,PP00,TT00
    REAL(8):: PX,PY,SX,SY,PX1,PY1,SX1,SY1,CD1,CL1,CD2,CL2,TORP,TORS,TOR1,PI
    REAL(8):: TOR_PIPE_TAU,TOR_PIPE_SIG

    INTEGER::PointCount,Index,Lo_flag,ForceLoadType
    Real(8)::lo_xi,Lo_X0,Force_Int(3,1),Force_Int_Tail(3,1),X_ex,Y_ex
    Real(8),Allocatable::Force_Tail(:,:)
    INTEGER::ElementIndex,NKINDSN,ml,mr
    REAL(8)::Flo_S(4),Flo_SI(AFP%Dim,AFP%NumE),Flo_SI_Tran(AFP%NumE,AFP%Dim),Fe(AFP%NumE,AFP%NE),Q_ORIG
    Real(8)::PFlexi(NsbN_ORIG,6),XFlexi_ORIG,YFlexi_ORIG,xFtemp,xxFtemp,Flo_xi
    Real(8)::M_Neu_X,M_Neu_Y,FaFx,FaFy,MXt,MYt
    Real(8),Allocatable :: PFlexi_1D(:,:),Mon2For(:,:)

    !!!Local!!!

    !PFlexi Elements
    !| 1 |  2 | 3  | 4 |    5   |   6   |
    !| X | Fx | Fy | y |  Dis_x | Dis_y |

    PI=DACOS(-1.0D0)


    !PP00=1/2*Rho_w*U^2*D (U:characteristical Velocity D:characteristical Length)
    PP00=1.0D0 !(0.5D0*DENSITY*UUXX00**2.0D0)
    !TT00=rho_w*J_total/rho_Structure*U^2/2/(Aera_Cylinder+Area_Plate)

    ! This is Dimensionless coeficience for Circular Cylinder With Plate
    ! TT00=((PI/64.0D0*LENGTH_LMM**2.0D0+0.13D0/12.0D0*LENGTH_LMM**2.0D0)/(PI/4.0D0+0.02D0)*DENSITY*UUXX00**2.0D0)
    ! This is Dimensionless coeficience for Rounded Square Cylinder With Plate
    TT00=DENSITY*(1.0D0/6.0D0-2.0D0*eta_r**4.0D0/3.0D0-2.0D0*eta_r**2.0D0*(1-eta_r)**2.0D0 &
        &+pi*eta_r**4.0D0/2.0D0-32.0D0*Eta_r**4.0D0/9.0D0/pi &
        &+2.0D0*pi*eta_r**2.0D0*(1.0D0/2.0D0-eta_r+4.0D0*eta_r/3.0D0/pi)**2.0D0+ETA_T*ETA_L**3.0D0/12.0D0+ETA_T*ETA_L*(1.0D0+ETA_L)**2.0D0/4.0D0)*UUXX00**2.0d0 &
        & /(2.0D0*(1.0D0-(4.0D0-PI)*ETA_R**2.0D0+eta_t*eta_l))

    PX=0.0D0
    PY=0.0D0

    SX=0.0D0
    SY=0.0D0

    TOR_PIPE_TAU=0.0D0
    TOR_PIPE_SIG=0.0D0

    IF ((IROTA==5102).Or.((IROTA==520)))Then
        DO I=1,NPBE
            IE=PBE(I,1)
            NL=PBE(I,2)
            IP=IPE(IE)

            IF (IP/=4) THEN
                WRITE(*,*)' FDVG9'
                PAUSE
                STOP
            END IF
            NET(1:4)=NODEELE(IE,1:4)

            IF(NL<=IP-1) THEN
                NODE1=NODEELE(IE,NL)
                NODE2=NODEELE(IE,NL+1)
            ELSE
                NODE1=NODEELE(IE,IP)
                NODE2=NODEELE(IE,1)
            END IF

            P1=PN(NODE1)
            P2=PN(NODE2)

            X1=X(NODE1)
            Y1=Y(NODE1)

            X2=X(NODE2)
            Y2=Y(NODE2)

            DX=X2-X1
            DY=Y2-Y1
            DS=DSQRT(DX*DX+DY*DY)

            CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)

            PM=0.5D0*(P1+P2)*DS!*1.0D0 B=1.0
            PX=PX+PM*FNX!正应力的x方向分力
            PY=PY+PM*FNY!正应力的y方向分力

            TOR_PIPE_SIG=TOR_PIPE_SIG-PM*FNX*(0.5D0*(Y1+Y2)-DIS0_Y)+PM*FNY*(0.5D0*(X1+X2)-DIS0_X)!正应力对圆柱的力矩（恒为0）

            CALL TXTY(X1,Y1,X2,Y2,FTX,FTY)

            CALL GAUSS_JACOB(5,5,IE)
            FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
            Q_ORIG=FD*WEI(3)*WEI(3)
            DVDX =0.0D0
            DUDY =0.0D0

            DO K=1,4
                FPDX(K)=(YPDITA22*FPDKSI(5,5,K)-YPDKSI12*FPDITA(5,5,K))/FD
                FPDY(K)=(XPDKSI11*FPDITA(5,5,K)-XPDITA21*FPDKSI(5,5,K))/FD
                NK=NET(K)
                UT=U(NK)
                VT=V(NK)
                DUDY=DUDY+FPDY(K)*UT
                DVDX=DVDX+FPDX(K)*VT
            END DO

            TAU=VISCOSITY*(DVDX-DUDY)*DS
            SX=SX+(-TAU*FTX)!切应力的x方向分力
            SY=SY+(-TAU*FTY)!切应力的y方向分力

            TOR_PIPE_TAU=TOR_PIPE_TAU-(-TAU*FTX)*(0.5D0*(Y1+Y2)-DIS0_Y)+(-TAU*FTY)*(0.5D0*(X1+X2)-DIS0_X) !切应力对圆柱的力矩
        END DO
        ForceXFlexi = (PX + SX)/PP00
        ForceYFlexi = (PY + SY)/PP00
        MomentFlexi = (TOR_PIPE_SIG + TOR_PIPE_TAU)/TT00
    End IF


    PX1 =   0.0D0
    PY1 =   0.0D0
    SX1 =   0.0D0
    SY1 =   0.0D0

    PFlexi = 0.0D0

    TORP    =   0.0D0
    TORS    =   0.0D0
    TOR     =   0.0D0

    MaxIndex = MaxVal(PlatePointIndex(:,1))
    Allocate(Force_Tail(Size(PlatePointIndex_Tail),4))
    Force_Tail = 0.0D0
    Allocate(PFlexi_1D(MaxIndex,3),Mon2For(MaxIndex,2))
    PFlexi_1D = 0.0D0
    Mon2For = 0.0D0

    Index = 0
    DO I=1,NSBE
        IE=SBE(I,1)!单元号
        NL=SBE(I,2)!起始的那个点号
        IP=IPE(IE) !节点数

        IF (IP/=4) THEN
            WRITE(*,*) 'FDVG10'
            PAUSE
            STOP
        END IF

        NET(1:4)=NODEELE(IE,1:4)

        IF(NL<=IP-1) THEN
            NODE1=NODEELE(IE,NL)
            NODE2=NODEELE(IE,NL+1)
        ELSE
            NODE1=NODEELE(IE,IP)
            NODE2=NODEELE(IE,1)
        END IF

        P1=PN(NODE1)
        P2=PN(NODE2)

        X1=X(NODE1)
        Y1=Y(NODE1)

        X2=X(NODE2)
        Y2=Y(NODE2)

        xFtemp  = X_ALE(NODE1)
        xxFtemp = X_ALE(NODE2)

        If(xFtemp < xxFtemp)Then
            XFlexi_ORIG=xFtemp
            YFlexi_ORIG=Y_ALE(NODE1)
        Else
            XFlexi_ORIG=xxFtemp
            YFlexi_ORIG=Y_ALE(NODE2)
        End If

        !IF(DABS(XFlexi_ORIG-ETA_L-0.5D0)<=1E-5)THEN
        !    FFLGG=1
        !END IF

        DX=X2-X1
        DY=Y2-Y1

        DS=DSQRT(DX*DX+DY*DY)

        CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)

        PM=0.5D0*(P1+P2)*DS
        !!!!!!!!!!!!!!!!
        !PX1=PX1+PM*FNX!
        !PY1=PY1+PM*FNY!
        !!!!!!!!!!!!!!!!
        PX1=PM*FNX
        PY1=PM*FNY

        !JY!尾板压力对刚（圆）心的力矩=====
        TORP=TORP-PM*FNX*(0.5D0*(Y1+Y2)-DIS0_Y)+PM*FNY*(0.5D0*(X1+X2)-DIS0_X)              !板上法向压力对圆心的力矩  !JY!认为圆心是刚心

        CALL TXTY(X1,Y1,X2,Y2,FTX,FTY)

        CALL GAUSS_JACOB(5,5,IE)
        FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
        Q_ORIG=FD*WEI(3)*WEI(3)
        DVDX =0.0D0
        DUDY =0.0D0

        DO K=1,4
            FPDX(K)=(YPDITA22*FPDKSI(5,5,K)-YPDKSI12*FPDITA(5,5,K))/FD
            FPDY(K)=(XPDKSI11*FPDITA(5,5,K)-XPDITA21*FPDKSI(5,5,K))/FD
            NK=NET(K)
            UT=U(NK)
            VT=V(NK)
            DUDY=DUDY+FPDY(K)*UT
            DVDX=DVDX+FPDX(K)*VT
        END DO

        TAU=VISCOSITY*(DVDX-DUDY)*DS
        SX1=(-TAU*FTX)
        SY1=(-TAU*FTY)

        PFlexi(I,1)= XFlexi_ORIG
        PFlexi(I,2)=(PX1+SX1)/PP00 !Fx
        PFlexi(I,3)=(PY1+SY1)/PP00 !Fy
        PFlexi(I,4)= YFlexi_ORIG

        J = SBE2SBN_Index(I)
        JJ = PlatePointIndex(J,1)

        PFlexi_1D(JJ,1) = PFlexi(I,1)
        PFlexi_1D(JJ,2) = PFlexi_1D(JJ,2) + PFlexi(I,2)
        PFlexi_1D(JJ,3) = PFlexi_1D(JJ,3) + PFlexi(I,3)

        If(JJ == MaxIndex) Then
            X_ex = X_Neutal(JJ,1)
            Y_ex = Y_Neutal(JJ,1)

            !下面这个主要是为了控制求出来的力矩的方向正确
            PFlexi(I,5) = XFlexi_ORIG - X_ex !DIS From Surface to Neutal Point,X
            PFlexi(I,6) = Y_ex - YFlexi_ORIG !DIS From Surface to Neutal Point,Y

            Index = Index + 1
            Force_Tail(Index,1) = YFlexi_ORIG
            Force_Tail(Index + 1,1) = Y_ALE(NODE1)!只有最后一个不会被覆盖

            Force_Tail(Index,2) = PFlexi(I,2) !Fx
            Force_Tail(Index,3) = PFlexi(I,3) !Fy
            Force_Tail(Index,4) = 0.0D0       !Fz = 0

        Else
            X_ex = X_Neutal(JJ,1)
            Y_ex = Y_Neutal(JJ,1)

            !下面这个主要是为了控制求出来的力矩的方向正确
            PFlexi(I,5) = XFlexi_ORIG - X_ex !DIS From Surface to Neutal Point,X
            PFlexi(I,6) = Y_ex - YFlexi_ORIG !DIS From Surface to Neutal Point,Y

        End If

        !注意！！这里其实不需要特殊处理，因为少记录的点其实是在尾板上，这里不用考虑，后面有处理。
        !而且这里的处理不对，加错点了，但是也没有影响，在第一个单元，最后也是置零了。
        !但是现在还是不要处理吧，尾板上面的力另做处理。

        !IF (I==NSBE) Then
        !    PFlexi(I+1,1)=X_ALE(NODE1)
        !    PFlexi(I+1,2)=(PX1+SX1)/PP00
        !    PFlexi(I+1,3)=(PY1+SY1)/PP00
        !    PFlexi(I+1,4)=Y_ALE(NODE2)
        !    J = SBE2SBN_Index(I+1)!??????????????20211220
        !    JJ = PlatePointIndex(J,1)
        !    PFlexi_1D(JJ,1) = PFlexi(I+1,1)
        !    PFlexi_1D(JJ,2) = PFlexi_1D(JJ,2) + PFlexi(I+1,2)
        !    PFlexi_1D(JJ,3) = PFlexi_1D(JJ,3) + PFlexi(I+1,3)
        !End if

        !尾板切应力对刚（圆）心的力矩,在柔性程序里面是不需要的
        !TORS=TORS-(-TAU*FTX)*(0.5D0*(Y1+Y2)-DIS0_Y)+(-TAU*FTY)*(0.5D0*(X1+X2)-DIS0_X)  !板上切向粘性剪切力对圆心的力矩
    END DO


    !!!!投机取巧的方法!!!!!
    PFlexi_1D(MaxIndex,2)= 0.0D0!PFlexi_1D(MaxIndex-1,2) + (PFlexi_1D(MaxIndex-1,2)-PFlexi_1D(MaxIndex-2,2))
    PFlexi_1D(MaxIndex,3)= 0.0D0!PFlexi_1D(MaxIndex-1,3) + (PFlexi_1D(MaxIndex-1,3)-PFlexi_1D(MaxIndex-2,3))

    !NGQ1 = 3
    !Allocate(GQ3xi(NGQ1),GQ3wei(NGQ1))
    !Call ANCF_Gauss_Jacob(GQ3xi,GQ3wei,NGQ1)

    !!!!!!!!!Cal Force On The Tail!!!!!!!!!!!!!!DO NOT DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !Way2!!!!
    !Force_Tail(Index+1,2) = 0.0D0 !Force_Tail(Index,2)!在同一个单元！
    !Force_Tail(Index+1,3) = 0.0D0 !Force_Tail(Index,3)
    !Force_Tail(Index+1,4) = 0.0D0 !Fz = 0
    !
    !Force_Int_Tail(1,1) = Sum(Force_Tail(:,2))
    !Force_Int_Tail(2,1) = Sum(Force_Tail(:,3))
    !
    !PFlexi_1D(MaxIndex,2) = 0.0D0!Force_Int_Tail(1,1)  !Fx
    !PFlexi_1D(MaxIndex,3) = 0.0D0!Force_Int_Tail(2,1)  !Fy
    !write(58845,'(*(E12.4,1x))') T,PFlexi_1D(MaxIndex,2),PFlexi_1D(MaxIndex,3),PFlexi_1D(MaxIndex-1,2),PFlexi_1D(MaxIndex-1,3)
    !!!!Way2!!!!

    !!!Force_Tail(Index+1,1)
    !Force_Tail(Index+1,2) = Force_Tail(Index,2)!在同一个单元
    !Force_Tail(Index+1,3) = Force_Tail(Index,3)
    !Force_Tail(Index+1,4) = 0.0D0 !Fz = 0
    !Force_Tail(:,1) = Force_Tail(:,1) + 0.5D0*Eta_T       !平移到0-Eta_T之间
    !
    !Force_Int_Tail = 0.0D0
    !Do I = 1, NGQ1
    !    lo_xi = (1.0D0 + GQ3xi(I)) * Eta_T / 2.0D0 / Eta_T!0-1
    !    Lo_X0 = lo_xi * Eta_T
    !
    !    Lo_flag = 0
    !    Findp: Do K= 1,Size(Force_Tail(:,1))
    !        IF (Force_Tail(K,1) == Lo_X0) Then!不像另外一个，这里必然有个0（0.5）是相等的
    !            Lo_flag = 1
    !            Exit Findp
    !        ElseIF (( Force_Tail(K,1) < Lo_X0).And.( Lo_X0< Force_Tail(K+1,1))) Then
    !            Lo_flag = 2
    !            Exit Findp
    !        End IF
    !    End Do Findp
    !
    !    IF(Lo_flag == 1)Then
    !        Force_Int(1,1) = Force_Tail(K,2)   !Fx
    !        Force_Int(2,1) = Force_Tail(K,3)   !Fy
    !        Force_Int(3,1) = Force_Tail(K,4)   !Fz = 0
    !
    !    ElseIF(Lo_flag == 2)Then
    !        Force_Int(1,1) = (Lo_X0 -Force_Tail(K + 1,1)) / (Force_Tail(K,1) - Force_Tail(K+1,1)) * Force_Tail(K,2) &
    !            & + (Lo_X0 -Force_Tail(K,1)) / (Force_Tail(K+1,1) - Force_Tail(K,1)) * Force_Tail(K+1,2)
    !        Force_Int(2,1) = (Lo_X0 -Force_Tail(K + 1,1)) / (Force_Tail(K,1) - Force_Tail(K+1,1)) * Force_Tail(K,3) &
    !            & + (Lo_X0 -Force_Tail(K,1)) / (Force_Tail(K+1,1) - Force_Tail(K,1)) * Force_Tail(K+1,3)
    !        Force_Int(3,1) = 0.0D0
    !    End If
    !    Force_Int_Tail = Force_Int_Tail + Force_Int * GQ3wei(I)
    !End Do
    !Force_Int_Tail = Force_Int_Tail * 1.0D0/2.0D0 *Eta_T
    !
    !PFlexi_1D(MaxIndex,2) = Force_Int_Tail(1,1)  !Fx
    !PFlexi_1D(MaxIndex,3) = Force_Int_Tail(2,1)  !Fy
    !write(58845,'(*(E12.4,1x))') PFlexi_1D(MaxIndex,2),PFlexi_1D(MaxIndex,3),PFlexi_1D(MaxIndex-1,2),PFlexi_1D(MaxIndex-1,3)
    !!!!!!!!Cal Force On The Tail!!!!!!!!!!!!!!DO NOT DELETE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF(Mod(NSTEP,50)==0)Then
        Write(5182,*)'Zone T="',NSTEP,'"'
        Do I =1,MaxIndex
            Write(5182,*) PFlexi_1D(I,1),PFlexi_1D(I,2)
        End Do
        Write(5182,*)
        Write(5183,*)'Zone T="',NSTEP,'"'
        Do I =1,MaxIndex
            Write(5183,*) PFlexi_1D(I,1),PFlexi_1D(I,3)
        End Do
        Write(5183,*)
    End IF

    ForceLoadType = 2
    IF(ForceLoadType==1)Then
        !Distributed Force!
        AF_Fex = 0.0d0
        Fe = 0.0d0
        Do I = 1 ,AFP%NE
            Force_Int = 0.0
            Do J = 1,NGQ1
                lo_xi = (1.0D0 + GQ3xi(J)) * AFP%LE / 2.0D0 / AFP%LE
                Lo_X0 = (I -1) * AFP%LE + 0.5D0 + lo_xi * AFP%LE

                Flo_S   = 0.0D0
                Flo_SI  = 0.0D0
                Flo_SI_Tran = 0.0D0
                Call ANCF_Shape_Fun(Flo_S,Flo_SI,lo_xi,0)
                Flo_SI_Tran = Transpose(Flo_SI)

                IF((Lo_X0 < PFlexi_1D(1,1)).Or.(Lo_X0 > PFlexi_1D(MaxIndex,1))) Then !比第一个还要小
                    Write(*,*) 'Mesh of Plate is Too Sparse'
                    Pause
                    Stop
                End If

                Lo_flag = 0
                Findp2: Do K= 1,MaxIndex  !查找在哪两个点之间，进行线性插值，可以优化到只循环一次，存在矩阵里面
                    IF (( PFlexi_1D(K,1) < Lo_X0).And.( Lo_X0< PFlexi_1D(K+1,1))) Then
                        Lo_flag = K
                        Exit Findp2
                    End IF
                End Do Findp2!!!!!!!!有点问题有点问题有点问题

                !f_Int = (x-x1)/(x0-x1)*y0 + (x-x0)/(x1-x0)*y1 线性插值
                Force_Int(1,1) = (Lo_X0 -PFlexi_1D(Lo_flag + 1,1)) / (PFlexi_1D(Lo_flag,1) - PFlexi_1D(Lo_flag+1,1)) * PFlexi_1D(Lo_flag,2) &
                    & + (Lo_X0 -PFlexi_1D(Lo_flag,1)) / (PFlexi_1D(Lo_flag+1,1) - PFlexi_1D(Lo_flag,1)) * PFlexi_1D(Lo_flag+1,2)
                Force_Int(2,1) = (Lo_X0 -PFlexi_1D(Lo_flag + 1,1)) / (PFlexi_1D(Lo_flag,1) - PFlexi_1D(Lo_flag+1,1)) * PFlexi_1D(Lo_flag,3) &
                    & + (Lo_X0 -PFlexi_1D(Lo_flag,1)) / (PFlexi_1D(Lo_flag+1,1) - PFlexi_1D(Lo_flag,1)) * PFlexi_1D(Lo_flag+1,3)
                Force_Int(3,1) = 0.0D0

                ! int S' * q
                Do K = 1,AFP%NumE
                    Fe(K,I)= Fe(K,I) + GQ3wei(J) * &
                        & Flo_SI_Tran(K,1)*Force_Int(1,1)+Flo_SI_Tran(K,2)*Force_Int(2,1)+Flo_SI_Tran(K,3)*Force_Int(3,1)
                End do
            End Do
            ! ds = AFP%LE, 1/2 from Gauss intergration
            Fe(:,I) = Fe(:,I) * AFP%LE * 1.0D0 / 2.0D0

            Index = AFP%NumD*(I-1)
            AF_Fex(Index+1:Index+AFP%NumE) = AF_Fex(Index+1:Index+AFP%NumE) + Fe(:,I)
        End Do
        !Distributed Force!
    Else
        !Concentrated Force!
        ElementIndex = 0.0d0
        Fe = 0.0d0
        AF_Fex = 0.0d0
        DO i = 1, MaxIndex
            ElementIndex=Floor((PFlexi_1D(i,1)-0.5d0)/AFP%LE)
            IF (ElementIndex<AFP%NE)THEN
                Flo_S   = 0.0D0
                Flo_SI  = 0.0D0
                Flo_SI_Tran = 0.0D0
                !(PFlexi_1D(i,1) + PFlexi_1D(i,1)) / 2.0D0  Move The Action Point To The Middle Of The Mesh
                Flo_xi = ((PFlexi_1D(I,1) + PFLEXI_1D(I+1,1)) / 2.0D0 -0.5d0-(ElementIndex)*AFP%LE)/AFP%LE
                Call ANCF_Shape_Fun(Flo_S,Flo_SI,Flo_xi,0)
                Flo_SI_Tran = Transpose(Flo_SI)

                Do j = 1,AFP%NumE
                    Fe(j,ElementIndex+1)=Flo_SI_Tran(j,1)*PFlexi_1D(i,2)+Flo_SI_Tran(j,2)*PFlexi_1D(i,3)!S' * q
                End do

                AF_Fex(AFP%NumD*ElementIndex+1:AFP%NumD*ElementIndex+AFP%NumE)&
                    &   =  AF_Fex(AFP%NumD*ElementIndex+1:AFP%NumD*ElementIndex+AFP%NumE)+Fe(:,ElementIndex+1)!&
                !&   * (PFlexi_1D(I+1,1) - PFlexi_1D(I,1)) / AFP%LE ! Single Mesh Width - Non-Dimentional
            ELSE
                Flo_S   = 0.0D0
                Flo_SI  = 0.0D0
                Flo_SI_Tran = 0.0D0
                Flo_xi = (PFlexi_1D(i,1)-0.5d0-(ElementIndex-1)*AFP%LE)/AFP%LE
                Call ANCF_Shape_Fun(Flo_S,Flo_SI,Flo_xi,0)
                Flo_SI_Tran = Transpose(Flo_SI)

                Do j = 1,AFP%NumE
                    Fe(j,ElementIndex)=Flo_SI_Tran(j,1)*PFlexi_1D(i,2)+Flo_SI_Tran(j,2)*PFlexi_1D(i,3)
                End do

                AF_Fex(AFP%NumD*(ElementIndex-1)+1:AFP%NumD*(ElementIndex-1)+AFP%NumE) &
                    & = AF_Fex(AFP%NumD*(ElementIndex-1)+1:AFP%NumD*(ElementIndex-1)+AFP%NumE)+Fe(:,ElementIndex) !&
                !& * (PFlexi_1D(I,1) - PFlexi_1D(I-1,1)) / AFP%LE
            End if
        End do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Concentrated Force!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    End IF

    !Bench2
    !AF_Fex = AF_Fex * 1.18e-3 * 1.0D0**2.0D0  / (AFP%RHO * AFP%A) !/ 10.0D0  !D = 1.0
    !Sun
    AF_Fex = AF_Fex * DENSITY * 1.0D0**2.0D0  / (AFP%RHO * AFP%A) !/ 10.0D0  !D = 1.0

    ml = AFP%fixed(1)
    mr = AFP%fixed(2)

    If  ((ml > 0).And.(IROTA /= 520).And.(IROTA /= 5102))Then
        AF_Fex(1:AFP%NumD) = 0.0
    Else IF (IROTA == 520)Then
        !AF_Fex(2) = 1.0D0 * AF_Fex(2)&
        !& + ForceYFlexi *DENSITY * 1.0D0**2.0D0 / (AFP%RHO * AFP%A)
        AF_Fex(2) = ForceYFlexi *DENSITY * 1.0D0**2.0D0 / (AFP%RHO * AFP%A)
        AF_Fex(1) = 0.0
        AF_Fex(3:AFP%NumD) = 0.0
    Else IF (IROTA == 5102)Then
        AFFex2 = AF_Fex(2)
        AF_Fex(1:AFP%NumD) = 0.0
    End If

    IF(MOD(NSTEP,5)==0)THEN
        WRITE(5201,'(4E18.10)') T,ForceYFlexi
    END IF

    DEAllocate(PFlexi_1D)

    End SubRoutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SubRoutine  ANCF_FORCE_EL
    Use CONSTANT,Only:UUXX00,DAMPING_Y,FRE_Y,IROTA
    Use VARIABLE,oNLY:Eta_T
    Use ANCF
    Use Blas95
    Use Local

    Implicit None
    Integer::ml,mr                          !Node of GaussQuadrature
    Real(8)::Lo_xi,Lo_S(1,4),Lo_S_V(4),Lo_S_Tran(4,1),Lo_SI(AFP%Dim,AFP%NumE),Lo_SI_Tran(AFP%NumE,AFP%Dim)
    Real(8)::Lo_SS(1,4),Lo_SS_V(4),Lo_SS_Tran(4,1),Lo_SSI(AFP%Dim,8),Lo_SSI_Tran(8,AFP%Dim)
    Integer::I,Is,Ie,J,K,Index
    Real(8)::FeL_e(AFP%NumE,1),FeL_e1(AFP%NumE,1),FeL_e2(AFP%NumE,1),e_ele(AFP%NumE),d_ele(4,AFP%Dim),d_ele_Tran(AFP%Dim,4),StrainD(1,AFP%NumE),StrainD_Tran(AFP%NumE,1),Strain(1,1),eee(AFP%NumE,1)
    Real(8)::rr(AFP%Dim),r_x(AFP%Dim,1),r_x_Tran(1,AFP%Dim),r_x_v(AFP%Dim),r_xx(AFP%Dim,1),r_xx_V(AFP%Dim),f1(AFP%Dim,1),f1_Tran(1,AFP%Dim),f1_v(AFP%Dim),f,g1,g,kappa,g_e(1,AFP%NumE),f_e(1,AFP%NumE)
    Real(8)::coeMA(AFP%Dim,AFP%Dim),coeMB(AFP%Dim,AFP%Dim),coeMAA(AFP%Dim,AFP%NumE),coeMBB(AFP%Dim,AFP%NumE),K_K_e(AFP%NumE,1)!注意这里的kke已经转置过了，不用再转置，kke这个形式就是正确的20211215
    Real(8)::coeMCC(AFP%Dim,AFP%NumE),coeMDD(AFP%Dim,AFP%NumE)
    Real(8)::coef(1),coe1g2(1),dcoef,dcoe1g2,ds
    Real(8)::e_ele_M(AFP%NumE,1),Matrix(AFP%NumE,AFP%NumE)
    Real(8)::r_xrep(AFP%Dim,AFP%NumE),r_xxrep(AFP%Dim,AFP%NumE),Ones(1,AFP%NumE),fel1(AFP%Dim,AFP%NumE),fel2(AFP%Dim,AFP%NumE)
    Real(8)::PI

    PI=DACOS(-1.0D0)


    AF_FeL = 0.0

    Do I = 1, AFP%NE
        Is = AFP%NumD*I - AFP%NumD + 1
        Ie = AFP%NumD*(I+1)

        !!!!!!!!!!!!!!!!!!!!
        e_ele = AF_e(Is:Ie)
        !!!!!!!!!!!!!!!!!!!!

        FeL_e1 = 0.0
        Do J = 1,NGQ2
            lo_xi = (1.0D0 + GQ5xi(J)) * AFP%LE / 2.0D0 / AFP%LE
            lo_S   = 0.0D0
            lo_SI  = 0.0D0
            lo_SI_Tran = 0.0D0
            Call ANCF_Shape_Fun(lo_S,lo_SI,lo_xi,1)
            !Lo_S_Tran = TransPose(lo_S)
            !Lo_SI_Tran = TransPose(lo_SI

            Lo_S_V =Reshape(Lo_S,(/4/))

            Call eps_eps_e(Lo_S_V,e_ele,eee)
            eee = eee * GQ5wei(J)
            FeL_e1 = FeL_e1 + eee

        End Do
        !FeL_e1 =  100.0D0 *  FeL_e1 *AFP%A * AFP%E* AFP%LE * 1.0D0 / 2.0D0
        FeL_e1 =  5.0D0 *  FeL_e1 *AFP%A * AFP%E / (AFP%RHO * AFP%A * UUXX00**2.0D0) * AFP%LE * 1.0D0 / 2.0D0  !zc
        !Bench 2
        !FeL_e1 =  FeL_e1 *AFP%A * AFP%E / (AFP%RHO * AFP%A * 51.3D0**2.0D0) * AFP%LE * 1.0D0 / 2.0D0  !zc
        !* AFP%LE * 1.0D0 / 2.0D0 -- Gauss Integeration
        !*AFP%A * AFP%E / (AFP%RHO * AFP%A * UUXX00**2.0D0) -- Dimensionaless coefficient
        !FeL_e1 =  100.0D0 * FeL_e1 * AFP%E / (AFP%RHO * U_inf**2.0D0) * AFP%LE * 1.0D0 / 2.0D0
        !FeL_e1 =  100.0D0 * FeL_e1  / AFP%RHO  * AFP%E  * AFP%LE * 1.0D0 / 2.0D0
        !FeL_e1 =  FeL_e1 *AFP%A * AFP%E / (AFP%RHO  * UUXX00**2.0D0) * AFP%LE * 1.0D0 / 2.0D0

        FeL_e2 = 0.0
        K_K_e = 0.0

        Do J = 1 , NGQ1
            lo_xi = (1.0D0 + GQ3xi(J)) * AFP%LE / 2.0D0 / AFP%LE

            lo_S   = 0.0D0
            lo_S_Tran = 0.0D0
            lo_SI  = 0.0D0
            lo_SI_Tran = 0.0D0
            Call ANCF_Shape_Fun(lo_S,lo_SI,lo_xi,1)
            !lo_S_Tran = Transpose(lo_S)
            !lo_SI_Tran = Transpose(lo_SI)

            lo_SS   = 0.0D0
            lo_SS_Tran = 0.0D0
            lo_SSI  = 0.0D0
            lo_SSI_Tran = 0.0D0
            Call ANCF_Shape_Fun(lo_SS,lo_SSI,lo_xi,2)
            !lo_SS_Tran = Transpose(lo_SS)
            !lo_SSI_Tran = Transpose(lo_SSI)

            Lo_S_V = Reshape(Lo_S,(/4/))
            Lo_SS_V = Reshape(Lo_SS,(/4/))

            Call S_e(Lo_S_V,e_ele,r_x)
            Call S_e(Lo_SS_V,e_ele,r_xx)

            r_x_V = Reshape(r_x,(/AFP%Dim/))
            r_xx_V = Reshape(r_xx,(/AFP%Dim/))

            r_x_Tran = Transpose(r_x)

            Call Cross(r_x,r_xx,f1)
            f1_Tran = Transpose(f1)
            f1_V = Reshape(f1,(/AFP%Dim/))
            coef = 3.0D0 * reshape(Matmul(f1_Tran,f1),(/1/)) / reshape(Matmul(r_x_Tran,r_x),(/1/))
            coe1g2 = 1.0D0 / reshape(Matmul(r_x_Tran,r_x),(/1/))**3.0D0
            !dcoef = dsqrt(Sum(f1(1:3,1)**2.0D0))
            !dcoe1g2 = dsqrt(Sum(r_x(1:3,1)**2.0D0))**2.0D0
            !dcoe1g2 = dsqrt(Sum(r_x(1:3,1)**2.0D0))**3.0D0
            dcoef = coef(1)
            dcoe1g2 = coe1g2(1)

            Call kappa_kappa_e(Lo_S_V,Lo_SS_V,dcoef,dcoe1g2,r_x_V,r_xx_V,f1_V,K_K_e)
            FeL_e2 = FeL_e2 + K_K_e * GQ3wei(J)

        End Do
        !FeL_e2 = FeL_e2 * AFP%E * AFP%I * AFP%LE *1.0D0 / 2.0D0
        FeL_e2 = FeL_e2 * AFP%E * AFP%I  /(1.0D0-0.35D0**2.0D0)/ (AFP%RHO * AFP%A * UUXX00**2.0D0 * 1.0D0**2.0D0 ) * AFP%LE *1.0D0 / 2.0D0   !zc
        !Bench 2
        !FeL_e2 = FeL_e2 * AFP%E * AFP%I / (AFP%RHO * AFP%A * 51.3D0**2.0D0 * 1.0D0**2.0D0 ) * AFP%LE *1.0D0 / 2.0D0   !zc
        ! Nondimentional Equation Form Goes as Follows
        !FeL_e2 = FeL_e2 / 0.84D0* AFP%E * Eta_T**2.0D0 /(12.0D0 *AFP%RHO * U_inf**2.0D0 * 1.0D0**2.0D0 ) * AFP%LE * 1.0D0 / 2.0D0 !D = 1.0
        !FeL_e2 = FeL_e2 * AFP%E * AFP%I / (AFP%RHO  * UUXX00**2.0D0 * 1.0D0**2.0D0 ) * AFP%LE * 1.0D0 / 2.0D0 !D = 1.0
        !FeL_e2 = FeL_e2 / 0.84D0  * Eta_T**2.0D0 /(12.0D0 * 1.0D0 **2.0D0) *  1.0D0 / AFP%RHO  * AFP%E * AFP%LE * 1.0D0 / 2.0D0

        Index = 0
        Do K = Is,Ie
            Index = Index +1
            AF_FeL(K) = AF_FeL(K) + FeL_e1(Index,1) + FeL_e2(Index,1)
        End Do

    End Do

    ml = AFP%fixed(1)
    mr = AFP%fixed(2)

    IF  ((ml > 0).And.(IROTA /= 520).And.(IROTA /= 5102))Then
        AF_FeL(1:AFP%NumD) = 0.0
    Else IF (IROTA == 520) Then
        AF_FeL(1) = 0.0D0
        AF_FeL(3:AFP%NumD) = 0.0D0
        !AF_FeL(2) =  1.0D0 * AF_FeL(2)  &
        !    & + 1.0D0**2.0D0 / AFP%A * PI**2.0D0 * DAMPING_Y * FRE_Y * AF_ed(2) &
        !    & + 1.0D0**2.0D0 / AFP%A * PI**3.0D0 * FRE_Y** 2.0D0 * AF_e(2)
        !AF_FeL(2) =  -1.0D0 * AF_FeL(2)  &
        !    & + 1.0D0**2.0D0 / AFP%A * PI**2.0D0 * DAMPING_Y * FRE_Y * AF_ed(2) &
        !    & +1.0D0**2.0D0 / AFP%A * PI**3.0D0 * FRE_Y** 2.0D0 * AF_e(2)
        !AF_FeL(2) =  1.0D0 * AF_FeL(2) * (Eta_T *1.0D0 / (PI * 1.0D0**2.0D0/4.0D0))  &
        !    & + 4.0D0 * PI * DAMPING_Y * FRE_Y * AF_ed(2) &
        !    & +(2.0D0 * PI * FRE_Y) ** 2.0D0 * AF_e(2)
        !AF_FeL(2) = 4.0D0 * PI * DAMPING_Y * FRE_Y * AF_ed(2) &
        !& +(2.0D0 * PI * FRE_Y) ** 2.0D0 * AF_e(2)
        AF_FeL(2) =  1.0D0 * AF_FeL(2)  &
            & + PI / 4.0D0 * (  4.0D0 * PI * DAMPING_Y * FRE_Y * AF_ed(2) &
            &                 + 4.0D0 * PI **2.0D0 * FRE_Y** 2.0D0 * AF_e(2))
    Else IF (IROTA == 5102) Then
        AFFel2 = AF_FeL(2)
        AF_FeL(1:AFP%NumD) = 0.0
    END IF

    End SubRoutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SubRoutine eps_eps_e(Sx,e,eee)
    Implicit None

    Real(8),Intent(IN)::Sx(4),e(12)
    Real(8),Intent(OUT)::eee(12,1)
    Real(8)::t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10, t11, t12, t13, t14, t17, t18, t19, t20, t21, t22, t23, t24
    Real(8)::t15, t29, t30, t31, t32, t33, t34, t35, t36, t16, t25, t26, t27, t28, t37, t38, t39

    t2 = Sx(1)
    t3 = e(1)
    t4 = t2*t3
    t5 = Sx(2)
    t6 = e(4)
    t7 = t5*t6
    t8 = Sx(3)
    t9 = e(7)
    t10 = t8*t9
    t11 = Sx(4)
    t12 = e(10)
    t13 = t11*t12
    t14 = t4+t7+t10+t13
    t17 = e(2)
    t18 = t2*t17
    t19 = e(5)
    t20 = t5*t19
    t21 = e(8)
    t22 = t8*t21
    t23 = e(11)
    t24 = t11*t23
    t15 = t18+t20+t22+t24
    t29 = e(3)
    t30 = t2*t29
    t31 = e(6)
    t32 = t5*t31
    t33 = e(9)
    t34 = t8*t33
    t35 = e(12)
    t36 = t11*t35
    t16 = t30+t32+t34+t36
    t25 = t14**2
    t26 = t25*(1.0/2.0)
    t27 = t15**2
    t28 = t27*(1.0/2.0)
    t37 = t16**2
    t38 = t37*(1.0/2.0)
    t39 = t26+t28+t38-1.0/2.0
    eee(:,1) = [t2*t14*t39,&
        & t2*t15*t39,&
        & t2*t16*t39,&
        & t5*t14*t39,&
        & t5*t15*t39,&
        & t5*t16*t39,&
        & t8*t14*t39,&
        & t8*t15*t39,&
        & t8*t16*t39,&
        & t11*t14*t39,&
        & t11*t15*t39,&
        & t11*t16*t39]
    End SubRoutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SubRoutine kappa_kappa_e(Sx,Sxx,coef,inv_rx_rx_cubed,rx,rxx,v,k_k_e)
    Implicit None

    Real(8)::t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t20, t21, t22, t23
    Real(8)::t24, t25, t26, t27, t28, t29, t30, t31, t32, t33, t34, t35, t36, t37, t38, t39, t40, t41, t42
    Real(8)::Sx(4),Sxx(4),coef,inv_rx_rx_cubed,rx(3),rxx(3),v(3),k_k_e(12,1)

    t2 = Sx(1)
    t3 = Sxx(1)
    t4 = rx(1)
    t5 = v(3)
    t6 = rxx(3)
    t7 = t2*t6
    t8 = rx(3)
    t9 = t7-t3*t8
    t10 = rx(2)
    t11 = rxx(1)
    t12 = t3*t4
    t13 = v(2)
    t14 = rxx(2)
    t15 = t3*t10
    t16 = v(1)
    t17 = t15-t2*t14
    t18 = t12-t2*t11
    t19 = Sx(2)
    t20 = Sxx(2)
    t21 = t6*t19
    t22 = t21-t8*t20
    t23 = t4*t20
    t24 = t10*t20
    t25 = t24-t14*t19
    t26 = t23-t11*t19
    t27 = Sx(3)
    t28 = Sxx(3)
    t29 = t6*t27
    t30 = t29-t8*t28
    t31 = t4*t28
    t32 = t10*t28
    t33 = t32-t14*t27
    t34 = t31-t11*t27
    t35 = Sx(4)
    t36 = Sxx(4)
    t37 = t6*t35
    t38 = t37-t8*t36
    t39 = t4*t36
    t40 = t10*t36
    t41 = t40-t14*t35
    t42 = t39-t11*t35
    k_k_e = reshape([-inv_rx_rx_cubed*(t5*t17+t9*t13+coef*t2*t4),&
        &inv_rx_rx_cubed*(t5*t18+t9*t16-coef*t2*t10),&
        & -inv_rx_rx_cubed*(t13*t18-t16*t17+coef*t2*t8),&
        &-inv_rx_rx_cubed*(t5*t25+t13*t22+coef*t4*t19),&
        &inv_rx_rx_cubed*(t5*t26+t16*t22-coef*t10*t19),&
        &-inv_rx_rx_cubed*(t13*t26-t16*t25+coef*t8*t19),&
        &-inv_rx_rx_cubed*(t5*t33+t13*t30+coef*t4*t27),&
        &inv_rx_rx_cubed*(t5*t34+t16*t30-coef*t10*t27),&
        &-inv_rx_rx_cubed*(t13*t34-t16*t33+coef*t8*t27),&
        &-inv_rx_rx_cubed*(t5*t41+t13*t38+coef*t4*t35),&
        &inv_rx_rx_cubed*(t5*t42+t16*t38-coef*t10*t35),&
        &-inv_rx_rx_cubed*(t13*t42-t16*t41+coef*t8*t35)],(/12,1/))

    End SubRoutine

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SubRoutine S_e(S,e,Calout)
    Implicit None

    Real(8),Intent(In)::S(4),e(12)
    Real(8),Intent(Out)::Calout(3,1)

    CalOut = Reshape([  (S(1)*e(1) + S(2)*e(4) + S(3)*e(7) + S(4)*e(10)),&
        & (S(1)*e(2) + S(2)*e(5) + S(3)*e(8) + S(4)*e(11)),&
        & (S(1)*e(3) + S(2)*e(6) + S(3)*e(9) + S(4)*e(12))],(/3,1/))
    End SubRoutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SubRoutine ANCF_Gauss_Jacob(XXI,Fweight,Ng)
    Implicit None

    Integer::Ng
    Real(8)::XXI(Ng),Fweight(Ng)

    Select Case(Ng)
    case(2)
        XXI(1)=-dsqrt(1.0D0 /3.0D0 )
        XXI(2)=dsqrt(1.0D0 /3.0D0 )

        Fweight(1)=   1.0D0
        Fweight(2)=   1.0D0
    Case(3)
        XXI(1)=-dsqrt(3.0D0 /5.0D0 )
        XXI(2)=0.0D0
        XXI(3)=dsqrt(3.0D0 /5.0D0 )

        Fweight(1)=5.0D0/9.0D0
        Fweight(2)=8.0D0/9.0D0
        Fweight(3)=5.0D0/9.0D0
    Case(4)
        XXI(1)=-dsqrt(3.0D0 /7.0D0  + dsqrt(120.0D0 )/35.0D0 )
        XXI(2)=-dsqrt(3.0D0 /7.0D0  - dsqrt(120.0D0 )/35.0D0 )
        XXI(3)= dsqrt(3.0D0 /7.0D0  - dsqrt(120.0D0 )/35.0D0 )
        XXI(4)= dsqrt(3.0D0 /7.0D0  + dsqrt(120.0D0 )/35.0D0 )

        Fweight(1)=   1.0D0 /2.0D0  - 5.0D0 /(3.0D0 *dsqrt(120.0D0 ))
        Fweight(2)=   1.0D0 /2.0D0  + 5.0D0 /(3.0D0 *dsqrt(120.0D0 ))
        Fweight(3)=   1.0D0 /2.0D0  + 5.0D0 /(3.0D0 *dsqrt(120.0D0 ))
        Fweight(4)=   1.0D0 /2.0D0  - 5.0D0 /(3.0D0 *dsqrt(120.0D0 ))
    Case(5)
        XXI(1)=-(dsqrt(5.0D0  + 2.0D0 *dsqrt(10.0D0 / 7.0D0 )))/3.0D0
        XXI(2)=-(dsqrt(5.0D0  - 2.0D0 *dsqrt(10.0D0  / 7.0D0 )))/3.0D0
        XXI(3)=0.0D0
        XXI(4)=dsqrt(5.0D0  - 2.0D0 *dsqrt(10.0D0 /7.0D0 ))/3.0D0
        XXI(5)=dsqrt(5.0D0 +2.0D0 *dsqrt(10.0D0 /7.0D0 ))/3.0D0

        Fweight(1)=(322.0D0 - 13.0D0 *dsqrt(70.0D0 ))/900.0D0
        Fweight(2)=(322.0D0  + 13.0D0 *dsqrt(70.0D0 ))/900.0D0
        Fweight(3)=128.0D0  / 225.0D0
        Fweight(4)=(322.0D0 + 13.0D0 *dsqrt(70.0D0 ))/900.0D0
        Fweight(5)=(322.0D0 - 13.0D0 *dsqrt(70.0D0 ))/900.0D0
    End Select

    End SubRoutine

    SubRoutine ANCF_Shape_Fun(S,SI,XI,N)
    Use ANCF,Only:AFP

    Implicit  None
    Integer::N,I
    Real(8)::S(4), XI,Ieye(AFP%DIM,AFP%DIM),SI(AFP%DIM,AFP%NumE)
    Real(8),Allocatable::SIM(:)

    Allocate(SIM(AFP%DIM*AFP%NumE))

    s   = 0.0d0
    SIM = 0.0d0
    SI  = 0.0d0

    Select case(N)
    case(0)
        s(1) = 1.0D0 - 3.0D0*xi**2.0D0 + 2.0D0*xi**3.0D0
        s(2) = AFP%LE * (xi - 2.0D0*xi**2.0D0 + xi**3.0D0)
        s(3) = 3.0D0*xi**2.0D0 - 2.0D0*xi**3.0D0
        s(4) = AFP%LE * (-xi**2.0D0 + xi**3.0D0)
    case(1)
        s(1) = (6.0D0*xi**2.0D0 - 6.0D0*xi)/AFP%LE
        s(2) = 1.0D0 - 4.0D0*xi + 3.0D0*xi**2.0D0
        s(3) = (-6.0D0*xi**2.0D0 + 6.0D0*xi)/AFP%LE
        s(4) = -2.0D0*xi + 3.0D0*xi**2.0D0
    case(2)
        s(1) = (12.0D0*xi-6.0D0)/AFP%LE**2.0D0
        s(2) = (-4.0D0+6.0D0*xi)/AFP%LE
        s(3) = (6.0D0-12.0D0*xi)/AFP%LE**2.0D0
        S(4) = (-2.0D0+6.0D0*xi)/AFP%LE
    End Select

    Ieye = 0.0
    Forall(I=1:AFP%DIM) Ieye(I,I) = 1.0
    SIM =[s(1)*Ieye,s(2)*Ieye,s(3)*Ieye,s(4)*Ieye]
    SI = Reshape(SIM,(/AFP%DIM,AFP%NumE/))

    DEAllocate(SIM)
    End SubRoutine

    SubRoutine cross(a, b,c)
    Real(8), DIMENSION(3) :: c
    Real(8), DIMENSION(3), INTENT(IN) :: a, b

    c(1) = a(2) * b(3) - a(3) * b(2)
    c(2) = a(3) * b(1) - a(1) * b(3)
    c(3) = a(1) * b(2) - a(2) * b(1)
    END SubRoutine cross


    SubRoutine CALSBE2SBN_Index
    Use Ancf
    Use VARIABLE,Only:NSBE,NSBN,SBE,SBN,NODEELE,IPE,X_ALE
    Implicit None
    Integer::I,J,K,INo,IE,NL,IP,NODE1,NODE2

    Allocate(SBE2SBN_Index(NSBN))

    Do I = 1,NSBE
        IE=SBE(I,1)
        NL=SBE(I,2)
        IP=IPE(IE)

        IF (IP/=4) THEN
            WRITE(*,*)' FDVG10'
            PAUSE
            STOP
        END IF

        IF(NL<=IP-1) THEN
            NODE1=NODEELE(IE,NL)
            NODE2=NODEELE(IE,NL+1)
        ELSE
            NODE1=NODEELE(IE,IP)
            NODE2=NODEELE(IE,1)
        END IF

        If(X_ALE(Node1)<=X_ALE(Node2))Then!使用x小的点作为单元的代表值
            Do J = 1,NSBN
                INo = SBN(J)
                If (NODE1==INo) Then
                    SBE2SBN_Index(I) = J
                End If
            End Do
        Else
            Do J = 1,NSBN
                INo = SBN(J)
                If (NODE2==INo) Then
                    SBE2SBN_Index(I) = J
                End If
            End Do
        End If

        !注意，下面的处理方式不对，漏掉的那个点不在最后一个单元上，在最尾端那里。
        !If (I==NSBE) Then
        !    Do J = 1,NSBN
        !        INo = SBN(J)
        !        If (NODE2==INo) Then
        !            SBE2SBN_Index(I+1) = J
        !        End If
        !    End Do
        !End IF

        SBE2SBN_Index( I + 1 ) = 0!暂时不做处理,暂时不需要

    End Do


    End SubRoutine

    SubRoutine CAlNumEMeshPoint!Not in Use
    Use ANCF
    Use LOCAL
    Use VARIABLE,Only:NSBN_ALE,SBN,X_ALE

    Implicit None
    Integer::I,J,K,Num,ElementIndex,IIA(AFP%NumE)
    Real(8)::X0
    Real(8)::CopyPlatePointIndex(NSBN_ALE,2)

    Num = Ceiling(MaxVal(PlatePointIndex(:,1)) / AFP%NumE) + 1
    Allocate(NumEMeshPoint(Num,AFP%NumE))
    NumEMeshPoint = 0.0D0

    CopyPlatePointIndex = PlatePointIndex
    !!CopyPlatePointIdex 一开始是以LabelID升序排列的
    !|----------|-------|
    !|EleIndex  |LabelId|
    !|  56      |   1   |
    !|  72      |   2   |
    !|  63      |   3   |
    !......

    IIA = 0
    Do I = 1,Size(PlatePointIndex_UDF(:,1))
        J = PlatePointIndex_UDF(I,1) !这个I就是LabelID
        X0 = X_ALE(J)
        ElementIndex=Ceiling((X0-0.5d0)/AFP%LE)
        IIA(ElementIndex) = IIA(ElementIndex) + 1!避免数据覆盖
        If(IIA(ElementIndex)>Num)Then
            Write(*,*) 'ERR#84511'
            Pause
            Stop
        Else
            NumEMeshPoint(IIA(ElementIndex),ElementIndex) = PlatePointIndex(I,1)!现在要存的是中轴线的Index
        End IF
    End Do

    End SubRoutine


    !SubRoutine MatrxiOutPut(FileName,Matrix,Rows,Cloums)
    !Implicit None
    !
    !Integer::Unit,Rows,Cloums,I
    !Real(8)::Seed,Matrix(Rows,Cloums)
    !Character(12)::FileName
    !Character(16)::FileNameL
    !
    !Unit = 125062
    !
    !Write(FileNameL,'(*(A))') FileName,'.txt'
    !Unit = ABS(Floor(Seed * 1500))
    !Open(Unit,File = FileNameL)
    !Do I = 1, Rows
    !    Write(Unit,'(*(f12.4,2X))') Matrix(I,:)
    !End Do
    !Close(Unit)
    !
    !End SubRoutine MatrxiOutPut
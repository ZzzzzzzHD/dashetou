    SubRoutine FlexibleMain
    !use lapack95
    Use VARIABLE,Only:dt,Dis_y,Dis_ang,Vel_ang,Vel_Y
    Use FlexibleOnly

    Implicit None

    !!!!!!!Local!!!!!!!!!
    Integer::i
    Integer,Allocatable::ipiv(:)
    Real(8),Allocatable::ORKeff(:,:)
    Real(8)::gam,beta,h,A1,A2,A3,A4,A5,A6,C_DAMPING(FL_Params_n,FL_Params_n)
    
    Real(8)::FX0,R_FX0,ElementIndex,FLO_XI,FLO_S(4)
    !!!!!!!!!!!!!!!!!!!!!

    !!New-Mark beta!!!!
    gam = 0.6D0
    beta            =   0.3025d0
    h               =   dt
    A1              =   1.0d0/(beta * h**2)
    A2              =   gam / (beta * h)
    A3              =   1 / (beta * h);
    A4              =   1/(2*beta)-1;
    A5              =   (gam/(2*beta)-1)*h;
    A6              =   gam/beta-1;
    !!New-Mark beta!!!!


    !!等效刚度矩阵!!
    C_DAMPING = DampingCoeiff*M_MASSM
    
    Keff = K_StiffM + A1*M_MASSM + A2 * C_DAMPING
    Allocate(ORKeff(FL_Params_n,FL_Params_n))
    ORKeff =  Keff

    Call Flexible_Force

    Foq     = Fq !上一步的位移速度和加速度数据
    Foqd    = Fqd
    Foqdd   = Fqdd

    Feff = FEULERBER + matmul(M_MASSM , (A1*Foq+A3*Foqd+A4*Foqdd)) + matmul(C_DAMPING , (A2*Foq+A6*Foqd+A5*Foqdd))

    Feff(1) = 0
    Feff(2) = 0
    
    !Allocate(G(FL_Params_n,FL_Params_n))  !这个是用来储存A矩阵巧利斯基分解后形成的上三角矩阵
    !G    = 0.0d0
    
    Allocate(Ipiv(FL_Params_n))
    Ipiv = 0 
    
    !Call cholesky_d(FL_Params_n, Keff, G, Feff)
    
    !Call lu(KEFF,IPIV)
    
    Call matrixinv(ORKEFF,KEFF,FL_Params_n)
    !call getrf(Keff,ipiv)
    !call getri(Keff,ipiv)
    
    Fq  = matmul(Keff,Feff)
    Fqdd =   A1*(Fq-Foq) - A3*Foqd- (1/(2*beta)-1)*Foqdd
    Fqd  =   Foqd + (1-gam)*h*Foqdd + gam*h*Fqdd

    Dis_y   =Fq(1)
    Dis_ang =Fq(2) !这一个刚好可以作为圆柱的y和angle
    Vel_ang =Fqd(1) 
    Vel_Y   =Fqd(2)

    
    Fdis0_Y      =Fdis_Y
    Fdis0_ang    =Fdis_ang
    FVel0_Y      =FVel_Y
    FVel0_ang    =FVel_ang
    !write(*,*) Dis_y
    Do i = 1,FL_Params_ne+1!加速度就没有管了，似乎后面没有用
        Fdis_Y(i)        = Fq(2*i-1)
        Fdis_ang(i)      = Fq(2*i)  !这个数组可以在后面全部储存尾板上的数据，就是需要处理第一个元素就行，利用几何关系
        FVel_Y(i)        = Fqd(2*i-1)
        FVel_ang(i)      = Fqd(2*i)
    End do
    
    !Fdis_Y(1)=0
    !Fdis_ang(1)=0
    !FVel_Y(1)=0
    !FVel_ang(1)=0
    !write(*,*) Dis_y,Dis_ang
    !write(*,'(*(f20.17,3X))') Fdis_Y
    !write(*,'(*(f20.17,3X))') Fdis0_Y
    !write(*,'(*(f12.6,3X))') FVel_Y
    !write(*,'(*(f12.6,3X))') Fqdd
    
    !Do I = 1 ,size(X_CA_ORR(:,1))
    !    FX0 = X_CA_ORR(I,1) !这是中轴线的x坐标
    !    R_FX0 = FX0 - 0.5
    !    ElementIndex = FLOOR(R_FX0/FL_PARAMS_L)
    !    IF (ElementIndex < FL_PARAMS_NE) Then
    !        FLO_XI=(R_FX0-ElementIndex*FL_PARAMS_L)/FL_PARAMS_L
    !        FLO_S = 0.0D0
    !        CALL SHAPE_FUN(FLO_S,4,FLO_XI,0)
    !        X_CA_ORR(I,2) = FLO_S(1)*FDIS_Y(ElementIndex+1)+FLO_S(2)*FDIS_ANG(ElementIndex+1)+ &
    !            & FLO_S(3)*FDIS_Y(ElementIndex+2)+FLO_S(4)*FDIS_ANG(ElementIndex+2)
    !    Else
    !        FLO_XI=(R_FX0-(ELEMENTINDEX-1)*FL_PARAMS_L)/FL_PARAMS_L
    !        FLO_S=0.0D0
    !        CALL SHAPE_FUN(FLO_S,4,FLO_XI,0)
    !        X_CA_ORR(I,2) = FLO_S(1)*FDIS_Y(ELEMENTINDEX)+FLO_S(2)*FDIS_ANG(ELEMENTINDEX)+ &
    !            & FLO_S(3)*FDIS_Y(ELEMENTINDEX+1)+FLO_S(4)*FDIS_ANG(ELEMENTINDEX+1)
    !    End If
    !End Do
    

    DEAllocate(ORKeff)
    End SubRoutine


    
    
    
    
    SubRoutine Make_Mass_Matrix
    Use FlexibleOnly

    Implicit None
    !!!Local!!!
    Integer::Ng,Ie,istart,iend,kin,ml,mr,il,ir,i,j
    Real(8)::Me(4,4),mm(4,4),FLo_xi,FLo_S(1,4),FLo_Tran_S(4,1),PI
    Real(8),Allocatable::XXI(:),Fweight(:)
    !!!

    Ng = 3
    Allocate(XXI(Ng))
    Allocate(Fweight(Ng))

    Call Gauss_Jacob_Fexible1D(XXI,Fweight,Ng)

    Do Ie =1,FL_Params_ne
        istart = 2 * ie -1
        iend   = 2 * ie +2

        Me     = 0.0d0

        Do kin = 1, Ng
            !xi=(1 + XXI(k)) * Params_L / 2.0d0
            FLo_xi=(1 + XXI(kin)) / 2.0d0
            Call Shape_Fun(FLo_S,4, FLo_xi, 0)

            FLo_Tran_S = transpose(FLo_S)
            mm = matmul(FLo_Tran_S, FLo_S)
            Me = Me + Fweight(kin) * mm;
        End do

        Me = Me *0.5D0
        !Me = Me * FL_params_rho * FL_params_A * FL_params_L !注意这是有量纲方程这样弄，无量纲的不一样

        M_MASSM(istart:iend, istart:iend) = M_MASSM(istart:iend, istart:iend) + Me;




    End do
    
    ml = Fixed(1)
    mr = Fixed(2)

    if ((ml - 1 <= 1e-8 ).Or.(ml - 2 <= 1e-8 )) then
        il = 0;
        M_MASSM(il+1:il+ml,:) = 0.0D0
        M_MASSM(:,il+1:il+ml) = 0.0D0
        DO I = il+1,il+ml
            DO J = il+1,il+ml
                If (I==J) then
                    M_MASSM(I,J) = 1.0d0;
                end if
            end Do
        end Do
    Else if(ml - 3 <= 1e-8) then
        M_MASSM(2,:) = 0.0d0
        M_MASSM(:,2) = 0.0d0
        M_MASSM(2,2) = 1.0d0
    End if

    
    if ((mr - 1 <= 1e-8 ).Or.(mr - 2 <= 1e-8 )) then
        ir = FL_params_n - 2;
        M_MASSM(ir+1:ir+mr,:) = 0.0D0
        M_MASSM(:,ir+1:ir+mr) = 0.0D0
        DO I = ir+1,ir+mr
            DO J = ir+1,ir+mr
                If (I==J) then
                    M_MASSM(I,J) = 1.0d0;
                end if
            end Do
        end Do

    Else if(ml - 3 <= 1e-8) then
        M_MASSM(FL_params_n,:) = 0.0d0
        M_MASSM(:,FL_params_n) = 0.0d0
        M_MASSM(FL_params_n,FL_params_n) = 1
    End if

    DEAllocate(XXI)
    DEAllocate(Fweight)

    open(517,file='MassMatrix.txt')
    do i = 1, FL_Params_n
        write(517,'(*(f12.4,1X))') M_MASSM(i,:)
    End do
    close(517)

    End SubRoutine

    SubRoutine Gauss_Jacob_Fexible1D(XXI,Fweight,Ng)
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

    SubRoutine Shape_Fun(S,Elements, XI, N)
    USE FlexibleOnly,ONLY:FL_Params_L
    Implicit  None


    Integer::N,Elements
    Real(8)::S(Elements), XI

    s = 0.0d0
    Select case(N)
    case(0)
        s(1) = 1.0D0 - 3.0D0*xi**2.0D0 + 2.0D0*xi**3.0D0
        s(2) = FL_Params_L * (xi - 2.0D0*xi**2.0D0 + xi**3.0D0)
        s(3) = 3.0D0*xi**2.0D0 - 2.0D0*xi**3.0D0
        s(4) = FL_Params_L * (-xi**2.0D0 + xi**3.0D0)

    case(1)
        s(1) = (6.0D0*xi**2.0D0 - 6.0D0*xi)/FL_Params_L
        s(2) = 1.0D0 - 4.0D0*xi + 3.0D0*xi**2.0D0
        s(3) = (-6.0D0*xi**2.0D0 + 6.0D0*xi)/FL_Params_L
        s(4) = -2.0D0*xi + 3.0D0*xi**2.0D0

    case(2)
        s(1) = (12.0D0*xi-6.0D0)/FL_Params_L**2.0D0
        s(2) = (-4.0D0+6.0D0*xi)/FL_Params_L
        s(3) = (6.0D0-12.0D0*xi)/FL_Params_L**2.0D0
        S(4) = (-2.0D0+6.0D0*xi)/FL_Params_L
    End Select



    End SubRoutine
    
    SubRoutine Make_Stiff_Matrix
    Use FlexibleOnly
    Use Constant 
    Implicit None
    !!!Local!!!
    Integer::Ng,ie,istart,iend,kin,il,ir,i,j,ml,mr
    Real(8)::PI,Ke(4,4),kk(4,4),Fxi,SXX(1,4),SXX_TRAN(4,1)
    Real(8),Allocatable::XXI(:),Fweight(:)
    !!!Local!!!
    
    pi=dacos(-1.0d0)
    Ng = 3
    Allocate(XXI(Ng))
    Allocate(Fweight(Ng))
    XXI(Ng)   = 0.0d0
    Fweight(Ng) = 0.0d0
    Call Gauss_Jacob_Fexible1D(XXI,Fweight,Ng)
    
    SXX=0.0d0
    Fxi=0.0d0
    
    Do ie = 1,FL_Params_ne
        istart = 2*ie - 1;
        iend = 2*ie + 2;
        Ke = 0.0d0;

        Do kin = 1,3
            Fxi = (1 + XXI(kin)) / 2
            call shape_fun(SXX,4, Fxi, 2)

            SXX_TRAN = transpose(SXX)
            kk = matmul(SXX_TRAN, SXX)
            Ke = Ke + Fweight(kin) * kk
        end do
        Ke = Ke * FL_Params_E * FL_Params_I/(MASSR * DENSITY * FL_params_A * FL_params_L * uuxx00**2 ) * 0.5D0

        K_StiffM(istart:iend,istart:iend) = K_StiffM(istart:iend, istart:iend) + Ke  
    End do
    
    ml = Fixed(1)
    mr = Fixed(2)
    
    if ((ml - 1 <= 1e-8 ).Or.(ml - 2 <= 1e-8 )) then
        il = 0;
        K_StiffM(il+1:il+ml,:) = 0.0D0
        K_StiffM(:,il+1:il+ml) = 0.0D0
        DO I = il+1,il+ml
            DO J = il+1,il+ml
                If (I==J) then
                    K_StiffM(I,J) = 1.0d0;
                end if
            end Do 
        end Do
    Else if(ml - 3 <= 1e-8) then
        K_StiffM(2,:) = 0.0d0
        K_StiffM(:,2) = 0.0d0
        K_StiffM(2,2) = 1.0d0
    End if
    
      
     if ((mr - 1 <= 1e-8 ).Or.(mr - 2 <= 1e-8 )) then
        ir = FL_Params_n - 2;
        K_StiffM(ir+1:ir+mr,:) = 0.0D0
        K_StiffM(:,ir+1:ir+mr) = 0.0D0
         DO I = ir+1,ir+mr
            DO J = ir+1,ir+mr
                If (I==J) then
                    K_StiffM(I,J) = 1.0d0;
                end if
            end Do 
         end Do

    Else if(ml - 3 <= 1e-8) then
        K_StiffM(FL_Params_n,:) = 0.0d0
        K_StiffM(:,FL_Params_n) = 0.0d0
        K_StiffM(FL_Params_n,FL_Params_n) = 1
    End if
    
    DEAllocate(XXI)
    DEAllocate(Fweight)    
    
    
    End SubRoutine
    
    
    SubRoutine Flexible_Force
    Use FlexibleOnly
    USE VARIABLE,ONLY: NPBE,PBE,NSBE,SBE,IPE,PN,X,Y,NODEELE,U,V,WEI,T,&
                    XPDKSI11,YPDITA22,XPDITA21,YPDKSI12,FPDX,FPDKSI,FPDITA,FPDY,TOR,M1,M2,M3,M4,&
                    DT,VEL0_ANG, VEL_ANG,&
                    WORK,WORK_TOR,WORK_F_Y,WORK_F_X,CL,CD,TOR0,CL0,CD0,&
                    VEL0_X,VEL_X,VEL0_Y,VEL_Y,DIS_X,DIS0_X,DIS_Y,DIS0_Y,DIS_ANG,DIS0_ANG,&
                    FORCE_X_CENTER_0,&
                    CL_CUT,ETA_R,ETA_T,ETA_L,X_ALE,Y_ALE,Nsbn_ORIG
    USE CONSTANT,ONLY: VISCOSITY,DENSITY,UUXX00,IROTA,LENGTH_LMM,MASSR,LL0

    Implicit None
    !!!Local!!!
	INTEGER:: I,K,IE,NL,NK,IP,NODE1,NODE2,NET(4)
	REAL(8):: P1,P2,PM,X1,X2,Y1,Y2,DX,DY,DS,FNX,FNY,FTX,FTY,FD,Q,DVDX,DUDY,UT,VT,TAU,PP00,TT00
	REAL(8):: PX,PY,SX,SY,PX1,PY1,SX1,SY1,CD1,CL1,CD2,CL2,TORP,TORS,TOR1,PI    
    REAL(8):: TOR_PIPE_TAU,TOR_PIPE_SIG
    
    INTEGER::ElementIndex,J,NKINDSN
    REAL(8)::Flo_S(4),Fe(4,FL_Params_ne),Q_ORIG,Flo_xi
    Real(8)::PFlexi(NsbN_ORIG,3),XFlexi_ORIG,YFlexi_ORIG,xFtemp,xxFtemp
    Real(8),Allocatable :: PFlexi_1D(:,:)

    !Real(8),Allocatable::
    !!!Local!!!
    
    
    PI=DACOS(-1.0D0)

    !这两个无量纲参数可能需要修改
    !PP00=1/2*Rho_w*U^2*D (U:characteristical Velocity D:characteristical Length)
    PP00=(0.5D0*DENSITY*UUXX00**2.0D0)
    !TT00=rho_w*J_total/rho_Structure*U^2/2/(Aera_Cylinder+Area_Plate)

    ! This is Dimensionless coeficience for Circular Cylinder With Plate
    !TT00=((PI/64.0D0*LENGTH_LMM**2.0D0+0.13D0/12.0D0*LENGTH_LMM**2.0D0)/(PI/4.0D0+0.02D0)*DENSITY*UUXX00**2.0D0)
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

        PM=0.5D0*(P1+P2)*DS
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
    

    PX1 =   0.0D0
    PY1 =   0.0D0
    SX1 =   0.0D0
    SY1 =   0.0D0

    PFlexi(Nsbn_ORIG,3) = 0.0D0
    
    TORP    =   0.0D0
    TORS    =   0.0D0
    TOR     =   0.0D0

    
    DO I=1,NSBE
        IE=SBE(I,1)
        NL=SBE(I,2)
        IP=IPE(IE)

        IF (IP/=4) THEN
            WRITE(*,*)' FDVG10'
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

        XFlexi_ORIG=X_ALE(NODE1)
        YFlexi_ORIG=Y_ALE(NODE1)
        
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
        SX1=SX1+(-TAU*FTX)
        SY1=SY1+(-TAU*FTY)

        
        PFlexi(I,1)=XFlexi_ORIG
        PFlexi(I,2)=(PX1+SX1)/PP00
        PFlexi(I,3)=(PY1+SY1)/PP00
        
        IF (I==NSBE) Then
            PFlexi(I+1,1)=X_ALE(NODE2)
            PFlexi(I+1,2)=(PX1+SX1)/PP00
            PFlexi(I+1,3)=(PY1+SY1)/PP00
        End if
        
        
        !write(*,*) PFlexi(I,3)
        !JY!尾板切应力对刚（圆）心的力矩,在柔性程序里面是不需要的
        !TORS=TORS-(-TAU*FTX)*(0.5D0*(Y1+Y2)-DIS0_Y)+(-TAU*FTY)*(0.5D0*(X1+X2)-DIS0_X)         !板上切向粘性剪切力对圆心的力矩

    END DO
    
    
    !这里主要是找出到底有多少个x坐标的种类，将相同x坐标的力进行相加
    DO I =1,NSBN_ORIG
        
            xFtemp=PFlexi(I,1)
            Do J = I+1,NSBN_ORIG
                xxFtemp =PFlexi(J,1)
                IF (xFtemp==xxFtemp) THEN
                    PFlexi(J,1)=0
                    PFlexi(I,2)=PFlexi(I,2)+PFlexi(J,2)
                    PFlexi(I,3)=PFlexi(I,3)+PFlexi(J,3)
                END IF
            End DO
        End do
        
        NKindSN=0
        DO I =1,NSBN_ORIG
            xFtemp=PFlexi(I,1)
            IF(xFtemp /= 0) then
                NKindSN =  NKindSN + 1
            End if
        End do
        
        Allocate(PFlexi_1D(NKindSN,3))
        PFlexi_1D = 0.0D0
        J = 1
        DO I = 1, NSBN_ORIG
            IF (PFlexi(I,1)/=0)THEN
                PFlexi_1D(J,1)=PFlexi(I,1)
                PFlexi_1D(J,2)=PFlexi(I,2)
                PFlexi_1D(J,3)=PFlexi(I,3)
                J=J+1
            END IF
        END DO
        
        Call Sort_Ascend_Multi(PFlexi_1D,NKindSN,3)

        !do i = 1, NKindSN
        !write(*,'(*(f12.6,3X))') PFlexi_1D(i,:)
        !End do
        
        ElementIndex=0.0d0
        FEULERBER = 0.0d0

        Fe = 0.0d0
        
        DO i = 1, NKindSN
            ElementIndex=Floor((PFlexi_1D(i,1)-0.5d0)/FL_Params_L)

            IF (ElementIndex<FL_Params_ne)THEN
                Flo_S=0
                Flo_xi = (PFlexi_1D(i,1)-0.5d0-(ElementIndex)*FL_Params_L)/FL_Params_L
                Call Shape_Fun(Flo_S,4,Flo_xi,0)

                !write(*,*) Flo_S,Flo_xi
                !PAUSE

                Do j = 1,4
                    Fe(j,ElementIndex+1)=Fe(j,ElementIndex+1)+PFlexi_1D(i,3)*Flo_S(j)
                End do
                Fe = Fe * 0.5D0
                
                FEULERBER(2*(ElementIndex+1)-1:2*(ElementIndex+1)+2)=FEULERBER(2*(ElementIndex+1)-1:2*(ElementIndex+1)+2)+Fe(1:4,ElementIndex+1)
            ELSE
                Flo_S=0
                Flo_xi = (PFlexi_1D(i,1)-0.5d0-(ElementIndex-1)*FL_Params_L)/FL_Params_L
                !write(*,*) 'sssssssssssssssssss',Flo_xi
                
                Call Shape_Fun(Flo_S,4,Flo_xi,0)
                Do j = 1,4
                    Fe(j,ElementIndex)=Fe(j,ElementIndex)+PFlexi_1D(i,3)*Flo_S(j)
                End do
                Fe = Fe * 0.5D0
                
                FEULERBER(2*ElementIndex-1:2*ElementIndex+2)=FEULERBER(2*ElementIndex-1:2*ElementIndex+2)+Fe(1:4,ElementIndex)
            End if
        End do

        !WRITE(*,'(*(f12.4,3x))') FEULERBER
        FEULERBER = FEULERBER / 2.0 / MassR
        DEAllocate(PFlexi_1D)

    End SubRoutine
    

    SubRoutine Sort_Ascend_Multi(XXSORT,Rows,Columns)
    IMPLICIT NONE

    INTEGER::I,J,Rows,Columns
    REAL*8::XXSORT(Rows,Columns)
    REAL*8::TMP(1,Columns)

    DO I=1,Rows-1
        DO J=I+1,Rows
            IF(XXSORT(J,1) < XXSORT(I,1))THEN
                TMP(1,:)=XXSORT(I,:)
                XXSORT(I,:)=XXSORT(J,:)
                XXSORT(J,:)=TMP(1,:)
            END IF
        END DO
    END DO

    End SubRoutine
    
    SubRoutine Sort_Descend_Multi(XXSORT,Rows,Columns)
    IMPLICIT NONE

    INTEGER::I,J,Rows,Columns
    REAL*8::XXSORT(Rows,Columns)
    REAL*8::TMP(1,Columns)

    DO I=1,Rows-1
        DO J=I+1,Rows
            IF(XXSORT(J,1) > XXSORT(I,1))THEN
                TMP(1,:)=XXSORT(I,:)
                XXSORT(I,:)=XXSORT(J,:)
                XXSORT(J,:)=TMP(1,:)
            END IF
        END DO
    END DO

    End SubRoutine
    
    SubRoutine Sort_Ascend_Multi_V2(XXSORT,Rows,Columns,PColums)
    IMPLICIT NONE

    INTEGER::I,J,Rows,Columns,PColums
    REAL*8::XXSORT(Rows,Columns)
    REAL*8::TMP(1,Columns)

    DO I=1,Rows-1
        DO J=I+1,Rows
            IF(XXSORT(J,PColums) < XXSORT(I,PColums))THEN
                TMP(1,:)=XXSORT(I,:)
                XXSORT(I,:)=XXSORT(J,:)
                XXSORT(J,:)=TMP(1,:)
            END IF
        END DO
    END DO

    End SubRoutine
    

    subroutine matrixinv(a,b,n)
    ! subroutine to calculate the inverse of a matrix using Gauss-Jordan elimination
    ! the inverse of matrix a(n,n) is calculated and stored in the matrix b(n,n)
    integer :: i,j,k,l,m,n,irow
    real(8):: big,a(n,n),b(n,n),dum

    !build the identity matrix
    do i = 1,n
        do j = 1,n
            b(i,j) = 0.0
        end do
        b(i,i) = 1.0
    end do

    do i = 1,n ! this is the big loop over all the columns of a(n,n)
        ! in case the entry a(i,i) is zero, we need to find a good pivot; this pivot
        ! is chosen as the largest value on the column i from a(j,i) with j = 1,n
        big = a(i,i)
        do j = i,n
            if (a(j,i).gt.big) then
                big = a(j,i)
                irow = j
            end if
        end do
        ! interchange lines i with irow for both a() and b() matrices
        if (big.gt.a(i,i)) then
            do k = 1,n
                dum = a(i,k) ! matrix a()
                a(i,k) = a(irow,k)
                a(irow,k) = dum
                dum = b(i,k) ! matrix b()
                b(i,k) = b(irow,k)
                b(irow,k) = dum
            end do
        end if
        ! divide all entries in line i from a(i,j) by the value a(i,i);
        ! same operation for the identity matrix
        dum = a(i,i)
        do j = 1,n
            a(i,j) = a(i,j)/dum
            b(i,j) = b(i,j)/dum
        end do
        ! make zero all entries in the column a(j,i); same operation for indent()
        do j = i+1,n
            dum = a(j,i)
            do k = 1,n
                a(j,k) = a(j,k) - dum*a(i,k)
                b(j,k) = b(j,k) - dum*b(i,k)
            end do
        end do
    end do

    ! substract appropiate multiple of row j from row j-1
    do i = 1,n-1
        do j = i+1,n
            dum = a(i,j)
            do l = 1,n
                a(i,l) = a(i,l)-dum*a(j,l)
                b(i,l) = b(i,l)-dum*b(j,l)
            end do
        end do
    end do

    end subroutine
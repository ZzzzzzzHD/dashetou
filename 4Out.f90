
    !
    !
    SUBROUTINE OUTPUT
    USE VARIABLE
    USE CONSTANT

    CHARACTER(15):: TITLE

    WRITE(FILENAME,'(I7.6,A)') NSTEP,'.dat'
    !		OPEN(UNIT=10,FILE='F:\RESULT\'//FILENAME)
    !		OPEN(UNIT=10,FILE=FILENAME, FORM='UNFORMATTED')
    OPEN(UNIT=10,FILE=FILENAME)

    WRITE(TITLE,'(F10.5,A)') T,' s'
    WRITE(10,*) 'TITLE=" ' // TITLE// ' " '
    WRITE(10,*) 'VARIABLES="X","Y","U","V" ,"P"  , "W" '
    WRITE(10,'(1X,A7,I6,A4,I6,A31)')  'ZONE N=',NN, ', E=',NE,',F=FEPOINT, ET=QUADRILATERAL'
    WRITE(10,*) 'DT=(SINGLE SINGLE SINGLE SINGLE SINGLE  SINGLE)'

    !		allocate(r8a(nn))
    !		r8a=0.0d0
    !		call ffnodes
    !		CALL vortices

    CALL VOR2

    DO I=1,NN
        WRITE(10,*)	X(I),Y(I),U(I),V(I),PN1(I) ,VOR(I)
    END DO
    DO I=1,NE
        WRITE(10,*) (NODEELE(I,J),J=1,4)
    END DO
    CLOSE(10)

    !		DEALLOCATE(r8a)

    END SUBROUTINE OUTPUT
    !
    !
    !
    !
    !
    SUBROUTINE vortices
    USE variable, only: nn,ne,x,y,nodeele,i4array,vor,u0,v0,XPDKSI11,XPDITA21,YPDKSI12,YPDITA22,FPDX,FPDY,FPDKSI,FPDITA

    IMPLICIT NONE

    INTEGER:: JJ,NET(4),K,NK,I,J
    REAL(8):: FD,Q,QR,UU,VV,DUDY,DVDX,ww
    REAL(8):: VT,UT


    integer:: me,mn,node
    allocate(i4array(nn))
    i4array=0
    vor=0.0d0

    do me=1,ne
        net(:)=nodeele(me,:)

        do mn=1,4
            node=net(mn)

            IF (mn==1) THEN
                i=4; j=4
            ELSE IF(mn==2) THEN
                i=3; j=4
            ELSE IF(mn==3) THEN
                i=3; j=3
            ELSE
                i=4; j=3
            END IF

            CALL GAUSS_JACOB(I,J,me)
            FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12


            UU  =0.0d0;		DUDY=0.0d0
            VV  =0.0d0;		DVDX=0.0d0
            DO K=1,4
                FPDX(K)=            (YPDITA22*FPDKSI(I,J,K)-YPDKSI12*FPDITA(I,J,K))/FD
                FPDY(K)=            (XPDKSI11*FPDITA(I,J,K)-XPDITA21*FPDKSI(I,J,K))/FD

                NK=NET(K)
                UT=U0(NK)
                VT=V0(NK)
                DUDY=DUDY+FPDY(K)*UT
                DVDX=DVDX+FPDX(K)*VT
            END DO
            ww=dvdx-dudy

            vor(node)=vor(node)+ww
            i4array(node)=i4array(node)+1

        end do
    end do

    do i=1,nn
        j=i4array(i)
        vor(i)=vor(i)/float(j)
    end do

    deallocate(i4array)

    END SUBROUTINE vortices
    !
    !
    !
    !
    !
    SUBROUTINE VOR2
    USE VARIABLE,ONLY:  NODEELE,NE,XPDKSI11,YPDITA22,XPDITA21,YPDKSI12,FPDX,FPDKSI,FPDITA,FPDY,U0,V0,&
        VOR,SELE,maxnum_node_ele,NODE_ELE,NN,vore
    IMPLICIT NONE

    INTEGER:: IE,NET(4),K,NK,I,NUM,M
    REAL(8):: FD,DUDY,DVDX,VT,UT,SALL

    allocate(vore(ne))

    Do IE=1,Ne

        NET=NODEELE(IE,:)
        CALL GAUSS_JACOB(5,5,IE)
        FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12

        DUDY =0.0d0;	DVDX =0.0d0

        DO K=1,4
            FPDX(K)=(YPDITA22*FPDKSI(5,5,K)-YPDKSI12*FPDITA(5,5,K))/FD
            FPDY(K)=(XPDKSI11*FPDITA(5,5,K)-XPDITA21*FPDKSI(5,5,K))/FD

            NK=NET(K)
            UT=U0(NK)
            VT=V0(NK)

            DUDY=DUDY+FPDY(K)*UT
            DVDX=DVDX+FPDX(K)*VT;
        END DO

        VORE(IE)=DVDX-DUDY
    END DO

    VOR=0.0D0

    DO I=1,NN
        NUM=0
        SALL=0.0D0

        DO  K=1,maxnum_node_ele
            M=NODE_ELE(I,K)

            IF (M.NE.0) THEN
                SALL=SALL+SELE(M)
                VOR(I)=VOR(I)+VORE(M)*SELE(M)
                NUM=NUM+1
            ELSE
                EXIT
            END IF

            IF(NUM==0) THEN
                PRINT*, 'HJP99'
                PAUSE
                STOP
            END IF

        END DO
        VOR(I)=VOR(I)/SALL

    END DO
    deallocate(vore)

    END SUBROUTINE VOR2
    !
    !
    !
    !从圆角方柱那里来的力计算程序



    SUBROUTINE FORCE

    USE VARIABLE,ONLY: NPBE,PBE,NSBE,SBE,IPE,PN,X,Y,NODEELE,U,V,WEI,T,&
        XPDKSI11,YPDITA22,XPDITA21,YPDKSI12,FPDX,FPDKSI,FPDITA,FPDY,TOR,M1,M2,M3,M4,&
        DT,VEL0_ANG, VEL_ANG,&
        WORK,WORK_TOR,WORK_F_Y,WORK_F_X,CL,CD,TOR0,CL0,CD0,CD1,CL1,CD2,CL2,&
        VEL0_X,VEL_X,VEL0_Y,VEL_Y,DIS_X,DIS0_X,DIS_Y,DIS0_Y,DIS_ANG,DIS0_ANG,&
        FORCE_X_CENTER_0,&
        CL_CUT,ETA_R,ETA_T,ETA_L
    USE CONSTANT,ONLY: VISCOSITY,DENSITY,UUXX00,IROTA,LENGTH_LMM,MASSR!,LL0

    IMPLICIT NONE
    INTEGER:: I,K,IE,NL,NK,IP,NODE1,NODE2,NET(4)
    REAL(8):: P1,P2,PM,X1,X2,Y1,Y2,DX,DY,DS,FNX,FNY,FTX,FTY,FD,Q,DVDX,DUDY,UT,VT,TAU,PP00
    REAL(8):: PX,PY,SX,SY,PX1,PY1,SX1,SY1,TORP,TORS,TOR1,PI            !!!!
    !JY!REAL(8)::CL,CD
    Real(8):: LL0                   !!!!LL0=0.0

    !JY!=========================================================================================================================
    REAL(8):: TOR_PIPE_TAU,TOR_PIPE_SIG!JY!圆柱上由切应力引起的转矩
    REAL(8):: DW_TOR,DW_F_Y,DW_F_X!JY!DW微功，WORK累计做功

    !!JY!!============计算气动中心需要注意=============================!!JY!!
    !!!
    !由于气动中心可能会在原点徘徊，因此直接计算的话会出现0/0的现象，计算机处理这种类型的计算会出现问题，比如很高的毛刺。
    !因此使用加权平均的方法来求解，即X=(CL1*X1+CL2*X2+...+CLN*XN)/CL的方法行不通,当CL接近0时，暂时强迫XC=XC_OLD，其余后期处理时再考虑

    !!JY!!============计算气动中心需要注意=============================!!JY!!

    REAL(8):: FORCE_X_CENTER,FORCE_Y_CENTER!JY!绝对坐标!X方向气动中心，Y方向气动中心
    REAL(8):: FORCE_X_CENTER_RELA,FORCE_Y_CENTER_RELA!JY!相对坐标!X方向气动中心，Y方向气动中心
    REAL(8):: DIS_K_M!JY!刚心与质心间距，这里暂时定为常值 为2/(25PI+2)
    REAL(8):: DIS_F_K!JY!气动中心到刚心的距离
    REAL(8):: DIS_ANG_F_K_RELA!JY!气动中心与刚心连线与尾板方向之间的相对角度
    REAL(8):: TT00!JY!对应PP00，弯矩无量纲化用的系数
    !JY!=========================================================================================================================

    LL0 = 0.0D0

    PI=DACOS(-1.0D0)
    DIS_K_M=2.0D0/(25.0D0*PI+2.0D0)

    PX=0.0D0
    PY=0.0D0
    SX=0.0D0
    SY=0.0D0


    TOR_PIPE_TAU=0.0D0
    TOR_PIPE_SIG=0.0D0

    FORCE_X_CENTER=0.0D0
    FORCE_Y_CENTER=0.0D0


    DO I=1,NPBE
        IE=PBE(I,1)          !单元
        NL=PBE(I,2)          !单元节点  1-4
        IP=IPE(IE)           !四边形

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
        PX=PX+PM*FNX
        PY=PY+PM*FNY

        TOR_PIPE_SIG=TOR_PIPE_SIG-PM*FNX*(0.5D0*(Y1+Y2)-DIS0_Y)+PM*FNY*(0.5D0*(X1+X2)-DIS0_X)
        !==========================================================

        CALL TXTY(X1,Y1,X2,Y2,FTX,FTY)

        CALL GAUSS_JACOB(5,5,IE)
        FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
        Q=FD*WEI(3)*WEI(3)
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
        SX=SX+(-TAU*FTX)
        SY=SY+(-TAU*FTY)

        !JY!圆柱切应力对刚（圆）心的力矩=====
        TOR_PIPE_TAU=TOR_PIPE_TAU-(-TAU*FTX)*(0.5D0*(Y1+Y2)-DIS0_Y)+(-TAU*FTY)*(0.5D0*(X1+X2)-DIS0_X)         !板上切向粘性剪切力对圆心的力矩

        !JY!=================================

        !!JY!计算气动中心======加权求解，无法进行
        FORCE_X_CENTER=FORCE_X_CENTER+(PM*FNY-TAU*FTY)*0.5D0*(X1+X2)
        FORCE_Y_CENTER=FORCE_Y_CENTER+(PM*FNX-TAU*FTX)*0.5D0*(Y1+Y2)
        !!JY!计算气动中心******

    END DO
    !*
    !*

    PX1=0.0D0
    PY1=0.0D0
    SX1=0.0D0
    SY1=0.0D0

    TORP=0.0D0
    TORS=0.0D0
    TOR =0.0D0

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

        DX=X2-X1
        DY=Y2-Y1
        DS=DSQRT(DX*DX+DY*DY)

        CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)

        PM=0.5D0*(P1+P2)*DS
        PX1=PX1+PM*FNX
        PY1=PY1+PM*FNY
        !==========================================================
        !JY!尾板压力对刚（圆）心的力矩=====
        TORP=TORP-PM*FNX*(0.5D0*(Y1+Y2)-DIS0_Y)+PM*FNY*(0.5D0*(X1+X2)-DIS0_X)              !板上法向压力对圆心的力矩  !JY!认为圆心是刚心

        !JY!暂留，加入刚心概念===========
        !!!!!    TORP=TORP-PM*FNX*(0.5D0*(Y1+Y2)-K_Y-DIS0_Y)+PM*FNY*(0.5D0*(X1+X2)-K_X-DIS0_X)  !JY!  (K_X,K_Y)刚心坐标
        !JY!=============================


        !==========================================================

        CALL TXTY(X1,Y1,X2,Y2,FTX,FTY)

        CALL GAUSS_JACOB(5,5,IE)
        FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
        Q=FD*WEI(3)*WEI(3)
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
        !============================================================
        !JY!尾板切应力对刚（圆）心的力矩=====
        TORS=TORS-(-TAU*FTX)*(0.5D0*(Y1+Y2)-DIS0_Y)+(-TAU*FTY)*(0.5D0*(X1+X2)-DIS0_X)         !板上切向粘性剪切力对圆心的力矩

        !JY!=================================
        !============================================================

        !!JY!计算气动中心======
        FORCE_X_CENTER=FORCE_X_CENTER+(PM*FNY-TAU*FTY)*0.5D0*(X1+X2)
        FORCE_Y_CENTER=FORCE_Y_CENTER+(PM*FNX-TAU*FTX)*0.5D0*(Y1+Y2)
        !!JY!计算气动中心******

    END DO

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

    PX=PX/PP00
    PY=PY/PP00
    SX=SX/PP00
    SY=SY/PP00

    PX1=PX1/PP00
    PY1=PY1/PP00
    SX1=SX1/PP00
    SY1=SY1/PP00


    CD=(PX+SX+PX1+SX1)            !系统受力
    CL=(PY+SY+PY1+SY1)

    CD1=(PX+SX)                   !圆柱受力
    CL1=(PY+SY)

    CD2=(PX1+SX1)                 !隔板受力
    CL2=(PY1+SY1)

    !JY!松弛尝试==========
    CL=LL0*CL0+(1-LL0)*CL
    CD=LL0*CD0+(1-LL0)*CD
    !=====================

    WRITE(999,'(7E18.8)') T,CD,CL,CD1,CL1,CD2,CL2
    !=============================================================!为了得到与平动结构运动方程相同的无量纲形式，对力矩进行了无量纲化

    TOR1=TORP+TORS+TOR_PIPE_TAU+TOR_PIPE_SIG
    TOR=TOR1/TT00

    TOR=LL0*TOR0+(1-LL0)*TOR
    !WRITE(777,'(6E18.8)') T,TOR,TOR1,TORP,TORS,TOR_PIPE_TAU  !zhang cheng 2023 12 06



    !JY!计算气动中心绝对坐标====
    IF(DABS(CL)<CL_CUT)THEN
        FORCE_X_CENTER=FORCE_X_CENTER_0!JY!CL会掠过零点，因此只能暂时强迫其为上一个值
    ELSE
        FORCE_X_CENTER=FORCE_X_CENTER/(CL*PP00)
        !JY!松弛一下
        FORCE_X_CENTER=LL0*FORCE_X_CENTER_0+(1-LL0)*FORCE_X_CENTER
    END IF



    FORCE_Y_CENTER=FORCE_Y_CENTER/(CD*PP00)

    !JY!计算气动中心************

    !!JY!计算气动中心相对（圆柱中心为原点，且结构转角为零）坐标
    !DIS_F_K=DSQRT((FORCE_X_CENTER-DIS_X)**2+(FORCE_Y_CENTER-DIS_Y)**2)
    !
    !IF(DIS_F_K==0)THEN
    !    FORCE_X_CENTER_RELA=0
    !    FORCE_Y_CENTER_RELA=0
    !ELSE
    !    IF(FORCE_Y_CENTER-DIS_Y<0)THEN
    !        DIS_ANG_F_K_RELA=-DACOS((FORCE_X_CENTER-DIS_X)/DIS_F_K)-DIS_ANG
    !    ELSE
    !        DIS_ANG_F_K_RELA=DACOS((FORCE_X_CENTER-DIS_X)/DIS_F_K)-DIS_ANG
    !    END IF
    !
    !    FORCE_X_CENTER_RELA=DIS_F_K*DCOS(DIS_ANG_F_K_RELA)
    !    FORCE_Y_CENTER_RELA=DIS_F_K*DSIN(DIS_ANG_F_K_RELA)
    !END IF

    !JY!计算气动中心相对（圆柱中心为原点，且结构转角为即时转角）坐标
    FORCE_X_CENTER_RELA=FORCE_X_CENTER-DIS_X
    FORCE_Y_CENTER_RELA=FORCE_Y_CENTER-DIS_Y


    WRITE(2331,'(5E18.8)') T,FORCE_X_CENTER,FORCE_Y_CENTER
    WRITE(2332,'(5E18.8)') T,FORCE_X_CENTER_RELA,FORCE_Y_CENTER_RELA

    FORCE_X_CENTER_0=FORCE_X_CENTER
    !JY!计算气动中心************


    !JY!做功输出
    DW_TOR=0.5D0*DT*(TOR0*VEL0_ANG+TOR*VEL_ANG)*TT00
    DW_F_Y=0.5D0*DT*(CL0*(VEL0_Y+DIS_K_M*VEL0_ANG*DCOS(DIS0_ANG))+CL*(VEL_Y+DIS_K_M*VEL_ANG*DCOS(DIS_ANG)))*PP00!JY!因为刚心与中心不重合
    DW_F_X=0.5D0*DT*(CD0*(VEL0_X-DIS_K_M*VEL0_ANG*DSIN(DIS0_ANG))+CD*(VEL_X-DIS_K_M*VEL_ANG*DSIN(DIS_ANG)))*PP00!JY!因为刚心与中心不重合
    WORK_TOR=WORK_TOR+DW_TOR
    WORK_F_Y=WORK_F_Y+DW_F_Y
    WORK_F_X=WORK_F_X+DW_F_X
    WORK=WORK+DW_TOR+DW_F_Y+DW_F_X
    WRITE(1025,'(5E18.8)') T-0.5D0*DT,WORK,DW_TOR+DW_F_Y+DW_F_X
    WRITE(1026,'(5E18.8)') T-0.5D0*DT,WORK_TOR,DW_TOR
    WRITE(1027,'(5E18.8)') T-0.5D0*DT,WORK_F_Y,DW_F_Y
    WRITE(1028,'(5E18.8)') T-0.5D0*DT,WORK_F_X,DW_F_X

    !JY!=================================================================================================================================================

    END SUBROUTINE FORCE


    !JY!查看ale右边界的速度
    SUBROUTINE OUTPUT_ale_right
    USE VARIABLE
    USE CONSTANT

    CHARACTER(15):: TITLE
    integer j

    WRITE(FILENAME,*) 'ale_right_',NSTEP,'.dat'

    OPEN(UNIT=10,FILE=FILENAME)

    WRITE(TITLE,'(F10.5,A)') T,' s'
    WRITE(10,*) 'TITLE=" ' // TITLE// ' " '
    WRITE(10,*) 'VARIABLES="Y","U","V","X"'

    DO I=1,Nrbn
        j=rbN_ALE(I)
        WRITE(10,*) Y(J),U(J),V(J),X(J)
    END Do

    CLOSE(10)


    END SUBROUTINE OUTPUT_ale_right





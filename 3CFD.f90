    SUBROUTINE C_ISOPE
    USE VARIABLE

    IMPLICIT NONE
    INTEGER:: I,J
    REAL(8):: KSI,ITA

    XI(1) = -DSQRT(3.0D0)/3.0D0
    XI(2) = -XI(1)
    XI(3) =  1.0D0
    XI(4) = -1.0D0
    XI(5) =  0.0D0

    WEI(1)=  1.0D0
    WEI(2)=  1.0D0
    WEI(3)=  2.0D0

    DO I=1,5
        KSI=XI(I)
        DO J=1,5
            ITA=XI(J)
            F(I,J,1)=0.25D0*(1.0D0-KSI)*(1.0D0-ITA)
            F(I,J,2)=0.25D0*(1.0D0+KSI)*(1.0D0-ITA)
            F(I,J,3)=0.25D0*(1.0D0+KSI)*(1.0D0+ITA)
            F(I,J,4)=0.25D0*(1.0D0-KSI)*(1.0D0+ITA)

            FPDKSI(I,J,1)=-0.25D0*(1.0D0-ITA)
            FPDKSI(I,J,2)= 0.25D0*(1.0D0-ITA)
            FPDKSI(I,J,3)= 0.25D0*(1.0D0+ITA)
            FPDKSI(I,J,4)=-0.25D0*(1.0D0+ITA)

            FPDITA(I,J,1)=-0.25D0*(1.0D0-KSI)
            FPDITA(I,J,2)=-0.25D0*(1.0D0+KSI)
            FPDITA(I,J,3)= 0.25D0*(1.0D0+KSI)
            FPDITA(I,J,4)= 0.25D0*(1.0D0-KSI)
        END DO
    END DO

    END SUBROUTINE C_ISOPE


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    SUBROUTINE GAUSS_JACOB(I,J,ELENUM)
    USE VARIABLE
    USE CONSTANT

    IMPLICIT NONE
    INTEGER:: I,J,K,L,NET(4),ELENUM
    REAL(8):: FPDKSIK,FPDITAK,XNEW(4),YNEW(4)

    XPDKSI11=0.0D0
    XPDITA21=0.0D0
    YPDKSI12=0.0D0
    YPDITA22=0.0D0

    DO K=1,4
        NET(K)=NODEELE( ELENUM,K )
    END DO

    DO K=1,4
        L=NET(K)
        XNEW(K)=X(L)
        YNEW(K)=Y(L)
    END DO

    DO K=1,4
        FPDKSIK=FPDKSI(I,J,K)
        FPDITAK=FPDITA(I,J,K)
        XPDKSI11=XPDKSI11+FPDKSI(I,J,K)*XNEW(K)
        XPDITA21=XPDITA21+FPDITA(I,J,K)*XNEW(K)
        YPDKSI12=YPDKSI12+FPDKSI(I,J,K)*YNEW(K)
        YPDITA22=YPDITA22+FPDITA(I,J,K)*YNEW(K)
    END DO

    END SUBROUTINE GAUSS_JACOB



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    SUBROUTINE INDEX00
    USE VARIABLE
    USE CONSTANT
    IMPLICIT NONE

    INTEGER:: I,J,M,K,MM,NK,L,L1,NODEOLD,NODE,NAA,IG,IP
    INTEGER:: IN,MI,II,JEMP,NNN,ERR,NET(4),NQ,KK
    INTEGER,ALLOCATABLE:: JMEN(:), NEN(:,:), J2(:),NENQ(:,:)


    ALLOCATE(JMEN(NN),J2(NN),IA(NN+1),I4ARRAY(NN*20))
    ALLOCATE(NENQ(10,NN))

    N1DA=0
    NNN =0
    JMEN=0
    NEN =0

    DO I=1,NE
        NET(:)=NODEELE(I,:)  !结点
        DO J=1,4
            MM=NET(J)
            JMEN(MM)=JMEN(MM)+1      !MM节点出现的个数=NUMBER OF SURROUNDING ELEMENTS
            NK=JMEN(MM)
            IF(NK>10) THEN
                WRITE(*,*)' EM9'
                WRITE(*,*) NK,MM,X(MM),Y(MM)
                PAUSE
                STOP
            END IF
            NENQ(NK,MM)=I	         !第NK个节点MM属于单元I
        END DO
    END DO

    IA(1)=1
    DO I=1,NN
        IG=0
        NQ=JMEN(I)                  ! NUMBER OF SURROUNDING ELEMENTS  1   2   4
        DO J=1,NQ
            NK=NENQ(J,I)            ! ELEMENT INDEX单元
            NET(:)=NODEELE(NK,:)
            DO K=1,4
                L=NET(K)            !结点
                IF (L<=I) THEN
                    IG=IG+1
                    !				IF(IG>=20) THEN; WRITE(*,*) 'IG20'; STOP; END IF
                    J2(IG)=L
                END IF
            END DO
        END DO


        DO M=1,IG-1
            DO L=M+1,IG
                IF(J2(M)==J2(L)) THEN
                    J2(M)=0    !去掉重复
                END IF
            END DO
        END DO

        IN=0
        DO J=1,IG
            IF(J2(J).NE.0) THEN
                IN=IN+1
                J2(IN)=J2(J)
            END IF
        END DO

        M=IN-1
        DO K=1,M
            DO KK=K+1,IN
                IF( J2(K) .GE. J2(KK) ) THEN
                    JEMP=J2(K)
                    J2(K)=J2(KK)
                    J2(KK)=JEMP       !按顺序从小到大排列
                END IF
            END DO
        END DO

        IA(I+1)=IA(I)+IN          !--------XHL---------IA:节点周围不大于自身的节点个数
        MI=IA(I)-1

        DO II=1,IN
            NAA=MI+II
            N1DA=N1DA+1
            I4ARRAY(NAA)=J2(II)
            IF(NAA/=N1DA) THEN
                WRITE(*,*)' RC5'
                PAUSE
                STOP
            END IF
        END DO
    END DO

    IF(NAA+1/=IA(NN+1)) THEN      !--------XHL---------？
        WRITE(*,*) NAA,N1DA,IA(NN+1)
        WRITE(*,*) ' BFG'
        PAUSE
        STOP
    END IF

    ALLOCATE(JA(N1DA))
    JA(1:N1DA)=I4ARRAY(1:N1DA)
    DEALLOCATE(I4ARRAY)

    WRITE(*,'(A10,I12)')'1D ARRAY =',N1DA
    END SUBROUTINE INDEX00
    !
    !
    SUBROUTINE Clock2AntiClock_Total
    USE VARIABLE
    USE CONSTANT


    INTEGER:: I,J,K,NODE1,NODE2,IP,K0
    REAL(8):: S,XX2,YY2,XX1,YY1
    REAL(8):: NODEELE_INVERSE(4)


    DO I=1,NE_ORIG
        S=0.0D0
        K0=3
        DO J=1,4
            IF (J<=K0) THEN
                NODE1=NODEELE_ORIG(I,J)
                NODE2=NODEELE_ORIG(I,J+1)
            ELSE
                NODE1=NODEELE_ORIG(I,4)
                NODE2=NODEELE_ORIG(I,1)
            END IF
            XX1=X_ORIG(NODE1)
            YY1=Y_ORIG(NODE1)
            XX2=X_ORIG(NODE2)
            YY2=Y_ORIG(NODE2)
            S=S+(XX1*YY2-XX2*YY1)!
        END DO

        IF (S>0.0D0) THEN

        ELSE
            !JY!将顺时针改为逆时针================================

            S=0
            DO III=1,4
                NODEELE_INVERSE(III)=NODEELE_ORIG(I,III)
            END DO
            NODEELE_ORIG(I,4)=NODEELE_INVERSE(1)
            NODEELE_ORIG(I,3)=NODEELE_INVERSE(2)
            NODEELE_ORIG(I,2)=NODEELE_INVERSE(3)
            NODEELE_ORIG(I,1)=NODEELE_INVERSE(4)

            DO J=1,4
                IF (J<=K0) THEN
                    NODE1=NODEELE_ORIG(I,J)
                    NODE2=NODEELE_ORIG(I,J+1)
                ELSE
                    NODE1=NODEELE_ORIG(I,4)
                    NODE2=NODEELE_ORIG(I,1)
                END IF
                XX1=X_ORIG(NODE1)
                YY1=Y_ORIG(NODE1)
                XX2=X_ORIG(NODE2)
                YY2=Y_ORIG(NODE2)
                S=S+(XX1*YY2-XX2*YY1)
            END DO
            !===========================================================

            IF (S>0.0D0) THEN

            ELSE
                WRITE(*,*) ' PMCCY'
                WRITE(*,*) ' THE ELEMENT IS:', I

                WRITE(*,*) S
                PAUSE
            END IF

        END IF
    END DO


    END SUBROUTINE Clock2AntiClock_Total!
    !
    !
    !


    SUBROUTINE Clock2AntiClock_ALE
    USE VARIABLE
    USE CONSTANT


    INTEGER:: I,J,K,NODE1,NODE2,IP,K0
    REAL(8):: S,XX2,YY2,XX1,YY1
    REAL(8):: NODEELE_INVERSE(4)

    Allocate(SELE_Ale(Ne_ale))
    DO I=1,NE_Ale
        S=0.0D0
        K0=3
        DO J=1,4
            IF (J<=K0) THEN
                NODE1=NODEELE_Ale(I,J)
                NODE2=NODEELE_Ale(I,J+1)
            ELSE
                NODE1=NODEELE_Ale(I,4)
                NODE2=NODEELE_Ale(I,1)
            END IF
            XX1=X_Ale(NODE1)
            YY1=Y_Ale(NODE1)
            XX2=X_Ale(NODE2)
            YY2=Y_Ale(NODE2)
            S=S+(XX1*YY2-XX2*YY1)!
        END DO

        IF (S>0.0D0) THEN
            SELE_Ale(I)=S*0.5D0
        ELSE
            !JY!将顺时针改为逆时针================================

            S=0
            DO III=1,4
                NODEELE_INVERSE(III)=NODEELE_Ale(I,III)
            END DO
            NODEELE_Ale(I,4)=NODEELE_INVERSE(1)
            NODEELE_Ale(I,3)=NODEELE_INVERSE(2)
            NODEELE_Ale(I,2)=NODEELE_INVERSE(3)
            NODEELE_Ale(I,1)=NODEELE_INVERSE(4)

            DO J=1,4
                IF (J<=K0) THEN
                    NODE1=NODEELE_Ale(I,J)
                    NODE2=NODEELE_Ale(I,J+1)
                ELSE
                    NODE1=NODEELE_Ale(I,4)
                    NODE2=NODEELE_Ale(I,1)
                END IF
                XX1=X_Ale(NODE1)
                YY1=Y_Ale(NODE1)
                XX2=X_Ale(NODE2)
                YY2=Y_Ale(NODE2)
                S=S+(XX1*YY2-XX2*YY1)
            END DO

            IF (S>0.0D0) THEN
                SELE_Ale(I)=S*0.5D0
            ELSE
                WRITE(*,*) ' PMCCY'
                WRITE(*,*) ' THE ELEMENT IS:', I

                WRITE(*,*) S
                PAUSE
            END IF

        END IF
    END DO

    DeAllocate(SELE_Ale)
    End SubRoutine



    SUBROUTINE CHARACTERISTIC_LENGTH
    USE VARIABLE
    USE CONSTANT


    INTEGER:: I,J,K,NODE1,NODE2,IP,K0
    REAL(8):: S,XX2,YY2,XX1,YY1
    REAL(8):: NODEELE_INVERSE(4)

    IF (ALLOCATED(SELE)) DEALLOCATE(SELE)

    ALLOCATE(SELE(NE))

    DO I=1,NE
        S=0.0D0
        IP=IPE(I)
        K0=IP-1
        DO J=1,IP
            IF (J<=K0) THEN
                NODE1=NODEELE(I,J)
                NODE2=NODEELE(I,J+1)
            ELSE
                NODE1=NODEELE(I,4)
                NODE2=NODEELE(I,1)
            END IF
            XX1=X(NODE1)
            YY1=Y(NODE1)
            XX2=X(NODE2)
            YY2=Y(NODE2)
            S=S+(XX1*YY2-XX2*YY1)!
        END DO

        IF (S>0.0D0) THEN
            SELE(I)=S*0.5D0
        ELSE
            !JY!将顺时针改为逆时针================================
            S=0
            DO III=1,4
                NODEELE_INVERSE(III)=NODEELE(I,III)
            END DO
            NODEELE(I,4)=NODEELE_INVERSE(1)
            NODEELE(I,3)=NODEELE_INVERSE(2)
            NODEELE(I,2)=NODEELE_INVERSE(3)
            NODEELE(I,1)=NODEELE_INVERSE(4)

            DO J=1,IP
                IF (J<=K0) THEN
                    NODE1=NODEELE(I,J)
                    NODE2=NODEELE(I,J+1)
                ELSE
                    NODE1=NODEELE(I,4)
                    NODE2=NODEELE(I,1)
                END IF
                XX1=X(NODE1)
                YY1=Y(NODE1)
                XX2=X(NODE2)
                YY2=Y(NODE2)
                S=S+(XX1*YY2-XX2*YY1)
            END DO
            !===========================================================

            IF (S>0.0D0) THEN
                SELE(I)=S*0.5D0
            ELSE
                WRITE(*,*) ' PMCCY_Length'
                WRITE(*,*) ' THE ELEMENT IS A :', I

                Call Output
                WRITE(*,*) S
                PAUSE
            END IF

        END IF
    END DO

    IF (ALLOCATED(CLE)) DEALLOCATE(CLE)
    ALLOCATE(CLE(NE))

    CLE=DSQRT(SELE)

    END SUBROUTINE CHARACTERISTIC_LENGTH
    !
    !
    !
    !
    !
    SUBROUTINE RHME0
    USE VARIABLE
    USE CONSTANT
    IMPLICIT NONE

    INTEGER:: IE,IEOLD,IET,NET(4),K,NK,I,J
    REAL(8):: FD,Q,QR,UU,VV,DUDX,DUDY,DVDX,DVDY,PT,VT,UT,PG,DPDX,DPDY

    REAL(8):: VISE, VISM, VIST, F0, DTS ,ROU

    DTS=DT/FLOAT(NSUB)
    ROU=DENSITY
    VISM=VISCOSITY/DENSITY

    DO  IE=1,NE
        NET=NODEELE(IE,:)

        DO  I=1,2
            DO J=1,2

                CALL GAUSS_JACOB(I,J,IE)
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                Q=FD*WEI(I)*WEI(J)

                UU   =0.0D0;  	VV   =0.0D0;
                DPDX =0.0D0; 	DPDY =0.0D0;
                DUDX =0.0D0;	DUDY =0.0D0;	DVDX =0.0D0;	DVDY =0.0D0
                PG   =0.0D0;


                DO  K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(I,J,K)-YPDKSI12*FPDITA(I,J,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(I,J,K)-XPDITA21*FPDKSI(I,J,K))/FD

                    NK=NET(K)
                    UT=U(NK)
                    VT=V(NK)
                    PT=PN(NK)

                    UU=UU+ F(I,J,K)* ((UT) -UM_ALE(NK))
                    VV=VV+ F(I,J,K)* ((VT) -VM_ALE(NK))

                    PG=PG+ F(I,J,K)*  PT

                    DUDX=DUDX+FPDX(K)*UT; 		DUDY=DUDY+FPDY(K)*UT
                    DVDX=DVDX+FPDX(K)*VT; 		DVDY=DVDY+FPDY(K)*VT
                    DPDX=DPDX+FPDX(K)*PT; 		DPDY=DPDY+FPDY(K)*PT
                END DO

                !
                VIST=(CS*CLE(IE))**2.0D0* DSQRT( 2.0D0*(DUDX**2.0D0 +DVDY**2.0D0 )+(DUDY+DVDX)**2.0D0 )              !-----------XHL------------？？？？？？？？？？？？
                VISE=VIST+VISM

                DO  K=1,4
                    NK=NET(K)
                    UB(NK)= UB(NK) &
                        + (Q*DTS)* (-F(I,J,K)*DPDX/ROU   - F(I,J,K)*(UU*DUDX+VV*DUDY) - VISE*(FPDX(K)*DUDX+FPDY(K)*DUDY) )
                    !        					+ (Q*DTS)* ( FPDX(K)*PG/ROU      - F(I,J,K)*(UU*DUDX+VV*DUDY) - VISE*(FPDX(K)*DUDX+FPDY(K)*DUDY) )

                    VB(NK)= VB(NK) &
                        + (Q*DTS)* (-F(I,J,K)*DPDY/ROU   - F(I,J,K)*(UU*DVDX+VV*DVDY) - VISE*(FPDX(K)*DVDX+FPDY(K)*DVDY) )
                    !         					+ (Q*DTS)* ( FPDY(K)*PG/ROU      - F(I,J,K)*(UU*DVDX+VV*DVDY) - VISE*(FPDX(K)*DVDX+FPDY(K)*DVDY) )

                END DO

            END DO
        END DO
    END DO

    END SUBROUTINE RHME0
    !
    !
    !	=========================================================================================================================
    !
    !
    SUBROUTINE ME_SOLUTION
    USE CONSTANT
    USE VARIABLE

    IMPLICIT NONE
    INTEGER:: I,J,K

    DO I=1,NN
        U(I)=UB(I)*(1.0D0/AM(I)) +U0(I)
        V(I)=VB(I)*(1.0D0/AM(I)) +V0(I)
    END DO


    END SUBROUTINE ME_SOLUTION
    !
    !
    !
    !
    !
    SUBROUTINE IMPOSE_VELOCITY_BC
    Use FlexibleOnly
    Use ANCF
    USE CONSTANT
    USE VARIABLE
    Use Blas95
    Use Local
    IMPLICIT NONE

    INTEGER:: I,J,K, NODE,ITY,II,ElementIndex,Is,Ie
    REAL(8):: Y000 , Y111,XA,XB,YA,YB,X0,Y0,DD,PI,THETA2,R00,DTHETA,THETA1,FLO_S(4),FLO_SI(AFP%Dim,AFP%NumE),FLO_xi,V_interp(3)
    Real(8)::edele(12),W_Ome(3,1),W_Ome_Num,Index,Distancd_X ,Distancd_Y

    PI=DACOS(-1.0D0)
    U(PBN(1:NPBN))=0.0D0


    IF(IROTA==5107) THEN
        IF (T>T_START) THEN
            
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                X0=X_ALE(J)
                Y0=Y_ALE(J)

                DD=DSQRT(X0**2+Y0**2)!到原点的距离     !JY!注意是X_ALE,这个数组保持最初始的数据不动
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF
                
                !U(J)=-DD*DSIN(DIS_ANG+THETA2)*VEL_ANG
                !V(J)=DD*DCOS(DIS_ANG+THETA2)*VEL_ANG+VEL_Y
                U(J)=0.0D0
                V(J)=0.0D0

            END DO

        ELSE
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF
    ELSE IF((IROTA==510).OR.(IROTA==5101)) THEN
        DO I=1,NPBN_ALE
            J=PBN_ALE(I)
            U(J)=0.0D0
            V(J)=0.0D0
        END DO
    ELSE IF((IROTA==5102).OR.(IROTA==520)) THEN
        IF (T>T_START) THEN
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                U(J)= 0.0D0
                V(J)= AF_ed(2)
            END DO
        ELSE
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO
        END IF
    ELSE IF(IROTA==200.or.IROTA==201) THEN
        IF (T>T_START) THEN

            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                X0=X_ALE(J)
                Y0=Y_ALE(J)

                DD=DSQRT(X0**2+Y0**2)!到原点的距离     !JY!注意是X_ALE,这个数组保持最初始的数据不动
                THETA2=DASIN(Y0/DD)
                IF (X0<0.0D0)   THEN
                    THETA2=PI-THETA2
                END IF

                U(J)=-DD*DSIN(DIS_ANG+THETA2)*VEL_ANG
                V(J)= DD*DCOS(DIS_ANG+THETA2)*VEL_ANG+VEL_Y

            END DO

        ELSE
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF

    ELSE IF(IROTA==2201) THEN
        IF (T>T_START) THEN

            DO I=1,NPBN_ALE
                J=PBN_ALE(I)

                U(J)= 0.0D0
                V(J)= VEL_Y_UP

            END DO

        ELSE
            DO I=1,NPBN_ALE
                J=PBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF

    ELSE IF(IROTA==0 .OR. IROTA==2101) THEN
        DO I=1,NPBN_ALE
            J=PBN_ALE(I)
            U(J)=0.0D0
            V(J)=0.0D0
        END DO
        
    END IF

    !===============PLATE or DC=====================
    IF(IROTA==5107) THEN
        IF (TT0>T_START) THEN
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=UM_ALE(J)
                V(J)=vM_ALE(J)
            END DO
        ELSE
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO
        END IF
    ELSE IF ((IROTA==510).OR.(IROTA==5102).OR.(IROTA==5101).OR.(IROTA==520)) THEN
        IF (TT0>T_START) THEN
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                X0=X_ALE(J)
                Y0=Y_ALE(J)
                X0=X0-0.5D0
                ElementIndex=FLOOR(X0/AFP%LE)
                IF (ElementIndex<AFP%NE) THEN
                    Is = AFP%NumD*ELEMENTINDEX + 1
                    Ie = AFP%NumD*ELEMENTINDEX + AFP%NumE
                    edele = AF_ed(Is:Ie)
                    FLO_XI=(X0-ElementIndex*AFP%LE)/AFP%LE
                    FLO_S=0.0D0
                    Call ANCF_Shape_Fun(FLo_S,FLo_SI,FLo_xi,0)
                    Call S_e(FLo_S,edele,V_interp)
                    U(J)= V_interp(1)
                    V(J)= V_interp(2)
                    !
                    !K = PlatePointIndex(I,1)
                    !U(J) = U(J) + (Delta_x510(K)-Delta_x0510(K)) / DT
                    !V(J) = V(J) + (Delta_y510(K)-Delta_y0510(K)) / DT
                    
                    W_Ome(:,1) = (Sec_Rota(I,:) - Sec_Rota0(I,:)) / DT
                    IF(DAbs(Dsqrt(Sum(Sec_Rota(I,1:3)**2))) /= 0.0D0 )Then
                        W_Ome_Num = Dsqrt(W_Ome(1,1)**2.0D0 + W_Ome(2,1)**2.0D0 + W_Ome(3,1)**2.0D0) &
                            & / Dsqrt(Sum(Sec_Rota(I,1:3)**2))
                    Else
                        W_Ome_Num = 0.0D0
                    End If
                    
                    If(Sec_Rota(I,2)<0.0D0)Then
                        W_Ome_Num = -W_Ome_Num
                    End IF
                    
                    IF(T<=T_Start + 0.02)Then
                        W_Ome_Num = 0.0
                    End IF
                    
                    !Index = PlatePointIndex(I,1)
                    !Distancd_X = X_Neutal(Index,1) - X0_ALe(J)
                    !Distancd_Y = Y_Neutal(Index,1) - Y0_ALe(J)
                    !
                    !U(J) = U(J) - W_Ome_Num * Distancd_Y
                    !V(J) = V(J) + W_Ome_Num * Distancd_X

                ELSE
                    Is = AFP%NumD*(ELEMENTINDEX-1) + 1
                    Ie = AFP%NumD*(ELEMENTINDEX-1) + AFP%NumE
                    edele = AF_ed(Is:Ie)

                    FLO_XI=(X0-(ElementIndex-1)*AFP%LE)/AFP%LE!!20211214 -1
                    FLO_S=0.0D0
                    Call ANCF_Shape_Fun(FLo_S,FLo_SI,FLo_xi,0)
                    Call S_e(FLo_S,edele,V_interp)
                    U(J)= V_interp(1)
                    V(J)= V_interp(2)
                    
                    !K = PlatePointIndex(I,1)
                    !U(J) = U(J) + (Delta_x510(K)-Delta_x0510(K)) / DT
                    !V(J) = V(J) + (Delta_y510(K)-Delta_y0510(K)) / DT

                    !这一部分得弄个说明
                    W_Ome(:,1) = (Sec_Rota(I,:) - Sec_Rota0(I,:)) / DT
                    IF(DAbs(Dsqrt(Sum(Sec_Rota(I,1:3)**2))) /= 0.0D0 )Then
                        W_Ome_Num = Dsqrt(W_Ome(1,1)**2.0D0 + W_Ome(2,1)**2.0D0 + W_Ome(3,1)**2.0D0) &
                            & / Dsqrt(Sum(Sec_Rota(I,1:3)**2))
                    Else
                        W_Ome_Num = 0.0D0
                    End If

                    If(Sec_Rota(I,2)<0.0D0)Then
                        W_Ome_Num = -W_Ome_Num
                    End IF
                    
                    IF(T<=T_Start + 0.02)Then
                        W_Ome_Num = 0.0
                    End IF

                END IF
            END DO
        ELSE
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO
        END IF


    ELSE IF (IROTA==200.or.IROTA==201) THEN
        IF (TT0>T_START) THEN

            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                X0=X_ALE(J)
                Y0=Y_ALE(J)

                DD=DSQRT(X0**2+Y0**2)
                THETA2=DASIN(Y0/DD)

                U(J)=-DD*DSIN(DIS_ANG+THETA2)*VEL_ANG
                V(J)= DD*DCOS(DIS_ANG+THETA2)*VEL_ANG+VEL_Y
            END DO

        ELSE
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF
    ELSE IF (IROTA==2201) THEN
        IF (TT0>T_START) THEN

            DO I=1,NSBN_ALE
                J=SBN_ALE(I)

                U(J)= 0.0D0
                V(J)= VEL_Y
            END DO

        ELSE
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF
        
    ELSE IF (IROTA==2101) THEN
        IF (TT0>T_START) THEN

            DO I=1,NSBN_ALE
                J=SBN_ALE(I)

                U(J)= 0.0D0
                V(J)= VEL_Y
            END DO

        ELSE
            DO I=1,NSBN_ALE
                J=SBN_ALE(I)
                U(J)=0.0D0
                V(J)=0.0D0
            END DO

        END IF
        
    ELSE IF (IROTA==0) THEN

        DO I=1,NSBN_ALE
            J=SBN_ALE(I)
            U(J)=0.0D0
            V(J)=0.0D0
        END DO

    ELSE

        WRITE(*,*)"ERROR IN THE IROTA"
        PAUSE
        STOP

    END IF

    !========================================
    If(Irota==510)Then
        DO I=1,NLBN
            J=LBN(I)
            Y000 = Y(J)
            YA = MInval(Y)
            YB = MaxVal(Y)
            U(J)=1.50D0*UUXX00*1.0D0/DABS(YA * YB) *(Y000-YA)*(YB-Y000)
            V(J)=0.0D0
        END DO
        DO I=1,NDBN
            J=DBN(I)
            U(J)=0.0D0
            V(J)=0.0D0
        END DO
        DO I=1,NUBN
            J=UBN(I)
            U(J)=0.0D0
            V(J)=0.0D0
        END DO
    ELSEIF((Irota==5101).OR.(IROTA==520).OR.(IROTA==5102))Then
        DO I=1,NLBN
            J=LBN(I)
            U(J)=1.00D0
            V(J)=0.0D0
        END DO
        DO I=1,NDBN
            J=DBN(I)
            !U(J)=0.0D0
            V(J)=0.0D0
        END DO
        DO I=1,NUBN
            J=UBN(I)
            !U(J)=0.0D0
            V(J)=0.0D0
        END DO
    Else
        DO I=1,NLBN
            J=LBN(I)
            Y000 = Y(J)
            U(J)=1.00D0
            V(J)=0.0D0
        END DO
        DO I=1,NDBN
            J=DBN(I)
            !U(J)=0.0D0
            V(J)=0.0D0
        END DO
        DO I=1,NUBN
            J=UBN(I)
            !U(J)=0.0D0
            V(J)=0.0D0
        END DO
    END IF


    DO  I=1,NRBN
        NODE=RBN(I )
        Y000=Y(NODE)
        ITY=0
        DO  J=1,MAXNUM_NODE_NODE
            K=NODE_NODE(NODE,J)
            Y111=Y(K)
            IF ( (Y111==Y000).AND.(K/=NODE) ) THEN
                ITY=1
                EXIT
            END IF
        END DO

        IF  (ITY==0) THEN
            PRINT*, 'TRY56'
            PAUSE
            STOP
        ELSE
            U(NODE)=U(K)
            V(NODE)=V(K)
        END IF
    END DO


    END SUBROUTINE IMPOSE_VELOCITY_BC
    !
    !
    !
    !
    !
    !
    SUBROUTINE PRESSURE_1BC
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: I,J,K,M,KK,KU,NODE,IE,NET(4),NK,NUM,MM,KKSD
    REAL(8):: UU,Y111,Y000,ITY,VISE,VisM,VisT,ROU,FD,VT,UT,SALL,mDUDX,mDVDY,PUPX
    REAL(8),ALLOCATABLE::DUDX(:),DVDY(:),DUDXE(:),DVDYE(:),PHIE(:),PHI(:)

    ALLOCATE(DUDX(NN))
    ALLOCATE(DVDY(NN))
    ALLOCATE(DUDXE(NE))
    ALLOCATE(DVDYE(NE))
    ALLOCATE(PHIE(NE))
    ALLOCATE(PHI(NN))

    DUDXE=0.0D0
    DVDYE=0.0D0
    DUDX=0.0D0
    DVDY=0.0D0
    PHI=0.0d0
    PHIE=0.0d0

    mDUDX=0.0d0
    mDVDY=0.0d0

    ROU=DENSITY
    VisM=VISCOSITY/density

    VisT=0.0D0
    VisE=VisM+VisT
    !
    !Do IE=1,Ne
    !    NET=NODEELE(IE,:)
    !    CALL GAUSS_JACOB(5,5,IE)
    !    FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
    !
    !    mDUDX=0.0d0
    !    mDVDY=0.0d0
    !    DO K=1,4
    !        FPDX(K)=(YPDITA22*FPDKSI(5,5,K)-YPDKSI12*FPDITA(5,5,K))/FD
    !        FPDY(K)=(XPDKSI11*FPDITA(5,5,K)-XPDITA21*FPDKSI(5,5,K))/FD
    !
    !        NK=NET(K)
    !        UT=U(NK)
    !        VT=V(NK)
    !
    !        mDUDX=mDUDX+FPDX(K)*UT
    !        mDVDY=mDVDY+FPDY(K)*VT
    !    END DO
    !
    !    PHIE(IE)=mDUDX+mDVDY
    !    DUDXE(IE)=mDUDX
    !    DVDYE(IE)=mDVDY
    !END DO
    !
    !PHI=0.0d0
    !DUDX=0.0D0
    !DVDY=0.0D0
    !
    !DO I=1,NN
    !    NUM=0
    !    SALL=0.0D0
    !    DO  K=1,maxnum_node_ele
    !        MM=NODE_ELE(I,K)
    !        IF (MM.NE.0) THEN
    !            SALL=SALL+SELE(MM)
    !            PHI(I)=PHI(I)+PHIE(MM)*SELE(MM)
    !            DUDX(I)=DUDX(I)+DUDXE(MM)*SELE(MM)
    !            DVDY(I)=DVDY(I)+DVDYE(MM)*SELE(MM)
    !            NUM=NUM+1
    !        ELSE
    !            EXIT
    !        END IF
    !        IF(NUM==0) THEN
    !            PRINT*, 'HJP99'
    !            PAUSE
    !            STOP
    !        END IF
    !    END DO
    !    PHI(I)=PHI(I)/SALL
    !    DUDX(I)=DUDX(I)/SALL
    !    DVDY(I)=DVDY(I)/SALL
    !END DO

    DO  M=1,NRBN  !Only the OutFlow Boundary needs this modified condition
        KU=RBN(M)
        !UU=VisE*DUDX(KU)-1/4*(u(KU)**2+v(KU)**2)*(1-tanh(u(KU)/(1*1/20)))  !Prescribe Pressure Boundary Condition
        UU=0.0D0
        I=IA(KU+1)-1
        AP(I)=1
        BP(KU)=UU
        DO J=IA(KU),IA(KU+1)-2
            K=JA(J)
            BP(K)=BP(K)-AP(J)*UU
            AP(J)=0.0D0
        END DO
        DO  I=KU+1,NN
            DO  J=IA(I),IA(I+1)-2
                K=JA(J)
                IF(K==KU) THEN
                    BP(I)=BP(I)-AP(J)*UU
                    AP(J)=0.0D0
                END IF
            END  DO
        END DO
    END DO

    END SUBROUTINE PRESSURE_1BC
    !
    !
    !
    !
    !
    !
    SUBROUTINE E_ET_Q
    !	PRE_CONDITIONED: AP=>E=>E<T>=>Q=E*E<T>
    USE  CONSTANT
    USE  VARIABLE

    IMPLICIT NONE
    INTEGER:: I,J,K,J1,J2,M
    REAL(8):: COOO

    IF(ALLOCATED(EAP  )) DEALLOCATE(EAP  )
    IF(ALLOCATED(EAPT )) DEALLOCATE(EAPT )
    IF(ALLOCATED(JEAPT)) DEALLOCATE(JEAPT)
    IF(ALLOCATED(IEAPT)) DEALLOCATE(IEAPT)


    ALLOCATE (EAP  (N1DA))
    ALLOCATE (EAPT (N1DA))
    ALLOCATE (JEAPT(N1DA))
    ALLOCATE (IEAPT (NN+1)   )

    ALLOCATE (R8ARRAY(NN))
    ALLOCATE (I4ARRAY(NN))

    DO I=1,NN
        K=IA(I+1)-1
        COOO=AP(K)

        IF(COOO<0.0D0) THEN
            WRITE(*,*) 'DXS0'
            PAUSE
            STOP
        ELSE
            EAP(K)=DSQRT(COOO)   !将AP开根号
        END IF

        R8ARRAY(I)=EAP(K)
    END DO

    I4ARRAY=0
    DO I=1,NN
        J1=IA(I)
        J2=IA(I+1)-1
        DO J=J1,J2
            K=JA(J)
            I4ARRAY(K)=I4ARRAY(K)+1
            IF(J/=J2) THEN
                EAP(J)=WMGA*AP(J)/R8ARRAY(K)
            END IF
        END DO
    END DO

    IEAPT(1)=1
    DO I=1,NN
        IEAPT(I+1)=IEAPT(I)+I4ARRAY(I)
        I4ARRAY(I)=0
    END DO

    DO I=1,NN
        DO J=IA(I),IA(I+1)-1
            K=JA(J)
            I4ARRAY(K)=I4ARRAY(K)+1
            M=IEAPT(K)+I4ARRAY(K)-1
            EAPT(M)=EAP(J)
            JEAPT(M)=I
        END DO
    END DO
    DEALLOCATE(R8ARRAY,I4ARRAY)

    END SUBROUTINE E_ET_Q
    !
    !
    !
    !
    !
    SUBROUTINE SOLUTION_BICGS (AAA,BBB,XXX)                 !没看
    !*	SOLUTION FOR LINEAR SYSTEM WITH PRECONDITIONED
    !*	BI-CONJUGATE-GRADIENT STABLE(BICSTAB).
    !*	AAA: ONE DIMENSION COEFFICIENT ARRAY;
    !*	BBB: RIGHT HAND TERM; XXX: UNKNOWN VARIABLE.
    USE CONSTANT
    USE VARIABLE

    IMPLICIT NONE
    INTEGER:: KK,ITERA_MAX,I,J,K
    REAL(8):: AAA(N1DA),BBB(NN),XXX(NN)
    REAL(8):: ROL,RAQ,ALF,W,ERROR,EPS,ROLI,BATE,WU,WD

    ALLOCATE( BICGS_PI(NN) )
    ALLOCATE( BICGS_R(NN)  )
    ALLOCATE( DRO(NN)      )
    ALLOCATE( BICGS_Y(NN)  )
    ALLOCATE( BICGS_S(NN)  )
    ALLOCATE( BICGS_T(NN)  )
    ALLOCATE( BICGS_Z(NN)  )
    ALLOCATE( NIU(NN)      )
    ALLOCATE( INVE_T(NN)   )
    ALLOCATE( INVE_S(NN)   )
    ALLOCATE( BICG_TEMPARRAY(NN) )

    BICGS_PI=0.0D0
    BICGS_R=0.0D0
    DRO=0.0D0
    BICGS_Y=0.0D0
    BICGS_S=0.0D0
    BICGS_T=0.0D0
    BICGS_Z=0.0D0
    NIU=0.0D0
    BICG_TEMPARRAY=0.0D0
    INVE_T=0.0D0
    INVE_S=0.0D0

    KK=0
    ITERA_MAX=1000

    BICG_TEMPARRAY =0.0D0
    DO I=1,NN
        DO J=IA(I),IA(I+1)-2
            K=JA(J)
            BICG_TEMPARRAY(I)=BICG_TEMPARRAY(I)+AAA(J)*XXX(K)
            BICG_TEMPARRAY(K)=BICG_TEMPARRAY(K)+AAA(J)*XXX(I)
        END DO
        BICG_TEMPARRAY(I)=BICG_TEMPARRAY(I)+AAA(IA(I+1)-1)*XXX(I)
    END DO

    BICGS_R     = BBB-BICG_TEMPARRAY
    DRO         = BICGS_R
    BICGS_PI    = 0.00D0
    NIU         = 0.00D0
    ROL         = 1.00D0
    ALF         = 1.00D0
    W           = 1.00D0
    ERROR       = 10.0D0
    KK          = 0

    DO KK=1,ITERA_MAX

        ROLI=0.0D0
        DO I=1,NN
            ROLI=ROLI+BICGS_R(I)*DRO(I)
        END DO

        BATE=ROLI/ROL*ALF/W
        ROL=ROLI

        DO I=1,NN
            BICGS_PI(I)=BICGS_R(I)+BATE*(BICGS_PI(I)-W*NIU(I))
        END DO

        CALL CHOLESKY(BICGS_Y,BICGS_PI,0)

        NIU=0.0D0
        DO I=1,NN
            DO J=IA(I),IA(I+1)-2
                K=JA(J)
                NIU(I)=NIU(I)+AAA(J)*BICGS_Y(K)
                NIU(K)=NIU(K)+AAA(J)*BICGS_Y(I)
            END DO
            NIU(I)=NIU(I)+AAA(IA(I+1)-1)*BICGS_Y(I)
        END DO

        RAQ=0.0D0
        DO I=1,NN
            RAQ=RAQ+DRO(I)*NIU(I)
        END DO

        ALF=ROLI/RAQ

        DO I=1,NN
            BICGS_S(I)=BICGS_R(I)-ALF*NIU(I)
        END DO

        CALL CHOLESKY(BICGS_Z,BICGS_S,0)

        BICGS_T=0.0D0
        DO I=1,NN
            DO J=IA(I),IA(I+1)-2
                K=JA(J)
                BICGS_T(I)=BICGS_T(I)+AAA(J)*BICGS_Z(K)
                BICGS_T(K)=BICGS_T(K)+AAA(J)*BICGS_Z(I)
            END DO
            BICGS_T(I)=BICGS_T(I)+AAA(IA(I+1)-1)*BICGS_Z(I)
        END DO

        CALL CHOLESKY( INVE_T,BICGS_T,1 )
        CALL CHOLESKY( INVE_S,BICGS_S,1 )

        WU=0.0D0
        WD=0.0D0

        DO I=1,NN
            WU=WU+INVE_T(I)*INVE_S(I)
            WD=WD+INVE_T(I)*INVE_T(I)
        END DO
        W=WU/WD

        DO I=1,NN
            XXX(I)=XXX(I)+ALF*BICGS_Y(I)+W*BICGS_Z(I)
            BICGS_R(I)=BICGS_S(I)-W*BICGS_T(I)
        END DO




        ERROR=0.0D0
        DO I=1,NN
            ERROR=ERROR+DABS(BICGS_R(I))
        END DO
        ERROR=ERROR/FLOAT(NN)

        IF (ERROR <= EPS_BICGSTAB) THEN
            WRITE(*,'(A10,I12,10X,A10,E12.4)') 'INTR NUM =' ,KK, 'ERROR    =' ,ERROR
            EXIT
        END IF


    END DO

    IF (ERROR > EPS_BICGSTAB) THEN
        WRITE(*,*)
        WRITE(*,*) ' NUM-INTERATION=',KK-1
        WRITE(*,*) ' ERROR         =',ERROR
        WRITE(*,*) ' EPS_BICGSTAB  =',EPS_BICGSTAB
        WRITE(*,*) ' SOLUTION IS DIVERGED! FAILURE IN SOLUTION! '
        PAUSE
        STOP
    END IF

    DEALLOCATE( BICGS_PI,BICG_TEMPARRAY,BICGS_S,BICGS_T,DRO)
    DEALLOCATE( BICGS_Z,INVE_T,INVE_S,NIU,BICGS_Y,BICGS_R)

    END SUBROUTINE SOLUTION_BICGS
    !
    !
    !
    !
    !
    SUBROUTINE CHOLESKY (XXX,BBB,MARKER)
    !*	SOLUTION FOR TRIGONAL MATRIX ADOPTING CHOLESKY METHOD.
    USE CONSTANT
    USE VARIABLE

    IMPLICIT NONE
    INTEGER:: MARKER,I,J1,J2,J,K
    REAL(8):: XXX(NN), BBB(NN),A
    REAL(8),ALLOCATABLE:: YYY(:)

    ALLOCATE(YYY(NN))
    YYY(1)=BBB(1)/EAP(IA(1))

    DO I=2,NN
        A=0.0D0
        J1=IA(I)
        J2=IA(I+1)-2
        DO J=J1,J2
            K=JA(J)
            A=A+EAP(J)*YYY(K)
        END DO
        YYY(I)=(BBB(I)-A)/EAP(J2+1)
    END DO

    IF (MARKER==1) THEN
        XXX=YYY
        RETURN
    END IF

    XXX(NN)=YYY(NN)/EAPT(IEAPT(NN))

    DO I=NN-1,1,-1
        A=0.0D0
        J1=IEAPT(I)+1
        J2=IEAPT(I+1)-1
        DO J=J1,J2
            K=JEAPT(J)
            A=A+EAPT(J)*XXX(K)
        END DO
        XXX(I)=(YYY(I)-A)/EAPT(J1-1)
    END DO

    END SUBROUTINE CHOLESKY
    !
    !
    !
    !
    !
    SUBROUTINE RHPOISSON0
    USE CONSTANT
    USE VARIABLE
    IMPLICIT NONE

    INTEGER:: ELE,NET(4),I,J,K,NK,IE
    REAL(8):: FD,Q, DUDX0, DVDY0, VIST, DUDX,DUDY,DVDX,DVDY,UU,VV,FYG


    DO ELE=1,NE
        NET(:)=NODEELE(ELE,:)

        DO I=1,2
            DO J=1,2

                CALL GAUSS_JACOB(I,J,ELE)
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                Q=FD*WEI(I)*WEI(J)

                UU   =0.0D0;	VV   =0.0D0;	DUDX =0.0D0;	DUDY =0.0D0
                DVDX =0.0D0;	DVDY =0.0D0;
                DUDX0=0.0D0;	DVDY0=0.0D0
                !				FYG  =0.0D0

                DO K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(I,J,K)-YPDKSI12*FPDITA(I,J,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(I,J,K)-XPDITA21*FPDKSI(I,J,K))/FD

                    NK=NET(K)

                    UU=UU+ (U(NK)-UM_ALE(NK))*F(I,J,K)	!!!!!!!!!!!!!!!!!!!!!
                    VV=VV+ (V(NK)-VM_ALE(NK))*F(I,J,K)	!!!!!!!!!!!!!!!!!!!!!

                    DUDX= DUDX + FPDX(K)*U(NK)
                    DUDY= DUDY + FPDY(K)*U(NK)
                    DVDX= DVDX + FPDX(K)*V(NK)
                    DVDY= DVDY + FPDY(K)*V(NK)
                    !					FYG = FYG  + 0.0D0 !FALE_ALL(NK)*F(I,J,K)

                    DUDX0=DUDX0+FPDX(K)*U0(NK)
                    DVDY0=DVDY0+FPDY(K)*V0(NK)
                END DO
                !					PRINT*, FYG

                DO K=1,4
                    NK=NET(K)
                    BP(NK)= BP(NK) + F(I,J,K) * ( - (DUDX0+DVDY0))* Q/DT &
                        - ( FPDX(K)*(UU*DUDX+VV*DUDY) + FPDY(K)*(UU*DVDX+VV*DVDY) )*Q !+ ( FPDY(K)*FYG*Q )
                END DO

            END DO
        END DO
    END DO

    END SUBROUTINE RHPOISSON0



    SUBROUTINE DYNAMIC_DT
    USE CONSTANT,ONLY: DT000
    USE VARIABLE,ONLY: NE,IPE,NODEELE,U,V,X,Y,DT
    IMPLICIT NONE

    INTEGER:: I,J,K,NODE,NK
    REAL(8):: UU, VV, AUV,  R8HUGE,  DDTT, DTMIN,XX,YY,XP(4),YP(4),X1,X2,X3,Y1,Y2,Y3,SS,LL(4),LMIN

    R8HUGE=HUGE(R8HUGE)

    DO I=1,NE

        UU=0.0D0
        VV=0.0D0
        XX=0.0D0
        YY=0.0D0

        NK=IPE(I)

        IF(NK/=4) THEN
            WRITE(*,*)' HJYR55'
            PAUSE
            STOP
        END IF

        DO K=1,NK
            NODE=NODEELE(I,K)
            UU=UU+U(NODE)
            VV=VV+V(NODE)
            XP(K)=X(NODE)
            YP(K)=Y(NODE)
            XX=XX+XP(K)
            YY=YY+YP(K)
        END DO
        

        UU=UU/ FLOAT(NK)
        VV=VV/ FLOAT(NK)
        XX=XX/FLOAT(NK)
        YY=YY/FLOAT(NK)
        AUV=DSQRT(UU*UU+VV*VV)

        DO K=1,NK
            IF(K<=NK-1) THEN
                X1=XP(K)
                X2=XP(K+1)
                Y1=YP(K)
                Y2=YP(K+1)
            ELSE
                X1=XP(K)
                X2=XP(1)
                Y1=YP(K)
                Y2=YP(1)
            END IF
            X3=XX
            Y3=YY

            SS   = 0.5D0*((X1*Y2-X2*Y1)+(X2*Y3-X3*Y2)+(X3*Y1-X1*Y3))
            IF(SS<0.0D0) THEN
                WRITE(*,*)' FVTR7','      ELEMENT:',I
                call OUTPUT
                PAUSE
                STOP
            END IF

            LL(K)= (2.0D0*SS)/(DSQRT((X2-X1)**2.0D0+(Y2-Y1)**2.0D0))
        END DO

        LMIN=MINVAL(LL)
        DDTT=LMIN/AUV


        IF(DDTT<=R8HUGE) THEN
            DTMIN=DDTT
            R8HUGE=DDTT
        END IF
    END DO


    DTMIN=DTMIN*0.5D0		! JASMIN

    IF (DTMIN > DT000) THEN
        DT= DT000
    ELSE
        DT=DTMIN
    END IF

    WRITE(*,'(A10,F12.6,10X,A10,F12.6)') 'DT       =',DT, 'DTMIN    =',DTMIN

    END SUBROUTINE DYNAMIC_DT

    SUBROUTINE BOUNDARY_INTEGRAL_PREEQU
    USE VARIABLE
    IMPLICIT NONE
    INTEGER:: N1,N2,IWB,I,K,NK,IE,NUM,IW,NET(4)
    REAL(8):: X1,Y1,X2,Y2,FNX,FNY,S,DELTAV,DELTAU

    DO NUM=1,NLBE_ORIG       !入流边界边界积分
        IE=LBE_ORIG(NUM,1)   !单元
        IW=LBE_ORIG(NUM,2)   !节点
        NET(1:4)=NODEELE(IE,1:4)

        IF(IW==1) THEN
            N1=NET(1)
            N2=NET(2)
            IWB=4
        ELSE IF(IW==2) THEN
            N1=NET(2)
            N2=NET(3)
            IWB=3
        ELSE IF(IW==3) THEN
            N1=NET(3)
            N2=NET(4)
            IWB=3
        ELSE
            N1=NET(4)
            N2=NET(1)
            IWB=4
        END IF

        X1=X(N1); Y1=Y(N1)
        X2=X(N2); Y2=Y(N2)

        CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)                  !FNX=(Y2-Y1)/DS      FNY=(X1-X2)/DS

        DO I=1,2
            IF((IW==1).OR.(IW==3)) THEN
                CALL GAUSS_JACOB(I,IWB,IE)
                S=WEI(I)*DSQRT(XPDKSI11**2.0D0+YPDKSI12**2.0D0)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    !考虑对偏P偏X进行近似出来后的结果，DELTA U 约掉了
                    !BP(NK)=BP(NK)-(U0(NK)*FNX+V0(NK)*FNY)*F(I,IWB,K)*S

                    !原版本
                    !BP(NK)=BP(NK)-((THETA(1)*DELTAU+U0(NK))*FNX+(THETA(1)*DELTAV+V0(NK))*FNY)*F(I,IWB,K)*S

                    !第二步修改
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(I,IWB,K)*S
                END DO
            ELSE
                CALL GAUSS_JACOB(IWB,I,IE)
                S=WEI(I)*DSQRT(XPDITA21*XPDITA21+YPDITA22*YPDITA22)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(IWB,I,K)*S
                END DO
            END IF
        END DO
    END DO


    DO NUM=1,NPBE           !圆柱表面边界积分
        IE=PBE(NUM,1)   !单元
        IW=PBE(NUM,2)   !节点
        NET(1:4)=NODEELE(IE,1:4)

        IF(IW==1) THEN
            N1=NET(1)
            N2=NET(2)
            IWB=4
        ELSE IF(IW==2) THEN
            N1=NET(2)
            N2=NET(3)
            IWB=3
        ELSE IF(IW==3) THEN
            N1=NET(3)
            N2=NET(4)
            IWB=3
        ELSE
            N1=NET(4)
            N2=NET(1)
            IWB=4
        END IF

        X1=X(N1); Y1=Y(N1)
        X2=X(N2); Y2=Y(N2)

        CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)                  !FNX=(Y2-Y1)/DS      FNY=(X1-X2)/DS

        DO I=1,2
            IF((IW==1).OR.(IW==3)) THEN
                CALL GAUSS_JACOB(I,IWB,IE)
                S=WEI(I)*DSQRT(XPDKSI11**2.0D0+YPDKSI12**2.0D0)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(I,IWB,K)*S
                END DO
            ELSE
                CALL GAUSS_JACOB(IWB,I,IE)
                S=WEI(I)*DSQRT(XPDITA21*XPDITA21+YPDITA22*YPDITA22)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(IWB,I,K)*S
                END DO
            END IF
        END DO
    END DO

    DO NUM=1,NSBE         !尾板表面边界积分
        IE=SBE(NUM,1)   !单元
        IW=SBE(NUM,2)   !节点
        NET(1:4)=NODEELE(IE,1:4)

        IF(IW==1) THEN
            N1=NET(1)
            N2=NET(2)
            IWB=4
        ELSE IF(IW==2) THEN
            N1=NET(2)
            N2=NET(3)
            IWB=3
        ELSE IF(IW==3) THEN
            N1=NET(3)
            N2=NET(4)
            IWB=3
        ELSE
            N1=NET(4)
            N2=NET(1)
            IWB=4
        END IF

        X1=X(N1); Y1=Y(N1)
        X2=X(N2); Y2=Y(N2)

        CALL NXNY(X1,Y1,X2,Y2,FNX,FNY)                  !FNX=(Y2-Y1)/DS      FNY=(X1-X2)/DS

        DO I=1,2
            IF((IW==1).OR.(IW==3)) THEN
                CALL GAUSS_JACOB(I,IWB,IE)
                S=WEI(I)*DSQRT(XPDKSI11**2.0D0+YPDKSI12**2.0D0)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(I,IWB,K)*S
                END DO
            ELSE
                CALL GAUSS_JACOB(IWB,I,IE)
                S=WEI(I)*DSQRT(XPDITA21*XPDITA21+YPDITA22*YPDITA22)
                DO K=1,4
                    NK=NET(K)
                    DELTAU=U(NK)-U0(NK)
                    DELTAV=V(NK)-V0(NK)
                    BP(NK)=BP(NK)-((THETA(1)*DELTAU)*FNX+(THETA(1)*DELTAV)*FNY)*F(IWB,I,K)*S
                END DO
            END IF
        END DO
    END DO


    END SUBROUTINE BOUNDARY_INTEGRAL_PREEQU


    SUBROUTINE MIN_X0_X00
    USE VARIABLE
    USE CONSTANT
    USE FlexibleOnly
    IMPLICIT NONE
    REAL(8):: I,J,X0,FX_temp(NSBN_ALE,1),xFtemp,xxFtemp
    Integer::NKindSN


    X00=HUGE(X00)
    DO I=1,NSBN_ALE
        J=SBN_ALE(I)
        X0=X_ALE(J)
        IF (X00>X0) THEN
            X00=X0
        END IF
    END DO

    WRITE(*,*) 'X00=',X00

    DO I = 1, NSBN_ALE
        J=SBN_ALE(I)
        FX_temp(I,1)=X_ALE(J)
    END DO

    DO I =1,NSBN_ALE
        xFtemp=FX_temp(I,1)
        Do J = I+1,NSBN_ALE
            xxFtemp =FX_temp(J,1)
            IF (xFtemp==xxFtemp) THEN
                FX_temp(J,1)=0
            END IF
        End DO
    End do

    NKindSN=0
    DO I =1,NSBN_ALE
        xFtemp=FX_temp(I,1)
        IF(xFtemp /= 0) then
            NKindSN =  NKindSN + 1
        End if
    End do


    Allocate(X_Plate_Orig(NKindSN,1))
    X_Plate_Orig = 0.0

    Allocate(FA_THETA(NKindSN),FA_THETA0(NKindSN))
    FA_THETA = 0.0
    FA_THETA0 = 0.0

    X_Plate_Orig = 0.0D0
    J = 1
    DO I = 1, NSBN_ALE
        IF (FX_temp(I,1)/=0)THEN
            X_Plate_Orig(J,1)=FX_temp(I,1)
            J=J+1
        END IF
    END DO

    Call Sort_Ascend_Multi(X_Plate_Orig,NKindSN,1)

    Allocate(X_CA_ORR(NKindSN,2))
    X_CA_ORR(:,1) = X_Plate_Orig(:,1)    !x
    X_CA_ORR(:,2) = 0.0             !y


    END SUBROUTINE MIN_X0_X00

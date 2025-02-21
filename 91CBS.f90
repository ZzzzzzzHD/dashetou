    SubRoutine STEP1
    !This SubRoutine is aimed to assemble RHS vetor of step 1
    USE CONSTANT
    USE VARIABLE

    Implicit None
    !Local variables
    Real(8):: LdeltaT_half
    INTEGER:: IE,IEOLD,IET,NET(4),K,NK,I,J
    REAL(8):: FD,Q,QR,UU,VV,DUDX,DUDY,DVDX,DVDY,PT,VT,UT,PG,DPDX,DPDY
    REAL(8):: VisE, VisM, VisT, F0, Dts ,ROU,umean,vmean
    !Modified data/variables 
    !=============================================================
    !Attention1: temporally theta2=1, further i will modify theta2!
    !Attention2: temporally there isn't mesh movement!
    !=============================================================

    ROU=DENSITY
    VisM=VISCOSITY/density

    Do ie =1,ne
        LdeltaT_half=DT*0.5d0
        net=NodeEle(ie,1:4)
        Do i=1,2
            Do j=1,2
                CALL GAUSS_JACOB(I,J,IE)
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                Q=FD*WEI(I)*WEI(J)

                UU   =0.0d0;  	VV   =0.0d0;
                DPDX =0.0d0; 	DPDY =0.0d0;
                DUDX =0.0d0;	DUDY =0.0d0;	DVDX =0.0d0;	DVDY =0.0d0
                PG   =0.0d0;

                DO  K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(I,J,K)-YPDKSI12*FPDITA(I,J,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(I,J,K)-XPDITA21*FPDKSI(I,J,K))/FD
                    NK=NET(K)
                    UT=U0(NK)
                    VT=V0(NK)
                    PT=PN(NK)

                    UU=UU+ F(I,J,K)* ((UT) -uM_Ale(nk))
                    VV=VV+ F(I,J,K)* ((VT) -vM_Ale(nk))

                    !UU=UU+ F(I,J,K)* UT
                    !VV=VV+ F(I,J,K)* VT

                    PG=PG+ F(I,J,K)*  PT

                    DUDX=DUDX+FPDX(K)*UT; 		DUDY=DUDY+FPDY(K)*UT
                    DVDX=DVDX+FPDX(K)*VT; 		DVDY=DVDY+FPDY(K)*VT
                    DPDX=DPDX+FPDX(K)*PT; 		DPDY=DPDY+FPDY(K)*PT
                END DO
                umean=uu/float(4)
                vmean=vv/float(4)

                VisT=(Cs*CLE(IE))**2.0D0* DSQRT( 2.0D0*(DUDX**2.0D0 +DVDY**2.0D0 )+(DUDY+DVDX)**2.0D0 )
                VisE=VisM+VisT

                Do k=1,4
                    !!Attention UB is now delta u^*, differs from UB in plate program!!
                    nk=net(k)
                    UB(nk)=UB(nk) &
                        + Q  * dt*( -F(i,j,k) * (uu * dudx + vv * dudy ) &
                        !-   LdeltaT_half *  FPDX(K) * umean * ( uu * dudx + vv * dudy ) &
                        -   LdeltaT_half * ( FPDX(K)*UU+FPDY(K)*VV)*( uu * dudx + vv * dudy)  &
                        -           VisE * ( FPDX(K) * DUDX + FPDY(K) * DUDY ))

                    VB(nk)=VB(nk) &
                        +Q  *dt* ( -F(I,J,K) *( uu * dvdx + vv * dvdy ) &
                        ! -LdeltaT_half    *  FPDY(K) * vmean * ( uu * dvdx + vv * dvdy ) &
                        - LdeltaT_half * ( FPDX(K)*UU+FPDY(K)*VV)*( uu * dvdx + vv * dvdy ) &
                        - VisE            * ( FPDX(K) * DVDX + FPDY(K) * DVDY ))

                End Do
            End Do
        End Do
    End Do
    End SubRoutine STEP1


    SubRoutine STEP1_Solution
    USE CONSTANT
    USE VARIABLE
    integer::i

    Do i=1,nn
        Ub(I)=U(I)+Ub(I)*(1.0D0/AM(I)) !这里计算出来的是delta u
        Vb(I)=V(I)+Vb(I)*(1.0D0/AM(I))
    End Do
    End SubRoutine


    SubRoutine STEP2
    USE CONSTANT
    USE VARIABLE
    Implicit None
    !local variables
    INTEGER:: ELE,NET(4),I,J,K,NK,ie
    REAL(8):: FD,Q, DUDX0, DVDY0, vist, dudx,dudy,dvdx,dvdy,uu,vv,DeltaU,Deltav,ddudx,ddudy,ddvdx,ddvdy
    !=======================step 2 begins===================
    do ele=1,ne
        net(:)=nodeele(ele,:)
        do i=1,2
            do j=1,2
                CALL gauss_jacob(i,j,ele)
                fd=xpdksi11*ypdita22-xpdita21*ypdksi12
                q=fd*wei(i)*wei(j)
                uu   =0.0d0;	vv   =0.0d0;	dudx =0.0d0;	dudy =0.0d0
                dvdx =0.0d0;	dvdy =0.0d0;
                dudx0=0.0d0;	dvdy0=0.0d0;    DeltaU=0.0d0;   DeltaV=0.0d0
                ddudx=0.0d0;    ddudy=0.0d0;    ddvdx=0.0d0;    ddvdy=0.0d0
                do k=1,4
                    fpdx(k)=(ypdita22*fpdksi(i,j,k)-ypdksi12*fpdita(i,j,k))/fd
                    fpdy(k)=(xpdksi11*fpdita(i,j,k)-xpdita21*fpdksi(i,j,k))/fd

                    nk=net(k)

                    uu=uu+ U0(nk)*f(i,j,k)
                    vv=vv+ V0(nk)*f(i,j,k)

                    dudx= dudx + fpdx(k)*u0(nk)
                    dudy= dudy + fpdy(k)*u0(nk)
                    dvdx= dvdx + fpdx(k)*v0(nk)
                    dvdy= dvdy + fpdy(k)*v0(nk)

                    DeltaU=DeltaU+(U(NK)-U0(nk))*f(i,j,k)
                    DeltaV=DeltaV+(v(nk)-v0(nk))*f(i,j,k)

                    ddudx= ddudx + fpdx(k)*ub(nk)
                    !ddudy= ddudy + fpdy(k)*ub(nk)
                    !ddvdx= ddvdx + fpdx(k)*vb(nk)
                    ddvdy= ddvdy + fpdy(k)*vb(nk)

                end do
                DO K=1,4
                    NK=NET(K)
                    !BP(NK)= BP(NK) + F(i,j,k) * ( - (DUDX0+DVDY0))* Q/DT &
                    ! - ( FPDX(K)*(uu*dudx+vv*dudy) + FPDY(K)*(uu*dvdx+vv*dvdy) )*Q
                    !!这里可能有问题!!
                    BP(NK)= BP(NK)&
                        !+ Q* (( FPDX(K)*UU+FPDY(K)*VV ) &
                        + Q*(-F(I,J,K)*(DUDX+DVDY) &
                        + theta(1) * ( FPDX(K)* DeltaU + FPDY(K) * DeltaV ))
                    !+ theta(1) * ( -f(i,j,k)* (ddvdx+ddvdy)))
                END DO
            End Do
        End Do
    End Do

    End Subroutine STEP2



    Subroutine STEP3
    USE CONSTANT
    USE VARIABLE

    Implicit None

    !local

    INTEGER:: IE,IEOLD,IET,NET(4),K,NK,I,J
    REAL(8):: FD,Q,QR,UU,VV,DUDX,DUDY,DVDX,DVDY,PT,VT,UT,PG,DPDX,DPDY,Umean,Vmean,PTOLD,DPODX,DPODY

    Do ie=1,ne
        net=NodeEle(ie,1:4)
        Do i=1,2
            Do j=1,2
                CALL GAUSS_JACOB(I,J,IE)
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                Q=FD*WEI(I)*WEI(J)

                DPDX =0.0d0; 	DPDY =0.0d0;
                PT   =0.0d0;    UU=0.0d0;      VV=0.0d0;
                Umean=0.0d0;    Vmean=0.0d0;    UT=0.0D0;   VT=0.0D0
                DPODX=0.0D0;    DPODY=0.0D0

                DO  K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(I,J,K)-YPDKSI12*FPDITA(I,J,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(I,J,K)-XPDITA21*FPDKSI(I,J,K))/FD

                    NK=NET(K)
                    UT=U0(NK)
                    VT=V0(NK)

                    UU = UU + F(I,J,K)* ((UT) -uM_Ale(nk))
                    VV = VV + F(I,J,K)* ((VT) -vM_Ale(nk))

                    PT=PN1(NK)
                    DPDX=DPDX+FPDX(K)*PT; 		DPDY=DPDY+FPDY(K)*PT

                    PTOLD=PN(NK)
                    DPODX=DPODX+FPDX(K)*PTOLD; 		DPODY=DPODY+FPDY(K)*PTOLD

                END DO

                Do k=1,4
                    NK=NET(K)
                    !UB2(nk)= UB2(NK) + (Q * ( -DT * F(I,J,K)*DPDX) ) / density
                    !VB2(nk)= VB2(NK) + (Q * ( -DT * F(I,J,K)*DPDY) ) / density
                    UB2(nk)= UB2(NK) + (Q * (( -DT * F(I,J,K)*DPDX) + DT**2.0D0 / 2.0D0 *(FPDX(K)*UU+FPDY(K)*VV) * DPODX))/density
                    VB2(nk)= VB2(NK) + (Q * (( -DT * F(I,J,K)*DPDY) + DT**2.0D0 / 2.0D0 *(FPDX(K)*UU+FPDY(K)*VV) * DPODY))/density
                END DO
            END DO
        END DO
    END DO
    End SubRoutine STEP3



    SubRoutine STEP3_Solution
    USE CONSTANT
    USE VARIABLE
    integer::i

    Do i=1,nn
        UB2(I)=UB2(I)*(1.0D0/AM(I)) !这里计算出来的是delta u
        VB2(I)=VB2(I)*(1.0D0/AM(I))
    End Do

    End SubRoutine


    SubRoutine OutFlow_Boundary
    Use CONSTANT
    Use VARIABLE

    Implicit None
    INTEGER:: N1,N2,IWB,I,J,K,NK,IE,NUM,IW,NET(4),KK,KM,ITY,II,JJ
    REAL(8):: X1,Y1,X2,Y2,FNX,FNY,S,DELTAV,DELTAU,Y111,Y000,NODE,DUDX,DVDY,FD,Q,UT,VT
    REAL(8):: VisE, VisM, VisT, ROU
    REAL(8):: ThetaOB,AlphaOB1,AlphaOB2

    ThetaOB=1
    AlphaOB1=0
    AlphaOB2=0

    ROU=DENSITY
    VisM=VISCOSITY/density
    VisT=0.0D0
    VisE=VisM+VisT

    DO NUM=1,NRBE       !In fact, Its a Boundary Integration of variable U added by DONG(2014)
        IE=RBE(NUM,1)   !单元
        IW=RBE(NUM,2)   !节点
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
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                DVDY =0.0d0
                DUDX =0.0d0

                DO K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(I,IWB,K)-YPDKSI12*FPDITA(I,IWB,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(I,IWB,K)-XPDITA21*FPDKSI(I,IWB,K))/FD
                    NK=NET(K)
                    UT=U(NK)
                    VT=V(NK)
                    DUDX=DUDX+FPDX(K)*UT
                    DVDY=DVDY+FPDY(K)*VT
                END DO

                CALL GAUSS_JACOB(I,IWB,IE)
                S=WEI(I)*DSQRT(XPDKSI11**2.0D0+YPDKSI12**2.0D0)
                DO K=1,4
                    NK=NET(K)
                    !UB(NK)=UB(NK)+dt*(PN(nk)*FNX+0.5*(U(NK)**2+V(NK)**2)*0.5*(1-TANH((FNX*U(NK)+FNY*V(NK))/1/(1/20)))-VisE*(PUPX+PVPY)*FNX)*F(I,IWB,K)*S
                    UB(NK)=UB(NK)+dt*(PN(nk)*FNX+((ThetaOb+AlphaOB2)*0.5*(U(NK)**2+V(NK)**2)*FNX&
                        &+(1-thetaOB+AlphaOB1)*0.5*(FNX*U(NK)+FNY*V(NK))*U(NK))&
                        &*0.5*(1-TANH((FNX*U(NK)+FNY*V(NK))/1/(1/20)))-VisE*(DUDX+DVDY))*F(I,IWB,K)*S

                    VB(NK)=VB(NK)+dt*(PN(nk)*FNY+((ThetaOb+AlphaOB2)*0.5*(U(NK)**2+V(NK)**2)*FNY&
                        &+(1-thetaOB+AlphaOB1)*0.5*(FNX*U(NK)+FNY*V(NK))*V(NK))&
                        &*0.5*(1-TANH((FNX*U(NK)+FNY*V(NK))/1/(1/20)))-VisE*(DUDX+DVDY))*F(I,IWB,K)*S
                END DO
            ELSE
                CALL GAUSS_JACOB(IWB,I,IE)
                FD=XPDKSI11*YPDITA22-XPDITA21*YPDKSI12
                DVDY =0.0d0
                DUDX =0.0d0

                DO K=1,4
                    FPDX(K)=(YPDITA22*FPDKSI(IWB,I,K)-YPDKSI12*FPDITA(IWB,I,K))/FD
                    FPDY(K)=(XPDKSI11*FPDITA(IWB,I,K)-XPDITA21*FPDKSI(IWB,I,K))/FD
                    NK=NET(K)
                    UT=U(NK)
                    VT=V(NK)
                    DUDX=DUDX+FPDX(K)*UT
                    DVDY=DVDY+FPDY(K)*VT
                END DO

                CALL GAUSS_JACOB(IWB,I,IE)
                S=WEI(I)*DSQRT(XPDITA21*XPDITA21+YPDITA22*YPDITA22)
                DO K=1,4
                    NK=NET(K)
                    UB(NK)=UB(NK)+dt*(PN(nk)*FNX+((ThetaOb+AlphaOB2)*0.5*(U(NK)**2+V(NK)**2)*FNX&
                        &+(1-thetaOB+AlphaOB1)*0.5*(FNX*U(NK)+FNY*V(NK))*U(NK))&
                        &*0.5*(1-TANH((FNX*U(NK)+FNY*V(NK))/1/(1/20)))-VisE*(DUDX+DVDY))*F(IWB,I,K)*S
                    VB(NK)=VB(NK)+dt*(PN(nk)*FNY+((ThetaOb+AlphaOB2)*0.5*(U(NK)**2+V(NK)**2)*FNY&
                        &+(1-thetaOB+AlphaOB1)*0.5*(FNX*U(NK)+FNY*V(NK))*V(NK))&
                        &*0.5*(1-TANH((FNX*U(NK)+FNY*V(NK))/1/(1/20)))-VisE*(DUDX+DVDY))*F(IWB,I,K)*S
                END DO
            END IF
        END DO
    END DO

    End SubRoutine OutFlow_Boundary
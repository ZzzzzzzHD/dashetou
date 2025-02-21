
    MODULE CONSTANT
    IMPLICIT NONE

    INTEGER:: NSTEP=00, MAXSTEP=6000000, NW

    REAL(8):: DT000, Tout=0.000d0, dtout

    REAL(8):: GRAVITY=-9.81d0 , DENSITY=1.0D0

    REAL(8):: WMGA   =1.6d0   , EPS_BICGSTAB=1.0d-07
    REAL(8):: UUXX00=1.0d0, length_lmm,mass_add,F_0,F_1,T_00


    INTEGER:: rightlayer, leftlayer, uplayer, downlayer
    INTEGER:: ncc,ndr


    Real(8):: VISCOSITY,Cs

    REAL(8)::MassR ,T_LMM

    REAL(8):: fre_R, damping_R
    REAL(8):: FRE_X,DAMPING_X
    REAL(8):: FRE_Y,DAMPING_Y

    REAL(8):: omg_Ale,amc_ale,hh,t_start,cfl

    integer:: re_start,	mesh_transform, irota,CSD_METHOD
    character(12):: file_name

    real(8):: LL0,epsilon_ang

    real(8)::time_stop,TIME_CALL_FORCE

    real(8)spring_out,spring_in
    real(8)::DS_R,R_pg,DS_L,L_pg,DS_U,U_pg,DS_D,D_pg
    
    END MODULE CONSTANT
    !
    !
    !
    !
    !
    MODULE VARIABLE
    IMPLICIT NONE


    REAL(8),ALLOCATABLE::THETA(:),UB2(:),VB2(:)
    REAL(8)            ::NW_START,NW_END,ETA_R,ETA_T,ETA_L

    INTEGER             :: NE_ORIG, NN_ORIG   , NN,NE
    INTEGER,ALLOCATABLE:: NODEELE_ORIG(:,:)  , NODEELE(:,:)
    REAL(8),ALLOCATABLE:: X_ORIG(:),Y_ORIG(:), X(:),Y(:), cxye(:,:),NODE_NODE_LMM_Angle(:,:,:),SELE_Iter(:),DELTX(:),DELTY(:)
    INTEGER,ALLOCATABLE:: NODE_NODE_LMM(:,:),NODE_NODE_LMM_ELE(:,:,:)
    INTEGER,ALLOCATABLE:: IPE(:)


    INTEGER,ALLOCATABLE:: ELE_ADJ(:,:)  , NODE_ELE(:,:)  , NODE_NODE(:,:)
    INTEGER            :: MAXNUM_ELE_ADJ, maxnum_node_ele, maxnum_node_node


    REAL(8),SAVE :: XI(5),WEI(3),FPDX(4),FPDY(4),F(5,5,5),FPDKSI(5,5,5),FPDITA(5,5,5)
    REAL(8),SAVE :: XPDKSI11,XPDITA21,YPDKSI12,YPDITA22


    REAL(8):: T,TSUB,DT,TT0,x00,fa_spr,fusai_spr
    INTEGER:: NSUB,ITER_NUM

    Real(8):: Dis0_ang, Dis_ang, Vel0_ang, Vel_ang, Acc0_ang, Acc_ang, tor, tor0, dxRk, dxRf
    Real(8):: DIS0_X,DIS_X,VEL0_X,VEL_X,ACC0_X,ACC_X,DXXK,DXXF
    Real(8):: DIS0_Y,DIS_Y,VEL0_Y,VEL_Y,ACC0_Y,ACC_Y,DXYK,DXYF
    Real(8):: DIS0_Y_UP,DIS_Y_UP,VEL0_Y_UP,VEL_Y_UP,ACC0_Y_UP,ACC_Y_UP

    Real(8)::m1,m2,m3,m4

    REAL(8),ALLOCATABLE:: U(:),V(:),U0(:),V0(:),PN(:), PN1(:),vor(:)


    REAL(8),ALLOCATABLE:: UB(:),VB(:),BP(:),AM(:),AP(:),AM1(:),SELE(:),SELE_ALE(:),CLE(:)


    INTEGER,ALLOCATABLE:: IA(:),JA(:),IEAPT(:),JEAPT(:)
    INTEGER:: n1da

    REAL(8),ALLOCATABLE:: EAPT(:),EAP(:),BICGS_PI(:),BICGS_R(:),INVE_T(:),DRO(:),BICGS_Y(:)
    REAL(8),ALLOCATABLE:: INVE_S(:),BICGS_S(:),BICGS_T(:),BICGS_Z(:),NIU(:),BICG_TEMPARRAY(:)


    CHARACTER(LEN=80  ):: FILENAME,MESHFILENAME


    INTEGER,ALLOCATABLE:: I4ARRAY(:),I4MATRIX(:,:),ARRAY_LARGE(:),nodetype(:)
    REAL(8),ALLOCATABLE:: r8a(:),r8matrix(:,:),R8ARRAY(:)

    REAL(8),ALLOCATABLE:: VORE(:)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER:: nbe_orig, NLBE_ORIG, NUBE_ORIG, NRBE_ORIG, NDBE_ORIG, npbe_orig, nsbe_orig
    INTEGER:: nlbn_orig, nrbn_orig, ndbn_orig, nubn_orig, npbn_orig, nsbn_orig
    INTEGER,allocatable:: BE_ORIG(:,:), LBE_ORIG(:,:), UBE_ORIG(:,:), RBE_ORIG(:,:), DBE_ORIG(:,:), PBE_ORIG(:,:), sbe_orig(:,:)
    INTEGER,allocatable:: pbn_orig(:), dbn_orig(:), ubn_orig(:), lbn_orig(:), rbn_orig(:), sbn_orig(:)

    INTEGER:: nlbn_ALE, nrbn_ALE, ndbn_ALE, nubn_ALE, npbn_ALE,	nsbn_ale,						NN_ALE,NE_ALE
    INTEGER,allocatable:: pbn_ALE(:), dbn_ALE(:), ubn_ALE(:), lbn_ALE(:), rbn_ALE(:), sbn_ale(:),	NODEELE_ALE(:,:)
    REAL(8),allocatable:: x_ale(:),y_ale(:)
    INTEGER:: N1DA_ALE
    INTEGER,allocatable:: ia_ale(:), ja_ale(:)
    REAL(8),allocatable:: ap_ale(:),bp_ale(:)



    INTEGER,allocatable:: nodeele_temp(:,:)
    REAL(8),allocatable:: x_temp(:),y_temp(:)

    INTEGER:: nlbn, nrbn, ndbn, nubn, npbn, nsbn
    INTEGER,allocatable:: pbn(:), dbn(:), ubn(:), lbn(:), rbn(:), sbn(:)
    REAL(8),allocatable:: xbc_ale(:),ybc_ale(:)
    REAL(8),allocatable:: x0_ale(:),y0_ale(:),x1_ale(:),y1_ale(:)

    INTEGER,allocatable:: bn_ale(:)
    INTEGER:: nbn_ale

    REAL(8),allocatable:: apx_ale(:),apy_ale(:),um_ale(:),vm_ale(:)

    INTEGER:: npbe,nrbe,nsbe
    INTEGER,allocatable:: pbe(:,:),rbe(:,:),sbe(:,:)

    REAL(8)::WORK,WORK_TOR,WORK_F_Y,WORK_F_X!JY!w累计做功,TOR累计做功,CL累计做功,CD累计做功
    REAL(8)::CD0,CL0
    REAL(8)::CD,CL,CD1,CL1,CD2,CL2
    REAL(8)::FORCE_X_CENTER_0
    REAL(8)::CL_CUT!JY!计算升力中心时用

    REAL(8):: INT_DIS_ANG,INT_DIS_Y!JY!初始位置  1219：ANG_0.5_Y_0

    END MODULE VARIABLE
    !
    !
    !
    !
    !
    !
    MODULE FlexibleOnly
    Implicit None

    REAL(8)::PLATE_DS
    REAL(8),ALLOCATABLE::FA_THETA(:),FA_THETA0(:)

    REAL(8)::FL_L           ! 暂时给定，以后改在ReadMe里面,就是eta_l
    INTEGER::FL_Params_ne
    INTEGER::FL_Params_n,Fixed(2)
    REAL(8)::FL_Params_L
    REAL(8)::FL_Params_rho
    REAL(8)::FL_Params_A   =  0.01
    REAL(8)::FL_Params_E
    REAL(8)::FL_Params_I   =  1e-7

    Real(8)::ForceXFlexi,ForceYFlexi,MomentFlexi,DampingCoeiff
    REAL(8),Allocatable::Fdis_ang(:),Fdis_Y(:),Fdis0_ang(:),Fdis0_Y(:)
    REAL(8),Allocatable::FVel_ang(:),FVel_Y(:),FVel0_ang(:),FVel0_Y(:)
    REAL(8),Allocatable::FAcc_ang(:),FAcc_Y(:),FAcc0_ang(:),FAcc0_Y(:)

    REAL(8),Allocatable::M_MASSM(:,:),K_StiffM(:,:),FEULERBER(:),Keff(:,:),Feff(:)
    REAL(8),Allocatable::Fq(:),Fqd(:),Fqdd(:),Foq(:),Foqd(:),Foqdd(:)
    Real(8),Allocatable::NSE_Alea(:,:)
    Real(8),Allocatable::X_Plate_Orig(:,:)
    Real(8),Allocatable::X_CA_ORR(:,:)

    Integer::NBN_ALE_KS
    Integer,Allocatable::BN_ALE_KS(:)
    Contains
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine lu(a,p)
    !   in situ decomposition, corresponds to LAPACK's dgebtrf
    real(8), intent(inout) :: a(:,:)
    integer, intent(out  ) :: p(:)
    integer                :: n, i,j,k,kmax
    n = size(a,1)
    p = [ ( i, i=1,n ) ]
    do k = 1,n-1
        kmax = maxloc(abs(a(p(k:),k)),1) + k-1
        if (kmax /= k ) then
            p([k, kmax]) = p([kmax, k])
            a([k, kmax],:) = a([kmax, k],:)
        end if
        a(k+1:,k) = a(k+1:,k) / a(k,k)
        forall (j=k+1:n) a(k+1:,j) = a(k+1:,j) - a(k,j)*a(k+1:,k)
    end do
    end subroutine


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    End MODULE FlexibleOnly


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine cholesky_d(n, A, G, b)
    implicit none
    integer,    intent(in)    :: n
    real*8,     intent(in)    :: A(n,n)
    real*8,     intent(out)   :: G(n,n)
    real*8, intent(inout) :: b(n)
    real*8 :: tmp
    integer    :: i,j

    ! Light check of positive definite
    do i = 1, n
        if (A(i,i).le.0.0d0) then
            b(:) = 0.0d0
            return
        end if
    end do

    ! [1]
    G(:,:)=0.0d0
    do j = 1, n
        G(j,j) = sqrt( A(j,j) - dot_product(G(j,1:j-1),G(j,1:j-1)) )
        do i = j+1, n
            G(i,j)  = ( A(i,j) - dot_product(G(i,1:j-1),G(j,1:j-1)) ) / G(j,j)
        end do
    end do
    ! write(6,'(3f10.5)') (G(i,:), i=1,n)

    ! [2]
    do i = 1, n
        tmp = 0.0d0
        do j = 1, i-1
            tmp = tmp + G(i,j)*b(j)
        end do
        b(i) = (b(i)-tmp)/G(i,i)
    end do

    ! [3]
    do i = n, 1, -1
        tmp = 0.0d0
        do j = i+1, n
            tmp = tmp + b(j)*G(j,i)
        end do
        b(i) = (b(i)-tmp)/G(i,i)
    end do
    end subroutine cholesky_d

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





    MODULE DEFINE_INTERFACE
    INTEGER::  DEBUG,L,III,NPTS,NELM,DISDOUBLE,VISDOUBLE,NUMSEGPTS(1),IMAX,JMAX,KMAX,NM(4,12),I,J,K,POSITIONCOORDSYS
    INTEGER::  ATTACHTOZONE,SCOPE,FONTTYPE,HEIGHTUNITS,ANCHOR,BOXTYPE,ISFILLED,GEOMTYPE,LINEPATTERN,NUMELLIPSEPTS
    INTEGER::  ZONE,NUMSEGMENT,BOXCOLOR,BOXFILLCOLOR,FILLCOLOR,COLOR,TEXTCOLOR,ARROWHEADSTYLE,ARROWHEADATTACHMENT
    REAL(8)::  XP,YP,ZP,FH,LINESPACING,PATTERNLENGTH,ARROWHEADANGLE,BOXMARGIN,BOXLINETHICKNESS
    REAL(8)::  ZGEOMDATA(1),TEXTANGLE,LINETHICKNESS,ARROWHEADSIZE,XGEOMDATA(1),YGEOMDATA(1)
    CHARACTER*1 NULCHAR

    INTERFACE
    INTEGER FUNCTION Tecini(Title,Variables,FName,ScratchDir,Debug,VIsDouble)

    !MS$ATTRIBUTES STDCALL :: tecini
    !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
    !MS$ATTRIBUTES REFERENCE :: ScratchDir,Debug,VIsDouble
    CHARACTER*(*) Title,Variables,FName,ScratchDir
    INTEGER::     Debug,VIsDouble
    END FUNCTION Tecini

    INTEGER FUNCTION Teczne(ZoneTitle,IMx,JMx,KMx,ZFormat,DupList)
    !MS$ATTRIBUTES STDCALL :: teczne
    !MS$ATTRIBUTES REFERENCE :: ZoneTitle,IMx,JMx,KMx
    !MS$ATTRIBUTES REFERENCE :: ZFormat,DupList
    CHARACTER*(*) ZoneTitle,ZFormat,DupList
    INTEGER    :: IMx,JMx,KMx
    END FUNCTION Teczne

    INTEGER FUNCTION Tecdat (N,FieldData,IsDouble)
    !MS$ATTRIBUTES STDCALL :: tecdat
    !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
    INTEGER:: N,IsDouble
    REAL(8):: FieldData(*)
    END FUNCTION Tecdat

    INTEGER FUNCTION Tecnod(NData)
    !MS$ATTRIBUTES STDCALL :: tecnod
    !MS$ATTRIBUTES REFERENCE :: NData
    INTEGER NData(*)
    END FUNCTION Tecnod

    INTEGER FUNCTION  Tecend
    !MS$ATTRIBUTES STDCALL :: tecend
    END FUNCTION Tecend

    END INTERFACE

    END MODULE DEFINE_INTERFACE




    SUBROUTINE Polygon_Area (NP,X,Y,S)
    IMPLICIT NONE

    INTEGER:: Np
    REAL(8):: X(NP+1),Y(NP+1),S

    INTEGER:: J
    REAL(8):: XX2,YY2,XX1,YY1

    S=0.0D0
    DO J=1,NP
        XX1=X(J)
        YY1=Y(J)
        XX2=X(J+1)
        YY2=Y(J+1)
        S=S+(XX1*YY2-XX2*YY1)
    END DO

    IF (S>0.0D0) THEN
        S=S*0.5D0
    ELSE
        WRITE(*,*) ' pmccy'
        write(*,*) ' Err Subroutine_Polygon_Area'
        PAUSE
        STOP
    END IF

    END SUBROUTINE Polygon_Area



    subroutine Dsvrgp(n,ra,rb,iperm)
    implicit none
    integer:: n,iperm(n)
    real(8):: ra(n),rb(n)

    integer:: i,j,ntp
    real(8):: temp

    do i=n-1,1,-1
        do j=1,i
            if(ra(j)>ra(j+1)) then
                temp=ra(j)
                ra(j)=ra(j+1)
                ra(j+1)=temp

                ntp=iperm(j)
                iperm(j)=iperm(j+1)
                iperm(j+1)=ntp
            end if
        end do
    end do
    rb=ra       ! Caution
    end subroutine dsvrgp

    !
    !
    !
    !
    !
    !
    SUBROUTINE NXNY(X1,Y1,X2,Y2,FNX,FNY)
    IMPLICIT NONE
    REAL(8):: X1,Y1,X2,Y2,DS,FNX,FNY

    DS=DSQRT((X2-X1)**2.0D0+(Y2-Y1)**2.0D0) !两点间的距离
    FNX=(Y2-Y1)/DS
    FNY=(X1-X2)/DS

    END SUBROUTINE NXNY

    SUBROUTINE KNXNY(X1,Y1,X2,Y2,FNX,FNY)
    IMPLICIT NONE
    REAL(8):: X1,Y1,X2,Y2,DS,FNX,FNY

    DS=DSQRT((X2-X1)**2.0D0+(Y2-Y1)**2.0D0) !两点间的距离
    FNX=(Y1-Y2)/DS
    FNY=(X1-X2)/DS

    END SUBROUTINE KNXNY

    !
    !
    !
    !
    !
    SUBROUTINE TXTY(X1,Y1,X2,Y2,FTX,FTY)
    IMPLICIT NONE
    REAL(8):: X1,Y1,X2,Y2,DS,FTX,FTY

    DS=DSQRT((X2-X1)**2.0D0+(Y2-Y1)**2.0D0)
    FTX=(X2-X1)/DS
    FTY=(Y2-Y1)/DS

    END SUBROUTINE TXTY





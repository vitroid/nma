      program nmae
      implicit real*8(a-h,o-z)
      CSINTH=1.0d0
      IPOT=4
      !read(5,*) maxit
      !read(5,*) in,in
      io=in+1
      call nma(csinth,ipot,in,io)
      rewind in
      rewind io
      stop
      end      
      subroutine NMA(csinth,ipot,in,io)
C ----------------------------------------------------------------------
C
C       NORMAL MODE ANALYSIS FOR WATER MOLECULES
C       TO THREE TRANSLATIONAL MODES AND THREE ROTATINAL MODES
C                          CODED BY H. TANAKA      1988, 01, 08
C                               SEE J.CHEM.PHYS. 87, 6070 (1987)
C       FIRST AND SECOND DERIVATIVES ARE EXPRESSED IN ANALYTICAL FORM
C       EWALD SUM IS NOT TAKEN INTO ACCOUNT
C ----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LPq=576,NPq=Lpq)
      PARAMETER (NNPq=NPq*(NPq-1)/2,Mq=6,NDq=Mq*NPq,NRq=3)
      PARAMETER (NRRq=NRq*(NRq+1)/2,NR2q=2*NRq,ND2q=NDq*2)
      PARAMETER (NDPq=NDq*(NDq+1)/2)
      COMMON /MT/ E(NRq),EV(NRq,NRq),VW(NDq,2),T(NRq,NRq),
     *            EKV(NDq),EKM(NRq,NRq,NPq),S(NRq,NRq)
      COMMON /AR/ X(LPq),Y(LPq),Z(LPq),
     *            EA(LPq),EB(LPq),EC(LPq),
     *            XX1(NPq),YY1(NPq),ZZ1(NPq),
     *            XX2(NPq),YY2(NPq),ZZ2(NPq),
     *            XX3(NPq),YY3(NPq),ZZ3(NPq),
     *            XX4(NPq),YY4(NPq),ZZ4(NPq)
      COMMON /BR/ AX1(NPq),AX2(NPq),AX3(NPq),AX4(NPq),
     *          AY1(NPq),AY2(NPq),AY3(NPq),AY4(NPq),
     *          AZ1(NPq),AZ2(NPq),AZ3(NPq),AZ4(NPq),
     *          BX1(NPq),BX2(NPq),BX3(NPq),BX4(NPq),
     *          BY1(NPq),BY2(NPq),BY3(NPq),BY4(NPq),
     *          CX1(NPq),CX2(NPq),CX3(NPq),CX4(NPq),
     *          CY1(NPq),CY2(NPq),CY3(NPq),CY4(NPq),
     *          CZ1(NPq),CZ2(NPq),CZ3(NPq),CZ4(NPq)
      COMMON /BR/ AAX1(NPq),AAX2(NPq),AAX3(NPq),AAX4(NPq),
     *          ABX1(NPq),ABX2(NPq),ABX3(NPq),ABX4(NPq),
     *          ACX1(NPq),ACX2(NPq),ACX3(NPq),ACX4(NPq),
     *          BBX1(NPq),BBX2(NPq),BBX3(NPq),BBX4(NPq),
     *          BCX1(NPq),BCX2(NPq),BCX3(NPq),BCX4(NPq),
     *          CCX1(NPq),CCX2(NPq),CCX3(NPq),CCX4(NPq)
      COMMON /BR/ AAY1(NPq),AAY2(NPq),AAY3(NPq),AAY4(NPq),
     *          ABY1(NPq),ABY2(NPq),ABY3(NPq),ABY4(NPq),
     *          ACY1(NPq),ACY2(NPq),ACY3(NPq),ACY4(NPq),
     *          BBY1(NPq),BBY2(NPq),BBY3(NPq),BBY4(NPq),
     *          BCY1(NPq),BCY2(NPq),BCY3(NPq),BCY4(NPq),
     *          CCY1(NPq),CCY2(NPq),CCY3(NPq),CCY4(NPq)
      COMMON /BR/ AAZ1(NPq),AAZ2(NPq),AAZ3(NPq),AAZ4(NPq),
     *          ACZ1(NPq),ACZ2(NPq),ACZ3(NPq),ACZ4(NPq),
     *          CCZ1(NPq),CCZ2(NPq),CCZ3(NPq),CCZ4(NPq)
      COMMON /CR/ FDI(NDq),SD(NDq,NDq)
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6
      COMMON /DR/ BYL,BZL
      COMMON /IN/ LC(NPq),N,NN,N3,N6
      DATA BK,ANUM/1.380658d-16,6.0221367d+23/
      DATA ANGLE,HOLEN,COLEN/104.52D0,0.9572D0,0.1546D0/
      DATA QE/332.17752D0/
      DATA UJ/4.184D3/,HCP/6.62618D-34/
      DATA SG,WM/1.0D-8,18.0D0/
      ioi=100
      EPS=2.0D-16
      QQ=0.5564d0
      AD1=7.313907049d5
      AD2=7.360955117d2
      PI=4.0D0*ATAN(1.0D0)
      NTHREE=NRq
      QS=QQ**2*QE
      SW=SQRT(QS*UJ*1.0D3/18.0D0)/(6.0D0*PI)
      P1=AD1/QS
      P2=-AD2/QS
      P3=-P1*12.0D0
      print *,"p1 p2",p1,p2
      P4=-P2*6.0D0
      P5=-P3*13.0D0
      P6=-P4*7.0D0
      RANGLE=PI*ANGLE/360.0D0
      OHZ=HOLEN*COS(RANGLE)
      HYL=HOLEN*SIN(RANGLE)
      HZL=16.0D0*OHZ/18.0D0
      OL=-OHZ+HZL
      CL=COLEN+OL
      wnm=15.0d0
      read(5,*) n,dtp,temp,bxl,byl,bzl,rc,ndata   
      write(io) n,dtp,temp,bxl,byl,bzl,rc,ndata   
      WMM=WM/ANUM
      ENU=QS/ANUM*UJ*1.0D+7
      ETEMP=ENU/BK
      NN=N-1
      N3=N*3
      N6=N3+N3
      N63=N3+N3-3
      RIX=(16.0D0*OL*OL+2.0D0*(HYL*HYL+HZL*HZL))/WM
      RIY=(16.0D0*OL*OL+2.0D0*HZL*HZL)/WM
      RIZ=2.0D0*HYL*HYL/WM
      GC=ANUM*BK/1.0D10
      PN=FLOAT(N)
      PNI=1.0D0/PN
      EMU=1.0D-10*ENU*ANUM*PNI
      EMUPID=EMU*PN
      RC2=RC*RC
      RL=RC-2.0D0
      RQ=1.0D0/(RL-RC)**5
      RQ3=30.0D0/(RL-RC)**5
      RQ6=60.0D0/(RL-RC)**5
      PID=2.0D0*PI
      WRITE(6,601) mmm,N,volm
601   FORMAT(5X,i7,1x,'NORMAL MODE ANALYSIS FOR WATER '/
     *          5X,'NUMBER OF MOLECULES = ',I5,1x,f10.5)
      do lll=1,ndata
         do i=1,n
 1015       READ(5,*) X(i),Y(i),Z(i),EA(i),EB(i),EC(i)
            print *, lll,i,X(i),y(i),z(i)
         enddo
      vol=bxl*byl*bzl
      volm=vol*anum/(pn*1.0d24)
      DO 48 I=1,N
      TH=EA(I)
      PH=EB(I)
      PS=EC(I)
c     a(i)=ea(i)
c     b(i)=eb(i)
c     c(i)=ec(i)
      SINA=SIN(TH)
      SINB=SIN(PH)
      SINC=SIN(PS)
      COSA=COS(TH)
      COSB=COS(PH)
      COSC=COS(PS)
      COSAP=SINA*SINC
      IF(ABS(SINA).GT.CSINTH)  THEN
      LC(I)=1
      CALL LC1(SINA,SINB,SINC,COSA,COSB,COSC,HYL,HZL,OL,CL,I)
      ELSE
      LC(I)=2
      AP=ACOS(COSAP)
      SINAP=SIN(AP)
      COSCP=COSA/SINAP
      CP=ACOS(COSCP)
      SINCP=SINA*COSC/SINAP
      IF(SINCP.LT.0.0D0) THEN
      CP=PID-CP
      END IF
      COSBP=-(COSC*SINB+COSA*COSB*SINC)/SINAP
      BP=ACOS(COSBP)
      SINBP=(COSC*COSB-COSA*SINB*SINC)/SINAP
      IF(SINBP.LT.0.0D0) THEN
      BP=PID-BP
      END IF
      EA(I)=AP
      EB(I)=BP
      EC(I)=CP
      CALL LC2(SINAP,SINBP,SINCP,COSAP,COSBP,COSCP,HYL,HZL,OL,CL,I)
      END IF
48    CONTINUE
c      WRITE(6,1122)(LC(I),I=1,N)
1122  FORMAT(65I2)
      DO 148 I=1,N
      TH=EA(I)
      PH=EB(I)
      PS=EC(I)
      SINA=SIN(TH)
      SINB=SIN(PH)
      SINC=SIN(PS)
      COSA=COS(TH)
      COSB=COS(PH)
      COSC=COS(PS)
      IF(LC(I).EQ.1) THEN
      EV(1,1)=RIX*COSC**2+RIY*SINC**2
      EV(1,2)=(RIX-RIY)*SINA*SINC*COSC
      EV(1,3)=0.0D0
      EV(2,2)=(RIX*SINC**2+RIY*COSC**2)*SINA**2+RIZ*COSA**2
      EV(2,3)=RIZ*COSA
      EV(3,3)=RIZ
      EV(2,1)=EV(1,2)
      EV(3,1)=EV(1,3)
      EV(3,2)=EV(2,3)
      ELSE
      EV(1,1)=RIY*COSC**2+RIZ*SINC**2
      EV(1,2)=(RIY-RIZ)*SINA*SINC*COSC
      EV(1,3)=0.0D0
      EV(2,2)=(RIY*SINC**2+RIZ*COSC**2)*SINA**2+RIX*COSA**2
      EV(2,3)=RIX*COSA
      EV(3,3)=RIX
      EV(2,1)=EV(1,2)
      EV(3,1)=EV(1,3)
      EV(3,2)=EV(2,3)
      END IF
      ICON=1
      CALL HOQRVW(EV,NRq,NTHREE,E,EPS,VW,ICON)
      J1=3*I-2
      J2=J1+1
      J3=J2+1
      EKV(J1)=1.0D0/SQRT(E(1))
      EKV(J2)=1.0D0/SQRT(E(2))
      EKV(J3)=1.0D0/SQRT(E(3))
      EKM(1,1,I)=EV(1,1)
      EKM(2,1,I)=EV(2,1)
      EKM(3,1,I)=EV(3,1)
      EKM(1,2,I)=EV(1,2)
      EKM(2,2,I)=EV(2,2)
      EKM(3,2,I)=EV(3,2)
      EKM(1,3,I)=EV(1,3)
      EKM(2,3,I)=EV(2,3)
      EKM(3,3,I)=EV(3,3)
  148 CONTINUE
      DO 58 I=1,NDq
      DO 58 J=1,NDq
      SD(J,I)=0.0D0
58    CONTINUE
      CALL DRVTV
      EP=EP/FLOAT(N)
      WRITE(6,*) lll,' ENERGY = ',EP,epc
      DO 777 I=1,N
      I1=3*I-2
      I2=I1+1
      I3=I2+1
      K1=N3+I1
      K2=K1+1
      K3=K2+1
      DO 777 J=1,N
      J1=3*J-2
      J2=J1+1
      J3=J2+1
      L1=N3+J1
      L2=L1+1
      L3=L2+1
      T(1,1)=SD(L1,K1)
      T(2,1)=SD(L2,K1)
      T(3,1)=SD(L3,K1)
      T(1,2)=SD(L1,K2)
      T(2,2)=SD(L2,K2)
      T(3,2)=SD(L3,K2)
      T(1,3)=SD(L1,K3)
      T(2,3)=SD(L2,K3)
      T(3,3)=SD(L3,K3)
      DO 779 K=1,NRq
      DO 779 L=1,NRq
      SS=0.0D0
      DO 778 M=1,NRq
      SS=SS+T(L,M)*EKM(M,K,I)
778   CONTINUE
      S(L,K)=SS
779   CONTINUE
      DO 769 K=1,NRq
      DO 769 L=1,NRq
      SS=0.0D0
      DO 768 M=1,NRq
      SS=SS+S(M,K)*EKM(M,L,J)
768   CONTINUE
      T(L,K)=SS
769   CONTINUE
      SD(L1,K1)=T(1,1)*EKV(J1)*EKV(I1)
      SD(L2,K1)=T(2,1)*EKV(J2)*EKV(I1)
      SD(L3,K1)=T(3,1)*EKV(J3)*EKV(I1)
      SD(L1,K2)=T(1,2)*EKV(J1)*EKV(I2)
      SD(L2,K2)=T(2,2)*EKV(J2)*EKV(I2)
      SD(L3,K2)=T(3,2)*EKV(J3)*EKV(I2)
      SD(L1,K3)=T(1,3)*EKV(J1)*EKV(I3)
      SD(L2,K3)=T(2,3)*EKV(J2)*EKV(I3)
      SD(L3,K3)=T(3,3)*EKV(J3)*EKV(I3)
777   CONTINUE
      DO 877 I=1,N
      I1=3*I-2
      I2=I1+1
      I3=I2+1
      DO 877 J=1,N
      J1=3*J-2
      J2=J1+1
      J3=J2+1
      K1=N3+J1
      K2=K1+1
      K3=K2+1
      T(1,1)=SD(K1,I1)
      T(2,1)=SD(K2,I1)
      T(3,1)=SD(K3,I1)
      T(1,2)=SD(K1,I2)
      T(2,2)=SD(K2,I2)
      T(3,2)=SD(K3,I2)
      T(1,3)=SD(K1,I3)
      T(2,3)=SD(K2,I3)
      T(3,3)=SD(K3,I3)
      DO 869 K=1,NRq
      DO 869 L=1,NRq
      SS=0.0D0
      DO 868 M=1,NRq
      SS=SS+EKM(M,L,J)*T(M,K)
868   CONTINUE
      S(L,K)=SS
869   CONTINUE
       SD(K1,I1)=S(1,1)*EKV(J1)
       SD(K2,I1)=S(2,1)*EKV(J2)
       SD(K3,I1)=S(3,1)*EKV(J3)
       SD(K1,I2)=S(1,2)*EKV(J1)
       SD(K2,I2)=S(2,2)*EKV(J2)
       SD(K3,I2)=S(3,2)*EKV(J3)
       SD(K1,I3)=S(1,3)*EKV(J1)
       SD(K2,I3)=S(2,3)*EKV(J2)
       SD(K3,I3)=S(3,3)*EKV(J3)
877   CONTINUE
      DO 888 I=1,N6
      DO 889 J=1,I
      SD(J,I)=SD(I,J)
889   CONTINUE
888   CONTINUE
C     *        CALL HOQRVW(A,KA,N,E,EPS,W,IND)                         *
C     *          A   ....GIVEN REAL SYMMETRIC MATRIX.                  *
C     *          KA  ....GIVEN ADJUSTABLE DIMENSION OF A.              *
C     *          N   ....GIVEN ORDER OF THE MATRIX A.                  *
C     *          E   ....RESULTANT EIGENVALUES.                        *
C     *          EPS ....GIVEN CRITERION FOR THE CONVERGENCE TEST.     *
C     *          W   ....WORK AREA.                                    *
C     *          IND ....GIVEN INDEX.                                  *
C     *             0          ....DO NOT COMPUTE EIGENVECTORS.        *
C     *             1          ....COMPUTE ALL THE EIGENVECTORS.       *
C     *              ....RESULTANT CONDITION CODE.                     *
C     *             0          ....NORMAL TERMINATION.                 *
C     *             30000      ....PARAMETER ERROR.                    *
      ICON=0
      CALL HOQRVW(SD,NDq,N6,FDI,EPS,VW,ICON)
C     CALL DSEIG1(SDS,N6,FDI,SD,NDq,NE,VW,ICON)
      IF(ICON.NE.0) THEN
      WRITE(6,*) 'ERROR  IN NORMAL MODE ANALYSIS ERROR CODE = ',ICON
      ELSE
      WRITE(6,610)
610   FORMAT(2X,'EIGEN VALUE IN CM**(-1) ')
      idd=0
      DO 772 I=1,N6
      IF(FDI(I).GT.0.0D0) THEN
      FDI(I)=SQRT(FDI(I))*SW
      ELSE
      FDI(I)=-SQRT(-FDI(I))*SW
      END IF
     
772   CONTINUE
      if(fdi(n63).lt.wnm) then
      idd=1
      write(9,901) idd,lll,mmm,fdi(n63)
      write(8,901) idd,lll,mmm,fdi(n63)
901   format(1x,i2,2(1x,i5),1x,f10.5)
      else
      write(9,901) idd,lll,mmm,fdi(n63)
      end if
      do i=1,N6
         print *, i, FDI(i)
      enddo
      ! WRITE(6,611)(I,FDI(I),I=1,N6)
CCCC  WRITE(6,611)((I,FDI(I)),I=N3+1,N6)
611   FORMAT(5(1X,I5,1x,f7.2))
      write(io) (fdi(i),i=1,n6),ep,epc,volm
c      N63=N6-3
c      WRITE(io) FDI,SD,LC
      END IF
	end do
      return
      END
      SUBROUTINE LC1(SINA,SINB,SINC,COSA,COSB,COSC,HYL,HZL,OL,CL,I)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LPq=576,NPq=Lpq)
      PARAMETER (NNPq=NPq*(NPq-1)/2,Mq=6,NDq=Mq*NPq,NRq=3)
      PARAMETER (NRRq=NRq*(NRq+1)/2,NR2q=2*NRq,ND2q=NDq*2)
      PARAMETER (NDPq=NDq*(NDq+1)/2)
      COMMON /MT/ E(NRq),EV(NRq,NRq),VW(NDq,2),T(NRq,NRq),
     *            EKV(NDq),EKM(NRq,NRq,NPq),S(NRq,NRq)
      COMMON /AR/ X(LPq),Y(LPq),Z(LPq),
     *            EA(LPq),EB(LPq),EC(LPq),
     *            XX1(NPq),YY1(NPq),ZZ1(NPq),
     *            XX2(NPq),YY2(NPq),ZZ2(NPq),
     *            XX3(NPq),YY3(NPq),ZZ3(NPq),
     *            XX4(NPq),YY4(NPq),ZZ4(NPq)
      COMMON /BR/ AX1(NPq),AX2(NPq),AX3(NPq),AX4(NPq),
     *          AY1(NPq),AY2(NPq),AY3(NPq),AY4(NPq),
     *          AZ1(NPq),AZ2(NPq),AZ3(NPq),AZ4(NPq),
     *          BX1(NPq),BX2(NPq),BX3(NPq),BX4(NPq),
     *          BY1(NPq),BY2(NPq),BY3(NPq),BY4(NPq),
     *          CX1(NPq),CX2(NPq),CX3(NPq),CX4(NPq),
     *          CY1(NPq),CY2(NPq),CY3(NPq),CY4(NPq),
     *          CZ1(NPq),CZ2(NPq),CZ3(NPq),CZ4(NPq)
      COMMON /BR/ AAX1(NPq),AAX2(NPq),AAX3(NPq),AAX4(NPq),
     *          ABX1(NPq),ABX2(NPq),ABX3(NPq),ABX4(NPq),
     *          ACX1(NPq),ACX2(NPq),ACX3(NPq),ACX4(NPq),
     *          BBX1(NPq),BBX2(NPq),BBX3(NPq),BBX4(NPq),
     *          BCX1(NPq),BCX2(NPq),BCX3(NPq),BCX4(NPq),
     *          CCX1(NPq),CCX2(NPq),CCX3(NPq),CCX4(NPq)
      COMMON /BR/ AAY1(NPq),AAY2(NPq),AAY3(NPq),AAY4(NPq),
     *          ABY1(NPq),ABY2(NPq),ABY3(NPq),ABY4(NPq),
     *          ACY1(NPq),ACY2(NPq),ACY3(NPq),ACY4(NPq),
     *          BBY1(NPq),BBY2(NPq),BBY3(NPq),BBY4(NPq),
     *          BCY1(NPq),BCY2(NPq),BCY3(NPq),BCY4(NPq),
     *          CCY1(NPq),CCY2(NPq),CCY3(NPq),CCY4(NPq)
      COMMON /BR/ AAZ1(NPq),AAZ2(NPq),AAZ3(NPq),AAZ4(NPq),
C    *          ABZ1(NPq),ABZ2(NPq),ABZ3(NPq),ABZ4(NPq),
     *          ACZ1(NPq),ACZ2(NPq),ACZ3(NPq),ACZ4(NPq),
C    *          BBZ1(NPq),BBZ2(NPq),BBZ3(NPq),BBZ4(NPq),
C    *          BCZ1(NPq),BCZ2(NPq),BCZ3(NPq),BCZ4(NPq),
     *          CCZ1(NPq),CCZ2(NPq),CCZ3(NPq),CCZ4(NPq)
      COMMON /CR/ FDI(NDq),SD(NDq,NDq)
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6
      COMMON /DR/ BYL,BZL
      COMMON /IN/ LC(NPq),N,NN,N3,N6
      SP12=-(SINC*COSB+COSA*SINB*COSC)
      DDA12=SINA*SINB*COSC
      DDB12=SINC*SINB-COSA*COSB*COSC
      DDC12=-COSC*COSB+COSA*SINB*SINC
      DAA12=COSA*SINB*COSC
      DAB12=SINA*COSB*COSC
      DAC12=-SINA*SINB*SINC
      DBB12=-SP12
      DBC12=COSC*SINB+COSA*COSB*SINC
      DCC12=-SP12
      SP13=SINA*SINB
      DDA13=COSA*SINB
      DDB13=SINA*COSB
      DAA13=-SP13
      DAB13=COSA*COSB
      DBB13=-SP13
      SP22=-SINC*SINB+COSA*COSB*COSC
      DDA22=-SINA*COSB*COSC
      DDB22=-SINC*COSB-COSA*SINB*COSC
      DDC22=-COSC*SINB-COSA*COSB*SINC
      DAA22=-COSA*COSB*COSC
      DAB22=SINA*SINB*COSC
      DAC22=SINA*COSB*SINC
      DBB22=-SP22
      DBC22=-COSC*COSB+COSA*SINB*SINC
      DCC22=-SP22
      SP23=-SINA*COSB
      DDA23=-COSA*COSB
      DDB23=SINA*SINB
      DAA23=-SP23
      DAB23=COSA*SINB
      DBB23=-SP23
      SP32=SINA*COSC
      DDA32=COSA*COSC
      DDC32=-SINA*SINC
      DAA32=-SP32
      DAC32=-COSA*SINC
      DCC32=-SP32
      SP33=COSA
      DDA33=-SINA
      DAA33=-COSA
      SY12=SP12*HYL
      SY22=SP22*HYL
      SY32=SP32*HYL
      SZ13=SP13*HZL
      SZ23=SP23*HZL
      SZ33=SP33*HZL
      XX1(I)=SY12+SZ13
      YY1(I)=SY22+SZ23
      ZZ1(I)=SY32+SZ33
      XX2(I)=-SY12+SZ13
      YY2(I)=-SY22+SZ23
      ZZ2(I)=-SY32+SZ33
      XX3(I)=SP13*OL
      YY3(I)=SP23*OL
      ZZ3(I)=SP33*OL
      XX4(I)=SP13*CL
      YY4(I)=SP23*CL
      ZZ4(I)=SP33*CL
      A12H=DDA12*HYL
      A22H=DDA22*HYL
      A32H=DDA32*HYL
      A13H=DDA13*HZL
      A23H=DDA23*HZL
      A33H=DDA33*HZL
      B12H=DDB12*HYL
      B22H=DDB22*HYL
      B13H=DDB13*HZL
      B23H=DDB23*HZL
      C12H=DDC12*HYL
      C22H=DDC22*HYL
      C32H=DDC32*HYL
      A13C=DDA13*CL
      A23C=DDA23*CL
      A33C=DDA33*CL
      B13C=DDB13*CL
      B23C=DDB23*CL
      A13O=DDA13*OL
      A23O=DDA23*OL
      A33O=DDA33*OL
      B13O=DDB13*OL
      B23O=DDB23*OL
      AA12H=DAA12*HYL
      AA22H=DAA22*HYL
      AA32H=DAA32*HYL
      AA13H=DAA13*HZL
      AA23H=DAA23*HZL
      AA33H=DAA33*HZL
      AB12H=DAB12*HYL
      AB22H=DAB22*HYL
      AB13H=DAB13*HZL
      AB23H=DAB23*HZL
      AC12H=DAC12*HYL
      AC22H=DAC22*HYL
      AC32H=DAC32*HYL
      BB12H=DBB12*HYL
      BB22H=DBB22*HYL
      BB13H=DBB13*HZL
      BB23H=DBB23*HZL
      BC12H=DBC12*HYL
      BC22H=DBC22*HYL
      CC12H=DCC12*HYL
      CC22H=DCC22*HYL
      CC32H=DCC32*HYL
      AA13C=DAA13*CL
      AA23C=DAA23*CL
      AA33C=DAA33*CL
      AB13C=DAB13*CL
      AB23C=DAB23*CL
      BB13C=DBB13*CL
      BB23C=DBB23*CL
      AA13O=DAA13*OL
      AA23O=DAA23*OL
      AA33O=DAA33*OL
      AB13O=DAB13*OL
      AB23O=DAB23*OL
      BB13O=DBB13*OL
      BB23O=DBB23*OL
      AX1(I)=A12H+A13H
      AY1(I)=A22H+A23H
      AZ1(I)=A32H+A33H
      AX2(I)=-A12H+A13H
      AY2(I)=-A22H+A23H
      AZ2(I)=-A32H+A33H
      AX3(I)=A13O
      AY3(I)=A23O
      AZ3(I)=A33O
      AX4(I)=A13C
      AY4(I)=A23C
      AZ4(I)=A33C
      BX1(I)=B12H+B13H
      BY1(I)=B22H+B23H
      BX2(I)=-B12H+B13H
      BY2(I)=-B22H+B23H
      BX3(I)=B13O
      BY3(I)=B23O
      BX4(I)=B13C
      BY4(I)=B23C
      CX1(I)=C12H
      CY1(I)=C22H
      CZ1(I)=C32H
      CX2(I)=-C12H
      CY2(I)=-C22H
      CZ2(I)=-C32H
      CX3(I)=0.0D0
      CY3(I)=0.0D0
      CZ3(I)=0.0D0
      CX4(I)=0.0D0
      CY4(I)=0.0D0
      CZ4(I)=0.0D0
      AAX1(I)=AA12H+AA13H
      AAY1(I)=AA22H+AA23H
      AAZ1(I)=AA32H+AA33H
      AAX2(I)=-AA12H+AA13H
      AAY2(I)=-AA22H+AA23H
      AAZ2(I)=-AA32H+AA33H
      AAX3(I)=AA13O
      AAY3(I)=AA23O
      AAZ3(I)=AA33O
      AAX4(I)=AA13C
      AAY4(I)=AA23C
      AAZ4(I)=AA33C
      BBX1(I)=BB12H+BB13H
      BBY1(I)=BB22H+BB23H
      BBX2(I)=-BB12H+BB13H
      BBY2(I)=-BB22H+BB23H
      BBX3(I)=BB13O
      BBY3(I)=BB23O
      BBX4(I)=BB13C
      BBY4(I)=BB23C
      CCX1(I)=CC12H
      CCY1(I)=CC22H
      CCZ1(I)=CC32H
      CCX2(I)=-CC12H
      CCY2(I)=-CC22H
      CCZ2(I)=-CC32H
      CCX3(I)=0.0D0
      CCY3(I)=0.0D0
      CCZ3(I)=0.0D0
      CCX4(I)=0.0D0
      CCY4(I)=0.0D0
      CCZ4(I)=0.0D0
      ABX1(I)=AB12H+AB13H
      ABY1(I)=AB22H+AB23H
      ABX2(I)=-AB12H+AB13H
      ABY2(I)=-AB22H+AB23H
      ABX3(I)=AB13O
      ABY3(I)=AB23O
      ABX4(I)=AB13C
      ABY4(I)=AB23C
      ACX1(I)=AC12H
      ACY1(I)=AC22H
      ACZ1(I)=AC32H
      ACX2(I)=-AC12H
      ACY2(I)=-AC22H
      ACZ2(I)=-AC32H
      ACX3(I)=0.0D0
      ACY3(I)=0.0D0
      ACZ3(I)=0.0D0
      ACX4(I)=0.0D0
      ACY4(I)=0.0D0
      ACZ4(I)=0.0D0
      BCX1(I)=BC12H
      BCY1(I)=BC22H
      BCX2(I)=-BC12H
      BCY2(I)=-BC22H
      BCX3(I)=0.0D0
      BCY3(I)=0.0D0
      BCX4(I)=0.0D0
      BCY4(I)=0.0D0
      RETURN
      END
      SUBROUTINE LC2(SINAP,SINBP,SINCP,COSAP,COSBP,COSCP,
     *               HYL,HZL,OL,CL,I)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LPq=576,NPq=Lpq)
      PARAMETER (NNPq=NPq*(NPq-1)/2,Mq=6,NDq=Mq*NPq,NRq=3)
      PARAMETER (NRRq=NRq*(NRq+1)/2,NR2q=2*NRq,ND2q=NDq*2)
      PARAMETER (NDPq=NDq*(NDq+1)/2)
      COMMON /MT/ E(NRq),EV(NRq,NRq),VW(NDq,2),T(NRq,NRq),
     *            EKV(NDq),EKM(NRq,NRq,NPq),S(NRq,NRq)
      COMMON /AR/ X(LPq),Y(LPq),Z(LPq),
     *            EA(LPq),EB(LPq),EC(LPq),
     *            XX1(NPq),YY1(NPq),ZZ1(NPq),
     *            XX2(NPq),YY2(NPq),ZZ2(NPq),
     *            XX3(NPq),YY3(NPq),ZZ3(NPq),
     *            XX4(NPq),YY4(NPq),ZZ4(NPq)
      COMMON /BR/ AX1(NPq),AX2(NPq),AX3(NPq),AX4(NPq),
     *          AY1(NPq),AY2(NPq),AY3(NPq),AY4(NPq),
     *          AZ1(NPq),AZ2(NPq),AZ3(NPq),AZ4(NPq),
     *          BX1(NPq),BX2(NPq),BX3(NPq),BX4(NPq),
     *          BY1(NPq),BY2(NPq),BY3(NPq),BY4(NPq),
     *          CX1(NPq),CX2(NPq),CX3(NPq),CX4(NPq),
     *          CY1(NPq),CY2(NPq),CY3(NPq),CY4(NPq),
     *          CZ1(NPq),CZ2(NPq),CZ3(NPq),CZ4(NPq)
      COMMON /BR/ AAX1(NPq),AAX2(NPq),AAX3(NPq),AAX4(NPq),
     *          ABX1(NPq),ABX2(NPq),ABX3(NPq),ABX4(NPq),
     *          ACX1(NPq),ACX2(NPq),ACX3(NPq),ACX4(NPq),
     *          BBX1(NPq),BBX2(NPq),BBX3(NPq),BBX4(NPq),
     *          BCX1(NPq),BCX2(NPq),BCX3(NPq),BCX4(NPq),
     *          CCX1(NPq),CCX2(NPq),CCX3(NPq),CCX4(NPq)
      COMMON /BR/ AAY1(NPq),AAY2(NPq),AAY3(NPq),AAY4(NPq),
     *          ABY1(NPq),ABY2(NPq),ABY3(NPq),ABY4(NPq),
     *          ACY1(NPq),ACY2(NPq),ACY3(NPq),ACY4(NPq),
     *          BBY1(NPq),BBY2(NPq),BBY3(NPq),BBY4(NPq),
     *          BCY1(NPq),BCY2(NPq),BCY3(NPq),BCY4(NPq),
     *          CCY1(NPq),CCY2(NPq),CCY3(NPq),CCY4(NPq)
      COMMON /BR/ AAZ1(NPq),AAZ2(NPq),AAZ3(NPq),AAZ4(NPq),
C    *          ABZ1(NPq),ABZ2(NPq),ABZ3(NPq),ABZ4(NPq),
     *          ACZ1(NPq),ACZ2(NPq),ACZ3(NPq),ACZ4(NPq),
C    *          BBZ1(NPq),BBZ2(NPq),BBZ3(NPq),BBZ4(NPq),
C    *          BCZ1(NPq),BCZ2(NPq),BCZ3(NPq),BCZ4(NPq),
     *          CCZ1(NPq),CCZ2(NPq),CCZ3(NPq),CCZ4(NPq)
      COMMON /CR/ FDI(NDq),SD(NDq,NDq)
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6
      COMMON /DR/ BYL,BZL
      COMMON /IN/ LC(NPq),N,NN,N3,N6
      SP12=COSCP*COSBP-COSAP*SINBP*SINCP
      DDA12=SINAP*SINBP*SINCP
      DDB12=-COSCP*SINBP-COSAP*COSBP*SINCP
      DDC12=-SINCP*COSBP-COSAP*SINBP*COSCP
      DAA12=COSAP*SINBP*SINCP
      DAB12=SINAP*COSBP*SINCP
      DBB12=-SP12
      DAC12=SINAP*SINBP*COSCP
      DBC12=SINCP*SINBP-COSAP*COSBP*COSCP
      DCC12=-SP12
      SP13=-(SINCP*COSBP+COSAP*SINBP*COSCP)
      DDA13=SINAP*SINBP*COSCP
      DDB13=SINCP*SINBP-COSAP*COSBP*COSCP
      DDC13=-COSCP*COSBP+COSAP*SINBP*SINCP
      DAA13=COSAP*SINBP*COSCP
      DAB13=SINAP*COSBP*COSCP
      DAC13=-SINAP*SINBP*SINCP
      DBB13=-SP13
      DBC13=COSCP*SINBP+COSAP*COSBP*SINCP
      DCC13=-SP13
      SP22=COSCP*SINBP+COSAP*COSBP*SINCP
      DDA22=-SINAP*COSBP*SINCP
      DDB22=COSCP*COSBP-COSAP*SINBP*SINCP
      DDC22=-SINCP*SINBP+COSAP*COSBP*COSCP
      DAA22=-COSAP*COSBP*SINCP
      DAB22=SINAP*SINBP*SINCP
      DAC22=-SINAP*COSBP*COSCP
      DBB22=-SP22
      DCC22=-SP22
      DBC22=-SINCP*COSBP-COSAP*SINBP*COSCP
      SP23=-SINCP*SINBP+COSAP*COSBP*COSCP
      DDA23=-SINAP*COSBP*COSCP
      DDB23=-SINCP*COSBP-COSAP*SINBP*COSCP
      DDC23=-COSCP*SINBP-COSAP*COSBP*SINCP
      DAA23=-COSAP*COSBP*COSCP
      DAB23=SINAP*SINBP*COSCP
      DAC23=SINAP*COSBP*SINCP
      DBB23=-SP23
      DBC23=-COSCP*COSBP+COSAP*SINBP*SINCP
      DCC23=-SP23
      SP32=SINAP*SINCP
      DDA32=COSAP*SINCP
      DDC32=SINAP*COSCP
      DAA32=-SP32
      DAC32=COSAP*COSCP
      DCC32=-SP32
      SP33=SINAP*COSCP
      DDA33=COSAP*COSCP
      DDC33=-SINAP*SINCP
      DAA33=-SP33
      DAC33=-COSAP*SINCP
      DCC33=-SP33
      SY12=SP12*HYL
      SY22=SP22*HYL
      SY32=SP32*HYL
      SZ13=SP13*HZL
      SZ23=SP23*HZL
      SZ33=SP33*HZL
      XX1(I)=SY12+SZ13
      YY1(I)=SY22+SZ23
      ZZ1(I)=SY32+SZ33
      XX2(I)=-SY12+SZ13
      YY2(I)=-SY22+SZ23
      ZZ2(I)=-SY32+SZ33
      XX3(I)=SP13*OL
      YY3(I)=SP23*OL
      ZZ3(I)=SP33*OL
      XX4(I)=SP13*CL
      YY4(I)=SP23*CL
      ZZ4(I)=SP33*CL
      A12H=DDA12*HYL
      A22H=DDA22*HYL
      A32H=DDA32*HYL
      A13H=DDA13*HZL
      A23H=DDA23*HZL
      A33H=DDA33*HZL
      B12H=DDB12*HYL
      B22H=DDB22*HYL
      B13H=DDB13*HZL
      B23H=DDB23*HZL
      C12H=DDC12*HYL
      C22H=DDC22*HYL
      C32H=DDC32*HYL
      C13H=DDC13*HZL
      C23H=DDC23*HZL
      C33H=DDC33*HZL
      A13C=DDA13*CL
      A23C=DDA23*CL
      A33C=DDA33*CL
      B13C=DDB13*CL
      B23C=DDB23*CL
      C13C=DDC13*CL
      C23C=DDC23*CL
      C33C=DDC33*CL
      A13O=DDA13*OL
      A23O=DDA23*OL
      A33O=DDA33*OL
      B13O=DDB13*OL
      B23O=DDB23*OL
      C13O=DDC13*OL
      C23O=DDC23*OL
      C33O=DDC33*OL
      AA12H=DAA12*HYL
      AA22H=DAA22*HYL
      AA32H=DAA32*HYL
      AA13H=DAA13*HZL
      AA23H=DAA23*HZL
      AA33H=DAA33*HZL
      AB12H=DAB12*HYL
      AB22H=DAB22*HYL
      AB13H=DAB13*HZL
      AB23H=DAB23*HZL
      AC12H=DAC12*HYL
      AC22H=DAC22*HYL
      AC32H=DAC32*HYL
      AC13H=DAC13*HZL
      AC23H=DAC23*HZL
      AC33H=DAC33*HZL
      BB12H=DBB12*HYL
      BB22H=DBB22*HYL
      BB13H=DBB13*HZL
      BB23H=DBB23*HZL
      BC12H=DBC12*HYL
      BC22H=DBC22*HYL
      BC13H=DBC13*HZL
      BC23H=DBC23*HZL
      CC12H=DCC12*HYL
      CC22H=DCC22*HYL
      CC32H=DCC32*HYL
      CC13H=DCC13*HZL
      CC23H=DCC23*HZL
      CC33H=DCC33*HZL
      AA13C=DAA13*CL
      AA23C=DAA23*CL
      AA33C=DAA33*CL
      AB13C=DAB13*CL
      AB23C=DAB23*CL
      AC13C=DAC13*CL
      AC23C=DAC23*CL
      AC33C=DAC33*CL
      BB13C=DBB13*CL
      BB23C=DBB23*CL
      BC13C=DBC13*CL
      BC23C=DBC23*CL
      CC13C=DCC13*CL
      CC23C=DCC23*CL
      CC33C=DCC33*CL
      AA13O=DAA13*OL
      AA23O=DAA23*OL
      AA33O=DAA33*OL
      AB13O=DAB13*OL
      AB23O=DAB23*OL
      AC13O=DAC13*OL
      AC23O=DAC23*OL
      AC33O=DAC33*OL
      BB13O=DBB13*OL
      BB23O=DBB23*OL
      BC13O=DBC13*OL
      BC23O=DBC23*OL
      CC13O=DCC13*OL
      CC23O=DCC23*OL
      CC33O=DCC33*OL
      AX1(I)=A12H+A13H
      AY1(I)=A22H+A23H
      AZ1(I)=A32H+A33H
      AX2(I)=-A12H+A13H
      AY2(I)=-A22H+A23H
      AZ2(I)=-A32H+A33H
      AX3(I)=A13O
      AY3(I)=A23O
      AZ3(I)=A33O
      AX4(I)=A13C
      AY4(I)=A23C
      AZ4(I)=A33C
      BX1(I)=B12H+B13H
      BY1(I)=B22H+B23H
      BX2(I)=-B12H+B13H
      BY2(I)=-B22H+B23H
      BX3(I)=B13O
      BY3(I)=B23O
      BX4(I)=B13C
      BY4(I)=B23C
      CX1(I)=C12H+C13H
      CY1(I)=C22H+C23H
      CZ1(I)=C32H+C33H
      CX2(I)=-C12H+C13H
      CY2(I)=-C22H+C23H
      CZ2(I)=-C32H+C33H
      CX3(I)=C13O
      CY3(I)=C23O
      CZ3(I)=C33O
      CX4(I)=C13C
      CY4(I)=C23C
      CZ4(I)=C33C
      AAX1(I)=AA12H+AA13H
      AAY1(I)=AA22H+AA23H
      AAZ1(I)=AA32H+AA33H
      AAX2(I)=-AA12H+AA13H
      AAY2(I)=-AA22H+AA23H
      AAZ2(I)=-AA32H+AA33H
      AAX3(I)=AA13O
      AAY3(I)=AA23O
      AAZ3(I)=AA33O
      AAX4(I)=AA13C
      AAY4(I)=AA23C
      AAZ4(I)=AA33C
      BBX1(I)=BB12H+BB13H
      BBY1(I)=BB22H+BB23H
      BBX2(I)=-BB12H+BB13H
      BBY2(I)=-BB22H+BB23H
      BBX3(I)=BB13O
      BBY3(I)=BB23O
      BBX4(I)=BB13C
      BBY4(I)=BB23C
      CCX1(I)=CC12H+CC13H
      CCY1(I)=CC22H+CC23H
      CCZ1(I)=CC32H+CC33H
      CCX2(I)=-CC12H+CC13H
      CCY2(I)=-CC22H+CC23H
      CCZ2(I)=-CC32H+CC33H
      CCX3(I)=CC13O
      CCY3(I)=CC23O
      CCZ3(I)=CC33O
      CCX4(I)=CC13C
      CCY4(I)=CC23C
      CCZ4(I)=CC33C
      ABX1(I)=AB12H+AB13H
      ABY1(I)=AB22H+AB23H
      ABX2(I)=-AB12H+AB13H
      ABY2(I)=-AB22H+AB23H
      ABX3(I)=AB13O
      ABY3(I)=AB23O
      ABX4(I)=AB13C
      ABY4(I)=AB23C
      ACX1(I)=AC12H+AC13H
      ACY1(I)=AC22H+AC23H
      ACZ1(I)=AC32H+AC33H
      ACX2(I)=-AC12H+AC13H
      ACY2(I)=-AC22H+AC23H
      ACZ2(I)=-AC32H+AC33H
      ACX3(I)=AC13O
      ACY3(I)=AC23O
      ACZ3(I)=AC33O
      ACX4(I)=AC13C
      ACY4(I)=AC23C
      ACZ4(I)=AC33C
      BCX1(I)=BC12H+BC13H
      BCY1(I)=BC22H+BC23H
      BCX2(I)=-BC12H+BC13H
      BCY2(I)=-BC22H+BC23H
      BCX3(I)=BC13O
      BCY3(I)=BC23O
      BCX4(I)=BC13C
      BCY4(I)=BC23C
      RETURN
      END
      SUBROUTINE DRVTV
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (LPq=576,NPq=Lpq)
      PARAMETER (NNPq=NPq*(NPq-1)/2,Mq=6,NDq=Mq*NPq,NRq=3)
      PARAMETER (NRRq=NRq*(NRq+1)/2,NR2q=2*NRq,ND2q=NDq*2)
      PARAMETER (NDPq=NDq*(NDq+1)/2)
      COMMON /MT/ E(NRq),EV(NRq,NRq),VW(NDq,2),T(NRq,NRq),
     *            EKV(NDq),EKM(NRq,NRq,NPq),S(NRq,NRq)
      COMMON /AR/ X(LPq),Y(LPq),Z(LPq),
     *            EA(LPq),EB(LPq),EC(LPq),
     *            XX1(NPq),YY1(NPq),ZZ1(NPq),
     *            XX2(NPq),YY2(NPq),ZZ2(NPq),
     *            XX3(NPq),YY3(NPq),ZZ3(NPq),
     *            XX4(NPq),YY4(NPq),ZZ4(NPq)
      COMMON /BR/ AX1(NPq),AX2(NPq),AX3(NPq),AX4(NPq),
     *          AY1(NPq),AY2(NPq),AY3(NPq),AY4(NPq),
     *          AZ1(NPq),AZ2(NPq),AZ3(NPq),AZ4(NPq),
     *          BX1(NPq),BX2(NPq),BX3(NPq),BX4(NPq),
     *          BY1(NPq),BY2(NPq),BY3(NPq),BY4(NPq),
     *          CX1(NPq),CX2(NPq),CX3(NPq),CX4(NPq),
     *          CY1(NPq),CY2(NPq),CY3(NPq),CY4(NPq),
     *          CZ1(NPq),CZ2(NPq),CZ3(NPq),CZ4(NPq)
      COMMON /BR/ AAX1(NPq),AAX2(NPq),AAX3(NPq),AAX4(NPq),
     *          ABX1(NPq),ABX2(NPq),ABX3(NPq),ABX4(NPq),
     *          ACX1(NPq),ACX2(NPq),ACX3(NPq),ACX4(NPq),
     *          BBX1(NPq),BBX2(NPq),BBX3(NPq),BBX4(NPq),
     *          BCX1(NPq),BCX2(NPq),BCX3(NPq),BCX4(NPq),
     *          CCX1(NPq),CCX2(NPq),CCX3(NPq),CCX4(NPq)
      COMMON /BR/ AAY1(NPq),AAY2(NPq),AAY3(NPq),AAY4(NPq),
     *          ABY1(NPq),ABY2(NPq),ABY3(NPq),ABY4(NPq),
     *          ACY1(NPq),ACY2(NPq),ACY3(NPq),ACY4(NPq),
     *          BBY1(NPq),BBY2(NPq),BBY3(NPq),BBY4(NPq),
     *          BCY1(NPq),BCY2(NPq),BCY3(NPq),BCY4(NPq),
     *          CCY1(NPq),CCY2(NPq),CCY3(NPq),CCY4(NPq)
      COMMON /BR/ AAZ1(NPq),AAZ2(NPq),AAZ3(NPq),AAZ4(NPq),
C    *          ABZ1(NPq),ABZ2(NPq),ABZ3(NPq),ABZ4(NPq),
     *          ACZ1(NPq),ACZ2(NPq),ACZ3(NPq),ACZ4(NPq),
C    *          BBZ1(NPq),BBZ2(NPq),BBZ3(NPq),BBZ4(NPq),
C    *          BCZ1(NPq),BCZ2(NPq),BCZ3(NPq),BCZ4(NPq),
     *          CCZ1(NPq),CCZ2(NPq),CCZ3(NPq),CCZ4(NPq)
      COMMON /CR/ FDI(NDq),SD(NDq,NDq)
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6
      COMMON /DR/ BYL,BZL
      COMMON /IN/ LC(NPq),N,NN,N3,N6
      print *,"DRVTV",N
      EP=0.0D0
      DO 60 I=1,N-1
      XX=X(I)
      YY=Y(I)
      ZZ=Z(I)
      XL1=XX1(I)
      YL1=YY1(I)
      ZL1=ZZ1(I)
      XL2=XX2(I)
      YL2=YY2(I)
      ZL2=ZZ2(I)
      XL3=XX3(I)
      YL3=YY3(I)
      ZL3=ZZ3(I)
      XL4=XX4(I)
      YL4=YY4(I)
      ZL4=ZZ4(I)
      NTI1=3*I-2
      NRI1=N3+NTI1
      NRI2=NRI1+1
      NRI3=NRI2+1
      NTI2=NTI1+1
      NTI3=NTI2+1
      DO 70 J=I+1,N
      DXT=XX-X(J)
      DYT=YY-Y(J)
      DZT=ZZ-Z(J)
      do ix=0,0
         do iy=0,0
            do iz=0,0
!      do ix=-1,1
!      do iy=-1,1
!      do iz=-1,1
      dx=dxt+real(ix)*bxl
      dy=dyt+real(iy)*byl
      dz=dzt+real(iz)*bzl
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      R2=SX+SY+SZ
      print *,"R",i,j,R2**0.5
      IF(R2.LE.RC2) THEN 
      XK1=XX1(J)
      YK1=YY1(J)
      ZK1=ZZ1(J)
      XK2=XX2(J)
      YK2=YY2(J)
      ZK2=ZZ2(J)
      XK3=XX3(J)
      YK3=YY3(J)
      ZK3=ZZ3(J)
      XK4=XX4(J)
      YK4=YY4(J)
      ZK4=ZZ4(J)
      DX1=DX+XL1
      DY1=DY+YL1
      DZ1=DZ+ZL1
      DX2=DX+XL2
      DY2=DY+YL2
      DZ2=DZ+ZL2
      DX3=DX+XL3
      DY3=DY+YL3
      DZ3=DZ+ZL3
      DX4=DX+XL4
      DY4=DY+YL4
      DZ4=DZ+ZL4
      RR=SQRT(R2)
      IF(RR.GT.RL) THEN
      RSI=1.0D0/R2
      RI=1.0D0/RR
      RRL=RR-RC
      RRS=RR-RL
      RRL2=RRL*RRL
      RRS2=RRS*RRS
      SF=RRL2*RRL*RQ*(10.0D0*RRS2-5.0D0*RRL*RRS+RRL2)
      DF=RQ3*RRL2*RRS2
      DS=RQ6*RRL*RRS*(RRL+RRS)
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      DFRI=DF*RI
      SDX=DFRI*DX
      SDY=DFRI*DY
      SDZ=DFRI*DZ
      DSF=DS-DFRI
      DSFRSI=DSF*RSI
      SXX=DFRI+DSFRSI*SX
      SYY=DFRI+DSFRSI*SY
      SZZ=DFRI+DSFRSI*SZ
      SXY=DSFRSI*XY
      SYZ=DSFRSI*YZ
      SZX=DSFRSI*ZX
      END IF
C ::::::: FIRST AND SECOND DERIVATIVES FOR EACH SITE
      DX=DX1-XK1
      DY=DY1-YK1
      DZ=DZ1-ZK1
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=-RSI
      RTI=RI*RSI
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX1(I)*DX+BY1(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX1(I))*RTI
      DXAJ=(RI3X*RAJ-AX1(J))*RTI
      DXBI=(RI3X*RBI-BX1(I))*RTI
      DXBJ=(RI3X*RBJ-BX1(J))*RTI
      DXCI=(RI3X*RCI-CX1(I))*RTI
      DXCJ=(RI3X*RCJ-CX1(J))*RTI
      DYAI=(RI3Y*RAI-AY1(I))*RTI
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI
      DYBI=(RI3Y*RBI-BY1(I))*RTI
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI
      DYCI=(RI3Y*RCI-CY1(I))*RTI
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI
      DZAI=(RI3Z*RAI-AZ1(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ1(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I)
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J)
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I)
     *         +DX*ABX1(I)+DY*ABY1(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J)
     *         -DX*ABX1(J)-DY*ABY1(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I)
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J)
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I)
     *         +DX*BBX1(I)+DY*BBY1(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J)
     *         -DX*BBX1(J)-DY*BBY1(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I)
     *         +DX*BCX1(I)+DY*BCY1(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J)
     *         -DX*BCX1(J)-DY*BCY1(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I)
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J)
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX1(I)*AX1(J)+AY1(I)*AY1(J)+AZ1(I)*AZ1(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX1(I)*BX1(J)+AY1(I)*BY1(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX1(I)*AX1(J)+BY1(I)*AY1(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX1(I)*CX1(J)+AY1(I)*CY1(J)+AZ1(I)*CZ1(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX1(I)*AX1(J)+CY1(I)*AY1(J)+CZ1(I)*AZ1(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX1(I)*BX1(J)+BY1(I)*BY1(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX1(I)*CX1(J)+BY1(I)*CY1(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX1(I)*BX1(J)+CY1(I)*BY1(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX1(I)*CX1(J)+CY1(I)*CY1(J)+CZ1(I)*CZ1(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX
      FDY=DDY
      FDZ=DDZ
      FDAI=DDAI
      FDBI=DDBI
      FDCI=DDCI
      FDAJ=DDAJ
      FDBJ=DDBJ
      FDCJ=DDCJ
      SDXX=DXX
      SDYY=DYY
      SDZZ=DZZ
      SDXY=DXY
      SDYZ=DYZ
      SDZX=DZX
      SDXAI=DXAI
      SDYAI=DYAI
      SDZAI=DZAI
      SDXBI=DXBI
      SDYBI=DYBI
      SDZBI=DZBI
      SDXCI=DXCI
      SDYCI=DYCI
      SDZCI=DZCI
      SDXAJ=DXAJ
      SDYAJ=DYAJ
      SDZAJ=DZAJ
      SDXBJ=DXBJ
      SDYBJ=DYBJ
      SDZBJ=DZBJ
      SDXCJ=DXCJ
      SDYCJ=DYCJ
      SDZCJ=DZCJ
      SDAAI=DAAI
      SDABI=DABI
      SDACI=DACI
      SDBBI=DBBI
      SDBCI=DBCI
      SDCCI=DCCI
      SDAAJ=DAAJ
      SDABJ=DABJ
      SDACJ=DACJ
      SDBBJ=DBBJ
      SDBCJ=DBCJ
      SDCCJ=DCCJ
      SDAAK=DAAK
      SDABK=DABK
      SDBAK=DBAK
      SDACK=DACK
      SDCAK=DCAK
      SDBBK=DBBK
      SDBCK=DBCK
      SDCBK=DCBK
      SDCCK=DCCK
      EPI=RI
      print *,"ri11",ri
C ::::: @@@ 11 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX1-XK2
      DY=DY1-YK2
      DZ=DZ1-ZK2
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=-RSI
      RTI=RI*RSI
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX1(I)*DX+BY1(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX1(I))*RTI
      DXAJ=(RI3X*RAJ-AX2(J))*RTI
      DXBI=(RI3X*RBI-BX1(I))*RTI
      DXBJ=(RI3X*RBJ-BX2(J))*RTI
      DXCI=(RI3X*RCI-CX1(I))*RTI
      DXCJ=(RI3X*RCJ-CX2(J))*RTI
      DYAI=(RI3Y*RAI-AY1(I))*RTI
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI
      DYBI=(RI3Y*RBI-BY1(I))*RTI
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI
      DYCI=(RI3Y*RCI-CY1(I))*RTI
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI
      DZAI=(RI3Z*RAI-AZ1(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ1(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I)
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J)
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I)
     *         +DX*ABX1(I)+DY*ABY1(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J)
     *         -DX*ABX2(J)-DY*ABY2(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I)
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J)
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I)
     *         +DX*BBX1(I)+DY*BBY1(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J)
     *         -DX*BBX2(J)-DY*BBY2(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I)
     *         +DX*BCX1(I)+DY*BCY1(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J)
     *         -DX*BCX2(J)-DY*BCY2(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I)
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J)
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX1(I)*AX2(J)+AY1(I)*AY2(J)+AZ1(I)*AZ2(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX1(I)*BX2(J)+AY1(I)*BY2(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX1(I)*AX2(J)+BY1(I)*AY2(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX1(I)*CX2(J)+AY1(I)*CY2(J)+AZ1(I)*CZ2(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX1(I)*AX2(J)+CY1(I)*AY2(J)+CZ1(I)*AZ2(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX1(I)*BX2(J)+BY1(I)*BY2(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX1(I)*CX2(J)+BY1(I)*CY2(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX1(I)*BX2(J)+CY1(I)*BY2(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX1(I)*CX2(J)+CY1(I)*CY2(J)+CZ1(I)*CZ2(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI+RI
      print *,"ri12",ri
C ::::: @@@ 12 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX2-XK1
      DY=DY2-YK1
      DZ=DZ2-ZK1
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=-RSI
      RTI=RI*RSI
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX2(I)*DX+BY2(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX2(I))*RTI
      DXAJ=(RI3X*RAJ-AX1(J))*RTI
      DXBI=(RI3X*RBI-BX2(I))*RTI
      DXBJ=(RI3X*RBJ-BX1(J))*RTI
      DXCI=(RI3X*RCI-CX2(I))*RTI
      DXCJ=(RI3X*RCJ-CX1(J))*RTI
      DYAI=(RI3Y*RAI-AY2(I))*RTI
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI
      DYBI=(RI3Y*RBI-BY2(I))*RTI
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI
      DYCI=(RI3Y*RCI-CY2(I))*RTI
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI
      DZAI=(RI3Z*RAI-AZ2(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ2(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I)
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J)
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I)
     *         +DX*ABX2(I)+DY*ABY2(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J)
     *         -DX*ABX1(J)-DY*ABY1(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I)
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J)
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I)
     *         +DX*BBX2(I)+DY*BBY2(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J)
     *         -DX*BBX1(J)-DY*BBY1(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I)
     *         +DX*BCX2(I)+DY*BCY2(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J)
     *         -DX*BCX1(J)-DY*BCY1(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I)
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J)
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX2(I)*AX1(J)+AY2(I)*AY1(J)+AZ2(I)*AZ1(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX2(I)*BX1(J)+AY2(I)*BY1(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX2(I)*AX1(J)+BY2(I)*AY1(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX2(I)*CX1(J)+AY2(I)*CY1(J)+AZ2(I)*CZ1(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX2(I)*AX1(J)+CY2(I)*AY1(J)+CZ2(I)*AZ1(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX2(I)*BX1(J)+BY2(I)*BY1(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX2(I)*CX1(J)+BY2(I)*CY1(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX2(I)*BX1(J)+CY2(I)*BY1(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX2(I)*CX1(J)+CY2(I)*CY1(J)+CZ2(I)*CZ1(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI+RI
      print *,"ri21",ri
C ::::: @@@ 21 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX2-XK2
      DY=DY2-YK2
      DZ=DZ2-ZK2
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=-RSI
      RTI=RI*RSI
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX2(I)*DX+BY2(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX2(I))*RTI
      DXAJ=(RI3X*RAJ-AX2(J))*RTI
      DXBI=(RI3X*RBI-BX2(I))*RTI
      DXBJ=(RI3X*RBJ-BX2(J))*RTI
      DXCI=(RI3X*RCI-CX2(I))*RTI
      DXCJ=(RI3X*RCJ-CX2(J))*RTI
      DYAI=(RI3Y*RAI-AY2(I))*RTI
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI
      DYBI=(RI3Y*RBI-BY2(I))*RTI
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI
      DYCI=(RI3Y*RCI-CY2(I))*RTI
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI
      DZAI=(RI3Z*RAI-AZ2(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ2(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I)
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J)
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I)
     *         +DX*ABX2(I)+DY*ABY2(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J)
     *         -DX*ABX2(J)-DY*ABY2(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I)
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J)
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I)
     *         +DX*BBX2(I)+DY*BBY2(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J)
     *         -DX*BBX2(J)-DY*BBY2(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I)
     *         +DX*BCX2(I)+DY*BCY2(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J)
     *         -DX*BCX2(J)-DY*BCY2(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I)
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J)
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX2(I)*AX2(J)+AY2(I)*AY2(J)+AZ2(I)*AZ2(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX2(I)*BX2(J)+AY2(I)*BY2(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX2(I)*AX2(J)+BY2(I)*AY2(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX2(I)*CX2(J)+AY2(I)*CY2(J)+AZ2(I)*CZ2(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX2(I)*AX2(J)+CY2(I)*AY2(J)+CZ2(I)*AZ2(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX2(I)*BX2(J)+BY2(I)*BY2(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX2(I)*CX2(J)+BY2(I)*CY2(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX2(I)*BX2(J)+CY2(I)*BY2(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX2(I)*CX2(J)+CY2(I)*CY2(J)+CZ2(I)*CZ2(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI+RI
      print *,"ri22",ri
C ::::: @@@ 22 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX1-XK4
      DY=DY1-YK4
      DZ=DZ1-ZK4
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=2.0D0*RSI
      RTI=-RI*RSIM
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX1(I)*DX+BY1(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX1(I))*RTI
      DXAJ=(RI3X*RAJ-AX4(J))*RTI
      DXBI=(RI3X*RBI-BX1(I))*RTI
      DXBJ=(RI3X*RBJ-BX4(J))*RTI
      DXCI=(RI3X*RCI-CX1(I))*RTI
      DXCJ=(RI3X*RCJ-CX4(J))*RTI
      DYAI=(RI3Y*RAI-AY1(I))*RTI
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI
      DYBI=(RI3Y*RBI-BY1(I))*RTI
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI
      DYCI=(RI3Y*RCI-CY1(I))*RTI
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI
      DZAI=(RI3Z*RAI-AZ1(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ1(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I)
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J)
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I)
     *         +DX*ABX1(I)+DY*ABY1(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J)
     *         -DX*ABX4(J)-DY*ABY4(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I)
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J)
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I)
     *         +DX*BBX1(I)+DY*BBY1(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J)
     *         -DX*BBX4(J)-DY*BBY4(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I)
     *         +DX*BCX1(I)+DY*BCY1(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J)
     *         -DX*BCX4(J)-DY*BCY4(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I)
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J)
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX1(I)*AX4(J)+AY1(I)*AY4(J)+AZ1(I)*AZ4(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX1(I)*BX4(J)+AY1(I)*BY4(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX1(I)*AX4(J)+BY1(I)*AY4(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX1(I)*CX4(J)+AY1(I)*CY4(J)+AZ1(I)*CZ4(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX1(I)*AX4(J)+CY1(I)*AY4(J)+CZ1(I)*AZ4(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX1(I)*BX4(J)+BY1(I)*BY4(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX1(I)*CX4(J)+BY1(I)*CY4(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX1(I)*BX4(J)+CY1(I)*BY4(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX1(I)*CX4(J)+CY1(I)*CY4(J)+CZ1(I)*CZ4(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI-2.0D0*RI
      print *,"ri14",ri
C ::::: @@@ 14 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX2-XK4
      DY=DY2-YK4
      DZ=DZ2-ZK4
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=2.0D0*RSI
      RTI=-RI*RSIM
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX2(I)*DX+BY2(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX2(I))*RTI
      DXAJ=(RI3X*RAJ-AX4(J))*RTI
      DXBI=(RI3X*RBI-BX2(I))*RTI
      DXBJ=(RI3X*RBJ-BX4(J))*RTI
      DXCI=(RI3X*RCI-CX2(I))*RTI
      DXCJ=(RI3X*RCJ-CX4(J))*RTI
      DYAI=(RI3Y*RAI-AY2(I))*RTI
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI
      DYBI=(RI3Y*RBI-BY2(I))*RTI
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI
      DYCI=(RI3Y*RCI-CY2(I))*RTI
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI
      DZAI=(RI3Z*RAI-AZ2(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ2(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I)
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J)
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I)
     *         +DX*ABX2(I)+DY*ABY2(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J)
     *         -DX*ABX4(J)-DY*ABY4(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I)
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J)
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I)
     *         +DX*BBX2(I)+DY*BBY2(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J)
     *         -DX*BBX4(J)-DY*BBY4(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I)
     *         +DX*BCX2(I)+DY*BCY2(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J)
     *         -DX*BCX4(J)-DY*BCY4(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I)
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J)
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX2(I)*AX4(J)+AY2(I)*AY4(J)+AZ2(I)*AZ4(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX2(I)*BX4(J)+AY2(I)*BY4(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX2(I)*AX4(J)+BY2(I)*AY4(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX2(I)*CX4(J)+AY2(I)*CY4(J)+AZ2(I)*CZ4(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX2(I)*AX4(J)+CY2(I)*AY4(J)+CZ2(I)*AZ4(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX2(I)*BX4(J)+BY2(I)*BY4(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX2(I)*CX4(J)+BY2(I)*CY4(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX2(I)*BX4(J)+CY2(I)*BY4(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX2(I)*CX4(J)+CY2(I)*CY4(J)+CZ2(I)*CZ4(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI-2.0D0*RI
      print *,"ri24",ri
C ::::: @@@ 24 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX4-XK1
      DY=DY4-YK1
      DZ=DZ4-ZK1
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=2.0D0*RSI
      RTI=-RI*RSIM
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX4(I)*DX+BY4(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX4(I))*RTI
      DXAJ=(RI3X*RAJ-AX1(J))*RTI
      DXBI=(RI3X*RBI-BX4(I))*RTI
      DXBJ=(RI3X*RBJ-BX1(J))*RTI
      DXCI=(RI3X*RCI-CX4(I))*RTI
      DXCJ=(RI3X*RCJ-CX1(J))*RTI
      DYAI=(RI3Y*RAI-AY4(I))*RTI
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI
      DYBI=(RI3Y*RBI-BY4(I))*RTI
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI
      DYCI=(RI3Y*RCI-CY4(I))*RTI
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI
      DZAI=(RI3Z*RAI-AZ4(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ4(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I)
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J)
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I)
     *         +DX*ABX4(I)+DY*ABY4(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J)
     *         -DX*ABX1(J)-DY*ABY1(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I)
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J)
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I)
     *         +DX*BBX4(I)+DY*BBY4(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J)
     *         -DX*BBX1(J)-DY*BBY1(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I)
     *         +DX*BCX4(I)+DY*BCY4(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J)
     *         -DX*BCX1(J)-DY*BCY1(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I)
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J)
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX4(I)*AX1(J)+AY4(I)*AY1(J)+AZ4(I)*AZ1(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX4(I)*BX1(J)+AY4(I)*BY1(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX4(I)*AX1(J)+BY4(I)*AY1(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX4(I)*CX1(J)+AY4(I)*CY1(J)+AZ4(I)*CZ1(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX4(I)*AX1(J)+CY4(I)*AY1(J)+CZ4(I)*AZ1(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX4(I)*BX1(J)+BY4(I)*BY1(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX4(I)*CX1(J)+BY4(I)*CY1(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX4(I)*BX1(J)+CY4(I)*BY1(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX4(I)*CX1(J)+CY4(I)*CY1(J)+CZ4(I)*CZ1(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI-2.0D0*RI
      print *,"ri41",ri
C ::::: @@@ 41 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX4-XK2
      DY=DY4-YK2
      DZ=DZ4-ZK2
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=2.0D0*RSI
      RTI=-RI*RSIM
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX4(I)*DX+BY4(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX4(I))*RTI
      DXAJ=(RI3X*RAJ-AX2(J))*RTI
      DXBI=(RI3X*RBI-BX4(I))*RTI
      DXBJ=(RI3X*RBJ-BX2(J))*RTI
      DXCI=(RI3X*RCI-CX4(I))*RTI
      DXCJ=(RI3X*RCJ-CX2(J))*RTI
      DYAI=(RI3Y*RAI-AY4(I))*RTI
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI
      DYBI=(RI3Y*RBI-BY4(I))*RTI
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI
      DYCI=(RI3Y*RCI-CY4(I))*RTI
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI
      DZAI=(RI3Z*RAI-AZ4(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ4(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I)
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J)
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I)
     *         +DX*ABX4(I)+DY*ABY4(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J)
     *         -DX*ABX2(J)-DY*ABY2(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I)
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J)
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I)
     *         +DX*BBX4(I)+DY*BBY4(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J)
     *         -DX*BBX2(J)-DY*BBY2(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I)
     *         +DX*BCX4(I)+DY*BCY4(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J)
     *         -DX*BCX2(J)-DY*BCY2(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I)
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J)
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX4(I)*AX2(J)+AY4(I)*AY2(J)+AZ4(I)*AZ2(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX4(I)*BX2(J)+AY4(I)*BY2(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX4(I)*AX2(J)+BY4(I)*AY2(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX4(I)*CX2(J)+AY4(I)*CY2(J)+AZ4(I)*CZ2(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX4(I)*AX2(J)+CY4(I)*AY2(J)+CZ4(I)*AZ2(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX4(I)*BX2(J)+BY4(I)*BY2(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX4(I)*CX2(J)+BY4(I)*CY2(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX4(I)*BX2(J)+CY4(I)*BY2(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX4(I)*CX2(J)+CY4(I)*CY2(J)+CZ4(I)*CZ2(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI-2.0D0*RI
      print *,"ri42",ri
C ::::: @@@ 42 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX4-XK4
      DY=DY4-YK4
      DZ=DZ4-ZK4
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RSIM=-4.0D0*RSI
      RTI=-RI*RSIM
      RFI=RSI*RTI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      RTIM=-RTI
      DDX=RTIM*DX
      DDY=RTIM*DY
      DDZ=RTIM*DZ
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ)
      DDAI=RSIM*RAI
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ)
      DDAJ=RSIM*RAJ
      RBI=RI*(BX4(I)*DX+BY4(I)*DY)
      DDBI=RSIM*RBI
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY)
      DDBJ=RSIM*RBJ
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ)
      DDCI=RSIM*RCI
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ)
      DDCJ=RSIM*RCJ
      RFI3=RFI*3.0D0
      DXX=RFI3*SX+RTIM
      DYY=RFI3*SY+RTIM
      DZZ=RFI3*SZ+RTIM
      DXY=RFI3*XY
      DYZ=RFI3*YZ
      DZX=RFI3*ZX
      RI3=3.0D0*RI
      RI3X=DX*RI3
      RI3Y=DY*RI3
      RI3Z=DZ*RI3
      DXAI=(RI3X*RAI-AX4(I))*RTI
      DXAJ=(RI3X*RAJ-AX4(J))*RTI
      DXBI=(RI3X*RBI-BX4(I))*RTI
      DXBJ=(RI3X*RBJ-BX4(J))*RTI
      DXCI=(RI3X*RCI-CX4(I))*RTI
      DXCJ=(RI3X*RCJ-CX4(J))*RTI
      DYAI=(RI3Y*RAI-AY4(I))*RTI
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI
      DYBI=(RI3Y*RBI-BY4(I))*RTI
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI
      DYCI=(RI3Y*RCI-CY4(I))*RTI
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI
      DZAI=(RI3Z*RAI-AZ4(I))*RTI
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI
      DZBI=RI3Z*RBI*RTI
      DZBJ=RI3Z*RBJ*RTI
      DZCI=(RI3Z*RCI-CZ4(I))*RTI
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI
      RTI3=3.0D0*RTI
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I)
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I))
     *    +RTI3*RAI*RAI
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J)
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J))
     *    +RTI3*RAJ*RAJ
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I)
     *         +DX*ABX4(I)+DY*ABY4(I))
     *    +RTI3*RAI*RBI
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J)
     *         -DX*ABX4(J)-DY*ABY4(J))
     *    +RTI3*RAJ*RBJ
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I)
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I))
     *    +RTI3*RAI*RCI
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J)
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J))
     *    +RTI3*RAJ*RCJ
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I)
     *         +DX*BBX4(I)+DY*BBY4(I))
     *    +RTI3*RBI*RBI
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J)
     *         -DX*BBX4(J)-DY*BBY4(J))
     *    +RTI3*RBJ*RBJ
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I)
     *         +DX*BCX4(I)+DY*BCY4(I))
     *    +RTI3*RBI*RCI
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J)
     *         -DX*BCX4(J)-DY*BCY4(J))
     *    +RTI3*RBJ*RCJ
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I)
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I))
     *    +RTI3*RCI*RCI
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J)
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J))
     *    +RTI3*RCJ*RCJ
      DAAK=RTIM*(AX4(I)*AX4(J)+AY4(I)*AY4(J)+AZ4(I)*AZ4(J))
     *    +RTI3*RAI*RAJ
      DABK=RTIM*(AX4(I)*BX4(J)+AY4(I)*BY4(J))
     *    +RTI3*RAI*RBJ
      DBAK=RTIM*(BX4(I)*AX4(J)+BY4(I)*AY4(J))
     *    +RTI3*RAJ*RBI
      DACK=RTIM*(AX4(I)*CX4(J)+AY4(I)*CY4(J)+AZ4(I)*CZ4(J))
     *    +RTI3*RAI*RCJ
      DCAK=RTIM*(CX4(I)*AX4(J)+CY4(I)*AY4(J)+CZ4(I)*AZ4(J))
     *    +RTI3*RAJ*RCI
      DBBK=RTIM*(BX4(I)*BX4(J)+BY4(I)*BY4(J))
     *    +RTI3*RBI*RBJ
      DBCK=RTIM*(BX4(I)*CX4(J)+BY4(I)*CY4(J))
     *    +RTI3*RBI*RCJ
      DCBK=RTIM*(CX4(I)*BX4(J)+CY4(I)*BY4(J))
     *    +RTI3*RBJ*RCI
      DCCK=RTIM*(CX4(I)*CX4(J)+CY4(I)*CY4(J)+CZ4(I)*CZ4(J))
     *    +RTI3*RCI*RCJ
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI+4.0D0*RI
      print *,"ri44",ri
C ::::: @@@ 44 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      DX=DX3-XK3
      DY=DY3-YK3
      DZ=DZ3-ZK3
      SX=DX*DX
      SY=DY*DY
      SZ=DZ*DZ
      RSI=1.0D0/(SX+SY+SZ)
      RI=SQRT(RSI)
      RHI=RSI*RSI*RSI
      RDI=RHI*RHI
      DE=P1*RDI+P2*RHI
      DF=(P3*RDI+P4*RHI)*RI
      DS=(P5*RDI+P6*RHI)*RSI
      XY=DX*DY
      YZ=DY*DZ
      ZX=DZ*DX
      DFRI=DF*RI
      DDX=DFRI*DX
      DDY=DFRI*DY
      DDZ=DFRI*DZ
      RAI=RI*(AX3(I)*DX+AY3(I)*DY+AZ3(I)*DZ)
      DDAI=DF*RAI
      RAJ=RI*(AX3(J)*DX+AY3(J)*DY+AZ3(J)*DZ)
      DDAJ=DF*RAJ
      RBI=RI*(BX3(I)*DX+BY3(I)*DY)
      DDBI=DF*RBI
      RBJ=RI*(BX3(J)*DX+BY3(J)*DY)
      DDBJ=DF*RBJ
      RCI=RI*(CX3(I)*DX+CY3(I)*DY+CZ3(I)*DZ)
      DDCI=DF*RCI
      RCJ=RI*(CX3(J)*DX+CY3(J)*DY+CZ3(J)*DZ)
      DDCJ=DF*RCJ
      DSF=DS-DFRI
      DSFRSI=DSF*RSI
      DSFRI=DSF*RI
      DXX=DFRI+DSFRSI*SX
      DYY=DFRI+DSFRSI*SY
      DZZ=DFRI+DSFRSI*SZ
      DXY=DSFRSI*XY
      DYZ=DSFRSI*YZ
      DZX=DSFRSI*ZX
      DXDSF=DX*DSFRI
      DYDSF=DY*DSFRI
      DZDSF=DZ*DSFRI
      DXAI=DFRI*AX3(I)+RAI*DXDSF
      DXAJ=DFRI*AX3(J)+RAJ*DXDSF
      DXBI=DFRI*BX3(I)+RBI*DXDSF
      DXBJ=DFRI*BX3(J)+RBJ*DXDSF
      DXCI=DFRI*CX3(I)+RCI*DXDSF
      DXCJ=DFRI*CX3(J)+RCJ*DXDSF
      DYAI=DFRI*AY3(I)+RAI*DYDSF
      DYAJ=DFRI*AY3(J)+RAJ*DYDSF
      DYBI=DFRI*BY3(I)+RBI*DYDSF
      DYBJ=DFRI*BY3(J)+RBJ*DYDSF
      DYCI=DFRI*CY3(I)+RCI*DYDSF
      DYCJ=DFRI*CY3(J)+RCJ*DYDSF
      DZAI=DFRI*AZ3(I)+RAI*DZDSF
      DZAJ=DFRI*AZ3(J)+RAJ*DZDSF
      DZBI=RBI*DZDSF
      DZBJ=RBJ*DZDSF
      DZCI=DFRI*CZ3(I)+RCI*DZDSF
      DZCJ=DFRI*CZ3(J)+RCJ*DZDSF
      DAAI=DFRI*(AX3(I)*AX3(I)+AY3(I)*AY3(I)+AZ3(I)*AZ3(I)
     *         +DX*AAX3(I)+DY*AAY3(I)+DZ*AAZ3(I))
     *    +DSF*RAI*RAI
      DAAJ=DFRI*(AX3(J)*AX3(J)+AY3(J)*AY3(J)+AZ3(J)*AZ3(J)
     *         -DX*AAX3(J)-DY*AAY3(J)-DZ*AAZ3(J))
     *    +DSF*RAJ*RAJ
      DABI=DFRI*(AX3(I)*BX3(I)+AY3(I)*BY3(I)
     *         +DX*ABX3(I)+DY*ABY3(I))
     *    +DSF*RAI*RBI
      DABJ=DFRI*(AX3(J)*BX3(J)+AY3(J)*BY3(J)
     *         -DX*ABX3(J)-DY*ABY3(J))
     *    +DSF*RAJ*RBJ
      DBBI=DFRI*(BX3(I)*BX3(I)+BY3(I)*BY3(I)
     *         +DX*BBX3(I)+DY*BBY3(I))
     *    +DSF*RBI*RBI
      DBBJ=DFRI*(BX3(J)*BX3(J)+BY3(J)*BY3(J)
     *         -DX*BBX3(J)-DY*BBY3(J))
     *    +DSF*RBJ*RBJ
      DCCI=DFRI*(CX3(I)*CX3(I)+CY3(I)*CY3(I)+CZ3(I)*CZ3(I)
     *         +DX*CCX3(I)+DY*CCY3(I)+DZ*CCZ3(I))
     *    +DSF*RCI*RCI
      DCCJ=DFRI*(CX3(J)*CX3(J)+CY3(J)*CY3(J)+CZ3(J)*CZ3(J)
     *         -DX*CCX3(J)-DY*CCY3(J)-DZ*CCZ3(J))
     *    +DSF*RCJ*RCJ
      DACI=DFRI*(AX3(I)*CX3(I)+AY3(I)*CY3(I)+AZ3(I)*CZ3(I)
     *         +DX*ACX3(I)+DY*ACY3(I)+DZ*ACZ3(I))
     *    +DSF*RAI*RCI
      DACJ=DFRI*(AX3(J)*CX3(J)+AY3(J)*CY3(J)+AZ3(J)*CZ3(J)
     *         -DX*ACX3(J)-DY*ACY3(J)-DZ*ACZ3(J))
     *    +DSF*RAJ*RCJ
      DBCI=DFRI*(BX3(I)*CX3(I)+BY3(I)*CY3(I)
     *         +DX*BCX3(I)+DY*BCY3(I))
     *    +DSF*RBI*RCI
      DBCJ=DFRI*(BX3(J)*CX3(J)+BY3(J)*CY3(J)
     *         -DX*BCX3(J)-DY*BCY3(J))
     *    +DSF*RBJ*RCJ
      DAAK=DFRI*(AX3(I)*AX3(J)+AY3(I)*AY3(J)+AZ3(I)*AZ3(J))
     *    +DSF*RAI*RAJ
      DABK=DFRI*(AX3(I)*BX3(J)+AY3(I)*BY3(J))
     *    +DSF*RAI*RBJ
      DBAK=DFRI*(BX3(I)*AX3(J)+BY3(I)*AY3(J))
     *    +DSF*RAJ*RBI
      DBBK=DFRI*(BX3(I)*BX3(J)+BY3(I)*BY3(J))
     *    +DSF*RBI*RBJ
      DCCK=DFRI*(CX3(I)*CX3(J)+CY3(I)*CY3(J)+CZ3(I)*CZ3(J))
     *    +DSF*RCI*RCJ
      DACK=DFRI*(AX3(I)*CX3(J)+AY3(I)*CY3(J)+AZ3(I)*CZ3(J))
     *    +DSF*RAI*RCJ
      DCAK=DFRI*(CX3(I)*AX3(J)+CY3(I)*AY3(J)+CZ3(I)*AZ3(J))
     *    +DSF*RAJ*RCI
      DBCK=DFRI*(BX3(I)*CX3(J)+BY3(I)*CY3(J))
     *    +DSF*RBI*RCJ
      DCBK=DFRI*(CX3(I)*BX3(J)+CY3(I)*BY3(J))
     *    +DSF*RBJ*RCI
      FDX=DDX+FDX
      FDY=DDY+FDY
      FDZ=DDZ+FDZ
      FDAI=FDAI+DDAI
      FDBI=FDBI+DDBI
      FDCI=FDCI+DDCI
      FDAJ=FDAJ+DDAJ
      FDBJ=FDBJ+DDBJ
      FDCJ=FDCJ+DDCJ
      SDXX=DXX+SDXX
      SDYY=DYY+SDYY
      SDZZ=DZZ+SDZZ
      SDXY=DXY+SDXY
      SDYZ=DYZ+SDYZ
      SDZX=DZX+SDZX
      SDXAI=DXAI+SDXAI
      SDYAI=DYAI+SDYAI
      SDZAI=DZAI+SDZAI
      SDXBI=DXBI+SDXBI
      SDYBI=DYBI+SDYBI
      SDZBI=DZBI+SDZBI
      SDXCI=DXCI+SDXCI
      SDYCI=DYCI+SDYCI
      SDZCI=DZCI+SDZCI
      SDXAJ=DXAJ+SDXAJ
      SDYAJ=DYAJ+SDYAJ
      SDZAJ=DZAJ+SDZAJ
      SDXBJ=DXBJ+SDXBJ
      SDYBJ=DYBJ+SDYBJ
      SDZBJ=DZBJ+SDZBJ
      SDXCJ=DXCJ+SDXCJ
      SDYCJ=DYCJ+SDYCJ
      SDZCJ=DZCJ+SDZCJ
      SDAAI=DAAI+SDAAI
      SDABI=DABI+SDABI
      SDACI=DACI+SDACI
      SDBBI=DBBI+SDBBI
      SDBCI=DBCI+SDBCI
      SDCCI=DCCI+SDCCI
      SDAAJ=DAAJ+SDAAJ
      SDABJ=DABJ+SDABJ
      SDACJ=DACJ+SDACJ
      SDBBJ=DBBJ+SDBBJ
      SDBCJ=DBCJ+SDBCJ
      SDCCJ=DCCJ+SDCCJ
      SDAAK=DAAK+SDAAK
      SDABK=DABK+SDABK
      SDBAK=DBAK+SDBAK
      SDACK=DACK+SDACK
      SDCAK=DCAK+SDCAK
      SDBBK=DBBK+SDBBK
      SDBCK=DBCK+SDBCK
      SDCBK=DCBK+SDCBK
      SDCCK=DCCK+SDCCK
      EPI=EPI+DE
      print *,"de33",de
C ::::: @@@ 33 ::::::::::::::::::::::::::::::::::::::::::::::::::::
      NTJ1=3*J-2
      NRJ1=N3+NTJ1
      NRJ2=NRJ1+1
      NRJ3=NRJ2+1
      NTJ2=NTJ1+1
      NTJ3=NTJ2+1
      IF(RR.GT.RL) THEN
      U=EPI
      VDX=FDX
      VDY=FDY
      VDZ=FDZ
      EPI=EPI*SF
      FDX=VDX*SF+SDX*U
      FDY=VDY*SF+SDY*U
      FDZ=VDZ*SF+SDZ*U
      SDXX=SDXX*SF+2.0D0*VDX*SDX+SXX*U
      SDYY=SDYY*SF+2.0D0*VDY*SDY+SYY*U
      SDZZ=SDZZ*SF+2.0D0*VDZ*SDZ+SZZ*U
      SDXY=SDXY*SF+VDX*SDY+VDY*SDX+SXY*U
      SDYZ=SDYZ*SF+VDY*SDZ+VDZ*SDY+SYZ*U
      SDZX=SDZX*SF+VDZ*SDX+VDX*SDZ+SZX*U
      SDXAI=SDX*FDAI+SDXAI*SF
      SDXBI=SDX*FDBI+SDXBI*SF
      SDXCI=SDX*FDCI+SDXCI*SF
      SDYAI=SDY*FDAI+SDYAI*SF
      SDYBI=SDY*FDBI+SDYBI*SF
      SDYCI=SDY*FDCI+SDYCI*SF
      SDZAI=SDZ*FDAI+SDZAI*SF
      SDZBI=SDZ*FDBI+SDZBI*SF
      SDZCI=SDZ*FDCI+SDZCI*SF
      SDXAJ=SDX*FDAJ+SDXAJ*SF
      SDXBJ=SDX*FDBJ+SDXBJ*SF
      SDXCJ=SDX*FDCJ+SDXCJ*SF
      SDYAJ=SDY*FDAJ+SDYAJ*SF
      SDYBJ=SDY*FDBJ+SDYBJ*SF
      SDYCJ=SDY*FDCJ+SDYCJ*SF
      SDZAJ=SDZ*FDAJ+SDZAJ*SF
      SDZBJ=SDZ*FDBJ+SDZBJ*SF
      SDZCJ=SDZ*FDCJ+SDZCJ*SF
      FDAI=FDAI*SF
      FDBI=FDBI*SF
      FDCI=FDCI*SF
      FDAJ=FDAJ*SF
      FDBJ=FDBJ*SF
      FDCJ=FDCJ*SF
      SDAAI=SDAAI*SF
      SDABI=SDABI*SF
      SDACI=SDACI*SF
      SDBBI=SDBBI*SF
      SDBCI=SDBCI*SF
      SDCCI=SDCCI*SF
      SDAAJ=SDAAJ*SF
      SDABJ=SDABJ*SF
      SDACJ=SDACJ*SF
      SDBBJ=SDBBJ*SF
      SDBCJ=SDBCJ*SF
      SDCCJ=SDCCJ*SF
      SDAAK=SDAAK*SF
      SDABK=SDABK*SF
      SDBAK=SDBAK*SF
      SDACK=SDACK*SF
      SDCAK=SDCAK*SF
      SDBBK=SDBBK*SF
      SDBCK=SDBCK*SF
      SDCBK=SDCBK*SF
      SDCCK=SDCCK*SF
      END IF
      if(i.eq.j) epi=epi/2.0d0
      EPID=EPI*EMUPID
      EP=EP+EPID
      SD(NTI1,NTI1)=SDXX+SD(NTI1,NTI1)
      SD(NTI1,NTI2)=SDXY+SD(NTI1,NTI2)
      SD(NTI2,NTI1)=SDXY+SD(NTI2,NTI1)
      SD(NTI1,NTI3)=SDZX+SD(NTI1,NTI3)
      SD(NTI3,NTI1)=SDZX+SD(NTI3,NTI1)
      SD(NTI2,NTI2)=SDYY+SD(NTI2,NTI2)
      SD(NTI2,NTI3)=SDYZ+SD(NTI2,NTI3)
      SD(NTI3,NTI2)=SDYZ+SD(NTI3,NTI2)
      SD(NTI3,NTI3)=SDZZ+SD(NTI3,NTI3)
      SD(NTJ1,NTJ1)=SDXX+SD(NTJ1,NTJ1)
      SD(NTJ1,NTJ2)=SDXY+SD(NTJ1,NTJ2)
      SD(NTJ2,NTJ1)=SDXY+SD(NTJ2,NTJ1)
      SD(NTJ1,NTJ3)=SDZX+SD(NTJ1,NTJ3)
      SD(NTJ3,NTJ1)=SDZX+SD(NTJ3,NTJ1)
      SD(NTJ2,NTJ2)=SDYY+SD(NTJ2,NTJ2)
      SD(NTJ2,NTJ3)=SDYZ+SD(NTJ2,NTJ3)
      SD(NTJ3,NTJ2)=SDYZ+SD(NTJ3,NTJ2)
      SD(NTJ3,NTJ3)=SDZZ+SD(NTJ3,NTJ3)
      SD(NTJ1,NTI1)=-SDXX+SD(NTJ1,NTI1)
      SD(NTJ1,NTI2)=-SDXY+SD(NTJ1,NTI2)
      SD(NTJ2,NTI1)=-SDXY+SD(NTJ2,NTI1)
      SD(NTJ1,NTI3)=-SDZX+SD(NTJ1,NTI3)
      SD(NTJ3,NTI1)=-SDZX+SD(NTJ3,NTI1)
      SD(NTJ2,NTI2)=-SDYY+SD(NTJ2,NTI2)
      SD(NTJ2,NTI3)=-SDYZ+SD(NTJ2,NTI3)
      SD(NTJ3,NTI2)=-SDYZ+SD(NTJ3,NTI2)
      SD(NTJ3,NTI3)=-SDZZ+SD(NTJ3,NTI3)
      SD(NTI1,NTJ1)=-SDXX+SD(NTI1,NTJ1)
      SD(NTI1,NTJ2)=-SDXY+SD(NTI1,NTJ2)
      SD(NTI2,NTJ1)=-SDXY+SD(NTI2,NTJ1)
      SD(NTI1,NTJ3)=-SDZX+SD(NTI1,NTJ3)
      SD(NTI3,NTJ1)=-SDZX+SD(NTI3,NTJ1)
      SD(NTI2,NTJ2)=-SDYY+SD(NTI2,NTJ2)
      SD(NTI2,NTJ3)=-SDYZ+SD(NTI2,NTJ3)
      SD(NTI3,NTJ2)=-SDYZ+SD(NTI3,NTJ2)
      SD(NTI3,NTJ3)=-SDZZ+SD(NTI3,NTJ3)
      SD(NRI1,NRI1)=SDAAI+SD(NRI1,NRI1)
      SD(NRI1,NRI2)=SDABI+SD(NRI1,NRI2)
      SD(NRI2,NRI1)=SDABI+SD(NRI2,NRI1)
      SD(NRI1,NRI3)=SDACI+SD(NRI1,NRI3)
      SD(NRI3,NRI1)=SDACI+SD(NRI3,NRI1)
      SD(NRI2,NRI2)=SDBBI+SD(NRI2,NRI2)
      SD(NRI2,NRI3)=SDBCI+SD(NRI2,NRI3)
      SD(NRI3,NRI2)=SDBCI+SD(NRI3,NRI2)
      SD(NRI3,NRI3)=SDCCI+SD(NRI3,NRI3)
      SD(NRJ1,NRJ1)=SDAAJ+SD(NRJ1,NRJ1)
      SD(NRJ1,NRJ2)=SDABJ+SD(NRJ1,NRJ2)
      SD(NRJ2,NRJ1)=SDABJ+SD(NRJ2,NRJ1)
      SD(NRJ1,NRJ3)=SDACJ+SD(NRJ1,NRJ3)
      SD(NRJ3,NRJ1)=SDACJ+SD(NRJ3,NRJ1)
      SD(NRJ2,NRJ2)=SDBBJ+SD(NRJ2,NRJ2)
      SD(NRJ2,NRJ3)=SDBCJ+SD(NRJ2,NRJ3)
      SD(NRJ3,NRJ2)=SDBCJ+SD(NRJ3,NRJ2)
      SD(NRJ3,NRJ3)=SDCCJ+SD(NRJ3,NRJ3)
      SD(NRI1,NRJ1)=-SDAAK+SD(NRI1,NRJ1)
      SD(NRJ1,NRI1)=-SDAAK+SD(NRJ1,NRI1)
      SD(NRI1,NRJ2)=-SDABK+SD(NRI1,NRJ2)
      SD(NRJ2,NRI1)=-SDABK+SD(NRJ2,NRI1)
      SD(NRI2,NRJ1)=-SDBAK+SD(NRI2,NRJ1)
      SD(NRJ1,NRI2)=-SDBAK+SD(NRJ1,NRI2)
      SD(NRI1,NRJ3)=-SDACK+SD(NRI1,NRJ3)
      SD(NRJ3,NRI1)=-SDACK+SD(NRJ3,NRI1)
      SD(NRI3,NRJ1)=-SDCAK+SD(NRI3,NRJ1)
      SD(NRJ1,NRI3)=-SDCAK+SD(NRJ1,NRI3)
      SD(NRI2,NRJ2)=-SDBBK+SD(NRI2,NRJ2)
      SD(NRJ2,NRI2)=-SDBBK+SD(NRJ2,NRI2)
      SD(NRI2,NRJ3)=-SDBCK+SD(NRI2,NRJ3)
      SD(NRJ3,NRI2)=-SDBCK+SD(NRJ3,NRI2)
      SD(NRI3,NRJ2)=-SDCBK+SD(NRI3,NRJ2)
      SD(NRJ2,NRI3)=-SDCBK+SD(NRJ2,NRI3)
      SD(NRI3,NRJ3)=-SDCCK+SD(NRI3,NRJ3)
      SD(NRJ3,NRI3)=-SDCCK+SD(NRJ3,NRI3)
      SD(NTI1,NRI1)=SDXAI+SD(NTI1,NRI1)
      SD(NRI1,NTI1)=SDXAI+SD(NRI1,NTI1)
      SD(NTI1,NRI2)=SDXBI+SD(NTI1,NRI2)
      SD(NRI2,NTI1)=SDXBI+SD(NRI2,NTI1)
      SD(NTI1,NRI3)=SDXCI+SD(NTI1,NRI3)
      SD(NRI3,NTI1)=SDXCI+SD(NRI3,NTI1)
      SD(NTI2,NRI1)=SDYAI+SD(NTI2,NRI1)
      SD(NRI1,NTI2)=SDYAI+SD(NRI1,NTI2)
      SD(NTI2,NRI2)=SDYBI+SD(NTI2,NRI2)
      SD(NRI2,NTI2)=SDYBI+SD(NRI2,NTI2)
      SD(NTI2,NRI3)=SDYCI+SD(NTI2,NRI3)
      SD(NRI3,NTI2)=SDYCI+SD(NRI3,NTI2)
      SD(NTI3,NRI1)=SDZAI+SD(NTI3,NRI1)
      SD(NRI1,NTI3)=SDZAI+SD(NRI1,NTI3)
      SD(NTI3,NRI2)=SDZBI+SD(NTI3,NRI2)
      SD(NRI2,NTI3)=SDZBI+SD(NRI2,NTI3)
      SD(NTI3,NRI3)=SDZCI+SD(NTI3,NRI3)
      SD(NRI3,NTI3)=SDZCI+SD(NRI3,NTI3)
      SD(NTJ1,NRJ1)=SDXAJ+SD(NTJ1,NRJ1)
      SD(NRJ1,NTJ1)=SDXAJ+SD(NRJ1,NTJ1)
      SD(NTJ1,NRJ2)=SDXBJ+SD(NTJ1,NRJ2)
      SD(NRJ2,NTJ1)=SDXBJ+SD(NRJ2,NTJ1)
      SD(NTJ1,NRJ3)=SDXCJ+SD(NTJ1,NRJ3)
      SD(NRJ3,NTJ1)=SDXCJ+SD(NRJ3,NTJ1)
      SD(NTJ2,NRJ1)=SDYAJ+SD(NTJ2,NRJ1)
      SD(NRJ1,NTJ2)=SDYAJ+SD(NRJ1,NTJ2)
      SD(NTJ2,NRJ2)=SDYBJ+SD(NTJ2,NRJ2)
      SD(NRJ2,NTJ2)=SDYBJ+SD(NRJ2,NTJ2)
      SD(NTJ2,NRJ3)=SDYCJ+SD(NTJ2,NRJ3)
      SD(NRJ3,NTJ2)=SDYCJ+SD(NRJ3,NTJ2)
      SD(NTJ3,NRJ1)=SDZAJ+SD(NTJ3,NRJ1)
      SD(NRJ1,NTJ3)=SDZAJ+SD(NRJ1,NTJ3)
      SD(NTJ3,NRJ2)=SDZBJ+SD(NTJ3,NRJ2)
      SD(NRJ2,NTJ3)=SDZBJ+SD(NRJ2,NTJ3)
      SD(NTJ3,NRJ3)=SDZCJ+SD(NTJ3,NRJ3)
      SD(NRJ3,NTJ3)=SDZCJ+SD(NRJ3,NTJ3)
      SD(NTJ1,NRI1)=-SDXAI+SD(NTJ1,NRI1)
      SD(NRI1,NTJ1)=-SDXAI+SD(NRI1,NTJ1)
      SD(NTJ1,NRI2)=-SDXBI+SD(NTJ1,NRI2)
      SD(NRI2,NTJ1)=-SDXBI+SD(NRI2,NTJ1)
      SD(NTJ1,NRI3)=-SDXCI+SD(NTJ1,NRI3)
      SD(NRI3,NTJ1)=-SDXCI+SD(NRI3,NTJ1)
      SD(NTJ2,NRI1)=-SDYAI+SD(NTJ2,NRI1)
      SD(NRI1,NTJ2)=-SDYAI+SD(NRI1,NTJ2)
      SD(NTJ2,NRI2)=-SDYBI+SD(NTJ2,NRI2)
      SD(NRI2,NTJ2)=-SDYBI+SD(NRI2,NTJ2)
      SD(NTJ2,NRI3)=-SDYCI+SD(NTJ2,NRI3)
      SD(NRI3,NTJ2)=-SDYCI+SD(NRI3,NTJ2)
      SD(NTJ3,NRI1)=-SDZAI+SD(NTJ3,NRI1)
      SD(NRI1,NTJ3)=-SDZAI+SD(NRI1,NTJ3)
      SD(NTJ3,NRI2)=-SDZBI+SD(NTJ3,NRI2)
      SD(NRI2,NTJ3)=-SDZBI+SD(NRI2,NTJ3)
      SD(NTJ3,NRI3)=-SDZCI+SD(NTJ3,NRI3)
      SD(NRI3,NTJ3)=-SDZCI+SD(NRI3,NTJ3)
      SD(NTI1,NRJ1)=-SDXAJ+SD(NTI1,NRJ1)
      SD(NRJ1,NTI1)=-SDXAJ+SD(NRJ1,NTI1)
      SD(NTI1,NRJ2)=-SDXBJ+SD(NTI1,NRJ2)
      SD(NRJ2,NTI1)=-SDXBJ+SD(NRJ2,NTI1)
      SD(NTI1,NRJ3)=-SDXCJ+SD(NTI1,NRJ3)
      SD(NRJ3,NTI1)=-SDXCJ+SD(NRJ3,NTI1)
      SD(NTI2,NRJ1)=-SDYAJ+SD(NTI2,NRJ1)
      SD(NRJ1,NTI2)=-SDYAJ+SD(NRJ1,NTI2)
      SD(NTI2,NRJ2)=-SDYBJ+SD(NTI2,NRJ2)
      SD(NRJ2,NTI2)=-SDYBJ+SD(NRJ2,NTI2)
      SD(NTI2,NRJ3)=-SDYCJ+SD(NTI2,NRJ3)
      SD(NRJ3,NTI2)=-SDYCJ+SD(NRJ3,NTI2)
      SD(NTI3,NRJ1)=-SDZAJ+SD(NTI3,NRJ1)
      SD(NRJ1,NTI3)=-SDZAJ+SD(NRJ1,NTI3)
      SD(NTI3,NRJ2)=-SDZBJ+SD(NTI3,NRJ2)
      SD(NRJ2,NTI3)=-SDZBJ+SD(NRJ2,NTI3)
      SD(NTI3,NRJ3)=-SDZCJ+SD(NTI3,NRJ3)
      SD(NRJ3,NTI3)=-SDZCJ+SD(NRJ3,NTI3)
      END IF
      end do
      end do
      end do
   70 CONTINUE
   80 CONTINUE
60    CONTINUE
      RETURN
      END

      subroutine hoqrvw(a,ka,n,e,eps,w,ind)
!     固有値計算。numpacからの流用。
      implicit real*8 (a-h,o-z)
      dimension a(ka,n),e(n),w(n,2)
      data dmach/2.d-16/
      logical ll
      if(n.lt.1.or.n.gt.ka.or.eps.le.0.d0) go to 540
      ll=ind.eq.0
      do 210 k=1,n-2
      md=mod(n-k,4)+1
      e(k)=a(k,k)
      s=0.d0
      do 10 j=k+1,n
      w(j,1)=0.d0
      e(j)=a(k,j)
   10 s=e(j)*e(j)+s
      s=dsign(dsqrt(s),e(k+1))
      w(k,1)=-s
      e(k+1)=e(k+1)+s
      a(k,k+1)=e(k+1)
      h=e(k+1)*s
      if(h.le.0.d0) go to 210
      h=1.d0/h
      do 30 j=k+1,n-3,4
      do 20 i=k+1,n
   20 w(i,1)=e(j)*a(i,j)+e(j+1)*a(i,j+1)
     *+e(j+2)*a(i,j+2)+e(j+3)*a(i,j+3)+w(i,1)
   30 continue
      go to (100,40,60,80),md
   40 do 50 i=k+1,n
   50 w(i,1)=e(n)*a(i,n)+w(i,1)
      go to 100
   60 do 70 i=k+1,n
   70 w(i,1)=e(n-1)*a(i,n-1)+e(n)*a(i,n)+w(i,1)
      go to 100
   80 do 90 i=k+1,n
   90 w(i,1)=e(n-2)*a(i,n-2)+e(n-1)*a(i,n-1)+e(n)*a(i,n)+w(i,1)
  100 r=0.d0
      do 110 i=k+1,n
  110 r=e(i)*w(i,1)+r
      g=r*0.5d0*h
      do 120 j=k+1,n
  120 w(j,1)=(e(j)*g-w(j,1))*h
      do 140 j=k+1,n-3,4
      do 130 i=k+1,n
      a(i,j)=e(j)*w(i,1)+e(i)*w(j,1)+a(i,j)
      a(i,j+1)=e(j+1)*w(i,1)+e(i)*w(j+1,1)+a(i,j+1)
      a(i,j+2)=e(j+2)*w(i,1)+e(i)*w(j+2,1)+a(i,j+2)
  130 a(i,j+3)=e(j+3)*w(i,1)+e(i)*w(j+3,1)+a(i,j+3)
  140 continue
      go to (210,150,170,190),md
  150 do 160 i=k+1,n
  160 a(i,n)=e(n)*w(i,1)+e(i)*w(n,1)+a(i,n)
      go to 210
  170 do 180 i=k+1,n
      a(i,n-1)=e(n-1)*w(i,1)+e(i)*w(n-1,1)+a(i,n-1)
  180 a(i,n)=e(n)*w(i,1)+e(i)*w(n,1)+a(i,n)
      go to 210
  190 do 200 i=k+1,n
      a(i,n-2)=e(n-2)*w(i,1)+e(i)*w(n-2,1)+a(i,n-2)
      a(i,n-1)=e(n-1)*w(i,1)+e(i)*w(n-1,1)+a(i,n-1)
  200 a(i,n)=e(n)*w(i,1)+e(i)*w(n,1)+a(i,n)
  210 a(k,k)=h
      e(n)=a(n,n)
      a(n,n)=1.d0
      if(n.eq.1) go to 530
      e(n-1)=a(n-1,n-1)
      w(n-1,1)=a(n-1,n)
      a(n-1,n-1)=1.d0
      a(n-1,n)=0.d0
      a(n,n-1)=0.d0
      if(ll) go to 430
      do 420 k=n-2,1,-1
      h=-a(k,k)
      a(k,k)=1.d0
      if(h.ge.0.d0) go to 400
      do 220 i=k+1,n
      w(i,2)=0.d0
  220 a(i,k)=a(k,i)*h
      do 240 i=k+1,n-3,4
      do 230 j=k+1,n
  230 w(j,2)=a(i,j)*a(k,i)+a(i+1,j)*a(k,i+1)
     *+a(i+2,j)*a(k,i+2)+a(i+3,j)*a(k,i+3)+w(j,2)
  240 continue
      md=mod(n-k,4)+1
      go to (310,250,270,290),md
  250 do 260 j=k+1,n
  260 w(j,2)=a(n,j)*a(k,n)+w(j,2)
      go to 310
  270 do 280 j=k+1,n
  280 w(j,2)=a(n-1,j)*a(k,n-1)+a(n,j)*a(k,n)+w(j,2)
      go to 310
  290 do 300 j=k+1,n
  300 w(j,2)=a(n-2,j)*a(k,n-2)+a(n-1,j)*a(k,n-1)+a(n,j)*a(k,n)+w(j,2)
  310 do 330 j=k+1,n-3,4
*voption vec
      do 320 i=k+1,n
      a(i,j)=w(j,2)*a(i,k)+a(i,j)
      a(i,j+1)=w(j+1,2)*a(i,k)+a(i,j+1)
      a(i,j+2)=w(j+2,2)*a(i,k)+a(i,j+2)
  320 a(i,j+3)=w(j+3,2)*a(i,k)+a(i,j+3)
  330 continue
      go to (400,340,360,380),md
*voption vec
  340 do 350 i=k+1,n
  350 a(i,n)=w(n,2)*a(i,k)+a(i,n)
      go to 400
*voption vec
  360 do 370 i=k+1,n
      a(i,n-1)=w(n-1,2)*a(i,k)+a(i,n-1)
  370 a(i,n)=w(n,2)*a(i,k)+a(i,n)
      go to 400
*voption vec
  380 do 390 i=k+1,n
      a(i,n-2)=w(n-2,2)*a(i,k)+a(i,n-2)
      a(i,n-1)=w(n-1,2)*a(i,k)+a(i,n-1)
  390 a(i,n)=w(n,2)*a(i,k)+a(i,n)
  400 do 410 i=k+1,n
      a(i,k)=0.d0
  410 a(k,i)=0.d0
  420 continue
  430 gn=dabs(e(1))+dabs(w(1,1))
      w(n,1)=0.d0
      do 440 i=1,n-1
      w(i+1,2)=w(i,1)
  440 gn=dmax1(dabs(w(i,1))+dabs(e(i+1))+dabs(w(i+1,1)),gn)
      if(gn.eq.0.d0) go to 530
      del=gn*dmax1(eps,dmach)
      do 520 k=n,2,-1
  450 do 460 l=k,2,-1
      if(dabs(w(l,2)).lt.del) go to 470
  460 continue
      l=1
  470 if(l.eq.k) go to 520
      g=(e(k-1)+e(k))*0.5d0
      f=e(k)-g
      ww0=w(k,2)*w(k,2)
      z=dsign(dsqrt(f*f+ww0),f)+g
      if(l+1.eq.k) go to 490
      ww1=w(k-1,2)*w(k-1,2)
      do 480 itime=1,2
      f2=(e(k-1)-z)*(e(k)-z)-ww0
      f3=(e(k-2)-z)*f2+(z-e(k))*ww1
      df3=(e(k-2)-z)*(z-e(k)+z-e(k-1))-f2+ww1
      if(dabs(df3).gt.dabs(f3)) z=z-f3/df3
  480 continue
  490 e(l)=e(l)-z
      c=1.d0
      s=0.d0
      do 510 j=l,k-1
      r=dsqrt(e(j)*e(j)+w(j+1,2)*w(j+1,2))
      w(j,2)=s*r
      ee=e(j)*c
      ff=w(j+1,2)*c
      c=e(j)/r
      s=w(j+1,2)/r
      e(j)=((e(j+1)-z)*s+ff*c)*s+ee+z
      e(j+1)=(e(j+1)-z)*c-ff*s
      if(ll) go to 510
      do 500 i=1,n
      g=a(i,j+1)
      a(i,j+1)=g*c-a(i,j)*s
  500 a(i,j)=a(i,j)*c+g*s
  510 continue
      w(k,2)=e(k)*s
      e(k)=e(k)*c+z
      go to 450
  520 continue
      index=1
      if(ll) call sortdk(n,e,index)
      if(.not.ll) call srtvdd(n,e,a,ka,n,w(1,1),w(1,2),index)
  530 ind=0
      return
  540 ind=30000
      return
      end
      subroutine sortds(n,ak,a,ind)
      integer r,rr
      real*8 ak,bk
      dimension a(n),ak(n),rr(30),ll(30)
      entry      sortdi(n,ak,a,ind)
      if(n.lt.1) go to 200
      if(n.eq.1) go to 190
      if(ind.eq.0) go to 20
      do 10 i=1,n
   10 ak(i)=-ak(i)
   20 isp=0
      l=1
      r=n
   30 if(r-l.lt.16) go to 120
      m=(r+l)/2
      max=r
      if(ak(m).gt.ak(r)) max=m
      if(ak(l).gt.ak(max)) max=l
      if(max.eq.r) go to 40
      b=a(max)
      bk=ak(max)
      a(max)=a(r)
      ak(max)=ak(r)
      a(r)=b
      ak(r)=bk
   40 if(ak(l).ge.ak(m)) go to 50
      b=a(l)
      bk=ak(l)
      a(l)=a(m)
      ak(l)=ak(m)
      a(m)=b
      ak(m)=bk
   50 b=a(l)
      bk=ak(l)
      i=l
      j=r
      go to 80
   60 a(j)=a(i)
      ak(j)=ak(i)
   70 j=j-1
   80 if(bk.lt.ak(j)) go to 70
      if(j.le.i) go to 100
      a(i)=a(j)
      ak(i)=ak(j)
   90 i=i+1
      if(ak(i).lt.bk) go to 90
      if(j.gt.i) go to 60
      i=j
  100 a(i)=b
      ak(i)=bk
      isp=isp+1
      if(r-i.ge.i-l) go to 110
      ll(isp)=l
      rr(isp)=i-1
      l=i+1
      go to 30
  110 ll(isp)=i+1
      rr(isp)=r
      r=i-1
      go to 30
  120 if(r-l.lt.1) go to 160
      j=r
  130 b=a(j-1)
      bk=ak(j-1)
      i=j
  140 if(ak(i).ge.bk) go to 150
      a(i-1)=a(i)
      ak(i-1)=ak(i)
      i=i+1
      if(i.le.r) go to 140
  150 a(i-1)=b
      ak(i-1)=bk
      j=j-1
      if(j.gt.l) go to 130
  160 if(isp.eq.0) go to 170
      l=ll(isp)
      r=rr(isp)
      isp=isp-1
      go to 30
  170 if(ind.eq.0) return
      do 180 i=1,n
  180 ak(i)=-ak(i)
      ind=0
      return
  190 ind=0
      return
  200 ind=30000
      return
      end
      subroutine sortdk(n,ak,ind)
      integer*4 rr(30),ll(30),r
      real*8 ak(n),bk
      if(n.lt.1) go to 200
      if(n.eq.1) go to 190
      if(ind.eq.0) go to 20
      do 10 i=1,n
   10 ak(i)=-ak(i)
   20 isp=0
      l=1
      r=n
   30 if(r-l.lt.16) go to 120
      m=(r+l)/2
      max=r
      if(ak(m).gt.ak(r)) max=m
      if(ak(l).gt.ak(max)) max=l
      if(max.eq.r) go to 40
      bk=ak(max)
      ak(max)=ak(r)
      ak(r)=bk
   40 if(ak(l).ge.ak(m)) go to 50
      bk=ak(l)
      ak(l)=ak(m)
      ak(m)=bk
   50 bk=ak(l)
      i=l
      j=r
      go to 80
   60 ak(j)=ak(i)
   70 j=j-1
   80 if(bk.lt.ak(j)) go to 70
      if(j.le.i) go to 100
      ak(i)=ak(j)
   90 i=i+1
      if(ak(i).lt.bk) go to 90
      if(j.gt.i) go to 60
      i=j
  100 ak(i)=bk
      isp=isp+1
      if(r-i.ge.i-l) go to 110
      ll(isp)=l
      rr(isp)=i-1
      l=i+1
      go to 30
  110 ll(isp)=i+1
      rr(isp)=r
      r=i-1
      go to 30
  120 if(r-l.lt.1) go to 160
      j=r
  130 bk=ak(j-1)
      i=j
  140 if(ak(i).ge.bk) go to 150
      ak(i-1)=ak(i)
      i=i+1
      if(i.le.r) go to 140
  150 ak(i-1)=bk
      j=j-1
      if(j.gt.l) go to 130
  160 if(isp.eq.0) go to 170
      l=ll(isp)
      r=rr(isp)
      isp=isp-1
      go to 30
  170 if(ind.eq.0) return
      do 180 i=1,n
  180 ak(i)=-ak(i)
      ind=0
      return
  190 ind=0
      return
  200 ind=30000
      return
      end
      subroutine srtvdd(n,ak,v,kv,l,iw,w,ind)
      real*8 ak(n),v(kv,n),w(l)
      dimension iw(n)
      if(l.lt.1.or.n.lt.1.or.l.gt.kv) go to 90
      do 10 i=1,n
   10 iw(i)=i
      call sortdi(n,ak,iw,ind)
      do 70 k=1,n
      if(iw(k).eq.k) go to 70
      do 20 i=1,l
   20 w(i)=v(i,k)
      j=k
   30 jj=iw(j)
      iw(j)=j
      if(jj.eq.k) go to 50
      do 40 i=1,l
   40 v(i,j)=v(i,jj)
      j=jj
      go to 30
   50 do 60 i=1,l
   60 v(i,j)=w(i)
   70 continue
   80 ind=0
      return
   90 ind=30000
      return
      end

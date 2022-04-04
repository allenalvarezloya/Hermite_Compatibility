      subroutine abc(N,NE,LENC,IRN,A,IP,DW,CPERM)
      INTEGER MAXN, MAXNE, LIW, LDW
      PARAMETER( MAXN=100, MAXNE=MAXN*MAXN,
     &           LIW=10*MAXN+MAXNE, LDW=3*MAXN+MAXNE )
      INTEGER ICNTL(10), INFO(10), I, J, JOB, K, N, NE, NUM,
     &        IRN(MAXNE), LENC(MAXN), IP(MAXN),
     &        RPERM(MAXN), CPERM(MAXN), IW(LIW)
      DOUBLE PRECISION A(MAXNE), DW(LDW)
      EXTERNAL MC64AD, MC22AD
      INTRINSIC EXP
C Read matrix and set column pointers
C       READ(5,*) N, NE
C       IP(1) = 1
C       DO 10 J=1,N
C         READ(5,*) LENC(J),
C      &            (IRN(K),A(K),K=IP(J),IP(J)+LENC(J)-1)
C         IP(J+1) = IP(J) + LENC(J)
C    10 CONTINUE
C Set default values for ICNTL
      CALL MC64ID(ICNTL)
C Suppress error and warning messages
      ICNTL(1) = -1
      ICNTL(2) = -1
C Compute matchings
       JOB = 5
       CALL MC64AD(JOB,N,NE,IP,IRN,A,NUM,CPERM,LIW,IW,LDW,DW,ICNTL,INFO)
        IF (INFO(1).LT.0) WRITE(6,'(A,I2)')
     &    ' Error return from MC64A/AD, INFO(1) = ',INFO(1)
C         WRITE(6,'(A,I1,A,4I2)')
C      &    ' For JOB = ',JOB,' the array CPERM() is: ',(CPERM(J),J=1,N)
C DW(1:N) contains row scaling, DW(N+1:2N) contains column scaling
C Scale the matrix
      DO 25 J = 1,N
        DO 24 K = IP(J),IP(J+1)-1
          I = IRN(K)
          A(K) = A(K) * EXP(DW(I) + DW(N+J))
   24   CONTINUE
   25 CONTINUE
C Put identity permutation in RPERM
      DO 30 J = 1,N
        RPERM(J) = J
   30 CONTINUE
C RPERM contains row permutation, CPERM contains column permutation
C Permute the matrix with MC22A/AD
      CALL MC22AD(N,IRN,A,NE,LENC,CPERM,RPERM,IW(1),IW(2*N+1))
C Adjust column pointers according to the permutation
      IP(1) = 1
      DO 40 J=1,N
        IP(J+1) = IP(J) + LENC(J)
   40 CONTINUE
C Write matrix
C       WRITE(6,'(A)') ' The scaled and permuted matrix (JOB=5) is:'
C       WRITE(6,'(2I4)') N, NE
C       DO 50 J=1,N
C         WRITE(6,'(I4,4(I6,F6.3))') LENC(J),
C      &           (IRN(K),A(K),K=IP(J),IP(J)+LENC(J)-1)
C    50 CONTINUE
      end subroutine


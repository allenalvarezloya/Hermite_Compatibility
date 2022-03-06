subroutine eq(A,N,NE,LENC_number,ROW,COL,PERM) 
	implicit none 
	INTEGER MAXN, MAXNE, LIW, LDW, LENC_number
      PARAMETER( MAXN=100, MAXNE=MAXN*MAXN,&
                LIW=10*MAXN+MAXNE, LDW=3*MAXN+MAXNE )
      INTEGER ICNTL(10), INFO(10), I, J, JOB, K, N, NE, NUM,&
             IRN(MAXNE), LENC(MAXN), IP(MAXN), &
             RPERM(MAXN), CPERM(MAXN), IW(LIW)
      DOUBLE PRECISION A(MAXNE), DW(LDW), ROW(N), COL(N), PERM(N)
      INTEGER COUNT
    COUNT = 1
    LENC = LENC_number
    !  Autofill IRN 
    DO I = 1,N 
      DO J = 1,N 
        IRN(COUNT) = J 
        COUNT = COUNT + 1
      end DO 
    end DO
    COUNT = 1
    DO J = 1,N+1 
      IP(J) = COUNT
      COUNT = COUNT + N 
    end DO 
	call abc(N,NE,LENC,IRN,A,IP,DW,CPERM)
	! write(*,*) " "
	! write(*,*) " Output in main "
 !      WRITE(*,'(A)') ' The scaled and permuted matrix (JOB=5) is:'
 !      WRITE(*,'(2I4)') N, NE
 !      DO J=1,N
 !        WRITE(*,'(I4,4(I6,F6.3))') LENC(J), &
 !                (IRN(K),A(K),K=IP(J),IP(J)+LENC(J)-1)
	!   end DO
    DO J=1,N
      ROW(J) = DW(J)
      COL(J) = DW(N+J)
      PERM(J) = CPERM(J)
    END DO
end subroutine
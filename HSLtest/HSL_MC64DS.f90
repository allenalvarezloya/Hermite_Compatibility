PROGRAM HSL_MC64DS
	USE HSL_MC64_DOUBLE
	IMPLICIT NONE
	TYPE(MC64_CONTROL) CONTROL
	TYPE(MC64_INFO) INFO
	INTEGER :: MATRIX_TYPE
	INTEGER, ALLOCATABLE :: PTR(:), ROW(:)
	DOUBLE PRECISION, ALLOCATABLE :: VAL(:)
	INTEGER, ALLOCATABLE :: PERM(:)
	INTEGER EYE, I, J, JAY, JOB, K, KASE, M, N, NE
	INTEGER, PARAMETER :: WP = KIND(0.0D0)
	REAL(WP), ALLOCATABLE :: SCALE(:),A(:,:)
	REAL(WP), PARAMETER :: ZERO=0.0D0
	DO KASE = 1,3
		! KASE = 1 ... square matrix
		! KASE = 2 ... rectangular matrix
		! KASE = 3 ... symmetric matrix
		IF (KASE == 1) THEN
			WRITE(6,’(A)’) ’Square matrix’
			MATRIX_TYPE = 2
		ENDIF
		IF (KASE == 2) THEN
			WRITE(6,’(//A)’) ’Rectangular matrix’
			MATRIX_TYPE = 1
		ENDIF
		IF (KASE == 3) THEN
			WRITE(6,’(//A)’) ’Symmetric matrix’
			MATRIX_TYPE = 4
		ENDIF
		! Read matrix order and number of entries
		READ (5,*) M,N,NE
		! Allocate arrays of appropriate sizes
		ALLOCATE(PTR(N+1), VAL(NE), ROW(NE))
		ALLOCATE(PERM(M+N),SCALE(M+N))
		ALLOCATE(A(M,N))
		! Read matrix and right-hand side
		READ (5,*) (PTR(I),I=1,N+1)
		READ (5,*) (ROW(I),I=1,NE)
		READ (5,*) (VAL(I),I=1,NE)
		! Get matching
		JOB = 5
		CALL MC64_MATCHING(JOB,MATRIX_TYPE,M,N,PTR,ROW,VAL,CONTROL,INFO,PERM,SCALE)
		IF(INFO%FLAG<0) THEN
			WRITE(6,’(A,I2)’) &
			’ Failure of MC64_MATCHING with INFO%FLAG=’, INFO%FLAG
			STOP
		END IF
		! Print permutations
		WRITE(6,’(A/(10I8))’) ’Row permutation’,PERM(1:M)
		WRITE(6,’(A/(10I8))’) ’Column permutation’,PERM(M+1:M+N)
		IF (MATRIX_TYPE.EQ.4) THEN
			! Print symmetric scaling
			WRITE(6,’(A/(5D12.4))’) ’Symmetric scaling’,EXP(SCALE(1:N))
		ELSE
			! Print row scaling
			WRITE(6,’(A/(5D12.4))’) ’Row scaling’,EXP(SCALE(1:M))
			! Print column scaling
			WRITE(6,’(A/(5D12.4))’) ’Column scaling’,EXP(SCALE(M+1:M+N))
		ENDIF
		! Scale matrix and put scaled and permuted entries in full array
		A(1:M,1:N) = ZERO
		DO J = 1,N
			! EYE and JAY are indices in permuted matrix
			JAY = ABS(PERM(M+J))
			DO K = PTR(J), PTR(J+1)-1
				I = ROW(K)
				EYE = ABS(PERM(I))
				! Scale and permute the matrix
				A(EYE,JAY) = EXP(SCALE(I))*VAL(K)*EXP(SCALE(M+J))
				! Set symmetric counterpart if appropriate
				If (MATRIX_TYPE.GE.3) A(JAY,EYE) = A(EYE,JAY)
			ENDDO
		ENDDO
		! Write out scaled matrix
		WRITE(6,’(A)’) ’Scaled matrix’
		DO I = 1,M
			WRITE(6,’(3D12.4)’) (A(I,J),J=1,N)
		ENDDO
	ENDDO
END PROGRAM HSL_MC64DS
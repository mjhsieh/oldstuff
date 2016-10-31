! I intended to do fasta reader in an hour, but this reminds me of hell of fortran i/o
module blah
use, intrinsic :: ISO_FORTRAN_ENV, only: IOSTAT_EOR
contains
  subroutine GetLine(unit, line, stat, iomsg)
    ! Reads a complete line (end-of-record terminated) from a file.
    !---------------------------------------------------------------------------
    ! Passed variables
    INTEGER, INTENT(IN) :: unit
    !! Logical unit connected for formatted input to the file.
    CHARACTER(:), INTENT(OUT), ALLOCATABLE :: line
    !! The line read.
    INTEGER, INTENT(OUT) :: stat
    !! Error code, positive on error, IOSTAT_END (which is negative) on end of
    !! file.
    CHARACTER(*), INTENT(OUT) :: iomsg
    !! Error message - only defined if iostat is non-zero.
    !---------------------------------------------------------------------------
    ! Local variables
    ! Buffer to read the line (or partial line).
    CHARACTER(256) :: buffer
    INTEGER :: size           ! Number of characters read from the file.
    !***************************************************************************
    line = ''
    do
       read (unit, "(A)", ADVANCE='NO', IOSTAT=stat, IOMSG=iomsg, SIZE=size)  &
            buffer
       if (stat > 0) RETURN      ! Some sort of error.
       line = line // buffer(:size)
       if (stat < 0) THEN
          if (stat == IOSTAT_EOR) stat = 0
          return
       end if
    end do
  end subroutine GetLine
end module blah

program main
  use blah
  implicit none
  integer             :: IArgC,i,iostat
  integer             :: mynarg, FASTA_IO
  character(len=1000) :: str, iomsg
  CHARACTER(:), ALLOCATABLE :: line
  LOGICAL doesExist
  str = ''
  mynarg = IArgC()
  if ( mynarg .ne. 1 ) then
     print *,mynarg
     stop
  endif
  do i = 1, mynarg
     call GetArg(i,str)
     ! write(6,*) "test", len(trim(str))
  enddo
  inquire( FILE=trim(str), EXIST=doesExist )
  if ( .not. doesExist ) then
     write(6,*) "file not found"
     stop(1)
  endif
  open(unit=FASTA_IO, file=trim(str), status='old', action='read')
  do
     call GetLine(FASTA_IO,line,iostat,iomsg)
     if ( iostat .lt. 0 ) exit
     write(*,"(A)") line
  enddo
  close(FASTA_IO)
  stop
1000 continue
  write(6,'(/'' FASTA file is not found: '',a)') trim(str)
  close(FASTA_IO)
  stop
end program main

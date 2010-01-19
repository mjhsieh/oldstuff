C
C   This utility convert CHARMM parameter file into MDynamix 
C   parameter file (.ff) (understandable for makemol utility)
C
      program char2mdx
      character*80 str,aux,comment
      character*4 N1,N2,N3,N4
      ICONV=0
 10   read(*,'(a80)',end=200)str
      l=leng(str,80)
      if(str(1:1).eq.'*'.or.str(1:1).eq.'!'.or.l.le.1)then
         if(l.le.1)l=1
         str(1:1)='#'
         write(*,'(a)')str(1:l)
         go to 100
      end if
      if(str(1:5).eq.'BONDS')then
         write(*,'(a5)')'BONDS'
         ICONV=1
         go to 100
      end if
      if(str(1:5).eq.'ANGLE')then
         write(*,'(a6)')'ANGLES'
         ICONV=2
         go to 100
      end if
      if(str(1:5).eq.'DIHED')then
         write(*,'(a8)')'TORSIONS'
         ICONV=3
         go to 100
      end if
      if(str(1:5).eq.'IMPRO')then
         write(*,'(a)')'#   Improper angles not implemented in Makemol'
         ICONV=4
         go to 100
      end if
      if(str(1:5).eq.'NONBO')then
         write(*,'(a9)')'NONBONDED'
         ICONV=5
         go to 100
      end if
      aux=comment(str,l)
      la=leng(aux,80)
      if(ICONV.eq.1)then
         read(str,*,end=100,err=100)N1,N2,FB,RB
         FB=4.187*FB
         write(*,'(a4,2x,a4,2x,f7.4,2x,f9.1,2x,a)')N1,N2,RB,FB,aux(1:la)
         go to 100
      end if
      if(ICONV.eq.2)then
         read(str,*,end=50,err=50)N1,N2,N3,FA,RA,FU,RU
         FU=4.187*FU
         FA=4.187*FA
      write(*,'(a4,2x,a4,2x,a4,2x,f7.2,1x,f8.1,1x,f8.4,1x,f8.2,2x,a)')
     +   N1,N2,N3,RA,FA,RU,FU,aux(1:la)
         go to 100
 50      read(str,*,end=100,err=100)N1,N2,N3,FA,RA
         FA=4.187*FA
         write(*,'(a4,2x,a4,2x,a4,2x,f7.2,2x,f9.1,2x,a)')
     +   N1,N2,N3,RA,FA,aux(1:la)
         go to 100
      end if
      if(ICONV.eq.3)then
         read(str,*,end=100,err=100)N1,N2,N3,N4,FT,M,RT
         FT=4.187*FT
         write(*,'(a4,2x,a4,2x,a4,2x,a4,2x,f8.0,2x,f8.1,2x,i2,2x,a)')
     +   N1,N2,N3,N4,RT,FT,M,aux(1:la)
         go to 100
      end if
      if(ICONV.eq.4)go to 100
      if(ICONV.eq.5)then
         read(str,*,end=60,err=60)N1,DUM,EP,SI,DUMM,EP14,SI14
         SI=1.782*SI
         SI14=1.782*SI14
         EP=-4.187*EP
         EP14=-4.187*EP14
         write(*,'(a4,2x,4(f8.4,2x),a)')
     +   N1,SI,EP,SI14,EP14,aux(1:la)
         go to 100
 60      read(str,*,end=100,err=100)N1,DUM,EP,SI
         SI=1.782*SI
         EP=4.187*EP
         write(*,'(a4,2x,2(f8.4,2x),a)')
     +   N1,SI,EP,aux(1:la)
         go to 100
      end if
 100  do I=1,80
         str(I:I)=' '
      end do
      go to 10
 200  write(*,'(a3)')'END'
      stop
      end
*
*================ LENG ===============================================
*
	function LENG(STR,LS)
	character*1 STR(LS)
	IL=0 
	do I=1,LS
	  ISTR=ichar(STR(I))
	  if(ISTR.le.15)then
	    LENG=I-1
	    return
	  end if
        if(STR(I).ne.' ')IL=I
	end do 
	LENG=IL
	return
      end    
*
*============== COMMENT =================================
*
      character*80 function comment(str,L)
      character*80 str
      do J=1,80
         comment(J:J)=' '
      end do
      do I=1,L
         if(str(I:I).eq.'!')then
            LC=L-I
            comment(1:LC)=str(I:L)
            return
         end if
      end do
      return
      end

MODULE aloha_scattering_matrix_utils

    USE aloha_constants

    Implicit none

CONTAINS

    ! ***
    ! Surotine de CONTrole de la SYMetrie du module d une mat. de scattering
    ! ***
    Subroutine ContSym1(Dim,Spr,Statsym)

    Integer          :: i,j,Dim
    Integer          :: Statsym(Dim,Dim)
    !Real(kind=wp)    :: A21,A12!,Diff
    Real(kind=wp)    :: A21r,A12r,Diffr
    Real(kind=wp)    :: A21i,A12i,Diffi
    Complex(kind=wp) :: Spr(Dim,Dim)

    Do j=1,Dim
        Do i=1,Dim
            If (i /= j) Then
    !               A12=CDABS(Spr(i,j))
    !               A21=CDABS(Spr(j,i))
    !               Diff=DABS(A12-A21)
                A12r=real(Spr(i,j))
                A21r=real(Spr(j,i))
                A12i=aimag(Spr(i,j))
                A21i=aimag(Spr(j,i))
                Diffr=DABS(A12r-A21r)
                Diffi=DABS(A12i-A21i)
    !               If (Diff.LT.0.001) Then
                If ((Diffr < 0.001).AND.(Diffi < 0.001)) Then
                    Statsym(i,j)=0
                Else 
    !               If ((Diffr.GT.0.001).OR.(Diffi.GT.0.001)) Then
    !               If (Diff.GT.0.001) Then
                    Statsym(i,j)=-1
    !                  Stop 'Probleme de Symetrie de la matice '
                EndIf
            EndIf 
        end do
    end do

    End subroutine


! ***
    Subroutine ContSymgg(Dim,S,Statsymgg)

    Integer          :: i,j,Dim
    Integer          :: Statsymgg(Dim,Dim)
    Real(kind=wp)    :: A21r,A12r,Diffr
    Real(kind=wp)    :: A21i,A12i,Diffi
    Complex(kind=wp) :: S(Dim,Dim)

    Do j=1,Dim
        Do i=1,Dim
            If (i /= j) Then
                A12r=real(S(i,j))
                A21r=real(S(j,i))
                A12i=aimag(S(i,j))
                A21i=aimag(S(j,i))
                Diffr=DABS(A12r-A21r)
                Diffi=DABS(A12i-A21i)
                If ((Diffr < 0.001).AND.(Diffi < 0.001)) Then
                    Statsymgg(i,j)=0
                Else 
    !               If ((Diffr.GT.0.001).OR.(Diffi.GT.0.001)) Then
                    Statsymgg(i,j)=-1
    !                  Stop 'Probleme de Symetrie de la matice '
                EndIf
            EndIf 
        end do
    end do

    End subroutine




! ***
! Subroutine de CONTrole de la SYMetrie du module d une mat. de scattering
! ***
      Subroutine ContSymg(Dim,Sg,Statsymgg)


!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!cccc
!c modif altair.partenaires

      Integer *4 i,j,Dim
      Integer *4 Statsymgg(2*10*60*5,2*10*60*5)
      Real *8 A21r,A12r,Diffr
      Real *8 A21i,A12i,Diffi
      Complex *16 Sg(2*10*60*5,2*10*60*5)
      !Character *1 rep

!&&&      Integer *4 i,j,Dim
!&&&      Integer *4 Statsymgg(2*3*15*5,2*3*15*5)
!&&&      Real *8 A21r,A12r,Diffr
!&&&      Real *8 A21i,A12i,Diffi
!&&&      Complex *16 Sg(2*3*15*5,2*3*15*5)
!&&&      Character *1 rep

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!cccc

      Do 100 j=1,Dim
         Do 110 i=1,Dim
            If (i.NE.j) Then
               A12r=DBLE(Sg(i,j))
               A21r=DBLE(Sg(j,i))
               A12i=DIMAG(Sg(i,j))
               A21i=DIMAG(Sg(j,i))
               Diffr=DABS(A12r-A21r)
               Diffi=DABS(A12i-A21i)
               If ((Diffr.LT.0.001).AND.(Diffi.LT.0.001)) Then
                  Statsymgg(i,j)=0
               Else 
!               If ((Diffr.GT.0.001).OR.(Diffi.GT.0.001)) Then
                  Statsymgg(i,j)=-1
!                  Stop 'Probleme de Symetrie de la matice '
               EndIf
            EndIf 
 110     Continue
 100  Continue
  
      Return
      End subroutine

! ***
! Suroutine de CONTrole de l UNItarite d une matrice de scattering
! ***
      Subroutine ContUnig(Dim,Sprg,Statuni)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,Dim
      Integer *4 Statuni(10*60*5)
      Real *8 Somm,Som(10*60*5)!,Diff
      Complex *16 Sprg(10*60*5,10*60*5)
      !Character *1 rep

!&&&      Integer *4 i,j,Dim
!&&&      Integer *4 Statuni(3*15*5)
!&&&      Real *8 Diff,Somm,Som(3*15*5)
!&&&      Complex *16 Sprg(3*15*5,3*15*5)
!&&&      Character *1 rep

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc

      Do 100 j=1,Dim
         Somm=0.
         Do 110 i=1,Dim
            Somm=Somm+(CDABS(Sprg(i,j)))**2
 110     Continue
         Som(j)=Somm
!         Diff=DABS(Somm-1)
!         If (Diff.LT.0.01) Then
         If (Somm.LT.1.) Then
            Statuni(i)=0
         Else 
             Statuni(i)=-1
!             Write (*,*) 'Probleme d Unitarite '
         EndIf
 100  Continue
  
!      Write (*,1) (Som(j),j=1,Dim)
!   1  Format(F8.4) 
      If (Dim.EQ.1) Then
         Write (*,1) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.2) Then
         Write (*,2) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.3) Then
         Write (*,3) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.4) Then
         Write (*,4) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.5) Then
         Write (*,5) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.6) Then
         Write (*,6) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.7) Then
         Write (*,7) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.8) Then
         Write (*,8) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.9) Then
         Write (*,9) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.10) Then
         Write (*,10) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.11) Then
         Write (*,11) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.12) Then
         Write (*,12) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.13) Then
         Write (*,13) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.14) Then
         Write (*,14) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.15) Then
         Write (*,15) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.16) Then
         Write (*,16) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.17) Then
         Write (*,17) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.18) Then
         Write (*,18) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.19) Then
         Write (*,19) (Som(j),j=1,Dim)
      ElseIf (Dim.EQ.20) Then
         Write (*,20) (Som(j),j=1,Dim)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4))
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4))  
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 

      Return
      End subroutine

! ***
! Subroutine de CALcul de la MATrice d AMPLItude d une matrice complexe
! ***
      Subroutine CalMatAmplig(prhe,Sprg,Matamp)


!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 Matamp(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

!&&&      Integer *4 i,j,prhe
!&&&      Real *8 Matamp(2*3/2*15*5,2*3/2*15*5)
!&&&      Complex *16 Sprg(2*3/2*15*5,2*3/2*15*5)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!cccc

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matamp(i,j)=CDABS(Sprg(i,j))
 110     Continue
 100  Continue

      Return
      End subroutine



! ***
      Subroutine CalMatAmpligg(prhe,Sg,Matampg)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 Matampg(2*10*60*5,2*10*60*5)
      Complex *16 Sg(2*10*60*5,2*10*60*5)

!&&&      Integer *4 i,j,prhe
!&&&      Real *8 Matampg(2*3*15*5,2*3*15*5)
!&&&      Complex *16 Sg(2*3*15*5,2*3*15*5)


!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            Matampg(i,j)=CDABS(Sg(i,j))
 110     Continue
 100  Continue

      Return
      End subroutine

! ***
! Subroutine de CALcul de la MATrice de PHASe d une matrice complexe

      Subroutine CalMatPhasg(prhe,Sprg,Matpha)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,prhe
      Real *8 pReel,pIma,phi,pi
      Real *8 Matpha(2*10/2*60*5,2*10/2*60*5)
      Complex *16 Sprg(2*10/2*60*5,2*10/2*60*5)

!&&&      Integer *4 i,j,prhe
!&&&      Real *8 pReel,pIma,phi,pi
!&&&      Real *8 Matpha(2*3/2*15*5,2*3/2*15*5)
!&&&      Complex *16 Sprg(2*3/2*15*5,2*3/2*15*5)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc

      pi=4.*ATAN(1.)

      Do 100 j=1,prhe
         Do 110 i=1,prhe
            pReel=DBLE(Sprg(i,j))
            pIma=DIMAG(Sprg(i,j))
            If (pReel.EQ.0) Then
               If (pIma.GT.0) Then
                  Matpha(i,j)=90.
               ElseIf (pIma.LT.0) Then
                  Matpha(i,j)=-90.
               Else
                  Matpha(i,j)=0.
               EndIf
            Else
               phi=pIma/pReel
               If (pReel.GT.0) Then
                  Matpha(i,j)=(180./pi)*DATAN(phi)
               Else
                  Matpha(i,j)=180.+(180./pi)*DATAN(phi)
                  If (Matpha(i,j).GT.180.) Then
                     Matpha(i,j)=Matpha(i,j)-360.
                  EndIf
               EndIf
            EndIf 
 110     Continue
 100  Continue

      Return
      End subroutine

! *** 
! Subroutine pour IMPrimer une MATrice Reelle d Amplitude
! ***
      Subroutine ImpMatRAg(lprhe,cprhe,mat)


!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

!&&&      Integer *4 i,j,cprhe,lprhe
!&&&      Real *8 mat(2*3/2*15*5,2*3/2*15*5)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccccc

      If (cprhe.EQ.1) Then
         Write (*,1) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4)) 
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4)) 
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 
  21  Format(21(F8.4)) 
  22  Format(22(F8.4)) 
  23  Format(23(F8.4)) 
  24  Format(24(F8.4)) 
  25  Format(25(F8.4)) 
  26  Format(26(F8.4)) 
  27  Format(27(F8.4)) 
  28  Format(28(F8.4)) 
  29  Format(29(F8.4)) 
  30  Format(30(F8.4)) 

      Return
      End subroutine
! ***
      Subroutine ImpMatRAgg(lprhe,cprhe,matg)


!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 matg(2*10*60*5,2*10*60*5)

!&&&      Integer *4 i,j,cprhe,lprhe
!&&&      Real *8 matg(2*3*15*5,2*3*15*5)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!cccc

      If (cprhe.EQ.1) Then
         Write (*,1) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((matg(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.4)) 
   2  Format(2(F8.4)) 
   3  Format(3(F8.4)) 
   4  Format(4(F8.4)) 
   5  Format(5(F8.4)) 
   6  Format(6(F8.4)) 
   7  Format(7(F8.4)) 
   8  Format(8(F8.4)) 
   9  Format(9(F8.4)) 
  10  Format(10(F8.4)) 
  11  Format(11(F8.4)) 
  12  Format(12(F8.4)) 
  13  Format(13(F8.4)) 
  14  Format(14(F8.4)) 
  15  Format(15(F8.4)) 
  16  Format(16(F8.4)) 
  17  Format(17(F8.4)) 
  18  Format(18(F8.4)) 
  19  Format(19(F8.4)) 
  20  Format(20(F8.4)) 
  21  Format(21(F8.4)) 
  22  Format(22(F8.4)) 
  23  Format(23(F8.4)) 
  24  Format(24(F8.4)) 
  25  Format(25(F8.4)) 
  26  Format(26(F8.4)) 
  27  Format(27(F8.4)) 
  28  Format(28(F8.4)) 
  29  Format(29(F8.4)) 
  30  Format(30(F8.4)) 

      Return
      End subroutine

! Subroutine pour IMPrimer une MATrice Reelle de Phase
! ***
      Subroutine ImpMatRPg(lprhe,cprhe,mat)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc
!c modif altair.partenaires

      Integer *4 i,j,cprhe,lprhe
      Real *8 mat(2*10/2*60*5,2*10/2*60*5)

!&&&      Integer *4 i,j,cprhe,lprhe
!&&&      Real *8 mat(2*3/2*15*5,2*3/2*15*5)

!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc!ccc

      If (cprhe.EQ.1) Then
         Write (*,1) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.2) Then
         Write (*,2) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.3) Then
         Write (*,3) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.4) Then
         Write (*,4) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.5) Then
         Write (*,5) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.6) Then
         Write (*,6) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.7) Then
         Write (*,7) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.8) Then
         Write (*,8) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.9) Then
         Write (*,9) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.10) Then
         Write (*,10) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.11) Then
         Write (*,11) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.12) Then
         Write (*,12) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.13) Then
         Write (*,13) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.14) Then
         Write (*,14) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.15) Then
         Write (*,15) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.16) Then
         Write (*,16) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.17) Then
         Write (*,17) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.18) Then
         Write (*,18) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.19) Then
         Write (*,19) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.20) Then
         Write (*,20) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.21) Then
         Write (*,21) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.22) Then
         Write (*,22) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.23) Then
         Write (*,23) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.24) Then
         Write (*,24) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.25) Then
         Write (*,25) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.26) Then
         Write (*,26) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.27) Then
         Write (*,27) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.28) Then
         Write (*,28) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.29) Then
         Write (*,29) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      ElseIf (cprhe.EQ.30) Then
         Write (*,30) ((mat(i,j),j=1,cprhe),i=1,lprhe)
      EndIf

   1  Format(1(F8.2)) 
   2  Format(2(F8.2)) 
   3  Format(3(F8.2)) 
   4  Format(4(F8.2)) 
   5  Format(5(F8.2)) 
   6  Format(6(F8.2)) 
   7  Format(7(F8.2)) 
   8  Format(8(F8.2)) 
   9  Format(9(F8.2)) 
  10  Format(10(F8.2)) 
  11  Format(11(F8.2)) 
  12  Format(12(F8.2)) 
  13  Format(13(F8.2)) 
  14  Format(14(F8.2)) 
  15  Format(15(F8.2)) 
  16  Format(16(F8.2)) 
  17  Format(17(F8.2)) 
  18  Format(18(F8.2)) 
  19  Format(19(F8.2)) 
  20  Format(20(F8.2)) 
  21  Format(21(F8.2)) 
  22  Format(22(F8.2)) 
  23  Format(23(F8.2)) 
  24  Format(24(F8.2)) 
  25  Format(25(F8.2)) 
  26  Format(26(F8.2)) 
  27  Format(27(F8.2)) 
  28  Format(28(F8.2)) 
  29  Format(29(F8.2)) 
  30  Format(30(F8.2)) 

      Return
      End subroutine

END MODULE aloha_scattering_matrix_utils
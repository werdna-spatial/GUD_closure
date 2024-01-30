#ifdef ALLOW_GUD

CBOP
C    !ROUTINE: GUD_DIAGS.h
C    !INTERFACE:
C #include GUD_DIAGS.h

C    !DESCRIPTION:
C Contains indices into diagnostics array

      integer iPP
      integer iNfix
      integer iDenit
      integer iDenitN
      integer iPPplank
      integer iGRplank
      integer iConsDIN
      integer iConsPO4
      integer iConsSi
      integer iConsFe
      integer gud_nDiag
      integer iUTK_holl
      integer iUTK_2ZPP
      integer iUTK_C
      integer iUTK_D
      integer iUTK_E
      PARAMETER(iPP=     1)
      PARAMETER(iNfix=   2)
      PARAMETER(iDenit=  3)
      PARAMETER(iDenitN= 4)
      PARAMETER(iConsPO4=5)
      PARAMETER(iConsSi= 6)
      PARAMETER(iConsFe= 7)
      PARAMETER(iConsDIN=8)
      PARAMETER(iPPplank=9)
      PARAMETER(iUTK_holl=10)
      PARAMETER(iUTK_2ZPP=11)
      PARAMETER(iUTK_C=12)
      PARAMETER(iUTK_D=13)
      PARAMETER(iUTK_E=14)
      PARAMETER(iGRplank=iUTK_E+nPPplank)
      PARAMETER(gud_nDiag=iGRplank+nGRplank-1)

CEOP
#endif /* ALLOW_GUD */

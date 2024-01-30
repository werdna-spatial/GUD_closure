#include "GUD_OPTIONS.h"

CBOP
C     !ROUTINE: GUD_DIAGNOSTICS_INIT

C     !INTERFACE:
      SUBROUTINE GUD_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Routine to initialize Gud diagnostics

C     !USES:
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "GUD_SIZE.h"
#include "GUD_INDICES.h"

C     !INPUT PARAMETERS:
C     myThid :: Number of this instance of GUD_DIAGNOSTICS_INIT
      INTEGER myThid
CEOP

#ifdef ALLOW_GUD
#ifdef ALLOW_DIAGNOSTICS

C     !FUNCTIONS:
      INTEGER     ILNBLNK
      EXTERNAL    ILNBLNK
      CHARACTER*16 DIAGS_MK_UNITS
      EXTERNAL     DIAGS_MK_UNITS

C     !LOCAL VARIABLES:
C     msgBuf :: Informational/error meesage buffer
c     CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER j
      INTEGER        diagNum
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle
      INTEGER      ilnb
      CHARACTER*7  trcUnits
      CHARACTER*30 locName
      CHARACTER*4  diagSufx

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      WRITE(diagName,'(A)')'PP'
      WRITE(diagTitle,'(A)')'Primary Production'
      diagUnits = 'mmol C/m^3/s    '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'Nfix'
      WRITE(diagTitle,'(A)')'N fixation'
      diagUnits = 'mmol N/m^3/s    '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'Denit'
      WRITE(diagTitle,'(A)')'Denitrification'
      diagUnits = 'mmol N/m^3/s    '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef GUD_ALLOW_RADTRANS
      WRITE(diagName,'(A)')'rmud'
      WRITE(diagTitle,'(A)')'inverse cosine of solar zenith angle'
      diagUnits = '1               '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      DO j = 1, nlam
        WRITE(diagName,'(A,I3.3)')'Rirr', j
        WRITE(diagTitle,'(A, I3)')
     &        'irradiance reflectance for waveband ', j
        diagUnits = '1               '
        diagCode  = 'SM P    L1      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'Ed',j
        WRITE(diagTitle,'(A,I3)')'direct irradiance for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'WM P    LR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'Es',j
        WRITE(diagTitle,'(A,I3)')
     &        'diffuse downwelling irradiance for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'WM P    LR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'Eu',j
        WRITE(diagTitle,'(A,I3)')
     &        'diffuse upwelling irradiance for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'WM P    LR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'PAR',j
        WRITE(diagTitle,'(A,I3)')'PAR waveband ',j
        diagUnits = 'uEin/m^2/s      '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'PARF',j
        WRITE(diagTitle,'(A,I3)')'PAR at W point, waveband ',j
        diagUnits = 'uEin/m^2/s      '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef GUD_DIAG_IOP
        WRITE(diagName,'(A,I3.3)')'a',j
        WRITE(diagTitle,'(A,I3)')
     &        'total absorption for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'bt',j
        WRITE(diagTitle,'(A,I3)')
     &        'total scattering for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'bb',j
        WRITE(diagTitle,'(A,I3)')
     &        'total backscattering for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'aplk',j
        WRITE(diagTitle,'(A,I3)')
     &        'absorption by plankton for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'btplk',j
        WRITE(diagTitle,'(A,I3)')
     &        'scattering by plankton for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'bbplk',j
        WRITE(diagTitle,'(A,I3)')
     &        'backscattering by plankton for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'aprt',j
        WRITE(diagTitle,'(A,I3)')
     &        'absorption by particles for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'btprt',j
        WRITE(diagTitle,'(A,I3)')
     &        'scattering by particles for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'bbprt',j
        WRITE(diagTitle,'(A,I3)')
     &        'backscattering by particles for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'aCDOM',j
        WRITE(diagTitle,'(A,I3)')
     &        'absorption by CDOM for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

#ifdef GUD_DIAG_RADTRANS_SOLUTION
        WRITE(diagName,'(A,I3.3)')'Estop',j
        WRITE(diagTitle,'(A,I3)')
     &        'diffuse downwelling irradiance for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'WM P    LR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'Eubot',j
        WRITE(diagTitle,'(A,I3)')
     &        'diffuse upwelling irradiance for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'WM P    LR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'amp1_',j
        WRITE(diagTitle,'(A,I3)')
     &        'amplitude of d-w increasing mode for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'amp2_',j
        WRITE(diagTitle,'(A,I3)')
     &        'amplitude of d-w decreasing mode for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'x_',j
        WRITE(diagTitle,'(A,I3)')
     &        'amplitude of direct in diffuse for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'y_',j
        WRITE(diagTitle,'(A,I3)')
     &        'amplitude of direct in upward for waveband ',j
        diagUnits = 'W/m^2           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'r1_',j
        WRITE(diagTitle,'(A,I3)')
     &        'R coefficient of d-w increasing mode for waveband ',j
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'r2_',j
        WRITE(diagTitle,'(A,I3)')
     &        'R coefficient of d-w decreasing mode for waveband ',j
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'att1_',j
        WRITE(diagTitle,'(A,I3)')
     &        'attenuation of d-w increasing mode for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        WRITE(diagName,'(A,I3.3)')'att2_',j
        WRITE(diagTitle,'(A,I3)')
     &        'attenuation of d-w decreasing mode for waveband ',j
        diagUnits = '1/m             '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif
      ENDDO
#endif /* GUD_ALLOW_RADTRANS */

      WRITE(diagName,'(A)')'atten'
      WRITE(diagTitle,'(A)')'attenuation in layer'
      diagUnits = '1               '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'exfPAR'
      WRITE(diagTitle,'(A)')'PAR at surface'
      diagUnits = 'uEin/m^2/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'surfPAR'
      WRITE(diagTitle,'(A)')'PAR at surface'
      diagUnits = 'uEin/m^2/s      '
      diagCode  = 'SM P    L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'PARF'
      WRITE(diagTitle,'(A)')'PAR at top of layer'
      diagUnits = 'uEin/m^2/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'PAR'
      WRITE(diagTitle,'(A,I3)')'total PAR at layer center'
      diagUnits = 'uEin/m^2/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c nutrient consumption
      WRITE(diagName,'(A)')'C_DIN'
      WRITE(diagTitle,'(A)')'consumption of DIN'
      diagUnits = 'mmolN/m3/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'C_PO4'
      WRITE(diagTitle,'(A)')'consumption of PO4'
      diagUnits = 'mmolP/m3/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'C_Si'
      WRITE(diagTitle,'(A)')'consumption of Si'
      diagUnits = 'mmolSi/m3/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'C_Fe'
      WRITE(diagTitle,'(A)')'consumption of Fe'
      diagUnits = 'mmolFe/m3/s      '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

c additional source terms of nutrients
      WRITE(diagName,'(A)')'S_DIN'
      WRITE(diagTitle,'(A)')'non-transport source of DIN'
      diagUnits = 'mmolN/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'S_PO4'
      WRITE(diagTitle,'(A)')'non-transport source of PO4'
      diagUnits = 'mmolP/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'S_Si'
      WRITE(diagTitle,'(A)')'non-transport source of Si'
      diagUnits = 'mmolSi/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'S_Fe'
      WRITE(diagTitle,'(A)')'non-transport source of Fe'
      diagUnits = 'mmolFe/m3/s      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      DO j = 1, nGud
        diagName = 'gGUD'//PTRACERS_ioLabel(j)
C--     Set default name & tracer Units:
        WRITE(locName,'(A,I4.4,A)') 'Tracer ', j
        trcUnits  = 'kg/kg  '
C-      use units & names from data.ptracers :
        ilnb = ILNBLNK(PTRACERS_units(j))
        IF ( ilnb.GE.1 ) THEN
         ilnb = LEN(trcUnits)
         trcUnits = PTRACERS_units(j)(1:ilnb)
        ENDIF
        ilnb = ILNBLNK(PTRACERS_names(j))
        IF ( ilnb.GE.1 ) THEN
         ilnb = MIN(LEN(locName),ilnb)
         WRITE(locName,'(A)') PTRACERS_names(j)(1:ilnb)
        ENDIF
        ilnb = MAX(ILNBLNK(locName),1)
        WRITE(diagTitle,'(2A)') locName(1:ilnb), ' tendency from GUD'
        diagUnits = DIAGS_MK_UNITS( trcUnits//'/s', myThid )
        diagCode  = 'SMR     MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO

      DO j = 1, nPPplank
        WRITE(diagName,'(A,I4.4)')'PP',j
        WRITE(diagTitle,'(A,I4)')'Primary Production plankton ',j
        diagUnits = 'mmol C/m^3/s    '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO

      DO j = 1, nGRplank
        WRITE(diagName,'(A,I4.4)')'GR',j
        WRITE(diagTitle,'(A,I4)')'Grazing loss of plankton ',j
        diagUnits = 'mmol C/m^3/s    '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I             diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      ENDDO

#ifdef GUD_ALLOW_CARBON
      WRITE(diagName,'(A)')'pH'
      WRITE(diagTitle,'(A)')'pH'
      diagUnits = '1               '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'pCO2'
      WRITE(diagTitle,'(A)')'Partial Pressure of CO2'
      diagUnits = 'atm             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'fluxCO2'
      WRITE(diagTitle,'(A)')'Flux of CO2 - air-sea exch'
      diagUnits = 'mmol C/m2/s     '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'gDICsurf'
      WRITE(diagTitle,'(A)')
     &            'Tendency of DIC due to air-sea exch + oldvirt.flux'
      diagUnits = 'mmol C/m3/s     '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'gALKsurf'
      WRITE(diagTitle,'(A)')'Tendency of ALK due to oldvirt.flux'
      diagUnits = 'mmol eq/m3/s    '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'gO2surf'
      WRITE(diagTitle,'(A)')'Tendency of O2 due to air-sea exch'
      diagUnits = 'mmol O/m3/s     '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* GUD_ALLOW_CARBON */

      WRITE(diagName,'(A)')'UTK_holl'
      WRITE(diagTitle,'(A)')'holling'
      diagUnits = '?             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'UTK_2ZPP'
      WRITE(diagTitle,'(A)')'holling'
      diagUnits = '?             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'UTK_C'
      WRITE(diagTitle,'(A)')'C'
      diagUnits = '?             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      WRITE(diagName,'(A)')'UTK_D'
      WRITE(diagTitle,'(A)')'D'
      diagUnits = '?             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        
      WRITE(diagName,'(A)')'UTK_E'
      WRITE(diagTitle,'(A)')'E'
      diagUnits = '?             '
      diagCode  = 'SM P    MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I           diagName, diagCode, diagUnits, diagTitle, 0, myThid )


#endif /* ALLOW_DIAGNOSTICS */
#endif /* ALLOW_GUD */

      RETURN
      END


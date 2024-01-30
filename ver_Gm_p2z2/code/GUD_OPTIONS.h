#ifndef GUD_OPTIONS_H
#define GUD_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_GUD

#include "CPP_OPTIONS.h"

CBOP
C    !ROUTINE: GUD_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C options for gud package
CEOP

C tracer selection

#undef  GUD_ALLOW_NQUOTA
#undef  GUD_ALLOW_PQUOTA
#undef  GUD_ALLOW_FEQUOTA
#undef  GUD_ALLOW_SIQUOTA
#define GUD_ALLOW_CHLQUOTA
#undef  GUD_ALLOW_CDOM
#define GUD_ALLOW_CARBON

C optional bits

#define GUD_ALLOW_DENIT
#undef  GUD_ALLOW_EXUDE
#define ALLOW_OLD_VIRTUALFLUX

C light

#define GUD_READ_PAR
#undef  GUD_USE_QSW
#undef  GUD_AVPAR
#define GUD_ALLOW_GEIDER
#undef  GUD_ALLOW_RADTRANS

C initialize chl with radtrans as in darwin2
#undef  GUD_CHL_INIT_LEGACY

#define GUD_GEIDER_RHO_SYNTH

C grazing

C for quadratic grazing a la quota
#undef  GUD_GRAZING_SWITCH
C#define  GUD_GRAZING_SWITCH
C compute palat from size ratios
#undef  GUD_ALLOMETRIC_PALAT

C turn off grazing temperature dependence
#undef  GUD_NOZOOTEMP

#undef  GUD_TIME_GRAZING

C temperature

#undef  GUD_NOTEMP
#define GUD_TEMP_VERSION 2
#undef  GUD_TEMP_RANGE

C iron

#define GUD_MINFE
#define GUD_PART_SCAV
#define GUD_IRON_SED_SOURCE_VARIABLE

C debugging

#undef  GUD_DEBUG

#define GUD_ALLOW_CONS

#define GUD_UNUSED 0

C deprecated

C these are for gud_generate_random
#undef  GUD_RANDOM_TRAITS
#undef  GUD_TWO_SPECIES_SETUP
#undef  GUD_NINE_SPECIES_SETUP
#undef  GUD_ALLOW_DIAZ

#endif /* ALLOW_GUD */
#endif /* GUD_OPTIONS_H */


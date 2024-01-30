Implementation of the grazing switch occures at the code level within GUD_OPTIONS.h

The model was run with base data.gud files.

Comment out or change to the behavior you wish to use at ~ line 47 of GUD_OPTIONS.h

C for quadratic grazing a la quota
C#undef  GUD_GRAZING_SWITCH
#define  GUD_GRAZING_SWITCH


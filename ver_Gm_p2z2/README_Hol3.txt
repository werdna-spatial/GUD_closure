Implementation of the Holling 3 requires a code change in gud_grazing.F and values in Data.gud.  

This code change becomes a factor as the hollexp is changed from 1 to 2 in data.gud.

Change the commented code in gud_grazingF ~ line 211


        tmp = grazemax(jz)*grazTempFunc(jz)*X(jz)*
C  Holling vs Hill 
C This  makes a difference for hollexp != 1
C Darwin Default which is a HIll
C     &    (sumprey**hollexp/(sumprey**hollexp+kgrazesat(jz)**hollexp))*
C Holling USE THIS FOR HOL3
     &    (sumprey**hollexp/(sumprey**hollexp+kgrazesat(jz)))*
C
     &    (1.0 - EXP(-inhib_graz*sumprey))**inhib_graz_exp
C

undef("eosstat_surface_rho")
function eosstat_1000_rho(t,s)
begin
;========================================================================
;=========================================================================

       c0     =    0.0
       c1     =    1.0
       c2     =    2.0
       c3     =    3.0
       c4     =    4.0
       c5     =    5.0
       c8     =    8.0
       c10    =   10.0
       c16    =   16.0
       c1000  = 1000.0
       c10000 =10000.0
       c1p5   =    1.5
       p33    = c1/c3
       p5     = 0.500
       p25    = 0.250
       p125   = 0.125
       p001   = 0.001
       eps    = 1.0e-10
       eps2   = 1.0e-20
       bignum = 1.0e+30

      tmin = -5.0
      tmax = 50.0
      smin = 0.0
      smax = 50.0

;-----------------------------------------------------------------------
;  UNESCO EOS constants and JMcD bulk modulus constants
;-----------------------------------------------------------------------

;------ for density of fresh water (standard UNESCO)

      unt0 =   999.842594
      unt1 =  6.793952e-2
      unt2 = -9.095290e-3
      unt3 =  1.001685e-4
      unt4 = -1.120083e-6
      unt5 =  6.536332e-9

;------ for dependence of surface density on salinity (UNESCO)

      uns1t0 =  0.824493
      uns1t1 = -4.0899e-3
      uns1t2 =  7.6438e-5
      uns1t3 = -8.2467e-7
      uns1t4 =  5.3875e-9
      unsqt0 = -5.72466e-3
      unsqt1 =  1.0227e-4
      unsqt2 = -1.6546e-6
      uns2t0 =  4.8314e-4

;------ from Table A1 of Jackett and McDougall

      bup0s0t0 =  1.965933e+4
      bup0s0t1 =  1.444304e+2
      bup0s0t2 = -1.706103
      bup0s0t3 =  9.648704e-3
      bup0s0t4 = -4.190253e-5

      bup0s1t0 =  5.284855e+1
      bup0s1t1 = -3.101089e-1
      bup0s1t2 =  6.283263e-3
      bup0s1t3 = -5.084188e-5

      bup0sqt0 =  3.886640e-1
      bup0sqt1 =  9.085835e-3
      bup0sqt2 = -4.619924e-4

      bup1s0t0 =  3.186519
      bup1s0t1 =  2.212276e-2
      bup1s0t2 = -2.984642e-4
      bup1s0t3 =  1.956415e-6

      bup1s1t0 =  6.704388e-3
      bup1s1t1 = -1.847318e-4
      bup1s1t2 =  2.059331e-7
      bup1sqt0 =  1.480266e-4

      bup2s0t0 =  2.102898e-4
      bup2s0t1 = -1.202016e-5
      bup2s0t2 =  1.394680e-7
      bup2s1t0 = -2.040237e-6
      bup2s1t1 =  6.128773e-8
      bup2s1t2 =  6.207323e-10

;=================================================
; To prevent problems with garbage on land points or ghost cells
;=================================================

     TQ = where(t.lt.tmax,t,tmax)
     TQ = where(TQ.gt.tmin,TQ,tmin)

     SQ = where(s.lt.smax,s,smax)
     SQ = where(SQ.gt.smin,SQ,smin)

     p   = 100  ;1000 dbar
     p2  = p*p
     SQR = sqrt(SQ)
     T2  = TQ*TQ


;=================================================
; Calculate surface (p=0) values from UNESCO eqns.
;=================================================

    WORK1 = uns1t0 + uns1t1*TQ + \
           (uns1t2 + uns1t3*TQ + uns1t4*T2)*T2
    WORK2 = SQR*(unsqt0 + unsqt1*TQ + unsqt2*T2)

    RHO_S = unt1*TQ + (unt2 + unt3*TQ + (unt4 + unt5*TQ)*T2)*T2 \
                    + (uns2t0*SQ + WORK1 + WORK2)*SQ

    WORK3 = bup0s1t0 + bup0s1t1*TQ +                    \
           (bup0s1t2 + bup0s1t3*TQ)*T2 +                \
            p *(bup1s1t0 + bup1s1t1*TQ + bup1s1t2*T2) + \
            p2*(bup2s1t0 + bup2s1t1*TQ + bup2s1t2*T2)
    WORK4 = SQR*(bup0sqt0 + bup0sqt1*TQ + bup0sqt2*T2 + \
                 bup1sqt0*p)

    BULK_MOD  = bup0s0t0 + bup0s0t1*TQ +                    \
                (bup0s0t2 + bup0s0t3*TQ + bup0s0t4*T2)*T2 + \
                p *(bup1s0t0 + bup1s0t1*TQ +                \
                   (bup1s0t2 + bup1s0t3*TQ)*T2) +           \
                p2*(bup2s0t0 + bup2s0t1*TQ + bup2s0t2*T2) + \
                SQ*(WORK3 + WORK4)

    DENOMK = c1/(BULK_MOD - p)

    rho    = (unt0 + RHO_S)*BULK_MOD*DENOMK

    drhodt =  unt1 + c2*unt2*TQ +                      \
             (c3*unt3 + c4*unt4*TQ + c5*unt5*T2)*T2 + \
             (uns1t1 + c2*uns1t2*TQ +                 \
              (c3*uns1t3 + c4*uns1t4*TQ)*T2 +         \
              (unsqt1 + c2*unsqt2*TQ)*SQR )*SQ

    drhods  = c2*uns2t0*SQ + WORK1 + c1p5*WORK2

 return(rho)

end

## Main Section 2 # F-[--D--]-F *18
#puts "====== This is Main Section 2 ======"
#puts "  Main Section 2 starts from Z2 = 0."
set LD 0.09
set LC 0.55
set LF [expr {$LQ*2+$LD*2+$LC*2}]
set fQ [expr {$LF/4.0/abs(sin($u*$PI/180.0/2.0))}]
set KQ [expr {1.0/$fQ/$LQ}]
set NL 18
set P0 $P0_S2

set iL 1
while {$iL <= $NL} {

  set P $P0
  set SQ [expr {$KQ*$P/$qe*$LQ}]
  Quadrupole -length $LQ -strength $SQ
  #set IQ [expr {$IQ+1}]
  #set Z2 [expr {$Z2+$LQ}]
  #puts "  Quadrupole $IQ, S = $SQ, Z2 = $Z2"
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  Cavity -length $LC -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
  set P0 [expr {$P0+$GC*$LC*cos($PHI)}]
  #set Z2 [expr {$Z2+$LC}]
  
  set PMin $P0
  set PMax [expr {$PMin+$GC*$LQ*cos($PHI)}]
  set P [expr {$PMin+$pp*($PMax-$PMin)}]
  set SQ [expr {-1.0*$KQ*$P/$qe*$LQ}]
  set SQ_S [expr {$SQ/$NS}]
  set i 1
  while {$i <= $NS} {
    Quadrupole -length 0 -strength $SQ_S
    Cavity -length $LQ_S -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
    set i [expr {$i + 1}]
  }
  set P0 $PMax
  #set IQ [expr {$IQ+1}]
  #set Z2 [expr {$Z2+$LQ}]
  #puts "  Quadrupole $IQ, S = $SQ, Z2 = $Z2"

  Cavity -length $LC -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
  set P0 [expr {$P0+$GC*$LC*cos($PHI)}]
  #set Z2 [expr {$Z2+$LC}]
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  set iL [expr {$iL + 1}]
}

set P $P0
set SQ [expr {$KQ*$P/$qe*$LQ}]
Quadrupole -length $LQ -strength $SQ
#set IQ [expr {$IQ+1}]
#set Z2 [expr {$Z2+$LQ}]
#puts "  Quadrupole $IQ, S = $SQ, Z2 = $Z2"

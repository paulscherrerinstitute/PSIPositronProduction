## Main Section 4 # F-[-]-[-]-D-[-]-[-]-F *7
#puts "====== This is Main Section 4 ======"
#puts "  Main Section 4 starts from Z2 = 0."
set LD 0.05
set LC 1.5
set LF [expr {$LQ*2+$LD*8+$LC*4}]
set fQ [expr {$LF/4.0/abs(sin($u*$PI/180.0/2.0))}]
set KQ [expr {1.0/$fQ/$LQ}]
set NL 7
set P0 $P0_S4

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
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  Cavity -length $LC -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
  set P0 [expr {$P0+$GC*$LC*cos($PHI)}]
  #set Z2 [expr {$Z2+$LC}]
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  set P $P0
  set SQ [expr {-1.0*$KQ*$P/$qe*$LQ}]
  Quadrupole -length $LQ -strength $SQ
  #set IQ [expr {$IQ+1}]
  #set Z2 [expr {$Z2+$LQ}]
  #puts "  Quadrupole $IQ, S = $SQ, Z2 = $Z2"
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  Cavity -length $LC -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
  set P0 [expr {$P0+$GC*$LC*cos($PHI)}]
  #set Z2 [expr {$Z2+$LC}]
  
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]
  Drift -length $LD
  #set Z2 [expr {$Z2+$LD}]

  Cavity -length $LC -gradient $GC -phase $PHID -six_dim 1 -frequency $fC
  #set P0 [expr {$P0+$GC*$LC*cos($PHI)}]
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

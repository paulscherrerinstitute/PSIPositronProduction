set LQ 0.4 ;# Length of quadrupole
set qe 1
set PI [expr {acos(-1)}]
set PhA 90 ;# Phase advance in degree, default: 90
set u $PhA ;# Phase advance
set IQ 0 ;# Quadrupole number
set NS 400 ;# Number of quadrupole slices
set LQ_S [expr {$LQ/$NS}]

set NC_S1 8 ;# Number of structures
set NC_S2 18
set NC_S3 28
set NC_S4 28
set NC_S5 36
set NC [expr {$NC_S1+$NC_S2+$NC_S3+$NC_S4+$NC_S5}]
set LC 1.5 ;# Length of structure
set fC 2 ;# Frequency

#set PHID 0 ;# Phase in degree
set PHI [expr {$PHID * $PI / 180.0}]

# Input mean energy
Octave {
  A = load("$input_file");
  P0_S1 = mean(A(:,1));
  Tcl_SetVar("P0_S1", P0_S1);
} 
set P1_Exp 2.86 ;# Expected mean energy

#set P1_MakeUp 0.0125 ;# Energy make up (CLIC380GeV_Nov2020)
#set DEF_GEV 0.0125 ;# Energy make up (CLIC380GeV_Nov2020)

#set GC [expr ($P1_Exp - $P0_S1 + $P1_MakeUp) / ($NC * $LC * cos($PHI)) ]
set GC [expr ($P1_Exp - $P0_S1 + $DEF_GEV) / ($NC * $LC * cos($PHI)) ]

set P0_S2 [expr {$P0_S1+$GC*$NC_S1*$LC*cos($PHI)}]
set P0_S3 [expr {$P0_S2+$GC*$NC_S2*$LC*cos($PHI)}]
set P0_S4 [expr {$P0_S3+$GC*$NC_S3*$LC*cos($PHI)}]
set P0_S5 [expr {$P0_S4+$GC*$NC_S4*$LC*cos($PHI)}]
set P0_S6 [expr {$P0_S5+$GC*$NC_S5*$LC*cos($PHI)}]
#set Z2 0 ;# Exit Z position
set pp 0.00 ;# P0 percentage in quadrupole, [0,1]

#puts "INFO::"
#puts "  P0_S1 = $P0_S1 GeV"
#puts "  P0_S2 = $P0_S2 GeV"
#puts "  P0_S3 = $P0_S3 GeV"
#puts "  P0_S4 = $P0_S4 GeV"
#puts "  P0_S5 = $P0_S5 GeV"
#puts "  P0_S6 = $P0_S6 GeV"


set ch {1.0}

#
# Choice of Structure
#

#array set structure {a 1.9e-3 g 2.325e-2 l 3.05e-3}
#array set structure {a 1.74e-3 g 2.29452380952381e-3 l 3.05952380952381e-3}
array set structure {a 2.0e-3 g 2.325e-3 l 3.05e-3}
array set structure {a 2e-3 g 1.5e-3 l 2.5e-3}
#array set structure {a 1.78e-3 g 1.17e-3 l 1.67e-3}

#
# transverse wakefield
# s is given in micro metres
# return value is in V/pCm^2
#

proc w_transv {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set s0 [expr 0.169*pow($a,1.79)*pow($g,0.38)*pow($l,-1.17)]
    return [expr 4.0*377.0*3e8*$s0*1e-12/(acos(-1.0)*pow($a,4))*(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]
}

proc w_transv {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set tmp [expr $g/$l]
    set alpha [expr 1.0-0.4648*sqrt($tmp)-(1.0-2.0*0.4648)*$tmp]
    set s0 [expr $g/8.0*pow($a/($l*$alpha),2)]
    return [expr 4.0*377.0*3e8*$s0*1e-12/(acos(-1.0)*pow($a,4))*(1.0-(1.0+sqrt($s*1e-6/$s0))*exp(-sqrt($s*1e-6/$s0)))]
}

#
# longitudinal wakefield
# s is given in micro metres
# return value is in V/pCm
#

proc w_long {s} {
    global structure
    set a $structure(a)
    set g $structure(g)
    set l $structure(l)
    set tmp [expr $g/$l]
    set alpha [expr 1.0-0.4648*sqrt($tmp)-(1.0-2.0*0.4648)*$tmp]
    set s0 [expr $g/8.0*pow($a/($l*$alpha),2)]
    return [expr 377.0*3e8*1e-12/(acos(-1.0)*$a*$a)*exp(-sqrt($s*1e-6/$s0))]
}

# remove wakefields

proc w_transv {s} {
 return 0.0
}

proc w_long {s} {
 return 0.0
}

source $common_script_dir/wake_calc.tcl

#
# Define the beams
#

source $common_script_dir/make_beam.tcl

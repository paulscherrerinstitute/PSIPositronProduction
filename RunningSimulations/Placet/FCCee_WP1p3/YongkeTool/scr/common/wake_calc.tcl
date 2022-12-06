#
# generic routine to calculate the wakefields
# it creates a file with name f_name that can be read when a beam is created
# z_list is of longitudinal position and weight of the slices
#

proc calc_new_transv {charge z_list} {
    set tmp {}
#
# unit to store in the file is MV/pCm^2
#

#
# multiply charge with factor for pC and MV
#
    set z_l {}
    foreach l $z_list {
	lappend z_l $l
	set z0 [lindex $l 0]
	foreach j $z_l {
	    set z [lindex $j 0]
	    lappend tmp [expr 1.6e-7*1e-6*$charge*[w_transv [expr $z0-$z]]]
	}
    }
    return $tmp
}

proc calc_new_long {charge z_list} {
    set z_l {}
    set n [llength $z_list]
    set tmp {}
    foreach l $z_list {
	set z0 [lindex $l 0]
	set wgt0 [expr $charge*[lindex $l 1]]
	set sum 0.0
	foreach j $z_l {
	    set z [lindex $j 0]
	    set wgt [expr [lindex $j 1]*$charge]
	    set sum [expr $sum+$wgt*[w_long [expr $z0-$z]]]
	}
	set sum [expr $sum+0.5*$wgt0*[w_long 0.0]]
	lappend z_l $l
#
# multiply with factor for pC and MV
#
	set sum [expr -($sum*1.6e-7*1e-6)]
	lappend tmp "$z0 $sum"
    }
    return $tmp
}

proc calc {f_name charge a b sigma n} {
    set file [open $f_name w]
    puts $file $n
    set long [calc_new_long $charge [GaussList -min $a -max $b -sigma $sigma -charge 1.0 -n_slices [expr 5*$n]]]
    set z_list [GaussList -min $a -max $b -sigma $sigma -charge 1.0 -n_slices [expr $n]]
    set transv [calc_new_transv $charge $z_list]
    set x [lindex $long [expr int($n/2)*5+2]]
    puts $file [lindex $x 1]
    for {set i 0} {$i<$n} {incr i} {
	set x [lindex $long [expr 2+5*$i]]
	puts $file "$x [lindex [lindex $z_list $i] 1]"
    }
    foreach x $transv {
	puts $file $x
    }
    close $file
}

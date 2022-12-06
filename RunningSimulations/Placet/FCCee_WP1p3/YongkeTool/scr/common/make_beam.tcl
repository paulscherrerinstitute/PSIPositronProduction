proc make_beam_particles {e0 e_spread n} {
    global match
    set sum_e 0.0
    set f [open particles.tmp w]

    Random g -type gaussian
    Random g2 -type linear

    set ex [expr $match(emitt_x)*1e-7*0.511e-3/$e0]
    set ey [expr $match(emitt_y)*1e-7*0.511e-3/$e0]
    set bx $match(beta_x)
    set by $match(beta_y)
    set alphax $match(alpha_x)
    set alphay $match(alpha_y)
    set sx [expr sqrt($ex*$bx)*1e6]
    set spx [expr sqrt($ex/$bx)*1e6]
    set sy [expr sqrt($ey*$by)*1e6]
    set spy [expr sqrt($ey/$by)*1e6]
    set st $match(sigma_z)
    set se [expr 0.01*abs($e_spread)]

    if ($e_spread<0.0) {
	for {set i 0} {$i<$n} {incr i} {
	    set e [expr (1.0+$se*([g2]-0.5))*$e0]
	    set z [g]
	    while {abs($z)>=3.0} {
		set z [g]
	    }
	    set x [expr [g]*$sx]
	    set xp [expr [g]*$spx-$alphax*$x/$sx*$spx]
	    set y [expr [g]*$sy]
	    set yp [expr [g]*$spy-$alphay*$y/$sy*$spy]
	    puts $f "$e $x $y [expr $z*$st] $xp $yp"
	}
    } {
	for {set i 0} {$i<$n} {incr i} {
	    set e [g]
	    while {abs($e)>=3.0} {
		set e [g]
	    }
	    set e [expr (1.0+$e*$se)*$e0]
	    set z [g]
	    while {abs($z)>=3.0} {
		set z [g]
	    }
	    set x [expr [g]*$sx]
	    set xp [expr [g]*$spx-$alphax*$x/$sx*$spx]
	    set y [expr [g]*$sy]
	    set yp [expr [g]*$spy-$alphay*$y/$sy*$spy]
#	    set y [expr $y+$sy]
	    puts $f "$e $x $y [expr $z*$st] $xp $yp"
	}
    }
    close $f
    exec sort -b -g -k 4,4 particles.tmp > particles.in
    exec rm -f particles.tmp
}

proc make_beam_many {name nslice n} {
    global charge e_initial match n_total
    set ch {1.0}
#    set n_total 50000
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -particles $n_total \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l {1.0 0.0 0.0}
    SetRfGradientSingle $name 0 $l    
    make_beam_particles $e_initial $match(e_spread) [expr $nslice*$n]
    BeamRead -file particles.in -beam $name
    #exec rm -f particles.in
    #exec rm -f beam.dat
}

proc make_beam_many_gradient {name nslice n de} {
    global charge e_initial match
    set ch {1.0}
#    calc_old beam.dat [GaussList -charge $charge -min -3 -max 3 -sigma $match(sigma_z) -n_slices $nslice]
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread 0.0 \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "[expr 1.0+0.01*$de] 0.0 0.0"
    SetRfGradientSingle $name 0 $l    
    make_beam_particles $e_initial $match(e_spread) [expr $nslice*$n]
    BeamRead -file particles.in -beam $name
}

proc make_beam_many_phase {name nslice n dph} {
    global charge e_initial match
    set ch {1.0}
#    calc_old beam.dat [GaussList -charge $charge -min -3 -max 3 -sigma $match(sigma_z) -n_slices $nslice]
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread 0.0 \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase $dph \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "[expr 1.0] 0.0 0.0"
    SetRfGradientSingle $name 0 $l    
    make_beam_particles $e_initial $match(e_spread) [expr $nslice*$n]
    BeamRead -file particles.in -beam $name
}

proc make_beam_slice {name nslice n} {
    global charge e_initial match n_total
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)\
	    -particles $n_total

    set l {}
    lappend l {1.0 0.0 0.0}
    SetRfGradientSingle $name 0 $l
}

proc make_beam_train {name nb nslice n} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches $nb \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge $charge \
	    -phase 0.0 \
	    -last_wgt [expr 155-$nb] \
	    -overlapp [expr -0.08] \
	    -distance [expr 0.08] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l {1.0 0.0 0.0}
    SetRfGradientSingle $name 0 $l    
}

proc make_beam_train_main {name nb nslice n} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
#    set wl {}
#    for {set i 0} {$i<$nb} {incr i} {
#	lappend wl [wakelong [expr 0.2*$i]]
#    }
    set wl 0.0
    if {$nb>=2} {
	lappend wl [wakelong 0.08]
    }
    for {set i 2} {$i<$nb} {incr i} {
	lappend wl 0.0
    }
    MainBeam $name -bunches $nb \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -longrange_list $wl \
	    -e0 $e_initial \
	    -file beam.dat \
	    -charge $charge \
	    -phase 0.0 \
	    -last_wgt [expr 155-$nb] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l {1.0 0.0 0.0}
    SetRfGradientSingle $name 0 $l    
}

proc make_beam_slice_gradient {name nslice n grad} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "$grad 0.0 0.0"
    SetRfGradientSingle $name 0 $l  
}

proc make_beam_slice_energy {name nslice n eng} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 [expr $eng*$e_initial] \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "1.0 0.0 0.0"
    SetRfGradientSingle $name 0 $l  
}

proc make_beam_slice_energy_gradient {name nslice n eng grad} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3 3 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.0 \
	    -e0 [expr $eng*$e_initial] \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "$grad 0.0 0.0"
    SetRfGradientSingle $name 0 $l  
}

proc make_beam_slice35 {name nslice n} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3.5 3.5 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.5 \
	    -e0 $e_initial \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l {1.0 0.0 0.0}
    SetRfGradientSingle $name 0 $l    
}

proc make_beam_slice35_energy_gradient {name nslice n eng grad} {
    global charge e_initial match
    set ch {1.0}
    calc beam.dat $charge -3.5 3.5 $match(sigma_z) $nslice
    InjectorBeam $name -bunches 1 \
	    -macroparticles [expr $n] \
	    -energyspread [expr 0.01*$match(e_spread)*$e_initial] \
	    -ecut 3.5 \
	    -e0 [expr $eng*$e_initial] \
	    -file beam.dat \
	    -chargelist $ch \
	    -charge 1.0 \
	    -phase 0.0 \
	    -overlapp [expr -390*0.3/1.3] \
	    -distance [expr 0.3/1.3] \
	    -alpha_y $match(alpha_y) \
	    -beta_y $match(beta_y) \
	    -emitt_y $match(emitt_y) \
	    -alpha_x $match(alpha_x) \
	    -beta_x $match(beta_x) \
	    -emitt_x $match(emitt_x)

    set l {}
    lappend l "$grad 0.0 0.0"
    SetRfGradientSingle $name 0 $l    
}

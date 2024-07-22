# Wang-Landau algorithm v1.1 beta
# by Haohao Fu fhh2626_at_gmail.com
#
# Usage:
# (in NAMD config file)
#    set st_outputPrefix st
#    set st_outputFreq 100
#    set st_restartFreq 20000
#    set st_exchangeFreq 2000
#    set st_loop 20000
#    set st_langevinPiston   0
#    
#    # for Wang-Landau simulated tempering:
#    set st_temp {300 320 340 360 380 400 420 440 460 480 500}
#
#    # for Wang-Landau simulated solute scaling
#    set st_factors {1.0 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5}
#    set st_temp 300
#    set st_indicationFile st.pdb
#    set st_indicationCol B
#
#    # for Wang-Landau Lambding
#    set st_factors {0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0}
#    set st_temp 300
#    alchFile   st.pdb
#    alchCol    B
#    alchOutFreq   1000
#    alchOutFile   st_alch.out
#    alchVdwShiftCoeff   5
#    alchElecLambdaStart   0.5
#    alchVdwLambdaEnd      0.7
#    alchBondDecouple  off
#
# IMPORTANT:
#    
#    one must let outputEnergies == st_outputFreq in the NAMD config file!
#    otherwise NAMD and the call_back function will become very buggy!
#
#
# Optional:
#    # type of Wang-Landau method (0 or 1, default 0)
#    # 0: classic simulated tempering
#    # 1: simulated solute scaling
#    # 2: lambding
#    set st_type 0
#    
#    # set restart file name (0 for a new simulation, otherwise the restart file name, default 0)
#    set st_restart 0
#
#    # whether Colvars module is activated
#    # this feature is very flexible, as the Colvars input can be different in different temperature
#    set st_colvars 0
#    
#    # the name of colvars config files for different temperature 
#    #  should be colvars_300.in, colvars_320.in, etc.
#    set st_colvarsPrefix colvars
#
#    # the weight of first temperature
#    # if this is 2, then the probability of the first temperature is twice as that of others
#    set st_firstWeight 1
#
#    # if one wants to run a classical simulated tempering/ simulated solute scaling/ simulated lambding
#    set st_weight {0 10 30 60 100 200 400 800 1600 2400 3200}
#;#########################################################################

if {$st_exchangeFreq % $st_outputFreq != 0} {
    puts "Error! st_exchangeFreq must be a multiple of st_outputFreq!"
    exit 1
}
if {$st_restartFreq % $st_exchangeFreq != 0} {
    puts "Error! st_restartFreq must be a multiple of st_exchangeFreq!"
    exit 1
}

set st_Boltzmann 0.001987191
set st_weightedFactor 1.0

# type of simulation
if {[info exist st_type] == 0} {
    set st_type 0
}

# whether Colvars is activated
if {[info exist st_colvars] == 0} {
    set st_colvars 0
}

# the weight of the first temperature
if {[info exist st_firstWeight] == 0} {
    set st_firstWeight 1
}

# initialize solute scaling
if {$st_type == 1} {
    soluteScaling  on
    soluteScalingFactor  1.0
    soluteScalingFile   $st_indicationFile
    soluteScalingCol    $st_indicationCol
}

# initialize lambding
if {$st_type == 2} {
    alch on
    alchType fep
    alchLambda   0.0
    alchLambda2  0.0
    alchEquilSteps  0
}

# beta for simulated tempering
if {$st_type == 0} {
    set st_beta {}
}
if {$st_type == 1 || $st_type == 2} {
    set st_beta [expr 1 / ($st_Boltzmann * $st_temp)]
}
# weight of simulated tempering
# auto updated by the Wang-Landau method
if {[info exist st_weight] == 0} {
    set st_weight {}
    set st_updateWeight 1
} else {
    if {$st_type == 0} {
        if {[llength $st_weight] != [llength $st_temp]} {
            puts "Error! The length of st_weight must equal to st_temp for a classical simulated tempering!"
            exit 1
        }
    } 
    if {$st_type == 1} {
        if {[llength $st_weight] != [llength $st_factors]} {
            puts "Error! The length of st_weight must equal to st_factors for a classical simulated solute scaling!"
            exit 1
        }
    }
    if {$st_type == 2} {
        if {[llength $st_weight] != [llength $st_lambdas]} {
            puts "Error! The length of st_weight must equal to st_lambdas for a classical simulated solute scaling!"
            exit 1
        }
    }
    set st_updateWeight 0
}
# histogram for temperatures or factors
set st_histogram {}
# current temperature
if {$st_type == 0} {
    set st_currentTemp [lindex $st_temp 0]
    if {$st_colvars == 1} {
        cv reset
        cv configfile ${st_colvarsPrefix}_${st_currentTemp}.in
        cv targettemperature $st_currentTemp
    }
}
# current factor
if {$st_type == 1 || $st_type == 2} {
    set st_currentFactor [lindex $st_factors 0]
    if {$st_colvars == 1} {
        cv reset
        cv configfile ${st_colvarsPrefix}_${st_currentFactor}.in
    }
}
# whether changed temperature or factors once
set st_changed 0
# probabilities of changing temperatures or factors
set st_probability {}
# number of temperatures or factors
if {$st_type == 0} {
    set st_numTemp [llength $st_temp]
}
if {$st_type == 1 || $st_type == 2} {
    set st_numTemp [llength $st_factors]
}

# colvars biases
if {$st_colvars == 1} {
    set st_colvarsBiases [cv list biases]
}


# energy for simulated tempering
# and current energy for solute scaling
set st_energy 0
# energies for each factor in solute scaling
set st_energies {}

if {[info exist st_restart] == 0} {
    set st_restart 0
}

#; split a line when more than one space between two string
#; stupid tcl!
proc st_splitline {line} {
	set line2 [string trim $line]
	regsub -all {[[:blank:]]+} $line2 " " line3
	return [split $line3]
}

# initialization
if {$st_type == 1 || $st_type == 2} {
    for {set i 0} {$i < $st_numTemp} {incr i} {
        lappend st_energies 0
    }
}
if {$st_restart == 0} {

    # step
    set st_step 0
    
    for {set i 0} {$i < $st_numTemp} {incr i} {
        if {$st_type == 0} {
            lappend st_beta [expr 1 / ($st_Boltzmann * [lindex $st_temp $i])]
        }
        if {$st_updateWeight == 1} {
            lappend st_weight 0
        }
        lappend st_histogram 0
        lappend st_probability 0
    }
} else {
    # restart simulation from a restart file
    set fp [open $st_restart r]
    gets $fp temp_data
    set st_step $temp_data
    gets $fp temp_data
    
    if {$st_type == 0} {
        set st_currentTemp $temp_data
    }
    if {$st_type == 1 || $st_type == 2} {
        set st_currentFactor $temp_data
    }
    
    gets $fp temp_data
    set st_weightedFactor $temp_data
    gets $fp temp_data
    set st_changed $temp_data
    gets $fp temp_data
    set splitedFactors [st_splitline $temp_data]
    gets $fp temp_data
    set splitedWeight [st_splitline $temp_data]
    gets $fp temp_data
    set splitedHist [st_splitline $temp_data]
    gets $fp temp_data
    set splitedprobabilities [st_splitline $temp_data]
    set st_numTemp [llength $splitedFactors]
    
    if {$st_type == 1 || $st_type == 2} {
        set st_factors {}
    }
    
    for {set i 0} {$i < $st_numTemp} {incr i} {
        if {$st_type == 0} {
            lappend st_beta [lindex $splitedFactors $i]
        }
        if {$st_type == 1 || $st_type == 2} {
            lappend st_factors [lindex $splitedFactors $i]
        }
        lappend st_weight [lindex $splitedWeight $i]
        lappend st_histogram [lindex $splitedHist $i]
        lappend st_probability [lindex $splitedprobabilities $i]
    }
    close $fp
}

# return the min value of a list
proc st_min {inputList} {
    set min [lindex $inputList 0]
    for {set i 0} {$i < [llength $inputList]} {incr i} {
        if {[lindex $inputList $i] < $min} {
            set min [lindex $inputList $i]
        }
    }
    return $min
}

# return the max value of a list
proc st_max {inputList} {
    set max [lindex $inputList 0]
    for {set i 0} {$i < [llength $inputList]} {incr i} {
        if {[lindex $inputList $i] > $max} {
            set max [lindex $inputList $i]
        }
    }
    return $max
}

# return the sum value of a list
proc st_sum {inputList} {
    set sum 0
    for {set i 0} {$i < [llength $inputList]} {incr i} {
        set sum [expr $sum + [lindex $inputList $i]]
    }
    return $sum
}

# write restart file
proc st_writeRestart {} {
    global st_outputPrefix st_currentTemp st_weightedFactor st_beta st_weight st_factors
    global st_histogram st_probability st_changed st_numTemp st_step st_type st_currentFactor
    set fp [open "${st_outputPrefix}.st.restart" w]
    puts $fp $st_step
    
    if {$st_type == 0} {
        puts $fp $st_currentTemp
    }
    if {$st_type == 1 || $st_type == 2} {
        puts $fp $st_currentFactor
    }
    
    puts $fp $st_weightedFactor
    puts $fp $st_changed
    set factors ""
    set weight ""
    set hist ""
    set probabilities ""
    for {set i 0} {$i < $st_numTemp} {incr i} {
        if {$st_type == 0} {
            append factors " " [lindex $st_beta $i]
        }
        if {$st_type == 1 || $st_type == 2} {
            append factors " " [lindex $st_factors $i]
        }
        append weight " " [lindex $st_weight $i]
        append hist " " [lindex $st_histogram $i]
        append probabilities " " [lindex $st_probability $i]
    }
    puts $fp $factors
    puts $fp $weight
    puts $fp $hist
    puts $fp $probabilities
    close $fp
}

# get potential of this structure
# write outputs 
proc st_callback {labels values} {
    global st_energy st_outputFreq st_step
    # the 10th term is potential
    set value [st_splitline $values]
    set st_energy [lindex $value 10]
}
callback st_callback


# get all the energies corresponding to each lambda
proc st_getAllEnergy {} {
    global st_type st_numTemp st_factors st_energy st_energies
    
    if {$st_type == 1} {
        for {set i 0} {$i < $st_numTemp} {incr i} {
            soluteScalingFactor [lindex $st_factors $i]
            run 0
        
            set st_energies [lreplace $st_energies $i $i $st_energy]
        }
    }
    
    if {$st_type == 2} {
        for {set i 0} {$i < $st_numTemp} {incr i} {
            alchLambda [lindex $st_factors $i]
            run 0
        
            set st_energies [lreplace $st_energies $i $i $st_energy]
        }
    }
}

# update the probabilities of changing temperatures/scaling factors
proc st_updateProbabilities {potential} {
    global st_weight st_probability st_numTemp st_beta st_energies st_type
    set lnPN {}
    set PN {}
    set offset 0
    
    # Compute the probability for each temperature.  This is done in log space to avoid overflow.
    if {$st_type == 0} {
        for {set i 0} {$i < $st_numTemp} {incr i} {
            lappend lnPN [expr [lindex $st_weight $i] - $potential * [lindex $st_beta $i]]
        }
    }
    if {$st_type == 1 || $st_type == 2} {
        for {set i 0} {$i < $st_numTemp} {incr i} {
            lappend lnPN [expr [lindex $st_weight $i] - [lindex $st_energies $i] * $st_beta]
        }
    }
    
    set maxLogProb [st_max $lnPN]
        for {set i 0} {$i < $st_numTemp} {incr i} {
            set offset [expr $offset + exp([lindex $lnPN $i] - $maxLogProb)]
        }
        set offset [expr $maxLogProb + log($offset)]
        for {set i 0} {$i < $st_numTemp} {incr i} {
            lappend PN [expr exp([lindex $lnPN $i] - $offset)]
        }
    
    set st_probability $PN
}

# change temperatures based on $st_probability
# also update weights and histograms
proc st_changeTemperatures {} {
    global st_weight st_updateWeight st_probability st_numTemp st_langevinPiston st_temp 
    global st_currentTemp st_weightedFactor st_histogram st_changeTemperatures st_changed
    global st_type st_factors st_currentFactor st_colvarsPrefix st_colvars st_colvarsBiases
    global st_firstWeight
    
    set r [expr rand()]
    for {set i 0} {$i < $st_numTemp} {incr i} {
        if {$r < [lindex $st_probability $i]} {
            # for simulated tempering, change temperature
            if {$st_type == 0} {
                if {[lindex $st_temp $i] != $st_currentTemp} {
                    # change temperature
                    langevinTemp [lindex $st_temp $i]

                    # if langevin piston is used to control pressure
                    if {$st_langevinPiston != 0} {
                        langevinPistonTemp [lindex $st_temp $i]
                    }

                    rescalevels [expr sqrt([lindex $st_temp $i] / $st_currentTemp)]
                    print "ST: Changing temperature from $st_currentTemp K to [lindex $st_temp $i] K!"
                    set st_currentTemp [lindex $st_temp $i]
                    
                    # if Colvars is used
                    if {$st_colvars == 1} {
                        cv reset
                        cv configfile ${st_colvarsPrefix}_${st_currentTemp}.in
                        
                        foreach bias $st_colvarsBiases {
                            if {[file exists ${st_colvarsPrefix}_${st_currentTemp}_${bias}.colvars.state]} {
                                cv bias $bias load ${st_colvarsPrefix}_${st_currentTemp}_${bias}
                            }
                        }
                        cv targettemperature $st_currentTemp
                    }
                    
                    set st_changed 1
                }
            }
            # for solute scaling, change scaling factor
            if {$st_type == 1} {
                if {[lindex $st_factors $i] != $st_currentFactor} {
                    # change lambda
                    soluteScalingFactor [lindex $st_factors $i]

                    print "ST: Changing soluteScalingFactor from $st_currentFactor to [lindex $st_factors $i] !"
                    set st_currentFactor [lindex $st_factors $i]
                    
                    # if Colvars is used
                    if {$st_colvars == 1} {
                        cv reset
                        cv configfile ${st_colvarsPrefix}_${st_currentFactor}.in
                        
                        foreach bias $st_colvarsBiases {
                            if {[file exists ${st_colvarsPrefix}_${st_currentFactor}_${bias}.colvars.state]} {
                                cv bias $bias load ${st_colvarsPrefix}_${st_currentFactor}_${bias}
                            }
                        }
                    }
                    
                    set st_changed 1
                }
            }
            # for lambding, change alchLambda
            if {$st_type == 2} {
                if {[lindex $st_factors $i] != $st_currentFactor} {
                    # change lambda
                    alchLambda [lindex $st_factors $i]
                    print "ST: Changing alchLambda from $st_currentFactor to [lindex $st_factors $i] !"
                    set st_currentFactor [lindex $st_factors $i]
                    set st_changed 1
                }
            }
            
            # update weights
            if {$st_updateWeight == 1} {
                
                # make sampling of the first replica more
                if {$i == 0} {
                        set st_weight [lreplace $st_weight $i $i [expr [lindex $st_weight $i] - ($st_weightedFactor / $st_firstWeight)]]
                } else {
                        set st_weight [lreplace $st_weight $i $i [expr [lindex $st_weight $i] - $st_weightedFactor]]
                }
                set st_histogram [lreplace $st_histogram $i $i [expr [lindex $st_histogram $i] + 1]]
                
                set st_fixedHistogram $st_histogram
                set st_fixedHistogram [lreplace $st_fixedHistogram 0 0 [expr [lindex $st_fixedHistogram 0] / $st_firstWeight]]
                
                set minCount [st_min $st_fixedHistogram]
                if {$minCount > 20 && $minCount >= [expr 0.2 * [st_sum $st_fixedHistogram] / $st_numTemp]} {
                    set st_weightedFactor [expr $st_weightedFactor * 0.5]
                    set weight0 [lindex $st_weight 0]
                    set st_histogram {}
                    for {set j 0} {$j < $st_numTemp} {incr j} {
                        lappend st_histogram 0
                        set st_weight [lreplace $st_weight $j $j [expr [lindex $st_weight $j] - $weight0]]
                    }  
                } elseif {[lindex $st_probability $i] > 0.99 && $st_changed == 0} {
                    set st_weightedFactor [expr $st_weightedFactor * 2]
                    set st_histogram {}
                    for {set j 0} {$j < $st_numTemp} {incr j} {
                        lappend st_histogram 0
                    }
                }
            }
            
            return
        } else {
            set r [expr $r - [lindex $st_probability $i]]
        }
    }
}

# the result string
proc st_resultString {} {
    global st_energy st_currentTemp st_numTemp st_weight st_step st_outputFreq st_type st_currentFactor
    set st_step [expr $st_step + $st_outputFreq]
    if {$st_type == 0} {
        set str "ST: Step: $st_step Temp: $st_currentTemp Energy: $st_energy Weight:"
    }
    if {$st_type == 1 || $st_type == 2} {
        set str "ST: Step: $st_step Factor: $st_currentFactor Energy: $st_energy Weight:"
    }
    for {set i 0} {$i < $st_numTemp} {incr i} {
        append str " " [lindex $st_weight $i]
    }
    return $str
}
        
# main loop
set innerLoopTime [expr $st_exchangeFreq / $st_outputFreq]
for {set loop 0} {$loop < $st_loop} {incr loop} {
    for {set innerLoop 0} {$innerLoop < $innerLoopTime} {incr innerLoop} {
        run $st_outputFreq
        print [st_resultString]
    }
    if {$st_type == 1 || $st_type == 2} {
        st_getAllEnergy
    }
    
    if {$st_colvars == 1} {
        if {$st_type == 0} {
            foreach bias $st_colvarsBiases {
                cv bias $bias save ${st_colvarsPrefix}_${st_currentTemp}_${bias}
            }
        }
        if {$st_type == 1} {
            foreach bias $st_colvarsBiases {
                cv bias $bias save ${st_colvarsPrefix}_${st_currentFactor}_${bias}
            }
        }
    }
    
    st_updateProbabilities $st_energy
    st_changeTemperatures
    # only write restart file when running the Wang-Landau algorithm
    if {$st_step % $st_restartFreq == 0 && $st_updateWeight == 1} {
        st_writeRestart
    }
}
"""
    Data type for cells, machines and states.
"""
struct MachineState end

"""
    Data type for cells, machines and positions on a cell.
"""
struct MachineData end

"""
    Struct for the so-called machinestate with field names
        time
        location_index
        cell_id
        cell_type_label
        number_machines_state_one
        number_machines_state_two
        number_machines_state_three
        number_neighbours_different_label
"""
struct MachineStateCellProperties 
    time
    location_index
    cell_id
    cell_type_label
    number_machines_state_one
    number_machines_state_two
    number_machines_state_three
    number_neighbours_different_label
end


"""
    Cell age data writer
    Params:

    time
    location_index
    cell_centre_coord_x  
    cell_centre_coord_y  
    cell_centre_coord_z  
    cell_age
"""
struct CellAgeProperties 
    time
    location_index
    cell_centre_coord_x  
    cell_centre_coord_y  
    cell_centre_coord_z  
    cell_age
end 

struct CellAgeData end  
"""
    Struct for the so-called machinedata
        time             
        cell_id           
        cell_label       
        machine_id        
        machine_state     
        machine_coord_x  
        machine_coord_y  
        machine_coord_z  
"""
struct MachineDataProperties
    time             
    cell_id           
    cell_label       
    machine_id        
    machine_state     
    machine_coord_x  
    machine_coord_y  
    machine_coord_z  
end

"""
    Get time series and compeition between cell_labels with field names
        time
        label₀
        label₁
        competition
"""
struct MachineStateTimeSeries
    time
    target
    attacker
end

"""
    Struct for population counts per unique time with field names
        time_to_match
        count_label_0
        count_label_1
        count_total
        count_machines_state_one
        count_machines_state_two
        count_machines_state_three
        count_machines_total
        count_machine_state_one_to_pop_01
        count_machine_state_two_to_pop_01
        count_machine_state_three_to_pop_01
        count_machine_total_to_pop_01
        count_neighbours_diff_to_cell_0
        count_neighbours_diff_to_cell_1
"""
struct PopulationCountsPerTimeStep
    time_to_match
    count_target
    count_attacker
    count_total
    count_machines_state_one
    count_machines_state_two
    count_machines_state_three
    count_machines_total
    count_machine_state_one_to_attacker
    count_machine_state_two_to_attacker
    count_machine_state_three_to_attacker
    count_machine_total_to_attacker
    count_neighbours_diff_to_target
    count_neighbours_diff_to_attacker
end


"""
    Struct for different times with field names
        time_vec
        unique_time_vec
        shifted_time_vec
        unique_shifted_time_vec
"""
struct TimeVecVariants 
    time_vec
    unique_time_vec
    shifted_time_vec
    unique_shifted_time_vec
end


"""
    Time from start of apoptosis to death of a cell
"""
struct ApoptosisTime
    tₐ::Float64
end


"""
    The radius of a circle
"""
struct Radius 
    r::Float64
end


"""
    The length of a cylinder
"""
struct Cylinder
    length
end

"""
 The structure of a capsule 
"""
struct Capsule
    hemisphere_radius::Radius
    cylinder_length::Cylinder
 end


 struct CellLabelsTimeSeries
    neighbour_to_target
    neighbour_to_attacker
end

struct CellTimeSeries
    MachineStateTimeSeries
    CellLabelsTimeSeries
end


struct TimeForDeathRate
    t::Float64
end

struct DeathRate
    d::Float64
end

struct DeathRateTime
    d::Float64
end

struct N₀AlreadySet
    n₀
end


struct K₁   
    k₁
end


struct T6SSKineticsPositionTime
    t
end

struct CellDivisionTime
    t
end

struct Local end
struct Global end
struct DeathRateMichaelisMentenForm end
struct CellBoundaryNormal end
struct DataAllCellLabels end



# cell scaling
# time, location_index, cell_id, 
struct GenericDatToMat end
struct CellScalingProperties
    time
    location_index             
    cell_id               
    cell_centre_coord_x  
    cell_centre_coord_y  
    cell_centre_coord_z 
    cell_radius
    cell_length 
    unknown_value
end

struct CellOrientationProperties
    time
    location_index             
    cell_id               
    cell_centre_coord_x  
    cell_centre_coord_y  
    cell_centre_coord_z 
    cell_orient_direction_x  
    cell_orient_direction_y  
    cell_orient_direction_z 
end
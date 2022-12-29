"""
    Per each time step extract the cell label for each cell and count. 

    The function returns the count for each time per each time step.
"""
function get_type_count_per_time(::MachineState,cell_data,time_to_match)
    cell_type_vector = [0,1]
    cells_matching_time = filter(x -> x.time == time_to_match,cell_data)
    count_each_type = map(y -> count(x -> x.cell_type_label == y,cells_matching_time),cell_type_vector)
    return count_each_type
end


"""
    This function gets the population counts per each cell label and the counts of machines. Returns the following time series: 
    1. `time_to_match`
    2. `count_target`
    3. `count_attacker`
    4. `count_total`
    5. `count_machines_state_one`
    6. `count_machines_state_two`
    7. `count_machines_state_three`
    8. `count_machines_total`
    9. `count_machine_state_one_to_attacker`
    10. `count_machine_state_two_to_attacker`
    11. `count_machine_state_three_to_attacker`
    12. `count_machine_total_to_attacker`
    13. `count_neighbours_diff_to_target`
    14. `count_neighbours_diff_to_attacker`
"""
function get_population_counts_per_time_step(cell_data,time_to_match)

    cells_matching_time = filter(x -> x.time == time_to_match,cell_data)
    subset_cell_target = filter(x->x.cell_type_label == 0,cells_matching_time)
    subset_cell_attacker = filter(x->x.cell_type_label == 1,cells_matching_time)


    count_target = count(x -> x.cell_type_label == 0,cells_matching_time)
    count_attacker = count(x -> x.cell_type_label == 1,cells_matching_time)
    count_total = count_target + count_attacker


    count_machines_state_one = sum([c.number_machines_state_one for c in cells_matching_time])
    count_machines_state_two = sum([c.number_machines_state_two for c in cells_matching_time])
    count_machines_state_three = sum([c.number_machines_state_three for c in cells_matching_time])
    count_machines_total = count_machines_state_one + count_machines_state_two + count_machines_state_three

    count_machine_state_one_to_attacker = count_machines_state_one / count_attacker
    count_machine_state_two_to_attacker = count_machines_state_two / count_attacker
    count_machine_state_three_to_attacker = count_machines_state_three / count_attacker
    count_machine_total_to_attacker = count_machines_total / count_attacker
    # Intend to get count of all cells, but when no more cell label = 0, just get all cells that are label 1.
    count_neighbours_diff_to_target = 
        length(subset_cell_target) != 0 ? 
        sum([c.number_neighbours_different_label for c in subset_cell_target]) : 
        count_attacker
    count_neighbours_diff_to_attacker = sum([c.number_neighbours_different_label for c in subset_cell_attacker])

    return PopulationCountsPerTimeStep(
        time_to_match,
        count_target,
        count_attacker,
        count_total,
        count_machines_state_one,
        count_machines_state_two,
        count_machines_state_three,
        count_machines_total,
        count_machine_state_one_to_attacker,
        count_machine_state_two_to_attacker,
        count_machine_state_three_to_attacker,
        count_machine_total_to_attacker,
        count_neighbours_diff_to_target,
        count_neighbours_diff_to_attacker)
end

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
    2. `count_label_0`
    3. `count_label_1`
    4. `count_total`
    5. `count_machines_state_one`
    6. `count_machines_state_two`
    7. `count_machines_state_three`
    8. `count_machines_total`
    9. `count_machine_state_one_to_pop_01`
    10. `count_machine_state_two_to_pop_01`
    11. `count_machine_state_three_to_pop_01`
    12. `count_machine_total_to_pop_01`
    13. `count_neighbours_diff_to_cell_0`
    14. `count_neighbours_diff_to_cell_1`
"""
function get_population_counts_per_time_step(cell_data,time_to_match)

    cells_matching_time = filter(x -> x.time == time_to_match,cell_data)
    subset_cell_label_0 = filter(x->x.cell_type_label == 0,cells_matching_time)
    subset_cell_label_1 = filter(x->x.cell_type_label == 1,cells_matching_time)


    count_label_0 = count(x -> x.cell_type_label == 0,cells_matching_time)
    count_label_1 = count(x -> x.cell_type_label == 1,cells_matching_time)
    count_total = count_label_0 + count_label_1


    count_machines_state_one = sum([c.number_machines_state_one for c in cells_matching_time])
    count_machines_state_two = sum([c.number_machines_state_two for c in cells_matching_time])
    count_machines_state_three = sum([c.number_machines_state_three for c in cells_matching_time])
    count_machines_total = count_machines_state_one + count_machines_state_two + count_machines_state_three

    count_machine_state_one_to_pop_01 = count_machines_state_one / count_label_1
    count_machine_state_two_to_pop_01 = count_machines_state_two / count_label_1
    count_machine_state_three_to_pop_01 = count_machines_state_three / count_label_1
    count_machine_total_to_pop_01 = count_machines_total / count_label_1
    count_neighbours_diff_to_cell_0 = sum([c.number_neighbours_different_label for c in subset_cell_label_0])
    count_neighbours_diff_to_cell_1 = sum([c.number_neighbours_different_label for c in subset_cell_label_1])

    return PopulationCountsPerTimeStep(
        time_to_match,
        count_label_0,
        count_label_1,
        count_total,
        count_machines_state_one,
        count_machines_state_two,
        count_machines_state_three,
        count_machines_total,
        count_machine_state_one_to_pop_01,
        count_machine_state_two_to_pop_01,
        count_machine_state_three_to_pop_01,
        count_machine_total_to_pop_01,
        count_neighbours_diff_to_cell_0,
        count_neighbours_diff_to_cell_1)
end

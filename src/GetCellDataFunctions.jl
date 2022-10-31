"""
    Function to get cell data when the data type is about the cell populations, called `MachineState`. This name represents the saved file name from Chaste.

    The purpose of the function is to take a raw matrix, which was read into memory through the `readdlm` function from `DelimitedFiles` package. 

    Each time step in the Chaste simulation saves to file $N$ elements. These elements are represented in the struct `MachineStateCellProperties` and contains fieldnames:
    1. `time`: the time value.
    2. `location_index`: Index of cell related to its location in the simulated domain.
    3. `cell_id`: Identification label of the cell.
    4. `cell_type_label`: Cell type, can be $0$ or $1$.
    5. `number_machines_state_one`: Count of total machines in all cells in state $1$.
    6. `number_machines_state_two`: Count of total machines in all cells in state $2$.
    7. `number_machines_state_three`: Count of total machines in all cells in state $3$.
    8. `number_neighbours_different_label`: Count of neighbouring cells that have a different label to the identified cell.
    

    This function returns a vector, each row is a unique time step and cell to encapsulate the time series of the population.
"""
function get_cell_data(::MachineState,
    cell_population_data::Matrix{Any})
    
    # Gather properties of cells 
    number_cell_properties = 7
    total_number_cells = Int((size(cell_population_data,2)-1)/number_cell_properties)
    total_number_time_steps = size(cell_population_data,1)
    
    # Separate time vector and state properties
    simulation_time_only = cell_population_data[:,1]
    simulation_no_time = cell_population_data[:,2:end]


    # Separate each cells across each rows
    # Get index of first and last cell properties
    cell_block_per_row = []
    for i in 1:total_number_cells
        first = number_cell_properties*i - 6
        push!(cell_block_per_row,first)
        second = number_cell_properties*i
        push!(cell_block_per_row,second)
    end


    
    # Separate each cell at each time step into a vector
    cells_data = []
    for j = 1:total_number_time_steps, i = 1:total_number_cells
        first = 2*i - 1
        second = 2*i
        cell_data = simulation_no_time[j,cell_block_per_row[first]:cell_block_per_row[second]]
        cell_time_data = vcat(simulation_time_only[j],cell_data)
        push!(cells_data,cell_time_data)
    end

    # Get cell data size counting at each time step
    size_cell_data = size(cells_data,1)

    # Get cell properties without time present
    non_time_cells = [i[2:end] for i in cells_data]
    
    # Assert that the intersection of checking empty rows is not empty
    # Essentially testing the data read is not empty as a whole
    @assert any(map(x->all(isempty.(x)),non_time_cells)) "Population starts smaller than it finishes, there should be empty cells"

    # Find the indices where without time we find empty vectors
    # Remove empty vectors
    # Assert no empty cell remain
    empty_indx = findall(map(x->all(isempty.(x)),non_time_cells))
    non_zero_elements = [cells_data[i] for i = 1:size_cell_data if i ∉ empty_indx]
    @assert all([i[2] == i[3] for i in non_zero_elements]) "The second and third data point of the cell is the location index and cell id, these should ne identical"
    @assert all([all([!isempty(j) for j in i]) for i in non_zero_elements]) "No one of the elements per time should be empty"

    # Map each vectory to the MachineStateCellProperties struct
    cell_properties = map(x -> MachineStateCellProperties(x...),non_zero_elements)
    return cell_properties
end




"""
    Function to get cell machine data when the data type is about the cell populations, called `MachineData`. This name represents the saved file name from Chaste.

    The purpose of the function is to take a raw matrix, which was read into memory through the `readdlm` function from `DelimitedFiles` package. 

    Each time step in the Chaste simulation saves to file $N$ elements. These elements are represented in the struct `MachineDataProperties` and contains fieldnames:
    1. `time`: the time value.
    2. `cell_id`: Identification label of the cell.
    3. `cell_label`: Cell type, can be $0$ or $1$.
    4. `machine_id`: Identification label of the machine.
    5. `machine_state`: State of machine identified.
    6. `machine_coord_x`: Coordinate of $x$ location.
    7. `machine_coord_y`: Coordinate of $y$ location.
    8. `machine_coord_z`: Coordinate of $z$ location.
    


    This function returns a vector, each row is a unique time step, machine id and cell to encapsulate the time series of the population.
"""
function get_cell_data(::MachineData,
    machine_population_data::Matrix{Any})
    
    # Gather properties of cells 
    number_machine_properties = 7
    total_number_machines = Int((size(machine_population_data,2)-1)/number_machine_properties)
    total_number_time_steps = size(machine_population_data,1)
    
    # Separate time vector and state properties
    simulation_time_only = machine_population_data[:,1]
    simulation_no_time = machine_population_data[:,2:end]


    # Separate each cells across each rows
    # Get index of first and last cell properties
    machine_block_per_row = []
    for i in 1:total_number_machines
        first = number_machine_properties*i - 6
        push!(machine_block_per_row,first)
        second = number_machine_properties*i
        push!(machine_block_per_row,second)
    end


    
    # Separate each cell at each time step into a vector
    machines_data = []
    for j = 1:total_number_time_steps, i = 1:total_number_machines
        first = 2*i - 1
        second = 2*i
        machine_data = simulation_no_time[j,machine_block_per_row[first]:machine_block_per_row[second]]
        machine_time_data = vcat(simulation_time_only[j],machine_data)
        push!(machines_data,machine_time_data)
    end

    # Get cell data size counting at each time step
    size_machine_data = size(machines_data,1)

    # Get cell properties without time present
    non_time_machines = [i[2:end] for i in machines_data]
    
    # Assert that the intersection of checking empty rows is not empty
    # Essentially testing the data read is not empty as a whole
    @assert any(map(x->all(isempty.(x)),non_time_machines)) "Population starts smaller than it finishes, there should be empty cells"

    # Find the indices where without time we find empty vectors
    # Remove empty vectors
    # Assert no empty cell remain
    empty_indx = findall(map(x->all(isempty.(x)),non_time_machines))
    non_zero_elements = [machines_data[i] for i = 1:size_machine_data if i ∉ empty_indx]
    #@assert all([i[2] == i[3] for i in non_zero_elements]) "The second and third data point of the cell is the location index and cell id, these should ne identical"
    @assert all([all([!isempty(j) for j in i]) for i in non_zero_elements]) "No one of the elements per time should be empty"

    # Map each vectory to the MachineStateCellProperties struct
    machine_properties = map(x -> MachineDataProperties(x...),non_zero_elements)
    return machine_properties
end
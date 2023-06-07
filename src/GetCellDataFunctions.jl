"""
    Function to get cell data when the data type is about the cell populations, called `MachineState`. This name represents the saved file name from Chaste.

    The purpose of the function is to take a raw matrix, which was read into memory through the `readdlm` function from `DelimitedFiles` package. 

    Each time step in the Chaste simulation saves to file `N` elements. These elements are represented in the struct `MachineStateCellProperties` and contains fieldnames:
    1. `time`: the time value.
    2. `location_index`: Index of cell related to its location in the simulated domain.
    3. `cell_id`: Identification label of the cell.
    4. `cell_type_label`: Cell type, can be `0` or `1`.
    5. `number_machines_state_one`: Count of total machines in all cells in state `1`.
    6. `number_machines_state_two`: Count of total machines in all cells in state `2`.
    7. `number_machines_state_three`: Count of total machines in all cells in state `3`.
    8. `number_neighbours_different_label`: Count of neighbouring cells that have a different label to the identified cell.
    

    This function returns a vector, each row is a unique time step and cell to encapsulate the time series of the population.
"""
function get_cell_data(::MachineState,
    cell_population_data::Matrix)
    
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
    empty_machine_vec = map(x->!all(isempty.(x)),non_time_cells)
     
     @assert any(empty_machine_vec) "Population starts smaller than it finishes, there should be empty cells"

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

    Each time step in the Chaste simulation saves to file `N` elements. These elements are represented in the struct `MachineDataProperties` and contains fieldnames:
    1. `time`: the time value.
    2. `cell_id`: Identification label of the cell.
    3. `cell_label`: Cell type, can be `0` or `1`.
    4. `machine_id`: Identification label of the machine.
    5. `machine_state`: State of machine identified.
    6. `machine_coord_x`: Coordinate of `x` location.
    7. `machine_coord_y`: Coordinate of `y` location.
    8. `machine_coord_z`: Coordinate of `z` location.
    


    This function returns a vector, each row is a unique time step, machine id and cell to encapsulate the time series of the population.
"""
function get_cell_data(::MachineData,
    machine_population_data::Matrix)
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
     empty_machine_vec = map(x->!all(isempty.(x)),non_time_machines)
     
     @assert any(empty_machine_vec) "Population starts smaller than it finishes, there should be empty cells"
 
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

"""
    Get cell deata from the CellAgeData struct data - 
"""
 function get_cell_data(::CellAgeData,
    cell_population_data::Matrix)
    
    # Gather properties of cells 
    number_cell_properties = 5
    cell_property_index_space = number_cell_properties - 1
    total_number_cells = Int((size(cell_population_data,2)-1)/number_cell_properties)
    total_number_time_steps = size(cell_population_data,1)
    
    # Separate time vector and state properties
    simulation_time_only = cell_population_data[:,1]
    simulation_no_time = cell_population_data[:,2:end]


    # Separate each cells across each rows
    # Get index of first and last cell properties
    cell_block_per_row = []
    for i in 1:total_number_cells
        first = number_cell_properties*i - cell_property_index_space
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
    empty_machine_vec = map(x->!all(isempty.(x)),non_time_cells)
     
     @assert any(empty_machine_vec) "Population starts smaller than it finishes, there should be empty cells"

    # Find the indices where without time we find empty vectors
    # Remove empty vectors
    # Assert no empty cell remain
    empty_indx = findall(map(x->all(isempty.(x)),non_time_cells))
    non_zero_elements = [cells_data[i] for i = 1:size_cell_data if i ∉ empty_indx]
    @assert all([all([!isempty(j) for j in i]) for i in non_zero_elements]) "No one of the elements per time should be empty"

    # Map each vectory to the MachineStateCellProperties struct
    cell_properties = map(x -> CellAgeProperties(x...),non_zero_elements)
    return cell_properties
end









 """
    Extract the data frame from the cell state data
"""
 


function extract_dataframe_names(::MachineState, data)
    @chain data begin
        DataFrame(CellData = _) 
        @rtransform :time = :CellData.time 
        @rtransform :location_index = :CellData.location_index 
        @rtransform :cell_id = :CellData.cell_id
        @rtransform :cell_type_label = :CellData.cell_type_label
        @rtransform :number_machines_state_one = :CellData.number_machines_state_one 
        @rtransform :number_machines_state_two = :CellData.number_machines_state_two
        @rtransform :number_machines_state_three = :CellData.number_machines_state_three
        @rtransform :number_neighbours_different_label = :CellData.number_neighbours_different_label
        select(_,Not(:CellData))
        @transform :time = :time .- :time[1]
    end
end

"""
    Extract the data frame from the machine data
"""
 function extract_dataframe_names(::MachineData, data)
    @chain data begin
        DataFrame(CellData = _) 
        @rtransform :time = :CellData.time
        @rtransform :cell_id = :CellData.cell_id
        @rtransform :cell_type_label = :CellData.cell_label
        @rtransform :machine_id = :CellData.machine_id
        @rtransform :machine_state = :CellData.machine_state
        @rtransform :machine_coord_x = :CellData.machine_coord_x
        @rtransform :machine_coord_y = :CellData.machine_coord_y
        @rtransform :machine_coord_z = :CellData.machine_coord_z
        select(_,Not(:CellData))
        @transform :time = :time .- :time[1]
    end
end



"""
    Extract data from the cell age data into the corresponding data frame
"""
function extract_dataframe_names(::CellAgeData, data)
    @chain data begin
        DataFrame(CellData = _) 
        @rtransform :time = :CellData.time 
        @rtransform :location_index = :CellData.location_index 
        @rtransform :cell_centre_coord_x = :CellData.cell_centre_coord_x
        @rtransform :cell_centre_coord_y = :CellData.cell_centre_coord_y
        @rtransform :cell_centre_coord_z = :CellData.cell_centre_coord_z
        @rtransform :cell_age = :CellData.cell_age        
        select(_,Not(:CellData))
        @transform :time = :time .- :time[1]
        @transform :cell_age = :cell_age .- :cell_age[1]
    end
end



"""
    From machine path to data frame and use the fields of MachineData
    for each column
"""
function get_cell_dataframe(T,path)
    df = @chain path begin
        readdlm(_) 
        get_cell_data(T,_) 
        extract_dataframe_names(T,_)
    end
    return df
end


""" 
    From String -> DataFrame    
    String: path to the chaste .dat file from 
    simulation of the cell populations
    Output dataframe will have the time series
    time, target, attacker, m₁, m₂, m₃, m, neigh_diff_to_attacker, neigh_diff_to_target
"""
function get_cell_dataframe_TS(T,chaste_path_dat_file)
    
    # Get cell dataframe and add attacker and target column
    cell_df = @chain chaste_path_dat_file begin
        get_cell_dataframe(T,_)
        @rtransform :cell_type = :cell_type_label == 1 ? "attacker" : "target"
    end

    # extract attacker and target time series
    attacker_target_TS = @chain cell_df begin
        @by([:cell_type,:time],
            :N = size(:cell_type,1))
        unstack(_,
        :time,:cell_type,:N,fill=0)
    end

    # extract machines time series
    machines_TS = @chain cell_df begin
        @rsubset :cell_type_label == 1
        @by([:cell_type_label,:time],
            :N = size(:cell_type_label,1),
            :m₁ = sum(:number_machines_state_one),
            :m₂ = sum(:number_machines_state_two),
            :m₃ = sum(:number_machines_state_three))
        @rtransform :mₜ = :m₁ + :m₂ + :m₃ 
        @select $(Not([:cell_type_label,:N]))
    end


    # extract target time series
    neighbours_TS = @chain cell_df begin
        @select :time :cell_type :number_neighbours_different_label 
        @by([:time, :cell_type],
            :neigh = sum(:number_neighbours_different_label))
        unstack(_,
            :time,:cell_type,:neigh,
            fill=0,
            renamecols=x->Symbol(:neigh_diff_to_, x))
    end


    # join all time series together
    all_TS = @chain attacker_target_TS begin
        innerjoin(_,neighbours_TS,on=:time)
        innerjoin(_,machines_TS,on=:time)
        @select :time :target :attacker :m₁ :m₂ :m₃ :mₜ :neigh_diff_to_attacker :neigh_diff_to_target
    end

    return all_TS
end






""" 
    From String -> DataFrame    
    String: path to the chaste .dat file from 
    simulation of the cell populations
    Output dataframe will have the time series
    time, target, attacker, m₁, m₂, m₃, m, neigh_diff_to_attacker, neigh_diff_to_target
"""
function get_cell_dataframe_TS(::DataAllCellLabels,chaste_path_dat_file)
    label_dict = Dict(0 => "target", 1 => "attacker", 3 => "infected")
    map_label(x) = label_dict[x]
    # Get cell dataframe and add attacker and target column
    cell_df = @chain chaste_path_dat_file begin
        get_cell_dataframe(MachineState(),_)
        @rtransform :cell_type = map_label(:cell_type_label)
    end
    


    # extract attacker and target time series
    attacker_target_TS = @chain cell_df begin
        @by([:cell_type,:time],
            :N = size(:cell_type,1))
        unstack(_,
        :time,:cell_type,:N,fill=0)
    end

    # extract machines time series
    machines_TS = @chain cell_df begin
        @rsubset :cell_type_label == 1
        @by([:cell_type_label,:time],
            :N = size(:cell_type_label,1),
            :m₁ = sum(:number_machines_state_one),
            :m₂ = sum(:number_machines_state_two),
            :m₃ = sum(:number_machines_state_three))
        @rtransform :mₜ = :m₁ + :m₂ + :m₃ 
        @select $(Not([:cell_type_label,:N]))
    end


    # extract target time series
    neighbours_TS = @chain cell_df begin
        @select :time :cell_type :number_neighbours_different_label 
        @by([:time, :cell_type],
            :neigh = sum(:number_neighbours_different_label))
        unstack(_,
            :time,:cell_type,:neigh,
            fill=0,
            renamecols=x->Symbol(:neigh_diff_to_, x))
    end


    # join all time series together
    all_TS = @chain attacker_target_TS begin
        innerjoin(_,neighbours_TS,on=:time)
        innerjoin(_,machines_TS,on=:time)
        #@select :time :target :attacker :infected :m₁ :m₂ :m₃ :mₜ :neigh_diff_to_attacker :neigh_diff_to_target
    end
#=

    try
        @chain all_TS begin
            name(_)
            filter(x -> occursin("infected",x),_)
        end
    catch
        all_TS.infected .= 0
    end
    
    all_TS = @chain all_TS begin
        @select :time :target :attacker :infected :m₁ :m₂ :m₃ :mₜ :neigh_diff_to_attacker :neigh_diff_to_target
    end
=#
    return all_TS
end










"""
    The area a capule with a cylinder length l and radius r
"""
cₐ(r,l) = π*r^2 + 2*r*l


"""
    The length a capule with a cylinder length l and radius r
"""
cₗ(r,l)  = 2*r + l


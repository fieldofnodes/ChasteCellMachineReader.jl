


today_nog_gap() = replace(today() |> string,"-" => "")

"""
    Normalise a vector, x such that x∈[0,1]
"""
function get_norm_vec(x) 
    return (x .- minimum(x)) ./ (maximum(x) .- minimum(x))
end

"""
    Shift a U(0,1) to U(a,b)
"""
function shift_scale_uniform(min_value,max_value)
    range_dif = abs(min_value - max_value)
    return range_dif*rand()+min_value
end

"""
    Get mean values from a U(a,b)
"""
function get_mean_cycle_time(cycle_min_time,cycle_max_time,runs)
    return mean([shift_scale_uniform(cycle_min_time,cycle_max_time) for i in 1:runs])
end




"""
    Convert a matrix to a dataframe.
    Use the `:auto` variable.
    `:x1` is the radius of the target boundary
    `:x2`` is the competition value
    Then the variables from `:x3` onwards is shifted
    to `:x1` and so on.
"""
function convert_mat_to_df(mat)
    sz_mat = size(mat,2)
    df = @chain mat begin
        DataFrame(_,:auto) 
        rename(_,:x1 => :radius)
        rename(_,:x2 => :comp)
        rename(_,["x$i" => "x$(i-2)" for i in 3:sz_mat])
    end
    return df
end


"""
    NOT GENERALISED
    Specific to path such that folder names separated by `/` or `\\` are removed by two times,then the remaining path removed all folder names until a single folder name remains. Here the realisation number is displayed.
"""
function get_realisation_number(data_path)
    realisation = @chain data_path begin
        dirname(_)
        dirname(_)
        basename(_)
        replace(_,"Realisation_" => "")
        parse(Float64,_)
    end
    return realisation
end


"""
    Take a dataframe and the number of realisations and includes this as a 
    column in the dataframe.
"""
function add_realisations(
    data::DataFrame,
    realisations::Int64)
    data.realisations .= realisations
    return data
end




"""
    Compute the standard error from the standard deviation
    and the count of observartions.
"""
function std_error(σ,n)
    return σ/√n
end



"""
    If a value is missing, then output a zero
    Missing -> 0
"""
missing_to_zero(m::Missing) = 0


"""
    Function that takes an input and if missing than outputs 0 otherwise
    return what that number is.
"""
function if_missing_then_zero(x) 
    @assert isreal(x)|ismissing(x) "Input should be a Real number of missing"
    x̂ = ismissing(x) ? missing_to_zero(x) : x
    return x̂
end



"""
    Functions specifice to the chaste cell reader where we convert
    values per cell type label (0,1,3) and outputs to zero
    Function needs to be generalised.
"""
function convert_missing_to_zero(data) 
    df = @chain data begin 
            unstack(_,:cell_type_label,:N,renamecols=x->Symbol(:label_, x))
            @rtransform :label_0 = if_missing_then_zero(:label_0)
            @rtransform :label_3 = if_missing_then_zero(:label_3) 
            @rtransform :label_target = :label_0 + :label_3
            @select :time :label_1 :label_target
            stack(_,[:label_1, :label_target],:time)
            @rtransform :variable = :variable == "label_1" ? 1 : 0
            rename(_,:variable => :cell_type_label) 
            rename(_,:value => :N)
            @select :cell_type_label :time :N
        end
        return df
end




"""
    A once-off look behind patter when we want to change the patter
    of the look behind but still search for [0-9.]{1.5} after
"""
function input_neg_lookbehind(pattern)
    rx = Regex(string("(?<=",pattern,")[0-9.]{1,5}"))
    return rx
end

"""
    Take the regex patter, the string and the regex function to parse a 
    character into a Float64 type.
"""
function parse_regex_float_match(
    string_for_match,
    pattern_to_match,
    regex_match_function)

    rx = regex_match_function(pattern_to_match)
    m = match(rx,string_for_match,).match
    
    return parse(Float64,m)
end




"""
    Take cell path and return the folder name
    which comes from the basename of the filename
"""
function new_folder_name(cell_path)
    folder_name = @chain cell_path begin
        basename(_)    
        splitext(_)
        _[1]
        replace(_,"_machinestate"=>"")
    end
    return folder_name
end


"""
    Generate path from cell path to save to file
    input the cell path and an output descriptor 
"""
function generate_cell_output_path(cell_path,output_desc,file_ext)
    base = basename(cell_path)
    df_path = replace(base,".dat" => ".$(file_ext)")
    output = string(output_desc,"_",df_path)
    return output
end

"""
    From cell path, output description and new extension
    output the new path
"""
function make_output_path(base_path,cell_path,output_desc,file_ext)
    output_folder = new_folder_name(cell_path)
    output_file = generate_cell_output_path(cell_path,output_desc,file_ext)
    return  base_path *"/"*
            output_folder *"/"* 
            output_file
end



function convert_vecs_to_mat(data)
    return reduce(hcat,data.N_vec)
end

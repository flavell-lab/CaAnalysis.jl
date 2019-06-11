function highlight_stim(stim, α_highlight=0.1)
    for i = 1:size(stim, 1)
        gca().axvspan(stim[i,1], stim[i,2],
            alpha=α_highlight, facecolor="red", lw=0)
    end
end

"""
    get_idx_unit(data_dict::Dict, idx_unit)

Returns indicies of units to use.

Arguments
---------
* `data_dict`: data_dictionary
* `idx_unit`: indicies of units to use. Options:
  * `:all`: use all units
  * `:ok`: use pre-selected units "idx_ok"
  * `UnitRange{T}`: range. e.g. `2:50` uses units from 2 to 50
  * `Array{T,1}`: list. e.g. `[1, 5, 11, 32]`
"""
function get_idx_unit(data_dict::Dict, idx_unit)
    if idx_unit == :ok
        if !haskey(data_dict, "idx_good")
            error("Missing key \"idx_good\" in data_dict")
        end
        return data_dict["idx_good"]
    elseif idx_unit == :all
        return 1:size(data_dict["f"], 1)
    elseif isa(idx_unit, Union{UnitRange{T}, Array{T,1}} where T <: Signed)
        return idx_unit
    else
        error("Incorrect idx_unit. It must be UnitRange{T}, Array{T,1}, " *
        ":ok, or :all")
    end
end

"""
    get_idx_t(data_dict::Dict, idx_t)

Returns indicies of time points to use.

Arguments
---------
 * `data_dict`: data_dictionary
 * `idx_t`: indices of time points to use. Options:
   * `:all`: use all time points
   * `Tuple{T,T}`: number of time points to trim. e.g. `(50, 0)` removes the
   * first 50 points. `(50,100)` removes first 50 and last 100
   * `UnitRange{T}`: range. e.g. `50:500` uses time points from 50 to 500
   * `Array{T,1}`: list. e.g. `[1, 5, 11, 32]` to be used.
 """
function get_idx_t(data_dict::Dict, idx_t)
    if idx_t == :all
        return 1:size(data_dict["f"], 2)
    elseif isa(idx_t, Union{UnitRange{T}, Array{T,1}} where T <: Signed)
        return idx_t
    elseif isa(idx_t, Tuple{T,T} where T <: Signed)
        return 1 + idx_t[1]:(size(data_dict["f"], 2) - idx_t[2])
    else
        error("Incorrect idx_t. It must be UnitRange{T}, Array{T,1}, " *
        "Tuple{T, T} or :all where T <: Signed")
    end

end

"""
    get_data(data_dict::Dict; data_key="f_denoised", idx_unit=:ok, idx_t=:all)

Returns data given data_key, units to use, and time points.

Arguments
---------
* `data_dict`: data_dictionary
* `data_key`: data key to use
* `idx_unit`: see [`get_idx_unit()`](@ref)
* `idx_t`: see [`get_idx_t()`](@ref)
"""
function get_data(data_dict::Dict; data_key="f_denoised", idx_unit=:ok,
    idx_t=:all)
    unit_use = get_idx_unit(data_dict, idx_unit)
    t_use = get_idx_t(data_dict, idx_t)

    data_dict[data_key][unit_use, t_use]
end
function highlight_stim(stim, α_highlight=0.1)
    for i = 1:size(stim, 1)
        gca().axvspan(stim[i,1], stim[i,2],
            alpha=α_highlight, facecolor="red", lw=0)
    end
end

function get_unit_idx(data_dict, idx_unit)
    if idx_unit == :ok
        return data_dict["idx_good"]
    elseif idx_unit == :all
        return 1:size(data_dict["f"], 1)
    elseif isa(idx_unit, Union{UnitRange{T}, Array{T,1}} where T <: Signed)
        return idx_unit
    else
        return nothing
    end
end

"""
Divides the activity traces `traces::Dict` by the marker signal `marker_traces::Dict`
"""
function divide_by_marker_signal(traces::Dict, marker_traces::Dict)
    new_traces = Dict()
    for neuron in keys(traces)
        new_traces[neuron] = Dict()
        for frame in keys(traces[neuron])
            if traces[neuron][frame] > 0 && marker_traces[neuron][frame] > 0
                new_traces[neuron][frame] = traces[neuron][frame] / marker_traces[neuron][frame]
            end
        end
    end
    return new_traces
end


"""
Gets the background camera intensity from the MHD files to be subtracted.

# Arguments
- `timepts`: Time points to get background of
- `get_basename::Function`: Function that outputs the name of the MHD file from the time and channel
- `mhd_path::String`: Path to MHD files
- `ch::Int`: Channel
"""
function get_background(timepts, get_basename::Function, mhd_path::String, ch::Int)
    frame_bkg = Dict()
    @showprogress for t in timepts
        frame_bkg[t] = median(read_img(MHD(joinpath(mhd_path, get_basename(t, ch))*".mhd")))
    end
    return frame_bkg
end

"""
Subtracts from `traces` the corresponding background `frame_bkg`.
"""
function bkg_subtract(traces::Dict, frame_bkg)
    new_traces = Dict()
    warning = false
    for neuron in keys(traces)
        new_traces[neuron] = Dict()
        for t in keys(traces[neuron])
            new_traces[neuron][t] = traces[neuron][t] - frame_bkg[t]
            if new_traces[neuron][t] < 0
                warning = true
                new_traces[neuron][t] = 0
            end
        end
    end
    if warning
        @warn "Some neurons had lower intensity than the background."
    end
    return new_traces
end

"""
Normalizes traces.

# Arguments
- `traces::Dict`: Traces to normalize

# Optional keyword arguments
- `zero::Bool`: Whether to normalize to zero (so activities that are less than average would be negative). Default true.
- `fn::Function`: Function that determines average neuron activity
"""
function normalize_traces(traces::Dict; zero::Bool=true, fn::Function=mean)
    new_traces = Dict()
    for neuron in keys(traces)
        if length(keys(traces[neuron])) == 0
            continue
        end
        new_traces[neuron] = Dict()
        s = fn([x for x in values(traces[neuron])])
        for t in keys(traces[neuron])
            new_traces[neuron][t] = traces[neuron][t] / s
            if zero
                new_traces[neuron][t] -= 1
            end
        end
    end
    return new_traces
end

"""
Interpolates data at missing time points.

# Arguments
 - `traces::Dict`: Dictionary containing for each neuron a set of time points and activity values at those time points
 - `t_range`: Time points to interpolate to. Time points outside this range will be set to `fill_val`
 - `itp_method` (optional, default `Linear()`): Interpolation method.
 - `extrap_method` (optional, default `Interpolations.Flat()`): Extrapolation method (for interpolating data points outside all time points the neuron was detected in)
 - `fill_val` (optional, default `NaN`): Value to fill time points outside the interpolation range.
"""
function interpolate_traces(traces::Dict, t_range; itp_method=Linear(), extrap_method=Interpolations.Flat(), fill_val=NaN)
    new_traces = Dict()
    max_t = maximum(t_range)
    for neuron in keys(traces)
        t_vals = sort(collect(keys(traces[neuron])))
        activity_vals = [traces[neuron][t] for t in t_vals]
        itp = extrapolate(interpolate((t_vals,), activity_vals, Gridded(itp_method)), extrap_method)
        new_traces[neuron] = [(t in t_range) ? Float64(itp(t), RoundNearest) : fill_val for t in 1:max_t]
    end
    return new_traces
end


"""
Z-scores traces.
"""
function zscore_traces(traces::Dict)
    new_traces = Dict()
    for neuron in keys(traces)
        m = mean(collect(values(traces[neuron])))
        s = std(collect(values(traces[neuron])))
        new_traces[neuron] = Dict()
        for t in keys(traces[neuron])
            new_traces[neuron][t] = (traces[neuron][t] - m) / s
        end
    end
    return new_traces
end

"""
Applies multiple data processing steps to the traces. The order of processing steps is:

- Background-subtraction
- Delete low S/N neurons
- Delete neurons detected in too few time points
- Interpolate traces at missing data points
- Denoise
- Bleach-correct
- Divide activity by marker channel
- Normalize
- Zscore

# Arguments
- `activity_traces::Dict`: traces in the activity channel
- `marker_traces::Dict`: traces in the marker channel
- `threshold::Real`: Number of timepoint detections necessary to include a neuron in the analysis

# Keyword arguments
- `activity_bkg`: Background in activity channel. If left blank, background will not be subtracted.
- `marker_bkg`: Background in marker channel. If left blank, background will not be subtracted.
- `min_intensity::Real`: Minimum average intensity in the activity channel for a neuron (after background subtraction).
    Neurons with less than this much signal will be removed. Default 0. 
- `interpolate_t_range`: Time points to interpolate to (default `nothing` which skips data interpolation)
- `denoise::Bool`: Whether to apply a total variation denoising step.
- `bleach_corr::Bool`: Whether to bleach-correct the traces.
- `divide::Bool`: Whether to divide the activity channel traces by the marker channel traces.
- `normalize::Bool`: Whether to normalize the traces.
- `normalize_fn::Function`: Function to use to get "average" activity when normalizing traces.
- `zscore::Bool`: Whether to z-score the traces.
"""
function process_traces(activity_traces::Dict, marker_traces::Dict, threshold::Real; activity_bkg=nothing, marker_bkg=nothing,
        min_intensity::Real=0, interpolate_t_range=nothing, denoise::Bool=false, bleach_corr::Bool=false, divide::Bool=false, normalize::Bool=false, normalize_fn::Function=mean,
        zscore::Bool=false)

    interpolate = !isnothing(interpolate_t_range)

    activity_traces = copy(activity_traces)
    marker_traces = copy(marker_traces)
    # background subtract activity
    if !isnothing(activity_bkg)
        activity_traces = bkg_subtract(activity_traces, activity_bkg)
    end
    # background subtract marker
    if !isnothing(marker_bkg)
        marker_traces = bkg_subtract(marker_traces, marker_bkg)
    end
    # delete neurons with too low S/N
    neuron_rois = [roi for roi in keys(activity_traces) if
            length(values(activity_traces[roi])) > 0 && mean(values(activity_traces[roi])) > min_intensity]
    for roi in keys(activity_traces)
        if !(roi in neuron_rois) || length(activity_traces[roi]) < threshold
            delete!(activity_traces, roi)
        end
    end
    for roi in keys(marker_traces)
        if !(roi in neuron_rois) || length(marker_traces[roi]) < threshold
            delete!(marker_traces, roi)
        end
    end


    # interpolate traces
    if interpolate
        activity_traces = interpolate_traces(activity_traces, interpolate_t_range)
        marker_traces = interpolate_traces(marker_traces, interpolate_t_range)
    end

    # make traces array for futher processing
    all_traces = [activity_traces, marker_traces]

    for (idx, traces) in enumerate(all_traces)
        traces_arr, hmap, valid_rois = make_traces_array(traces, threshold=threshold, replace_blank=true)

        processed_traces_arr = traces_arr
        data_dict = Dict()
        data_dict["f"] = traces_arr
        data_dict["idx_ok"] = 1:size(data_dict["f"])[1]

        # denoise
        if denoise
            denoiser = DenoiserGSTV(100, 10.)
            denoise!(data_dict, denoiser.f)
            processed_traces_arr = data_dict["f_denoised"]
        else
            data_dict["f_denoised"] = data_dict["f"]
        end

        # bleach-correct
        if bleach_corr
            if interpolate
                fit_bleach!(data_dict, idx_t=interpolate_t_range)
            else
                fit_bleach!(data_dict)
            end
            processed_traces_arr = data_dict["f_bleach"]
        end

        # convert back to dictionary
        all_traces[idx] = Dict()
        for i in 1:length(valid_rois)
            all_traces[idx][valid_rois[i]] = Dict()
            for t in keys(activity_traces[valid_rois[i]])
                all_traces[idx][valid_rois[i]][t] = processed_traces_arr[i,t]
            end
        end
    end

    # divide activity by marker
    if divide
        all_traces[1] = divide_by_marker_signal(all_traces[1], all_traces[2])
    end

    if normalize
        for idx=1:2
            all_traces[idx] = normalize_traces(all_traces[idx], fn=normalize_fn)
        end
    end

    if zscore
        for idx=1:2
            all_traces[idx] = zscore_traces(all_traces[idx])
        end
    end


    return all_traces
end


"""
Gets all values from a set of processed traces. Ensure that the traces have been normalized prior to using this function!
"""
function get_all_values(processed_traces)
    trace_vals = []
    for neuron in keys(processed_traces)
        append!(trace_vals, collect(values(processed_traces[neuron])))
    end
    return trace_vals
end

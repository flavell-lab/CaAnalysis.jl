module CaAnalysis

using Statistics, PyPlot, HDF5, Dierckx, ProgressMeter, Distances

include("init.jl")
include("data.jl")
include("denoise.jl")
include("bleach.jl")
include("heatmap.jl")
include("single_unit.jl")
include("multivar.jl")
include("util/util.jl")
include("util/processing.jl")
include("util/plot.jl")

export import_data,

    # denoise.jl
    denoise_gstv,
    denoise_gstv!,
    denoise_tf,
    denoise_tf!,

    # bleach.jl
    fit_bleach,
    fit_bleach!,

    # heatmap.jl
    plot_cluster_cost,
    plot_heatmap,

    # unit.jl
    plot_unit_cor,
    get_idx_unit,
    get_idx_stim,
    get_idx_t,
    get_data,
    get_stim,

    # single_unit.jl
    compute_unit_cor,
    plot_unit_cor,

    # multivar.jl
    pca,
    plot_pca_var,
    plot_statespace_component,
    plot_statespace_3d,
    plot_statespace_2d,

    # util/processing.jl
    grad,
    integrate,
    standardize,
    discrete_diff,
    chain_process,

    # util/plot.jl
    highlight_stim

end # module

module CaAnalysis

using Statistics, PyPlot, HDF5, Dierckx, ProgressMeter, MultivariateStats

include("init.jl")
include("data.jl")
include("bleach.jl")
include("heatmap.jl")
include("single_unit.jl")
include("multivar.jl")
include("util/util.jl")
include("util/processing.jl")
include("util/plot.jl")
include("util/denoise.jl")

export import_data,
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
    plot_statespace_2d_stim,

    # util/processing.jl
    derivative,
    integrate,
    standardize,
    discrete_diff,
    chain_process,

    # util/plot.jl
    highlight_stim,


    # util/denoise.jl
    preview_denoise,
    denoise,
    denoise!,
    DenoiserTrendfilter,
    DenoiserGSTV

end # module

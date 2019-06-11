module CaAnalysis

using Statistics, PyPlot, HDF5, Dierckx

include("data.jl")
include("denoise.jl")
include("bleach.jl")
include("heatmap.jl")
include("util.jl")
include("single_unit.jl")
include("multivar.jl")

export import_data,

    # denoise.jl
    denoise_gstv,
    denoise_tf,

    # bleach.jl
    fit_bleach,

    # heatmap.jl
    plot_cluster_cost,
    plot_heatmap,

    # unit.jl
    plot_unit_cor,
    get_idx_unit,
    get_idx_t,
    get_data,
    highlight_stim,

    # single_unit.jl
    compute_unit_cor,
    plot_unit_cor,

    # multivar.jl
    pca,
    plot_statespace_component,
    plot_statespace_3d,
    plot_pca_var

end # module

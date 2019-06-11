using Clustering

"""
    plot_cluster_cost(args)

documentation
"""
function plot_cluster_cost(f::Array{T,2}, n=20) where T
    cost_list = zeros(n)

    for i = 1:n
        clustered=kmeans(Array(f'), i+1, init=:kmcen, maxiter=200,
            tol=1.0e-8)
        cost_list[i] = clustered.totalcost
    end

    subplot(2,1,1)
    plot(2:n+1, cost_list, "ko-")
    ylabel("Total cost")
    xlabel("Number of clusters")
    xlim(0,n+2)
    xticks(1:2:n+1)

    subplot(2,1,2)
    plot(3:n+1, diff(cost_list), "ko-")
    ylabel("Δ(Total cost)")
    xlabel("Number of clusters")
    xlim(0,n+2)
    xticks(1:2:n+1)

    tight_layout()

    nothing
end

function plot_cluster_cost(data_dict::Dict, n=20; data_key="f_bleach",
    idx_unit=:ok, idx_t=:all)
    f = get_data(data_dict, data_key=data_key, idx_unit=idx_unit,
        idx_t=idx_t)

    plot_cluster_cost(f, n)

    nothing
end

function plot_heatmap(f::Array{T,2}, n=10; vmin=0.5, vmax=1.5,
        cmap="magma", vol_rate=0.75) where T

    f_plot = f ./ mean(f, dims=2)

    clustered = kmeans(Array(f_plot'), n, init=:kmcen, maxiter=200,
        tol=1.0e-8)
    kmeans_order = sortperm(clustered.assignments);


    imshow(f_plot[kmeans_order,:], aspect=7, vmin=vmin, vmax=vmax,
        cmap=cmap)
    xlabel("Time (min)")
    ylabel("Unit")

    Δ_minute = 5 # minute
    n_vol = size(f, 2)
    t_rg = 0:floor(Int, n_vol * vol_rate / 60 / Δ_minute)
    x_rg = collect(t_rg) * 60 * Δ_minute / vol_rate
    xticks(x_rg, Δ_minute * collect(t_rg))

    nothing
end

function plot_heatmap(data_dict, n=10; cmap="magma", data_key="f_bleach",
    idx_unit=:ok, idx_t=:all, vmin=0.5, vmax=1.5, vol_rate=0.75)
    f = get_data(data_dict, data_key=data_key, idx_unit=idx_unit,
        idx_t=idx_t)

    plot_heatmap(f, n; vmin=vmin, vmax=vmax, cmap=cmap, vol_rate=vol_rate)

    nothing
end

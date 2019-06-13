function pca(f)
    X = f

    M = fit(PCA, X)
    Yt = transform(M, X)
    (M, X, Yt)
end

function plot_pca_var(M, n=10)
    var_list = (M.prinvars / M.tprinvar)[1:n]

    PyPlot.bar(1:n, var_list, color="k", alpha=0.5)
    plot(1:n, cumsum(var_list),"ko-")
    yticks(0.2:0.2:1.0)
    ylim(0,1.0)
    xlabel("PC #")
    ylabel("Var exp %")

    nothing
end

function plot_statespace_component(Y, n=10; idx_stim=nothing, α_highlight=0.1)
    @assert(n % 2 == 0)

    for i = 1:n
        subplot(Int(n/2), 2, i)
        plot(Y[i, :])
        title("Component $i")

        if !isnothing(idx_stim)
            highlight_stim(idx_stim, α_highlight)
        end
    end
    tight_layout()

    nothing
end

function plot_statespace_3d(Y, prjax=[1,2,3])
    if length(prjax) != 3 || !(eltype(prjax) <: Signed)
        error("prjax should be list of 3 integers")
    end

    ax = gcf().add_subplot(111, projection="3d")
    ax.plot(Y[prjax[1],:], Y[prjax[2],:], Y[prjax[3],:], color="k", alpha=0.5)

    ax.view_init(45,15)
    xlabel("Axis 1")
    ylabel("Axis 2")
    zlabel("Axis 3")

    nothing
end

function plot_statespace_2d(Y, prjax=[1,2])
    if length(prjax) != 2 || !(eltype(prjax) <: Signed)
        error("prjax should be list of 2 integers")
    end

    plot(Y[prjax[1],:], Y[prjax[2],:], color="k", alpha=0.5)

    xlabel("Axis 1")
    ylabel("Axis 2")

    nothing
end

# function calc_statespace_dist(Y, prjax=[1,2,3], dist=Euclidean())
#     colwise(dist, Y[prjax,1:end-1], Y[prjax,2:end])
# end

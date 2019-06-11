function grad_1d(x, y)
    spl = Spline1D(x, y, k=3) # spline order 3
    derivative(spl, x)
end

function pca(f, gradient=true, standardize=true)
    X = f

    n_unit = size(f, 1)
    t_max = size(f, 2)

    if gradient
        for i = n_unit
            X[i, :] = grad_1d(1:t_max, X[i, :])
        end
    end

    if standardize
        σ = std(f, dims=2)
        μ = mean(f, dims=2)
        X = (f .- μ) ./ σ
    end

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
end

function plot_statespace_component(Y, n=10; idx_stim=nothing, α_highlight=0.1)
    @assert(n % 2 == 0)

    for i = 1:n
        subplot(Int(n/2), 2, i)
        plot(Y[i, :])
        title("Component $i"))

        if !isnothing(idx_stim)
            highlight_stim(idx_stim, α_highlight)
        end
    end
    tight_layout()
end

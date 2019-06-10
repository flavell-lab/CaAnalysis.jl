function pca(f, gradient=true, standardize=true)
    X = f

    if gradient
        X = diff(X, dims=2)
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

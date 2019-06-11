function ∇(x::Array{T,1}, y::Array{T,1}, nu=1) where T
    spl = Spline1D(x, y, k=3) # spline order 3
    derivative(spl, x, nu=nu)
end

function ∇(y::Array{T,1}, nu=1) where T
    x = 1:length(y)
    spl = Spline1D(x, y, k=3)
    derivative(spl, x, nu=nu)
end

function standardize(f::Array{T,2}) where T
    σ = std(f, dims=2)
    μ = mean(f, dims=2)
    (f .- μ) ./ σ
end

function grad(f::Array{T,2}, nu=1) where T
    f_grad = zero(f)
    n_unit = size(f, 1)

    for i = 1:n_unit
        f_grad[i, :] = ∇(f[i, :], nu)
    end

    f_grad
end

function discrete_diff(f::Array{T,2}) where T
    diff(f, dims=2)
end

function chain_process(f::Array, g_list)
    f_proc = f

    for g = g_list
        f_proc = g(f_proc)
    end

    f_proc
end

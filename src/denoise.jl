using TotalVariation, Lasso

function denoise_gstv(f::Array{T, 2}; n=50, λ=5.0) where T
    f_denoised = zero(f)

    for roi_ = 1:size(f, 1)
        signal = f[roi_, :]
        f_denoised[roi_, :] = TotalVariation.gstv(signal, n, λ)
    end

    f_denoised
end

function denoise_gstv(data_dict::Dict; n=50, λ=5.0)
    data_dict["f_denoised"] = denoise_gstv(data_dict["f"], n=n, λ=λ)

    nothing
end

function denoise_tf(f::Array{T, 2}; order=2, λ=2.0) where T
    f_denoised = zero(f)

    for roi_ = 1:size(f, 1)
        signal = f[roi_, :]
        f_denoised[roi_, :] = fit(TrendFilter, signal, order, λ).β
    end

    f_denoised
end

function denoise_tf(data_dict::Dict; order=2, λ=2.0)
    data_dict["f_denoised"] = denoised_tf(data_dict["f"], order=order, λ=λ)

    nothing
end

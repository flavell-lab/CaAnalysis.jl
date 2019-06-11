using LsqFit

@. bleach_model(x, p) = p[1] * exp(-x * p[2]) + p[3] * exp(-x * p[4])

"""
    fit_bleach(f, p0, plot_fit=true)

Fits double exponential bleaching model:
bleach_model(x, p) = p[1] * exp(-x * p[2]) + p[3] * exp(-x * p[4])

Arguments
---------
* `f`: data to fit the bleaching model
* `p0`: initial parameters (4 element array)
* `plot_fit`: plot fit result if true
"""

function fit_bleach(f, p0, plot_fit=true)
    @assert(length(p0) == 4)

    y = f ./ f[1]
    x = 0:length(y)-1

    fitted = curve_fit(bleach_model, x, y, p0)

    if plot_fit
        subplot(2,1,1)
        title("Fitted")
        plot(y, "k", alpha=0.75)
        plot(x, bleach_model(x, fitted.param), "r")

        subplot(2,1,2)
        title("Residual")
        plot(fitted.resid, "k", alpha=0.75)

        tight_layout()
        println("Fitted parameter:", fitted.param)
    end
    bleach_curve = bleach_model(x, fitted.param)
    fitted, bleach_curve ./ bleach_curve[1]
end

"""
    fit_bleach(f::Array{T,2}, p0, plot_fit=true)

Calculates mean activity across the units and fits the double exponential
bleaching model:
bleach_model(x, p) = p[1] * exp(-x * p[2]) + p[3] * exp(-x * p[4])

Arguments
---------
* `f`: N x T data array. N: number of units, T: number of time points.
* `p0`: initial parameters (4 element array)
* `plot_fit`: plot fit result if true
"""
function fit_bleach(f::Array{T,2}, p0, plot_fit=true) where T
    y = dropdims(mean(f, dims=1), dims=1)
    fit_bleach(y, p0, plot_fit)
end

"""
    fit_bleach(data_dict::Dict, p0, plot_fit=true; idx_unit=:ok,
        data_key="f_denoised")

Calculates mean activity across the units and fits the double exponential
bleaching model:
bleach_model(x, p) = p[1] * exp(-x * p[2]) + p[3] * exp(-x * p[4])

Arguments
---------
* `f`: N x T data array. N: number of units, T: number of time points.
* `p0`: initial parameters (4 element array)
* `plot_fit`: plot fit result if true
* `idx_unit`: `:ok`: only using pre-filtered units or array/range of indicies.
* `data_key`: key of data_dict to be used for fitting the model
"""

function fit_bleach(data_dict::Dict, p0, plot_fit=true; idx_unit=:ok,
        data_key="f_denoised")
    f = data_dict[data_key][get_unit_idx(data_dict, idx_unit), :]

    fitted, bleach_curve = fit_bleach(f, p0, plot_fit)

    data_dict["bleach_param"] = fitted.param
    data_dict["bleach_resid"] = fitted.resid
    data_dict["bleach_curve"] = bleach_curve
    data_dict["f_bleach"] = Array((data_dict["f_denoised"]' ./ bleach_curve)')

    nothing
end

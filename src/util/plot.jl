function highlight_stim(idx_stim, α_highlight=0.1)
    ax = gca()
    for i = 1:size(idx_stim, 1)
        ax.axvspan(idx_stim[i,1], idx_stim[i,2],
            alpha=α_highlight, facecolor="red", lw=0)
    end

    nothing
end

function highlight_stim(x, y, z, stim, cmap)
    stim = Float64.(stim)
    ax = gca()
    ax.scatter(x, y, z, c=cmap(stim / maximum(stim)), edgecolor="none")

    nothing
end

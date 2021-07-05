@recipe function f(df::DataFrame,
                   titlefont  = font(8, "Courier"),
                   xtickfont  = font(3, "Courier"),
                   ytickfont  = font(3, "Courier"),
                   yguidefont = font(5, "Courier"),
                   xguidefont = font(5, "Courier"),
                   legendfont = font(5, "Courier"),
                   markershape = :circle,
                   markersize = 2,
                   markerstrokewidth = 0,
                   legend = nothing)

    # Default values
    titlefont  := titlefont
    xtickfont  := xtickfont
    ytickfont  := ytickfont
    yguidefont := yguidefont
    xguidefont := xguidefont
    legendfont := legendfont
    markershape := markershape
    markersize := markersize
    markerstrokewidth := markerstrokewidth 
    
    dfts = tots(df)
    names_x = names(dfts)[2:end]
    x = Array(dfts[:,2:end])
    ts = dfts[:,1]
    name_ts = names(dfts)[1]

    legend := isnothing(legend) ? (size(x,2) <= 4 ? :outertop : :outerright) : legend

    @series begin
        label := permutedims(names_x)
        color_palette := :tab10
        seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
        seriestype := :line

        ts,x 
    end
end 

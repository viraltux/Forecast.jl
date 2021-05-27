@recipe function f(fc::FORECAST;
                   titlefont  = font(8, "Courier"),
                   xtickfont  = font(3, "Courier"),
                   ytickfont  = font(3, "Courier"),
                   yguidefont = font(5, "Courier"),
                   xguidefont = font(5, "Courier"),
                   legendfont = font(5, "Courier"),
                   xlims = nothing,
                   ylims = nothing,
                   markershape = :circle,
                   markersize = 2,
                   markerstrokewidth = 0,
                   title = fc.call,
                   xguide = "Time",
                   legend = nothing)

    # Default values
    titlefont  := titlefont
    xtickfont  := xtickfont
    ytickfont  := ytickfont
    yguidefont := yguidefont
    xguidefont := xguidefont
    legendfont := legendfont
    title      := title
    xguide     := xguide
    markershape := markershape
    markersize := markersize
    markerstrokewidth := markerstrokewidth 

    # Making the forecast to occupy the last 1/4 of xlims
    sx = size(fc.model.x,1)
    sm = size(fc.mean,1)
    ts = vcat(fc.model.x[:,1],fc.mean[:,1])
    xr1 = max(1,sx-3*size(fc.mean,1))
    xr2 = sx+sm
    xr = eltype(ts) in [Date, DateTime, Time] ?
        (ts[xr1], ts[xr2]) : (Dates.value(ts[xr1]), Dates.value(ts[xr2]))
    xlims := isnothing(xlims) ? xr : xlims
    xend = Array(fc.model.x[xr1:end,2:end])
    my = minimum([minimum(xend),minimum(Array(fc.lower[:,2:end]))])
    My = maximum([maximum(xend),maximum(Array(fc.upper[:,2:end]))])
    ypad = 0.1 * maximum([abs(My),abs(my)])
    ylims := isnothing(ylims) ? (my-ypad,My+ypad) : ylims

    xts = fc.model.x
    names_x = names(xts)[2:end]
    x = Array(xts[:,2:end])
    ts = xts[:,1]
    name_ts = names(xts)[1]

    legend := isnothing(legend) ? (size(x,2) <= 4 ? :top : :outerright) : legend

    @series begin
        label := permutedims(fc.model.varnames)
        color_palette := :tab10
        seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
        seriestype := :line

        ts,x
    end

    alpha_dimming = 2
    if size(x,2) == 1
        @series begin
            label := nothing
            color_palette := :tab10
            seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
            linestyle := :dot
            seriesalpha := 0
            fillalpha := 1/(alpha_dimming*size(x,2))

            rmean = Array(fc.mean[:,2:end])
            lower1 = Array(fc.lower[:,2:div(end-1,2)+1])
            upper1 = Array(fc.upper[:,2:div(end-1,2)+1])

            ribbon :=  (rmean - lower1, upper1 - rmean)

            fc.mean[:,1], Array(fc.mean[:,2:end])
            
        end
    end

    @series begin
        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
        seriestype := :line
        linestyle := :dot
        seriesalpha := 0
        fillalpha := 1/(alpha_dimming*size(x,2))
        
        rmean = Array(fc.mean[:,2:end])
        lower2 = Array(fc.lower[:,div(end-1,2)+2:end])
        upper2 = Array(fc.upper[:,div(end-1,2)+2:end])

        ribbon :=  (rmean - lower2, upper2 - rmean)
        
        fc.mean[:,1], Array(fc.mean[:,2:end])
    end

    @series begin

        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:size(x,2)),1,size(x,2))
        seriestype := :line
        linestyle := :dot
        seriesalpha := 1

        fc.mean[:,1], Array(fc.mean[:,2:end])
    end

end 

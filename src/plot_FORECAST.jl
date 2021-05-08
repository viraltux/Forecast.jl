@recipe function f(fc::FORECAST)

    # plot configuration
    titlefont:= font(8, "Courier")
    xtickfont:= font(3, "Courier")
    ytickfont:= font(3, "Courier")
    yguidefont:= font(5, "Courier")
    xguidefont:= font(5, "Courier")
    legendfont:= font(5, "Courier")

    x = fc.model.x isa TimeArray ? values(fc.model.x) : fc.model.x

    yguide := fc.model.call
    xguide := "Time"

    sox = ndims(x) <= 1 ? (length(fc.model.x),1) : size(fc.model.x)
    sx = ndims(x) <= 1 ? (length(x),1) : size(x)
    sm = ndims(fc.mean) <= 1 ? (length(fc.mean),1) : size(fc.mean)

    legend := sx[2] <= 4 ? :outertop : :outerright
    

    xr = (max(1,sx[1]-3*size(fc.mean,1)),sx[1] + sm[1])
    xend = x[xr[1]:end,:]
    xlims := xr

    ypad = 0.1
    my = minimum([minimum(xend),minimum(fc.lower)])
    My = maximum([maximum(xend),maximum(fc.upper)])
    ylims := (my-ypad*abs(my),My+ypad*abs(My))
    
    @series begin
        label := permutedims(fc.model.varnames)
        color_palette := :tab10
        seriescolor := reshape(collect(1:sx[2]),1,sx[2])
        seriestype := :line
        x
    end

    alpha_dimming = 2
    if size(x,2) ==1
        @series begin
            label := nothing
            color_palette := :tab10
            seriescolor := reshape(collect(1:sm[2]),1,sm[2])
            linestyle := :dot
            seriesalpha := 0
            fillalpha := 1/(alpha_dimming*sx[2])

            rmean = fc.mean
            lower1 = fc.lower[:,1:Int(end/2)]
            upper1 = fc.upper[:,1:Int(end/2)]

            ribbon :=  (rmean - lower1, upper1 - rmean)
            sx[1]+1:sx[1]+sm[1], fc.mean
        end
    end

    @series begin
        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:sm[2]),1,sm[2])
        seriestype := :line
        linestyle := :dot
        seriesalpha := 0
        fillalpha := 1/(alpha_dimming*sx[2])
        
        rmean = fc.mean
        lower2 = fc.lower[:,Int(end/2)+1:end]
        upper2 = fc.upper[:,Int(end/2)+1:end]

        ribbon :=  (rmean - lower2, upper2 - rmean)
        sx[1]+1:sx[1]+sm[1], fc.mean
    end

    @series begin

        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:sm[2]),1,sm[2])
        seriestype := :line
        linestyle := :dot
        seriesalpha := 1
        
        rmean = fc.mean

        sx[1]+1:sx[1]+sm[1], fc.mean
    end

    
end 

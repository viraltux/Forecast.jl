@recipe function f(fc::FORECAST)

    # plot configuration
    legend:= :outerright
    
    titlefont:= font(8, "Courier")
    xtickfont:= font(3, "Courier")
    ytickfont:= font(3, "Courier")
    yguidefont:= font(5, "Courier")
    xguidefont:= font(5, "Courier")

    yguide := fc.model.call
    xguide := "Time"

    ypad = 0.1
    xlims := (0,size(fc.model.x,1)[1] + size(fc.mean,1)[1])
    my = minimum([minimum(fc.model.x),minimum(fc.lower)])
    My = maximum([maximum(fc.model.x),maximum(fc.upper)])
    ylims := (my-ypad*my,My+ypad*My)
              

    sx = ndims(fc.model.x) <= 1 ? (length(fc.model.x),1) : size(fc.model.x)
    sm = ndims(fc.mean) <= 1 ? (length(fc.mean),1) : size(fc.mean)


    @series begin
        label := permutedims(fc.model.varnames)
        color_palette := :tab10
        seriescolor := reshape(collect(1:sx[2]),1,sx[2])
        seriestype := :line
        ymirror := false
        guide_position := :left
        fc.model.x
    end

    @series begin
        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:sm[2]),1,sm[2])
        linestyle := :dot
        ymirror := false
        guide_position := :left

        rmean = fc.mean
        lower1 = fc.lower[:,1:Int(end/2)]
        upper1 = fc.upper[:,1:Int(end/2)]

        ribbon :=  (rmean - lower1, upper1 - rmean)
        sx[1]+1:sx[1]+sm[1], fc.mean
    end

    @series begin
        label := nothing
        color_palette := :tab10
        seriescolor := reshape(collect(1:sm[2]),1,sm[2])
        seriestype := :line
        linestyle := :dot
        ymirror := false
        guide_position := :left

        rmean = fc.mean
        lower2 = fc.lower[:,Int(end/2)+1:end]
        upper2 = fc.upper[:,Int(end/2)+1:end]

        ribbon :=  (rmean - lower2, upper2 - rmean)
        sx[1]+1:sx[1]+sm[1], fc.mean
    end
end 

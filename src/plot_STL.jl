@recipe function f(fSTL::STL)

    # subplots configuration
    legend:= false
    grid:= false
    layout := @layout [Data
                       Seasonal
                       Trend
                       Remainder]
    
    xtickfont:= font(3, "Courier")
    ytickfont:= font(3, "Courier")
    yguidefont:= font(5, "Courier")

    # reference bar
    a,b = extrema(skipmissing(values(fSTL.decomposition[!,:Seasonal].+
                                     fSTL.decomposition[!,:Trend].+
                                     fSTL.decomposition[!,:Remainder]))); hd = b-a
    a,b = extrema(skipmissing(values(fSTL.decomposition[!,:Seasonal]))); hs = b-a
    a,b = extrema(skipmissing(values(fSTL.decomposition[!,:Trend]))); ht = b-a
    a,b = extrema(skipmissing(values(fSTL.decomposition[!,:Remainder]))); hr = b-a
    mh = min(hd,hs,ht,hr)

    inset_subplots := [(1, bbox(-0.012, 0, 0.01, 1.0, :left)),
                       (2, bbox(-0.012, 0, 0.01, 1.0, :right)),
                       (3, bbox(-0.012, 0, 0.01, 1.0, :left)),
                       (4, bbox(-0.012, 0, 0.01, 1.0, :right))]

    @series begin
        subplot := 1
        yguide := "Data"
        xaxis := nothing
        bottom_margin := -5Plots.mm    
        fSTL.decomposition[!,:Seasonal] .+
            fSTL.decomposition[!,:Trend] .+
            fSTL.decomposition[!,:Remainder]
    end

    @series begin
        subplot := 2
        yguide := "Trend"
        xaxis := nothing
        ymirror:=true
        guide_position:=:left
        bottom_margin := -5Plots.mm    
        fSTL.decomposition[!,:Trend]
    end

    @series begin
        subplot := 3
        yguide := "Seasonal"
        xaxis := nothing
        seriestype := :sticks
        bottom_margin := -5Plots.mm    
        fSTL.decomposition[!,:Seasonal]
    end

    @series begin
        subplot := 4
        yguide := "Remainder"
        seriestype := :sticks
        #bottom_margin:=0Plots.mm    
        ymirror:=true
        guide_position:=:left

        fSTL.decomposition[!,:Timestamp], fSTL.decomposition[!,:Remainder]
    end

    @series begin
        subplot := 5
        background_color_inside := nothing
        framestyle := :none
        seriestype := :bar
        ylims := (0,hd)
        mh:mh
    end

    @series begin
        subplot := 6
        background_color_inside := nothing
        framestyle := :none
        seriestype := :bar
        ylims := (0,ht)
        mh:mh
    end
    
    @series begin
        subplot := 7
        framestyle := :none
        seriestype := :bar
        ylims := (0,hs)
        mh:mh
    end

    @series begin
        subplot := 8
        framestyle := :none
        seriestype := :bar
        ylims := (0,hr)
        mh:mh
    end

end

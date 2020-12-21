@recipe function f(fSTL::Main.Forecast.STL)


    # subplots configuration
    legend --> false
    # link := :both
    grid --> false
    #framestyle := [:yaxis :yaxis :yaxis :axes]
    layout := @layout [Data{0.99w}  ebar
                       Seasonal     ebar
                       Trend        ebar
                       Remainder    ebar]

    xtickfont:= font(3, "Courier")
    ytickfont:= font(3, "Courier")
    yguidefont:= font(5, "Courier")

    # reference bar
    a,b = extrema(skipmissing(values(fSTL.ta[:Seasonal].+
                                     fSTL.ta[:Trend].+
                                     fSTL.ta[:Remainder]))); hd = b-a
    a,b = extrema(skipmissing(values(fSTL.ta[:Seasonal]))); hs = b-a
    a,b = extrema(skipmissing(values(fSTL.ta[:Trend]))); ht = b-a
    a,b = extrema(skipmissing(values(fSTL.ta[:Remainder]))); hr = b-a
    mh = min(hd,hs,ht,hr)

    @series begin
        subplot := 1
        yguide := "Data"
        xaxis := nothing
        fSTL.ta[:Seasonal] .+ fsts.ta[:Trend] .+ fsts.ta[:Remainder]
    end

    @series begin
        subplot := 2
        framestyle := :none
        seriestype := :bar
        ylims := (0,hd)
        mh:mh
    end

    @series begin
        subplot := 3
        yguide := "Trend"
        xaxis := nothing
        fSTL.ta[:Trend]
    end

    @series begin
        subplot := 4
        framestyle := :none
        seriestype := :bar
        ylims := (0,ht)
        mh:mh
    end

    
    @series begin
        subplot := 5
        yguide := "Seasonal"
        xaxis := nothing
        seriestype := :sticks
        fSTL.ta[:Seasonal]
    end

    @series begin
        subplot := 6
        framestyle := :none
        seriestype := :bar
        ylims := (0,hs)
        mh:mh
    end


    @series begin
        subplot := 7
        yguide := "Remainder"
        seriestype := :sticks
        fSTL.ta[:Remainder]
    end

    @series begin
        subplot := 8
        framestyle := :none
        seriestype := :bar
        ylims := (0,hr)
        mh:mh
    end

    
end

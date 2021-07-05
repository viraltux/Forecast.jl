@recipe function f(fCCF::CCF)

    # plot configuration
    legend:= false

    titlefont:= font(8, "Courier")
    xtickfont:= font(3, "Courier")
    ytickfont:= font(3, "Courier")
    yguidefont:= font(5, "Courier")
    xguidefont:= font(5, "Courier")

    N = fCCF.N
    if split(fCCF.type,'_')[1] == "pacf"
        type = (fCCF.type) == "pacf_step" ? "Stepwise PACF" : "Real PACF"
        type = (fCCF.type) == "pacf_stepwise-real" ? "Stepwise and Real PACF   |s|r|" : type
        auto_prefix = ""
    else
        type = (fCCF.type == "cor") ? "Correlation" : "Covariance"
        auto_prefix = fCCF.auto ? "Auto-" : "Cross-"
    end

    yguide := auto_prefix*type
    xguide := "Lag"
    
    # Center ticks
    lr = length(fCCF.ccf)
    gap = div(lr,11)
    mp = div(lr+1,2)
    rs = mp:gap:lr
    ls = mp-gap:-gap:1

    if fCCF.auto
        xticks := 1:gap:lr
    else
        xt = vcat(reverse(collect(ls)),collect(rs))
        xticks := (xt, xt .- xt[div(end+1,2)])
    end
    
    @series begin
        seriestype := :sticks
        ymirror := false
        guide_position := :left
        linewidth := 100/(fCCF.lag+1)
        
        collect(1:size(fCCF.ccf)[1]+1) .- ((ndims(fCCF.ccf) == 2) ? .1 : 0),
        ndims(fCCF.ccf) == 2 ? fCCF.ccf[:,1] : fCCF.ccf
    end

    if ndims(fCCF.ccf) == 2
        @series begin
            seriestype := :sticks
            ymirror := false
            guide_position := :left
            linewidth := 100/(fCCF.lag+1)

            collect(1:size(fCCF.ccf)[1]+1) .+ .1,
            fCCF.ccf[:,2]
        end
    end
    
    if (fCCF.type == "cor") | (split(fCCF.type,'_')[1] == "pacf")
        lg = fCCF.lag
        a1 = fCCF.alpha[1]; c1 = fCCF.ci[1]
        a2 = fCCF.alpha[2]; c2 = fCCF.ci[2]
        dc = c2-c1

        @series begin
            seriestype := :hline
            seriescolor := :black
            linealpha := 0.5
            linestyle := :dash
            [-c1,c1]
        end

        @series begin
            seriestype := :hline
            seriescolor := :black
            linealpha := 0.5
            linestyle := :dot
            ax = fCCF.auto ? lg : 2*lg + 1
            annotations := 
                [(ax,c1+dc/8,string("CI %",a1*100),font(3, "Courier")),
                 (ax,c2+dc/8,string("CI %",a2*100),font(3, "Courier"))]
            [-c2,c2]
        end
    end
end 

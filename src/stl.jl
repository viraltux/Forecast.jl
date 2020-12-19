module Forecast

export stl

include("sma.jl")
include("loess.jl")

using TimeSeries

"""
    sodd(x)

Return the smallest odd integer greater than or equal to `x`.        
"""
function sodd(x::Number)::Integer
    cx = Integer(ceil(x))
    mod(cx,2)==0 ? cx+1 : cx
end

"""
    stl(Yv, np; robust=false, 
                nl=sodd(np), 
                ns=10*length(Yv)+1,
                nt=sodd(1.5*np/(1-1.5/ns)), 
                ni=robust ? 1 : 2,
                no=robust ? 15 : 0,
                spm=false,
                qsmp=max(div(np,7),2))
    
Decompose a time series into trend, seasonal, and remainder components.

    "STL has s simple design that consists of a sequence of applications of the loess smoother; the simplicity allows analysis of the properties of the procedure and allows fast computation, even for very long time series and large amounts of trend and seasonal smoothing. Other features of STL  are specification of amounts of seasonal and trend smoothing that range, in a nearly continous way, from very small amount of smoothing to a very large amount; robust estimates of the trend and seasonal components that are not distorted by aberrant behavior in the data; specification of the period of the seasonal component to any intenger multiple of the time sampling interval greater than one; and the ability to decompose time series with missing values."

Excerpt from "STL: A Seasonal, Trend Decomposition Procedure Based on Loess"
              Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning.
              Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

    Args:
        `np`: Sesonality.
        `robust`: Robust stl.
        `nl`: Smoothing parameter of the low-pass filter.
        `ns`: Smoothing parameter for the seasonal component. It is chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; must always be odd and at least 7. The default value is not advised on the original paper but it is  the same chosen by the stl implementation in R.
        `nt`: Smoothing parameter for the trend decomposition.
        `ni`: Number of inner loop cycles.
        `no`: Number of outer loop cycles. The paper advises 5 ("safe value") or 10 ("near certainty of convergence") cycles when robustness is required, however the default is set at 15 following R implementation.
        `spm`: Seasonal post-smoothing.
        `qsmp`: Loess q window for Seasonal post-smoothing.
    Returns:
        An `stl` object with the seasonal, trend and remainder components if Yv is an Array and
        a TimeSeries object with the same components if Yv is a TimeSeries.

# Examples
```julia-repl
julia> Forecast.stl(co2_ts, 365; robust=true, spm=true)
[ Info: Seasonal and Trend corvengence achieved
4609×3 TimeArray{Union{Missing, Float64},2,Date,Array{Union{Missing, Float64},2}} 1974-05-17 to 1986-12-31
│            │ Seasonal │ Trend    │ Remainder │
├────────────┼──────────┼──────────┼───────────┤
│ 1974-05-17 │ 3.4255   │ 329.9863 │ -0.0318   │
│ 1974-05-18 │ 3.3813   │ 329.9891 │ -0.2604   │
│ 1974-05-19 │ 3.3359   │ 329.9919 │ 0.1322    │
│ 1974-05-20 │ 3.2894   │ 329.9947 │ 0.3559    │
  ...
```
"""


function stl(Yv::TimeArray,np::Integer;
             robust=false,
             nl=sodd(np),
             ns=10*length(Yv)+1,
             nt=sodd(1.5*np/(1-1.5/ns)),
             ni=robust ? 1 : 2,
             no=robust ? 15 : 0,
             spm=false,
             qsmp=max(div(np,7),2))

    str = stl(TimeSeries.values(Yv),np;
              robust=robust, nl=nl, ns=ns, nt=nt,
              ni=ni, no=no, spm=spm, qsmp=qsmp)
    stl_ts = (Timestamp = TimeSeries.timestamp(Yv),
              Seasonal = str.seasonal,
              Trend = str.trend,
              Remainder = str.remainder)
    TimeArray(stl_ts; timestamp = :Timestamp)
end

function stl(Yv, #::AbstractVector{T},
             np::Integer;
             robust=false,
             nl=sodd(np),
             ns=10*length(Yv)+1,
             nt=sodd(1.5*np/(1-1.5/ns)),
             ni=robust ? 1 : 2,
             no=robust ? 15 : 0,
             spm=false,
             qsmp=max(div(np,7),2)) #where T<:Union{Missing, Number}

    @assert mod(ns,2)==1 & (ns>=7) "ns is chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; must always be odd and at least 7"

    function B(u)
        u < 1 ? (1.0-u^2)^2 : 0.0
    end

    N = length(Yv)
    # initial robustness weights
    rhov = ones(N)
    # intial trend
    Tv = Tv0 = zeros(N)
    Sv = Sv0 = Array{Float64,1}(undef,N)
    Rv = Array{Float64,1}(undef,N)
    Cv = Array{Float64,1}(undef,N+2*np)
    cnv = false # convergence flag
    for o in 0:no
        for k in 1:ni
            # Updating sesonal and trend components
            ## 1. Detrending (Yv = Tv + Sv)
            Sv = Yv - Tv

            ### Seasonal convergence criterion
            Md = maximum(abs.(skipmissing(Sv-Sv0)))
            M0 = maximum(skipmissing(Sv0))
            m0 = minimum(skipmissing(Sv0))
            cnv = (Md/(M0-m0) < 0.01)
            println(Md/(M0-m0))
            Sv0 = Sv

            # Sesonal Smoothing 2,3,4
            ## 2. Cycle-subseries Smoothing
            for csi in 1:np
                Cv[csi:np:N+2*np] = loess(csi:np:N,
                                          Sv[csi:np:N];
                                          q=ns,d=1,rho=rhov[csi:np:N],
                                          predict=(1.0*csi-np):np:(N+np))
            end
            ## 3. Low-Pass Filtering of Smoothed Cycle-Subseries
            ### centered support instead 1:N to balance out machine error
            Lv = loess(1-ceil(N/2):N-ceil(N/2),
                       collect(skipmissing(sma(sma(sma(Cv,np),np),3))),
                       d=1,q=nl,rho=rhov)
            ## 4. Detreending of Smoothed Cycle-Subseries

            ### Lv is substracted to prevent low-frenquency power
            ### from entering the seasonal component.
            Sv = Cv[np+1:end-np] - Lv

            ## 5. Deseasonalizing
            Dv = Yv - Sv

            # Trend Smoothing
            ## 6. Trend Smoothing
            ### centered support instead 1:N to balance out machine error
            ### (floor isntead ceil like in Lv to balance out even lengths)
            Tv = loess(1-floor(N/2):N-floor(N/2),Dv,q=nt,d=1,rho=rhov)

            ### Trend convergence criterion
            Md = maximum(abs.(skipmissing(Tv-Tv0)))
            M0 = maximum(skipmissing(Tv0))
            m0 = minimum(skipmissing(Tv0))
            cnv = cnv & (Md/(M0-m0) < 0.01)
            Tv0 = Tv

        end
        # Computation of robustness weights
        ## These new weights will reduce the influence of transient behaviour
        ## on the trend and seasonal components.

        # Tv and Sv defined everywhere but Rv is not defined
        # where Yv has missing values
        Rv = Yv - Tv - Sv
        if 0 < o <= no
            smRv = skipmissing(Rv)
            h = 6*median(abs.(smRv))
            rhov = B.(abs.(smRv)/h)
        end
    end

    if spm
        # Using div(np,7) as default to approximate
        # the q=51 for a np=365 chosen in the original paper
        Sv = loess(1-ceil(N/2):N-ceil(N/2),Sv,q=qsmp)
        Rv = Yv - Tv - Sv
    end
    
    if cnv
        @info "Seasonal and Trend corvengence achieved"
    else
        @warn "Seasonal and/or Trend convergence not achieved, consider increasing the number of inner and/or outer cycles"
    end

     (seasonal = Sv,
      trend = Tv,
      remainder = Rv)
end


end

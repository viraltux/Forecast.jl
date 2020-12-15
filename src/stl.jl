module Forecast

export stl

include("sma.jl")
include("loess.jl")

"""
    sodd(x)

Return the smallest odd integer greater than or equal to `x`.
    
# Examples
```julia-repl
julia> sood(3.6)
5
```
"""
function sodd(x::Number)::Integer
    cx = Integer(ceil(x))
    mod(cx,2)==0 ? cx+1 : cx
end

"""
    stl(Yv, np; robust=false, nl=sodd(np), ns=10*length(Yv)+1,
                nt=sodd(1.5*np/(1-1.5/ns)), 
                ni=robust ? 1 : 2,
                no=robust ? 15 : 0,
                spm=false)
    
Decompose a time series into trend, seasonal, and remainder components.

"STL has s simple desing that consists of a sequence of applications of the loess smoother; the simplicity allows analysis of the properties of the procedure and allows fast computation, even for very long time series and large amounts of trend and seasonal smoothing. Other features of STL  are specification of amounts of seasonal and trend smoothing that range, in a nearly continous way, from very small amount of smoothing to a very large amount; roubst estimatews of the trend and seasonal components that are not distorted by aberrant behavior in the data; specification of the period of the seasonal component to any intenger multiple of the time sampling interval greater than one; and the ability to decompose time series with missing values."

Excerpt from "STL: A Seasonal, Trend Decomposition Procedure Based on Loess
Robert B. Cleveland, William S. Cleveland, Jean E. McRae, and Irma Terpenning" published in the Journal of Official Statistics Vol. 6. No. 1, 1990, pp. 3-73 (c) Statistics Sweden.

Args:
    `np`: Sesonality.
    `robust`: Robust stl.
    `nl`: Smoothing parameter of the low-pass filter.
    `ns`: Smoothing parameter for the seasonal component. It is chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; must always be odd and at least 7. The default value is not advised on the original paper but it is  the same chosen by the stl implementation in R.
    `nt`: Smoothing parameter for the trend decomposition.
    `ni`: Number of inner loop cycles.
    `no`: Number of outer loop cycles. The paper advises 5 ("safe value") or 10 ("near certainty of convergence") cycles when robustness is required, however the default is set at 15 following R implementation.
    `spm`: Seasonal post-smoothing.

Returns:
    An `stl` object.
"""

function stl(Yv::,
             np;
             robust=false,
             nl=sodd(np),
             ns=10*length(Yv)+1,
             nt=sodd(1.5*np/(1-1.5/ns)),
             ni=robust ? 1 : 2,
             no=robust ? 15 : 0,
             spm=false)
    
    @assert mod(ns,2)==1 & (ns>=7) "ns is chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; must always be odd and at least 7"

    #TODO implement seasonal post-smoothing option?
    #TODO implement multiplicative decomposition Fv?
    #TODO implement convergence criteria? at least return its value
    # with a recomendation if not converge takes place
    # max(abs(Sv-Sv1))/(max(Sv)-min(Sv)) < 0.01 then convergence takes place
    
    function B(u)
        u < 1 ? (1-u^2)^2 : 0
    end

    N = length(Yv)
    # initial robustness weights
    rhov = ones(N)
    # intial trend
    Tv = zeros(N)
    Sv = Array{Float64,1}(undef,N)
    Rv = Array{Float64,1}(undef,N)
    for o in 1:(no+1)
        for k in 1:ni
            # Updating sesonal and trend components
            ## 1. Detrending (Yv = Tv + Sv)
            Sv = Yv - Tv

            # Sesonal Smoothing 2,3,4
            ## 2. Cycle-subseries Smoothing
            Cv = zeros(N+2*np)
            for csi in 1:np
                Cv[csi:np:N+2*np] = loess(csi:np:N,Sv[csi:np:N];
                                          q=ns,d=1,k=rhov[csi:np:N],
                                          predict=(1.0*csi-np):np:(N+np))
            end
            ## 3. Low-Pass Filtering of Smoothed Cycle-Subseries
            Lv = loess(1.:N,collect(skipmissing(sma(sma(sma(Cv,np),np),3))),d=1,q=nl)
            ## 4. Detreending of Smoothed Cycle-Subseries

            ### Lv is substracted to prevent low-frenquency power
            ### from entering the seasonal component.
            Sv = Cv[np+1:end-np] - Lv
            
            ## 5. Deseasonalizing
            Dv = Yv - Sv

            # Trend Smoothing
            ## 6. Trend Smoothing
            Tv = loess(1.0:N,Dv,q=nt,d=1,k=rhov) #Tvk1
        end
        # Computation of robustness weights
        ## These new weights will reduce the influence of transient behaviour
        ## on the trend and seasonal components.

        # Tv and Sv defined everywhere but Rv is not defined
        # where Yv has missing values
        Rv = Yv - Tv - Sv
        h = 6*median(abs.(Rv))
        rhov = B.(abs.(Rv)/h)
    end

    (Yv,Tv,Sv,Rv)
end

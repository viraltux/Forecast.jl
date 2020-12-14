include("sma.jl")
include("loess.jl")

function stl(Yv,np;no=1,ni=1,nl=365,nt=573,ns=35)
    # parameters
    # np #sesonality
    # no #outer loop cycles
    # ni #inner loop cycles
    # nl #smoothing parameter of the low-pass filter
    # nt #smoothing parameter for the trend decomposition
    # ns #smoothing parameter for the seasonal component

    
    function B(u)
        u < 1 ? (1-u^2)^2 : 0
    end

    N = length(Yv)
    # initial robustness weights
    rhov = ones(N)
    # intial trend
    Tv = zeros(N)
    for o in 1:no
        for k in 1:ni
            # Updating sesonal and trend components
            ## 1. Detrending (Yv = Tv + Sv)
            Sv = Yv - Tv

            # Sesonal Smoothing 2,3,4
            ## 2. Cycle-subseries Smoothing
            Cv = zeros(N+2*np)
            for csi in 1:np
                a = (csi-1)*(div(N,np)+2)+1
                b = csi*(div(N,np)+2)
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
            Dv = Yv - Svk1

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

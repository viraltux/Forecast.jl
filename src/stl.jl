include("sma.jl")
include("loess.jl")

function stl(Yv,np;
             robust=false,
             spm=false, #seasonal post-smoothing
             nl=(mod(ceil(np),2)==0) ? Integer(ceil(np)+1) : Integer(ceil(np)),
             ns=10*length(Yv)+1, 
             nt=(mod(ceil(1.5*np/(1-1.5/ns)),2)==0) ? Integer(ceil(1.5*np/(1-1.5/ns))+1) : Integer(ceil(1.5*np/(1-1.5/ns))),
             ni=robust ? 1 : 2,
             no=robust ? 15 : 0)

    # parameters
    # np #sesonality
    # no #outer loop cycles
    # ni #inner loop cycles
    # nl #smoothing parameter of the low-pass filter
    # nt #smoothing parameter for the trend decomposition
    # ns #smoothing parameter for the seasonal component
    # ns must be carefully tailored for each application
    
    #nodd(1.5*np/(1-1.5/ns)),
    function nodd(x)
        cx = ceil(x)
        (mod(cx,2)==0) ? cx+1 : cx
    end
    
    @assert mod(ns,2)==1 & (ns>=7) "ns is chosen of the basis of knowledge of the time series and on the basis of diagnostic methods; must always be odd and at least 7"

    
    # if robustness is not needed ni = 2 and no = 0 otherwise
    # ni = 1 and no = 5 (no = 10 near certainty convergence)
    # it can be coded a convergence criterion such an in section 3.3

    #TODO implement seasonal post-smoothing option
    #TODO implement multiplicative decomposition Fv
    #TODO implement convergence criteria or at least return its value
    # with a recomendation if not converge takes place
    # max(abs(Sv-Sv1))/(max(Sv)-min(Sv)) < 0.01 then convergence takes place

    #using R defaults for ns and no
    
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

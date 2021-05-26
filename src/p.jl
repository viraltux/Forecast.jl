"""
Package: Forecast

    function p(dx, x0)
               

Return reverse lagged differences of a given order for Vector, Array and DataFrame.

# Arguments
- `dx`: Array or DataFrame of data.
- `x0`: Initial constants the reverse difference. The default value represents an integration of order one and lag one with initial values at zero. The format for the initial values is Array{Real,3}(lag, variable, order)"

# Returns
Lagged differences Vector or Array of a given order.

# Examples
```julia-repl

# Order two with Lag two
julia> x = repeat(1:2,30);
julia> dx = d(x,2,2);
julia> x0 = zeros(2,1,2); # lag 2, 1 variable, order 1
julia> x0[1,:,:] = collect(1:2);
julia> p(dx,x0) ≈ x
true

# Calculation of π
julia> x = 0:0.001:1;
julia> y = sqrt.(1 .- x.^2);
julia> isapprox(4*p(y)[end]/1000 , π, atol = 0.01)
true
```
"""
function p(dx::AbstractArray,
           x0 = reshape(repeat([0],size(dx,2)),1,:,1))

    format_ol = "x0 format is Array{Real,3}(order, variable, lag)"
    @assert ndims(x0) <= 3 format_ol

    ndims(x0) == 0 && (x0 = reshape([x0],1,1,1))
    ndims(x0) == 1 && (x0 = reshape(x0,:,1,1))
    ndims(x0) == 2 && (x0 = reshape(x0,:,size(x0,2),1))

    or = size(x0,1)
    la = size(x0,3)
    n = size(dx,1)
    nv = size(dx,2)

    @assert nv == size(x0,2) format_ol
    @assert 0 <= la & la <= (n-1)
    
    if !isnothing(findfirst(ismissing, dx))
        @error "Missing values allowed at the start and end of `dx` but not within."
        if !isnothing(findfirst(isnan, dx))
            @error "NaN values not allowed; use `missing` at the start or end of `dx` but not within."
        end
    end

    px = similar(dx, nv == 1 ? (n+or*la,) : (n+or*la,nv))
    px[1:size(dx,1),:] = dx
    for i in or:-1:1
        for j in 1:la
            pxj  = px[j:la:end,:]
            npxj = size(pxj,1)
            px[j:la:end,:] = cumsum(vcat(x0[i:i,:,j], pxj), dims=1)[1:npxj,:]
        end
    end

    return px
    
end

function p(df::DataFrame,
           x0 = reshape(repeat([0],size(tots(df),2)-1),1,:,1))

    df = tots(df)
    ts = df[:,1]
    dfx = df[:,2:end]
    pnames = "p[" * string(size(x0)) * "]_" .* names(dfx)
    
    px = p(Array(dfx), x0)

    dfpx = DataFrame(reshape(px,size(px,1),size(px,2)),pnames)
    pdfts = DataFrame(vcat(nΔt(ts,-(size(px,1)-size(df,1))), ts), [names(df)[1]])
    
    return hcat(pdfts,dfpx)

end

function p(fc::FORECAST,
           x0::AbstractArray = [[0]])

    ts_x = tots(fc.model.x)
    names_x = names(ts_x)
    names_x[2:end] = "p$(size(x0))_" .* names_x[2:end]

    ts_mean = tots(fc.mean)
    ts_lower = tots(fc.lower)
    ts_upper = tots(fc.upper)
    names_mean = names(ts_mean)
    names_mean[2:end] = "p$(size(x0))_" .* names_mean[2:end]
    names_lower = names(ts_lower)
    names_lower[2:end] = "p$(size(x0))_" .* names_lower[2:end]
    names_upper = names(ts_upper)
    names_upper[2:end] = "p$(size(x0))_" .* names_upper[2:end]

    xts = ts_x[:,1]
    x = Array(ts_x[:,2:end])

    fmean = Array(ts_mean[:,2:end])
    flower = Array(ts_lower[:,2:end])
    fupper = Array(ts_upper[:,2:end])

    n = size(x,1)
    m = size(x,2)

    pfxmean = p(vcat(x,fmean), x0)
    pfx = convert(Matrix{Float64}, reshape(pfxmean,:,size(x,2))[1:end-size(fmean,1),:])
    pfmean = convert(Matrix{Float64}, reshape(pfxmean,:,size(x,2))[end-size(fmean,1)+1:end,:])

    fts_x = vcat(nΔt(xts,-(size(pfxmean,1)-size(x,1)-size(fmean,1))),xts)

    fts_mean = nΔt(xts,size(pfmean,1))

    pfc = deepcopy(fc)
    pfc.model.varnames = names_x[2:end]

    pfc.model.x = hcat(DataFrame(fts_x,[names_x[1]]),
                       DataFrame(pfx,names_x[2:end]))
    pfc.mean    = hcat(DataFrame(fts_mean[:,1:1],[names_mean[1]]),
                       DataFrame(pfmean,names_mean[2:end]))

    if size(fmean,2) > 1
        z = fupper .- repeat(fmean,1,2)
        pfmean = repeat(pfmean,1,2)
    else
        z = fupper .- fmean
    end
    
    pfc.upper = hcat(DataFrame(fts_mean[:,1:1],[names_upper[1]]),
                     DataFrame(pfmean .+ z,names_upper[2:end]))
    pfc.lower = hcat(DataFrame(fts_mean[:,1:1],[names_lower[1]]),
                     DataFrame(pfmean .- z,names_lower[2:end]))
    
    pfc.call = fc.call * "\nIntegrated with x0 size: " * string(size(x0))

    return(pfc)
    
end

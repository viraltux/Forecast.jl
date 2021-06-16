"""
Package: Forecast

Store results from the function stl

# Arguments
    `decomposition::DataFrame`    A time array with three time series from a fitted STL model
    `call::String`                method called to generate ta
""" 
mutable struct STL{T<:DataFrame}
    decomposition::T  # A DataFrame with three time series from a fitted STL model
    call::String      # method called t
end

function Base.show(io::IO, stlx::STL)
    printstyled("STL Object:",bold=true,color=:underline)
    print(" ", stlx.call)
end



"""
Package: Forecast

Store results from the function `summarize`

# Arguments
`quantiles::DataFrame`        DataFrame with the data quantiles for each column
`moment::DataFrame`           DataFrame with the first four moments for each column
`format::DataFrame`           DataFrame with number of types per column
"""
mutable struct SUMMARIZE
    quantiles::DataFrame
    moments::DataFrame
    format::DataFrame
end

function Base.show(io::IO, s::SUMMARIZE)
    pretty_table(s.quantiles, nosubheader = true, show_row_number=false)
    pretty_table(s.moments, nosubheader = true, show_row_number=false)
    pretty_table(s.format, nosubheader = true, show_row_number=false)
end


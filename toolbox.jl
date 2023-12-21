#---------------------------------------------------
# This file contains some self-defined functions
# that are used in the simulation
#---------------------------------------------------

function knnsearch(v1, v2)
    sorted_v1 = sort(v1)
    result = similar(v2, eltype(v1))
    for (i, val) in enumerate(v2)
        idx = searchsortedfirst(sorted_v1, val)
        if idx == 1
            result[i] = 1
        elseif idx == length(v1) + 1
            result[i] = length(v1)
        else
            dist1 = abs(val - sorted_v1[idx - 1])
            dist2 = abs(val - sorted_v1[idx])
            result[i] = dist1 < dist2 ? idx - 1 : idx
        end
    end
    return result
end

function index_to_value(A, B)
    return view(A, B)
end

function estimate_lambda1(x, FXP, emom)
    return (FXP.δ * (-1 + (1 + FXP.δ / x) * log(1 + x / FXP.δ)) - emom[1])^2
end 

nanmean(x)                      = mean(filter(!isnan,x));
nanmean(x,y)                    = mapslices(nanmean,x, dims=y);
nanstd(x)                       = std(filter(!isnan,x));
nanstd(x,y)                     = mapslices(nanstd,x, dims=y);
nanskw(x)                       = skewness(filter(!isnan,x));
nanskw(x,y)                     = mapslices(nanskw,x, dims=y);
nansum(x)                       = sum(filter(!isnan,x));
nansum(x,y)                     = mapslices(nansum,x,dims=y);
nanmax(x)                       = maximum(filter(!isnan,x));
nanmax(x,y)                     = mapslices(nanmax,x,dims=y);
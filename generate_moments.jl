#------------------------------------------------
# This file generates the moments to be matched
# given the simulated worker panel
#------------------------------------------------

function generate_moments(res)
    @unpack nw, nf, T = Set()

    # unemployment rate
    totemp                          = sum(res.e, dims=1);
    
    # transition rates
    j2j                             = sum(res.j2j[:,2:end],dims=1)./totemp[1:end-1];
    
    # mean log wage
    validw                          = res.w .> 0;
    logw                            = zeros(size(res.w));
    logw                           .= NaN;
    logw[validw]                   .= log.(res.w[validw]);
    mlw                             = nanmean(logw,1);

    # log wage variance
    sdlw                            = nanstd(logw,1);
    
    # log wage skewness
    sklw                            = nanskw(logw,1);

    # log wage growth
    sel                             = validw[:,2:end] .* validw[:,1:end-1];    # indicator of employment and nonneg wage in 2 consecutive periods
    dlogw                           = diff(logw, dims =2);
    dlw                             = nansum(sel .* dlogw, 1)./sum(sel, dims=1);
    sel                             = res.j2j[:,2:end] .* sel;
    dlwjj                           = nansum(sel .* dlogw, 1)./sum(sel, dims=1);
    
    # hiring (log) wages
    mhlw                            = nanmean(logw .* res.u2e, 1);
    sdhlw                           = nanstd(logw .* res.u2e, 1);
    
    # wage regressions
    lhs                             = logw[:];
    sel                             = validw[:];
    lhs                             = lhs[sel];
    rhs                             = hcat(ones(size(lhs)), res.ten[sel], res.eten[sel]);
    
    beta                            = (rhs' * rhs)\rhs' * lhs;
    b_ten                           = beta[2];
    b_eten                          = beta[3];
    
    # max log wage per worker
    nan_rows                        = findall(row -> all(isnan, row), eachrow(logw));
    logw_filtered                   = logw[setdiff(1:size(logw, 1), nan_rows), :];
    
    if size(logw_filtered, 1) == 0
        wrkmxlw = NaN * ones(48)
    else
        wrkmxlw = nanmax(logw_filtered, 2)
    end

    # mean and max log wage per firm
    jj                              = res.j[:];
    logw                            = logw[:];
    logw                            = logw[jj .!= 0];
    jj                              = jj[jj .!= 0];
    firms                           = unique(jj);
    frmmxlw                         = [maximum(logw[jj .== firm]) for firm in firms];
    frmmnlw                         = [mean(logw[jj .== firm]) for firm in firms];
  
    return Moments(mean(j2j), mean(mlw), mean(sdlw), mean(sklw), mean(dlw), mean(dlwjj), b_ten, b_eten, nanstd(wrkmxlw), nanstd(frmmxlw), nanstd(frmmnlw), mean(mhlw), mean(sdhlw))
end
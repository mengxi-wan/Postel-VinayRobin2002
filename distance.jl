function distance(px, emv, fxp)

    p         = par_convert(px, 2);
    FXP       = fxp;
    FXP.pj    = fxp.p_min .+ quantile(LogNormal(p.smpd[1], p.smpd[2]), fxp.Fj)
    FXP.iFb   = (FXP.pj .- FXP.p_min) .* (1 .- FXP.Fj) + exp(p.smpd[1] + p.smpd[2] * p.smpd[2]/2) * cdf(Normal(p.smpd[1] + p.smpd[2] * p.smpd[2], p.smpd[2]), log.(FXP.pj .- FXP.p_min))

    # simulate panel 
    sd        = simulate_panel(p, FXP)
    sm        = generate_moments(sd)

    # simulated moments
    smv       = [sm.mean_j2j; sm.mean_mlw; sm.mean_sdlw; sm.mean_sklw; sm.mean_dlw; sm.mean_dlwjj; sm.b_ten; sm.b_eten; sm.std_wrkmxlw; sm.std_frmmxlw; sm.std_frmmnlw; sm.mean_mhlw; sm.mean_sdhlw]

    # distance
    y         = smv - emv
    y         = y' * y
    
end

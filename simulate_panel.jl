#-----------------------------------------
# This function simulates a worker panel
# using the input parameters 
#-----------------------------------------

function simulate_panel(px, fxp)

    @unpack ρ, δ, λ1, λ0, b, p_min, Fj, pj, iFb = FixedParameters(fxp)
    @unpack sdeps, smpd = EstimatedParameters(px)
    @unpack nw, nf, T = Set()

    res = Initialize() 

    ####### simulate the first period #######
    # employment status
    urate                           = δ/(δ + λ0);
    dice                            = rand(nw,1);
    res.e[:,1]                      = (dice .>= urate);

    # draw initial productivity to the current firm and the productivity of last poaching firm
    dice                            = rand(nw,1); # draw from L(p)
    dice_Fp                         = (δ + λ1) * dice./(δ .+ λ1 * dice); # corresponding rank in F
    res.j[:,1]                      = knnsearch(Fj, dice_Fp);
    res.j[:,1]                      = res.e[:,1] .* res.j[:,1];
    dice                            = rand(nw,1); # draw from L(q|p)
    dice_Fq                         = (δ .+ λ1 .- (δ .+ λ1 * (1 .- dice_Fp))./sqrt.(dice))/λ1; # corresponding quantile in F
    res.jpoach[:,1]                 = knnsearch(Fj, dice_Fq);
    res.jpoach[:,1]                 = res.e[:,1] .* res.jpoach[:,1];

    # tenure: exponential distribution with intensity delta + lam1 * (1-dice_Fp)
    res.ten[:,1]                    = rand(nw,1);
    res.ten[:,1]                    = -log.(res.ten[:,1])./(δ .+ λ1 * (1 .- dice_Fp));
    res.ten[:,1]                    = floor.(res.ten[:,1] .* res.e[:,1]);

    # "employment" tenure (time since last unemployed)
    res.eten[:,1]                   = rand(nw,1);
    res.eten[:,1]                   = res.ten[:,1] - log.(res.eten[:,1])/δ;
    res.eten[:,1]                   = floor.(res.eten[:,1].*res.e[:,1]);

    ####### main simulation #######
    t = 2;
    while t <= T
        # event indicators:
        # indicator of unemployed job offer
        dice                        = rand(nw,1);
        u_offer                     = (dice .< λ0) .* (1 .- res.e[:,t-1]);
    
        # indicator of employed job offer or layoff
        dice                        = rand(nw,1);
        e_offer                     = (dice .< λ1) .* res.e[:,t-1];
        e_to_u                      = (dice .>= λ1) .* (dice .< λ1 + δ) .* res.e[:,t-1]; 
        e_to_u                      = e_to_u .> 0 
    
        # draw of the job quality of potential new offer and calculate potential surplus
        jq                          = zeros(nw,1);
        any_offer                   = (e_offer + u_offer .== 1);
        any_offer                   = (any_offer .> 0)
        dice                        = rand(sum(any_offer),1);
        jq[any_offer]               = knnsearch(Fj,dice);
    
        # update worker states:
        # unemployed hire
        u_to_e                      = u_offer .> 0;
        u_to_e                      = vec(u_to_e);
        res.e[u_to_e,t]            .= 1;
        res.j[u_to_e,t]            .= jq[u_to_e];
        res.jpoach[u_to_e,t]       .= 0;
        res.j2j[u_to_e,t]          .= 0;
        res.u2e[u_to_e,t]          .= 1;
        res.e2u[u_to_e,t]          .= 0;
        res.ten[u_to_e,t]          .= 0;
        res.eten[u_to_e,t]         .= 0;
    
        # employed layoff
        e_to_u                      = vec(e_to_u); 
        res.e[e_to_u,t]            .= 0;
        res.j[e_to_u,t]            .= 0;
        res.jpoach[e_to_u,t]       .= 0;
        res.j2j[e_to_u,t]          .= 0;
        res.u2e[e_to_u,t]          .= 0;
        res.e2u[e_to_u,t]          .= 1;
        res.ten[e_to_u,t]          .= 0;
        res.eten[e_to_u,t]         .= 0;
    
        # employed accepted offer
        e_to_e                      = e_offer .* (jq .> res.j[:,t-1]) .> 0;
        e_to_e                      = vec(e_to_e)
        res.e[e_to_e,t]            .= 1;
        res.j[e_to_e,t]             = jq[e_to_e];
        res.jpoach[e_to_e,t]        = res.j[e_to_e,t-1];
        res.j2j[e_to_e,t]          .= 1;
        res.u2e[e_to_e,t]          .= 0;
        res.e2u[e_to_e,t]          .= 0;
        res.ten[e_to_e,t]          .= 0;
        res.eten[e_to_e,t]          = res.eten[e_to_e,t-1] .+ 1;
    
        # employed renegotiation
        reneg                       = e_offer .* (jq .<= res.j[:,t-1]) .* (jq .> res.jpoach[:,t-1]) .> 0;
        reneg                       = vec(reneg)
        res.e[reneg,t]             .= 1;
        res.j[reneg,t]              = res.j[reneg,t-1];
        res.jpoach[reneg,t]         = jq[reneg];
        res.j2j[reneg,t]           .= 0;
        res.u2e[reneg,t]           .= 0;
        res.e2u[reneg,t]           .= 0;
        res.ten[reneg,t]            = res.ten[reneg,t-1] .+ 1;
        res.eten[reneg,t]           = res.eten[reneg,t-1] .+ 1;
    
        # employed stay put
        e_stay                      = res.e[:,t-1] .* (1 .- e_to_e) .* (1 .- e_to_u) .* (1 .- reneg) .> 0;
        e_stay                      = vec(e_stay);
        res.e[e_stay,t]            .= 1;
        res.j[e_stay,t]             = res.j[e_stay,t-1];
        res.jpoach[e_stay,t]        = res.jpoach[e_stay,t-1];
        res.j2j[e_stay,t]          .= 0;
        res.u2e[e_stay,t]          .= 0;
        res.e2u[e_stay,t]          .= 0;
        res.ten[e_stay,t]           = res.ten[e_stay,t-1] .+ 1;
        res.eten[e_stay,t]          = res.eten[e_stay,t-1] .+ 1;
     
        # move on to next period
        t = t+1;
    end

    # construct wages
    res.j[res.j .== 0]             .= 1;
    res.j[res.j .== 1001]          .= 1000;
    res.jpoach[res.jpoach .== 0]   .= 1;
    res.w                           = index_to_value(pj,res.jpoach) .- λ1 * (index_to_value(iFb,res.j) - index_to_value(iFb, res.jpoach))/(ρ + δ);
    res.w                           = res.w .* res.e; 
    res.j                           = res.j .* res.e; # change the job quality of unemployed workers back to zero
    res.jpoach                      = res.jpoach .* res.e; 

    # include epsilon
    ε                               = rand(LogNormal(-0.5*sdeps^2, sdeps), nw);
    res.w                           = res.w .* ε;

    # express tenure (and employment tenure) in years rather than months
    res.ten                         = res.ten/12;       
    res.eten                        = res.eten/12;

    return res
end


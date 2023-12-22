#--------------------------------------------------------
# Replication of Postel-Vinay and Robin (2002)
# Authors: Mengxi Wan
# Referenced matlab code from Postel-Vinay's website
# This is the main file that estimates the parameters
#--------------------------------------------------------

# Load packages and functions
using Parameters, Random, Statistics, Distances, NearestNeighbors, Optim, CSV, DataFrames, LineSearches, Distributions, NelderMead
include("par_convert.jl"); include("toolbox.jl"); include("simulate_panel.jl"); include("generate_moments.jl"); include("distance.jl")

# Load the empirical moments
emom = CSV.File("emom.csv", header=false) |> DataFrame;
emom = Vector(emom[:, 1]);

# Define the structs
@with_kw mutable struct FixedParameters         
    b::Float64                = 10                               
    j::Array{Float64, 1}      = collect(1:1:2000)
    ρ::Float64                = 0.008
    δ::Float64                = 0.025                            
    λ1::Float64               = 0.1                              
    λ0::Float64               = 0.3                                                                      
    p_min::Float64            = 10      
    Fj::Array{Float64, 1}     = collect(0:1/1999:1) 
    pj::Array{Float64, 1}     = collect(10:1/1999:12)
    iFb:: Array{Float64, 1}   = collect(0:1/1999:1) 
end

@with_kw mutable struct Moments
    mean_j2j::Float64
    mean_mlw::Float64
    mean_sdlw::Float64
    mean_sklw::Float64
    mean_dlw::Float64
    mean_dlwjj::Float64
    b_ten::Float64
    b_eten::Float64
    std_wrkmxlw::Float64
    std_frmmxlw::Float64
    std_frmmnlw::Float64
    mean_mhlw::Float64
    mean_sdhlw::Float64
end

@with_kw mutable struct Results
    e::Array{Float64, 2}
    j::Array{Int64, 2}
    jpoach::Array{Int64, 2}
    j2j::Array{Float64, 2}
    u2e::Array{Float64, 2}
    e2u::Array{Float64, 2}
    ten::Array{Float64, 2}
    eten::Array{Float64, 2}
    w::Array{Float64, 2}
end

@with_kw mutable struct Set
    seed                      = 88
    nw::Int64                 = 2000 
    T::Int64                  = 48                         
    nf::Int64                 = floor(nw/10)                    
end

@with_kw mutable struct EstimatedParameters   
    sdeps::Float64            
    smpd::Array{Float64, 1}   
end

function Initialize()
    @unpack nw, T             = Set()
    e                         = zeros(nw,T)
    j                         = zeros(nw,T)
    jpoach                    = zeros(nw,T)
    j2j                       = zeros(nw,T)
    u2e                       = zeros(nw,T)
    e2u                       = zeros(nw,T)
    ten                       = zeros(nw,T)
    eten                      = zeros(nw,T)
    w                         = zeros(nw,T)
    Results(e, j, jpoach, j2j, u2e, e2u, ten, eten, w)
end 

# λ1 can be estimated or fixed. Following are the codes to estimate λ1 
fxp                         = FixedParameters();
λ1                          = Optim.optimize(x -> estimate_lambda1(x[1], fxp, emom), [2 * emom[1]], method = LBFGS(), show_trace = true);
λ1                          = Optim.minimizer(λ1);
fxp.λ1                      = λ1[1];

# The productivity distribution is taken as given
Fj                          = CSV.File("productivity.csv", header=false) |> DataFrame;
fxp.Fj                      = Vector(Fj[:, 1]);

# We guess the unobserved worker heterogeneity accounts for 1/4 of wage variance
p0                          = EstimatedParameters(sqrt(emom[3]/4), [0,1]);
initial                     = [sqrt(p0.sdeps), p0.smpd[1], p0.smpd[2]];
initial                     = [0.5, 0.6, 1]

# Use Nelder-Mead to find the minimum of the distance between simulated and empirical moments
# NelderMead pacakge
result                      = NelderMead.optimise(x -> distance(x, emom, fxp), initial, ones(3) ./ 10)
pv_out, d_out, returncode, iters, simplex = result

# Optim package
# options                   = Optim.Options(show_trace = true, iterations = 100);
# result                    = Optim.optimize(x -> distance(x, emom, fxp), initial, NelderMead(), options)
# pv_out                    = Optim.minimizer(result)
# d_out                     = Optim.minimum(result)

# Calculate the simulated moments using the optimal parameters
ESP                         = par_convert(pv_out, 2);
SDT                         = simulate_panel(ESP, fxp);
SMO                         = generate_moments(SDT);
smom                        = [SMO.mean_j2j, SMO.mean_mlw, SMO.mean_sdlw, SMO.mean_sklw, SMO.mean_dlw, SMO.mean_dlwjj, SMO.b_ten, SMO.b_eten, SMO.std_wrkmxlw, SMO.std_frmmxlw, SMO.std_frmmnlw, SMO.mean_mhlw, SMO.mean_sdhlw];

# Show the model fit
println(emom')
println(smom')

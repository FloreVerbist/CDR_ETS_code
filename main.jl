#---------------------
# CODE INFO
# author: Flore Verbist 
# description: this code is used for the paper entitled "Carbon Removals Meet Emission Trading System Desings: A Precautionary Path towards Integration"

#----------------------
## MAIN FILE 
# Step 0: activating environment
using Pkg 
# cd("./Environment_CDR")
# pwd() # double check if you are in the right directory
# Pkg.activate(".") # activate a new environment
# Pkg.instantiate()
# Pkg.status() # double check if its activated

# STEP 1): data and package access
# A) input packages
using CSV
using DataFrames
using XLSX
using JuMP
using Gurobi
using Plots
using StatsPlots
using Statistics
using LaTeXStrings
using Plots; LaTeXStrings; pgfplotsx() # gr() --> pgfplotsx requires pdflatex viewer via Miktex or Texlive installer (make sure those programs are installed as well)
using Plots.PlotMeasures
using Vega
using NumericIO
using ProgressBars, Printf # progress bar
using StatsBase
using IterTools
using IterableTables
using QuadGK

# input general function files attached to the main file 
include("import_data.jl")                                         # processing of excel data 
include("case_processing.jl")                                     # initialising model 
include("results.jl")                                             # plots of results 

# Data   
# CAP
T_OPT = 2100-2024           # Time horizon 
T_OPT_plot = 2070           # Plotting horizon
Timeseries_data, Index_names, Parameters = import_data_CCUS("C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Scripts/Script_P2_CDR_ETS/Script_P2_CDR_ETS/Data_removals.xlsx")  # data reading 
Carbon_budget = sum(Timeseries_data[convert.(Int64,t+1), "ETS_cap [EUA/y]"] for t in 0:1:T_OPT )
Budget_range = 10000:2000:34000
EM_REF = EM_ref = 3000  # Initial emissions
DR = 0.06               # discountrate
DRc = 0.06
INFL = 0.0              # Inflation rate



# Marginal costs of Abatement 
# LINEAR 
BETA = 0.1 
GAMMA = 1 
# NON LINEAR 
BETA = 4.5*10^(-5)   #(-17)
GAMMA = 2

# Parameters for Marginal costs of Removals
# LINEAR 
A_r = 150 
B_r = 200
# NON LINEAR 
DISC_LOOP = 0.07    # endogenous part
U_R = 0.3 
M_R = 0.3 #0.3
D_R = 2500 #2500
K_R = - 0.1 # - 0.1
C_R = 0.88
Cum_remov_init = 0.0

m_R_add = 300   # exogenous part
r_R_add = 80.0


# Parameters for Marginal costs of abatement
# LINEAR 
A_a = 80.0  
B_a = 300.0
# NON LINEAR 
D_A = 100 #50 # endogenous part
C_A = 0.98
M_A = 0.1 # 0.05
K_A = -50
U_A = 1.0 
Cum_abat_init = 5000

m_A_add = 50.0 # exogenous part
r_A_add = 40.0



# ADMM parameters
ADMM_max_it = 800
RHO =  RHO2 = RHO1 =  10^(-2) 
RHO_vect = vcat(0.1, ones(9).*0.01, 0.03)
EAC_start = 0.0 #euro/tCO2
CDR_start = 0.0 
a = 1 
EPSILON = 1.0
n = 3


# Case specific parameters 
P_CAP = 300 #€/tCO2     # limit on removal integration
CAP_MIN = 300 #MtCO2    # limit of cap in the net-cap case
 
T1 = 20                 # Time at which case change
FRAC_eq =  0.85         # non-permanence fraction/obligation fraction
FRAC_R_vec = vcat(1.0*ones(7), FRAC_eq, 1.0)            # Fraction in reality 
# plotting information

Plot_title_ordered = [
    "Social Optimum", 
    "No Removals (ETS Endgame)", 
    "Disconnected Markets - Removal Target", 
    "Connected Markets - Net Cap", 
    "Net Cap with Removal Obligation (net-negative)", 
    "Net Cap with Equivalence Obligation (impermanence)",
    "Conditional Market Integration - Net Cap (price-based)",  
    "Connected Markets - Gross Cap", 
    "Gross Cap with Removal Target (net-zero)"]
Plot_title = [
    "Social Optimum", 
    "No Removals (ETS Endgame)", 
    "Disconnected Markets - Removal Target", 
    "Connected Markets - Net Cap", 
    "Connected Markets - Gross Cap", 
    "Conditional Market Integration - Net Cap (price-based)", 
    "Net Cap with Removal Obligation (net-negative)", 
    "Net Cap with Equivalence Obligation (impermanence)", 
    "Gross Cap with Removal Target (net-zero)"]
Plot_title = ["","","","","","","","",""]
CASE_NAME_ordered = [
    "C0-social_optimal",
    "C1-no_removals",
    "C2-disconnected", 
    "C3-full_integration_net_cap", 
    "C6-removal_obligation", 
    "C7-equivalence_obligation",
    "C5-conditional_integration",
    "C4-full_integration_gross_cap", 
    "C8-rem_obl_gross"
]
CASE_NAME = [
    "C0-social_optimal",
    "C1-no_removals",
    "C2-disconnected", 
    "C3-full_integration_net_cap", 
    "C4-full_integration_gross_cap", 
    "C5-conditional_integration",
    "C6-removal_obligation", 
    "C7-equivalence_obligation",
    "C8-rem_obl_gross"
]
Table_column_names = [
    "Total_net_emissions",
    "Total_emissions",
    "Fraction_removals",
    "Total_costs_per_net",
    "Total_costs_per_abat",
    "Total_budget_per_net",
    "Learning_gains_per_rem",
    "Learning_gains_per_em",
    "ETS_max",
    "ETS_avg",
    "RC_max",
    "ETS_price",
    "Fraction_costs_abat",
    "Fraction_costs_rem",
    "Emitter_cost_avg"
]

# Results pre-processing
Dict_results = Dict()
ETS_data = Dict() 
compile_dictionary(Dict_results, ETS_data)


ft1(x) = @.([10 (-2x+20)])
ft(x) = minimum(ft1(x); dims=2)
quadgk(var -> ft(var), 0, 10)[1] # this actually works correct integration 
# Preprocessing marginal cost curves 
LINEAR = false
LEARNING = true
if LINEAR == true 
    f(remov, t, Cum_remov) = @.(A_r + B_r*remov)
    #f2(abat, t, Cum_abat) = @.(BETA*(EM_REF*(abat)).^GAMMA)
    f2(abat, t, Cum_abat) = @.(A_a+ B_a*abat)
else 
# NOTE: chaning the parameters of this function, will automatically be taken into account
    if LEARNING== true
        f3(remov, t, Cum_remov) = @.[(1100 + remov*100) (D_R ./ (1 + Cum_remov[t-1]).^(M_R) .*(remov./((U_R/(1+exp(K_R*(t+1)))) .- remov)).^C_R + m_R_add*(1-t/T_OPT) + r_R_add)]
        f(remov, t, Cum_remov) = minimum(f3(remov, t,  Cum_remov); dims=2)
        fr(remov, t,  Cum_remov) = @.(D_R ./ (1 + Cum_remov[t-1]).^(M_R) .*(remov./((U_R/(1+exp(K_R*(t+1)))) .- remov)).^C_R + m_R_add*(1-t/T_OPT) + r_R_add)

        f1(abat, t,  Cum_abat) =  @.[(450+(abat)*150) (D_A ./ (1 +  Cum_abat[t-1]).^(M_A) .*( abat./((U_A/(1+exp(K_A*(t+1)))) .-  abat)).^C_A + m_A_add*(1-t/T_OPT) + r_A_add)]
        f2(abat, t,  Cum_abat) =  minimum(f1(abat, t,  Cum_abat); dims=2)
        fa(abat, t,  Cum_abat) = @.(D_A ./ (1 +  Cum_abat[t-1]).^(M_A) .*( abat./((U_A/(1+exp(K_A*(t+1)))) .-  abat)).^C_A + m_A_add*(1-t/T_OPT) + r_A_add)
        #g(abat, t, Cum_abat) = @.(minimum([400.0.+abat, @.(D_A ./ (1 + Cum_abat[t-1]).^(M_A) .*(abat./((U_A/(1+exp(K_A*(t)))) .- abat)).^C_A + F_A_add*(1-t/T_OPT))]))

    elseif LEARNING == false
        f(remov, t, Cum_remov) = @.(D_R ./ (1 + 12000).^(M_R) .*(remov./((U_R/(1+exp(K_R*(T_OPT/2+1)))) .- remov)).^C_R + m_R_add*(1-1/3) + r_R_add)
        f1(abat, t,  Cum_abat) =  @.[(290+(abat)*110) (D_A ./ (1 +  (Cum_abat_init+9000)).^(M_A) .*( abat./((U_A/(1+exp(K_A*(T_OPT/2+1)))) .-  abat)).^C_A + m_A_add*(1-1/3) + r_A_add)]
        f2(abat, t,  Cum_abat) =  minimum(f1(abat, t,  Cum_abat); dims=2)
        #g(abat, t, Cum_abat) = @.(minimum([400.0.+abat, @.(D_A ./ (1 + Cum_abat[t-1]).^(M_A) .*(abat./((U_A/(1+exp(K_A*(t)))) .- abat)).^C_A + F_A_add*(1-t/T_OPT))]))
    end
end
# preprocessing_marginal_cost_curves(LINEAR, LEARNING) # doesn't work yet, so use the lines above

# STEP 2): Running cases 

EPSILON = 1.0
DR = 0.06
for (C_nr, Case_name) in enumerate(CASE_NAME)
    Case_name = CASE_NAME[C_nr]
    RHO = RHO1 = RHO2 = RHO_vect[C_nr]
    FRAC_R = FRAC_R_vec[C_nr] # Reality on actual value of removals
    ###--------------------------------------------------------------
    results = Dict()
    ADMM = Dict()
    global BUDGET = Carbon_budget

    print("$(Case_name)  --  ")
    include("$(Case_name).jl")

    define_results!(Parameters,results,ADMM)

    model_e = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_r = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))
    model_ccb = Model(optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" =>  0))


    # Step 3b: process parameters
    process_sets_parameters!(model_e, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_r, Parameters,Timeseries_data, T_OPT)
    process_sets_parameters!(model_ccb, Parameters,Timeseries_data, T_OPT)


    # calculate equilibrium 
    ADMM["walltime"] = @elapsed Dict_results = ADMM_price_based_m!(Case_name,results, ADMM, model_e, model_r, model_ccb, Dict_results, ETS_data, LINEAR)  

    # Step: plotting
    plot_model!(Case_name, C_nr, model_e, model_r, model_ccb, ADMM, results, Dict_results,  true)
    ###--------------------------------------------------------------

end 

# STEP 3) a) Extracting and saving results to table
Dict_table = dict_table(Dict_results)
# b) Or loading the results from an existing table
# Dict_table = DataFrame(CSV.File("output_paper.csv")) 
#ordered_Cnr = [1, 2, 3, 4, 7, 8, 6, 5, 9]

# STEP 4) Making plots
#plot_pareto!(Dict_results,true) # use direct function --> doesnt work yet in this format - run directly from original script file
#learningcurves(results, ADMM)  
ccb_budget_plot(Dict_table)
emitter_budget_plot(Dict_table)
learning_plot(Dict_table)
carbon_lock_in_plot(Dict_table)
# deflation_effect()


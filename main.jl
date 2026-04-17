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

LINEAR = false
LEARNING = true

include("parameters.jl")
# Results pre-processing
Dict_results = Dict()
ETS_data = Dict() 
compile_dictionary(Dict_results, ETS_data)


# preprocessing_marginal_cost_curves(LINEAR, LEARNING) # doesn't work yet, so use the lines above

# STEP 2): Running cases 


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


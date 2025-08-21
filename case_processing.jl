
# RC(remov, t) = @. (D_r ./ (1 + Cum_remov[t-1]).^(M_r) .*(remov./((U_r/(1+exp(K_r*(t)))) .- remov)).^C_r)
# AC(abat, t) = @. (D_a ./ (1 + Cum_abat[t-1]).^(M_a) .*(abat./((U_a/(1+exp(K_a*(t+1)))) .- abat)).^C_a)

function compile_dictionary(Dict_results, ETS_data)


    for Case_name in CASE_NAME
        Dict_results[(Case_name,"Total_net_emissions")] = []
        Dict_results[(Case_name,"Total_emissions")] = []
        Dict_results[(Case_name,"Fraction_removals")] = []
        Dict_results[(Case_name,"Total_costs_per_net")] = []
        Dict_results[(Case_name,"Total_costs_per_abat")] = []
        Dict_results[(Case_name,"Total_budget_per_net")] = []
        Dict_results[(Case_name,"Learning_gains_per_rem")] = []
        Dict_results[(Case_name,"Learning_gains_per_em")] = []
        Dict_results[(Case_name,"ETS_max")] = []
        Dict_results[(Case_name,"ETS_avg")] = []
        Dict_results[(Case_name,"RC_max")] = []
        Dict_results[(Case_name,"ETS_price")] = []
        Dict_results[(Case_name,"Fraction_costs_abat")] = []
        Dict_results[(Case_name,"Fraction_costs_rem")]= []
        Dict_results[(Case_name,"Emitter_cost_avg")] = []
        for i in 1:T_OPT
            ETS_data[(Case_name, i)] = [] 
        end 
    end 
    return Dict_results, ETS_data
end

function preprocessing_marginal_cost_curves(LINEAR::Bool, LEARNING::Bool) # does not work yet
    # f denotes the function for the MRCC
    # f2 denotes the function for the MACC
    if LINEAR == true 
        f(remov, t, Cum_remov) = @.(A_r + B_r*remov)
        #f2(abat, t, Cum_abat) = @.(BETA*(EM_REF*(abat)).^GAMMA)
        f2(abat, t, Cum_abat) = @.(A_a+ B_a*abat)
    else 
    # NOTE: chaning the parameters of this function, will automatically be taken into account
        if LEARNING== true
            f3(remov, t, Cum_remov) = @.[(1100 + remov*100) (D_R ./ (1 + Cum_remov[t-1]).^(M_R) .*(remov./((U_R/(1+exp(K_R*(t+1)))) .- remov)).^C_R + m_R_add*(1-t/T_OPT) + r_R_add)]
            f(remov, t, Cum_remov) = minimum(f3(remov, t,  Cum_remov); dims=2)
            f1(abat, t,  Cum_abat) =  @.[(290+(abat)*110) (D_A ./ (1 +  Cum_abat[t-1]).^(M_A) .*( abat./((U_A/(1+exp(K_A*(t+1)))) .-  abat)).^C_A + m_A_add*(1-t/T_OPT) + r_A_add)]
            f2(abat, t,  Cum_abat) =  minimum(f1(abat, t,  Cum_abat); dims=2)
            #g(abat, t, Cum_abat) = @.(minimum([400.0.+abat, @.(D_A ./ (1 + Cum_abat[t-1]).^(M_A) .*(abat./((U_A/(1+exp(K_A*(t)))) .- abat)).^C_A + F_A_add*(1-t/T_OPT))]))
        elseif LEARNING == false
            f(remov, t, Cum_remov) = @.(D_R ./ (1 + 12000).^(M_R) .*(remov./((U_R/(1+exp(K_R*(T_OPT/2+1)))) .- remov)).^C_R + m_R_add*(1-1/3) + r_R_add)
            f1(abat, t,  Cum_abat) =  @.[(290+(abat)*110) (D_A ./ (1 +  (Cum_abat_init+9000)).^(M_A) .*( abat./((U_A/(1+exp(K_A*(T_OPT/2+1)))) .-  abat)).^C_A + m_A_add*(1-1/3) + r_A_add)]
            f2(abat, t,  Cum_abat) =  minimum(f1(abat, t,  Cum_abat); dims=2)
            #g(abat, t, Cum_abat) = @.(minimum([400.0.+abat, @.(D_A ./ (1 + Cum_abat[t-1]).^(M_A) .*(abat./((U_A/(1+exp(K_A*(t)))) .- abat)).^C_A + F_A_add*(1-t/T_OPT))]))
        end
    end
    return  f, f2
end



function process_sets_parameters!(model::Model, Parameters::DataFrame, Timeseries_data::DataFrame, T_OPT::Int64)
    "Function info: defining parameters and sets
    NOTE: find back the parameter explanation in the Parameters excel table "
    
    # extract the sets you need
    model.ext[:sets] = Dict()
    
    # create dictionary to store sets
    model.ext[:sets] = Dict()

    # define the sets
    T = model.ext[:sets][:T] = 0:1:T_OPT

    # generate a dictonary "parameters"
    model.ext[:parameters] = Dict()

    CAP =  model.ext[:parameters][:CAP] =  Dict(t =>  Timeseries_data[convert.(Int64,t+1), "ETS_cap [EUA/y]"] for t in T ) #.*0.48916
    CAPr = model.ext[:parameters][:CAPr] =  Dict(t =>  Timeseries_data[convert.(Int64,t+1), "RC_cap [EUA/y]"] for t in T )
    Budget = model.ext[:parameters][:Budget] = BUDGET
    #EUA_price = model.ext[:parameters][:EUA_price] = Dict(i => 100*(1 + Discountrate).^(i) for i in T) # just to check
    
    Beta =  model.ext[:parameters][:Beta] = BETA #Parameters[Index_names.=="Beta", "Value"][1] 
    Gamma = model.ext[:parameters][:Gamma] = GAMMA #Parameters[Index_names.=="Gamma", "Value"][1] 
    EM_ref =model.ext[:parameters][:EM_ref] = EM_REF #Parameters[Index_names.=="EM_ref", "Value"][1].*10^6 # FV: adjusted compared to original optimisation problem 
    FRAC_RO = model.ext[:parameters][:FRAC_RO] = Dict(t =>  Timeseries_data[convert.(Int64,t+1), "Frac_removal_obligation [%]"] for t in T )
    FRAC_RO_gross = model.ext[:parameters][:FRAC_RO_gross] = Dict(t =>  Timeseries_data[convert.(Int64,t+1), "Frac_removal_obligation_gross [%]"] for t in T )
    FRAC_EQ = model.ext[:parameters][:FRAC_EQ] = FRAC_eq

    EAC_price =  model.ext[:parameters][:EAC_price] =  Dict(t =>  EAC_start for t in T) # initial value of the EU ETS price = 0 
    CDR_price =  model.ext[:parameters][:CDR_price] =  Dict(t =>  CDR_start for t in T) # initial value of the EU ETS price = 0 
    #CCB_price = model.ext[:parameters][:CCB_price] = Dict(t =>  CCB_start for t in T) # initial value of the EU ETS price = 0 

    e = model.ext[:parameters][:e] = Dict(t => 0.0 for t in T) # number of emissions
    r = model.ext[:parameters][:r] = Dict(t =>  0.0  for t in T)  # number of removals
    Cum_remov = model.ext[:parameters][:Cum_remov] = Dict(t =>  0.0  for t in -1:1:T_OPT) 
    Cum_abat = model.ext[:parameters][:Cum_abat] = Dict(t =>  0.0  for t in -1:1:T_OPT) 
    Cum_abat[-1] = Cum_abat_init #MtCO2
    Cum_remov[-1] = Cum_remov_init #MtCO2


    x_bar = model.ext[:parameters][:x_bar] = Dict(t =>  0.0  for t in T) # total credits in EAC market
    x_e_bar = model.ext[:parameters][:x_e_bar] = Dict(t =>  0.0  for t in T) # total credits in EAC market
    ec_bar = model.ext[:parameters][:ec_bar] = Dict(t =>  0.0  for t in T) # total credits in EAC market
    rc_bar = model.ext[:parameters][:rc_bar] = Dict(t =>  0.0  for t in T) # total sold by CDR facility
    rc_e_bar = model.ext[:parameters][:rc_e_bar] = Dict(t =>  0.0  for t in T) # total emission credits bought by emitter
    rc_r_bar = model.ext[:parameters][:rc_r_bar] = Dict(t =>  0.0  for t in T) # total emission credits bought by emitter


    # ADMM rho value 
    ρ= model.ext[:parameters][:ρ] = RHO
    ρ1= model.ext[:parameters][:ρ1] = RHO
    ρ2= model.ext[:parameters][:ρ2] = RHO2
    # ADMM computational accuracy 
    
    epsilon = model.ext[:parameters][:epsilon] =  EPSILON
    DiscountFactor = model.ext[:parameters][:DiscountFactor] =  Dict(t =>  1/(1+DR)^t  for t in T) # Discount factor


    N = model.ext[:parameters][:N] = n
    # return model

    # MACC 


    # # high learning
    # D_a = model.ext[:parameters][:D_a]  = 1724
    # C_a = model.ext[:parameters][:C_a]  = 1.64
    # M_a = model.ext[:parameters][:M_a] = 0.621
    # K_a = model.ext[:parameters][:K_a] = -0.14
    # U_a = model.ext[:parameters][:U_a] = 0.99

    # # low learning
    # D_a = model.ext[:parameters][:D_a]  = 329
    # C_a = model.ext[:parameters][:C_a]  = 0.96
    # M_a = model.ext[:parameters][:M_a] = 0.075
    # K_a = model.ext[:parameters][:K_a] = -0.17
    # U_a = model.ext[:parameters][:U_a] = 0.98

    
    # customized
    # D_a = model.ext[:parameters][:D_a]  = 50
    # C_a = model.ext[:parameters][:C_a]  = 0.96
    # M_a = model.ext[:parameters][:M_a] = 0.00075
    # K_a = model.ext[:parameters][:K_a] = -20
    # U_a = model.ext[:parameters][:U_a] = 1.00

    D_a = model.ext[:parameters][:D_a]  = D_A # 50
    C_a = model.ext[:parameters][:C_a]  = C_A #0.96
    M_a = model.ext[:parameters][:M_a] = M_A  #0.05
    K_a = model.ext[:parameters][:K_a] = K_A #-5.0
    U_a = model.ext[:parameters][:U_a] = U_A #1.0

    # MRCC 
    # # high learning
    # D_r = model.ext[:parameters][:D_r]  = 7384 #7384
    # C_r = model.ext[:parameters][:C_r]  = 0.88
    # M_r = model.ext[:parameters][:M_r] = 0.621
    # K_r = model.ext[:parameters][:K_r] = -0.36
    # U_r = model.ext[:parameters][:U_r] =0.20

    # low learning 
    # D_r = model.ext[:parameters][:D_r]  = 1327 #7384
    # C_r = model.ext[:parameters][:C_r]  = 0.53
    # M_r = model.ext[:parameters][:M_r] = 0.075
    # K_r = model.ext[:parameters][:K_r] = -0.13
    # U_r = model.ext[:parameters][:U_r] =0.19

    # Custom learning
    D_r = model.ext[:parameters][:D_r]  = D_R #2500
    C_r = model.ext[:parameters][:C_r]  = C_R
    M_r = model.ext[:parameters][:M_r] = M_R #0.3
    K_r = model.ext[:parameters][:K_r] = K_R
    U_r = model.ext[:parameters][:U_r] =  U_R #0

    return model
end

function define_results!(Parameters::DataFrame,results::Dict,ADMM::Dict)


    Years = T_OPT + 1
    Iterations = ADMM_max_it #convert.(Int64, Parameters[Index_names.=="ADMM_max_it", "Value"][1])
    results["x"] = zeros(Iterations,Years) # number of credits which are no removals
    results["x_e"] = zeros(Iterations,Years) # number of credits which are no removals
    results["r"] = zeros(Iterations,Years) # number of removals
    results["e"] = zeros(Iterations,Years) # number of emission 
    results["ec"] = zeros(Iterations,Years) # number of emission certificates each year and each iteration
    results["rc"] = zeros(Iterations,Years) # number of removal certificates each year and each iteration (remover sold)
    results["rc_e"] = zeros(Iterations,Years) # number of removal certificates each year and each iteration (emitter bought)
    results["rc_r"] = zeros(Iterations,Years) # number of removal certificates each year and each iteration (emitter bought)
    results["Cum_remov"] = zeros(Iterations,Years) # number of removal certificates each year and each iteration (emitter bought)
    results["Cum_abat"] = zeros(Iterations,Years) # number of removal certificates each year and each iteration (emitter bought)


    results["EAC_price"] = EAC_start*ones(Iterations,Years) # allowance price
    results["CDR_price"] = CDR_start*ones(Iterations,Years)
 

    ADMM["Imbalances"] = Dict() 

    ADMM["Imbalances"]["EAC"] = zeros(Iterations,Years)
    ADMM["Imbalances"]["CDR"] = zeros(Iterations,Years)

    ADMM["Imbalances"]["Chi_EAC"] = zeros(Iterations,Years)
    ADMM["Imbalances"]["Chi_CDR"] = zeros(Iterations,Years)


    ADMM["Residuals"] = Dict()
    ADMM["Residuals"]["Primal_EAC"] =zeros(Iterations,1)
    ADMM["Residuals"]["Primal_CDR"] =zeros(Iterations,1)
    ADMM["Residuals"]["Dual_EC"] =zeros(Iterations,1)
    ADMM["Residuals"]["Dual_RC"] =zeros(Iterations,1)
    ADMM["Residuals"]["Dual_RC_E"]=zeros(Iterations,1)
    ADMM["Residuals"]["Dual_RC_R"]=zeros(Iterations,1)
    ADMM["Residuals"]["Dual_X"]=zeros(Iterations,1)
    ADMM["Residuals"]["Dual_X_E"]=zeros(Iterations,1)
    
    # ADMM["Tolerance"] = data["ADMM"]["epsilon"]*sqrt(2*data["nyears"])  
    ADMM["ρ"] = RHO*ones(ADMM_max_it,1) #Parameters[Index_names.=="rho", "Value"][1]*ones(Int64(Parameters[Index_names.=="ADMM_max_it", "Value"][1]),1)
    ADMM["ρ1"] = RHO*ones(ADMM_max_it,1) #Parameters[Index_names.=="rho", "Value"][1]*ones(Int64(Parameters[Index_names.=="ADMM_max_it", "Value"][1]),1)
    ADMM["ρ2"] = RHO2*ones(ADMM_max_it,1) #Parameters[Index_names.=="rho", "Value"][1]*ones(Int64(Parameters[Index_names.=="ADMM_max_it", "Value"][1]),1)
    ADMM["n_iter"] = 1 
    ADMM["walltime"] = 0

    return results, ADMM
end



function Area(f, t, Cum_remov, a, b, n)
    A = a+(b-a)/n/2
    xs = a:(b-a)/n:b
    deltas = diff(xs)
    cs = xs[1:end-1]
    sum(f(cs[i], t, Cum_remov) .* deltas[i] for i in 1:length(deltas))
end

function Area_abat(f2, t, Cum_abat, a, b, n)
    xs = a:(b-a)/n:b
    deltas = diff(xs)
    cs = xs[2:end]
    sum(f2(cs[i], t, Cum_abat) .* deltas[i] for i in 1:length(deltas))
end

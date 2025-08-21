function ADMM_price_based_m!(Case_name::String, results::Dict,ADMM::Dict, model_e::Model, model_r::Model, model_ccb::Model, dict::Dict, ETS_data::Dict, LINEAR::Bool) #, conditional_case::Bool)
    convergence = 0
    iterations = ProgressBar(1:ADMM_max_it-1)
    Cum_remov = model_r.ext[:parameters][:Cum_remov] 
    Cum_abat = model_e.ext[:parameters][:Cum_abat] 
    Beta = model_r.ext[:parameters][:Beta] 
    Gamma = model_r.ext[:parameters][:Gamma] 

    Remover!(model_r)
    CentralBank!(model_ccb)
    Emitter!(model_e)
    

    #AC(abat, t) = @. (D_a ./ (1 + Cum_abat[t-1]).^(M_a) .*(abat./((U_a/(1+exp(K_a*(t+1)))) .- abat)).^C_a)

    RC(remov, t, Cum_remov) = fr(remov, t, Cum_remov)
    AC(abat, t, Cum_abat) = fa(abat, t, Cum_abat)
 

    for iter in iterations
        if convergence == 0
            # Calculate penalty terms ADMM and update price to most recent value
            if iter > 1 # defining previous iteration values 

  
                model_e.ext[:parameters][:EAC_price] = Dict(t => results["EAC_price"][iter,t+1] for t in model_e.ext[:sets][:T])
                model_r.ext[:parameters][:EAC_price] = Dict(t => results["EAC_price"][iter,t+1] for t in model_e.ext[:sets][:T])
                model_ccb.ext[:parameters][:EAC_price] = Dict(t => results["EAC_price"][iter,t+1] for t in model_e.ext[:sets][:T]) #FV
                    
                model_e.ext[:parameters][:CDR_price] = Dict(t => results["CDR_price"][iter,t+1]  for t in model_e.ext[:sets][:T])
                model_r.ext[:parameters][:CDR_price] = Dict(t => results["CDR_price"][iter,t+1] for t in model_r.ext[:sets][:T])
                model_ccb.ext[:parameters][:CDR_price] = Dict(t => results["CDR_price"][iter,t+1] for t in model_r.ext[:sets][:T]) #FV

                model_e.ext[:parameters][:ec_bar] = Dict(t => results["ec"][iter-1,t+1] for t in model_e.ext[:sets][:T])
                model_r.ext[:parameters][:rc_bar] = Dict(t => results["rc"][iter-1,t+1] for t in model_e.ext[:sets][:T])
                model_ccb.ext[:parameters][:rc_e_bar] = Dict(t => results["rc_e"][iter-1,t+1] for t in model_e.ext[:sets][:T])
                model_ccb.ext[:parameters][:rc_r_bar] = Dict(t => results["rc_r"][iter-1,t+1] for t in model_e.ext[:sets][:T])


                model_e.ext[:parameters][:ρ] = ADMM["ρ1"][iter]
                model_r.ext[:parameters][:ρ] = ADMM["ρ2"][iter]
                model_ccb.ext[:parameters][:ρ1] = ADMM["ρ1"][iter]
                model_ccb.ext[:parameters][:ρ2] = ADMM["ρ2"][iter]
                    
                   
            end



            for t in  model_e.ext[:sets][:T]

                D_a = model_e.ext[:parameters][:D_a] 
                C_a = model_e.ext[:parameters][:C_a] 
                M_a = model_e.ext[:parameters][:M_a]
                K_a = model_e.ext[:parameters][:K_a] 
                U_a = model_e.ext[:parameters][:U_a] 
                

                #----------------------------

                    
                D_r = model_r.ext[:parameters][:D_r]  
                C_r = model_r.ext[:parameters][:C_r]  
                M_r = model_r.ext[:parameters][:M_r] 
                K_r = model_r.ext[:parameters][:K_r] 
                U_r = model_r.ext[:parameters][:U_r] 

                K_a = model_e.ext[:parameters][:K_a] 
                U_a = model_e.ext[:parameters][:U_a] 


                st = 0.0001

                if LINEAR == true
                    rem_max = 1.0 
                    ab_max = 1.0
                else 
                    if LEARNING == true
                        rem_max = U_r./(1+exp(K_r*(t)))
                        ab_max  = U_a./(1+exp(K_a*(t)))
                    elseif LEARNING == false
                        rem_max = U_r./(1+exp(K_r*(T_OPT/2)))
                        ab_max  = U_a./(1+exp(K_a*(T_OPT/2)))
                    end
                end
                #-----------------------------------
                abat = collect(0:st:ab_max)
                EAC_price_nom = maximum([0,model_e.ext[:parameters][:EAC_price][t]]) # FV: maximum maybe not correct here?
                LAMDA_ER = @.(EAC_price_nom)
                _, j = findmin(abs.(AC(abat,t, Cum_abat).-LAMDA_ER)) #FV_add: t+1
                Abat_t = maximum([0, abat[j]]).*EM_REF

                model_e.ext[:parameters][:e][t] = EM_REF - maximum([0,Abat_t])

                Cum_abat[t] = model_e.ext[:parameters][:Cum_abat][t] = model_e.ext[:parameters][:Cum_abat][t-1] + Abat_t
                
                results["Cum_abat"][iter, t+1] = Cum_abat[t]
                #----------------------------------------
                remov = collect(0:st:rem_max)
                CDR_price_nom = maximum([0,model_r.ext[:parameters][:CDR_price][t]]) # FV: maximum maybe not correct here?
                LAMDA_CDR = @.(CDR_price_nom)
                _, j = findmin(abs.(RC(remov,t+1, Cum_remov).-LAMDA_CDR)) #FV_add: t+1
                Remov_t = maximum([0, remov[j]]).*EM_REF

                model_r.ext[:parameters][:r][t] = maximum([0,Remov_t])
                
                Cum_remov[t] = model_r.ext[:parameters][:Cum_remov][t] = model_r.ext[:parameters][:Cum_remov][t-1] + model_r.ext[:parameters][:r][t]

                results["Cum_remov"][iter, t+1] = Cum_remov[t]
            end

            results["e"][iter,:] = [model_e.ext[:parameters][:e][t] for t in model_e.ext[:sets][:T]]
            results["r"][iter,:] = [model_r.ext[:parameters][:r][t] for t in model_e.ext[:sets][:T]] # NOTE: multiplied with FRAC_eq
            


            Update_emitter!(model_e)
            variables_extract(model_e)

            results["ec"][iter,:] = [value.(model_e.ext[:variables][:ec][t]) for t in model_e.ext[:sets][:T]]
            model_ccb.ext[:parameters][:ec_bar] = Dict(t => results["ec"][iter,t+1] for t in model_e.ext[:sets][:T] )

            Update_remover!(model_r)
            variables_extract(model_r)
           
            results["rc"][iter,:] = [value.(model_r.ext[:variables][:rc][t]) for t in model_r.ext[:sets][:T]]
            model_ccb.ext[:parameters][:rc_bar]  = Dict(t => results["rc"][iter,t+1] for t in model_r.ext[:sets][:T] )
            

            #Solve updated problem
            Update_CentralBank!(model_ccb)
            variables_extract(model_ccb)
            results["rc_e"][iter,:] = [value.(model_ccb.ext[:variables][:rc_e][t]) for t in model_ccb.ext[:sets][:T]]
            model_e.ext[:parameters][:rc_e_bar] = Dict(t => results["rc_e"][iter,t+1] for t in model_r.ext[:sets][:T] )
            results["rc_r"][iter,:] = [value.(model_ccb.ext[:variables][:rc_r][t]) for t in model_ccb.ext[:sets][:T]]
            model_r.ext[:parameters][:rc_r_bar] = Dict(t => results["rc_r"][iter,t+1] for t in model_r.ext[:sets][:T] )
            



              
            # # Imbalances
            ADMM["Imbalances"]["EAC"][iter,:] = [model_e.ext[:parameters][:CAP][t] for t in  model_e.ext[:sets][:T]]  + results["rc_e"][iter,:].*FRAC_EQ - results["ec"][iter,:] 
            ADMM["Imbalances"]["CDR"][iter,:] = [value.(model_r.ext[:variables][:rc][t]) - value.(model_ccb.ext[:variables][:rc_r][t]) for t in model_r.ext[:sets][:T]]

            # # Primal residuals
            ADMM["Residuals"]["Primal_EAC"][iter] = sqrt(sum(ADMM["Imbalances"]["EAC"][iter,:].^2))
            ADMM["Residuals"]["Primal_CDR"][iter] = sqrt(sum(ADMM["Imbalances"]["CDR"][iter,:].^2))

            if iter > 1 
                ADMM["Residuals"]["Dual_EC"][iter] = sqrt(sum(( (results["ec"][iter,:] - results["ec"][iter-1,:])).^2))
                ADMM["Residuals"]["Dual_RC"][iter] = sqrt(sum(( (results["rc"][iter,:] - results["rc"][iter-1,:])).^2))
                ADMM["Residuals"]["Dual_RC_E"][iter] = sqrt(sum(( (results["rc_e"][iter,:] - results["rc_e"][iter-1,:])).^2))
                ADMM["Residuals"]["Dual_RC_R"][iter] = sqrt(sum(( (results["rc_r"][iter,:] - results["rc_r"][iter-1,:])).^2))
            else
            end
           
            # Price updates 
            results["EAC_price"][iter+1,:] = results["EAC_price"][iter,:] - model_e.ext[:parameters][:ρ] .*ADMM["Imbalances"]["EAC"][iter,:]
            results["CDR_price"][iter+1,:] = results["CDR_price"][iter,:] - model_r.ext[:parameters][:ρ] .*ADMM["Imbalances"]["CDR"][iter,:] 

            # Progress bar
            set_description(iterations, string(@sprintf("Primal EAC %.3f -- Primal CDR %.3f -- Dual EC %.3f -- Dual RC %.3f -- Dual RC E %.3f -- Dual RC R %.3f",  ADMM["Residuals"]["Primal_EAC"][iter], ADMM["Residuals"]["Primal_CDR"][iter], ADMM["Residuals"]["Dual_EC"][iter], ADMM["Residuals"]["Dual_RC"][iter], ADMM["Residuals"]["Dual_RC_E"][iter], ADMM["Residuals"]["Dual_RC_R"][iter])))

            # Check convergence: primal and dual satisfy tolerance OR primal remains high and dual near zero. 
            # In case the latter happens, the solution needs to be inspected to ensure this can be interpreted as an equilibrium.
            if (ADMM["Residuals"]["Primal_EAC"][iter] <= model_e.ext[:parameters][:epsilon] && ADMM["Residuals"]["Primal_CDR"][iter] <= model_e.ext[:parameters][:epsilon] && ADMM["Residuals"]["Dual_EC"][iter] <= model_e.ext[:parameters][:epsilon] && ADMM["Residuals"]["Dual_RC"][iter] <= model_e.ext[:parameters][:epsilon] &&  ADMM["Residuals"]["Dual_RC_E"][iter] <= model_e.ext[:parameters][:epsilon] &&  ADMM["Residuals"]["Dual_RC_R"][iter] <= model_e.ext[:parameters][:epsilon])
                convergence = 1
            end
            ADMM["n_iter"] = copy(iter)
        end
    end
    # Actual amount of removals 
    results["r"][ADMM["n_iter"],:] = results["r"][ADMM["n_iter"],:].*FRAC_R #NOTE: multiplied with fraciton in reality in the general case we say R= EQ, not for the sensitivity analysis, EQ fixed at 0.85, R can vary

    Cum_remov_extract = Dict( i =>  results["Cum_remov"][ADMM["n_iter"],i+1] for i in 0:T_OPT)
    Cum_remov_extract[-1] = Cum_remov_init # creation of element in case that t = 0
    Cum_abat_extract = Dict( i =>  results["Cum_abat"][ADMM["n_iter"],i+1] for i in 0:T_OPT)
    Cum_abat_extract[-1] = Cum_abat_init # creation of element in case that t = 0
    
    Cum_NoLearning_remov = Dict( i =>  Cum_remov_init for i in -1:T_OPT)
    Cum_NoLearning_abat = Dict( i =>  Cum_abat_init for i in -1:T_OPT)


    for t in 0:T_OPT
        if value.(model_r.ext[:parameters][:r][t]) == 0.0 
            # model_r.ext[:parameters][:CDR_price][t] = convert(Vector{Float64}, RC(0.0,t, Cum_remov))
            results["CDR_price"][ADMM["n_iter"],t+1] = f(0.0, t,  Cum_remov_extract)[1]
        else
        end
    end


    parameters_extract(model_r)
    Total_costs_rem = sum(1/(1+DRc)^t*quadgk(var -> fr(var, t, Cum_remov_extract), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r.ext[:sets][:T])[1]
    Total_costs_rem_NL = sum(1/(1+DRc)^t*quadgk(var -> fr(var, t, Cum_NoLearning_remov), 0, r[t]./EM_REF)[1].*EM_REF for t in model_r.ext[:sets][:T])[1]
    parameters_extract(model_e)
    Total_costs_abat = sum(1/(1+DRc)^t*quadgk(var -> fa(var, t, Cum_abat_extract), 0,(EM_REF-e[t])./EM_REF)[1].*EM_REF for t in model_r.ext[:sets][:T])[1]
    Total_costs_abat_NL = sum(1/(1+DRc)^t*quadgk(var -> fa(var, t, Cum_NoLearning_abat), 0,(EM_REF-e[t])./EM_REF)[1].*EM_REF for t in model_r.ext[:sets][:T])[1]

    EAC_price = model_e.ext[:parameters][:EAC_price]
    CDR_price = model_r.ext[:parameters][:CDR_price] # FV: maybe this one should be changed as well. 


    Total_costs = Total_costs_rem + Total_costs_abat
    
   
    Total_budget =   sum([CAP[t].*values(EAC_price[t]) for t in 0:1:T_OPT]) #sum((results["ec"][ADMM["n_iter"],:].*[values(EAC_price[t]) for t in 0:1:T_OPT] -results["rc"][ADMM["n_iter"],:].*[values(CDR_price[t]) for t in 0:1:T_OPT])) 
    Total_net_emissions = sum(results["e"][ADMM["n_iter"],:] - results["r"][ADMM["n_iter"],:])
    Total_emissions  = sum(results["e"][ADMM["n_iter"],:]) 
    Total_removals = sum(results["r"][ADMM["n_iter"],:])
    Total_abated =  (T_OPT+1)*EM_REF - Total_emissions
    ETS_max = maximum([EAC_price[t] for t in 0:1:T_OPT])
    RC_max = maximum([results["CDR_price"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])
    ETS_avg = sum([EAC_price[t].*results["ec"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])./sum([results["e"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT])
    credit_costs_emitter_avg = sum([EAC_price[t].*results["ec"][ADMM["n_iter"],t+1]  for t in 0:1:T_OPT])./sum([results["e"][ADMM["n_iter"],t+1] for t in 0:1:T_OPT]) # NOTE: ec also includes rc 




    push!(dict[(Case_name,"Total_net_emissions")], Total_net_emissions/1000) #Gton
    push!(dict[(Case_name,"Total_emissions")], Total_emissions/1000) #Gton
    push!(dict[(Case_name,"Fraction_removals")], Total_net_emissions/Total_emissions)
    push!(dict[(Case_name,"Total_costs_per_net")], Total_costs/Total_net_emissions)
    push!(dict[(Case_name,"Total_costs_per_abat")], Total_costs/Total_abated)
    push!(dict[(Case_name,"Total_budget_per_net")], Total_budget/Total_net_emissions)
    push!(dict[(Case_name,"Learning_gains_per_rem")], (Total_costs_rem_NL -Total_costs_rem)/(Total_removals./FRAC_R))
    push!(dict[(Case_name,"Learning_gains_per_em")], (Total_costs_abat_NL- Total_costs_abat)/Total_abated)
    push!(dict[(Case_name,"ETS_max")], ETS_max)
    push!(dict[(Case_name,"ETS_avg")], ETS_avg)
    push!(dict[(Case_name,"RC_max")], RC_max)
    push!(dict[(Case_name,"Fraction_costs_abat")], Total_costs_abat/Total_costs)
    push!(dict[(Case_name,"Fraction_costs_rem")], Total_costs_rem/Total_costs)
    push!(dict[(Case_name,"Emitter_cost_avg")], credit_costs_emitter_avg + Total_costs_abat/Total_emissions)

    for t in 1:T_OPT
        push!(ETS_data[(Case_name, t)], values(EAC_price[t]))
    end

    return dict 
end

function Emitter!(model::Model)
    model.ext[:variables] = Dict()
    model.ext[:constraints] = Dict()

    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    e = model.ext[:parameters][:e]                           # Emissions emitted and not abated each year
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term

    ρ = model.ext[:parameters][:ρ] # Parameter controlling the price update step size in the ADMM algorithm
    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]
    FRAC_EQ = model.ext[:parameters][:FRAC_EQ]



    # Define variables
    #e = model.ext[:variables][:e] = @variable(model, e[T], lower_bound=0, base_name="Emissions bought in year t by industrial fringe")
    ec = model.ext[:variables][:ec] = @variable(model, ec[T], lower_bound=0, base_name="Emission allowances bought in year t by industrial fringe")

    # Constraints 
    Banking_e =  model.ext[:constraints][:Banking_e] = @constraint(model, [t=T],
    (sum(ec[i] for  i in 0:t)  >= sum(e[i]   for i in 0:t)))
    

    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    sum(DiscountFactor[t]*(ec[t].*EAC_price[t]) for t in T) +  
    # Penalty term ETS market
    ( ρ./2) * sum( DiscountFactor[t]*((ec[t] - (ec_bar[t] + 1/N*(CAP[t] + rc_e_bar[t].*FRAC_EQ -  ec[t])))^2) for t in T) .* a  )
    
    
    optimize!(model);





    return model
end



function Update_emitter!(model::Model)
    model.ext[:variables] = Dict()
    # Extract sets
    T = model.ext[:sets][:T] # set of years

    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    e = model.ext[:parameters][:e]                           # Emissions emitted and not abated each year
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term

    ρ = model.ext[:parameters][:ρ] # Parameter controlling the price update step size in the ADMM algorithm
    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]
    FRAC_EQ = model.ext[:parameters][:FRAC_EQ]

    # for t in T
    #     if EAC_price[t] < P_CAP
    #         rc_e_bar[t] = model.ext[:parameters][:rc_e_bar][t] = 0
    #     else 
    #         skip 
    #     end
    # end
    # Define variables
   # e = model.ext[:variables][:e] = @variable(model, [t=T], lower_bound=0, base_name="Emissions in year t by industrial fringe")
    ec = model.ext[:variables][:ec] = @variable(model, [t=T], lower_bound=0, base_name="Emission allowances bought in year t by industrial fringe")

    for t in T
        delete(model,model.ext[:constraints][:Banking_e][t])
    end




    # Constraints 
    Banking_e =  model.ext[:constraints][:Banking_e] = @constraint(model, [t=T],
    (sum(ec[i] for  i in 0:t)  >= sum(e[i]   for i in 0:t)))
    

    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    sum(DiscountFactor[t]*(ec[t].*EAC_price[t]) for t in T) +  # sum(DiscountFactor[t]*BETA*(EM_REF - e[t]).^(2)/2 for t in T) +
    # Penalty term ETS market
    ( ρ./2) * sum( DiscountFactor[t]*((ec[t] - (ec_bar[t] + 1/N*(CAP[t] + rc_e_bar[t].*FRAC_EQ -  ec[t])))^2) for t in T) .* a  )
    
    
    optimize!(model);

    return model 
end




function Remover!(model::Model)
    model.ext[:variables] = Dict()
    model.ext[:constraints] = Dict()


    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    r = model.ext[:parameters][:r]                           # Removals removed
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term

    ρ = model.ext[:parameters][:ρ] # Parameter controlling the price update step size in the ADMM algorithm
    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]

    # Define variables


    #r = model.ext[:variables][:r] = @variable(model, r[T], lower_bound=0, base_name="CDR removed in year t by remover")
    rc = model.ext[:variables][:rc] = @variable(model, rc[T], lower_bound=0, base_name="CDR removed sold in year t by remover")


    
    # Constraints 
    Banking_r =  model.ext[:constraints][:Banking_r] = @constraint(model, [t=T],
    (sum(rc[i] for  i in 0:t)  <= sum(r[i]   for i in 0:t)))
    
    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    - sum(DiscountFactor[t]*(rc[t].*CDR_price[t]) for t in T)  + #sum(DiscountFactor[t]*(A_r*r[t] + B_r*r[t]^(2)/2/3000) for t in T) +
    # Penalty term CDR market
    ( ρ./2) * sum( DiscountFactor[t]*(((rc[t]) -(rc_bar[t] + 1/N*(rc_r_bar[t] -  rc[t])))^2) for t in T) .* a )
    
    
    
    optimize!(model);




    return model
end






function Update_remover!(model::Model)
    model.ext[:variables] = Dict()


    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    r = model.ext[:parameters][:r]                           # Removals removed
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term

    ρ = model.ext[:parameters][:ρ] # Parameter controlling the price update step size in the ADMM algorithm
    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]

    # Define variables

    rc = model.ext[:variables][:rc] = @variable(model, [t=T], lower_bound=0, base_name="CDR allowances sold in year t by remover")
    #r = model.ext[:variables][:r] = @variable(model, [t=T], lower_bound=0, base_name="CDR removed sold in year t by remover")


    for t in T
        delete(model,model.ext[:constraints][:Banking_r][t])
     
    end



    # Constraints 
    Banking_r =  model.ext[:constraints][:Banking_r] = @constraint(model, [t=T],
    (sum(rc[i] for  i in 0:t)  <= sum(r[i]   for i in 0:t)))
    
    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    - sum(DiscountFactor[t]*(rc[t].*CDR_price[t] ) for t in T)  +  #sum(DiscountFactor[t]*(A_r*r[t] + B_r*r[t]^(2)/2/3000) for t in T) +
    # Penalty term CDR market
    ( ρ./2) * sum( DiscountFactor[t]*(((rc[t]) -(rc_bar[t] + 1/N*(rc_r_bar[t] -  rc[t])))^2) for t in T) .* a )
    
    
    optimize!(model);


    return model
end




function CentralBank!(model::Model)

    model.ext[:variables] = Dict()
    model.ext[:constraints] = Dict()


    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    r = model.ext[:parameters][:r]                           # Removals removed
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term

    ρ1= model.ext[:parameters][:ρ1] 
    ρ2= model.ext[:parameters][:ρ2] 

    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]
    FRAC_EQ = model.ext[:parameters][:FRAC_EQ]
    # Define variables


    rc_e = model.ext[:variables][:rc_e] = @variable(model, rc_e[T], lower_bound=0, base_name="CDR sold to emitter")
    rc_r = model.ext[:variables][:rc_r] = @variable(model, rc_r[T], lower_bound=0, base_name="CDR bought from remover")

    
    # Constraints 
    Banking_ccb =  model.ext[:constraints][:Banking_ccb] = @constraint(model, [t=T],
    (sum(rc_e[i] for  i in 0:t)  <= sum(rc_r[i]   for i in 0:t)))

    model.ext[:objective] = @objective(model, Min, 
    # MACC
    (- sum(DiscountFactor[t]*rc_e[t].*((EAC_price[t])) for t in T)) + sum(DiscountFactor[t]*(rc_r[t].*(CDR_price[t])) for t in T)  +
    # Penalty term ETS market
    ( ρ1./2) * sum( DiscountFactor[t]*((rc_e[t] - (rc_e_bar[t] - 1/N*(CAP[t] + rc_e[t].*FRAC_EQ -  ec_bar[t])))^2) for t in T) .* a   +
    ( ρ2./2) * sum( DiscountFactor[t]*(((rc_r[t]) -(rc_r_bar[t] - 1/N*(rc_r[t] -  rc_bar[t])))^2) for t in T) .* a 
    )
    
    
    optimize!(model);


    return model
end



function Update_CentralBank!(model::Model)

    model.ext[:variables] = Dict()


    # Extract sets
    T = model.ext[:sets][:T] # set of years
    
    # Extract parameters
    DiscountFactor = model.ext[:parameters][:DiscountFactor] # Discount factor
    r = model.ext[:parameters][:r]                           # Removals removed
    EAC_price = model.ext[:parameters][:EAC_price]           # EUA price 
    CDR_price = model.ext[:parameters][:CDR_price]           # CDR price



    ec_bar = model.ext[:parameters][:ec_bar]        # ADMM penalty term
    rc_bar = model.ext[:parameters][:rc_bar]        # ADMM penalty term
    rc_e_bar = model.ext[:parameters][:rc_e_bar]    # ADMM penalty term
    rc_r_bar = model.ext[:parameters][:rc_r_bar]    # ADMM penalty term
    ρ1= model.ext[:parameters][:ρ1] 
    ρ2= model.ext[:parameters][:ρ2] 
    CAP = model.ext[:parameters][:CAP] # cap 
    N = model.ext[:parameters][:N] # number of participants in the ETS market 
    Budget = model.ext[:parameters][:Budget]
    FRAC_EQ = model.ext[:parameters][:FRAC_EQ]
    # Define variables


    rc_e = model.ext[:variables][:rc_e] = @variable(model, [t=T], lower_bound=0, base_name="CDR sold to emitter")
    rc_r = model.ext[:variables][:rc_r] = @variable(model, [t=T], lower_bound=0, base_name="CDR bought from remover")

    for t in T
        delete(model, model.ext[:constraints][:Banking_ccb][t])
    end
    # delete(model, model.ext[:constraints][:End_ccb][T_OPT])
    # Constraints 
    Banking_ccb =  model.ext[:constraints][:Banking_ccb] = @constraint(model, [t=T],
    (sum(rc_e[i] for  i in 0:t)  <= sum(rc_r[i]   for i in 0:t)))


    
    # Objective 
    
    model.ext[:objective] = @objective(model, Min, 
    # MACC
    (- sum(DiscountFactor[t]*rc_e[t].*((EAC_price[t])) for t in T)) +  sum(DiscountFactor[t]*(rc_r[t].*(CDR_price[t])) for t in T)  +
    # Penalty term ETS market
    ( ρ1./2) * sum( DiscountFactor[t]*((rc_e[t] - (rc_e_bar[t] - 1/N*(CAP[t] + rc_e[t].*FRAC_EQ -  ec_bar[t])))^2) for t in T) .* a   +
    # ( ρ1./2) * sum( DiscountFactor[t]*(((1/N*(CAP[t] + rc_e[t] -  ec_bar[t])))^2) for t in T) .* a   +
    ( ρ2./2) * sum( DiscountFactor[t]*(((rc_r[t]) -(rc_r_bar[t] - 1/N*(rc_r[t] -  rc_bar[t])))^2) for t in T) .* a 
    )
    
    
    optimize!(model);

    
    return model
end







function plot_model!(Case_name::String, C_nr, model_e, model_r, model_ccb, ADMM, results, Dict_results::Dict, plotting::Bool)
    
    if !plotting
        skip
    else



        parameters_extract(model_e)
        parameters_extract(model_r)


        # Step 6: 
        EM = results["e"][ADMM["n_iter"],:]
        REM = results["r"][ADMM["n_iter"],:]
        hatch_vect = [x for p in zip(EM, zeros(length(EM))) for x in p].*[0; repeat([1, 1, 0, 0], convert(Int64, T_OPT/2)); 1]
        Doubled_time_vector =  [x for x in Timeseries_data[1:1:(T_OPT+1),"Year"] for _ in 1:2] + [repeat([0, 0, 0, -1], convert(Int64, T_OPT/2)); 0; 0]

        Peac = results["EAC_price"][ADMM["n_iter"],:]
        Pcdr = results["CDR_price"][ADMM["n_iter"],:]

        variables_extract(model_e)
        variables_extract(model_r)
        variables_extract(model_ccb)
        
        # parameters_extract(model)
        # variables_extract(model)


        EAC_surplus = [sum( ec[i]  - (EM[i+1])  for i in 0:t) for t in model_r.ext[:sets][:T]]
        CDR_surplus = [sum((REM[i+1] -rc[i])  for i in 0:t) for t in model_r.ext[:sets][:T]]
        CCB_surplus = [sum((rc_r[i] -rc_e[i])  for i in 0:t) for t in model_r.ext[:sets][:T]]


        ymax = 650.0


        Color = ["snow3" "snow2"  "slategray3"]
        Alpha = [1 0.4 0.5]

        pp = stackedarea(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM -REM],   ylims = (-200, 3000),  yticks =0:500:3000, labels = ["Emissions"  "Removals"], foreground_color_legend = nothing, title = "", ylabel = L"Emissions or Removals [MtCO$_2$/yr]", xlabel = "Years", legend = (0.02, 0.9), linecolor = Color, color = Color, alpha = Alpha, right_margin = 5mm)    
        plot!(Doubled_time_vector, [hatch_vect], ylims = (-200, 3000), labels = "", title = "",  linecolor = :white, alpha = 0.4, lw=1)  
        plot!(Timeseries_data[1:1:(T_OPT+1),"Year"], [(EM-REM)], color =:grey, labels="" )
        plot!(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM], color =:grey, labels="" )
        plot!(twiny(pp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], lw=2, ylims = (-200, 3000), label = ["ETS cap"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "navy" "black" "brown"], ls = [:solid :dash :dash :dash], margin = 0mm)
        plot!(twinx(pp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac[t]   for t in 1:1:(T_OPT+1)] [Pcdr[t]  for t in 1:1:(T_OPT+1)]], ylims = (0, ymax), yticks = 50:50:ymax, label = ["ETS price (MAC)" "RC price (MRC)"], foreground_color_legend = nothing, color = ["grey" "grey83"], ylabel = [L"Certificate price [€/tCO$_2$]"], legend = (0.58,0.98), lw = 2, margin = 0mm)  
        pprint = plot(pp, size  = (500, 400),  xlim = [2024,(2024+T_OPT)], titlefont = font(11),  title = Plot_title[C_nr], plot_titlefontsize = 12, topmargin = 2mm)
        annotate!(2032, 2570,text(L"\textbf{/////}", :white, :right, 6), linewidth =2)
        annotate!(2032, 2475,text("___", :grey, :right, 8), linewidth =2)
        annotate!(2032, 2390,text("___", :grey, :right, 8), linewidth =2)
        



        ppp = groupedbar(Timeseries_data[1:1:(T_OPT+1),"Year"], [EM -REM], bar_position =:stack, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), yticks =-500:250:2500, labels = ["Emissions"  "Removals"], foreground_color_legend = nothing, title = "", ylabel = L"Emissions or Removals [MtCO$_2$/yr]", xlabel = "Years", legend = (0.02, 0.85), lw = 0.75, color = Color, linecolor = [:gray31 :gray31], alpha = Alpha, right_margin = 5mm)
        plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [Timeseries_data[1:1:(T_OPT+1), "ETS_cap [EUA/y]"]./1], lw=[1 1], xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["ETS cap"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.98), color = ["black"  "black" "black" "brown"], ls = [:solid :dash :dash :dash], margin = 0mm)
        plot!(twiny(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [EM-REM], lw=0, xticks = (2025:5:(2024+T_OPT)), ylims = (-600, 2500), label = ["Net emissions"],foreground_color_legend = nothing, xformatter=_->"", yformatter=_->"", xlabel = "", ylabel="", legend = (0.02,0.91), color = ["black"  "black" "black" "brown"], ls = [:dash :dash :dash :dash], margin = 0mm, markershape = :xcross, markersize = 2)
        plot!(twinx(ppp), Timeseries_data[1:1:(T_OPT+1),"Year"], [[Peac[t]   for t in 1:1:(T_OPT+1)] [Pcdr[t].*FRAC_EQ  for t in 1:1:(T_OPT+1)] [Pcdr[t]  for t in 1:1:(T_OPT+1)]], ylims = (0, ymax), yticks = 0:50:ymax, background_color_legend = nothing, label = ["ETS price (MAC)" "RC price (net) (MRC)" "RC price (gross)"], foreground_color_legend = nothing, color = ["grey" "grey83" "grey83"], ls = [:solid :solid :dash ], ylabel = [L"Certificate price [€/tCO$_2$]"], legend = (0.5,0.98), lw = 2, margin = 0mm)  

        pprint = plot(ppp, size  = (500, 400),  xlim = [2023,(T_OPT_plot)],  titlefont = font(11),  title = Plot_title[C_nr], plot_titlefontsize = 12, topmargin = 2mm, gridlinewidth = 1, xtickfontsize=12,ytickfontsize=12, xguidefontsize=14,yguidefontsize=14,legendfontsize=12)
        
        filename_string  =  "C:/Users/VERBISTF/OneDrive - KU Leuven/PhD Flore/P4_negative_emissions/Figures/Plots/$(Case_name).svg"
        savefig(filename_string) 
        savefig(pprint, filename_string)


        return display(pprint)
    end
end

